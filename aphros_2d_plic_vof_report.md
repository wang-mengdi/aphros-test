# Aphros 2D PLIC VOF 代码结构分析报告

## 1. 项目概述

Aphros 是一个基于 C++17 开发的维度通用（1D-4D）计算流体力学代码库。PLIC（Piecewise Linear Interface Calculation）VOF（Volume of Fluid）是其中处理多相流界面捕捉的核心模块。

- **语言**: C++17，大量使用模板和 `.ipp` 文件分离模板实现
- **网格类型**: `MeshCartesian<double, dim>`，通过 `MULTIDIMX` 宏实例化 1-4 维

---

## 2. 编译方案

### 主构建系统：CMake
- 根 `CMakeLists.txt` 位于 `src/CMakeLists.txt`
- 使用 `add_object(T name)` 自定义函数，为每个模块创建 **OBJECT 库**（`name_obj`）和对应的 **INTERFACE 库**（`name`）
- 所有目标文件最终链接为共享库 `libaphros.so`（`src/CMakeLists.txt:361-365`）
- 对外 C 接口封装为 `libaphros_c.so`，并链接为可执行文件 `ap.mfer`（`src/CMakeLists.txt:382-388`）

### 传统构建系统：GNU Make
- `make/bootstrap` 脚本自动扫描 `src/` 下所有 `.c/.cpp/.h/.hpp/.inc/.ipp` 文件，生成 `.mk` 依赖规则
- 同时生成 Unix 和 Windows 两种规则集（`make/unix/` 和 `make/windows/`）
- `make/0*.mk` 和 `make/1*.mk` 文件控制功能开关（如 `0avx.mk` 留空表示禁用，对应 `1avx.mk` 启用）

### 编译选项
| 选项 | 功能 |
|------|------|
| `USE_AVX` (ON) | 启用 AVX2/FMA 向量化优化 Youngs 法向计算 |
| `USE_MPI` (ON) | 分布式并行计算 |
| `USE_OPENMP` (ON) | 共享内存并行 |
| `USE_DIM2` (ON) | 编译 2D 几何 |
| `USE_DIM3` (ON) | 编译 3D 几何 |
| `USE_DIM1` / `USE_DIM4` (OFF) | 1D / 4D 几何 |
| `USE_HDF` (ON) | HDF5 并行输出 |
| `USE_HYPRE` (ON) | Hypre 线性求解器 |
| `USE_BACKEND_CUBISM` (ON) | Cubism 分布式后端 |
| `USE_BACKEND_LOCAL` (ON) | 本地共享内存后端 |
| `USE_BACKEND_NATIVE` (ON) | Native 分布式后端 |

---

## 3. PLIC VOF 核心文件结构

### 3.1 重构数学核心：`Reconst<Scal>`

**文件**: `src/solver/reconst.h`

纯静态方法类，实现了 PLIC 重构的全部几何数学。维度通用（通过 `Vect` 模板参数实现 2D/3D 统一）。

| 函数 | 功能 |
|------|------|
| `GetLineU(n, a, h)` | 给定法向 `n` 和平面常数 `a`，计算体积分数 |
| `GetLineA(n, u, h)` | 逆运算：由法向和体积分数求平面常数（迭代法） |
| `GetLineU0` / `GetLineA0` | 有序法向版本，使用解析分段公式（涉及三次方程求解 `SolveCubic`） |
| `GetLineVol(n, a, h, dx, d)` | 计算 advection 位移 `dx` 后进入下游单元的体积盈余 |
| `GetLineVolStr(...)` | Lagrange Explicit 步骤的拉伸界面版本 |
| `GetLineFlux(n, a, h, q, dt, d)` | 计算 Eulerian 几何通量（面通量） |
| `GetLineFluxStr(...)` | 拉伸通量版本 |
| `GetCutPoly(xc, n, a, h)` | 返回 PLIC 平面切割单元的 polygon（可视化用） |
| `GetCenter(n, a, h)` | 切割 polygon 的质心 |
| `SolveCubic(a,b,c,d,k)` | 三次方程求解器 |
| `GetFitN(xx)` | 由点云拟合最优平面法向 |
| `GetInterPoly`, `GetInterLine` | Polygon 与平面/直线的相交工具 |

### 3.2 法向计算：`UNormal<M>`

**文件**: `src/solver/normal.h`, `src/solver/normal.ipp`, `src/solver/normal.cpp`

| 方法 | 算法 |
|------|------|
| `CalcNormalYoungs` / `CalcNormalYoungs1` / `CalcNormalYoungsAvx` | **Youngs 格式**：体积分数梯度（y 方向负梯度求法向）。3x3x3 模板，节点插值后平均回单元。`Avx` 版本用 `__m256d` 一次处理 4 个单元 |
| `CalcNormalHeight` / `CalcNormalHeight1` | **高度函数法**：沿各轴向累加 3 列单元求和得高度，选最陡方向，L1 归一化输出 |
| `CalcNormal` | **组合**：先 Youngs，再用高度函数覆盖（当高度函数在主导方向更陡时） |

### 3.3 单相 VOF 求解器：`Vof<EB>`

**文件**: `src/solver/vof.h`, `src/solver/vof.ipp`, `src/solver/vof.cpp`

继承自 `AdvectionSolver<M>`，后者继承自 `UnsteadyIterativeSolver` -> `UnsteadySolver`。

**时间步流程**（`Vof<EB>::Imp` -> `MakeIteration`）：

```
StartStep → 配置边界条件，t=0时初始重构
  │
  ▼
MakeIteration:
  ├─ ① 初始化: fcu_.iter_curr = fcu_.time_prev + dt·source·(1-u)
  │              fcuu_ = sharp(阈值0.5) 场 (Weymouth 散度项)
  │
  ├─ ② ReconstPlanes (重构):
  │     a. 检测界面单元 fci_: 0 < u < 1
  │     b. 计算法向 fcn_ = CalcNormal (Youngs + 高度函数)
  │     c. 计算平面常数 fca_ = GetLineA(n, u, h)
  │
  ├─ ③ 方向分裂 Advection Sweeps (方向顺序每步旋转):
  │     对每个 sweep 方向:
  │       ├─ interfacial 单元: GetLineFlux(n,α,h,flux,dt,d) 精确几何通量
  │       ├─ pure 单元: flux = mixture_flux × u
  │       └─ 传播颜色至下游
  │
  │     更新策略 (scheme 参数):
  │       • plain:   u += -dl
  │       • aulisa:  EI: u = (u-dl)/(1-ds); LE: u += u*ds - dl
  │       • weymouth: u += fcuu_*ds - dl  (默认)
  │     clip to [0,1], 边界通信
  │
  ├─ ④ 可选项: Sharpen (方向性 CFL=0.5 逆流锐化)
  ├─ ⑤ 可选项: FilterOrphan (清除体积<filterth 的碎片)
  ├─ ⑥ Recolor (连通分量标记) → 传播或并查集算法
  ├─ ⑦ BcClear (近边界快照到 0 或 1)
  │
  ▼
FinishStep: fcu_.time_curr = fcu_.iter_curr, PostStep 交换ghost数据
```

**关键数据字段**:

| 字段 | 含义 |
|------|------|
| `fcu_` (StepData) | 体积分数：time_prev, time_curr, iter_curr, iter_prev |
| `fcuu_` | Weymouth 散度项用的锐化体积分数 |
| `fccl_` | 颜色（连通分量标签） |
| `fcn_` | 法向向量（L1 归一化） |
| `fca_` | 平面常数 α |
| `fci_` | 界面掩码（bool） |
| `fcim_` | 图像向量（周期边界跨越计数） |
| `mebc_` | 对流传质边界条件 |

### 3.4 多相 VOF 求解器：`Vofm<EB>`

**文件**: `src/solver/vofm.h`, `src/solver/vofm.ipp`, `src/solver/vofm.cpp`

扩展单相 VOF 以处理多个不可混溶流体层：
- 使用 `Multi<FieldCell<Scal>>`（每层一个场的 `std::vector`）
- **法向平均**：同单元多层时按体积加权平均（`avgnorm0`/`avgnorm1` 参数）
- **Sweep 颜色追踪**：每层独立计算通量，颜色跨空单元传播
- **Recolor**：并查集或传播算法跨所有层
- 状态保存/加载：通过 HDF5

### 3.5 VOF 工具类：`UVof<M>`

**文件**: `src/util/vof.h`, `src/util/vof.ipp`, `src/util/vof.cpp`

| 功能 | 说明 |
|------|------|
| `DumpPoly` | 将 PLIC 界面 polygon 输出到 VTK |
| `DumpPolyMarch` | 用 marching squares/cubes 从体积分数提取等值面，节点插值 PLIC 法向计算平滑法向 |
| `Recolor` / `RecolorDirect` / `RecolorUnionFind` | 连通分量标记（迭代传播 / 并查集路径压缩） |
| `GetAdvectionBc` | 将对流传质 BC 转化为体积分数、颜色、图像、法向、α 的 BC |
| 模板工具 | `GetStencilValues`, `InterpolateToNodes`, `SumToNodesNan`, `GradientNodes` |

---

## 4. 支持文件

| 文件 | 说明 |
|------|------|
| `src/solver/advection.h` | 基类 `AdvectionSolver<M>`、`generic::Plic<Vect>` 结构体、`BCondAdvection` |
| `src/parse/vof.h` | 参数解析：`VofPar<M>` 结构体 + `ParseVofPar` |
| `src/solver/multi.h` | `Multi<T>` 模板——多相场容器 |
| `src/solver/solver.h` | `UnsteadySolver`、`UnsteadyIterativeSolver`、`StepData<T>` |
| `src/solver/trackerm.h` | `Trackerm<M>`——图像向量追踪（周期边界） |
| `src/solver/curv.h` | 曲率估算 `curvature::Particles` |
| `src/solver/partstrmeshm.h` | 粒子-弦方法曲率计算 `PartStrMeshM::Part()` |

---

## 5. 测试文件

| 文件 | 说明 |
|------|------|
| `src/test/advection/main.cpp` | 对流测试框架（使用 `Vof<M>` 或 `Vofm<M>`） |
| `src/test/advection/cmp.h` | 测试验证对比工具 |
| `src/test/reconst/plane.cpp` | Level-set vs PLIC 重构测试 |
| `src/test/reconst/main.cpp` | 重构测试入口 |

---

## 6. 总体算法架构图

```
┌──────────────────────────────────────────────────────┐
│              时间步 (Vof::MakeIteration)              │
├──────────────────────────────────────────────────────┤
│ 1. Setup: fcu_iter = fcu_prev + dt * source * (1-u)  │
│    fcuu = sharp(0.5) for Weymouth                     │
├──────────────────────────────────────────────────────┤
│ 2. Reconstruction (ReconstPlanes):                   │
│    a. 检测界面单元 (fci: 0<u<1)                       │
│    b. 计算法向 (Youngs + 高度函数)                     │
│    c. 计算平面常数 α = GetLineA(n,u,h)                │
├──────────────────────────────────────────────────────┤
│ 3. 方向分裂 Advection Sweeps (方向顺序每步旋转):      │
│    For each direction:                               │
│      a. 对 sweep 方向上的每个面:                       │
│         - 界面单元: GetLineFlux(n,α,h,flux,dt,d)      │
│         - 纯单元: flux = mixture_flux * u             │
│         - 传播颜色                                    │
│      b. 有限体积更新:                                  │
│         plain:       u += -dl                         │
│         aulisa:      EI + LE 交替                     │
│         weymouth:    u += fcuu*ds - dl                │
│      c. Clip [0,1], 通信/边界反射                      │
├──────────────────────────────────────────────────────┤
│ 4. 可选项: Sharpen, FilterOrphan                      │
│ 5. Recolor: 连通分量标记 (传播 或 并查集)              │
│ 6. BcClear: 近边界快照到 0 或 1                       │
└──────────────────────────────────────────────────────┘
```

## 7. 配置参数 (`VofPar<M>`)

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `scheme` | `weymouth` | 对流方案: plain / aulisa / weymouth |
| `dim` | 3 | 有效维度 |
| `clipth` | 1e-10 | 体积分数裁剪阈值 |
| `filterth` | 0 | 碎片过滤阈值 |
| `layers` | 4 | 多相层数 (仅 Vofm) |
| `sharpen` | false | 界面锐化开关 |
| `sharpen_cfl` | 0.5 | 锐化 CFL 数 |
| `recolor` | true | 连通分量标记 |
| `recolor_unionfind` | true | 使用并查集 |
| `recolor_reduce` | true | 缩减颜色空间 |
| `coalth` | 1e10 | 合并阈值 |
| `extrapolate_boundaries` | false | 外推边界 |
| `avgnorm0`/`avgnorm1` | 1.0 | 法向平均阈值 (仅 Vofm) |

---

## 8. Rising Bubble 测试算例 (Basilisk 基准)

基于 [Basilisk rising bubble test](http://basilisk.fr/src/test/rising.c)。模拟单个气泡在重力 + 表面张力下的浮力上升过程，用于验证终端速度和气泡形状。

### 8.1 文件结构

```
examples/208_rising/
├── std.conf              # 主配置文件
├── Makefile              # 构建/运行规则 (调用 ap.sim_base.makefile)
├── plot.py               # 可视化脚本 (Hausdorff 距离对比 + 界面线绘制)
├── README.md             # "based on basilisk.fr/src/test/rising.c"
├── ref/
│   └── c1g3l4s.txt       # MooNMD 参考界面线 (基准对比)
└── run_all               # 批量运行脚本 (完整精度 vs 简化粘度)

src/test/rising/
├── test                  # Python 测试脚本 (CTest 集成)
├── add.conf              # 测试用配置覆盖 (减少输出)
├── CMakeLists.txt         # CTest 注册
├── ref/
│   └── sm_0003.vtk       # 参考 VTK 输出 (Hausdorff 校验用)
└── .gitignore
```

### 8.2 物理参数 (std.conf)

| 参数 | 值 | 说明 |
|------|----|------|
| `dim` / `spacedim` | 2 | 2D 模拟 |
| `extent` | 2 | 域尺寸 |
| `backend` | native | 分布式后端 |
| `init_vf` | `list: sphere 0.5 0 0 0.25` | 初始圆形气泡，中心(0.5, 0)，半径 0.25 |
| `vel_init` | zero | 静止初始场 |
| `rho1` / `rho2` | 1000 / 100 | 液体 1000，气泡 100 (密度比 10:1) |
| `mu1` / `mu2` | 10 / 1 | 粘度比 10:1 |
| `sigma` | 24.5 | 表面张力系数 |
| `gravity` | (-0.98, 0, 0) | 向左的重力 (气泡向左浮升) |
| `tmax` | 3 | 模拟终点时间 |
| `cfl` / `cfla` | 0.8 | CFL 条件 |
| `fluid_solver` | proj (默认) | 投影法求解器 |
| `sharpen` | 1 | 启用界面锐化 |
| `vfsmooth` | 1 | VOF 平滑 |
| `tol` | 1e-4 | 压力/速度残差容限 |

### 8.3 边界条件 (std.conf:29-40)

```
wall 0 0 0 { box 0 0 0 10 }    # X0 面: 壁面
symm { box 0 0 0 10 }          # Y0 面: 对称 (默认)
symm { box 0 0.5 0 10 1e-6 }   # X1 面: 对称
```

- X 方向: 左端 wall (固壁)，右端 symm (对称)
- Y 方向: 两端 symm (对称)
- Z 方向 (2D): `hypre_periodic_z 1` (周期性)

### 8.4 网格配置

**Makefile** (`examples/208_rising/Makefile`):
```
m  = 256 64 1    # 全局网格 256×64
bs = 64 64 1     # 分块大小 64×64
np = 4            # 4 个 MPI 进程
```
→ 全局 256×64 cells，分 4 个 block，每 block 64×64。

**测试简化版** (`src/test/rising/test`):
```
m = 128 32 1
bs = 128 32 1
np = 1
```
→ 单进程 128×32 cells，用于快速验证。

### 8.5 运行方式

#### 方式一：Linux 完整运行 (生产级)

```bash
cd examples/208_rising

# ① 生成配置文件 (base.conf + mesh.conf + a.conf)
make cleanall
make conf

# ② 运行模拟
ap.mpirun ap.mfer --verbose

# ③ 或一步完成
make cleanrun
```

`make` 调用链：`Makefile` → `include $(shell ap.makesim)` → 读取 `ap.sim_base.makefile`
- `make conf`: 生成 `base.conf` (`ap.create_base_conf`)、`mesh.conf` (`ap.part`)、`a.conf` (`ap.create_a_conf`)
- `make run`: 执行 `ap.run ap.mfer --verbose`

关键辅助脚本 (`deploy/scripts/`):
| 脚本 | 功能 |
|------|------|
| `ap.makesim` | 输出 `ap.sim_base.makefile` 路径 |
| `ap.create_base_conf` | 生成 `base.conf` (复制 `sim_base.conf`) |
| `ap.part` | 根据 m, bs, np 生成 `mesh.conf` |
| `ap.run` | 启动模拟 (`ap.mfer`) |
| `ap.mpirun` | 带 MPI 启动 |

#### 方式二：CTest 测试 (快速验证)

测试在 `src/test/rising/` 下运行，通过 CMake/CTest 集成：

```bash
cd build
ctest -R rising
```

内部执行 `make -f sim.makefile cleanrun m='128 32 1' bs='128 32 1' np=1`，生成 `sm_0003.vtk`，与 `ref/sm_0003.vtk` 通过 Hausdorff 距离比较，阈值 0.001。

#### 方式三：直接编译运行 (Linux, 无 deploy 脚本)

```bash
cd src
make                                      # 编译 aphros
./ap.mfer std.conf mesh.conf              # 直接传配置文件
# 或
./ap.mfer --extra "include std.conf" --extra "include mesh.conf"
```

#### 方式四：Windows (MSVC + CMake, 当前环境)

编译 `ap.mfer` 完整可执行文件后运行：

```powershell
cd C:\Code\aphros\build\test\advection
.\ap.mfer.exe --verbose std.conf
```

需要预先准备:
1. `base.conf` ← 复制 `deploy/scripts/sim_base.conf`
2. `mesh.conf` ← 手动创建或运行 `ap.part` 等价逻辑

### 8.6 输出文件

| 文件模式 | 内容 | 格式 |
|----------|------|------|
| `sm_*.vtk` | Marching squares 提取的平滑界面线 | VTK PolyData |
| `s_*.vtk` | PLIC polygon 碎片 (dumppoly=1 时) | VTK PolyData |
| `vx_*.dat` / `vy_*.dat` / `p_*.dat` | 速度/压力/体积分数场 | raw binary |
| `*.h5` / `*.xmf` | HDF5 场输出 (dumpformat=hdf 时) | HDF5 |
| `stat.dat` | 统计摘要 (时间, 体积, CFL 等) | text |

输出帧命名: `<prefix>_<frame>.vtk`，`frame = t / dump_field_dt`。例如 `sm_0003.vtk` 对应 t=0.3 (dump_field_dt=0.1 时)。

### 8.7 可视化方法

#### 方法一: `plot.py` (Python + matplotlib)

`examples/208_rising/plot.py` 读取 `sm_*.vtk` 文件绘制界面轮廓线：

```bash
# 单帧: 叠加 aphros 结果 + MooNMD 参考
python plot.py 'sm_0030.vtk;Aphros' '../ref/c1g3l4s.txt;MooNMD;k;--' --output result.pdf

# 多文件对比
python plot.py --output comparison.png \
  'explviscous0.vtk;$\nabla \mathbf{u}$' \
  'explviscous1.vtk;$\nabla \mathbf{u} + \nabla \mathbf{u}^T$' \
  'ref/c1g3l4s.txt;ref (MooNMD);k;--'

# Hausdorff 距离计算
python plot.py 'sm_0030.vtk;a' 'ref/c1g3l4s.txt;ref' --dist 0 1
```

依赖: `aphros.plot` 模块 (位于 `deploy/scripts/aphros/plot.py`)，`matplotlib`, `numpy`。

**核心逻辑** (`plot.py`):
1. `read_lines_vtk(path)`: 读取 `sm_*.vtk` 中的 PLIC 界面 polygon 线
2. `read_lines_moonmd(path)`: 读取 MooNMD 参考 `.txt` 格式
3. `directed_hausdorff(a, b)`: 计算定向 Hausdorff 距离 (用于定量对比)
4. `ax.plot(*lines)`: matplotlib 绘制界面轮廓

#### 方法二: ParaView (VTK 3D 查看器)

`sm_*.vtk` 和 `s_*.vtk` 为标准 VTK PolyData 格式，可直接在 ParaView 中打开：
1. `File → Open` → 选择 `sm_0030.vtk`
2. 选择 "Surface" 或 "Wireframe" 表示
3. 可叠加 `vx_*.raw` / `vy_*.raw` 等场文件查看流场

#### 方法三: 场数据可视化 (raw/plain 格式)

如果设置 `dumpformat plain`，会输出 `.dat` 格式的体积分数、法向、平面常数文件：

```python
import aphros.plot
import matplotlib.pyplot as plt
import numpy as np

u = aphros.plot.ReadArray("u_0001.dat")    # 体积分数
a = aphros.plot.ReadArray("a_0001.dat")    # 平面常数 α
nx = aphros.plot.ReadArray("nx_0001.dat")  # 法向 x 分量
ny = aphros.plot.ReadArray("ny_0001.dat")  # 法向 y 分量

sx, sy = u.shape
hx, hy = 1.0/sx, 1.0/sy
x = (0.5 + np.arange(sx)) * hx
y = (0.5 + np.arange(sy)) * hy
xv, yv = np.meshgrid(x, y)
xv = xv.T; yv = yv.T

fig, ax = aphros.plot.PlotInitSq()
aphros.plot.PlotFieldGray(ax, u, vmin=0., vmax=1.)   # 灰度场
lines = aphros.plot.GetLines(xv, yv, a, nx, ny, hx, hy, u)  # PLIC 线段
aphros.plot.PlotLines(ax, *lines)
fig.savefig("result.png", dpi=150)
```

### 8.8 批量运行脚本 (`run_all`)

`examples/208_rising/run_all` 运行两个变体并对比：

```bash
# 变体 1: 默认粘度格式 (粘性应力 ∝ ∇u)
run 'explviscous1.vtk' ''

# 变体 2: 简化粘度 (粘性应力 ∝ ∇u + ∇u^T, 通过 --extra 注入)
run 'explviscous0.vtk' 'set int explviscous 0'

# 生成对比图
python plot.py --output 'explviscous.png' \
  'explviscous0.vtk;nop' \
  'explviscous1.vtk;yes' \
  'ref/c1g3l4s.txt;ref;k;--'
```

### 8.9 配置继承链

rising bubble 配置加载顺序：

```
sim_base.conf (deploy/scripts/，定义所有默认参数)
    ↓ include
a.conf (自动生成，include std.conf)
    ↓ include  
std.conf (examples/208_rising/，本例物理+数值参数)
    ↓ include (自动)
mesh.conf (ap.part 生成，网格拓扑)
    ↓ include (自动)
par.conf (可选，ap.pargen 生成的无量纲参数)
```

`--extra` 命令行参数可以在最后覆盖任何 `set` 变量。
