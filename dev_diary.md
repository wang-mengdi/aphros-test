# Windows/MSVC 适配开发日记

> 日期: 2026-06-25 ~ 2026-07-01
> 目标: 在 Windows 上编译运行 Aphros 2D PLIC VOF
> 环境: Windows 11, Visual Studio 2022 Community (MSVC 14.44), CMake 3.31.6

---

## 第一阶段: 探索与编译

### 1. 理解代码结构

用户要求分析 Aphros 的 2D PLIC VOF 代码结构。通过阅读源码，写成了 `aphros_2d_plic_vof_report.md`，涵盖：
- 编译方案 (CMake + Make 双系统)
- `Reconst<Scal>` 几何运算层
- `UNormal<M>` 法向计算 (Youngs + 高度函数)
- `Vof<EB>` / `Vofm<EB>` 单相/多相求解器
- 时间步流程 (重构 → 方向分裂 → 锐化 → Recolor)

### 2. 寻找测试算例

用户要求找 drop falling 类测试。发现：
- **最接近**: `examples/205_multivof/buoyancy_breakup/` (浮力碎裂)
- **标准 benchmark**: `examples/208_rising/` (Basilisk rising bubble, Hysing Case 1)
- **基础测试**: `src/test/advection/` (时间反演剪切流)

### 3. 首次编译尝试

环境检测:
- CMake: ✅ (VS 2022 自带, `C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe`)
- MSVC: ✅ (14.44.35207)
- WSL: ❌ (有二进制但未安装 Linux 发行版)
- Make/GCC/Cygwin: ❌

CMake 配置命令 (仅2D, 无MPI/HDF5/HYPRE):
```
cmake src/ -DUSE_MPI=OFF -DUSE_HDF=OFF -DUSE_HYPRE=OFF -DUSE_OPENCL=OFF
  -DUSE_DIM3=OFF -DUSE_DIM1=OFF -DUSE_DIM4=OFF -DUSE_DIM2=ON
  -DUSE_AVX=OFF -DUSE_OPENMP=OFF -DUSE_BACKEND_CUBISM=OFF
  -DUSE_BACKEND_NATIVE=OFF -DUSE_BACKEND_LOCAL=ON
  -DUSE_MFER=ON
```

### 4. 编译问题及修复

#### 问题 1: gitgen.cpp 缺失
- **原因**: `src/gitrev` 是 bash 脚本, Windows 无法执行
- **修复**: 手动创建 `src/util/gitgen.cpp`:
  ```cpp
  const char* kGitRev = "unknown";
  const char* kGitMsg = "";
  const char* kGitDiff = "";
  ```

#### 问题 2: -Wextra 非法参数
- **原因**: `src/client.cmake:91` 的 `add_compile_options(-Wall -pedantic -Wextra)` 是 GCC 风格, MSVC 不支持
- **修复**: 注释掉该行

#### 问题 3: M_PI 未定义
- **原因**: MSVC 需要 `_USE_MATH_DEFINES` 宏才定义 `M_PI`
- **修复**: `src/client.cmake` 添加:
  ```cmake
  if(MSVC)
    add_compile_definitions(_USE_MATH_DEFINES NOMINMAX)
  endif()
  ```

#### 问题 4: aphros.lib 找不到 (导入库问题)
- **原因**: MSVC 构建 SHARED DLL 时, 无 `__declspec(dllexport)` 不生成导入库
- **第一次尝试**: 添加 `WINDOWS_EXPORT_ALL_SYMBOLS` → 生成了 `.lib` 但仍有符号缺失
- **最终方案**: MSVC 下改用 STATIC 库 (`src/CMakeLists.txt:362-366`)

#### 问题 5: kForceLink_init_vel 未解析
- **原因**: 静态库链接时, 模块注册数组 (`kReg_conjugate[]`) 被 `/OPT:REF` 删除。该数组的值从不被读取, 仅作为副作用调用 `RegisterModule()`
- **修复**: `src/CMakeLists.txt` 添加 `target_link_options(aphros PRIVATE "/OPT:NOREF")`

### 5. 编译结果

成功编译:
- `ap.mfer.exe` — 完整求解器 (Hydro)
- `t.advection.exe` — 对流测试

---

## 第二阶段: 运行 2D Advection 测试

### 运行

创建最小配置:
```
set int hl 2
set int loc_maxcomm 16
set int dim 2
...
set string advection_solver vof
set string vof_scheme aulisa
```

结果: 32×32 网格, tmax=0.3, 40 步, 体积守恒 sumu=0.12617, 总耗时 0.058s

### 可视化

- Python `plot.py` 读取 `u/a/nx/ny_*.dat` 文件
- 修复 `deploy/scripts/aphros/vtk.py` 的 numpy 2.0 兼容问题 (`np.float` → `np.float64`, `np.int` → `np.int64`)
- 生成 `u_0001.png` (灰度 VF 场 + PLIC 界面线)

---

## 第三阶段: Rising Bubble 完整模拟

### 7. 理解配置系统

Rising bubble 配置继承链:
```
a.conf → include base.conf (sim_base.conf 全默认)
       → include mesh.conf (ap.part 生成)
       → include std.conf (用户提供, 物理参数)
       → include add.conf (可选覆盖)
```

### 8. 配置文件准备

创建了 `std.conf` (基于 `examples/208_rising/std.conf`, 改 `backend local`):
- 物理: ρ₁=1000, ρ₂=100, σ=24.5, gravity=(-0.98,0,0)
- 网格: extent=2, bs×by = 8 blocks for OpenMP
- 初始: sphere at (0.5,0,0), r=0.25
- BC: wall + 2× symm

### 9. 第一次运行: 闪退 (ACCESS_VIOLATION)

- `enable_fluid 1` + `enable_surftens 1` → 崩溃
- `enable_fluid 0` → 成功 (advection-only, 只跑 VOF)
- 原因定位: `sim_base.conf` 默认 `linsolver_gen hypre` / `linsolver_vort hypre`
  HYPRE 未编译 → 查找模块失败 → crash

### 10. conjugate 求解器可以工作

- 设置 `linsolver_gen conjugate`, `linsolver_symm conjugate`, `linsolver_vort conjugate`
- 加上 `/OPT:NOREF` 修复后, conjugate 正常工作了 (Jacobi 也可以)
- 全物理模拟开始正常运行

### 11. 输出文件

256×256 网格, tmax=3.0, dt≈0.0026:
- `sm_0000.vtk` ~ `sm_0030.vtk`: 31 帧 PLIC 界面 (marching squares)
- `stat.dat`: 时间步统计
- 可视化脚本 `plot_rising.py`: 30 帧面板图

---

## 第四阶段: SimLiquid 对比分析

### 12. 代码审查

分析了 `C:\Code\SimLiquid\`:
- xmake 构建, 依赖 Eigen/amgcl/spdlog/TBB
- `SimBuilder::BuildBubbleRising` 设置相同的 Hysing Case 1 参数
- **发现 SimLiquid 缺少粘度参数** (`m_ViscosityLiquid`/`m_ViscosityGas`)

### 13. 算法对比

| 模块 | Aphros | SimLiquid |
|------|--------|-----------|
| 法向 | Youngs + 高度函数 | MYC |
| α 偏移 | 解析三次方程 | 简化分段公式 |
| 通量 | GetLineFlux 精确 | ClipWithHalfPlane |
| 对流 | plain/aulisa/weymouth | Weymouth 式 |
| 线性求解 | HYPRE/CG/Jacobi | amgcl |
| 表面张力 | 粒子-弦法 | 中心差分 CSF |
| 粘性 | 隐式扩散 | 未实现 |
| BC | wall/symm/inlet 完整 | 全壁面零速 |

### 14. 对比方案讨论

1. **先降级 Aphros**: 禁用 HF + 改用简单表面张力公式, 匹配 SimLiquid 算法
2. **定量指标**: 质心轨迹, 终端速度, 体积守恒, Hausdorff 距离
3. **网格收敛性**: 多分辨率对比验证一致性

---

## 修改文件汇总

### 直接修改的源文件

| 文件 | 改动 | 说明 |
|------|------|------|
| `src/client.cmake` | 注释 `-Wextra` | MSVC 不支持 GCC 风格警告 |
| | 添加 `_USE_MATH_DEFINES NOMINMAX` | MSVC 需要此宏定义 M_PI |
| `src/CMakeLists.txt` | MSVC下 `aphros` 改为 STATIC | 无需导出符号 |
| | 添加 `/OPT:NOREF` | 防止链接器删除模块注册数组 |
| | MSVC下 `aphros_c` 改为 STATIC | 同上 |
| `src/util/CMakeLists.txt` | 删除 `add_dependencies(${T} gitrev)` | gitrev 是 bash 脚本无法运行 |
| `deploy/scripts/aphros/vtk.py` | `np.float` → `np.float64` | NumPy 2.0 兼容 |
| | `np.int` → `np.int64` | 同上 |
| `src/gitrev.bat` | 新建空文件 | Windows 占位 |

### 新创建的文件

| 文件 | 说明 |
|------|------|
| `src/util/gitgen.cpp` | Git 版本信息 (手动生成) |
| `aphros_2d_plic_vof_report.md` | PLIC VOF 代码结构报告 |
| `rising_bubble_ic_bc_report.md` | 初始条件与边界条件分析 |
| `aphros_vs_simliquid_comparison.md` | 两库算法对比 |
| `build/plot_rising.py` | rising bubble 可视化脚本 |
| `build/vis.py` | advection 可视化脚本 |
| `C:\Code\aphros\build/` | CMake 构建产物目录 |

---

## 后续计划

1. **版本控制**: fork cselab/aphros, 分支 `windows-msvc-support`, 提交所有修改
2. **算法降级**: 禁用高度函数法向 + 替换表面张力为简单中心差分 (匹配 SimLiquid)
3. **定量对比**: Hausdorff 距离, 质心轨迹, 终端速度, 网格收敛性
4. **添加 CI**: Windows MSVC + vcpkg 依赖, 自动构建和测试
