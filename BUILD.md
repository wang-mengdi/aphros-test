# Windows MSVC 编译指南

## 前置条件

- **Visual Studio 2022 Community**（或更高），安装 "使用 C++ 的桌面开发" 工作负载
  - 含 MSVC 编译器、CMake 3.31+、Windows 10 SDK
- **Git**
- **Python 3** + numpy + matplotlib（可视化用，非必须）

> ⚠️ 不需要：WSL、MinGW、GNU Make、MPI、HYPRE、HDF5、OpenCL。仅需 Visual Studio。

---

## 1. 克隆并准备

```powershell
git clone <你的fork地址> aphros
cd aphros
git checkout windows-msvc-support
```

### 1.1 手动创建 gitgen.cpp

`src/gitrev` 是 bash 脚本，Windows 无法执行。需手动创建占位文件：

```powershell
@"
// unknown
const char* kGitRev = "unknown";
const char* kGitMsg = "";
const char* kGitDiff = "";
"@ | Out-File -FilePath "src\util\gitgen.cpp" -Encoding ascii
```

---

## 2. 编译

### 2.1 初始化 MSVC 环境

每次新开 PowerShell 都要跑这句，否则 `cl.exe` / `nmake` 找不到：

```powershell
. "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
```

### 2.2 CMake 配置

```powershell
$env:APHROS_PREFIX = "$pwd\build\install"

mkdir build -Force
cd build

cmake ..\src `
  -G "NMake Makefiles" `
  -DCMAKE_BUILD_TYPE=Release `
  -DAPHROS_PREFIX="$env:APHROS_PREFIX" `
  -DUSE_MPI=OFF -DUSE_HDF=OFF -DUSE_HYPRE=OFF -DUSE_OPENCL=OFF `
  -DUSE_DIM3=OFF -DUSE_DIM1=OFF -DUSE_DIM4=OFF -DUSE_DIM2=ON `
  -DUSE_AVX=OFF -DUSE_OPENMP=OFF `
  -DUSE_BACKEND_CUBISM=OFF -DUSE_BACKEND_NATIVE=OFF -DUSE_BACKEND_LOCAL=ON `
  -DUSE_TESTS=ON -DUSE_MFER=ON `
  -DUSE_EXPLORER=OFF -DUSE_CONF2PY=OFF `
  -DUSE_WARNINGS=OFF
```

> **关键点**：
> - **必须用 `-G "NMake Makefiles"`**，默认 Visual Studio generator 会失败
> - **`APHROS_PREFIX` 是强制环境变量**，不设直接 fatal error
> - **`USE_WARNINGS=OFF`** 避免 GCC 风格的 `-Wextra`（MSVC 不认识）

### 2.3 构建

```powershell
# ① 构建 aphros 静态库（约 60-80 个 .obj）
nmake aphros

# ② 构建 ap.mfer.exe
nmake mfer

# ③ 可选：构建对流测试
nmake t.advection
```

产物位置：
| 文件 | 路径 |
|------|------|
| 静态库 | `build\aphros.lib` |
| 主程序 | `build\ap.mfer.exe` |
| 测试程序 | `build\test\advection\t.advection.exe` |

---

## 3. 运行 Rising Bubble（Hysing Case 1）

### 3.1 创建配置文件

在 `build\` 下创建 4 个文件：

**① `base.conf`** — 直接复制：

```powershell
Copy-Item ..\deploy\scripts\sim_base.conf build\base.conf
```

**② `mesh.conf`** — 256×64 单进程网格（域 2×0.5，匹配 Hysing 基准）：

```
set int px 1
set int py 1
set int pz 1
set int bx 4
set int by 1
set int bz 1
set int bsx 64
set int bsy 64
set int bsz 1
set double extent 2
```

全局尺寸 = `px·bx·bsx × py·by·bsy` = 4×64 × 1×64 = 256×64。
`extent=2` 设 x 方向域长为 2，y 方向按网格比例自动缩放：`64/256 × 2 = 0.5`，
得到域 [0,2]×[0,0.5]，与 Hysing 基准和 `examples/208_rising/Makefile`（`m = 256 64 1`）一致。
气泡在 y=0 对称面，仅模拟上半域。

**③ `std.conf`** — 基于 `examples/208_rising/std.conf`，改用 `backend local`：

```
set string backend local
set int openmp 0
set double extent 2
set int dim 2
set int spacedim 2
set int hypre_periodic_z 1

set string vel_init zero
set string init_vf list
set string list_path "inline
sphere 0.5 0 0 0.25
"

set int vfsmooth 1
set int sharpen 1
set double tol 1e-4
set double cfl 0.8
set double cfla 0.8
set int linsolver_gen_maxnorm 1
set string linsolver_symm conjugate
set int hypre_symm_maxiter 300
set int linreport 0
set int stat_step_every 10
set int report_step_every 10

set string bc_path "inline
wall 0 0 0 {
  box 0 0 0 10
}
symm {
  box 0 0 0 10
}
symm {
  box 0 0.5 0 10 1e-6
}
"

set double tmax 3
set double dump_field_dt 0.1

set double rho1 1000
set double mu1 10
set double rho2 100
set double mu2 1
set double sigma 24.5
set vect gravity -0.98 0 0

set string dumplist omz p
set int dumppoly 0
set int verbose_stages 0
```

**④ `a.conf`** — 入口文件，**必须覆盖 linsolver_gen 和 linsolver_vort**：

```
include base.conf
include mesh.conf
include std.conf
set string linsolver_gen conjugate
set string linsolver_vort conjugate
set int openmp 0
```

> ⚠️ `sim_base.conf` 默认 `linsolver_gen hypre` 和 `linsolver_vort hypre`。HYPRE 未编译，不覆盖会导致入口程序直接闪退（exit code -1073741819 / ACCESS_VIOLATION）。

### 3.2 运行

```powershell
cd build
.\ap.mfer.exe --verbose a.conf
```

预期输出：
```
Loading config from 'a.conf'
OpenMP num_threads: 1
Hostname ...
global mesh=(256,256)
surface tension dt=0.0026105
=====
STEP=0 t=0.00000000 dt=0.00000100 ...
Dump n=0 t=0.00000100 to sm_0000.vtk
STEP=10 t=0.02349548 dt=0.00261050 ...
...
```

产物：
- `sm_0000.vtk` ~ `sm_0030.vtk` — 31 帧 PLIC 界面（ParaView 可直接打开）
- `stat.dat` — 时间步统计

---

## 4. 可视化

### ParaView

```powershell
paraview                          # 启动 ParaView
File → Open → sm_..vtk            # 选中所有帧，自动识别为时间序列
```

### Python matplotlib

```powershell
cd build
python plot_rising.py             # 生成 30 帧面板图 rising_bubble_panels.png
```

### 与 MooNMD 参考比较

```powershell
cd build
Copy-Item ..\examples\208_rising\plot.py . -Force
python plot.py 'sm_0030.vtk;Aphros' '..\examples\208_rising\ref\c1g3l4s.txt;MooNMD;k;--' --output comparison.pdf
```

---

## 5. 源文件修改说明

`windows-msvc-support` 分支共修改 4 个源文件：

### `src/client.cmake` — MSVC 兼容性

```cmake
# 添加 Windows 宏定义
if(MSVC)
  add_compile_definitions(_USE_MATH_DEFINES NOMINMAX)
endif()

# 禁用 GCC 风格警告（MSVC 不认识 -Wextra）
# add_compile_options(-Wall -pedantic -Wextra)   ← 已注释
```

### `src/CMakeLists.txt` — 静态链接 + 防止模块注册被删除

```cmake
# aphros: MSVC 下用静态库 (SHARED 需要 dllexport)
if(MSVC)
  add_library(${T} STATIC ${ObjLibsTarget})
  target_link_options(${T} PRIVATE "/OPT:NOREF")  # 防止删除模块注册数组
else()
  add_library(${T} SHARED ${ObjLibsTarget})
endif()

# aphros_c: 同理
if(MSVC)
  add_library(${T} STATIC ...)
else()
  add_library(${T} SHARED ...)
endif()
```

### `src/util/CMakeLists.txt` — 移除 bash 脚本依赖

```cmake
# 删除: add_dependencies(${T} gitrev)
# gitrev 是 bash 脚本，Windows 无 bash
```

### `deploy/scripts/aphros/vtk.py` — NumPy 2.0 兼容

```python
# np.float → np.float64
# np.int   → np.int64
```

---

## 6. 常见故障

| 症状 | 原因 | 解决 |
|------|------|------|
| `cl : D8021: invalid numeric argument '/Wextra'` | MSVC 收到 GCC 标志 | 确认 cmake 用了 `-DUSE_WARNINGS=OFF` |
| `CMake Error: APHROS_PREFIX not set` | 环境变量缺失 | `$env:APHROS_PREFIX = "$pwd\build\install"` |
| `LINK : LNK1181: cannot open 'aphros.lib'` | aphros 未先构建 | 先跑 `nmake aphros` 再 `nmake mfer` |
| `LINK : LNK1181: cannot open 'aphros_c.lib'` | aphros_c 未构建或仍是 SHARED | 确认 CMakeLists 中 `if(MSVC)` 分支用了 `STATIC` |
| `variable 'hl' of type 'int' not found` | 缺少 `base.conf` | `a.conf` 需要有 `include base.conf` |
| 运行闪退，exit code `-1073741819` | linsolver_gen/vort 仍是 hypre | `a.conf` 中加 `set string linsolver_gen conjugate` |
| `nmake` 不可用 | 未初始化 MSVC 环境 | 重新跑 `. vcvars64.bat` |
| `cmake` 不可用 | 同上 | VS 自带的 cmake 也在 vcvars 初始化后加入 PATH |
| `MSB8066: Custom build ... gitrev` | gitrev.bat 占位缺失 | 创建空文件 `src/gitrev.bat` |
