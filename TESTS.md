# Aphros VOF Advection → SimViewer

用 aphros VOF 求解器运行对流测试，输出 SimViewer 可直接可视化的二进制数据。

## 快速开始

以下所有命令在 **PowerShell**（已执行 `vcvars64.bat`）中运行，当前目录为 `aphros/build`。

### 1. 编译

```powershell
# 0. 初始化 MSVC 环境（每次新开 PowerShell 都要执行）
. "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"

# 如果 vcvars 没把 nmake/cl 加入 PATH，手动补一下：
$env:PATH = "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.42.34433\bin\Hostx64\x64;$env:PATH"

# 1. cmake 配置（只需一次）
$env:APHROS_PREFIX = "$pwd\install"
cmake ..\src -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release `
  -DAPHROS_PREFIX="$env:APHROS_PREFIX" `
  -DUSE_MPI=OFF -DUSE_HDF=OFF -DUSE_HYPRE=OFF -DUSE_OPENCL=OFF `
  -DUSE_DIM3=OFF -DUSE_DIM1=OFF -DUSE_DIM4=OFF -DUSE_DIM2=ON `
  -DUSE_AVX=OFF -DUSE_OPENMP=OFF `
  -DUSE_BACKEND_CUBISM=OFF -DUSE_BACKEND_NATIVE=OFF -DUSE_BACKEND_LOCAL=ON `
  -DUSE_TESTS=ON -DUSE_MFER=ON -DUSE_WARNINGS=OFF

# 编译 simviewer 测试
nmake t.simviewer
```

产物：`build\test\simviewer\t.simviewer.exe`

### 2. 运行 Translation 测试

```powershell
.\test\simviewer\t.simviewer.exe translation 64
```

> 参数：`translation` 或 `vortex`，可选分辨率（默认 64）

运行后在当前目录生成 `simviewer_translation\` 文件夹。

### 3. 用 SimViewer 可视化

```powershell
# 启动 SimViewer（位于 C:\Code\SimLiquid）
cd C:\Code\SimLiquid
xmake run viewer
```

在 SimViewer 窗口中：**File → Open** → 选择 `C:\Code\aphros-test\build\simviewer_translation` 文件夹。

- 蓝色三角网格 = VOF 分数场
- 红色线段 = PLIC 重建界面
- 灰色网格 = 计算网格
- 空格键播放/暂停，← → 逐帧浏览

## 测试场景

### Translation（均匀平移）

圆形液滴在匀速场中平移，验证 VOF 对流精度：

```
t.simviewer.exe translation [分辨率]
```

| 参数 | 值 |
|------|-----|
| 初始形状 | 圆心 (0.25, 0.25)，半径 0.15 |
| 速度场 | u=1.0, v=0.4（恒定） |
| 模拟时间 | 0 → 0.5 |
| 输出帧间隔 | 0.025（~20 帧） |

### Single Vortex（单涡）

圆形液滴在周期性涡流场中变形后复原（对标 SimLiquid 单涡测试）：

```
t.simviewer.exe vortex [分辨率]
```

| 参数 | 值 |
|------|-----|
| 初始形状 | 圆心 (0.5, 0.75)，半径 0.15 |
| 速度场 | 单涡流场，周期 T=8.0 |
| 模拟时间 | 0 → 8.0 |
| 输出帧间隔 | 0.04（~200 帧，25fps） |

## 输出格式

```
simviewer_translation/
├── description.yaml      # SimViewer 场景描述
├── frame_count.txt        # 总帧数
└── results/
    ├── 0/
    │   ├── vof.out        # VOF fracmap 三角网格（含拓扑）
    │   ├── grid.out       # 网格线（仅第 0 帧）
    │   └── interface.out  # PLIC 界面线段
    ├── 1/
    │   ├── vof.out        # 无拓扑（TopoFixed）
    │   └── interface.out
    └── ...
```

数据格式遵循 [SimViewer Data Format](https://github.com/zhxxch/SimViewer/blob/main/docs/data-format.md)。

## 文件

```
src/test/simviewer/
├── CMakeLists.txt
├── main.cpp           # translation + vortex 两个测试场景
├── simviewer.h        # SimViewer 二进制导出（vof/grid/interface）
└── README.md
```
