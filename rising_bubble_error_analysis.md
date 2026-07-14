# Aphros Rising Bubble：误差分析与参考解来源

## 1. 问题

将 Aphros 模拟的 t=3 时刻气泡界面与 MooNMD 参考解对比，发现严重不吻合。
经排查，根本原因是 **网格配置错误导致域尺寸不符**。

## 2. 根本原因：域尺寸不匹配

### Hysing 基准测试的标准域

Basilisk `rising.c`（http://basilisk.fr/src/test/rising.c）明确写道：

> "The domain will span [0:2]×[0:0.5] and will be resolved with 256×64 grid points."

- 域：[0, 2] × [0, 0.5]（矩形，4:1 宽高比）
- 网格：256 × 64
- 气泡中心 (0.5, 0)，半径 0.25，y=0 为对称面（仅模拟上半域）
- 重力 -0.98 在 x 方向（x 为垂直方向）

### BUILD.md 原始错误配置

| 参数 | 正确值（Hysing/Basilisk/Makefile） | BUILD.md 原始错误值 |
|------|--------------------------------------|---------------------|
| by | 1 | 2 |
| bsy | 64 | 128 |
| 网格 | 256×64 | 256×256 |
| 域 | 2×0.5 | 2×2 |

`examples/208_rising/Makefile` 中 `m = 256 64 1` 确认正确网格为 256×64。

### 错误影响

域从 2×0.5 变为 2×2 后，壁面约束大幅减弱，气泡上升过快且变形不同：
- Aphros（错误域 2×2）：气泡质心从 0.5 升至 1.44（上升 0.94）
- MooNMD（正确域 2×0.5）：气泡质心从 0.5 升至 1.07（上升 0.57）

### 修复

`mesh.conf` 改为：

```
set int bx 4
set int by 1
set int bsx 64
set int bsy 64
set double extent 2
```

`extent=2` 设 x 方向域长 2，y 方向按比例缩放：`64/256 × 2 = 0.5`，得到正确的 [0,2]×[0,0.5] 域。

## 3. 误差计算方法

### 使用的脚本

`build/error_analysis.py`，基于 `examples/208_rising/plot.py` 的 VTK 读取函数。

### 误差指标

对两组界面点集 A（Aphros）和 B（MooNMD）计算：

| 指标 | 公式 | 含义 |
|------|------|------|
| 有向 Hausdorff 距离 A→B | `max_{a∈A} min_{b∈B} ‖a-b‖` | Aphros 上每个点到 MooNMD 的最大最近距离 |
| 对称 Hausdorff 距离 | `max(d(A→B), d(B→A))` | 两个方向取最大，整体最大偏差 |
| 均值最小距离 | `mean_{a∈A} min_{b∈B} ‖a-b‖` | 平均偏差 |
| RMS 最小距离 | `sqrt(mean(min‖a-b‖²))` | 均方根偏差 |
| 质心距离 | `‖centroid_A - centroid_B‖` | 气泡整体位置差 |

### 坐标对齐

MooNMD 数据 `c1g3l4s.txt` 的第一列为 y，第二列为 x（`plot.py:read_lines_moonmd` 中 `yy, xx = lines`）。
读取后执行 `yy -= 0.5` 将 y 坐标平移到以对称面为中心。
Aphros VTK 数据用 `mirror=True` 镜像 y（因仅模拟上半域）。

### 命令

```powershell
# 仅对比图
python plot.py 'sm_0030.vtk;Aphros' '..\examples\208_rising\ref\c1g3l4s.txt;MooNMD;k;--' --output comparison.pdf --figsize 4 4

# 对比图 + Hausdorff 距离（打印到 stdout）
python plot.py 'sm_0030.vtk;Aphros' '..\examples\208_rising\ref\c1g3l4s.txt;MooNMD;k;--' --output comparison.pdf --figsize 4 4 --dist 0 1

# 完整误差分析
python error_analysis.py
```

需设置 `$env:PYTHONPATH = "$pwd\..\deploy\scripts"`（aphros Python 模块所在）。

### 错误配置下的误差值（域 2×2 vs 2×0.5，仅供参考）

| 指标 | 值 |
|------|-----|
| 对称 Hausdorff 距离 | 0.3955 |
| 均值最小距离 (Aphros→MooNMD) | 0.2271 |
| RMS 最小距离 (Aphros→MooNMD) | 0.2522 |
| 质心距离 | 0.3734 |

修复域尺寸后应大幅改善。

## 4. 参考解来源

### 数据文件

`examples/208_rising/ref/c1g3l4s.txt` — MooNMD 求解器在 t=3 时刻的气泡界面轮廓（650 个点）。

### 来源链

1. **Hysing et al. (2009)** 定义了二维气泡上升基准测试（Case 1 和 Case 2）：
   > S. Hysing, S. Turek, D. Kuzmin, N. Parsin, G. Schafbach,
   > *"Quantitative benchmark computations of two-dimensional bubble dynamics"*,
   > Int. J. Numer. Methods Fluids, 60(11), 1259–1288, 2009.

2. **MooNMD**（Mathematics and object-oriented Numerics in Magdeburg）求解器
   — 德国马格德堡大学 Turek 课题组的有限元流体求解器，用于生成参考解。

3. 数据托管在 **FEATFLOW** 基准测试网站：
   http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble/bubble_verification.html

4. Basilisk `rising.c` 直接引用同一份数据文件（`../c1g3l4s.txt`），
   gnuplot 脚本 `u 2:($1-0.5)` 对应 `plot.py` 的 `yy -= 0.5` 坐标平移。

### 文件名编码

`c1g3l4s`：c1 = Case 1，g3 = grid level 3，l4s = 网格细化参数。

## 5. 开源情况

| 代码 | 许可证 | 来源 |
|------|--------|------|
| Aphros | MIT | github.com/cselab/aphros（ETH Zurich，当前仓库） |
| Basilisk | GPL | basilisk.fr |
| MooNMD | 公开可下载 | featflow.de（FEATFLOW 项目，Turek 课题组） |

## 6. 相关文件索引

| 文件 | 说明 |
|------|------|
| `BUILD.md` §3 | 运行指南（已修复 mesh.conf） |
| `examples/208_rising/std.conf` | 物理参数与边界条件配置 |
| `examples/208_rising/Makefile` | 正确网格参数 `m = 256 64 1` |
| `examples/208_rising/plot.py` | 对比与 Hausdorff 距离计算脚本 |
| `examples/208_rising/ref/c1g3l4s.txt` | MooNMD 参考解（650 点界面轮廓） |
| `build/error_analysis.py` | 完整误差分析脚本 |
| `build/mesh.conf` | 网格配置（已修复为 256×64） |
| `rising_bubble_ic_bc_report.md` | 初始条件与边界条件代码流程分析 |
