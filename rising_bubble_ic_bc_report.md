# Aphros Rising Bubble: 初始条件与边界条件设定总结

## 1. 配置入口

所有参数通过 `examples/208_rising/std.conf` 设定。解析链条：

```
distrsolver.cpp:86  Parser::ParseFile(std.conf)
  → parser.cpp:452  ParseFile() → ReadMultiline() → Cmd() → CmdSet()
    → vars.cpp:155  Vars::Map<string>::SetStr()  存入类型化参数表
```

`Vars` (`src/parse/vars.h:79-82`) 内含四个 Map：`String`、`Int`、`Double`、`Vect`，所有 `set` 命令最终落入对应 map。

---

## 2. 初始体积分数 (init_vf)

### 配置

```
set string init_vf list
set string list_path "inline
sphere 0.5 0 0 0.25
"
```

### 代码流程

```
kernel/hydro.ipp:1301  InitVf(fcvf, var, m)
  → init_u.h:131       分支 "list"
    → init.ipp:646      ReadPrimList() → 检测到 "inline"
      → init.ipp:22     UPrimList::GetPrimitives()
        → primlist.ipp:606  ParseSphere() 解析 "0.5 0 0 0.25"
          → primlist.ipp:149  构建 Primitive 结构体:
              cx=0.5, cy=0, cz=0, rx=ry=rz=0.25
              ls = (1 - ((x-c)/r)^2) * r^2   (level-set 函数)
```

### 体到每个 Cell

```
init.ipp:73  GetSphereOverlap(center, h, c, r)  对每个 cell 计算重叠体积
  → fcvf[c] = 0~1 的体积分数
```

### Primitive 支持类型 (`src/func/primlist.ipp:565-624`)

sphere, ring, box, roundbox, cylinder, polygon, polygon2, ruled, smooth_step

---

## 3. 初始速度 (vel_init)

### 配置

```
set string vel_init zero
```

### 代码流程

```
kernel/hydro.ipp:1354   fcvel.Reinit(m, Vect(0))  先全部置零
  → hydro.ipp:609       "zero" 分支: 直接 nop (已是零)
```

### 其他可选速度场 (Module 注册, `src/func/init_vel.cpp:104`)

| 名称 | 类 |
|------|----|
| `uniform` | `Uniform<M>` — 匀速场 |
| `taylor-green` | `TaylorGreen<M>` |
| `single_vortex` | `SingleVortex<M>` |
| `couette` | `Couette<M>` |
| `kelvin-helmholtz` | `KelvinHelmholtz<M>` |

---

## 4. 边界条件 (bc_path)

### 配置

```
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
```

语法：`<bc_type> [参数] { <primitive> }`。每个 `{}` 块定义一组面的集合（primitive box 覆盖到哪个面就应用该 BC）。

### 代码流程

```
kernel/hydro.ipp:1390  InitBc(var, m, known_keys)
  → hydro.ipp:769      var.String["bc_path"] → 读取 "inline" 剩余内容
    → codeblocks.cpp:110  ParseCodeBlocks() → 解析出 3 个 Block:

        Block 0: name="wall 0 0 0"   content="box 0 0 0 10"
        Block 1: name="symm"         content="box 0 0 0 10"
        Block 2: name="symm"         content="box 0 0.5 0 10 1e-6"

    → init_bc.h:191    ParseGroups():
        对每个面 face:
          对每个 Block (按顺序):
            lsmax(face_center) > 0 && is_boundary(face)
              → me_group[face] = block_index
              → me_nci[face] = neighbor_cell_id

    → init_bc.h:42/95   ParseBCondFluid():
        "wall 0 0 0" → BCondFluid{type=wall, velocity=(0,0,0)}
        "symm"       → BCondFluid{type=symm}

    → hydro.ipp:710    default_adv():
        wall  → halo=fill
        symm  → halo=reflect
```

### BC 覆盖规则

Block 按顺序处理，**后处理的覆盖前面的**。最终效果：

| 面 | 结果 |
|----|------|
| X0 (x=0) | Block0 wall → Block1 symm 覆盖 → `symm` |
| X1 (x=extent) | `symm` (Block1) |
| Y0 (y=0) | `symm` (Block1) |
| Y1 (y=1) | `symm` (Block1) |

Block2 的 box 中心在 y=0.5，薄层 1e-6 基本不覆盖有效面。

### BCondFluid 类型 (`src/solver/fluid.h:101-110`)

```cpp
enum class BCondFluidType {
  wall,          // 无滑移壁面
  slipwall,      // 滑移壁面
  inlet,         // 速度入口
  inletflux,     // 通量入口
  inletpressure, // 压力入口
  outlet,        // 出口
  outletpressure,// 压力出口
  symm,          // 对称
};
```

### BC 转换为求解器格式 (`src/solver/fluid.h:123-152`)

| 流体 BC | 速度 BC 类型 |
|---------|-------------|
| wall, inlet | Dirichlet (固定值) |
| slipwall, symm | Mixed (部分 Dirichlet + Neumann) |

---

## 5. 物理参数

### 配置

```
set double rho1 1000    # 液体密度
set double mu1 10       # 液体粘度
set double rho2 100     # 气泡密度
set double mu2 1        # 气泡粘度
set double sigma 24.5   # 表面张力系数
set vect gravity -0.98 0 0   # 向左的重力
```

### 代码流程 (`kernel/hydro.ipp:1810-1861`)

```
CalcMixture():
  fc_rho[c] = rho1*(1-a) + rho2*a   混合密度 (a=smoothed VF)
  fc_mu[c]  = mu1*(1-a) + mu2*a     混合粘度
  ff_rho[f] = rho1*(1-a) + rho2*a   面密度

  febp[f] += gravity · n · ff_rho[f]   重力体力加至面压力梯度

CalcSurfaceTension():
  fc_force ← 表面张力 (由 sigma 和曲率计算)
```

### 时间步约束

| 约束 | 来源 | 公式 |
|------|------|------|
| CFL(st) | 表面张力 | `dt < sqrt((rho1+rho2) * h^3 / (4π * sigma))` |
| CFL(visc) | 粘性 | `dt < min(mu1/rho1, mu2/rho2) * h^2` |
| CFL(conv) | 对流 | `dt < cfl * h / maxvel` |

---

## 6. 关键文件索引

| 功能 | 文件 | 行号 |
|------|------|------|
| 配置解析入口 | `src/distr/distrsolver.cpp` | 86 |
| Vars 参数存储 | `src/parse/vars.h` | 79-82 |
| Config 文件解析 | `src/parse/parser.cpp` | 452 |
| init_vf 调度 | `src/func/init_u.h` | 119-160 |
| init_vf list 实现 | `src/func/init.ipp` | 14-95 |
| Sphere primitive 解析 | `src/func/primlist.ipp` | 149-179, 606 |
| 体到体积分数 | `src/func/init.ipp` | 73 |
| init_vel 调度 | `src/util/hydro.ipp` | 30-614 |
| vel_init zero | `src/util/hydro.ipp` | 609 |
| BC 初始化 | `src/kernel/hydro.ipp` | 1390 |
| bc_path 解析与分组 | `src/util/hydro.ipp` | 725-779 |
| CodeBlocks 解析 | `src/parse/codeblocks.cpp` | 110-137 |
| Face 分组 | `src/func/init_bc.h` | 120-191 |
| BCondFluid 解析 | `src/func/init_bc.h` | 31-101 |
| BCondFluid 枚举 | `src/solver/fluid.h` | 101-110 |
| 物理参数混合 | `src/kernel/hydro.ipp` | 1810-1861 |
| 表面张力 | `src/kernel/hydro.ipp` | 1858-1861 |
| Module 注册基类 | `src/util/module.h` | 21-84 |
