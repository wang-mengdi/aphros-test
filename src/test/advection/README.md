# Advection test

Time-reversed advection of circle with stream function `sin(x)sin(y)`.

## Build and run

    ./conf
    ./t.advection

## Plot

    ./plot.py

## Optional recorded 3D transport inputs (2026-09-25)

The existing velocity modes and initialization defaults are unchanged.
The test driver supports these additional inputs:

- `set string init_vel deformation3d`: on the unit cube,
  `u=2*sin(pi*x)^2*sin(2*pi*y)*sin(2*pi*z)`,
  `v=-sin(pi*y)^2*sin(2*pi*z)*sin(2*pi*x)`, and
  `w=-sin(pi*z)^2*sin(2*pi*x)*sin(2*pi*y)`.
  Requires `spacedim=3` and `dim=3`. All components reverse at `t >= revt`.
  Choose a timestep that divides `revt` exactly for symmetric reversal.
- `set string init_vf raw` plus `set string init_vf_raw_path initial.raw`:
  use the existing `dump::Raw` reader to load headerless native Float64
  volume fractions, x contiguous, matching the full mesh size. This driver
  restricts that option to one whole-domain mesh block and checks file length
  and finite values in [0,1]. On this Windows build the format is little endian.
- `set int advection_walls 1`: build native reflecting scalar boundary
  conditions and set normal velocity flux to zero on nonperiodic faces.
  Set both `loc_periodic_x/y/z` and `hypre_periodic_x/y/z` consistently to 0
  on wall axes. The option defaults to 0. This does not add momentum physics.

These options change the prescribed inputs and test boundary setup, not the
VOF reconstruction, plane-volume formulas, sweep order, clipping or curvature.
The v14 comparison uses the exact initial cell averages and saved outputs from
SimLiquid v12 (64^3, dt=1/256, revt=0.5, tmax=1). Its interface stays away from
walls; it is not a contact-angle or wall-interface accuracy test.
