set_project("aphros-benchmarks")
set_version("0.1.0")
set_languages("c11", "cxx14")
add_rules("mode.debug", "mode.release")

-- Keep 2D available; opt in to 3D for the standard PLIC benchmarks.
-- Keep the existing CMake files as the upstream reference.
option("aphros_3d")
    set_default(false)
    set_showmenu(true)
    set_description("Build the 3D local backend in addition to 2D")
option_end()

add_includedirs("src")
add_defines(
    "_USE_MATH_DEFINES", "NOMINMAX",
    "_USE_DIM1_=0", "_USE_DIM2_=1", "_USE_DIM4_=0",
    "_USE_MPI_=0", "_USE_HDF_=0", "_USE_AVX_=0",
    "_USE_BACKEND_CUBISM_=0", "_USE_BACKEND_LOCAL_=1",
    "_USE_BACKEND_NATIVE_=0", "_USE_HYPRE_=0", "_USE_AMGX_=0",
    "_USE_OPENCL_=0"
)

add_defines("_USE_DIM3_=" .. (has_config("aphros_3d") and "1" or "0"))

-- Equivalent to the object modules collected by src/CMakeLists.txt for the
-- local backend in the enabled dimensions. Disabled optional backends are omitted.
local aphros_sources = {
    "src/color/color.c",
    "src/distr/report.cpp",
    "src/distr/distr.cpp",
    "src/distr/distr_particles.cpp",
    "src/distr/local.cpp",
    "src/distr/distrsolver.cpp",
    "src/distr/distrbasic.cpp",
    "src/dump/dumper.cpp",
    "src/dump/hdf.cpp",
    "src/dump/raw.cpp",
    "src/dump/xmf.cpp",
    "src/dump/vtk.cpp",
    "src/dump/dump.cpp",
    "src/func/primlist.cpp",
    "src/func/init.cpp",
    "src/func/init_vel.cpp",
    "src/func/init_contang.cpp",
    "src/geom/mesh.cpp",
    "src/inside/main.c",
    "src/linear/linear.cpp",
    "src/march/main.c",
    "src/overlap/overlap.cpp",
    "src/parse/parser.cpp",
    "src/parse/codeblocks.cpp",
    "src/parse/vars.cpp",
    "src/parse/evalexpr.cpp",
    "src/parse/template.cpp",
    "src/parse/argparse.cpp",
    "src/solver/solver.cpp",
    "src/solver/approx.cpp",
    "src/solver/embed.cpp",
    "src/solver/approx_eb.cpp",
    "src/solver/vofm.cpp",
    "src/solver/vof.cpp",
    "src/solver/tracer.cpp",
    "src/solver/particles.cpp",
    "src/solver/normal.cpp",
    "src/solver/curv.cpp",
    "src/solver/partstrmeshm.cpp",
    "src/solver/convdiffi.cpp",
    "src/solver/convdiffe.cpp",
    "src/solver/convdiffvg.cpp",
    "src/solver/simple.cpp",
    "src/solver/proj.cpp",
    "src/solver/proj_eb.cpp",
    "src/solver/fluid_dummy.cpp",
    "src/solver/electro.cpp",
    "src/util/suspender.cpp",
    "src/util/system.c",
    "src/util/sysinfo.cpp",
    "src/util/hydro.cpp",
    "src/util/linear.cpp",
    "src/util/fluid.cpp",
    "src/util/visual.cpp",
    "src/util/convdiff.cpp",
    "src/util/vof.cpp",
    "src/util/distr.cpp",
    "src/util/mpi.cpp",
    "src/util/events.cpp",
    "src/util/timer.cpp",
    "make/windows/gitgen.cpp",
    "src/util/filesystem.cpp",
    "src/util/git.cpp",
    "src/util/posthook_default.cpp",
    "src/util/subcomm_dummy.cpp",
    "src/util/histogram.cpp",
    "src/util/logger.cpp",
    "src/util/format.cpp",
    "src/util/fixed_allocator.cpp",
    "src/util/hydro_post.cpp",
    "src/young/main.c",
    "src/kernel/hydro.cpp"
}

target("t.simviewer")
    set_kind("binary")
    set_default(true)
    set_runtimes("MD")
    add_ldflags("/OPT:NOREF")
    for _, source in ipairs(aphros_sources) do
        add_files(source)
    end
    add_files("src/test/simviewer/main.cpp")

target("ap.mfer")
    set_kind("binary")
    set_default(false)
    set_runtimes("MD")
    add_ldflags("/OPT:NOREF")
    for _, source in ipairs(aphros_sources) do
        add_files(source)
    end
    add_files("src/aphros_c/git.cpp", "src/aphros_c/main.cpp",
              "src/aphros_c/parser.cpp", "src/main.c")

target("t.advection")
    set_kind("binary")
    set_default(false)
    set_runtimes("MD")
    add_ldflags("/OPT:NOREF")
    for _, source in ipairs(aphros_sources) do
        add_files(source)
    end
    add_files("src/test/advection/main.cpp")
