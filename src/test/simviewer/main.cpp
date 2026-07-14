// VOF advection tests with SimViewer binary output.
// Matching aphros staged-execution pattern from src/test/advection.
// Usage: t.simviewer [translation|vortex] [resolution]

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>

#include <distr/distrbasic.h>
#include <solver/vof.h>

#include "simviewer.h"

using M = MeshCartesian<double, 2>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

static Vect VortexVel(Vect x, Scal t, Scal period) {
  Scal c = std::cos(M_PI * t / period);
  Scal sx = std::sin(M_PI * x[0]), sy = std::sin(M_PI * x[1]);
  Scal cx = std::cos(M_PI * x[0]), cy = std::cos(M_PI * x[1]);
  return Vect(-2.0 * c * cy * sx * sx * sy, 2.0 * c * cx * sx * sy * sy);
}

void RunSim(M& m, Vars& var) {
  auto sem = m.GetSem("simviewer");

  struct {
    std::unique_ptr<Vof<M>> as;
    FieldCell<Scal> fcu, fccl, fc_src;
    FieldEmbed<Scal> fe_flux;
    MapEmbed<BCondAdvection<Scal>> mf_cond;
    std::string outdir, test_name;
    Scal tmax, frame_dt, vortex_period, next_dump;
    Vect uniform_vel;
    uint32_t frame_count;
    int hl;
  } * ctx(sem);
  auto& s = *ctx;

  // Init stages (matching advection test pattern)
  if (sem.Nested("init-field")) {
    s.test_name = var.String["simviewer_test"];
    s.tmax = var.Double["simviewer_tmax"];
    s.frame_dt = var.Double["simviewer_frame_dt"];
    s.vortex_period = var.Double["simviewer_vortex_period"];
    s.uniform_vel = Vect(var.Double["simviewer_vel_x"], var.Double["simviewer_vel_y"]);
    Vect circle_c(var.Double["circle_cx"], var.Double["circle_cy"]);
    Scal circle_r = var.Double["circle_r"];

    s.fcu.Reinit(m, 0);
    for (auto c : m.Cells()) {
      s.fcu[c] = (m.GetCenter(c).dist(circle_c) < circle_r) ? 1.0 : 0.0;
    }
    m.Comm(&s.fcu);
  }
  if (sem("init-create")) {
    s.fc_src.Reinit(m, 0);
    s.fe_flux.Reinit(m, 0);
    s.fccl.Reinit(m, 0);

    for (auto f : m.Faces()) {
      Vect x = m.GetCenter(f);
      Vect vel = (s.test_name == "vortex") ?
        VortexVel(x, 0.0, s.vortex_period) : s.uniform_vel;
      s.fe_flux[f] = vel.dot(m.GetSurface(f));
    }

    Scal dx = m.GetCellSize()[0];
    Scal maxvel = (s.test_name == "vortex") ? 2.0 : s.uniform_vel.norm();
    Scal dt = var.Double["cfl"] * dx / std::max(maxvel, 1e-10);
    typename Vof<M>::Par p;
    p.dim = 2;
    s.as.reset(new Vof<M>(m, m, s.fcu, s.fccl, s.mf_cond,
                          &s.fe_flux, &s.fc_src, 0., dt, p));

    s.outdir = "simviewer_" + s.test_name;
    s.hl = var.Int["hl"];
    simviewer::InitOutput(s.outdir, m, s.hl);
    s.next_dump = 0;
    s.frame_count = 0;
    if (m.IsRoot()) std::cout << "Test: " << s.test_name
      << " tmax=" << s.tmax << " dt=" << dt << std::endl;
  }

  // Main loop (matching advection test pattern)
  sem.LoopBegin();
  if (sem("empty")) {}
  if (sem("checkloop")) {
    if (s.as->GetTime() >= s.tmax) sem.LoopBreak();
  }
  if (sem("vel")) {
    if (s.test_name == "vortex") {
      Scal t = s.as->GetTime();
      for (auto f : m.Faces()) {
        Vect x = m.GetCenter(f);
        s.fe_flux[f] = VortexVel(x, t, s.vortex_period).dot(m.GetSurface(f));
      }
    }
  }
  if (sem.Nested("start"))  { s.as->StartStep(); }
  if (sem.Nested("iter"))   { s.as->MakeIteration(); }
  if (sem.Nested("finish")) { s.as->FinishStep(); }
  if (sem("dump")) {
    Scal t = s.as->GetTime();
    if (t >= s.next_dump - 1e-12) {
      bool init = (s.frame_count == 0);
      simviewer::ExportVof(s.outdir, s.frame_count, s.as->GetField(), m, init, s.hl);
      auto plic = s.as->GetPlic();
      simviewer::ExportInterface(s.outdir, s.frame_count,
                                 *plic.vfci[0], *plic.vfcn[0], *plic.vfca[0], m, s.hl);
      simviewer::ExportNormal(s.outdir, s.frame_count,
                              *plic.vfci[0], *plic.vfcn[0], m, s.hl);
      simviewer::ExportVelocity(s.outdir, s.frame_count, m,
                                s.fe_flux, s.as->GetTimeStep(), s.hl);
      if (init) simviewer::ExportGrid(s.outdir, m, s.hl);
      s.frame_count++;
      simviewer::UpdateFrameCount(s.outdir, s.frame_count);
      if (m.IsRoot()) std::cout << "Frame " << (s.frame_count-1)
        << " t=" << t << " dt=" << s.as->GetTimeStep() << std::endl;
      s.next_dump += s.frame_dt;
    }
  }
  sem.LoopEnd();

  if (sem()) {
    if (m.IsRoot()) std::cout << "Done. " << s.frame_count
      << " frames in " << s.outdir << std::endl;
  }
}

int main(int argc, const char** argv) {
  std::string test = "translation"; int res = 64;
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "vortex") == 0) test = "vortex";
    else if (std::strcmp(argv[i], "translation") == 0) test = "translation";
    else { int r = std::atoi(argv[i]); if (r > 0) res = r; }
  }
  Scal tmax = (test == "vortex") ? 8.0 : 0.5;
  Scal frame_dt = (test == "vortex") ? 0.2 : 0.025;
  Scal cx = (test == "vortex") ? 0.5 : 0.25;
  Scal cy = (test == "vortex") ? 0.75 : 0.25;

  std::string conf = R"EOF(
set int dim 2
set string backend local
set int openmp 0
set int bx 1
set int by 1
set int bz 1
set int bsx )EOF" + std::to_string(res) + R"EOF(
set int bsy )EOF" + std::to_string(res) + R"EOF(
set int bsz 1
set int px 1
set int py 1
set int pz 1
set double cfl 0.5
set string simviewer_test )EOF" + test + R"EOF(
set double simviewer_tmax )EOF" + std::to_string(tmax) + R"EOF(
set double simviewer_frame_dt )EOF" + std::to_string(frame_dt) + R"EOF(
set double simviewer_vortex_period 8.0
set double simviewer_vel_x 1.0
set double simviewer_vel_y 0.4
set double circle_cx )EOF" + std::to_string(cx) + R"EOF(
set double circle_cy )EOF" + std::to_string(cy) + R"EOF(
set double circle_r 0.15
set int hl 2
set int loc_periodic_x 0
set int loc_periodic_y 0
)EOF";

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, RunSim, conf);
}
