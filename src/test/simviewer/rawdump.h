// Minimal VOF-only wrapper that dumps raw float output matching SimLiquid's .raw format.
// Compile as part of the simviewer test to get raw VOF data for comparison.
//
// Usage: t.simviewer.exe rawof translation 64
// Output: simviewer_translation/vof_frame_XX.raw (resolution*resolution floats)

#pragma once
#include <fstream>
#include <string>
#include "geom/mesh.h"

namespace simviewer {

template <class M, class Scal>
void ExportVofRaw(const std::string& dirname, uint32_t frame,
                  const FieldCell<Scal>& vof, const M& m, int hl = 2) {
  auto h = m.GetCellSize();
  int nx = static_cast<int>(m.GetInBlockCells().GetSize()[0]);
  int ny = static_cast<int>(m.GetInBlockCells().GetSize()[1]);

  std::vector<float> buf(nx * ny, 0.0f);
  for (auto c : m.Cells()) {
    auto idx = m.GetIndexCells().GetMIdx(c);
    int i = static_cast<int>(idx[0]);
    int j = static_cast<int>(idx[1]);
    if (i >= 0 && i < nx && j >= 0 && j < ny) {
      buf[j * nx + i] = static_cast<float>(vof[c]);
    }
  }

  std::string fname = dirname + "/vof_frame_" + std::to_string(frame) + ".raw";
  std::ofstream f(fname, std::ios::binary);
  f.write(reinterpret_cast<const char *>(buf.data()), buf.size() * sizeof(float));
}

// Export normal (nx,ny) and alpha as raw floats, row-major j*nx+i
template <class M, class Scal>
void ExportNormalAlphaRaw(const std::string& dirname, uint32_t frame,
                          const FieldCell<bool>& fci,
                          const FieldCell<typename M::Vect>& fcn,
                          const FieldCell<Scal>& fca,
                          const M& m, int hl = 2) {
  int nx = static_cast<int>(m.GetInBlockCells().GetSize()[0]);
  int ny = static_cast<int>(m.GetInBlockCells().GetSize()[1]);
  // normal: 2 floats per cell
  std::vector<float> bufn(nx * ny * 2, 0.0f);
  // alpha: 1 float per cell
  std::vector<float> bufa(nx * ny, 0.0f);
  for (auto c : m.Cells()) {
    auto idx = m.GetIndexCells().GetMIdx(c);
    int i = static_cast<int>(idx[0]);
    int j = static_cast<int>(idx[1]);
    if (i >= 0 && i < nx && j >= 0 && j < ny) {
      int k = j * nx + i;
      if (fci[c]) {
        bufn[k * 2]     = static_cast<float>(fcn[c][0]);
        bufn[k * 2 + 1] = static_cast<float>(fcn[c][1]);
        bufa[k]         = static_cast<float>(fca[c]);
      }
    }
  }
  {
    std::string fn = dirname + "/normal_frame_" + std::to_string(frame) + ".raw";
    std::ofstream f(fn, std::ios::binary);
    f.write(reinterpret_cast<const char*>(bufn.data()), bufn.size() * sizeof(float));
  }
  {
    std::string fa = dirname + "/alpha_frame_" + std::to_string(frame) + ".raw";
    std::ofstream f(fa, std::ios::binary);
    f.write(reinterpret_cast<const char*>(bufa.data()), bufa.size() * sizeof(float));
  }
}

} // namespace simviewer
