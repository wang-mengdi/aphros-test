// SimViewer binary export utilities for aphros
// Writes the SimViewer data format for 2D VOF visualization.
// Based on SimViewer data-format.md:
//   vof: fracmap indexed triangles (TopoFixed=true)
//   grid: static lines (frame 0 only)
//   interface: animated lines (PLIC segments)
//
// Uses only C++14-compatible features (no std::filesystem).

#pragma once

#include <array>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#define mkdir_p(path) _mkdir(path)
#else
#include <sys/stat.h>
#define mkdir_p(path) mkdir(path, 0755)
#endif

#include "geom/mesh.h"
#include "solver/reconst.h"

namespace simviewer {

// --- Simple path helpers (no std::filesystem needed) ---
inline std::string path_join(const std::string& a, const std::string& b) {
  if (a.empty()) return b;
  if (b.empty()) return a;
  if (a.back() == '/' || a.back() == '\\') return a + b;
  return a + "/" + b;
}

inline void ensure_dir(const std::string& path) {
  // Create directory, ignoring "already exists" errors
  mkdir_p(path.c_str());
  // On Windows, mkdir returns -1 for existing dir, which is fine
}

inline void make_dirs(const std::string& path) {
  // Create all parent directories
  std::string cur;
  for (size_t i = 0; i < path.size(); ++i) {
    if (path[i] == '/' || path[i] == '\\') {
      if (!cur.empty()) ensure_dir(cur);
    }
    cur += path[i];
  }
  if (!cur.empty()) ensure_dir(cur);
}

// --- Binary write helpers (little-endian) ---
inline void WriteU32(std::ostream& out, uint32_t val) {
  out.write(reinterpret_cast<const char*>(&val), sizeof(val));
}
inline void WriteF32(std::ostream& out, float val) {
  out.write(reinterpret_cast<const char*>(&val), sizeof(val));
}
inline void WriteVec2(std::ostream& out, float x, float y) {
  WriteF32(out, x);
  WriteF32(out, y);
}

// --- Directory setup ---
inline void InitOutput(const std::string& dirname) {
  // Remove existing output if present
  std::string cmd = "rmdir /s /q \"" + dirname + "\" 2>nul";
#ifdef _WIN32
  std::system(cmd.c_str());
#else
  cmd = "rm -rf \"" + dirname + "\"";
  std::system(cmd.c_str());
#endif

  make_dirs(path_join(path_join(dirname, "results"), "0"));

  // Write description.yaml
  {
    std::ofstream f(path_join(dirname, "description.yaml"));
    f << "Dimension: 2\n";
    f << "Radius: 0.55\n";
    f << "Objects:\n";
    // vof: fracmap indexed triangles, topo fixed
    f << "  - Name: vof\n";
    f << "    Animated: true\n";
    f << "    Primitive: Triangles\n";
    f << "    Shader: fracmap\n";
    f << "    Indexed: true\n";
    f << "    TopoFixed: true\n";
    f << "    Material:\n";
    f << "      Albedo: [0, 0, 1, 1]\n";
    // grid: static lines
    f << "  - Name: grid\n";
    f << "    Primitive: Lines\n";
    f << "    Material:\n";
    f << "      Albedo: [0.2, 0.2, 0.2, 1]\n";
    // interface: animated lines
    f << "  - Name: interface\n";
    f << "    Animated: true\n";
    f << "    Primitive: Lines\n";
    f << "    Material:\n";
    f << "      Albedo: [1, 0, 0, 1]\n";
  }

  // Write initial frame_count.txt
  {
    std::ofstream f(path_join(dirname, "frame_count.txt"));
    f << "1\n";
  }
}

inline void UpdateFrameCount(const std::string& dirname, uint32_t n) {
  std::ofstream f(path_join(dirname, "frame_count.txt"));
  f << n << "\n";
}

inline std::string frame_dir(const std::string& dirname, uint32_t frame) {
  return path_join(path_join(dirname, "results"), std::to_string(frame));
}

// --- VOF field export (fracmap indexed triangles) ---
template <class M>
void ExportVof(
    const std::string& dirname, uint32_t frame,
    const FieldCell<typename M::Scal>& vof, const M& m, bool is_initial) {
  using Scal = typename M::Scal;
  auto fdir = frame_dir(dirname, frame);
  make_dirs(fdir);

  std::ofstream f(path_join(fdir, "vof.out"), std::ios::binary);

  uint32_t num_cells = static_cast<uint32_t>(m.GetInBlockCells().GetSize().prod());
  uint32_t vertex_count = num_cells * 4;

  // vertex_count
  WriteU32(f, vertex_count);

  // positions: vertex_count * vec2<float32>
  for (auto c : m.Cells()) {
    for (size_t q = 0; q < 4; ++q) {
      auto pos = m.GetNode(m.GetNode(c, q));
      WriteVec2(f, static_cast<float>(pos[0]), static_cast<float>(pos[1]));
    }
  }

  // heats: vertex_count * float32 (same vof per cell -> 4 copies)
  for (auto c : m.Cells()) {
    float val = static_cast<float>(vof[c]);
    for (int i = 0; i < 4; ++i) WriteF32(f, val);
  }

  // indices (only frame 0, TopoFixed=true)
  if (is_initial) {
    uint32_t index_count = num_cells * 6;
    WriteU32(f, index_count);
    for (uint32_t i = 0; i < num_cells; ++i) {
      uint32_t indices[6] = {
          i * 4 + 0, i * 4 + 1, i * 4 + 2,
          i * 4 + 2, i * 4 + 1, i * 4 + 3};
      f.write(reinterpret_cast<const char*>(indices), sizeof(indices));
    }
  }
}

// --- Grid lines export (static, frame 0 only) ---
template <class M>
void ExportGrid(const std::string& dirname, const M& m) {
  auto fdir = frame_dir(dirname, 0);
  make_dirs(fdir);
  std::ofstream f(path_join(fdir, "grid.out"), std::ios::binary);

  auto res = m.GetInBlockCells().GetSize();
  uint32_t nx = static_cast<uint32_t>(res[0]);
  uint32_t ny = static_cast<uint32_t>(res[1]);

  auto origin = m.GetNode(IdxNode(0));
  auto top_right = origin + m.GetCellSize() * generic::Vect<double, 2>(nx, ny);
  float ox = static_cast<float>(origin[0]);
  float oy = static_cast<float>(origin[1]);
  float tx = static_cast<float>(top_right[0]);
  float ty = static_cast<float>(top_right[1]);
  auto dx = static_cast<float>(m.GetCellSize()[0]);
  auto dy = static_cast<float>(m.GetCellSize()[1]);

  uint32_t cnt = 2 * (nx + 1) + 2 * (ny + 1);
  WriteU32(f, cnt);

  for (uint32_t i = 0; i <= nx; ++i) {
    float x = ox + i * dx;
    WriteVec2(f, x, oy);
    WriteVec2(f, x, ty);
  }
  for (uint32_t j = 0; j <= ny; ++j) {
    float y = oy + j * dy;
    WriteVec2(f, ox, y);
    WriteVec2(f, tx, y);
  }
}

// --- Interface (PLIC) line segments export ---
template <class M>
void ExportInterface(
    const std::string& dirname, uint32_t frame,
    const FieldCell<bool>& interface_mask,
    const FieldCell<typename M::Vect>& normal,
    const FieldCell<typename M::Scal>& alpha,
    const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  auto fdir = frame_dir(dirname, frame);
  make_dirs(fdir);

  std::vector<float> segments;
  for (auto c : m.Cells()) {
    if (!interface_mask[c]) continue;
    auto poly = Reconst<Scal>::GetCutPoly(
        m.GetCenter(c), normal[c], alpha[c], m.GetCellSize());
    if (poly.size() >= 2) {
      for (size_t i = 0; i + 1 < poly.size(); ++i) {
        segments.push_back(static_cast<float>(poly[i][0]));
        segments.push_back(static_cast<float>(poly[i][1]));
        segments.push_back(static_cast<float>(poly[i + 1][0]));
        segments.push_back(static_cast<float>(poly[i + 1][1]));
      }
      if (poly.size() > 2) {
        segments.push_back(static_cast<float>(poly.back()[0]));
        segments.push_back(static_cast<float>(poly.back()[1]));
        segments.push_back(static_cast<float>(poly[0][0]));
        segments.push_back(static_cast<float>(poly[0][1]));
      }
    }
  }

  std::ofstream f(path_join(fdir, "interface.out"), std::ios::binary);
  uint32_t vertex_count = static_cast<uint32_t>(segments.size() / 2);
  WriteU32(f, vertex_count);
  f.write(reinterpret_cast<const char*>(segments.data()),
          segments.size() * sizeof(float));
}

} // namespace simviewer
