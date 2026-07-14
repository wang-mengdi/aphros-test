// SimViewer binary export for 2D VOF — aphros edition
// Objects: vof (fracmap triangles), grid (static lines),
//          interface (PLIC segments), normal, velocity

#pragma once
#include <array>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#ifdef _WIN32
#include <direct.h>
#define sim_mkdir(p) _mkdir(p)
#else
#include <sys/stat.h>
#define sim_mkdir(p) mkdir(p, 0755)
#endif
#include "geom/mesh.h"
#include "solver/reconst.h"

namespace simviewer {

// --- helpers ---
inline std::string j(const std::string& a, const std::string& b) {
  if (a.empty()) return b; if (b.empty()) return a;
  if (a.back()=='/'||a.back()=='\\') return a+b; return a+"/"+b;
}
inline void mkdirs(const std::string& p) {
  std::string c; for (size_t i=0;i<p.size();++i){if(p[i]=='/'||p[i]=='\\'){if(!c.empty())sim_mkdir(c.c_str());}c+=p[i];}
  if(!c.empty())sim_mkdir(c.c_str());
}
inline std::string fdir(const std::string& root, uint32_t f) { return j(j(root,"results"),std::to_string(f)); }

inline void w32(std::ostream& o, uint32_t v){o.write((const char*)&v,4);}
inline void wf(std::ostream& o, float v)    {o.write((const char*)&v,4);}
inline void wv2(std::ostream& o, float x, float y){wf(o,x);wf(o,y);}

// --- halo helper: extract hl from block size difference ---
template <class M>
int GetHl(const M& m) {
  // block size = physical + 2*hl in each dimension
  // We don't have the physical size without config, but hl is stored
  // in the mesh's internal state. Simplest: just return 2 (matches aphros default).
  return 2;
}

// --- interior cell iteration (block already excludes halos) ---
template <class M>
void ForEachCell(const M& m, std::function<void(IdxCell)> body) {
  for (auto c : m.Cells()) body(c);
}

// --- Init ---
template <class M>
void InitOutput(const std::string& dirname, const M& m, int hl) {
  std::string cmd = "rmdir /s /q \"" + dirname + "\" 2>nul";
#ifdef _WIN32
  std::system(cmd.c_str());
#else
  std::system(("rm -rf \""+dirname+"\"").c_str());
#endif
  mkdirs(j(j(dirname,"results"),"0"));

  auto full = m.GetInBlockCells().GetSize();
  auto h = m.GetCellSize();
  // First physical node is at index hl (skip halo nodes before domain)
  int nx = full[0], ny = full[1];
  float ox = (float)(m.GetNode(IdxNode(0))[0] + hl * h[0]);
  float oy = (float)(m.GetNode(IdxNode(0))[1] + hl * h[1]);

  std::ofstream f(j(dirname,"description.yaml"));
  f << "Dimension: 2\nRadius: 0.55\n";
  f << "GridResolution: [" << nx << ", " << ny << "]\n";
  f << "GridSpacing: " << h[0] << "\n";
  f << "GridOrigin: [" << ox << ", " << oy << "]\n";
  f << "Objects:\n"
       "  - Name: vof\n    Animated: true\n    Primitive: Triangles\n"
       "    Shader: fracmap\n    Indexed: true\n    TopoFixed: true\n"
       "    Material:\n      Albedo: [0, 0, 1, 1]\n"
       "  - Name: grid\n    Primitive: Lines\n"
       "    Material:\n      Albedo: [0.2, 0.2, 0.2, 1]\n"
       "  - Name: interface\n    Animated: true\n    Primitive: Lines\n"
       "    Material:\n      Albedo: [1, 0, 0, 1]\n"
       "  - Name: normal\n    Animated: true\n    Primitive: Lines\n"
       "    Material:\n      Albedo: [0, 1, 0, 1]\n"
       "  - Name: velocity\n    Animated: true\n    Primitive: Lines\n"
       "    Material:\n      Albedo: [1, 1, 0, 1]\n";
  std::ofstream fc(j(dirname,"frame_count.txt")); fc << "1\n";
}
inline void UpdateFrameCount(const std::string& d, uint32_t n) {
  std::ofstream f(j(d,"frame_count.txt")); f << n << "\n";
}

// --- VOF ---
template <class M>
void ExportVof(const std::string& dir, uint32_t frame,
               const FieldCell<typename M::Scal>& vof, const M& m, bool init, int hl) {
  auto fd = fdir(dir,frame); mkdirs(fd);
  std::ofstream f(j(fd,"vof.out"), std::ios::binary);
  uint32_t nc = (uint32_t)m.GetInBlockCells().GetSize().prod();
  w32(f, nc*4);
  ForEachCell(m, [&](IdxCell c){
    for(size_t q=0;q<4;++q){auto p=m.GetNode(m.GetNode(c,q));wv2(f,(float)p[0],(float)p[1]);}
  });
  ForEachCell(m, [&](IdxCell c){
    float v=(float)vof[c]; for(int i=0;i<4;++i)wf(f,v);
  });
  if(init){w32(f,nc*6); for(uint32_t i=0;i<nc;++i){uint32_t x[6]={i*4+0,i*4+1,i*4+2,i*4+2,i*4+1,i*4+3};f.write((const char*)x,24);}}
}

// --- Grid ---
template <class M>
void ExportGrid(const std::string& dir, const M& m, int hl) {
  auto fd = fdir(dir,0); mkdirs(fd);
  std::ofstream f(j(fd,"grid.out"), std::ios::binary);
  auto full=m.GetInBlockCells().GetSize(); auto h=m.GetCellSize();
  uint32_t nx=(uint32_t)full[0], ny=(uint32_t)full[1];
  float ox=(float)(m.GetNode(IdxNode(0))[0]+hl*h[0]), oy=(float)(m.GetNode(IdxNode(0))[1]+hl*h[1]);
  float dx=(float)h[0], dy=(float)h[1];
  w32(f,2*(nx+1)+2*(ny+1));
  for(uint32_t i=0;i<=nx;++i){float x=ox+i*dx;wv2(f,x,oy);wv2(f,x,oy+ny*dy);}
  for(uint32_t j=0;j<=ny;++j){float y=oy+j*dy;wv2(f,ox,y);wv2(f,ox+nx*dx,y);}
}

// --- Interface ---
template <class M>
void ExportInterface(const std::string& dir, uint32_t frame,
                     const FieldCell<bool>& mask, const FieldCell<typename M::Vect>& n,
                     const FieldCell<typename M::Scal>& a, const M& m, int hl) {
  using S=typename M::Scal; auto fd=fdir(dir,frame); mkdirs(fd);
  std::vector<float> seg;
  ForEachCell(m, [&](IdxCell c){
    if(!mask[c])return; auto p=Reconst<S>::GetCutPoly(m.GetCenter(c),n[c],a[c],m.GetCellSize());
    if(p.size()<2)return;
    for(size_t i=0;i+1<p.size();++i){seg.push_back((float)p[i][0]);seg.push_back((float)p[i][1]);seg.push_back((float)p[i+1][0]);seg.push_back((float)p[i+1][1]);}
    if(p.size()>2){seg.push_back((float)p.back()[0]);seg.push_back((float)p.back()[1]);seg.push_back((float)p[0][0]);seg.push_back((float)p[0][1]);}
  });
  std::ofstream f(j(fd,"interface.out"),std::ios::binary); w32(f,(uint32_t)(seg.size()/2)); f.write((const char*)seg.data(),seg.size()*4);
}

// --- Normal ---
template <class M>
void ExportNormal(const std::string& dir, uint32_t frame,
                  const FieldCell<bool>& mask, const FieldCell<typename M::Vect>& n,
                  const M& m, int hl) {
  auto fd=fdir(dir,frame); mkdirs(fd);
  std::vector<float> seg; float len=(float)(m.GetCellSize()[0]*0.5);
  ForEachCell(m, [&](IdxCell c){
    if(!mask[c])return; auto cen=m.GetCenter(c); auto end=cen+n[c]*len;
    seg.push_back((float)cen[0]);seg.push_back((float)cen[1]);seg.push_back((float)end[0]);seg.push_back((float)end[1]);
  });
  std::ofstream f(j(fd,"normal.out"),std::ios::binary); w32(f,(uint32_t)(seg.size()/2)); f.write((const char*)seg.data(),seg.size()*4);
}

// --- Velocity ---
template <class M>
void ExportVelocity(const std::string& dir, uint32_t frame, const M& m,
                    const FieldEmbed<typename M::Scal>& fe, typename M::Scal dt, int hl) {
  using V=typename M::Vect; using S=typename M::Scal;
  auto fd=fdir(dir,frame); mkdirs(fd);
  S dx=m.GetCellSize()[0], dy=m.GetCellSize()[1];
  std::vector<float> seg;
  ForEachCell(m, [&](IdxCell c){
    V ctr=m.GetCenter(c);
    S u=(fe[m.GetFace(c,IdxNci(0))]+fe[m.GetFace(c,IdxNci(1))])*0.5f/dy;
    S v=(fe[m.GetFace(c,IdxNci(2))]+fe[m.GetFace(c,IdxNci(3))])*0.5f/dx;
    V end=ctr+V(u,v)*dt;
    seg.push_back((float)ctr[0]);seg.push_back((float)ctr[1]);seg.push_back((float)end[0]);seg.push_back((float)end[1]);
  });
  std::ofstream f(j(fd,"velocity.out"),std::ios::binary); w32(f,(uint32_t)(seg.size()/2)); f.write((const char*)seg.data(),seg.size()*4);
}

} // namespace simviewer
