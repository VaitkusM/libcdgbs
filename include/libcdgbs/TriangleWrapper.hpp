#pragma once

// ─── Pull in the C API ─────────────────────────────────
extern "C" {
#include "triangle_api.h"   // context*, triangle_mesh_create, etc. :contentReference[oaicite:0]{index=0}
#include "triangle.h"       // triangleio typedef :contentReference[oaicite:1]{index=1}
}

// Undefine Triangle’s macros so they won’t clash with Eigen
#ifdef dest
#undef dest
#endif
#ifdef mark
#undef mark
#endif

#include <vector>
#include <Eigen/Core>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

namespace libcdgbs {

  /// Your TriMesh type
  using Mesh = OpenMesh::TriMesh_ArrayKernelT<>;

  /// Wraps the wo80/C‐API of Triangle into a C++ class
  class TriangleWrapper {
  public:
    TriangleWrapper(char const* options = "pQzYq");
    ~TriangleWrapper();

    /// loops: each loop is a sequence of sub‑curves (each a polyline)
    Mesh triangulate(
      std::vector<std::vector<std::vector<Eigen::Vector3d>>> const& loops,
      double L_target
    );

  private:
    context* ctx_;  // the native Triangle context
  };

} // namespace libcdgbs
