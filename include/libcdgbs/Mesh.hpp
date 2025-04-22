#pragma once
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

namespace libcdgbs {
  struct MeshTraits : public OpenMesh::DefaultTraits {
    using Point = OpenMesh::Vec3d; // the default would be Vec3f
    using Normal = OpenMesh::Vec3d;
    using TexCoord2D = OpenMesh::Vec2d;
  };
  using Mesh = OpenMesh::TriMesh_ArrayKernelT<MeshTraits>;
}