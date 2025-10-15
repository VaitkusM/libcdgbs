#pragma once
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

namespace libcdgbs {
  struct MeshTraits : public OpenMesh::DefaultTraits {
    using Point = OpenMesh::Vec3d; // the default would be Vec3f
    using Normal = OpenMesh::Vec3d;
    using TexCoord2D = OpenMesh::Vec2d;
  };
  class Mesh : public OpenMesh::TriMesh_ArrayKernelT<MeshTraits> { 
  public:
    struct ClosestResult {
      FaceHandle face;
      Point   pt_xyz;
      Point   pt_bary;
      double  dist2    = std::numeric_limits<double>::infinity();
    };
    ClosestResult findClosestFace(Point pt) const;
    double dist2Face(Point pt, FaceHandle ff) const;

    Point interpolateInFace(FaceHandle f, Point bary) const;
    

    static void barycentricCoordinates(
      const Point& pt,
      const Point& p1,
      const Point& p2,
      const Point& p3,
      double& u, double& v, double& w
    );
  };
}