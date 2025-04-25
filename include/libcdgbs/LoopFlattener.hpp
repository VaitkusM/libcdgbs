#pragma once

#include <Eigen/Core>

namespace libcdgbs {
  class LoopFlattener {
  public:
    using Vec3 = Eigen::Vector3d;
    using Vec2 = Eigen::Vector2d;
    using SubCurve3D = std::vector<Vec3>;
    using Curve3D = std::vector<SubCurve3D>;
    using Curves3D = std::vector<Curve3D>;
    using DoubleVector = std::vector<double>;
    // Main method: takes input curves (each composed of sub-curves) and associated normals
    // Returns flattened curves into the XY plane (z = 0)
    static Curve3D developLoop(const Curve3D& curve_loop, const Curve3D& normals, bool flipped = false);
    static Curve3D normalizeLoopAngles(const Curve3D& curve_loop, const Curve3D& normals, bool flipped = false);
    static Curve3D closeLoop(const Curve3D& curve_loop);
    static void getGeodesicCurvatures(
      const SubCurve3D& curve_pts,
      const SubCurve3D& normals,
      DoubleVector& curvatures,
      bool                loop = false,
      bool                integrated = false,
      bool                flipped = false
    );
    static double getGeodesicCurvature(
      const SubCurve3D& curve_pts,
      const SubCurve3D& normals,
      size_t idx,
      DoubleVector& curvatures,
      bool                loop = false,
      bool                integrated = false,
      bool                flipped = false
    );
    static double integrateCurvature(
      const DoubleVector& lengths,
      const DoubleVector& curvatures,
      SubCurve3D& curve_2d,
      double starting_angle = 0.0,
      Vec3 starting_pos = Vec3(0.0, 0.0, 0.0)
    );

  };
}