#pragma once

#include <Eigen/Eigen>

namespace libcdgbs{
class LoopFlattener {
public:
  using Vec3 = Eigen::Vector3d;
  using Vec2 = Eigen::Vector2d;
  using SubCurve3D = std::vector<Vec3>;
  using Curve3D = std::vector<SubCurve3D>;
  using Curves3D = std::vector<Curve3D>;
  // Main method: takes input curves (each composed of sub-curves) and associated normals
  // Returns flattened curves into the XY plane (z = 0)
  static std::vector<Curve3D> developCurves(const std::vector<Curve3D>& curves, const std::vector<Curve3D>& normals);
};
}