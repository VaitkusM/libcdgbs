#include "libcdgbs/LoopFlattener.hpp"

using namespace libcdgbs;

using Vec3 = Eigen::Vector3d;
using Vec2 = Eigen::Vector2d;
using SubCurve3D = std::vector<Vec3>;
using Curve3D = std::vector<SubCurve3D>;
using Curves3D = std::vector<Curve3D>;

// Thx to ChatGPT
// Main method: takes input curves (each composed of sub-curves) and associated normals
// Returns flattened curves into the XY plane (z = 0)
Curves3D LoopFlattener::developCurves(const Curves3D& curves, const Curves3D& normals) {
  Curves3D developedCurves;
  developedCurves.reserve(curves.size());

  for (size_t curveIdx = 0; curveIdx < curves.size(); ++curveIdx) {
    const auto& subCurves = curves[curveIdx];
    const auto& subNormals = normals[curveIdx];

    if (subCurves.size() != subNormals.size()) {
      continue;
    }

    std::vector<Vec2> directions;
    std::vector<double> lengths;

    // First pass: extract edge directions and lengths
    for (size_t subIdx = 0; subIdx < subCurves.size(); ++subIdx) {
      const auto& points = subCurves[subIdx];
      const auto& norms = subNormals[subIdx];
      if (points.size() < 2 || points.size() != norms.size()) continue;

      for (size_t i = 0; i < points.size() - 1; ++i) {
        Vec3 edge = points[i + 1] - points[i];
        Vec3 normal = norms[i].normalized();

        Vec3 tangent1 = edge.normalized();
        Vec3 tangent2 = normal.cross(tangent1).normalized();
        tangent1 = tangent2.cross(normal).normalized();

        double edgeLength = edge.norm();
        Eigen::Matrix<double, 3, 2> basis;
        basis.col(0) = tangent1;
        basis.col(1) = tangent2;

        Vec2 projected = basis.transpose() * edge;
        projected = projected.normalized() * edgeLength;

        directions.push_back(projected.normalized());
        lengths.push_back(edgeLength);
      }
    }

    // Compute angles between directions
    std::vector<double> angles;
    for (size_t i = 0; i < directions.size(); ++i) {
      const Vec2& dir1 = directions[i];
      const Vec2& dir2 = directions[(i + 1) % directions.size()];
      double dot = std::clamp(dir1.dot(dir2), -1.0, 1.0);
      double angle = std::acos(dot);
      double sign = (dir1.x() * dir2.y() - dir1.y() * dir2.x()) > 0 ? 1.0 : -1.0;
      angles.push_back(sign * angle);
    }

    // Normalize angles to sum up to (n - 2) * PI for closed polygon with n vertices
    double targetSum = (directions.size() - 2) * M_PI;
    double actualSum = std::accumulate(angles.begin(), angles.end(), 0.0);
    double correctionFactor = targetSum / actualSum;
    for (auto& angle : angles) angle *= correctionFactor;

    // Second pass: construct the developed curve
    Curve3D developedCurve;
    Vec2 currentPos(0.0, 0.0);
    Vec2 currentDir(1.0, 0.0);
    size_t edgeIdx = 0;

    std::vector<Vec2> developedFlatPoints;
    developedFlatPoints.emplace_back(currentPos);

    for (size_t subIdx = 0; subIdx < subCurves.size(); ++subIdx) {
      const auto& points = subCurves[subIdx];
      const auto& norms = subNormals[subIdx];
      if (points.size() < 2 || points.size() != norms.size()) continue;

      SubCurve3D developedSubCurve;
      developedSubCurve.reserve(points.size());
      developedSubCurve.emplace_back(currentPos.x(), currentPos.y(), 0.0);

      for (size_t i = 0; i < points.size() - 1; ++i, ++edgeIdx) {
        double angle = angles[edgeIdx];
        Eigen::Rotation2D<double> rot(angle);
        currentDir = rot * currentDir;
        currentPos += currentDir * lengths[edgeIdx];
        developedSubCurve.emplace_back(currentPos.x(), currentPos.y(), 0.0);
        developedFlatPoints.emplace_back(currentPos);
      }

      developedCurve.push_back(std::move(developedSubCurve));
    }

    // Compute total arclength and displacement to close the loop
    double totalArc = std::accumulate(lengths.begin(), lengths.end(), 0.0);
    Vec2 loopDisplacement = developedFlatPoints.front() - developedFlatPoints.back();

    // Apply correction to each point by relative arclength
    double accumulatedLength = 0.0;
    size_t flatPointIdx = 0;
    for (size_t subIdx = 0; subIdx < developedCurve.size(); ++subIdx) {
      auto& sub = developedCurve[subIdx];
      for (size_t i = 0; i < sub.size(); ++i, ++flatPointIdx) {
        double weight = accumulatedLength / totalArc;
        Vec2 correction = weight * loopDisplacement;
        sub[i].x() += correction.x();
        sub[i].y() += correction.y();
        if (i < sub.size() - 1 && edgeIdx < lengths.size()) {
          accumulatedLength += lengths[flatPointIdx];
        }
      }
    }

    developedCurves.push_back(std::move(developedCurve));
  }

  return developedCurves;
}