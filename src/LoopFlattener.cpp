#include "libcdgbs/LoopFlattener.hpp"
#include <iostream>
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

  const Vec3 globalZ(0, 0, 1);
  for (size_t curveIdx = 0; curveIdx < curves.size(); ++curveIdx) {
    const auto& subCurves = curves[curveIdx];
    const auto& subNormals = normals[curveIdx];
    if (subCurves.size() != subNormals.size()) continue;

    // First pass: compute edge lengths and 2D directions by rotating each edge to align its normal with Z
    std::vector<Vec2> directions;
    std::vector<double> lengths;
    for (size_t subIdx = 0; subIdx < subCurves.size(); ++subIdx) {
      const auto& pts = subCurves[subIdx];
      const auto& norms = subNormals[subIdx];
      if (pts.size() < 2 || norms.size() != pts.size()) continue;
      for (size_t i = 0; i + 1 < pts.size(); ++i) {
        Vec3 edge = pts[i + 1] - pts[i];
        Vec3 normal = norms[i].normalized();
        double L3 = edge.norm();

        // rotation that maps 'normal' to Z
        Vec3 axis = normal.cross(globalZ);
        double axisLen = axis.norm();
        Vec3 rotated;
        if (axisLen < 1e-8) {
          // normal already near ±Z
          rotated = (normal.dot(globalZ) > 0 ? edge : -edge);
        }
        else {
          axis.normalize();
          double ang = std::acos(std::clamp(normal.dot(globalZ), -1.0, 1.0));
          Eigen::AngleAxisd R(ang, axis);
          rotated = R * edge;
        }

        // projected 2D direction in XY-plane
        Vec2 d2(rotated.x(), rotated.y());
        double L2 = d2.norm();
        Vec2 dir2d = (L2 > 1e-8 ? d2 / L2 : Vec2(1, 0));

        directions.push_back(dir2d);
        lengths.push_back(L3);
      }
    }

    size_t n = directions.size();
    if (n == 0) continue;

    // Compute signed exterior (turn) angles between successive segments
    std::vector<double> extAngles(n);
    for (size_t i = 0; i < n; ++i) {
      const Vec2& d1 = directions[i];
      const Vec2& d2 = directions[(i + 1) % n];
      double crossv = d1.x() * d2.y() - d1.y() * d2.x();
      double dotv = std::clamp(d1.dot(d2), -1.0, 1.0);
      extAngles[i] = std::atan2(crossv, dotv);
    }

    // Compute interior angles = π − exterior
    std::vector<double> intAngles(n);
    for (size_t i = 0; i < n; ++i) {
      intAngles[i] = M_PI - extAngles[i];
    }

    // Mark actual corners by detecting large exterior turns
    std::vector<bool> isCorner(n, false);
    const double thetaCorner = M_PI / 4; // >45° considered a true corner
    for (size_t i = 0; i < n; ++i) {
      if (std::abs(extAngles[i]) > thetaCorner)
        isCorner[i] = true;
    }

    // Normalize only non-corner interior angles so sum = (n-2)π
    double sumRes = 0, sumUnres = 0;
    for (size_t i = 0; i < n; ++i) {
      (isCorner[i] ? sumRes : sumUnres) += intAngles[i];
    }
    double target = (static_cast<double>(n) - 2.0) * M_PI;
    if (sumUnres > 1e-8) {
      double scale = (target - sumRes) / sumUnres;
      for (size_t i = 0; i < n; ++i)
        if (!isCorner[i]) intAngles[i] *= scale;
    }

    // Convert back to exterior angles
    for (size_t i = 0; i < n; ++i) {
      extAngles[i] = M_PI - intAngles[i];
    }

    // Second pass: build flattened loops in 2D
    Curve3D developedCurve;
    developedCurve.reserve(subCurves.size());
    Vec2 pos(0, 0);
    Vec2 dir = directions[0];
    size_t edgeIdx = 0;
    std::vector<Vec2> flatPts;
    flatPts.reserve(n + 1);
    flatPts.push_back(pos);

    for (size_t subIdx = 0; subIdx < subCurves.size(); ++subIdx) {
      const auto& pts = subCurves[subIdx];
      size_t m = pts.size();
      if (m < 2) continue;
      SubCurve3D devSub;
      devSub.reserve(m);
      devSub.emplace_back(pos.x(), pos.y(), 0);
      for (size_t i = 0; i + 1 < m; ++i, ++edgeIdx) {
        if (edgeIdx > 0) dir = Eigen::Rotation2D<double>(extAngles[edgeIdx - 1]) * dir;
        pos += dir * lengths[edgeIdx];
        devSub.emplace_back(pos.x(), pos.y(), 0);
        flatPts.push_back(pos);
      }
      developedCurve.push_back(std::move(devSub));
    }

    // Displacement correction to ensure closure
    double totalArc = std::accumulate(lengths.begin(), lengths.end(), 0.0);
    Vec2 delta = flatPts.front() - flatPts.back();
    double acc = 0; size_t idx = 0;
    for (auto& sub : developedCurve) {
      for (size_t i = 0;i < sub.size();++i, ++idx) {
        double w = (totalArc > 0 ? acc / totalArc : 0);
        Vec2 corr = delta * w;
        sub[i].x() += corr.x();
        sub[i].y() += corr.y();
        if (i + 1 < sub.size()) acc += lengths[idx];
      }
    }

    developedCurves.push_back(std::move(developedCurve));
  }

  std::cout << "Developing curves... DONE!" << std::endl;
  return developedCurves;
}