#include "libcdgbs/LoopFlattener.hpp"
#include "libcdgbs/GeomUtil.hpp"
#include <iostream>
using namespace libcdgbs;
using namespace GeomUtil;

using Vec3 = Eigen::Vector3d;
using Vec2 = Eigen::Vector2d;
using SubCurve3D = std::vector<Vec3>;
using Curve3D = std::vector<SubCurve3D>;
using Curves3D = std::vector<Curve3D>;

void LoopFlattener::getGeodesicCurvatures(
  const SubCurve3D& curve_pts,
  const SubCurve3D& normals,
  DoubleVector& curvatures,
  bool                loop,
  bool                integrated,
  bool                flipped
)
{
  size_t num_corners = curve_pts.size() - (loop ? 1 : 2);
  curvatures.clear();
  curvatures.resize(num_corners);

  for (size_t ii = 0; ii < num_corners; ++ii) {
    curvatures[ii] = getGeodesicCurvature(
      curve_pts,
      normals,
      ii + 1,
      curvatures,
      loop,
      integrated,
      flipped
    );
  }
}

double LoopFlattener::getGeodesicCurvature(
  const SubCurve3D& curve_pts,
  const SubCurve3D& normals,
  size_t idx,
  DoubleVector& curvatures,
  bool                loop,
  bool                integrated,
  bool                flipped
)
{
  size_t num_pts = curve_pts.size();
  if (num_pts != normals.size()) {
    return 0.0;
  }

  size_t idx_m1 = circular_index(idx, -1, num_pts);
  size_t idx_p1 = circular_index(idx, +1, num_pts);

  Vector v1 = curve_pts[idx_p1] - curve_pts[idx];
  Vector v2 = curve_pts[idx_m1] - curve_pts[idx];
  Vector nn = normals[idx];
  nn *= flipped ? -1.0 : 1.0;
  Vector v1_proj = projectToPlane(v1, nn);
  Vector v2_proj = projectToPlane(v2, nn);


  double angle = getAngle(
    v1_proj,
    v2_proj,
    nn,
    true
  );
  return angle / (integrated ? 1.0 : getVoronoiLength(curve_pts, idx));
}

double LoopFlattener::integrateCurvature(
  const DoubleVector& lengths,
  const DoubleVector& curvatures,
  SubCurve3D& curve_2d,
  double starting_angle,
  Vec3 starting_pos
)
{
  curve_2d.clear();
  if ((lengths.size() - 1) != curvatures.size()) {
    return starting_angle;
  }
  curve_2d.resize(lengths.size() + 1);

  double curr_angle = starting_angle;
  Vector curr_dir = Vector(cos(starting_angle), sin(starting_angle), 0.0);
  curve_2d[0] = starting_pos;
  curve_2d[1] = curve_2d[0] + lengths[0] * curr_dir;
  for (size_t ii = 2; ii < curve_2d.size(); ++ii) {
    curr_angle += curvatures[ii - 2];
    curr_dir = Vector(cos(curr_angle), sin(curr_angle), 0.0);
    curve_2d[ii] = curve_2d[ii - 1] + lengths[ii - 1] * curr_dir;
  }
  return curr_angle;
}



Curve3D LoopFlattener::developLoop(const Curve3D& curves, const Curve3D& normals) {
  Curve3D developedCurves(curves.size());

  double curr_angle = 0.0;
  Vector curr_pos(0.0, 0.0, 0.0);
  for (size_t side = 0; side < curves.size(); ++side) {
    size_t side_p1 = circular_index(side, 1, curves.size());
    const auto& curv = curves[side];
    const auto& norm = normals[side];
    const auto& curv_p1 = curves[side_p1];
    size_t num_pts = curv.size();

    DoubleVector edge_lengths(num_pts - 1, 1.0 / (num_pts - 1.0));
    DoubleVector curvatures(num_pts - 2, 0.0);

    bool flipped = false;

    getEdgeLengths(
      curv,
      edge_lengths,
      false//ribbon.periodic
    );

    getGeodesicCurvatures(
      curv,
      norm,
      curvatures,
      false,//ribbon.periodic,
      true,
      flipped
    );

    curr_angle = integrateCurvature(
      edge_lengths,
      curvatures,
      developedCurves[side],
      curr_angle,
      curr_pos
    );

    curr_pos = developedCurves[side].back();
    auto corner_angle = GeomUtil::getAngle(
      curv_p1[1] - curv_p1[0],
      curv[curv.size() - 2] - curv[curv.size() - 1],
      normals[side_p1][0],
      true
    );

    curr_angle += corner_angle;
  }

  return developedCurves;
}

Curve3D LoopFlattener::normalizeLoopAngles(
  const Curve3D& curve_loop,
  const Curve3D& normals,
  bool flipped
)
{
  size_t num_sides = curve_loop.size();
  Curve3D normalized_curv(num_sides);

  // Collecting angles
  DoubleVector interior_angles, normalized_interior_angles;
  double angle_sum = 0.0;
  for (size_t side = 0; side < num_sides; ++side) {
    size_t side_p1 = circular_index(side, 1, num_sides);
    const auto& curv = curve_loop[side];
    const auto& norm = normals[side];
    const auto& curv_p1 = curve_loop[side_p1];
    size_t num_pts = curv.size();
    DoubleVector interior_angles_for_side(num_pts - 2, M_PI);

    DoubleVector geodesic_curvatures;
    getGeodesicCurvatures(
      curv,
      norm,
      geodesic_curvatures,
      false,//ribbon.periodic,
      true,
      flipped
    );
    for (size_t ii = 0; ii < geodesic_curvatures.size(); ++ii) {
      interior_angles_for_side[ii] = unsignAngle(M_PI - geodesic_curvatures[ii]);
    }
    angle_sum += std::accumulate(
      interior_angles_for_side.begin(),
      interior_angles_for_side.end(),
      0.0
    );

    interior_angles.insert(
      interior_angles.end(),
      interior_angles_for_side.begin(),
      interior_angles_for_side.end()
    );

    auto corner_angle = GeomUtil::getAngle(
      curv_p1[1] - curv_p1[0],
      curv[curv.size() - 2] - curv[curv.size() - 1],
      normals[side_p1][0],
      false,
      false
    );
    interior_angles.push_back(corner_angle);
    angle_sum += corner_angle;
  }

  // Normalizing angles
  normalized_interior_angles = interior_angles;
  auto scaling_factor = M_PI * double(interior_angles.size() - 2) / angle_sum;
  size_t curr_idx = 0;
  for (size_t side = 0; side < num_sides; ++side) {
    const auto& curv = curve_loop[side];
    size_t num_pts = curv.size();

    for (size_t ii = 0; ii < num_pts - 2; ++ii) {
      normalized_interior_angles[curr_idx] *= scaling_factor;
      curr_idx++;
    }
    normalized_interior_angles[curr_idx] *= scaling_factor;
    curr_idx++;
  }

  // Integrating normalized angles
  curr_idx = 0;
  double curr_angle = 0.0;
  Vector curr_pos(0.0, 0.0, 0.0);
  for (size_t side = 0; side < num_sides; ++side) {
    const auto& curve = curve_loop[side];
    size_t num_pts = curve.size();
    DoubleVector lengths;

    GeomUtil::getEdgeLengths(curve, lengths);

    curr_angle = integrateCurvature(
      lengths,
      M_PI - DoubleVector(
        normalized_interior_angles.begin() + curr_idx,
        normalized_interior_angles.begin() + curr_idx + num_pts - 2),
      normalized_curv[side],
      curr_angle,
      curr_pos
    );
    curr_angle += M_PI - normalized_interior_angles[curr_idx + num_pts - 2];
    curr_pos = normalized_curv[side].back();
    curr_idx += num_pts - 1;
  }

  return normalized_curv;
}


Curve3D LoopFlattener::closeLoop(
  const Curve3D& curve_loop
)
{
  Curve3D closed_loop = curve_loop;
  size_t num_sides = curve_loop.size();
  double lengths_sum = 0.0;
  DoubleVector lengths_sum_cumm(1, 0.0);

  for (size_t side = 0; side < num_sides; ++side) {
    const auto& curv = curve_loop[side];
    size_t num_pts = curv.size();

    for (size_t ii = 0; ii < num_pts - 1; ++ii) {
      lengths_sum += getEdgeLength(curv, ii);
      lengths_sum_cumm.push_back(lengths_sum);
    }
  }
  lengths_sum_cumm /= lengths_sum;
  Vector diff = curve_loop.back().back() - curve_loop.front().front();

  size_t curr_idx = 0;
  for (size_t side = 0; side < num_sides; ++side) {
    const auto& curv = curve_loop[side];
    size_t num_pts = curv.size();

    for (size_t ii = 0; ii < num_pts; ++ii) {
      closed_loop[side][ii] -= lengths_sum_cumm[curr_idx + ii] * diff;
    }
    curr_idx += num_pts - 1;
  }

  return closed_loop;
}