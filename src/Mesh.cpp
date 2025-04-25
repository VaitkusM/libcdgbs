#include "libcdgbs/Mesh.hpp"

using namespace libcdgbs;

Mesh::FaceHandle Mesh::findClosestFace(Point pt) const
{
  FaceHandle closest_face;
  double min_dist = std::numeric_limits<double>::max();

  for (auto ff : faces()) {
    double dist = dist2Face(pt, ff);
    if (dist < min_dist) {
      min_dist = dist;
      closest_face = ff;
    }
  }

  return closest_face;
}

double Mesh::dist2Face(Point pt, FaceHandle ff) const
{
  return (pt - calc_face_centroid(ff)).norm();
}

void Mesh::barycentricCoordinates(
  const Point& pt,
  const Point& p1,
  const Point& p2,
  const Point& p3,
  double& u, double& v, double& w)
{
  auto v0 = p2 - p1;
  auto v1 = p3 - p1;
  auto v2 = pt - p1;

  auto d00 = v0.dot(v0);
  auto d01 = v0.dot(v1);
  auto d11 = v1.dot(v1);
  auto d20 = v2.dot(v0);
  auto d21 = v2.dot(v1);

  auto denom = d00 * d11 - d01 * d01;
  const double epsilon = 1e-10;
  if (std::abs(denom) < epsilon) {
    u = v = w = 0.0; // Degenerate or nearly degenerate case
    return;
  }

  u = (d11 * d20 - d01 * d21) / denom;
  v = (d00 * d21 - d01 * d20) / denom;
  w = 1.0 - u - v;
}