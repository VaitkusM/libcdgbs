#pragma once

#include <Eigen/Core>

namespace GeomUtil {
  using Vector = Eigen::Vector3d;
  using DoubleVector = std::vector<double>;

  double getAngle(
    Vector v1,
    Vector v2,
    bool   turning = false,
    bool   signed_angle = true
  );

  double getAngle(
    Vector v1,
    Vector v2,
    Vector nn,
    bool   turning = false,
    bool   signed_angle = true
  );

  double unsignAngle(double angle);

  Vector projectToPlane(Vector vv, Vector nn);

  double getVoronoiLength(
    const std::vector<Vector>& curve_pts,
    size_t          idx,
    bool            loop = false
  );

  double getEdgeLength(
    const std::vector<Vector>& curve_pts,
    size_t idx,
    bool loop = false
  );

  void getEdgeLengths(
    const std::vector<Vector>& curve_pts,
    DoubleVector& edge_lengths,
    bool loop = false
  );


  Vector vectorArea(const std::vector<Vector>& pts);
  Vector vectorArea(const std::vector<std::vector<Vector> >& loop);
  Vector centroid(const std::vector<Vector>& pts);
  Vector centroid(const std::vector<std::vector<Vector> >& loop);
  void shiftByVector(std::vector<Vector>& loop, Vector vv);
  void shiftByVector(std::vector<std::vector<Vector> >& loop, Vector vv);
  void alignPointSets(std::vector<Vector>& source, const std::vector<Vector>& target);
  void alignPointSets(std::vector<std::vector<Vector> >& source, const std::vector<std::vector<Vector> >& target);
  void alignCurveLoopPCA(const std::vector<std::vector<Vector> >& pts, std::vector<std::vector<Vector> >& pts_aligned, bool shift_cog);
  
  void completeOrthoFrame(
    Vector  v0,
    Vector  v1,
    Vector& v2,
    bool    right_handed = true
  );

  void getPrincipalAxes(
    const std::vector<std::vector<Vector> >& curves,
    Vector& xx,
    Vector& yy,
    Vector& zz,
    Vector& cog
  );

  inline size_t circular_index(size_t i, int offset, size_t n) {
    return (i + n + (offset % static_cast<int>(n))) % n;
  }

  DoubleVector operator+(const DoubleVector& array, double scalar);
  DoubleVector operator+(double scalar, const DoubleVector& array);
  DoubleVector& operator+=(DoubleVector& array, double scalar);
  DoubleVector operator-(const DoubleVector& array, double scalar);
  DoubleVector operator-(double scalar, const DoubleVector& array);
  DoubleVector& operator-=(DoubleVector& array, double scalar);
  DoubleVector operator*(const DoubleVector& array, double scalar);
  DoubleVector operator*(double scalar, const DoubleVector& array);
  DoubleVector& operator*=(DoubleVector& array, double scalar);
  DoubleVector operator/(const DoubleVector& array, double scalar);
  DoubleVector operator/(double scalar, const DoubleVector& array);
  DoubleVector& operator/=(DoubleVector& array, double scalar);
}