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