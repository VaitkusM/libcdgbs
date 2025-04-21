#include "libcdgbs/GeomUtil.hpp"
#include <Eigen/Geometry>

namespace GeomUtil {

  double getAngle(
    Vector v1,
    Vector v2,
    bool   turning,
    bool   signed_angle
  )
  {
    double s = (v1.cross(v2)).norm();
    double c = v1.dot(v2);
    double angle = atan2(s, c);
    angle = turning ? M_PI - angle : angle;
    return signed_angle ? angle : unsignAngle(angle);
  }

  double getAngle(Vector v1,
    Vector v2,
    Vector nn,
    bool   turning,
    bool   signed_angle
  )
  {
    double s = (v1.cross(v2)).norm();
    double c = v1.dot(v2);
    double angle = ((v1.cross(v2)).dot(nn)) >= 0 ? atan2(s, c) : -atan2(s, c);
    angle = turning ? M_PI - angle : angle;
    return signed_angle ? angle : unsignAngle(angle);
  }

  double unsignAngle(double angle) {
    while (angle < 0) {
      angle += 2.0 * M_PI;
    }
    return angle;
  }

  Vector projectToPlane(Vector vv, Vector nn)
  {
    nn.normalize();
    return vv - (vv.dot(nn)) * nn;
  }

  double getVoronoiLength(
    const std::vector<Vector>& curve_pts,
    size_t          idx,
    bool            loop
  )
  {
    size_t idx_m1 = circular_index(idx, -1, curve_pts.size());
    return (getEdgeLength(curve_pts, idx_m1, loop) +
      getEdgeLength(curve_pts, idx, loop)) / 2.0;
  }

  double getEdgeLength(
    const std::vector<Vector>& curve_pts,
    size_t idx,
    bool loop)
  {
    size_t num_pts = curve_pts.size();
    size_t idx_p1 = loop ? circular_index(idx, 1, num_pts) : (idx + 1);
    if (idx_p1 >= num_pts) {
      return 0.0;
    }
    return (curve_pts[idx_p1] - curve_pts[idx]).norm();
  }

  void getEdgeLengths(
    const std::vector <Vector>& curve_pts,
    DoubleVector& edge_lengths,
    bool            loop
  )
  {
    size_t num_edges = curve_pts.size() - (loop ? 0 : 1);

    edge_lengths.clear();
    edge_lengths.resize(num_edges);

    for (size_t ii = 0; ii < num_edges; ++ii) {
      edge_lengths[ii] = getEdgeLength(curve_pts, ii, loop);
    }
  }

  double integrateCurvature(
    const DoubleVector& lengths,
    const DoubleVector& curvatures,
    std::vector<Vector>& curve_2d,
    double              starting_angle,
    Vector              starting_pos
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



  // Convenience operators for DoubleVector

  DoubleVector operator+(const DoubleVector& array, double scalar)
  {
    DoubleVector new_array(array);
    for (auto& dd : new_array) {
      dd += scalar;
    }
    return new_array;
  }

  DoubleVector operator+(double scalar, const DoubleVector& array)
  {
    DoubleVector new_array(array);
    for (auto& dd : new_array) {
      dd += scalar;
    }
    return new_array;
  }

  DoubleVector& operator+=(DoubleVector& array, double scalar)
  {
    for (auto& dd : array) {
      dd += scalar;
    }
    return array;
  }

  DoubleVector operator-(const DoubleVector& array, double scalar)
  {
    DoubleVector new_array(array);
    for (auto& dd : new_array) {
      dd -= scalar;
    }
    return new_array;
  }

  DoubleVector operator-(double scalar, const DoubleVector& array)
  {
    DoubleVector new_array(array);
    for (auto& dd : new_array) {
      dd = scalar - dd;
    }
    return new_array;
  }

  DoubleVector& operator-=(DoubleVector& array, double scalar)
  {
    for (auto& dd : array) {
      dd -= scalar;
    }
    return array;
  }

  DoubleVector operator*(const DoubleVector& array, double scalar)
  {
    DoubleVector new_array(array);
    for (auto& dd : new_array) {
      dd *= scalar;
    }
    return new_array;
  }

  DoubleVector operator*(double scalar, const DoubleVector& array)
  {
    DoubleVector new_array(array);
    for (auto& dd : new_array) {
      dd *= scalar;
    }
    return new_array;
  }

  DoubleVector& operator*=(DoubleVector& array, double scalar)
  {
    for (auto& dd : array) {
      dd *= scalar;
    }
    return array;
  }

  DoubleVector operator/(const DoubleVector& array, double scalar)
  {
    DoubleVector new_array(array);
    for (auto& dd : new_array) {
      dd /= scalar;
    }
    return new_array;
  }

  DoubleVector operator/(double scalar, const DoubleVector& array)
  {
    DoubleVector new_array(array);
    for (auto& dd : new_array) {
      dd = scalar / dd;
    }
    return new_array;
  }

  DoubleVector& operator/=(DoubleVector& array, double scalar)
  {
    for (auto& dd : array) {
      dd /= scalar;
    }
    return array;
  }
}



