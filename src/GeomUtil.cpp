#include "libcdgbs/GeomUtil.hpp"
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>

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

  Vector vectorArea(const std::vector<Vector>& pts)
  {
    Vector va(0.0, 0.0, 0.0);
    const int num_pts = pts.size();
    for (int ii = 0; ii < num_pts; ++ii) {
      const auto v1 = pts[ii];
      const auto v2 = pts[circular_index(ii, 1, num_pts)];
      va += (v1.cross(v2));
    }
    return va;
  }

  Vector vectorArea(const std::vector<std::vector<Vector> >& loop)
  {
    Vector va(0.0, 0.0, 0.0);
    for (const auto& loop_pts : loop) {
      const int num_pts = loop_pts.size();
      for (int ii = 0; ii < num_pts; ++ii) {
        const auto v1 = loop_pts[ii];
        const auto v2 = loop_pts[circular_index(ii,1, num_pts)];
        va += (v1.cross(v2));
      }
    }
    return va;
  }


  Vector centroid(const std::vector<Vector>& pts)
  {
    Vector cog(0.0, 0.0, 0.0);
    for (auto pp : pts) {
      cog += pp;
    }
    cog /= pts.size();
    return cog;
  }

  Vector centroid(const std::vector<std::vector<Vector> >& loop)
  {
    Vector cog(0.0, 0.0, 0.0);
    size_t num = 0;
    for (const auto& loop_pts : loop) {
      for (auto pp : loop_pts) {
        cog += pp;
        num++;
      }
    }
    cog /= num;
    return cog;
  }

  void shiftByVector(std::vector<Vector>& loop, Vector vv)
  {
    for (auto& pp : loop) {
      pp += vv;
    }
  }

  void shiftByVector(std::vector<std::vector<Vector> >& loop, Vector vv)
  {
    for (auto& curve : loop) {
      for (auto& pp : curve) {
        pp += vv;
      }
    }
  }


  void alignPointSets(std::vector<Vector>& source, const std::vector<Vector>& target)
  {
    Eigen::Matrix2d cov, rot;
    cov = Eigen::Matrix2d::Zero();
    for (size_t ii = 0; ii < source.size(); ++ii) {
      Eigen::Vector2d ss, tt;
      ss(0) = source[ii][0];
      ss(1) = source[ii][1];
      tt(0) = target[ii][0];
      tt(1) = target[ii][1];
      cov += tt * ss.transpose();
    }
    Eigen::JacobiSVD<Eigen::Matrix2d> svd(cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
    rot = svd.matrixU() * svd.matrixV().transpose();
    for (size_t ii = 0; ii < source.size(); ++ii) {
      Eigen::Vector2d ss, rotss;
      ss(0) = source[ii][0];
      ss(1) = source[ii][1];
      rotss = rot * ss;
      source[ii][0] = rotss(0);
      source[ii][1] = rotss(1);
    }
  }

  void alignPointSets(std::vector<std::vector<Vector>>& source, const std::vector<std::vector<Vector>>& target)
  {
    Eigen::Matrix2d cov, rot;
    cov = Eigen::Matrix2d::Zero();
    for (size_t cc = 0; cc < source.size(); ++cc) {
      for (size_t pt = 0; pt < source[cc].size(); ++pt) {
        Eigen::Vector2d ss, tt;
        ss(0) = source[cc][pt][0];
        ss(1) = source[cc][pt][1];
        tt(0) = target[cc][pt][0];
        tt(1) = target[cc][pt][1];
        cov += tt * ss.transpose();
      }
    }
    Eigen::JacobiSVD<Eigen::Matrix2d> svd(cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
    rot = svd.matrixU() * svd.matrixV().transpose();
    for (size_t cc = 0; cc < source.size(); ++cc) {
      for (size_t pt = 0; pt < source[cc].size(); ++pt) {
        Eigen::Vector2d ss, rotss;
        ss(0) = source[cc][pt][0];
        ss(1) = source[cc][pt][1];
        rotss = rot * ss;
        source[cc][pt][0] = rotss(0);
        source[cc][pt][1] = rotss(1);
      }
    }
  }

  void completeOrthoFrame(
    Vector  v0,
    Vector  v1,
    Vector& v2,
    bool    right_handed
  )
  {
    v2 = (right_handed ? 1 : -1) * (v0.cross(v1)).normalized();
  }

  void getPrincipalAxes(
    const std::vector<std::vector<Vector> >& curves,
    Vector& xx,
    Vector& yy,
    Vector& zz,
    Vector& cog
  )
  {
    int num_pts = 0;
    for (const auto& pts : curves) {
      num_pts += pts.size();
    }

    cog = centroid(curves);
    Eigen::MatrixXd PP(num_pts, 3);
    size_t idx = 0;
    for (const auto& pts : curves) {

      for (size_t ii = 0; ii < pts.size(); ++ii) {
        PP(idx, 0) = pts[ii][0] - cog[0];
        PP(idx, 1) = pts[ii][1] - cog[1];
        PP(idx, 2) = pts[ii][2] - cog[2];
        idx++;
      }
    }
    Eigen::Matrix3d CC = PP.transpose() * PP;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver;
    solver.compute(CC);
    auto ev = solver.eigenvectors();

    xx[0] = ev(0, 2);
    xx[1] = ev(1, 2);
    xx[2] = ev(2, 2);

    yy[0] = ev(0, 1);
    yy[1] = ev(1, 1);
    yy[2] = ev(2, 1);

    zz[0] = ev(0, 0);
    zz[1] = ev(1, 0);
    zz[2] = ev(2, 0);


    xx = getAngle(xx, Vector(1.0, 0.0, 0.0)) >= 0 ? xx : -xx;

    completeOrthoFrame(Vector(0, 0, 1), xx, yy);
  }

  void alignCurveLoopPCA(const std::vector<std::vector<Vector> >& pts, std::vector<std::vector<Vector> >& pts_aligned, bool shift_cog)
  {

    Vector cog(0.0, 0.0, 0.0);
    Vector xx(0.0, 0.0, 0.0);
    Vector yy(0.0, 0.0, 0.0);
    Vector zz(0.0, 0.0, 0.0);

    getPrincipalAxes(pts, xx, yy, zz, cog);

    xx = -xx;
    yy = -yy;

    pts_aligned = pts;

    shiftByVector(pts_aligned, -cog);

    Eigen::Matrix3d rot, irot;
    rot(0, 0) = xx[0]; rot(0, 1) = yy[0]; rot(0, 2) = zz[0];
    rot(1, 0) = xx[1]; rot(1, 1) = yy[1]; rot(1, 2) = zz[1];
    rot(2, 0) = xx[2]; rot(2, 1) = yy[2]; rot(2, 2) = zz[2];
    irot = rot.transpose();

    for (size_t side = 0; side < pts.size(); ++side) {
      int num_pts = pts[side].size();
      for (int ii = 0; ii < num_pts; ++ii) {
        Eigen::Vector3d vv;
        vv(0) = pts_aligned[side][ii][0];
        vv(1) = pts_aligned[side][ii][1];
        vv(2) = pts_aligned[side][ii][2];
        auto vr = irot * vv;
        pts_aligned[side][ii][0] = vr(0);
        pts_aligned[side][ii][1] = vr(1);
        pts_aligned[side][ii][2] = vr(2);
      }
    }
    if (!shift_cog) {
      shiftByVector(pts_aligned, cog);
    }
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



