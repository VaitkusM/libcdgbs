#include "libcdgbs/LoopFlattener.hpp"
#include "libcdgbs/GeomUtil.hpp"
#include <iostream>
#include <Eigen/Cholesky>
using namespace libcdgbs;
using namespace GeomUtil;

using Vec3 = Eigen::Vector3d;
using Vec2 = Eigen::Vector2d;
using SubCurve3D = std::vector<Vec3>;
using Curve3D = std::vector<SubCurve3D>;
using Curves3D = std::vector<Curve3D>;
using Loop2D = LoopFlattener::Loop2D;
using DistConstraint = LoopFlattener::DistConstraint;

using Eigen::Vector2d;
using Eigen::VectorXd;
using Eigen::MatrixXd;


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



Curve3D LoopFlattener::developLoop(const Curve3D& curves, const Curve3D& normals, bool flipped) {
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
      // std::cout << "Interior angle (" << side << "," << ii << "):"  << interior_angles_for_side[ii] << std::endl;
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
      projectToPlane(curv_p1[1] - curv_p1[0], normals[side].back()).normalized(),
      projectToPlane(curv[curv.size() - 2] - curv[curv.size() - 1], normals[side].back()).normalized(),
      normals[side].back(),
      false,
      false
    );
    interior_angles.push_back(corner_angle);
    // std::cout << "Corner angle: " << corner_angle << std::endl;
    angle_sum += corner_angle;
  }

  // std::cout << "Angle sum: " << angle_sum << std::endl;

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


static inline Eigen::Matrix2d Rot(double th) 
{
  double c = std::cos(th), s = std::sin(th);
  Eigen::Matrix2d R;
  R << c, -s,
        s,  c;
  return R;
}

double LoopFlattener::sigmoid(double u) 
{
  // numerically stable sigmoid
  if (u >= 0.0) {
    double e = std::exp(-u);
    return 1.0 / (1.0 + e);
  } else {
    double e = std::exp(u);
    return e / (1.0 + e);
  }
}

Vector2d LoopFlattener::transformPoint(const Loop2D& L, int k) 
{
  const Vector2d& p = L.ref[k];
  if (L.fixed) return p;
  const double s = LoopFlattener::sigmoid(L.u);
  return s * Rot(L.theta) * p + Vector2d(L.tx, L.ty);
}

// 2x4 Jacobian of transformed point wrt [tx, ty, theta, u]
static inline Eigen::Matrix<double,2,4> dPoint_dParams(const Loop2D& L, int k) 
{
  Eigen::Matrix<double,2,4> J;
  J.setZero();
  if (L.fixed) return J;

  const Vector2d& p = L.ref[k];
  const double s = LoopFlattener::sigmoid(L.u);
  const double ds_du = s * (1.0 - s);
  const Eigen::Matrix2d R = Rot(L.theta);

  J.col(0) = Vector2d(1.0, 0.0); // d/dtx
  J.col(1) = Vector2d(0.0, 1.0); // d/dty

  Eigen::Matrix2d K; // 90deg rot operator
  K << 0.0, -1.0,
       1.0,  0.0;
  J.col(2) = s * (R * (K * p));   // d/dtheta
  J.col(3) = ds_du * (R * p);         // d/du

  return J;
}

struct VarIndexing {
    std::vector<int> offset; // per loop: -1 if fixed, else starting index
    int nvars = 0;
};

static VarIndexing buildIndexing(const std::vector<Loop2D>& loops) 
{
  VarIndexing VI;
  VI.offset.assign(loops.size(), -1);
  int off = 0;
  for (int i = 0; i < (int)loops.size(); ++i) {
    if (!loops[i].fixed) {
      VI.offset[i] = off;
      off += 4;
    }
  }
  VI.nvars = off;
  return VI;
}

static double cost(const std::vector<Loop2D>& loops,
                   const std::vector<DistConstraint>& C,
                   bool relative = true) 
{
  double sse = 0.0;
  for (const auto& c : C) {
    Vector2d a = LoopFlattener::transformPoint(loops[c.loopA], c.idxA);
    Vector2d b = LoopFlattener::transformPoint(loops[c.loopB], c.idxB);
    auto d_target = c.targetDist;
    if(d_target < 1e-8) d_target = 1e-8; // avoid singularity
    double r = (a - b).norm() - d_target;
    if(relative) {
      r /= d_target; // for clarity
    }
    sse += r * r;
  }
  return sse;
}

static void applyStep(std::vector<Loop2D>& loops, const VarIndexing& VI,
                      const VectorXd& dx, double alpha) 
{
  for (int i = 0; i < (int)loops.size(); ++i) {
    int off = VI.offset[i];
    if (off < 0) continue;
    loops[i].tx    += alpha * dx[off + 0];
    loops[i].ty    += alpha * dx[off + 1];
    loops[i].theta += 0; //alpha * dx[off + 2];
    loops[i].u     += alpha * dx[off + 3];

    // Optional: keep theta in a reasonable range (not required, but helps debugging)
    if (loops[i].theta >  M_PI) loops[i].theta -= 2.0 * M_PI;
    if (loops[i].theta < -M_PI) loops[i].theta += 2.0 * M_PI;
  }
}

// Damped Gauss–Newton (LM-style) with simple backtracking
bool LoopFlattener::solveHolePlacements(std::vector<Loop2D>& loops,
                                        const std::vector<DistConstraint>& constraints,
                                        int maxIters,
                                        double lambda0,
                                        double tolStep,
                                        double tolCostRel,
                                        bool relative) 
{
  VarIndexing VI = buildIndexing(loops);
  const int n = VI.nvars;
  if (n == 0) return true;

  double lambda = lambda0;
  double oldCost = cost(loops, constraints, relative);

  // std::cout << " Iter " << -1 << ": cost = " << oldCost << ", lambda = " << lambda << std::endl;
  // // std::cout << "  Step norm = " << dx.norm() << ", alpha = " << alpha << std::endl;
  // // std::cout << "  Params:";
  // for (int i = 0; i < (int)loops.size(); ++i) {
  //   if (!loops[i].fixed) {
  //     std::cout << "   Loop " << i << ": t = (" << loops[i].tx << ", " << loops[i].ty << "); theta = " << loops[i].theta << "; s = " << sigmoid(loops[i].u) << std::endl;
  //   }
  // }

  for (int it = 0; it < maxIters; ++it) {
    MatrixXd H = MatrixXd::Zero(n, n);
    VectorXd g = VectorXd::Zero(n);

    // Build normal equations H = J^T J, g = J^T r
    for (const auto& c : constraints) {
      const Loop2D& LA = loops[c.loopA];
      const Loop2D& LB = loops[c.loopB];

      Vector2d xa = transformPoint(LA, c.idxA);
      Vector2d xb = transformPoint(LB, c.idxB);

      Vector2d v = xa - xb;
      double dist = v.norm();
      if (dist < 1e-12) dist = 1e-12;
      Vector2d u = v / dist;

      double r = dist - c.targetDist;

      Eigen::RowVector4d JrA = Eigen::RowVector4d::Zero();
      Eigen::RowVector4d JrB = Eigen::RowVector4d::Zero();

      int offA = VI.offset[c.loopA];
      int offB = VI.offset[c.loopB];

      if (!LA.fixed) {
        auto Jxa = dPoint_dParams(LA, c.idxA);      // 2x4
        for (int j = 0; j < 4; ++j) JrA[j] = u.dot(Jxa.col(j));
      }
      if (!LB.fixed) {
        auto Jxb = dPoint_dParams(LB, c.idxB);      // 2x4
        for (int j = 0; j < 4; ++j) JrB[j] = -u.dot(Jxb.col(j));
      }

      // r_rel = (dist - di)/di  (same as dist/di - 1)
      if(relative) {
        const double dmin = 1e-8; // choose based on your scale; must be > 0
        double di  = std::max(c.targetDist, dmin);
        double inv = 1.0 / di;

        r   *= inv;
        JrA *= inv;
        JrB *= inv;
      }

      if (offA >= 0) {
        g.segment<4>(offA) += JrA.transpose() * r;
        H.block<4,4>(offA, offA) += JrA.transpose() * JrA;
      }
      if (offB >= 0) {
        g.segment<4>(offB) += JrB.transpose() * r;
        H.block<4,4>(offB, offB) += JrB.transpose() * JrB;
      }
      if (offA >= 0 && offB >= 0) {
        H.block<4,4>(offA, offB) += JrA.transpose() * JrB;
        H.block<4,4>(offB, offA) += JrB.transpose() * JrA;
      }
    }

    // Damping: (H + lambda I) dx = -g  (classic LM-style “damped least squares”) :contentReference[oaicite:1]{index=1}
    MatrixXd Hd = H;
    Hd.diagonal().array() += lambda;

    VectorXd dx = Hd.ldlt().solve(-g); // Eigen’s LDLT solver for symmetric systems :contentReference[oaicite:2]{index=2}
    if (!dx.allFinite()) {
      lambda *= 10.0;
      continue;
    }

    if (dx.norm() < tolStep) return true;

    // Backtracking line search on the step
    const std::vector<Loop2D> saved = loops;
    bool accepted = false;
    double alpha = 1.0;

    for (int ls = 0; ls < 12; ++ls) {
      loops = saved;
      applyStep(loops, VI, dx, alpha);

      double newCost = cost(loops, constraints, relative);
      if (newCost <= oldCost) {
        double rel = std::abs(oldCost - newCost) / std::max(1.0, oldCost);
        oldCost = newCost;
        lambda = std::max(lambda * 0.3, 1e-12);
        accepted = true;
        if (rel < tolCostRel) return true;
        break;
      }
      alpha *= 0.5;
    }

    if (!accepted) {
      loops = saved;
      lambda *= 10.0; // increase damping when we cannot find a decreasing step
    }
    // std::cout << " Iter " << it << ": cost = " << oldCost << ", lambda = " << lambda << std::endl;
    // // std::cout << "  Step norm = " << dx.norm() << ", alpha = " << alpha << std::endl;
    // // std::cout << "  Params:";
    // for (int i = 0; i < (int)loops.size(); ++i) {
    //   if (!loops[i].fixed) {
    //     std::cout << "   Loop " << i << ": t = (" << loops[i].tx << ", " << loops[i].ty << "); theta = " << loops[i].theta << "; s = " << sigmoid(loops[i].u) << std::endl;
    //   }
    // }
  }
  return false; // max iters reached
}