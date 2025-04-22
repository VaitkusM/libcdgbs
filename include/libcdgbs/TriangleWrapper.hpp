#pragma once
#include <vector>
#include <Eigen/Core>
#include <libcdgbs/Mesh.hpp>

namespace libcdgbs {

  class TriangleWrapper {
  public:
    /// loops: each loop is a sequence of sub‑curves (each a polyline)
    Mesh triangulate_loop(
      std::vector<std::vector<std::vector<Eigen::Vector3d>>> const& loops,
      double L_target
    );

  };

} // namespace libcdgbs
