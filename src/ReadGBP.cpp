#include "libcdgbs/SurfGBS.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

using namespace libcdgbs;

bool SurfGBS::readGBP(const std::string& filename, double target_length)
{
  std::ifstream in(filename);
  if (!in) {
    std::cerr << "Cannot open file\n";
    return false;
  }

  std::vector<std::vector<Ribbon> > ribbons;
  const size_t num_loops = 1;
  std::vector<size_t> num_sides;
  std::vector<std::vector<size_t>> num_rows;
  std::vector<std::vector<size_t>> num_cols;

  ribbons.clear();
  num_sides.clear();
  num_rows.clear();
  num_cols.clear();

  ribbons.resize(1);
  num_sides.resize(1);
  num_rows.resize(1);
  num_cols.resize(1);

  int du;
  in >> num_sides[0] >> du;                // <# of ribbons>
  num_rows[0].resize(num_sides[0]);
  num_cols[0].resize(num_sides[0]);

  Eigen::Vector3d ccp;
  in >> ccp(0) >> ccp(1) >> ccp(2);   // <ccp_x> <ccp_y> <ccp_z>

  const size_t nl = (du + 1) / 2;
  const size_t dv = nl - 1;
  const size_t ncp = num_sides[0] * (1 + du / 2) * nl;
  const size_t nc = du + 1;

  // build a uniform, clamped knot‐vector in U
  Geometry::DoubleVector knotsU;
  knotsU.insert(knotsU.end(), du + 1, 0.0);
  knotsU.insert(knotsU.end(), du + 1, 1.0);

  // build a uniform, clamped knot‐vector in V
  Geometry::DoubleVector knotsV;
  knotsV.insert(knotsV.end(), dv + 1, 0.0);
  knotsV.insert(knotsV.end(), dv + 1, 1.0);

  std::vector<Eigen::Vector3d> cps(ncp);

  for (size_t i = 0; i < ncp; ++i) {
    double x, y, z;
    in >> x >> y >> z;
    cps[i] = { x, y, z };
  }

  std::vector<std::vector<std::vector<Eigen::Vector3d> > > corners(num_sides[0]);
  for (size_t side = 0; side < num_sides[0]; ++side) {
    corners[side].resize(nc - nl);
    for (size_t column = 0; column < nc - nl; ++column) {
      corners[side][column].resize(nl);
    }
  }

  int idx = 0;
  for (size_t layer = 0; layer < nl; ++layer) {
    for (size_t side = 0; side < num_sides[0]; ++side) {
      size_t side_p1 = (side + 1) % num_sides[0];
      for (int column = layer; column < int(nc) - 1 - int(layer); ++column) {
        if (column < int(nc) - int(nl)) {
          corners[side][column][layer] = cps[idx++];
        }
        else {
          corners[side_p1][layer][nc - 1 - column] = cps[idx++];
        }
      }
    }
  }

  for (size_t side = 0; side < num_sides[0]; ++side) {
    num_cols[0][side] = du + 1;
    num_rows[0][side] = nl;
    const size_t side_p1 = (side + 1) % num_sides[0];
    Geometry::PointVector ctrl;
    ctrl.resize(num_cols[0][side] * num_rows[0][side]);
    for (size_t layer = 0; layer < nl; ++layer) {
      for (size_t column = 0; column < nc; ++column) {
        if (column < nc- nl) {
          auto cp = corners[side][column][layer];
          ctrl[column * num_rows[0][side] + layer] = {cp(0), cp(1), cp(2)};
        }
        else {
          auto cp = corners[side_p1][layer][nc - 1 - column];
          ctrl[column * num_rows[0][side] + layer] = { cp(0), cp(1), cp(2) };
        }
      }
    }
    // construct & store the B‑spline surface
    ribbons[0].emplace_back(du, dv, knotsU, knotsV, ctrl);
  }

  load_ribbons(ribbons, target_length);

  return true;
}