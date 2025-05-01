#include "libcdgbs/SurfGBS.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

using namespace libcdgbs;

bool SurfGBS::readCGB(const std::string& filename, double target_length)
{
  std::ifstream in(filename);
  if (!in) {
    std::cerr << "Cannot open file\n";
    return false;
  }

  int cw_type, dom_type;
  double par_dil, dom_tol;
  in >> cw_type >> dom_type >> par_dil >> dom_tol;

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

  in >> num_sides[0];                // <# of ribbons>
  num_rows[0].resize(num_sides[0]);
  num_cols[0].resize(num_sides[0]);

  Eigen::Vector3d ccp;
  in >> ccp(0) >> ccp(1) >> ccp(2);   // <ccp_x> <ccp_y> <ccp_z>

  for (size_t side = 0; side < num_sides[0]; ++side) {
    int degS, nLayers;
    in >> degS >> nLayers;
    // ensure consistency: degH == nLayers-1
    int degV = nLayers - 1;   // vertical degree
    int degU = degS;          // horizontal degree

    // build a uniform, clamped knot‐vector in U
    Geometry::DoubleVector knotsU;
    knotsU.insert(knotsU.end(), degU + 1, 0.0);
    knotsU.insert(knotsU.end(), degU + 1, 1.0);

    // build a uniform, clamped knot‐vector in V
    Geometry::DoubleVector knotsV;
    knotsV.insert(knotsV.end(), degV + 1, 0.0);
    knotsV.insert(knotsV.end(), degV + 1, 1.0);

    // compute control‐point counts
    num_cols[0][side] = degU + 1;
    num_rows[0][side] = nLayers;

    // read the (nCtrlV × nCtrlU) grid of 3D points
    Geometry::PointVector ctrl;
    ctrl.resize(num_cols[0][side] * num_rows[0][side]);
    for (size_t i = 0; i < num_rows[0][side]; ++i) {
      for (size_t j = 0; j < num_cols[0][side]; ++j) {
        double x, y, z;
        in >> x >> y >> z;
        ctrl[j * num_rows[0][side] + i] = { x, y, z };
      }
    }

    // construct & store the B‑spline surface
    ribbons[0].emplace_back(degU, degV, knotsU, knotsV, ctrl);
  }
  
  std::vector<double> multipliers(num_sides[0]);
  for (size_t i = 0; i < num_sides[0]; ++i) {
    in >> multipliers[i];
  }

  load_ribbons(ribbons, target_length);

  return true;
}