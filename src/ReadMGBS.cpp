#include "libcdgbs/SurfGBS.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

using namespace libcdgbs;

bool SurfGBS::readMGBS(const std::string& filename, double target_length)
{
  std::ifstream in(filename);
  if (!in) {
    std::cerr << "Cannot open file\n";
    return false;
  }

  std::vector<std::vector<Ribbon> > ribbons;
  size_t num_loops;
  std::vector<size_t> num_sides;
  std::vector<std::vector<size_t>> num_rows;
  std::vector<std::vector<size_t>> num_cols;

  ribbons.clear();
  num_sides.clear();
  num_rows.clear();
  num_cols.clear();

  in >> num_loops;                          // <# of loops>
  ribbons.resize(num_loops);
  num_sides.resize(num_loops);
  num_rows.resize(num_loops);
  num_cols.resize(num_loops);
  for (size_t loop = 0; loop < num_loops; ++loop) {
    in >> num_sides[loop];                // <# of ribbons in loop li>
    num_rows[loop].resize(num_sides[loop]);
    num_cols[loop].resize(num_sides[loop]);
    for (size_t side = 0; side < num_sides[loop]; ++side) {
      int degH, degS, nLayers, zero;
      in >> degH >> degS >> nLayers >> zero;
      // ensure consistency: degH == nLayers-1
      int degV = nLayers - 1;   // vertical degree
      int degU = degS;          // horizontal degree

      // read full knot‐vector line for U
      in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      std::string knotLine;
      std::getline(in, knotLine);
      std::istringstream ksi(knotLine);
      Geometry::DoubleVector knotsU{
          std::istream_iterator<double>(ksi),
          std::istream_iterator<double>()
      };

      // build a uniform, clamped knot‐vector in V
      Geometry::DoubleVector knotsV;
      knotsV.insert(knotsV.end(), degV + 1, 0.0);
      knotsV.insert(knotsV.end(), degV + 1, 1.0);

      // compute control‐point counts
      num_cols[loop][side] = knotsU.size() - degU - 1;
      num_rows[loop][side] = nLayers;

      // read the (nCtrlV × nCtrlU) grid of 3D points
      Geometry::PointVector ctrl;
      ctrl.resize(num_cols[loop][side] * num_rows[loop][side]);
      for (size_t i = 0; i < num_rows[loop][side]; ++i) {
        for (size_t j = 0; j < num_cols[loop][side]; ++j) {
          double x, y, z;
          in >> x >> y >> z;
          ctrl[j*num_rows[loop][side] + i] = {x, y, z};
        }
      }

      // construct & store the B‑spline surface
      ribbons[loop].emplace_back(degU, degV, knotsU, knotsV, ctrl);
    }
  }

  int sideRes, meshRes;
  // read the last two ints
  if (!(in >> sideRes >> meshRes)) {
    std::cerr << "Missing side/mesh resolution\n";
    return false;
  }

  load_ribbons(ribbons, target_length);

  return true;
}