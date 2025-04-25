#include "libcdgbs/SurfGBS.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

using namespace libcdgbs;

bool SurfGBS::readMLP(const std::string& filename, double target_length)
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
    num_rows[loop].resize(num_sides[loop], 2); // MLP format support only 2 rows
    num_cols[loop].resize(num_sides[loop]);
    for (size_t side = 0; side < num_sides[loop]; ++side) {
      in >> num_cols[loop][side];                // <# of ribbons in loop li>
      int degV = num_rows[loop][side] - 1;   // vertical degree
      int degU = num_cols[loop][side] - 1;          // horizontal degree

      // build a uniform, clamped knot‐vector in U
      Geometry::DoubleVector knotsU;
      knotsU.insert(knotsU.end(), degU + 1, 0.0);
      knotsU.insert(knotsU.end(), degU + 1, 1.0);

      // build a uniform, clamped knot‐vector in V
      Geometry::DoubleVector knotsV;
      knotsV.insert(knotsV.end(), degV + 1, 0.0);
      knotsV.insert(knotsV.end(), degV + 1, 1.0);

      // read the (nCtrlV × nCtrlU) grid of 3D points
      Geometry::PointVector ctrl;
      ctrl.resize(num_cols[loop][side] * num_rows[loop][side]);
      for (size_t i = 0; i < num_cols[loop][side]; ++i) {
        double dummy1, dummy2, x0, y0, z0, x1, y1, z1;
        in >> dummy1 >> dummy2 >> x0 >> y0 >> z0 >> x1 >> y1 >> z1;
        ctrl[i * num_rows[loop][side]] = { x0, y0, z0 };
        ctrl[i * num_rows[loop][side] + 1] = { x1, y1, z1 };
      }

      // construct & store the B‑spline surface
      ribbons[loop].emplace_back(degU, degV, knotsU, knotsV, ctrl);
    }
  }

  load_ribbons(ribbons, target_length);

  return true;
}