#include "libcdgbs/SurfGBS.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

using namespace libcdgbs;

bool SurfGBS::readNGBS(const std::string& filename, double target_length)
{
  std::ifstream in(filename);
  if (!in) {
    std::cerr << "Cannot open file\n";
    return false;
  }

  // Read in first line starting with #
  std::string line;
  std::getline(in, line);

  int dom_type;
  double dom_res, bd_res;
  in >> dom_type >> dom_res >> bd_res;

  int autohs;
  double ahsa, ahsr;
  in >> autohs >> ahsa >> ahsr;

  int mu, cvmode;
  in >> mu >> cvmode;

  in >> num_loops;
  ribbons.resize(num_loops);
  num_sides.resize(num_loops);
  num_rows.resize(num_loops);
  num_cols.resize(num_loops);
  for (size_t loop = 0; loop < num_loops; ++loop) {
    in >> num_sides[loop];
    num_rows[loop].resize(num_sides[loop]);
    num_cols[loop].resize(num_sides[loop]);
    for (size_t side = 0; side < num_sides[loop]; ++side) {
      int degH, degS, nLayers;
      in >> degS >> degH >> nLayers;
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

      // read h-splits
      double hsl, hsr;
      in >> hsl >> hsr;

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
          ctrl[j * num_rows[loop][side] + i] = { x, y, z };
        }
      }

      // construct & store the B‑spline surface
      ribbons[loop].emplace_back(degU, degV, knotsU, knotsV, ctrl);

      // Ribbon CPs
      int rcp_n, rcp_deg;
      double rcp_w;
      in >> rcp_n >> rcp_w >> rcp_deg;
      if(rcp_n > 0) {
        std::vector<Eigen::Vector3d> rcp(rcp_n);
        for (size_t i = 0; i < rcp_n; ++i) {
          double x, y, z;
          in >> x >> y >> z;
          rcp[i] = { x, y, z };
        }
      }

      // alpha-beta splines
      int alpha_deg, beta_deg;
      in >> alpha_deg; beta_deg = alpha_deg;
      // read full knot‐vector line for alpha
      in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      knotLine.clear();
      std::getline(in, knotLine);
      ksi = std::istringstream(knotLine);
      Geometry::DoubleVector alpha_knots{
          std::istream_iterator<double>(ksi),
          std::istream_iterator<double>()
      };

      auto alpha_ncp = alpha_knots.size() - alpha_deg - 1;
      std::vector<double> alpha_cps(alpha_ncp), beta_cps(alpha_ncp);
      for (size_t i = 0; i < alpha_ncp; ++i) {
        double alpha, beta, zero;
        in >> alpha >> beta >> zero;
        alpha_cps[i] = alpha;
        beta_cps[i] = beta;
      }

      // Control Vectors
      int cv_n;
      in >> cv_n;
      if (cv_n > 0) {
        std::vector<Eigen::Vector3d> cv_cps(cv_n);
        std::vector<double> cv_u(cv_n);
        for (size_t i = 0; i < cv_n; ++i) {
          double u, x, y, z;
          in >> u >> x >> y >> z;
          cv_u[i] = u;
          cv_cps[i] = { x, y, z };
        }
      }
    }
  }

  load_ribbons(ribbons, target_length);

  return true;
}