#include <fstream>
#include <iomanip>
#include "libcdgbs/SurfGBS.hpp"

using namespace libcdgbs;

bool SurfGBS::writeMGBS(const std::string& filename) const
{
  std::ofstream out(filename);
  if (!out) return false;

  out << ribbons.size() << '\n';                         // <# of loops>

  for (size_t loop = 0; loop < ribbons.size(); ++loop) {
    const auto& sides = ribbons[loop];
    out << sides.size() << '\n';                         // <# of ribbons in loop li>

    for (size_t side = 0; side < sides.size(); ++side) {
      const auto& rib = sides[side];

      const size_t degS = rib.basisU().degree();        // s-degree (U degree)
      const size_t nCols = rib.basisU().knots().size() - degS - 1;

      const size_t degV = rib.basisV().degree();        // actual V degree
      const size_t nRows = rib.basisV().knots().size() - degV - 1; // rows (layers)
      const size_t degH = deg_h[loop][side];  // h-degree (must equal nLayers-1 for your reader)

      // Header line for this side: <h-degree> <s-degree> <#layers> 0
      out << degH << ' ' << degS << ' ' << nRows << " 0\n";

      // Full U knot-vector on one line (your reader consumes the entire line into doubles)
      out << std::setprecision(17);
      const auto& knotsU = rib.basisU().knots();
      for (size_t k = 0; k < knotsU.size(); ++k) {
        if (k) out << ' ';
        out << knotsU[k];
      }
      out << '\n';

      // Control points: row-major [row (V)][col (U)], matching readMGBS' nested loops
      for (size_t row = 0; row < nRows; ++row) {
        for (size_t col = 0; col < nCols; ++col) {
          const auto cp = rib.controlPoint(col, row);
          out << cp[0] << ' ' << cp[1] << ' ' << cp[2] << '\n';
        }
      }
    }
  }

  // Tail ints (your reader expects two ints and doesn't otherwise use them)
  out << 400 << ' ' << 200 << '\n';

  return true;
}