#include "libcdgbs/SurfGBS.hpp"
#include <vector>
#include <cmath>
#include <cstddef>
#include <stdexcept>

using namespace libcdgbs;

/*
Periodic uniform cubic B-spline basis with C1 continuity achieved by double interior knots.

INPUT:
  u in [0,1) (treated periodically: u==1 maps to u==0)
  N = number of Bézier segments around the loop (N>=2 recommended)

OUTPUT:
  w: dense vector of length M=2N, containing basis values such that
       C(u) = sum_{i=0..M-1} w[i] * P_i
  where P_i are the periodic B-spline control points.

INDEXING / Bézier correspondence (critical):
  Bézier segment k has control points Q_{k,0..3}.
  Define B-spline control points P_i (i mod 2N) by
      P_{2k}   = Q_{k,1}
      P_{2k+1} = Q_{k,2}
  Then (implied by C1 consistency):
      Q_{k,0} = 0.5*(P_{2k-1} + P_{2k})
      Q_{k,3} = 0.5*(P_{2k+1} + P_{2k+2})

ACTIVE BASIS ON SEGMENT k:
  Indices {2k-1, 2k, 2k+1, 2k+2} (mod 2N).

LOCAL PARAMETER:
  k = floor(N*u) in {0,...,N-1}
  t = N*u - k in [0,1)

BASIS VALUES (in terms of cubic Bernstein b0..b3):
  w[2k-1] = 0.5*b0
  w[2k]   = 0.5*b0 + b1
  w[2k+1] = b2 + 0.5*b3
  w[2k+2] = 0.5*b3
*/
std::vector<double> SurfGBS::periodicC1CubicBSplineBasis(double u, size_t N)
{
    if (N == 0) throw std::invalid_argument("N must be >= 1.");
    const std::size_t M = 2 * N;

    // Wrap u periodically to [0,1)
    double uw = u - std::floor(u);
    if (uw >= 1.0) uw = 0.0; // defensive

    // Segment index and local coordinate
    const double s = uw * static_cast<double>(N);
    std::size_t k = static_cast<std::size_t>(std::floor(s));
    if (k >= N) k = 0;                 // if uw extremely close to 1
    const double t = s - static_cast<double>(k);

    // Cubic Bernstein polynomials
    const double omt = 1.0 - t;
    const double b0 = omt * omt * omt;
    const double b1 = 3.0 * t * omt * omt;
    const double b2 = 3.0 * t * t * omt;
    const double b3 = t * t * t;

    std::vector<double> w(M, 0.0);

    // Base index: P_{2k} corresponds to Bézier Q_{k,1}
    const std::size_t base = (2 * k) % M;

    auto add = [&](std::ptrdiff_t di, double val) {
        const std::size_t idx = static_cast<std::size_t>((static_cast<std::ptrdiff_t>(base) + di + static_cast<std::ptrdiff_t>(M)) % static_cast<std::ptrdiff_t>(M));
        w[idx] += val; // += to be safe even if indices wrap (e.g. very small N)
    };

    add(-1, 0.5 * b0);        // index 2k-1
    add( 0, 0.5 * b0 + b1);   // index 2k
    add( 1, b2 + 0.5 * b3);   // index 2k+1
    add( 2, 0.5 * b3);        // index 2k+2

    return w;
}
