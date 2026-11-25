// Regression test: automatic h-widths next to merged (multi-segment) sides.
//
// The domain is a planar 9 x 1 rectangle whose bottom side is given as
// three straight ribbons that join smoothly and are therefore merged:
//
//   side 0  left,   length 1                 (0,1) -> (0,0)
//   side 1  bottom, three merged segments    (0,0) -> (9,0)
//   side 2  right,  length 1                 (9,0) -> (9,1)
//   side 3  top,    length 9                 (9,1) -> (0,1)
//
// The automatic rule ends a neighbour's window where its arc length from
// the shared corner reaches twice the side's own length (here 2), at the
// first boundary sample past that point, and stores the window as a
// fraction of the neighbour's curve parameter. On a merged neighbour,
// restrict_params caps the window at the adjacent segment.

#include "libcdgbs/SurfGBS.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

namespace {

using libcdgbs::SurfGBS;

int failures = 0;

void check(bool ok, const std::string& what) {
    if (!ok) {
        ++failures;
        std::cerr << "FAILED: " << what << std::endl;
    }
}

void checkNear(double value, double expected, double tol, const std::string& what) {
    if (!(std::abs(value - expected) <= tol)) {
        ++failures;
        std::cerr << "FAILED: " << what << ": " << value << ", expected " << expected
                  << std::endl;
    }
}

// A straight cubic ribbon from (x0, y0) to (x1, y1) in the plane z = 0,
// uniformly parameterized, with the cross derivative pointing to the
// left of the direction of travel (into a counterclockwise loop).
SurfGBS::Ribbon straightRibbon(double x0, double y0, double x1, double y1) {
    const double dx = x1 - x0, dy = y1 - y0;
    const double len = std::sqrt(dx * dx + dy * dy);
    const double nx = -dy / len, ny = dx / len;
    const double depth = 0.3;
    Geometry::PointVector cpts;
    for (size_t c = 0; c < 4; ++c) {
        const double u = double(c) / 3.0;
        const double x = x0 + u * dx, y = y0 + u * dy;
        cpts.push_back({x, y, 0.0});
        cpts.push_back({x + depth * nx, y + depth * ny, 0.0});
    }
    const Geometry::DoubleVector knots_u{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
    const Geometry::DoubleVector knots_v{0.0, 0.0, 1.0, 1.0};
    return SurfGBS::Ribbon(3, 1, knots_u, knots_v, cpts);
}

// The rectangle with bottom segments of lengths a, b and c (a + b + c = 9).
bool build(SurfGBS& g, double a, double b, bool restrict, bool use_h_widths) {
    std::vector<std::vector<SurfGBS::Ribbon>> ribbons(1);
    ribbons[0].push_back(straightRibbon(0.0, 1.0, 0.0, 0.0));
    ribbons[0].push_back(straightRibbon(0.0, 0.0, a, 0.0));
    ribbons[0].push_back(straightRibbon(a, 0.0, a + b, 0.0));
    ribbons[0].push_back(straightRibbon(a + b, 0.0, 9.0, 0.0));
    ribbons[0].push_back(straightRibbon(9.0, 0.0, 9.0, 1.0));
    ribbons[0].push_back(straightRibbon(9.0, 1.0, 0.0, 1.0));
    SurfGBS::InputParams params;
    // Chosen so that no side length is a multiple of it: the sample
    // counts then do not depend on rounding in the measured lengths.
    params.target_length = 0.26;
    params.restrict_params = restrict;
    params.use_h_widths = use_h_widths;
    g.load_ribbons(ribbons, params);
    return g.compute_domain_boundary() && g.compute_domain_mesh() &&
           g.compute_local_parameters();
}

// Checks h of `side` along the merged side 1 against a linear ramp of
// width w from the shared corner (t = 0 for side 0, t = 1 for side 2),
// and that s of `side` is pinned at its end value inside the window.
void checkRamps(const SurfGBS& g, size_t side, double w, const std::string& label) {
    const auto& vertices = g.domain_boundary_vertices[0][1];
    const auto& params = g.domain_boundary_params[0][1];
    check(vertices.size() == params.size(), label + ": sample counts");
    double worst_h = 0.0, worst_s = 0.0;
    for (size_t i = 0; i < vertices.size() && i < params.size(); ++i) {
        const double t = params[i];
        const double d = (side == 0) ? t : 1.0 - t; // distance from the corner
        const size_t v = size_t(vertices[i].idx());
        const double expected_h = std::min(1.0, d / w);
        worst_h = std::max(worst_h, std::abs(g.h_coords[v][0][side] - expected_h));
        // Strictly inside the window (clear of rounding at its edge) and
        // away from the shared corner, which belongs to the side itself.
        if (d > 1e-9 && d < w - 1e-9) {
            const double expected_s = (side == 0) ? 1.0 : 0.0;
            worst_s = std::max(worst_s, std::abs(g.s_coords[v][0][side] - expected_s));
        }
    }
    checkNear(worst_h, 0.0, 1e-8, label + ": h ramp of side " + std::to_string(side));
    checkNear(worst_s, 0.0, 1e-8, label + ": s window of side " + std::to_string(side));
}

} // namespace

int main() {
    // With the target length 0.26 the bottom segments of lengths 1, 7, 1
    // get 5, 26 and 5 samples. From either corner, arc length 2 is first
    // reached at the fourth interval of the middle segment (1 + 4 * 7/25),
    // which is the curve parameter 1/3 + (4/25)/3 of the merged side.
    const double w_long = 1.0 / 3.0 + (4.0 / 25.0) / 3.0;
    // The top side (length 9) gets 34 samples; arc length 2 is first
    // reached at the eighth interval (8 * 9/33).
    const double w_top = 8.0 / 33.0;
    const double tol_w = 1e-5; // arc-length sampling bisects to about 1e-6

    { // The stored widths, which do not depend on restrict_params.
        SurfGBS g;
        check(build(g, 1.0, 7.0, true, true), "build (1, 7, 1), restrict");
        check(g.num_sides[0] == 4, "four sides after merging");
        check(g.num_segments[0][1] == 3, "side 1 has three segments");
        checkNear(g.h_widths[0][0][1], w_long, tol_w, "side 0, next (merged)");
        checkNear(g.h_widths[0][2][0], w_long, tol_w, "side 2, prev (merged)");
        checkNear(g.h_widths[0][0][0], w_top, tol_w, "side 0, prev");
        checkNear(g.h_widths[0][2][1], w_top, tol_w, "side 2, next");
        for (size_t side : {size_t(1), size_t(3)}) {
            checkNear(g.h_widths[0][side][0], 1.0, 0.0, "long side, prev");
            checkNear(g.h_widths[0][side][1], 1.0, 0.0, "long side, next");
        }
        // restrict_params caps the window at the adjacent segment.
        checkRamps(g, 0, 1.0 / 3.0, "(1, 7, 1) restrict");
        checkRamps(g, 2, 1.0 / 3.0, "(1, 7, 1) restrict");
    }
    { // Without restrict_params the automatic width applies.
        SurfGBS g;
        check(build(g, 1.0, 7.0, false, true), "build (1, 7, 1), no restrict");
        checkNear(g.h_widths[0][0][1], w_long, tol_w, "no restrict: side 0, next");
        checkRamps(g, 0, g.h_widths[0][0][1], "(1, 7, 1) no restrict");
        checkRamps(g, 2, g.h_widths[0][2][0], "(1, 7, 1) no restrict");
    }
    { // Without h-widths and without restrict_params: the whole side.
        SurfGBS g;
        check(build(g, 1.0, 7.0, false, false), "build (1, 7, 1), no widths");
        checkRamps(g, 0, 1.0, "(1, 7, 1) no widths");
        checkRamps(g, 2, 1.0, "(1, 7, 1) no widths");
    }
    { // Without h-widths, restrict_params alone gives the adjacent segment.
        SurfGBS g;
        check(build(g, 1.0, 7.0, true, false), "build (1, 7, 1), restrict only");
        checkRamps(g, 0, 1.0 / 3.0, "(1, 7, 1) restrict only");
        checkRamps(g, 2, 1.0 / 3.0, "(1, 7, 1) restrict only");
    }
    { // Segments of length 3 (11 samples each): arc length 2 is reached
      // inside the adjacent segment, at its seventh interval, so the
      // automatic width is narrower than the cap and applies as it is.
        const double w_short = (7.0 / 10.0) / 3.0;
        SurfGBS g;
        check(build(g, 3.0, 3.0, true, true), "build (3, 3, 3), restrict");
        check(g.num_segments[0][1] == 3, "(3, 3, 3): three segments");
        checkNear(g.h_widths[0][0][1], w_short, tol_w, "(3, 3, 3): side 0, next");
        checkNear(g.h_widths[0][2][0], w_short, tol_w, "(3, 3, 3): side 2, prev");
        checkRamps(g, 0, g.h_widths[0][0][1], "(3, 3, 3) restrict");
        checkRamps(g, 2, g.h_widths[0][2][0], "(3, 3, 3) restrict");
    }

    if (failures == 0) {
        std::cout << "h_widths_test: all checks passed" << std::endl;
        return 0;
    }
    std::cerr << "h_widths_test: " << failures << " check(s) failed" << std::endl;
    return 1;
}
