#include "libcdgbs/SurfGBS.hpp"
#include "libcdgbs/LoopFlattener.hpp"
#include "libcdgbs/TriangleWrapper.hpp"

using namespace libcdgbs;
using Vec3 = Eigen::Vector3d;
using Vec2 = Eigen::Vector2d;
using SubCurve3D = std::vector<Vec3>;
using Curve3D = std::vector<SubCurve3D>;
using Curves3D = std::vector<Curve3D>;

SurfGBS::SurfGBS()
{

}

int SurfGBS::load_ribbons(const std::vector<std::vector<Ribbon> >& ribbon_surfs)
{
  ribbons = ribbon_surfs;

  int ret_code = 0;

  compute_domain_boundary();
  compute_domain_mesh();
  compute_local_parameters();
  compute_blend_functions();
  evaluate_mesh();

  return true;
}

bool SurfGBS::compute_domain_boundary()
{
  Curves3D points(num_loops);
  Curves3D normals(num_loops);
  for (size_t loop = 0; loop < num_loops; loop) {
    points[loop].resize(num_sides[loop]);
    normals[loop].resize(num_sides[loop]);
    for (size_t side = 0; side < num_sides[loop]; ++side) {
      const auto& rib = ribbons[loop][side];
      const size_t res = side_res[loop][side];
      points[loop][side].resize(res);
      normals[loop][side].resize(res);
      for (size_t i = 0; i < res; ++i) {
        const double u = double(i) / res;

        Geometry::VectorMatrix duv;
        auto pt = rib.eval(u, 0.0, 1, duv);


        points[loop][side][i] = { pt[0], pt[1], pt[2] };
        auto du = duv[1][0];
        auto dv = duv[0][1];
        auto nn = (du ^ dv).normalized();
        normals[loop][side][i] = { nn[0], nn[1], nn[2] };
      }
    }
  }

  domain_boundary_curves = LoopFlattener::developCurves(points, normals);

  return true;
}

bool SurfGBS::compute_domain_mesh()
{
  auto triangle_wrapper = TriangleWrapper();
  meshDomain = triangle_wrapper.triangulate(domain_boundary_curves, side_res.front().front());

  domain_boundary_vertices.clear();
  domain_boundary_vertices.resize(num_loops);
  size_t idx_vtx = 0;
  for (size_t loop_idx = 0; loop_idx < domain_boundary_curves.size(); ++loop_idx) {
    domain_boundary_vertices[loop_idx].resize(num_sides[loop_idx]);
    auto const& loop = domain_boundary_curves[loop_idx];
    size_t first_idx = idx_vtx;
    for (size_t si = 0; si < loop.size(); ++si) {
      auto sub = loop[si];
      for (size_t i = 0; i < sub.size(); ++i) {
        auto vtx = meshDomain.vertex_handle((si == loop.size() - 1 && i == sub.size() - 1) ? first_idx : idx_vtx);
        domain_boundary_vertices[loop_idx][si].push_back(vtx);
        if (i < sub.size() - 1) {
          ++idx_vtx;
        }
      }
    }
  }

  return true;
}

bool SurfGBS::compute_local_parameters()
{
  if (!compute_harmonic_parameters()) {
    std::cout << "Error computing harmonic parameters" << std::endl;
    return false;
  }

  return true;
}

bool SurfGBS::compute_blend_functions()
{
  blend_functions.clear();
  blend_functions.resize(meshDomain.n_vertices());
  for (const auto v : meshDomain.vertices()) {
    auto& Bf = blend_functions[v.idx()];
    Bf.resize(num_loops);
    for (size_t loop = 0; loop < num_loops; ++loop) {
      Bf[loop].resize(num_sides[loop]);
      for (size_t side = 0; side < num_sides[loop]; ++side) {
        Bf[loop][side].resize(
          num_rows[loop][side], 
          std::vector<double>(num_cols[loop][side], 0.0)
        );

        const auto deg_s = ribbons[loop][side].basisU().degree();
        const auto deg_h = ribbons[loop][side].basisV().degree();

        const auto s = std::min(std::max(s_coords[v.idx()][loop][side], 0.0), 1.0);
        const auto h = std::min(std::max(h_coords[v.idx()][loop][side], 0.0), 1.0);

        const auto& Bu = ribbons[loop][side].basisU();
        const auto& Bv = Geometry::BSBasis(3, {0, 0, 0, 0, 1, 1, 1, 1});

        const size_t span_u = Bu.findSpan(s), span_v = Bv.findSpan(h);
        Geometry::DoubleVector Bh, Bs;
        Bu.basisFunctions(span_u, s, Bs);
        Bv.basisFunctions(span_v, h, Bh);

        for (size_t row = 0; row < num_rows[loop][side]; ++row) {
          for (size_t col = 0; col <= deg_s; ++col) {
            const size_t ri = row;
            const size_t ci = col + span_u - deg_s;
            const double mu = get_mu(v, loop, side, ri, ci);
            Bf[loop][side][ri][ci] =
              mu * Bs[col] * Bh[row];
          }
        }
      }
    }
  }

  return true;
}

bool SurfGBS::evaluate_mesh(bool reset)
{
  if (reset) {
    meshSurface = Mesh(meshDomain);
  }
  for (const auto v : meshDomain.vertices()) {

    OpenMesh::Vec3d pt(0.0, 0.0, 0.0);
    double sum = 0.0;
    for (size_t loop = 0; loop < num_loops; ++loop) {
      for (size_t side = 0; side < num_sides[loop]; ++side) {
        for (size_t row = 0; row < num_rows[loop][side]; ++row) {
          for (size_t col = 0; col < num_cols[loop][side]; ++col) {
            const auto Bf = blend_functions[v.idx()][loop][side][row][col];
            const auto cp = OpenMesh::Vec3d(ribbons[loop][side].controlPoint(col, row).data());
            pt += Bf * cp;
            sum += Bf;
          }
        }
      }
    }
    pt /= sum;
    meshSurface.point(meshSurface.vertex_handle(v.idx())) = pt;
  }

  return true;
}


double SurfGBS::get_mu(const Mesh::VertexHandle& vtx, size_t loop, size_t side, size_t row, size_t col) const
{
  const auto side_m1 = prev(loop, side);
  const auto side_p1 = next(loop, side);

  const auto h = h_coords[vtx.idx()][loop][side];
  const auto hm1 = h_coords[vtx.idx()][loop][side_m1];
  const auto hp1 = h_coords[vtx.idx()][loop][side_p1];

  const auto p = num_rows[loop][side];
  const auto alpha = pow(hm1, p) / (pow(hm1, p) + pow(h, p));
  const auto beta = pow(hp1, p) / (pow(hp1, p) + pow(h, p));

  double mu = 1.0;
  if (col < num_rows[loop][side_m1]) {
    mu = alpha;
  }
  if (col >= num_cols[loop][side] - num_rows[loop][side_p1]) {
    mu = beta;
  }

  return mu;

}

inline size_t circular_index(size_t i, int offset, size_t n) {
  return (i + n + (offset % static_cast<int>(n))) % n;
}

size_t SurfGBS::prev(size_t loop, size_t side) const {
  return circular_index(side, -1, num_sides[loop]);
}

size_t SurfGBS::next(size_t loop, size_t side) const {
  return circular_index(side, +1, num_sides[loop]);
}