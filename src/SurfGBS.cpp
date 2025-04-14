#include "libcdgbs/SurfGBS.hpp"

using namespace libcdgbs;

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
  //ToDo

  return true;
}

bool SurfGBS::compute_domain_mesh()
{
  //ToDo

  return true;
}

bool SurfGBS::compute_local_parameters()
{
  //ToDo

  return true;
}

bool SurfGBS::compute_blend_functions()
{
  for (const auto v : meshDomain.vertices()) {
    for (size_t loop = 0; loop < num_loops; ++loop) {
      for (size_t side = 0; side < num_sides[loop]; ++side) {
        const auto s = s_coords[v.idx()][loop][side];
        const auto h = h_coords[v.idx()][loop][side];

        Geometry::DoubleVector Bh, Bs;
        ribbons[loop][side].basisU().basisFunctions(0, s, Bs);
        ribbons[loop][side].basisV().basisFunctions(0, h, Bh);


        for (size_t row = 0; row < num_rows[loop][side]; ++row) {
          for (size_t col = 0; col < num_cols[loop][side]; ++col) {
            blend_functions[v.idx()][loop][side][row][col] =
              get_mu(v, loop, side, row, col) * Bs[col] * Bh[row];
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

    for (size_t loop = 0; loop < num_loops; ++loop) {
      for (size_t side = 0; side < num_sides[loop]; ++side) {
        for (size_t row = 0; row < num_rows[loop][side]; ++row) {
          for (size_t col = 0; col < num_cols[loop][side]; ++col) {
            const auto Bf = blend_functions[v.idx()][loop][side][row][col];
            const auto cp = OpenMesh::Vec3d(ribbons[loop][side].controlPoint(row, col).data());
            pt += Bf * cp;
          }
        }
      }
    }

    meshSurface.point(meshSurface.vertex_handle(v.idx())) = pt;
  }

  return true;
}


double SurfGBS::get_mu(const Mesh::VertexHandle& vtx, size_t loop, size_t side, size_t row, size_t col) const
{
  auto side_m1 = prev(loop, side);
  auto side_p1 = next(loop, side);

  double h = h_coords[vtx.idx()][loop][side];
  double hm1 = h_coords[vtx.idx()][loop][side_m1];
  double hp1 = h_coords[vtx.idx()][loop][side_p1];

  double mu = 1.0;
  double alpha = pow(hm1, num_rows[loop][side]) / (pow(hm1, num_rows[loop][side]) + pow(h, num_rows[loop][side]));
  double beta = pow(hp1, num_rows[loop][side]) / (pow(hp1, num_rows[loop][side]) + pow(h, num_rows[loop][side]));
  if(col < num_rows[loop][side_m1]) {
    mu = alpha;
  }
  if(col > num_cols[loop][side] - num_rows[loop][side_p1]) {
    mu = beta;
  }

  return mu;

}

inline size_t circular_index(size_t i, int offset, size_t n) {
  return (i + n + (offset % static_cast<int>(n))) % n;
}

size_t SurfGBS::prev(size_t loop, size_t side) const{
  return circular_index(side, -1, num_sides[loop]);
}

size_t SurfGBS::next(size_t loop, size_t side) const {
  return circular_index(side, +1, num_sides[loop]);
}