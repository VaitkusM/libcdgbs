#include "libcdgbs/SurfGBS.hpp"
#include "libcdgbs/MatrixUtil.hpp"

using namespace libcdgbs;
using namespace MatrixUtil;

template <typename T>
void concatenateVectors(
  const std::vector<T>& v1,
  const std::vector<T>& v2,
  std::vector<T>& vc,
  bool                  common_joint = false
)
{
  const size_t s1 = v1.size();
  const size_t s2 = v2.size();
  const size_t shift = (common_joint ? 1 : 0);

  if (s1 == 0 && s2 == 0) {
    return;
  }

  vc.clear();
  vc.resize(s1 + s2 - shift);

  for (size_t ii = 0; ii < s1; ++ii) {
    vc[ii] = v1[ii];
  }
  for (size_t ii = 0; ii < s2; ++ii) {
    if (ii >= shift) {
      vc[s1 + ii - shift] = v2[ii];
    }
  }
}

bool SurfGBS::compute_harmonic_parameters()
{
  //Initialize
  s_coords.clear();
  s_coords.resize(meshDomain.n_vertices());
  h_coords.clear();
  h_coords.resize(meshDomain.n_vertices());
  for (const auto v : meshDomain.vertices()) {
    auto& s = s_coords[v.idx()];
    auto& h = h_coords[v.idx()];
    s.resize(num_loops);
    h.resize(num_loops);
    for (size_t loop = 0; loop < num_loops; ++loop) {
      s[loop].resize(num_sides[loop]);
      h[loop].resize(num_sides[loop]);
    }
  }

  auto& mesh = meshDomain;
  size_t num_vert = mesh.n_vertices();

  SparseMatrix QQ(num_vert, num_vert);

  buildMatrixVertexLaplace(mesh, QQ);

  size_t num_sides_all = std::accumulate(num_sides.begin(), num_sides.end(), 0);

  { // Compute h-coordinates
    std::vector<VertexHandle> sides_pts;
    for (size_t li = 0; li < num_loops; ++li) {
      for (size_t si = 0; si < num_sides[li]; ++si) {
        for (size_t i = 0; i < domain_boundary_vertices[li][si].size() - 1; ++i) {
          auto vtx = domain_boundary_vertices[li][si][i];
          sides_pts.push_back(vtx);
        }
      }
    }

    size_t num_cons = sides_pts.size();
    SparseMatrix Ch(num_cons, num_vert), KKT;
    DenseMatrix pp = DenseMatrix::Zero(num_vert, num_sides_all);
    DenseMatrix dh = DenseMatrix::Ones(num_cons, num_sides_all);
    DenseMatrix rhs, x;

    addConstraint2Matrix(mesh, sides_pts, Ch);

    size_t idx = 0;
    for (size_t loop = 0; loop < num_loops; ++loop) {
      for (size_t side = 0; side < num_sides[loop]; ++side) {
        size_t side_m1 = prev(loop, side);
        size_t side_p1 = next(loop, side);
        std::vector<VertexHandle> side_pts;
        std::vector<VertexHandle> side_pts_m1;
        std::vector<VertexHandle> side_pts_p1;
        for (size_t i = 0; i < domain_boundary_vertices[loop][side].size(); ++i) {
          auto vtx = domain_boundary_vertices[loop][side][i];
          side_pts.push_back(vtx);
        }
        for (size_t i = 0; i < domain_boundary_vertices[loop][side_m1].size(); ++i) {
          auto vtx = domain_boundary_vertices[loop][side_m1][i];
          side_pts_m1.push_back(vtx);
        }
        for (size_t i = 0; i < domain_boundary_vertices[loop][side_p1].size(); ++i) {
          auto vtx = domain_boundary_vertices[loop][side_p1][i];
          side_pts_p1.push_back(vtx);
        }

        addConstraint2RHS(
          mesh,
          side_pts,
          dh,
          idx,
          false,
          0.0,
          0.0,
          0.0,
          true,
          true
        );

        addConstraint2RHS(
          mesh,
          side_pts_m1,
          dh,
          idx,
          true,
          0.0,
          1.0,
          0.0,
          true,
          true
        );

        addConstraint2RHS(
          mesh,
          side_pts_p1,
          dh,
          idx,
          true,
          0.0,
          0.0,
          1.0,
          true,
          true
        );

        ++idx;
      }
    }

    buildMatrixKKTSystem(QQ, Ch, KKT);
    KKT.makeCompressed();
    Ch.resize(0, 0);
    //CC.data().squeeze();
    buildMatrixKKTRHS(pp, dh, rhs);
    dh.resize(0, 0);

    if (!solveLinearSystem(KKT, rhs, x)) {
      return false;
    }

    for (auto v : mesh.vertices()) {
      size_t idx = 0;
      for (size_t loop = 0; loop < num_loops; ++loop) {
        for (size_t side = 0; side < num_sides[loop]; ++side) {
          h_coords[v.idx()][loop][side] = x(v.idx(), idx);
          ++idx;
        }
      }
    }

  }

  //Compute s-coordinates
  for (size_t loop = 0; loop < num_loops; ++loop) {
    for (size_t side = 0; side < num_sides[loop]; ++side) {
      // std::cout << "Computing s-coordinates for loop " << loop << " side " << side << std::endl;
      size_t side_m1 = prev(loop, side);
      size_t side_p1 = next(loop, side);
      std::vector<VertexHandle> side_pts;
      std::vector<VertexHandle> side_pts_m1;
      std::vector<VertexHandle> side_pts_p1;
      for (size_t i = 0; i < domain_boundary_vertices[loop][side].size() - 1; ++i) {
        auto vtx = domain_boundary_vertices[loop][side][i];
        side_pts.push_back(vtx);
      }
      for (size_t i = 0; i < domain_boundary_vertices[loop][side_m1].size() - 1; ++i) {
        auto vtx = domain_boundary_vertices[loop][side_m1][i];
        side_pts_m1.push_back(vtx);
      }
      for (size_t i = 0; i < domain_boundary_vertices[loop][side_p1].size(); ++i) {
        auto vtx = domain_boundary_vertices[loop][side_p1][i];
        side_pts_p1.push_back(vtx);
      }

      std::vector<VertexHandle> sides_pts;
      concatenateVectors(
        side_pts_m1,
        side_pts,
        sides_pts,
        false
      );
      concatenateVectors(
        std::vector<VertexHandle>(sides_pts),
        side_pts_p1,
        sides_pts,
        false
      );

      // std::cout << "sides_pts size: " << sides_pts.size() << std::endl;

      size_t num_cons = sides_pts.size();
      SparseMatrix CC(num_cons, num_vert), KKT;
      DenseMatrix pp = DenseMatrix::Zero(num_vert, 1);
      DenseMatrix dd;
      DenseMatrix rhs, x;

      addConstraint2Matrix(mesh, sides_pts, CC);

      // std::cout << "CC size: " << CC.rows() << " " << CC.cols() << std::endl;
      // //Printing non-zeroes of CC
      // for (int k = 0; k < CC.outerSize(); ++k) {
      //   for (SparseMatrix::InnerIterator it(CC, k); it; ++it) {
      //     std::cout << "CC[" << it.row() << "][" << it.col() << "] = " << it.value() << std::endl;
      //   }
      // }

      addConstraint2RHS(
        mesh,
        side_pts_m1,
        dd,
        0,
        false,
        0.0,
        0.0,
        0.0,
        false,
        false
      );

      addConstraint2RHS(
        mesh,
        side_pts,
        dd,
        0,
        true,
        0.0,
        0.0,
        1.0,
        false,
        false
      );

      addConstraint2RHS(
        mesh,
        side_pts_p1,
        dd,
        0,
        false,
        1.0,
        0.0,
        1.0,
        false,
        true
      );

      // //printing elements of dd
      // for (int k = 0; k < dd.rows(); ++k) {
      //   std::cout << "dd[" << k << "] = " << dd(k, 0) << std::endl;
      //   // draw dashes after every 100th element
      //   if (k % 100 == 99) {
      //     std::cout << "------------------------" << std::endl;
      //   }
      // }

      buildMatrixKKTSystem(QQ, CC, KKT);
      KKT.makeCompressed();
      CC.resize(0, 0);
      //CC.data().squeeze();
      buildMatrixKKTRHS(pp, dd, rhs);
      dd.resize(0, 0);

      if (!solveLinearSystem(KKT, rhs, x)) {
        return false;
      }

      for (auto v : mesh.vertices()) {
        s_coords[v.idx()][loop][side] = x(v.idx(), 0);
      }

    }
  }


  return true;
}