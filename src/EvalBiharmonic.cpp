#include "libcdgbs/Mesh.hpp"
#include "libcdgbs/SurfGBS.hpp"
#include "libcdgbs/GeomUtil.hpp"
#include "libcdgbs/MatrixUtil.hpp"

using namespace libcdgbs;
using namespace MatrixUtil;

bool SurfGBS::evaluate_mesh_biharmonic(Mesh& mesh, bool reset)
{
  if (reset) {
    mesh = Mesh(meshDomain);
  }

  size_t num_v = mesh.n_vertices();
  size_t num_b = 0;
  for(const auto& vv : mesh.vertices()) {
    if (mesh.is_boundary(vv)) {
      ++num_b;
    }
  }
  size_t num_i = num_v - num_b;


  // Computing coefficient matrix
  SparseMatrix KKT(num_v + num_i, num_v + num_i);

  SparseMatrix LL(num_v, num_v);

  buildMatrixVertexLaplace(meshDomain, LL);

  SparseMatrix L_all_bd(num_v, num_b);
  Eigen::VectorXi idx_i(num_i), idx_b(num_b);
  size_t i_i = 0, i_b = 0;
  for (const auto& vv : meshDomain.vertices()) {
    if (!meshDomain.is_boundary(vv)) {
      idx_i(i_i) = vv.idx();
      i_i++;
    }
    else {
      idx_b(i_b) = vv.idx();
      i_b++;
    }
  }

  std::vector<int> i_v(num_v);
  std::iota(i_v.begin(), i_v.begin() + num_v, 0);
  Eigen::Map<Eigen::VectorXi> idx_v(i_v.data(), i_v.size());

  SparseMatrix MM(num_v, num_v);
  SparseMatrix L_int_all(num_i, num_v);
  SparseMatrix L_int_all_T(num_v, num_i);

  buildMatrixVertexMass(meshDomain, MM);

  std::vector<int> Ri_00(num_v), Ri_01(num_v), Ri_10(num_i);
  std::vector<int> Ci_00(num_v), Ci_01(num_i), Ci_10(num_v);
  std::iota(Ri_00.begin(), Ri_00.begin() + num_v, 0);
  std::iota(Ri_01.begin(), Ri_01.begin() + num_v, 0);
  std::iota(Ri_10.begin(), Ri_10.begin() + num_i, num_v);
  std::iota(Ci_00.begin(), Ci_00.begin() + num_v, 0);
  std::iota(Ci_01.begin(), Ci_01.begin() + num_i, num_v);
  std::iota(Ci_10.begin(), Ci_10.begin() + num_v, 0);
  Eigen::Map<Eigen::VectorXi> R_00(Ri_00.data(), Ri_00.size()), R_01(Ri_01.data(), Ri_01.size()), R_10(Ri_10.data(), Ri_10.size());
  Eigen::Map<Eigen::VectorXi> C_00(Ci_00.data(), Ci_00.size()), C_01(Ci_01.data(), Ci_01.size()), C_10(Ci_10.data(), Ci_10.size());

  slice(LL, idx_i, idx_v, L_int_all);
  slice(LL, idx_v, idx_i, L_int_all_T);
  buildMatrixKKTSystem(-MM, -L_int_all, KKT);

  KKT.makeCompressed();
  
  //Computing rhs
  slice(LL, idx_v, idx_b, L_all_bd);

  Eigen::MatrixXd pos_bd = Eigen::MatrixXd::Zero(num_b, 3);
  Eigen::MatrixXd der_bd = Eigen::MatrixXd::Zero(num_v, 3);
  Eigen::MatrixXd der_scale = Eigen::MatrixXd::Zero(num_v, 1);
  for (size_t loop = 0; loop < num_loops; ++loop) {
    for (size_t side = 0; side < ribbons[loop].size(); ++side) {
      // Boundary cross derivative samples and tangent derivative samples
      const auto& cd = boundary_crossderivatives[loop][side];
      const auto& td = boundary_tangents_unnormalized[loop][side];
      std::vector<Mesh::Point> grad_h, grad_s;
      auto f_h = [this, loop, side](Mesh::VertexHandle vh) {
        const auto& hs = side_deformed[loop][side] ? h_coords_deformed : h_coords;
        return hs[vh.idx()][loop][side];
      };
      auto f_s = [this, loop, side](Mesh::VertexHandle vh) {
        return s_coords[vh.idx()][loop][side];
      };
      meshDomain.computeFaceGradientofFunction(f_h, grad_h);
      meshDomain.computeFaceGradientofFunction(f_s, grad_s);
      for (size_t ii = 0; ii < side_res[loop][side]; ++ii) {
        auto idx = domain_boundary_vertices[loop][side][ii].idx();
        auto pt = boundary_points[loop][side][ii];

        if (ii < side_res[loop][side] - 1) {
          pos_bd(idx, 0) = pt[0];
          pos_bd(idx, 1) = pt[1];
          pos_bd(idx, 2) = pt[2];
        }

        if (ii > 0 && ii < side_res[loop][side] - 1) {
          size_t idx_m1 = domain_boundary_vertices[loop][side][ii - 1].idx();
          size_t idx_p1 = domain_boundary_vertices[loop][side][ii + 1].idx();

          auto vim1 = meshDomain.vertex_handle(idx_m1);
          auto vi = meshDomain.vertex_handle(idx);
          auto vip1 = meshDomain.vertex_handle(idx_p1);
          auto pim1 = meshDomain.point(vim1);
          auto pi = meshDomain.point(vi);
          auto pip1 = meshDomain.point(vip1);
          auto tim1 = pim1 - pi;
          auto ti = pip1 - pi;
          double lim1 = tim1.norm();
          double li = ti.norm();
          double ll = (lim1 + li) / 2.0;

          Mesh::Point gh(0.0, 0.0, 0.0), gs(0.0, 0.0, 0.0);
          double area_sum = 0.0;
          for (const auto& ff : meshDomain.vf_range(vi)) {
            double area = meshDomain.calc_sector_area(meshDomain.halfedge_handle(ff));
            gh += grad_h[ff.idx()] * area;
            gs += grad_s[ff.idx()] * area;
            area_sum += area;
          }
          gh /= area_sum;
          gs /= area_sum;

          Mesh::Point tt = (ti.normalized() - tim1.normalized()).normalize_cond();
          Mesh::Point nn = (false && loop > 0 ? -1.0 : 1.0) * Mesh::Point(-tt[1], tt[0], 0.0);
          double dh = gh | nn;
          double ds = gs | nn;

          auto der = ll * (cd[ii] * dh + td[ii] * ds);
          der_bd(idx, 0) = der[0];
          der_bd(idx, 1) = der[1];
          der_bd(idx, 2) = der[2];
          der_scale(idx, 0) = ll * dh;
        }
        else if (ii == 0) {
          size_t idx_m1 = domain_boundary_vertices[loop][prev(loop ,side)][side_res[loop][side] - 2].idx();
          size_t idx_p1 = domain_boundary_vertices[loop][side][ii + 1].idx();
          auto vim1 = meshDomain.vertex_handle(idx_m1);
          auto vi = meshDomain.vertex_handle(idx);
          auto vip1 = meshDomain.vertex_handle(idx_p1);
          auto pim1 = meshDomain.point(vim1);
          auto pi = meshDomain.point(vi);
          auto pip1 = meshDomain.point(vip1);
          auto tim1 = pim1 - pi;
          auto ti = pip1 - pi;
          double lim1 = tim1.norm();
          double lip1 = ti.norm();
          double ll = (lim1 + lip1) / 2.0;

          Mesh::Point gg(0.0, 0.0, 0.0);
          double area_sum = 0.0;
          for (const auto& ff : meshDomain.vf_range(vi)) {
            double area = meshDomain.calc_sector_area(mesh.halfedge_handle(ff));
            gg += grad_h[ff.idx()] * area;
            area_sum += area;
          }
          gg /= area_sum;

          Mesh::Point tt = ti.normalize_cond();
          Mesh::Point nn = (false && loop > 0 ? -1.0 : 1.0) * Mesh::Point(-tt[1], tt[0], 0.0);
          double dh = gg | nn;

          der_scale(idx, 0) = ll * dh;
        }
        else if (ii == domain_boundary_vertices[loop][side].size() - 1) {
          size_t idx_m1 = domain_boundary_vertices[loop][side][ii - 1].idx();
          size_t idx_p1 = domain_boundary_vertices[loop][next(loop, side)][1].idx();

          auto vim1 = meshDomain.vertex_handle(idx_m1);
          auto vi = meshDomain.vertex_handle(idx);
          auto vip1 = meshDomain.vertex_handle(idx_p1);
          auto pim1 = meshDomain.point(vim1);
          auto pi = meshDomain.point(vi);
          auto pip1 = meshDomain.point(vip1);
          auto tim1 = pim1 - pi;
          auto ti = pip1 - pi;
          double lim1 = tim1.norm();
          double lip1 = ti.norm();
          double ll = (lim1 + lip1) / 2.0;

          Mesh::Point gg(0.0, 0.0, 0.0);
          double area_sum = 0.0;
          for (const auto& ff : meshDomain.vf_range(vi)) {
            double area = meshDomain.calc_sector_area(mesh.halfedge_handle(ff));
            gg += grad_h[ff.idx()] * area;
            area_sum += area;
          }
          gg /= area_sum;

          Mesh::Point tt = (-tim1).normalize_cond();
          Mesh::Point nn = (false && loop > 0 ? -1.0 : 1.0) * Mesh::Point(-tt[1], tt[0], 0.0);
          double dh = gg | nn;

          der_scale(idx, 0) = ll * dh;
        }
      }
    }
  }

  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(num_v + num_i, 3);
  rhs.block(0, 0, num_v, 3) = L_all_bd * pos_bd + der_bd;

  Eigen::MatrixXd xyz;

  auto solver_biharmonic = new Eigen::SparseLU<Eigen::SparseMatrix<double> >;
  solver_biharmonic->compute(KKT);

  if (solver_biharmonic->info() != Eigen::Success) {
    omerr() << "  !!!ERROR during factorization!!!" << std::endl;
  }

  xyz = solver_biharmonic->solve(rhs);

  if (solver_biharmonic->info() != Eigen::Success) {
    omerr() << "  !!!ERROR during backsolve!!!" << std::endl;
  }

  //solveLinearSystem(KKT, rhs, xyz);
  //omerr() << "Biharmonic system error: " << (KKT*xyz - rhs).norm()/rhs.norm() << std::endl;
  //Eigen::SparseMatrix<double> KKT_T(KKT.transpose());
  //omerr() << "KKT symmetry error: " << (KKT - KKT_T).norm() << std::endl;

  //omerr() << "KKT matrix:" << std::endl;
  //omerr() << KKT;
  //omerr() << std::endl;

  //omerr() << "rhs matrix:" << std::endl;
  //omerr() << rhs;
  //omerr() << std::endl;

  for (auto& vv : meshDomain.vertices()) {

    if (!meshDomain.is_boundary(vv)) {
      mesh.point(vv) = {
        xyz(num_v - num_b + vv.idx(), 0),
        xyz(num_v - num_b + vv.idx(), 1),
        xyz(num_v - num_b + vv.idx(), 2)
      };
    }
    else {
      mesh.point(vv) = {
        pos_bd(vv.idx(), 0),
        pos_bd(vv.idx(), 1),
        pos_bd(vv.idx(), 2)
      };
    }
  }

  mesh.request_face_normals();
  mesh.update_face_normals();
  mesh.request_vertex_normals();
  mesh.update_vertex_normals();

  return true;
}