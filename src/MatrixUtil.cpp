#include "libcdgbs/MatrixUtil.hpp"

using namespace libcdgbs;

double sampleInterval(
  size_t ii,
  size_t num_samples,
  double u_min,
  double u_max
)
{
  double tt = (double)ii / (num_samples - 1);
  return (1.0 - tt) * u_min + tt * u_max;
}

void sampleInterval(
  size_t        num_samples,
  std::vector<double>& samples,
  double        u_min,
  double        u_max
)
{
  samples.clear();
  samples.resize(num_samples);
  for (size_t ii = 0; ii < num_samples; ++ii) {
    samples[ii] = sampleInterval(ii, num_samples, u_min, u_max);
  }
}

void MatrixUtil::buildMatrixVertexLaplace(
  const Mesh& mesh,
  SparseMatrix& Lap,
  bool          preallocated
)
{
  int num_v = mesh.n_vertices();
  if (Lap.rows() < num_v || Lap.cols() < num_v) {
    Lap.resize(num_v, num_v);
  }
  if (!preallocated) {
    DenseVector nnz = DenseVector::Zero(num_v);
    for (auto vv : mesh.vertices()) {
      nnz(vv.idx()) = 1 + mesh.valence(vv);
    }
    Lap.makeCompressed();
    Lap.reserve(nnz);
  }

  for (auto ff : mesh.faces()) {
    VertexHandles v;
    for (auto fv : mesh.fv_range(ff)) {
      v.push_back(fv);
    }

    VertexHandle v1 = v[0];
    VertexHandle v2 = v[1];
    VertexHandle v3 = v[2];

    auto p1 = mesh.point(v1);
    auto p2 = mesh.point(v2);
    auto p3 = mesh.point(v3);

    double A_i = ((p2 - p1) % (p3 - p1)).norm();
    double v1_cot = 0.5 * ((p2 - p1) | (p3 - p1)) / A_i;
    double v2_cot = 0.5 * ((p1 - p2) | (p3 - p2)) / A_i;
    double v3_cot = 0.5 * ((p2 - p3) | (p1 - p3)) / A_i;


    Lap.coeffRef(v1.idx(), v1.idx()) += v3_cot + v2_cot;
    Lap.coeffRef(v1.idx(), v2.idx()) += -v3_cot;
    Lap.coeffRef(v1.idx(), v3.idx()) += -v2_cot;

    Lap.coeffRef(v2.idx(), v1.idx()) += -v3_cot;
    Lap.coeffRef(v2.idx(), v2.idx()) += v3_cot + v1_cot;
    Lap.coeffRef(v2.idx(), v3.idx()) += -v1_cot;

    Lap.coeffRef(v3.idx(), v1.idx()) += -v2_cot;
    Lap.coeffRef(v3.idx(), v2.idx()) += -v1_cot;
    Lap.coeffRef(v3.idx(), v3.idx()) += v2_cot + v1_cot;
  }


}


void MatrixUtil::addConstraint2Matrix(
  const Mesh& mesh,
  const VertexHandles& fixed_vertices,
  SparseMatrix& CC,
  bool                 preallocated
)
{
  int num_vert = mesh.n_vertices();
  int num_cons = fixed_vertices.size();
  if (CC.rows() < num_cons || CC.cols() < num_vert) {
    CC.resize(num_cons, num_vert);
  }
  if (!preallocated) {
    DenseVector nnz = DenseVector::Zero(num_vert);
    for (const auto& vv : fixed_vertices) {
      nnz(vv.idx()) = 1;
    }
    CC.makeCompressed();
    CC.reserve(nnz);
  }

  size_t cons_idx = 0;
  for (const auto& vv : fixed_vertices) {
    CC.coeffRef(cons_idx, vv.idx()) = 1;
    cons_idx++;
  }
}

void MatrixUtil::addConstraint2RHS(
  const Mesh& mesh,
  const VertexHandles& fixed_vertices,
  DenseMatrix& rhs,
  size_t               col_idx,
  bool                 ramp,
  double               const_value,
  double               ramp_begin,
  double               ramp_end,
  bool                 preallocated,
  bool                 last_included
)
{
  size_t num_cons = fixed_vertices.size();
  size_t num_rows_old = rhs.rows();
  size_t num_cols_old = rhs.cols();
  if (!preallocated || num_rows_old < num_cons) {
    rhs.conservativeResize(num_rows_old + num_cons, std::max(num_cols_old, (size_t)1));
  }

  SampledCurve pts(num_cons);
  for (size_t ii = 0; ii < num_cons; ++ii) {
    auto pt = mesh.point(fixed_vertices[ii]);
    pts[ii] = { pt[0], pt[1], pt[2] };
  }


  for (size_t ii = 0; ii < num_cons; ++ii) {
    double value;
    if (ramp) {
      value = sampleInterval(
        ii,
        num_cons + (last_included ? 0 : 1),
        ramp_begin,
        ramp_end
      );
    }
    else {
      value = const_value;
    }
    size_t cons_idx = preallocated ? fixed_vertices[ii].idx() : (num_rows_old + ii);
    rhs(cons_idx, col_idx) = value;
  }
}

void MatrixUtil::buildMatrixKKTSystem(
  const SparseMatrix& AA,
  const SparseMatrix& CC,
  SparseMatrix& LL
)
{
  int num_cols_A = AA.cols();
  int num_rows_A = AA.rows();
  int num_cols_C = CC.cols();
  int num_rows_C = CC.rows();

  LL.resize(num_rows_A + num_rows_C, num_cols_A + num_rows_C);
  DenseVector num_nnz = DenseVector(LL.cols());

  for (int kk = 0; kk < AA.outerSize(); ++kk) {
    for (SparseMatrix::InnerIterator it(AA, kk); it; ++it) {
      num_nnz(it.col())++;
    }
  }

  for (int kk = 0; kk < CC.outerSize(); ++kk) {
    for (SparseMatrix::InnerIterator it(CC, kk); it; ++it) {
      num_nnz(it.col())++;
      num_nnz(num_cols_A + it.row())++;
    }
  }

  LL.reserve(num_nnz);

  for (int kk = 0; kk < AA.outerSize(); ++kk) {
    for (SparseMatrix::InnerIterator it(AA, kk); it; ++it) {
      LL.coeffRef(it.row(), it.col()) = it.value();
    }
  }

  for (int kk = 0; kk < CC.outerSize(); ++kk) {
    for (SparseMatrix::InnerIterator it(CC, kk); it; ++it) {
      LL.coeffRef(num_rows_A + it.row(), it.col()) = it.value();
      LL.coeffRef(it.col(), num_rows_A + it.row()) = it.value();
    }
  }

  LL.makeCompressed();
}

void MatrixUtil::buildMatrixKKTRHS(
  const DenseMatrix& bb,
  const DenseMatrix& dd,
  DenseMatrix& ff
)
{
  ff.resize(bb.rows() + dd.rows(), bb.cols());
  ff.topRows(bb.rows()) = bb;
  ff.bottomRows(dd.rows()) = dd;
}

bool MatrixUtil::solveLinearSystem(
  const SparseMatrix& AA,
  const DenseMatrix& bb,
  DenseMatrix& xx
)
{
  Eigen::SparseLU<SparseMatrix> solver;
  std::cout << "Factorizing KKT matrix...";
  solver.compute(AA);
  std::cout << " DONE!" << std::endl;

  if (solver.info() != Eigen::Success) {
    std::cout << "  !!!ERROR during factorization!!!" << std::endl;
    return false;
  }

  std::cout << "Backsolving for KKT system...";
  xx = solver.solve(bb);
  std::cout << " DONE!" << std::endl;

  if (solver.info() != Eigen::Success) {
    std::cout << "  !!!ERROR during backsolve!!!" << std::endl;
    return false;
  }
  return true;
}