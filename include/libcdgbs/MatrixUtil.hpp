#pragma once

#include <Eigen/Eigen>
#include "libcdgbs/Mesh.hpp"

namespace MatrixUtil {
  using DenseMatrix = Eigen::MatrixXd;
  using DenseVector = Eigen::VectorXd;
  using SparseMatrix = Eigen::SparseMatrix<double>;
  using VertexHandle = libcdgbs::Mesh::VertexHandle;
  using VertexHandles = std::vector<VertexHandle>;
  using SampledCurve = std::vector<Eigen::Vector3d>;
  using Mesh = libcdgbs::Mesh;

  void buildMatrixVertexLaplace(
    const Mesh& mesh,
    SparseMatrix& Lap,
    bool preallocated = false
  );

  void buildMatrixVertexMass(
    const Mesh& mesh,
    SparseMatrix& Mass,
    bool preallocated = false
  );

  void addConstraint2Matrix(
    const Mesh&        mesh,
    const VertexHandles& fixed_vertices,
    SparseMatrix&        CC,
    bool                 preallocated   = false
  );

  void addConstraint2RHS(
    const Mesh&        mesh,
    const VertexHandles& fixed_vertices,
    DenseMatrix&     rhs,
    size_t               col_idx        = 0,
    bool                 ramp           = false,
    double               const_value    = 0.0, 
    double               ramp_begin     = 0.0,
    double               ramp_end       = 1.0,
    bool                 preallocated   = false,
    bool                 last_included  = true
  );

    void addConstraint2RHS(
    const Mesh&                mesh,
    const VertexHandles&       fixed_vertices,
    const std::vector<double>& values,
    DenseMatrix&               rhs,
    size_t                     col_idx        = 0,
    bool                       reversed       = false,
    bool                       preallocated   = false,
    bool                       last_included  = true
  );

  

  void buildMatrixKKTSystem(
    const SparseMatrix& QQ,
    const SparseMatrix& CC,
    SparseMatrix& KKT
  );

  void buildMatrixKKTRHS(
    const DenseMatrix& pp,
    const DenseMatrix& dd,
    DenseMatrix& rhs
  );

  bool solveLinearSystem(
    const SparseMatrix& AA,
    const DenseMatrix& bb,
    DenseMatrix& xx
  );

    // Act like the matlab X(row_indices,col_indices) operator, where
  // row_indices, col_indices are non-negative integer indices.
  //
  // Inputs:
  //   X  m by n matrix
  //   R  list of row indices
  //   C  list of column indices
  // Output:
  //   Y  #R by #C matrix
  void slice(
    const Eigen::SparseMatrix<double> &X,
    const Eigen::VectorXi             &R,
    const Eigen::VectorXi             &C,
    Eigen::SparseMatrix<double>       &Y 
  );

  // Act like the matlab Y(row_indices,col_indices) = X
  // 
  // Inputs:
  //   X  xm by xn rhs matrix
  //   R  list of row indices
  //   C  list of column indices
  //   Y  ym by yn lhs matrix
  // Output:
  //   Y  ym by yn lhs matrix, same as input but Y(R,C) = X
  void slice_into(
    const Eigen::SparseMatrix<double> &X,
    const Eigen::VectorXi             &R,
    const Eigen::VectorXi             &C,
    Eigen::SparseMatrix<double>       &Y
  );
}