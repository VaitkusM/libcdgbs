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
}