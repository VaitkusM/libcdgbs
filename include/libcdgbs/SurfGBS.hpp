#pragma once
#include <geometry.hh>
#include "Mesh.hpp"
#include <Eigen/Eigen>

namespace libcdgbs {
  class SurfGBS {
  public:
    using Ribbon = Geometry::BSSurface;
    using Point2D = Geometry::Point2D;
    using Point3D = Geometry::Point3D;
    using VertexHandle = Mesh::VertexHandle;

    enum DomainType {
      CURVED,
      POLYGONAL,
      REGULAR
    } domain_type;

    std::vector<std::vector<Ribbon> > ribbons;

    Mesh meshDomain;
    Mesh meshSurface;

    std::vector<std::vector<std::vector<Eigen::Vector3d> > > developed_boundary_curves;
    std::vector<std::vector<std::vector<Eigen::Vector3d> > > developed_boundary_curves_normalized;
    std::vector<std::vector<std::vector<Eigen::Vector3d> > > domain_boundary_curves;
    std::vector<std::vector<std::vector<VertexHandle> > > domain_boundary_vertices;

    size_t num_loops;
    std::vector<size_t> num_sides;
    std::vector<std::vector<size_t>> num_rows;
    std::vector<std::vector<size_t>> num_cols;
    std::vector<std::vector<size_t>> deg_h;

    double target_length = 3.0;
    std::vector<std::vector<size_t> > side_res;

    std::vector<std::vector<std::vector<double> > > s_coords;
    std::vector<std::vector<std::vector<double> > > h_coords;
    std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > blend_functions;


    SurfGBS();

    void load_ribbons(const std::vector<std::vector<Ribbon> >& ribbon_surfs, double target_length = 3.0);
    void load_ribbons_and_evaluate(const std::vector<std::vector<Ribbon> >& ribbon_surfs, double target_length, Mesh& mesh);

    bool readMGBS(const std::string& filename, double target_length = 3.0);
    bool writeOBJ(const Mesh& mesh, const std::string& filename);

    bool compute_domain_boundary();
    bool compute_domain_mesh();
    bool compute_local_parameters();
    bool compute_blend_functions();
    bool evaluate_mesh(Mesh& mesh, bool reset = true);

    bool compute_harmonic_parameters();



    size_t prev(size_t loop, size_t side) const;
    size_t next(size_t loop, size_t side) const;


    // Point3D evaluate(const Point2D& uv) const;
    // Point3D evaluate(const Mesh::VertexHandle& vtx) const;

    // double get_s_coord(const Point2D& uv, size_t loop, size_t side) const;
    // double get_s_coord(const Mesh::VertexHandle& vtx, size_t loop, size_t side) const;

    // double get_h_coord(const Point2D& uv, size_t loop, size_t side) const;
    // double get_h_coord(const Mesh::VertexHandle& vtx, size_t loop, size_t side) const;

    //double get_mu(const Point2D& uv, size_t loop, size_t side) const;
    double get_mu(const VertexHandle& vtx, size_t loop, size_t side, size_t row, size_t col) const;
  };
}