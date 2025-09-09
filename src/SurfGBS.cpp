#include <fstream>

#include "libcdgbs/SurfGBS.hpp"
#include "libcdgbs/LoopFlattener.hpp"
#include "libcdgbs/TriangleWrapper.hpp"
#include "libcdgbs/GeomUtil.hpp"

using namespace libcdgbs;
using Vec3 = Eigen::Vector3d;
using Vec2 = Eigen::Vector2d;
using SubCurve3D = std::vector<Vec3>;
using Curve3D = std::vector<SubCurve3D>;
using Curves3D = std::vector<Curve3D>;

//calculate arclength of each ribbon boundary B-spline curve
double getLength(const Geometry::BSSurface& rib, double u_start = 0.0, double u_end = 1.0) {
  double length = 0.0;
  for (size_t i = 0; i < 200; ++i) {
    double ua = (1.0 - (double(i) / 200.0)) * u_start + (double(i) / 200.0) * u_end;
    double ub = (1.0 - (double(i + 1) / 200.0)) * u_start + (double(i + 1) / 200.0) * u_end;;
    auto pa = rib.eval(ua, 0.0);
    auto pb = rib.eval(ub, 0.0);
    length += (pb - pa).norm();
  }
  return length;
}


SurfGBS::SurfGBS()
{

}

void SurfGBS::load_ribbons(const std::vector<std::vector<Ribbon> >& ribbon_surfs, double target_length, bool merge_corners, double deform_value, bool restrict_params, bool c1_merge)
{
  ribbons = ribbon_surfs;

  SurfGBS::target_length = target_length;
  SurfGBS::deform_value = deform_value;
  SurfGBS::restrict_params = restrict_params;
  SurfGBS::c1_merge = c1_merge;

  num_segments.clear();
  num_segments.resize(ribbons.size());
  for (size_t loop = 0; loop < ribbons.size(); ++loop) {
    num_segments[loop].resize(ribbons[loop].size(), 1);
  }
  init_data();

  if (merge_corners) {
    merge_smooth_corners();
  }
}

[[maybe_unused]]
static void writeLoops(std::vector<std::vector<std::vector<Eigen::Vector3d> > > loops,
  std::string filename = "/tmp/boundary.obj") {
  size_t index = 1;
  std::ofstream f(filename);
  for (const auto& loop : loops)
    for (const auto& curve : loop) {
      size_t start = index;
      for (const auto& p : curve) {
        f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
        index++;
      }
      f << 'l';
      for (size_t i = start; i < index; ++i)
        f << ' ' << i;
      f << std::endl;
    }
}

void SurfGBS::load_ribbons_and_evaluate(const std::vector<std::vector<Ribbon> >& ribbon_surfs, double target_length, Mesh& mesh, bool merge_smooth_corners, double deform_value, bool restrict_params, bool c1_merge)
{
  load_ribbons(ribbon_surfs, target_length, merge_smooth_corners, deform_value, restrict_params, c1_merge);
  compute_domain_boundary();
  // writeLoops(domain_boundary_curves, "/tmp/boundary.obj");
  compute_domain_mesh();
  compute_local_parameters();
  compute_blend_functions();
  evaluate_mesh(mesh);
}


void SurfGBS::init_data()
{
  num_sides.clear();
  num_rows.clear();
  num_cols.clear();
  side_res.clear();
  side_segment_res.clear();
  deg_h.clear();
  deform_splines.clear();
  side_deformed.clear();

  num_loops = ribbons.size();
  num_sides.resize(num_loops);
  num_rows.resize(num_loops);
  num_cols.resize(num_loops);
  side_res.resize(num_loops);
  side_segment_res.resize(num_loops);
  deg_h.resize(num_loops);
  deform_splines.resize(num_loops);
  side_deformed.resize(num_loops);
  for (size_t loop = 0; loop < num_loops; ++loop) {
    const size_t ns = ribbons[loop].size();
    num_sides[loop] = ns;
    deg_h[loop].resize(ns, 3);
    num_rows[loop].resize(ns);
    num_cols[loop].resize(ns);
    side_res[loop].resize(ns);
    side_segment_res[loop].resize(ns);
    for (size_t side = 0; side < ns; ++side) {
      //Nmber of B-spline control points in U and V directions
      num_rows[loop][side] = ribbons[loop][side].basisV().knots().size() - ribbons[loop][side].basisV().degree() - 1;
      num_cols[loop][side] = ribbons[loop][side].basisU().knots().size() - ribbons[loop][side].basisU().degree() - 1;
      side_res[loop][side] = 100; // default resolution
      side_segment_res[loop][side].resize(num_segments[loop][side], 100); // default resolution per segment
    }
  }

  // Set number of samples based on target edge length
  for (size_t loop = 0; loop < num_loops; ++loop) {
    for (size_t side = 0; side < num_sides[loop]; ++side) {
      auto rib = ribbons[loop][side];
      if(num_segments[loop][side] == 1) {
        auto curve_length = getLength(rib);

        auto num_samples = std::max(1.0, curve_length / target_length);
        side_res[loop][side] = std::max(static_cast<size_t>(num_samples), size_t(5));
        //side_res[loop][side] = 200; //forcing 200 for now
        side_segment_res[loop][side][0] = side_res[loop][side];
      }
      else {
        side_res[loop][side] = 0;
        for(size_t seg = 0; seg < num_segments[loop][side]; ++seg) {
          double u_start = double(seg) / (num_segments[loop][side]);
          double u_end = double(seg + 1) / (num_segments[loop][side]);
          auto curve_length = getLength(rib, u_start, u_end);

          auto num_samples = std::max(1.0, curve_length / target_length);
          const size_t seg_res = std::max(static_cast<size_t>(num_samples), size_t(5));
          side_res[loop][side] += seg_res;
          side_segment_res[loop][side][seg]= seg_res;
          if(seg > 0) {
            side_res[loop][side] -= 1; //remove duplicate point at segment boundary
          }
        }
      }
    }
  }

  // Initialize deformation splines (bi-quadratic B-splines with 3x3 control points)
  for (size_t loop = 0; loop < num_loops; ++loop) {
    for (size_t side = 0; side < num_sides[loop]; ++side) {
      const size_t deg_u = 2;
      const size_t deg_v = 2;
      const Geometry::DoubleVector knots {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};

      Geometry::VectorVector cps;
      cps.push_back({ 0.0, 0.0, 0.0 });
      cps.push_back({ 0.5, 0.5, 0.5 });
      cps.push_back({ 1.0, 1.0, 1.0 });
      cps.push_back({ 0.0, 0.0, 0.0 });
      cps.push_back({ 0.5, 0.5, 0.5 });
      cps.push_back({ 1.0, 1.0, 1.0 });
      cps.push_back({ 0.0, 0.0, 0.0 });
      cps.push_back({ 0.5, 0.5, 0.5 });
      cps.push_back({ 1.0, 1.0, 1.0 });
      
      deform_splines[loop].push_back(Geometry::BSSurface(deg_u, deg_v, knots, knots, cps));
      side_deformed[loop].push_back(false);
    }
  }


}

bool SurfGBS::compute_domain_boundary()
{
  developed_boundary_curves.clear();
  developed_boundary_curves_normalized.clear();
  domain_boundary_curves.clear();
  developed_boundary_curves.resize(num_loops);
  developed_boundary_curves_normalized.resize(num_loops);
  domain_boundary_curves.resize(num_loops);
  domain_boundary_params.clear();
  domain_boundary_params.resize(num_loops);

  Curves3D points(num_loops);
  Curves3D normals(num_loops);
  for (size_t loop = 0; loop < num_loops; ++loop) {
    points[loop].resize(num_sides[loop]);
    normals[loop].resize(num_sides[loop]);
    domain_boundary_params[loop].resize(num_sides[loop]);
    for (size_t side = 0; side < num_sides[loop]; ++side) {
      const auto& rib = ribbons[loop][side];
      const size_t res = side_res[loop][side];
      points[loop][side].reserve(res);
      normals[loop][side].reserve(res);
      domain_boundary_params[loop][side].reserve(res);
      if (num_segments[loop][side] > 1) {
        for (size_t i = 0; i < num_segments[loop][side]; ++i) {
          const double u_start = double(i) / (num_segments[loop][side]);
          const double u_end = double(i + 1) / (num_segments[loop][side]);
          const size_t seg_res = side_segment_res[loop][side][i];
          for(size_t j = 0; j < seg_res; ++j) {
            if(i > 0 && j == 0) {
              continue; //skip duplicate point at segment boundary
            }
            double u = (1.0 - (double(j) / (seg_res - 1))) * u_start + (double(j) / (seg_res - 1)) * u_end;
            Geometry::VectorMatrix duv;
            auto pt = rib.eval(u, 0.0, 1, duv);

            points[loop][side].push_back({ pt[0], pt[1], pt[2] });
            auto du = duv[1][0];
            auto dv = duv[0][1];
            auto nn = (du ^ dv).normalized();
            normals[loop][side].push_back({ nn[0], nn[1], nn[2] });
            domain_boundary_params[loop][side].push_back(u);
          }
        }
      }
      else {
        for (size_t i = 0; i < res; ++i) {
          double u = double(i) / (res - 1);
          Geometry::VectorMatrix duv;
          auto pt = rib.eval(u, 0.0, 1, duv);

          points[loop][side].push_back({ pt[0], pt[1], pt[2] });
          auto du = duv[1][0];
          auto dv = duv[0][1];
          auto nn = (du ^ dv).normalized();
          normals[loop][side].push_back({ nn[0], nn[1], nn[2] });
          domain_boundary_params[loop][side].push_back(u);
        }
      }
    }

    // Flip vector ?
    auto va = GeomUtil::vectorArea(points[loop]);
    Eigen::Vector3d avg_n(0.0, 0.0, 0.0);
    size_t num_pts = 0;
    for (size_t side = 0; side < normals[loop].size(); ++side) {
      for (size_t ii = 0; ii < normals[loop][side].size(); ++ii) {
        avg_n += normals[loop][side][ii];
        num_pts++;
      }
    }
    avg_n /= num_pts;
    bool flip = loop > 0 && (va.dot(avg_n)) < 0;
    if(flip) {
      std::cout << "Flipping loop " << loop << " with vector area " << va.transpose() << " and average normal " << avg_n.transpose() << std::endl;
    }

    developed_boundary_curves[loop] =
      LoopFlattener::developLoop(points[loop], normals[loop], flip);

    developed_boundary_curves_normalized[loop] =
      LoopFlattener::normalizeLoopAngles(points[loop], normals[loop], flip);

    domain_boundary_curves[loop] =
      LoopFlattener::closeLoop(developed_boundary_curves_normalized[loop]);

    auto cog = GeomUtil::centroid(domain_boundary_curves[loop]);
    GeomUtil::shiftByVector(domain_boundary_curves[loop], -cog);
    // GeomUtil::shiftByVector(developed_boundary_curves_normalized[loop], -cog);
    // GeomUtil::shiftByVector(developed_boundary_curves[loop], -cog);

    auto aligned = domain_boundary_curves[loop];
    GeomUtil::alignCurveLoopPCA(domain_boundary_curves[loop], aligned, false);
    domain_boundary_curves[loop] = aligned;
  }

  if (debug_outputs) {
    writeLoops(domain_boundary_curves, "boundary_uv_0.obj");
  }

  if (num_loops > 1) { // Handling inner loops
    //Computing perimeter surface
    std::vector<std::vector<Ribbon> > perimeter_ribbons(1, ribbons.front());
    SurfGBS perimeter_gbs;
    perimeter_gbs.load_ribbons(perimeter_ribbons, target_length);
    perimeter_gbs.compute_domain_boundary();
    perimeter_gbs.compute_domain_mesh();
    perimeter_gbs.compute_local_parameters();
    perimeter_gbs.compute_blend_functions();
    perimeter_gbs.evaluate_mesh(perimeter_gbs.meshSurface, true);

    if (debug_outputs) {
      writeOBJ(perimeter_gbs.meshSurface, "perimeter.obj");
    }

    //Projecting inner loops
    for (size_t loop = 1; loop < num_loops; ++loop) {
      auto curv_proj = domain_boundary_curves[loop];
      const auto& curv = points[loop];
      perimeter_gbs.projectCurves2Domain(curv, curv_proj);
      if (debug_outputs) {
        writeLoops({ curv_proj }, std::string("boundary_proj_") + std::to_string(loop) + ".obj");
      }
      auto cog = GeomUtil::centroid(curv_proj);
      GeomUtil::shiftByVector(curv_proj, -cog);
      GeomUtil::alignPointSets(domain_boundary_curves[loop], curv_proj);
      GeomUtil::shiftByVector(domain_boundary_curves[loop], cog);
    }
  }

  if (debug_outputs) {
    writeLoops(domain_boundary_curves, "boundary_uv_1.obj");
    writeLoops(developed_boundary_curves, "boundary_developed.obj");
    writeLoops(developed_boundary_curves_normalized, "boundary_normalized.obj");
    writeLoops(points, "boundary_xyz.obj");
  }

  return true;
}

bool SurfGBS::compute_domain_mesh()
{
  auto triangle_wrapper = TriangleWrapper();
  meshDomain = triangle_wrapper.triangulate_loop(domain_boundary_curves, target_length);

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
  if(!compute_deformed_parameters()) {
    std::cout << "Error computing deformed parameters" << std::endl;
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
        const auto deg_h = SurfGBS::deg_h[loop][side];

        const auto &hs = side_deformed[loop][side] ? h_coords_deformed : h_coords;

        const auto s = std::min(std::max(s_coords[v.idx()][loop][side], 0.0), 1.0);
        const auto h = std::min(std::max(hs[v.idx()][loop][side], 0.0), 1.0);

        const auto& Bu = ribbons[loop][side].basisU();
        const auto& Bv = Geometry::BSBasis(deg_h, [deg_h]() {
          Geometry::DoubleVector knots(deg_h + 1, 0.0);
          knots.insert(knots.end(), deg_h + 1, 1.0);
          return knots;
          }());

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

bool SurfGBS::evaluate_mesh(Mesh& mesh, bool reset)
{
  if (reset) {
    mesh = Mesh(meshDomain);
  }
  for (const auto v : mesh.vertices()) {

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
    mesh.point(mesh.vertex_handle(v.idx())) = pt;
  }

  mesh.request_face_normals();
  mesh.update_face_normals();
  mesh.request_vertex_normals();
  mesh.update_vertex_normals();

  return true;
}


double SurfGBS::get_mu(const Mesh::VertexHandle& vtx, size_t loop, size_t side, size_t row, size_t col) const
{
  const auto eps = 1e-10;

  const auto side_m1 = prev(loop, side);
  const auto side_p1 = next(loop, side);

  const auto& hs = side_deformed[loop][side] ? h_coords_deformed : h_coords;
  const auto& hs_m1 = side_deformed[loop][side_m1] ? h_coords_deformed : h_coords;
  const auto& hs_p1 = side_deformed[loop][side_p1] ? h_coords_deformed : h_coords;

  const auto h = hs[vtx.idx()][loop][side];
  const auto hm1 = hs_m1[vtx.idx()][loop][side_m1];
  const auto hp1 = hs_p1[vtx.idx()][loop][side_p1];

  const auto p = num_rows[loop][side];
  const auto alpha = pow(hm1, p) / (pow(hm1, p) + pow(h, p));
  const auto beta = pow(hp1, p) / (pow(hp1, p) + pow(h, p));

  double mu = 1.0;
  if ((hm1 < eps && h < eps) || (h < eps && hp1 < eps)) {
    mu = 0.5;
  }
  else {
    if (col < num_rows[loop][side_m1]) {
      mu = alpha;
    }
    if (col >= num_cols[loop][side] - num_rows[loop][side_p1]) {
      mu = beta;
    }
  }

  return mu;

}

void SurfGBS::projectCurves2Domain(
  const std::vector<std::vector<Eigen::Vector3d>>& curves_xyz,
  std::vector<std::vector<Eigen::Vector3d>>& curves_uv
) const {
  if (meshDomain.n_faces() != meshSurface.n_faces() || meshDomain.n_faces() == 0) {
    return;
  }
  curves_uv.clear();
  curves_uv.resize(curves_xyz.size());
  for (size_t side = 0; side < curves_xyz.size(); ++side) {
    size_t num_pts = curves_xyz[side].size();
    curves_uv[side].resize(num_pts);
    for (size_t ii = 0; ii < num_pts; ++ii) {
      Eigen::Vector3d pt = curves_xyz[side][ii];
      auto closest_f = meshSurface.findClosestFace({ pt[0], pt[1], pt[2] });
      curves_uv[side][ii] = project2Triangle_uv(pt, closest_f);
    }
  }
}

Eigen::Vector3d SurfGBS::project2Triangle_uv(Eigen::Vector3d pt, Mesh::FaceHandle ff) const
{
  auto f = meshSurface.make_smart(ff);
  auto v1 = f.halfedge().to();
  auto v2 = f.halfedge().next().to();
  auto v3 = f.halfedge().next().next().to();
  auto p1_xyz = meshSurface.point(v1);
  auto p2_xyz = meshSurface.point(v2);
  auto p3_xyz = meshSurface.point(v3);
  auto p1_uv = meshDomain.point(v1);
  auto p2_uv = meshDomain.point(v2);
  auto p3_uv = meshDomain.point(v3);

  // Project point to the plane of the triangle
  auto pp = Mesh::Point(pt[0], pt[1], pt[2]) - (meshSurface.normal(ff) | (Mesh::Point(pt[0], pt[1], pt[2]) - p1_xyz)) * meshSurface.normal(ff);

  // Barycentric coordinates
  double u, v, w;
  Mesh::barycentricCoordinates(pp, p1_xyz, p2_xyz, p3_xyz, u, v, w);

  auto pt_uv = u * p1_uv + v * p2_uv + w * p3_uv;
  return { pt_uv[0], pt_uv[1], 0.0 };
}

void SurfGBS::merge_smooth_corners() {
  auto new_ribbons = ribbons;
  for(size_t loop = 0; loop < 1; ++loop) { // Skip inner loops
    new_ribbons[loop].clear();
    new_ribbons[loop].push_back(ribbons[loop].front());
    std::vector<size_t> num_merged_sides(1, 1);
    for(size_t side = 1; side < num_sides[loop]; ++side) {
      const auto side_m1 = prev(loop, side);
      // const auto side_p1 = next(loop, side);
      const auto rib1 = new_ribbons[loop].back();
      const auto rib2 = ribbons[loop][side];
      // Check angle between sides
      Geometry::VectorMatrix duv1, duv2;
      rib1.eval(1.0, 0.0, 1.0, duv1);
      rib2.eval(0.0, 0.0, 1.0, duv2);
      const double angle = std::acos(std::clamp(duv1[1][0].normalized() * duv2[1][0].normalized(), -1.0, 1.0)) * 180.0 / M_PI;
      // std::cout << "Loop " << loop << " sides " << side_m1 << " and " << side << " angle: "
      //  << angle << std::endl;
      if (angle < 1e-3) {
        // Create new ribbon
        const size_t deg_u1 = rib1.basisU().degree();
        const size_t deg_u2 = rib2.basisU().degree();
        const size_t deg_v = rib1.basisV().degree();
        if(deg_u1 != deg_u2) {
          std::cerr << "Error: Cannot merge ribbons with different degrees!" << std::endl;
          continue;
        }
        else {
          std::cout << "Merging ribbons at loop " << loop << " sides " << side_m1 << " and " << side << std::endl;
        }
        size_t num_segments = num_merged_sides.back() + 1;
        Geometry::DoubleVector knots_u, knots_v; // ToDo: Merging non-uniform knotvectors
        knots_u.insert(knots_u.end(), deg_u1 + 1, 0.0);
        for (size_t i = 1; i < num_segments; ++i) {
          knots_u.insert(knots_u.end(), deg_u1 - (c1_merge ? 1 : 0), i / double(num_segments));
        }
        knots_u.insert(knots_u.end(), deg_u1 + 1, 1.0);
        knots_v = rib1.basisV().knots();
        Geometry::PointVector cpts;
        for(size_t c = 0; c < rib1.numControlPoints().at(0) - (c1_merge ? 1 : 0); ++c) {
          for(size_t r = 0; r < rib1.numControlPoints().at(1); ++r) {
            cpts.push_back(rib1.controlPoint(c, r));
          }
        }
        for (size_t c = 1; c < rib2.numControlPoints().at(0); ++c) {
          for (size_t r = 0; r < rib2.numControlPoints().at(1); ++r) {
            cpts.push_back(rib2.controlPoint(c, r));
          }
        }
        Geometry::BSSurface new_rib(deg_u1, deg_v, knots_u, knots_v, cpts);
        new_ribbons[loop].back() = new_rib;
        ++num_merged_sides.back();
      }
      else {
        new_ribbons[loop].push_back(ribbons[loop][side]);
        num_merged_sides.push_back(1);
      }
    }

    {// Merge last and first ribbons if needed
      size_t side = 0;
      const auto side_m1 = prev(loop, side);
      const auto rib1 = new_ribbons[loop].back();
      const auto rib2 = new_ribbons[loop].front();
      // Check angle between sides
      Geometry::VectorMatrix duv1, duv2;
      rib1.eval(1.0, 0.0, 1.0, duv1);
      rib2.eval(0.0, 0.0, 1.0, duv2);
      const double angle = std::acos(std::clamp(duv1[1][0].normalized() * duv2[1][0].normalized(), -1.0, 1.0)) * 180.0 / M_PI;
      // std::cout << "Loop " << loop << " sides " << side_m1 << " and " << side << " angle: "
      //  << angle << std::endl;
      if (angle < 1e-3) {
        // Create new ribbon
        const size_t deg_u1 = rib1.basisU().degree();
        const size_t deg_u2 = rib2.basisU().degree();
        const size_t deg_v = rib1.basisV().degree();
        if (deg_u1 != deg_u2) {
          std::cerr << "Error: Cannot merge ribbons with different degrees!" << std::endl;
          continue;
        }
        else {
          std::cout << "Merging ribbons at loop " << loop << " sides " << side_m1 << " and " << side << std::endl;
        }
        size_t num_segments = num_merged_sides.back() + num_merged_sides.front();
        Geometry::DoubleVector knots_u, knots_v; // ToDo: Merging non-uniform knotvectors
        knots_u.insert(knots_u.end(), deg_u1 + 1, 0.0);
        for (size_t i = 1; i < num_segments; ++i) {
          knots_u.insert(knots_u.end(), deg_u1 - (c1_merge ? 1 : 0), i / double(num_segments));
        }
        knots_u.insert(knots_u.end(), deg_u1 + 1, 1.0);
        knots_v = rib1.basisV().knots();
        Geometry::PointVector cpts;
        for (size_t c = 0; c < rib1.numControlPoints().at(0) - (c1_merge ? 1 : 0); ++c) {
          for (size_t r = 0; r < rib1.numControlPoints().at(1); ++r) {
            cpts.push_back(rib1.controlPoint(c, r));
          }
        }
        for (size_t c = 1; c < rib2.numControlPoints().at(0); ++c) {
          for (size_t r = 0; r < rib2.numControlPoints().at(1); ++r) {
            cpts.push_back(rib2.controlPoint(c, r));
          }
        }
        Geometry::BSSurface new_rib(deg_u1, deg_v, knots_u, knots_v, cpts);
        new_ribbons[loop].front() = new_rib;
        new_ribbons[loop].pop_back();
        num_merged_sides.front() += num_merged_sides.back();
        num_merged_sides.pop_back();
      }
    }
    num_segments[loop] = num_merged_sides;
  }
  ribbons = new_ribbons;
  init_data();
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
