#include "libcdgbs/TriangleWrapper.hpp"
#include <cstdlib>
#include <cstring>
#include <map>

namespace libcdgbs {

  TriangleWrapper::TriangleWrapper(char const* options) {
    ctx_ = triangle_context_create();
    triangle_context_options(ctx_, const_cast<char*>(options));
  }

  TriangleWrapper::~TriangleWrapper() {
    triangle_context_destroy(ctx_);
  }

  Mesh TriangleWrapper::triangulate(
    std::vector<std::vector<std::vector<Eigen::Vector3d>>> const& loops,
    double L_target
  ) {
    // Prepare the input struct
    triangleio in{};
    std::map<std::pair<double, double>, int> index;
    std::vector<std::array<double, 2>> pts2D;
    std::vector<std::pair<int, int>> segs;

    // Flatten loops → sub‑curves → edges
    for (auto const& loop : loops) {
      for (auto const& sub : loop) {
        for (size_t i = 0; i + 1 < sub.size(); ++i) {
          auto const& P = sub[i], & Q = sub[i + 1];
          std::pair<double, double> a{ P.x(),P.y() }, b{ Q.x(),Q.y() };
          // get or create index for a
          int ia = index.count(a) ? index[a] : (index[a] = pts2D.size(), pts2D.push_back({ a.first,a.second }), int(pts2D.size() - 1));
          // same for b
          int ib = index.count(b) ? index[b] : (index[b] = pts2D.size(), pts2D.push_back({ b.first,b.second }), int(pts2D.size() - 1));
          segs.emplace_back(ia, ib);
        }
      }
      // close the loop
      if (!loop.empty() && !loop.front().empty()) {
        auto const& F = loop.front().front();
        auto const& L = loop.back().back();
        auto fa = std::make_pair(F.x(), F.y());
        auto lb = std::make_pair(L.x(), L.y());
        int ifa = index[fa], ilb = index[lb];
        if (ifa != ilb) segs.emplace_back(ilb, ifa);
      }
    }

    // fill in.pointlist
    in.numberofpoints = int(pts2D.size());
    in.pointlist = (REAL*)std::malloc(in.numberofpoints * 2 * sizeof(REAL));
    for (int i = 0; i < in.numberofpoints; ++i) {
      in.pointlist[2 * i] = pts2D[i][0];
      in.pointlist[2 * i + 1] = pts2D[i][1];
    }

    // fill in.segmentlist
    in.numberofsegments = int(segs.size());
    in.segmentlist = (int*)std::malloc(in.numberofsegments * 2 * sizeof(int));
    for (int i = 0; i < in.numberofsegments; ++i) {
      in.segmentlist[2 * i] = segs[i].first;
      in.segmentlist[2 * i + 1] = segs[i].second;
    }

    // 1) triangulate using the wo80 API
    // Recompute options for this run: no‐Steiner on loops + quality + area limit
    double A = (std::sqrt(3.0) / 4.0) * L_target * L_target;
    std::ostringstream opts;
    opts << "pQzYq20a" << A;
    std::string s = opts.str();
    triangle_context_options(ctx_, const_cast<char*>(s.c_str()));
    triangle_mesh_create(ctx_, &in);

    // 2) Copy it out into a standard triangleio
    triangleio out{};
    triangle_mesh_copy(ctx_, &out, /*edges=*/1, /*neighbors=*/0);

    // 3) Build an OpenMesh TriMesh
    Mesh mesh;
    std::vector<Mesh::VertexHandle> vh;
    vh.reserve(out.numberofpoints);
    for (int i = 0; i < out.numberofpoints; ++i) {
      vh.push_back(mesh.add_vertex({
        float(out.pointlist[2 * i]),
        float(out.pointlist[2 * i + 1]),
        0.0f
        }));
    }
    for (int i = 0; i < out.numberoftriangles; ++i) {
      mesh.add_face({
        vh[out.trianglelist[3 * i]],
        vh[out.trianglelist[3 * i + 1]],
        vh[out.trianglelist[3 * i + 2]]
        });
    }

    // 4) Clean up
    triangle_free(in.pointlist);
    triangle_free(in.segmentlist);
    triangle_free(out.pointlist);
    triangle_free(out.trianglelist);

    return mesh;
  }

} // namespace libcdgbs
