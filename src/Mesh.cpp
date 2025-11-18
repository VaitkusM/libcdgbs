#include "libcdgbs/Mesh.hpp"

using namespace libcdgbs;

// struct ClosestResult {
//   int           triangle = -1;   // index in input
//   Mesh::Point   closest;         // closest point on that triangle
//   Mesh::Point   bary;            // (u,v,w) with u+v+w=1
//   double        dist2    = std::numeric_limits<double>::infinity();
// };

// Helper: clamp t to [0,1]
inline double clamp01(double t) { return t < 0 ? 0 : (t > 1 ? 1 : t); }

// Closest point on segment AB to P
inline Mesh::Point closestPointSegment(const Mesh::Point& A, const Mesh::Point& B, const Mesh::Point& P, double* t_out = nullptr) {
  Mesh::Point AB = B - A;
  double denom = AB.sqrnorm();
  double t = denom > 0 ? dot(P - A, AB) / denom : 0.0;
  t = clamp01(t);
  if (t_out) *t_out = t;
  return A + AB * t;
}

// Robust closest point on triangle ABC to P (Ericson-style; handles degeneracy)
inline void closestPointTriangle(const Mesh::Point& A, const Mesh::Point& B, const Mesh::Point& C,
  const Mesh::Point& P, Mesh::Point& outPoint, Mesh::Point& outBary)
{
  // Check if triangle is degenerate (area ~ 0). Fall back to best edge/vertex.
  Mesh::Point N = cross(B - A, C - A);
  double area2 = N.sqrnorm();
  const double EPS = 1e-18;
  if (area2 < EPS) {
    // Compare three segments
    double tAB, tBC, tCA;
    Mesh::Point pAB = closestPointSegment(A, B, P, &tAB);
    Mesh::Point pBC = closestPointSegment(B, C, P, &tBC);
    Mesh::Point pCA = closestPointSegment(C, A, P, &tCA);
    double dAB = (P - pAB).sqrnorm(), dBC = (P - pBC).sqrnorm(), dCA = (P - pCA).sqrnorm();

    if (dAB <= dBC && dAB <= dCA) {
      outPoint = pAB;
      outBary = Mesh::Point(1.0 - tAB, tAB, 0.0);
    }
    else if (dBC <= dCA) {
      outPoint = pBC;
      outBary = Mesh::Point(0.0, 1.0 - tBC, tBC);
    }
    else {
      outPoint = pCA;
      outBary = Mesh::Point(tCA, 0.0, 1.0 - tCA);
    }
    return;
  }

  // Non-degenerate: project P onto plane, then clamp to triangle via barycentric regions
  Mesh::Point AB = B - A, AC = C - A, AP = P - A;
  double d1 = dot(AB, AP);
  double d2 = dot(AC, AP);
  if (d1 <= 0.0 && d2 <= 0.0) { outPoint = A; outBary = Mesh::Point(1, 0, 0); return; }

  Mesh::Point BP = P - B;
  double d3 = dot(AB, BP);
  double d4 = dot(AC, BP);
  if (d3 >= 0.0 && d4 <= d3) { outPoint = B; outBary = Mesh::Point(0, 1, 0); return; }

  double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    double v = d1 / (d1 - d3);
    outPoint = A + AB * v;
    outBary = Mesh::Point(1.0 - v, v, 0.0);
    return;
  }

  Mesh::Point CP = P - C;
  double d5 = dot(AB, CP);
  double d6 = dot(AC, CP);
  if (d6 >= 0.0 && d5 <= d6) { outPoint = C; outBary = Mesh::Point(0, 0, 1); return; }

  double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    double w = d2 / (d2 - d6);
    outPoint = A + AC * w;
    outBary = Mesh::Point(1.0 - w, 0.0, w);
    return;
  }

  double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    Mesh::Point BC = C - B;
    double denom = dot(BC, BC);
    double w = denom > 0 ? (dot(P - B, BC) / denom) : 0.0;
    w = clamp01(w);
    outPoint = B + BC * w;
    outBary = Mesh::Point(0.0, 1.0 - w, w);
    return;
  }

  // Inside face region
  double denom = (va + vb + vc);
  double v = vb / denom;
  double w = vc / denom;
  double u = 1.0 - v - w;
  outPoint = A + AB * v + AC * w;
  outBary = Mesh::Point(u, v, w);
}


Mesh::ClosestResult Mesh::findClosestFace(Point pt) const
{
  ClosestResult closest;

  for (auto ff : faces()) {
    std::vector<Point> vtx; vtx.reserve(3);
    for(auto v : ff.vertices()) {
      vtx.push_back(point(v));
    }
    
    Point closest_xyz, closest_bary;
    closestPointTriangle(vtx[0], vtx[1], vtx[2], pt, closest_xyz, closest_bary);
    auto dist2 = (closest_xyz - pt).sqrnorm();
    if (dist2 < closest.dist2) {
      closest.face = ff;
      closest.dist2 = dist2;
      closest.pt_xyz = closest_xyz;
      closest.pt_bary = closest_bary;
    }
  }

  return closest;
}

Mesh::Point Mesh::interpolateInFace(FaceHandle ff, Point bary) const
{
  auto f = make_smart(ff);
  std::vector<Point> vtx; vtx.reserve(3);
  for(auto v : f.vertices()) {
    vtx.push_back(point(v));
  }

  auto pt = bary[0] * vtx[0] + bary[1] * vtx[1] + bary[2] * vtx[2];
  return { pt[0], pt[1], 0.0 };
}

double Mesh::dist2Face(Point pt, FaceHandle ff) const
{
  return (pt - calc_face_centroid(ff)).norm(); 
}

void Mesh::computeFaceGradientofFunction(
  std::function<double(VertexHandle)> func,
  std::vector<Point>& outGrads
) const {
  outGrads.clear();
  outGrads.resize(n_faces());

  for (auto ff : faces()) {
    Mesh::Point gradient(0.0, 0.0, 0.0);
    double A = calc_sector_area(halfedge_handle(ff));
    for (auto fh : fh_range(ff)) {
      auto ev = calc_edge_vector(fh);
      auto vv = to_vertex_handle(next_halfedge_handle(fh));
      double u = func(vv);
      gradient += u * Mesh::Point(-ev[1], ev[0], 0.0);
    }
    gradient /= 2.0 * A;

    outGrads[ff.idx()] = gradient;
  }
}

void Mesh::barycentricCoordinates(
  const Point& pt,
  const Point& p1,
  const Point& p2,
  const Point& p3,
  double& u, double& v, double& w)
{
  auto v0 = p2 - p1;
  auto v1 = p3 - p1;
  auto v2 = pt - p1;

  auto d00 = v0.dot(v0);
  auto d01 = v0.dot(v1);
  auto d11 = v1.dot(v1);
  auto d20 = v2.dot(v0);
  auto d21 = v2.dot(v1);

  auto denom = d00 * d11 - d01 * d01;
  const double epsilon = 1e-10;
  if (std::abs(denom) < epsilon) {
    u = v = w = 0.0; // Degenerate or nearly degenerate case
    return;
  }

  u = (d11 * d20 - d01 * d21) / denom;
  v = (d00 * d21 - d01 * d20) / denom;
  w = 1.0 - u - v;
}