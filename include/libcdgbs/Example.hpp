#pragma once

extern "C" {
#include <triangle.h>
#include <triangle_api.h>
}
#undef dest
#undef mark

#include <geometry.hh>
#include "Mesh.hpp"
#include "SurfGBS.hpp"
#include <Eigen/Dense>

namespace libcdgbs {
    class Example {
    public:
        Eigen::MatrixXd m;
        context* triangle_ctx;
        Geometry::BSSurface rib;
        Mesh mesh;
        SurfGBS gbs;
        void say_hello();
    };
}
