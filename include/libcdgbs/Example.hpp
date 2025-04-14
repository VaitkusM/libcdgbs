#pragma once
#include <Eigen/Dense>
#include <triangle.h>
#include <triangle_api.h>
#include <geometry.hh>
#include "Mesh.hpp"
#include "SurfGBS.hpp"

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
