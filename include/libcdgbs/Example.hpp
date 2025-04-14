#pragma once
#include <Eigen/Dense>
#include <triangle.h>
#include <triangle_api.h>
#include <geometry.hh>

namespace libcdgbs {
    Eigen::MatrixXd m;
    context* triangle_ctx;
    Geometry::BSSurface rib;
    void say_hello();
}
