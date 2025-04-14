#pragma once
#include <Eigen/Dense>
#include <triangle.h>
#include <triangle_api.h>

namespace libcdgbs {
    Eigen::MatrixXd m;
    context* triangle_ctx;

    void say_hello();
}
