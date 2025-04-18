#pragma once

#include <geometry.hh>
#include "Mesh.hpp"
#include "SurfGBS.hpp"
#include <Eigen/Dense>

namespace libcdgbs {
    class Example {
    public:
        Eigen::MatrixXd m;
        Geometry::BSSurface rib;
        Mesh mesh;
        SurfGBS gbs;
        void say_hello();
    };
}
