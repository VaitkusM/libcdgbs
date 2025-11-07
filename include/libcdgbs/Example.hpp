#pragma once
#include "SurfGBS.hpp"

namespace libcdgbs {
    class Example {
    public:
        SurfGBS gbs;
        void say_hello(const std::string& filename = "", double target_length = 3.0, bool merge_smooth_corners = true, double deform_value = 0.0, bool restrict_params = true, bool c1_merge = false, double global_inner_loop_scale = 1.0);
    };
}
