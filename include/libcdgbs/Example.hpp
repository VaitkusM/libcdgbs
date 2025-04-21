#pragma once
#include "SurfGBS.hpp"

namespace libcdgbs {
    class Example {
    public:
        SurfGBS gbs;
        void say_hello(const std::string& filename = "");
    };
}
