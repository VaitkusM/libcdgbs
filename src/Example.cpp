#include "libcdgbs/Example.hpp"
#include <iostream>

namespace libcdgbs {
  void Example::say_hello() {
    gbs = SurfGBS();
    std::cout << "Hello from libcdgbs!" << std::endl;
  }
}
