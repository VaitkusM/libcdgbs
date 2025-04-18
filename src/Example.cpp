#include "libcdgbs/Example.hpp"
#include <iostream>

namespace libcdgbs {
  void Example::say_hello() {
    gbs = SurfGBS();
    gbs.readMGBS("one_loop_flat.mgbs");
    std::cout << "Hello from libcdgbs!" << std::endl;
  }
}
