#include "libcdgbs/Example.hpp"
#include <iostream>

namespace libcdgbs {
  void Example::say_hello() {
    gbs = SurfGBS();
    std::string filename("two_loops_sharp_flat"); 
    gbs.readMGBS(filename + ".mgbs");
    gbs.compute_domain_mesh();
    gbs.writeOBJ(filename + ".obj");
    std::cout << "Hello from libcdgbs!" << std::endl;
  }
}
