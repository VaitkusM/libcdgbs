#include "libcdgbs/Example.hpp"

int main(int argc, char* argv[]) {
    std::string filename = "patch-1";
    // Check if the user provided a filename as a command-line argument
    if (argc >= 2) {
        filename = argv[1];
    }

    libcdgbs::Example example;
    example.say_hello(filename);
    return 0;
}
