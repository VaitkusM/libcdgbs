#include "libcdgbs/Example.hpp"

int main(int argc, char* argv[]) {
    std::string filename = "threeloops";
    double target_length = 3.0;
    // Check if the user provided a filename as a command-line argument
    if (argc >= 2) {
        filename = argv[1];
        if(argc > 2) {
            target_length = std::stod(argv[2]);
        }
    }

    libcdgbs::Example example;
    example.say_hello(filename, target_length);
    return 0;
}
