#include "libcdgbs/Example.hpp"

int main(int argc, char* argv[]) {
    std::string filename = "Concave";
    double target_length = 3.0;
    bool merge_smooth_corners = true;
    // Check if the user provided a filename as a command-line argument
    if (argc >= 2) {
        filename = argv[1];
        if(argc > 2) {
            target_length = std::stod(argv[2]);
            merge_smooth_corners = (argc > 3) ? !(std::string(argv[3]) == "--nomerge") : true;
        }
    }

    libcdgbs::Example example;
    example.say_hello(filename, target_length, merge_smooth_corners);
    return 0;
}
