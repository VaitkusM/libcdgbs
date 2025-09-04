#include "libcdgbs/Example.hpp"

int main(int argc, char* argv[]) {
    std::string filename = "Concave";
    double target_length = 3.0;
    bool merge_smooth_corners = true;
    double deform_value = 0.0;
    // Check if the user provided a filename as a command-line argument
    if (argc >= 2) {
        filename = argv[1];
        if(argc > 2) {
            target_length = std::stod(argv[2]);
            merge_smooth_corners = (argc > 3) ? !(std::string(argv[3]) == "--nomerge") : true;
            deform_value = (argc > 4) ? std::stod(argv[4]) : 0.0;
        }
    }

    libcdgbs::Example example;
    example.say_hello(filename, target_length, merge_smooth_corners, deform_value);
    return 0;
}
