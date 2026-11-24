#include "libcdgbs/Example.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {

void printUsage(const char* prog) {
    std::cout
        << "Usage: " << prog << " [example] [options]\n"
        << "\n"
        << "  example                 input file (default: Concave). The extension\n"
        << "                          selects the format: .mgbs, .gbs, .cgb, .ngbs,\n"
        << "                          .gbp or .mlp; without one, these are tried\n"
        << "                          in that order\n"
        << "\n"
        << "Options, in any order (defaults in brackets):\n"
        << "  --target-length <x>     target edge length of the domain mesh [3]\n"
        << "  --deform <x>            deformation value [0]\n"
        << "  --inner-scale <x>       global scale of the inner loops [1]\n"
        << "  --merge | --nomerge     merge ribbons at smooth corners [merge]\n"
        << "  --restrict | --norestrict\n"
        << "                          restrict parameters to the adjacent segment\n"
        << "                          of a merged side [restrict]\n"
        << "  --c1merge | --c0merge   C1 or C0 merge of ribbons [c0merge]\n"
        << "  --usehwidths | --nohwidths\n"
        << "                          use h-widths [usehwidths]\n"
        << "  --arclength | --noarclength\n"
        << "                          arc-length boundary sampling [arclength]\n"
        << "  --mergeinnerc1 | --nomergeinnerc1\n"
        << "                          C1 merge of the inner loops [nomergeinnerc1]\n"
        << "  -h, --help              print this message\n"
        << "\n"
        << "For compatibility, up to three bare numbers are also accepted and\n"
        << "mean, in order: target length, deformation value, inner-loop scale.\n";
}

// Parses a complete floating-point number; false on anything else.
bool parseNumber(const std::string& s, double& value) {
    if (s.empty()) { return false; }
    char* end = nullptr;
    value = std::strtod(s.c_str(), &end);
    return end == s.c_str() + s.size();
}

} // namespace

int main(int argc, char* argv[]) {
    std::string filename = "Concave";
    libcdgbs::SurfGBS::InputParams params;

    struct Flag {
        const char* on;
        const char* off;
        bool* target;
    };
    const Flag flags[] = {
        {"--merge", "--nomerge", &params.merge_smooth_corners},
        {"--restrict", "--norestrict", &params.restrict_params},
        {"--c1merge", "--c0merge", &params.c1_merge},
        {"--usehwidths", "--nohwidths", &params.use_h_widths},
        {"--arclength", "--noarclength", &params.arclength_sampling},
        {"--mergeinnerc1", "--nomergeinnerc1", &params.merge_inner_c1},
    };
    struct Value {
        const char* name;
        double* target;
    };
    const Value values[] = {
        {"--target-length", &params.target_length},
        {"--deform", &params.deform_value},
        {"--inner-scale", &params.global_inner_loop_scale},
    };
    const size_t num_values = sizeof(values) / sizeof(values[0]);

    bool have_filename = false;
    size_t bare_numbers = 0;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        }
        bool matched = false;
        for (const Flag& f : flags) {
            if (arg == f.on || arg == f.off) {
                *f.target = (arg == f.on);
                matched = true;
            }
        }
        for (const Value& v : values) {
            if (arg == v.name) {
                if (i + 1 >= argc || !parseNumber(argv[i + 1], *v.target)) {
                    std::cerr << "Error: " << arg << " needs a number" << std::endl;
                    return 1;
                }
                ++i;
                matched = true;
            }
        }
        if (matched) { continue; }
        double number = 0.0;
        if (parseNumber(arg, number)) {
            // Negative numbers also start with '-', so numbers are tried
            // before the unknown-option check.
            if (bare_numbers >= num_values) {
                std::cerr << "Error: unexpected number '" << arg << "'" << std::endl;
                return 1;
            }
            *values[bare_numbers++].target = number;
        }
        else if (arg[0] == '-') {
            std::cerr << "Error: unknown option '" << arg << "' (see --help)" << std::endl;
            return 1;
        }
        else if (!have_filename) {
            filename = arg;
            have_filename = true;
        }
        else {
            std::cerr << "Error: unexpected argument '" << arg << "' (see --help)" << std::endl;
            return 1;
        }
    }

    std::cout << "Example: " << filename << "\n"
              << "  target length  " << params.target_length << "\n"
              << "  deform value   " << params.deform_value << "\n"
              << "  inner scale    " << params.global_inner_loop_scale << "\n"
              << std::boolalpha
              << "  merge          " << params.merge_smooth_corners << "\n"
              << "  restrict       " << params.restrict_params << "\n"
              << "  c1 merge       " << params.c1_merge << "\n"
              << "  use h-widths   " << params.use_h_widths << "\n"
              << "  arclength      " << params.arclength_sampling << "\n"
              << "  merge inner c1 " << params.merge_inner_c1 << std::endl;

    libcdgbs::Example example;
    return example.say_hello(filename, params) ? 0 : 1;
}
