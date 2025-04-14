#include "libcdgbs/Example.hpp"

int main() {
    libcdgbs::Example example;
    example.say_hello();
    example.m.setIdentity(4,5);
    return 0;
}
