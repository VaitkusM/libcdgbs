# libcdgbs

**libcdgbs** is a cross-platform C++ library implementing Curved Domain Generalized B-Spline (CD-GBS) patches. It is built using CMake and supports external libraries including:

- [Eigen](https://gitlab.com/libeigen/eigen.git)
- [Triangle](https://github.com/wo80/Triangle)
- [OpenMesh](https://www.graphics.rwth-aachen.de:9000/OpenMesh/OpenMesh)
- [libgeom](https://github.com/salvipeter/libgeom)

## ✨ Features

- ???

## 🧰 Requirements

- C++17 or newer
- CMake 3.15+
- GCC / Clang / MSVC
- Git (for cloning submodules)


## 🧱 Build Instructions

### 1. Clone the repository with submodules

```bash
git clone --recurse-submodules https://github.com/YOUR_USERNAME/libcdgbs.git
cd libcdgbs
```

If you forgot `--recurse-submodules`, run:

```bash
git submodule update --init --recursive
```

### 2. Build with CMake

#### On Linux/macOS:

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

#### On Windows (MinGW):

```bash
mkdir build
cd build
cmake .. -G "MinGW Makefiles"
cmake --build .
```

### 3. Run the test

```bash
./libcdgbs_test
```

---

## 🧩 Project Structure

```
libcdgbs/
├── include/           # Public headers
├── src/               # Source files
├── tests/             # Test executables
├── external/          # Submodules: Triangle, Eigen, OpenMesh, libgeom
├── CMakeLists.txt     # Main build config
└── README.md
```

---

## 📜 License

MIT License. See [LICENSE](LICENSE) for full details.