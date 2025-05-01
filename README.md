# libcdgbs

**libcdgbs** is a cross-platform C++ library implementing [Curved Domain Generalized B-Spline (CD-GBS) patches](https://doi.org/10.1016/j.cagd.2021.102019). It is built using CMake and supports external libraries including:

- [Eigen](https://gitlab.com/libeigen/eigen.git)
- [Triangle](https://github.com/libigl/triangle)
- [OpenMesh](https://www.graphics.rwth-aachen.de:9000/OpenMesh/OpenMesh)
- [libgeom](https://github.com/salvipeter/libgeom)

## ✨ Features

- Usable via libgeom-based API & file import.
- Export into OpenMesh format.

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

## 📝 To-Do

- [x] Multi-loop
- [ ] Support filetypes:
  - [x] MGBS
  - [ ] NGBS
  - [x] CGB
  - [x] GBS
  - [x] MLP
  - [ ] GBP
- [ ] Interior control point
- [ ] Ribbon control points
- [ ] MAT-based interior control

---

## 📜 License

MIT License. See [LICENSE](LICENSE) for full details.

## 📚 References

```bibtex
@article{vaitkus2021multisided,
  author = {M{\'a}rton Vaitkus and Tam{\'a}s V{\'a}rady and P{\'e}ter Salvi and {\'A}goston Sipos},
  title = {Multi-sided B-spline surfaces over curved, multi-connected domains},
  journal = {Computer Aided Geometric Design},
  volume = {89},
  pages = {102019},
  year = {2021},
  publisher = {Elsevier},
  doi = {10.1016/j.cagd.2021.102019},
  url = {https://doi.org/10.1016/j.cagd.2021.102019}
}

@misc{libcdgbs,
  title = {libcdgbs: Cross-platform C++ library for Curved Domain Generalized B-Spline (CD-GBS) patches},
  howpublished = {\\url{https://github.com/VaitkusM/libcdgbs}},
  author={M{\'a}rton Vaitkus},
  year = {2025}
}
```
