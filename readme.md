# How to build metropolis

Make sure that following libs are installed:
* [GMP](https://gmplib.org/) with C++-library (gmp и gmpxx);
* OpenMP (included to GCC);
* CMake 3.9 or newer;
* C++17 compiler.

```
git clone --recurse-submodules https://github.com/uyras/metropolis.git
cd metropolis
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target metropolis
```

Run instructions:
```
./build/metropolis --help
```

# How to use

Just run `./metropolis --help`, and read.
