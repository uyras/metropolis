# How to build metropolis

First, make sure that GMP library is installed (https://gmplib.org/).

```
git clone --recurse-submodules https://github.com/uyras/metropolis.git
cd metropolis
mkdir build && cd build
cmake ..
cmake --build .
```

# How to use

Just run `./metropolis --help`, and read.
