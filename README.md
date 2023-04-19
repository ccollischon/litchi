
# Litchi
This repository provides tools for creating Minkowski maps of data in the HEALPix-format. For more info on what that means and how it can be used, please read the following paper: TODO

If you are using results from litchi in a publication or presentation, please also cite the aforementioned paper.


For details on how to use litchi, see doc/manual.pdf

# Clone and build
This project depends on HEALPix (https://healpix.sourceforge.io/), which needs to be installed by the user.
Apart from that, it uses pybind11 (https://github.com/pybind/pybind11.git) and eigen (https://gitlab.com/libeigen/eigen.git) as submodules and handles them automatically. Clone both litchi and the submodules with the following commands:

```
git clone http://www.sternwarte.uni-erlangen.de/gitlab/collischon/litchi.git
cd litchi; git submodule update --init --recursive
cd ..
```

Create a build directory and create the Makefile if your standard compiler can handle C++20 (e.g. g++-10 and above, can be checked with `g++ -v`):
```
mkdir litchi-build; cd litchi-build;
cmake ../litchi ./
```

Some setups may enable you to load a module that changes your compiler version beforehand, (e.g. `module load gcc/10` on Remeis). Ask your administrator about this.

If you need to specify a different compiler (such as g++-10), name it before calling cmake:
```
CC=gcc-10 CXX=g++-10 cmake ../litchi ./
```
Compile:
```
make
```
Litchi requires HEALPix to be installed. Its folder should either be specified in an environment variable called HEALPIX (HEALPix does this if you compile it yourself, check with `echo $HEALPIX`).
Alternatively, it should be in the system’s include/library paths (e.g. after installing with the package manager or adding it manually).

Optional: create code documentation with doxygen by going to litchi/doc and then calling `doxygen doxygen.config`
