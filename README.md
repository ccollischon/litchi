
# Litchi
This repository provides tools for creating Minkowski maps of data in the HEALPix-format. For more info on what that means and how it can be used, please read the following paper: https://www.nature.com/articles/s42005-024-01751-1

If you are using results from litchi in a publication or presentation, please also cite the aforementioned paper.


For details on how to use litchi, see doc/manual.pdf

litchi is archived at Zenodo:    [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11940174.svg)](https://doi.org/10.5281/zenodo.11940174)



# Clone and build
This project depends on HEALPix (https://healpix.sourceforge.io/), which needs to be installed by the user.
Apart from that, it uses pybind11 (https://github.com/pybind/pybind11.git) and eigen (https://gitlab.com/libeigen/eigen.git) as submodules and handles them automatically. Clone both litchi and the submodules with the following commands:

```
git clone https://github.com/ccollischon/litchi.git
cd litchi; git submodule update --init --recursive
cd ..
```

Create a build directory and create the Makefile if your standard compiler can handle C++20 (e.g. g++-10 and above, can be checked with `g++ -v`):
```
cmake -S litchi/ -B litchi-build/
```

Some setups may enable you to load a module that changes your compiler version beforehand, (e.g. `module load gcc/10` on Remeis). Ask your administrator about this.

If you need to specify a different compiler (such as g++-10), name it before calling cmake:
```
CC=gcc-10 CXX=g++-10 cmake -S litchi/ -B litchi-build/
```
Compile:
```
make
```
Litchi requires HEALPix to be installed. Its folder should either be specified in an environment variable called HEALPIX (HEALPix does this if you compile it yourself, check with `echo $HEALPIX`).
Alternatively, it should be in the system’s include/library paths (e.g. after installing with the package manager or adding it manually).

Optional: create code documentation with doxygen by going to litchi/doc and then calling `doxygen doxygen.config`

# Running tests
Other compile targets are `tests`, `testGeometry` and `litchi_debug`, the latter providing an O2 optimized version containing debug symbols and profiling information of litchi. For running the tests (e.g. after modifying the code), a file from the Planck legacy archive is required and can be obtained by 
```
wget -O COM_CMB_IQU-smica_2048_R3.00_hm1.fits \ "http://pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID=COM_CMB_IQU-smica_2048_R3.00_hm1.fits"
```
The tests can then be run by `./tests` and `./testGeometry`.
The python bindings can be tested with a little script (note that using params.Nside to degrade the input is ill-advised in many cases, see e.g. the aforementioned paper):
```
import numpy as np
import litchieat as li
import healpy as hp

filename = "COM_CMB_IQU-smica_2048_R3.00_hm1.fits"
m = hp.read_map(filename)

params = li.paramStruct()
params.Nside = 512
params.rankA = 0
params.rankB = 2
params.curvIndex = 1
params.mint = np.mean(m)
params.numt = 1
params.function = "irrAniso"

params.NsideOut = 16
params.smoothRad = 5. * np.pi / 180

li.makeMinkmap(filename, params, "testbindings.fits")
```