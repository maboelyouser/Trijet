# PYTHIA configuration file.
# Generated on 16 ينا, 2024 EET 10:44:13 م with the user supplied options:
# --with-root=/home/mg/workfolder/install/python3
# --with-lhapdf6=/home/mg/workfolder/install/python3
# --with-hepmc3=/home/mg/workfolder/install/python3
# --with-rivet=/home/mg/workfolder/install/python3
# --with-fastjet3=/home/mg/workfolder/install/python3
# --with-python=/usr/include/python3.10

# Install directory prefixes.
PREFIX_BIN=/home/mg/workfolder/src/pythia8308/bin
PREFIX_INCLUDE=/home/mg/workfolder/src/pythia8308/include
PREFIX_LIB=/home/mg/workfolder/src/pythia8308/lib
PREFIX_SHARE=/home/mg/workfolder/src/pythia8308/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=true
FASTJET3_CONFIG=fastjet-config
FASTJET3_BIN=/home/mg/workfolder/install/python3/bin/
FASTJET3_INCLUDE=-I/home/mg/workfolder/install/python3/include
FASTJET3_LIB=-L/home/mg/workfolder/install/python3/lib -Wl,-rpath,/home/mg/workfolder/install/python3/lib -lfastjet

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=true
HEPMC3_CONFIG=HepMC3-config
HEPMC3_BIN=/home/mg/workfolder/install/python3/bin/
HEPMC3_INCLUDE=-I/home/mg/workfolder/install/python3/include
HEPMC3_LIB=-L/home/mg/workfolder/install/python3/lib -Wl,-rpath,/home/mg/workfolder/install/python3/lib -lHepMC3

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=true
LHAPDF6_CONFIG=lhapdf-config
LHAPDF6_BIN=/home/mg/workfolder/install/python3/bin/
LHAPDF6_INCLUDE=-I/home/mg/workfolder/install/python3/include
LHAPDF6_LIB=-L/home/mg/workfolder/install/python3/lib -Wl,-rpath,/home/mg/workfolder/install/python3/lib -lLHAPDF

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=true
RIVET_CONFIG=rivet-config
RIVET_BIN=/home/mg/workfolder/install/python3/bin/
RIVET_INCLUDE=-I/home/mg/workfolder/install/python3/include
RIVET_LIB=-L/home/mg/workfolder/install/python3/lib -Wl,-rpath,/home/mg/workfolder/install/python3/lib -lRivet

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=/home/mg/workfolder/install/python3/bin/
ROOT_INCLUDE=-I/home/mg/workfolder/install/python3/include
ROOT_LIB=-L/home/mg/workfolder/install/python3/lib -Wl,-rpath,/home/mg/workfolder/install/python3/lib -lCore

GZIP_USE=false
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=true
PYTHON_CONFIG=python-config
PYTHON_BIN=/usr/include/python3.10/
PYTHON_INCLUDE=-I/usr/include/python3.10
PYTHON_LIB=-L/usr/include/python3.10 -Wl,-rpath,/usr/include/python3.10

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
