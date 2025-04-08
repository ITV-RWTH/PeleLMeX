# AMReX
DIM = 3
COMP = gnu
PRECISION = DOUBLE
USE_EB = FALSE
USE_HYPRE = FALSE

# Profiling
PROFILE = FALSE
TINY_PROFILE = TRUE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_SYCL = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PeleLMeX

# Visualisation
USE_CATALYST = TRUE
USE_CONDUIT = TRUE
CONDUIT_DIR = ${EBROOTCONDUIT}
CATALYST_DIR = ${EBROOTCATALYST}

# PelePhysics
Chemistry_Model = air
Eos_Model = Fuego
Transport_Model = Constant

PELE_HOME ?= ../../..
PELE_PHYSICS_HOME ?= /p/project/h2index/kim17/sk_inSitu_Visualisation/PelePhysics
AMREX_HYDRO_HOME ?= /p/project/h2index/kim17/sk_inSitu_Visualisation/AMReX-Hydro
include $(PELE_HOME)/Exec/Make.PeleLMeX


