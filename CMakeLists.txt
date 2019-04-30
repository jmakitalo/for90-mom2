cmake_minimum_required (VERSION 2.6)
project (for90-mom2)
enable_language (Fortran)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin/")

set(BUILD_TYPE "release")

if(BUILD_TYPE MATCHES "release")
 message("Building release")

 # Support OpenMP.
 # May want to add also "-vec-report1" to get vectorization reports.
 if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
	set(CMAKE_Fortran_FLAGS "-fpp -openmp -parallel -O2")
 elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	set(CMAKE_Fortran_FLAGS " -fopenmp -O2 -ffree-line-length-0")
 endif()
elseif(BUILD_TYPE MATCHES "debug")
 message("Building debug")

 if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
	set(CMAKE_Fortran_FLAGS "-g")
 elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	set(CMAKE_Fortran_FLAGS "-g -fcheck=all -fbounds-check -fbacktrace -ffree-line-length-0")
 endif()
elseif(BUILD_TYPE MATCHES "profile")
 message("Building profile")

 set(CMAKE_Fortran_FLAGS "-pg")
else()
	message(FATAL_ERROR "Invalid build type!")
endif()


#---- run libraries check ----
message("---- Checking installed libraries")
# check MKL root folder
if(DEFINED ENV{MKLROOT})
	# set bash environment variable to CMAKE variable
	set(MKLROOT "$ENV{MKLROOT}")
	message("---- -- The intel mkl root folder: ${MKLROOT}")
else()
	message("---- -- The MKLROOT is not set, use default folder!")
	set(MKLROOT "/worktmp/intel/mkl")
endif()

 if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
	set(CMAKE_EXE_LINKER_FLAGS "-mkl")
 elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	set(CMAKE_EXE_LINKER_FLAGS "-Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -L${MKLROOT}/../compiler/lib/intel64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -liomp5 -lpthread")
 endif()

add_executable(mom
f77toms644/zbsubs.f
f77toms644/zbesyh.f
f77toms644/machcon.f
aux.f90
constants.f90
ffields.f90
interface.f90
mesh.f90
nlsurf.f90
solver.f90
strat.f90
vie.f90
bc.f90
cs.f90
int.f90
nfields.f90
pec.f90
solver_vie.f90
symmetry.f90
bessel.f90
diffr.f90
green.f90
linalg.f90
nfpost.f90
quad.f90
source.f90
sysmat.f90
common.f90
dipole.f90
greenprd.f90
main.f90
nlbulk.f90
rwgf.f90
srcint.f90
time.f90
)
