cmake_minimum_required (VERSION 2.6)
project (for90-mom2)
enable_language (Fortran)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin/")

set(BUILD_TYPE "release")

if(BUILD_TYPE MATCHES "release")
 message("Building release")

 # Support OpenMP.
 # May want to add also "-vec-report1" to get vectorization reports.
 set(CMAKE_Fortran_FLAGS "-fpp -openmp -parallel -O2")
elseif(BUILD_TYPE MATCHES "debug")
 message("Building debug")

 set(CMAKE_Fortran_FLAGS "-g")
elseif(BUILD_TYPE MATCHES "profile")
 message("Building profile")

 set(CMAKE_Fortran_FLAGS "-pg")
else()
	message(FATAL_ERROR "Invalid build type!")
endif()

set(CMAKE_EXE_LINKER_FLAGS "-mkl")

add_executable(mom
amos.f
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
