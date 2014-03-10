#TARGETCPU=opteron
OBJECTS=constants.o common.o linalg.o int.o mesh.o main.o green.o rwgf.o quad.o source.o \
	time.o srcint.o bessel.o greenprd.o symmetry.o aux.o sysmat.o solver.o \
	interface.o nfields.o diffr.o bc.o ffields.o cs.o nlsurf.o dipole.o nlbulk.o \
	nfpost.o
COMPILER=ifort -fpp -openmp -parallel -diag-disable 8291 -O2 -vec-report1
#COMPILER=ifort -diag-disable 8291 -vec-report1 -pg

mom: $(OBJECTS)
	$(COMPILER) $(OBJECTS) -o bin/mom -mkl

main.o: main.f90 interface.f90 interface.o
	$(COMPILER) -c main.f90

constants.o: constants.f90
	$(COMPILER) -c constants.f90

common.o: common.f90 source.f90 source.o greenprd.f90 greenprd.o nlsurf.f90 nlsurf.o nlbulk.f90 \
	nlbulk.o
	$(COMPILER) -c common.f90

linalg.o: linalg.f90 constants.f90 constants.o
	$(COMPILER) -c linalg.f90

int.o: int.f90 linalg.f90 linalg.o mesh.f90 mesh.o
	$(COMPILER) -c int.f90

mesh.o: mesh.f90 linalg.f90 linalg.o
	$(COMPILER) -c mesh.f90

green.o: green.f90 linalg.f90 linalg.o
	$(COMPILER) -c green.f90

rwgf.o: rwgf.f90 mesh.f90 mesh.o quad.f90 quad.o
	$(COMPILER) -c rwgf.f90

quad.o: quad.f90 mesh.f90 mesh.o
	$(COMPILER) -c quad.f90

srcint.o: srcint.f90 rwgf.f90 rwgf.o quad.f90 quad.o green.f90 green.o int.f90 int.o greenprd.f90 greenprd.o symmetry.f90 symmetry.o
	$(COMPILER) -c srcint.f90

source.o: source.f90 aux.f90 aux.o bessel.f90 bessel.o mesh.f90 mesh.o symmetry.f90 symmetry.o quad.f90 quad.o rwgf.f90 rwgf.o dipole.f90 dipole.o
	$(COMPILER) -c source.f90

time.o: constants.f90 constants.o time.f90
	$(COMPILER) -c time.f90

bessel.o: bessel.f90
	$(COMPILER) -c bessel.f90

greenprd.o: greenprd.f90 linalg.f90 linalg.o
	$(COMPILER) -c greenprd.f90

symmetry.o: symmetry.f90 linalg.f90 linalg.o
	$(COMPILER) -c symmetry.f90

aux.o: aux.f90 constants.f90 constants.o
	$(COMPILER) -c aux.f90

sysmat.o: sysmat.f90 srcint.f90 srcint.o common.f90 common.o
	$(COMPILER) -c sysmat.f90

bc.o: bc.f90 symmetry.f90 symmetry.o mesh.f90 mesh.o
	$(COMPILER) -c bc.f90

solver.o: solver.f90 sysmat.f90 sysmat.o source.f90 source.o bc.f90 bc.o time.f90 time.o
	$(COMPILER) -c solver.f90

interface.o: interface.f90 solver.f90 solver.o diffr.f90 diffr.o ffields.f90 ffields.o cs.f90 cs.o nfpost.f90 nfpost.o
	$(COMPILER) -c interface.f90

nfields.o: nfields.f90 srcint.f90 srcint.o
	$(COMPILER) -c nfields.f90

diffr.o: diffr.f90 nfields.f90 nfields.o source.f90 source.o common.f90 common.o
	$(COMPILER) -c diffr.f90

ffields.o: ffields.f90 rwgf.f90 rwgf.o symmetry.f90 symmetry.o dipole.f90 dipole.o
	$(COMPILER) -c ffields.f90

cs.o: cs.f90 source.f90 source.o
	$(COMPILER) -c cs.f90

nlsurf.o: nlsurf.f90 rwgf.f90 rwgf.o symmetry.f90 symmetry.o bc.f90 bc.o
	$(COMPILER) -c nlsurf.f90

nlbulk.o: nlbulk.f90 nfields.f90 nfields.o
	$(COMPILER) -c nlbulk.f90

dipole.o: dipole.f90 srcint.f90 srcint.o nfields.f90 nfields.o
	$(COMPILER) -c dipole.f90

nfpost.o: nfpost.f90 nfields.f90 nfields.o source.f90 source.o
	$(COMPILER) -c nfpost.f90

clean:
	rm *.o *.mod
