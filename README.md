FOR90-MOM2 USER'S GUIDE
=======================

Author information
------------------

Jouni Mäkitalo (jmakitalo15@gmail.com)  
Researcher at Tampere University of Technology (2011-2015)  
Department of Physics  
Optics Laboratory

References
----------

The code was developed for my Ph.D. thesis:

"Boundary Integral Operators in Linear and Second-order Nonlinear Nano-optics"
http://URN.fi/URN:ISBN:978-952-15-3539-0 

Disclaimer
----------

Although the program has been used to successfully model several cases and the results compare well with the results for a spherical particle, the author does not take any responsibility for the correctness of results obtained by the end-user.

Preface
-------

This is the user's guide to for90-mom2 Method of Moments electromagnetic wave scattering solver implemented in Fortran 90.

This program is based on the theory and expressions derived in Master of Science Thesis (Jouni Mäkitalo, 2011).
The nonlinear surface response code was developed in the Ph.D. Thesis (Jouni Mäkitalo, 2015).

The program is made for solving linear and weakly nonlinear scattering problems, where the particle size is on the order of the wavelength and whose material is described by a piece-wise constant complex permittivity. The following list gives the main features and application targets:

- Computation of full-wave solutions to electromagnetic scattering problems
- PMCHWT boundary element formulation of dielectric and lossy piece-wise homogeneous scatterers
- Deduce scattering and absorption cross-sections, near-fields and radiation patterns
- Excitations: plane-wave, focused Gaussian beams, electric dipole
- Enforce symmetries and decompose solutions from their irreducible representations
- Second-order nonlinear scattering with surface electric dipolar source
- Second-order nonlinear scattering with bulk dipolar source
- Multidomain problems: a domain for each constant permittivity, the domains may
  be nested and may be touching
- Problems with 2-D periodicity in the xy-plane, deduce near-fields, reflectivity and transmittance
- Solve a problem for multiple sources quickly by reusing matrix factorization.
- Nonlinear beamscanning: deduce the reflected second-harmonic power when a focused beam is
  scanned over a rectangular domain
- Eigenmodes of flat PEC scatterers/antennas.
- Volume integral formulation of scattering from inhomogeneous objects. Uses SWG basis and Duffy transform.

Building
--------

The code has been successfully built with the Intel Fortran compiler (ifort 12.1.2). It should build with other compilers, but this has not been tested.

The code has a dependency on LAPACK library. During the development, Intel's Math Kernel Library packaging of LAPACK was used.

The file structure of the code base is

 .git  
 *.f90  
 amos.f  
 CMakeLists.txt  

The amos.f is a Fortran 77 file that is part of the Netlib package. It implements various Bessel functions of complex argument.

To build, go to the source directory and type:

mkdir bin  
mkdir build  
cd build  
FC=ifort cmake ..  
make  

The executeable named mom will be placed in bin/ directory. The most probable reason for the build to fail is that the Math Kernel Library may be inaccessable. Revise the CMakeLists.txt for appropriate modifications.

Basics usage
------------

Run the program interactively by typing

./mom

The program will display a prompt and will receive commands. In Unix based operating systems you can pass the commands via a text file with stdin redirection. Say you have a file "input" which includes the commands that describe the probelm to be solved. Then just type

./mom < input

The program supports symmetric multiprocessing via OpenMP. To exploit this, set the environment variable

OMP_NUM_THREADS=n

to get n threads. If you run the program in the Sun grid environment (as in MGRID Akaatti),
the you can use the following script (named "myscript"):

 #!/bin/sh  
 #$ -S /bin/bash  
 #$ -N batch-name  
 #$ -cwd  
 #$ -j y  
 #$ -o output  
 #$ -pe openmp 8  
 export OMP_NUM_THREADS=$NSLOTS  
 /home/mydirectory/bin/mom < input  

To submit a batch with this script, type

qsub myscript

In the Slurm queue system (MGRID Merope), the script is of the form

 #!/bin/bash  
 #SBATCH -J batch-name  
 #  
 # Let's redirect job's out some other file than default slurm-%jobid-out  
 #SBATCH --output=output  
 #   
 # Maximum time.  
 #SBATCH --time=2:00:00  
 #SBATCH --ntasks=12  
 #  
 # Memory per node in MiB  
 #SBATCH --mem=2000  
 export OMP_NUM_THREADS=12  
 /home/mydirectory/bin/mom < input  

To submit a batch with this script, type

sbatch myscript

Commands
--------

This section explains the commands implemented in the main interface of the program. Note that the usage of the program should not be limited to the use of this interface. The interface only provides actions for the most common use of the program as was required during the development phase. The user may add new functionality to the source code and if necessary, add a new interface command for this.

Wavelengths are used to fix the frequency of the time-harmonic problem. This is mostly because usually resonances are intuitively related to size of geometrical features. Thus the wavelengts referenced below always correspond to the wavelength in vacuum.

In the following, brackets/braces are used to denote formal arguments of commands. Do not include the brackets in the actual arguments. A bracketed argument is to be replaced by only single actual argument and a braced argument can be replaced by multiple actual arguments.

The implementation uses SI units, so all parameters should represent values in this unit system.

COMMAND: name [str]  
DESCRIPTION:  
Sets the name of the computation. The name is used mainly to create output filenames. Avoid spaces and special characters in the name as this might cause trouble depending on the filesystem in use.

Example: name 'myscatteringproblem'

COMMAND: mesh [filename] [scale]  
DESCRIPTION:  
Loads the surface mesh that is used to construct the basis functions. The [filename] should point to a valid mesh file of some of the following formats:

- Gmsh mesh (*.msh)
- Netgen neutral mesh format (*.nmf)

The [scale] is a real number that is used to scale the points in the meshfile into points in the computational domain, which has units of meters.

Example: mesh 'sphere.msh' 1E-9

COMMAND: quad [tri] [tetra]  
DESCRIPTION:  
Quadrature rule for triangle and tetrahedron elements. Possible choises for [tri] are
tri_gl1
tri_gl3
tri_gl4
tri_gl7
tri_gl13

and for [tetra]
tetra_gl1
tetra_gl4
tetra_gl11

The number refers to the number of points.

COMMAND: wlrg [nwl] [first] [last]  
DESCRIPTION:  
Allocates datastructures for amount of [nwl] computations, which will have wavelengths that are equally spaced from [first] to [last], including the endpoints.

Example: wlrg 10 400E-9 1200E-9

COMMAND: wlls [nwl] {wavelengths}  
DESCRIPTION:  
Allocated datastrctures for amount of [nwl] computations, which will have wavelengths given by the space separated list {wavelengths}.

Example: wlls 3 300E-9 1.7E-6 2.3E-6

COMMAND: nsrc [n]  
DESCRIPTION:  
Allocates space for [n] sources. Must be called prior to ssrc.

COMMAND: ssrc [index] [type] {parameters}  
DESCRIPTION:  
Define the source of index [index] in the scattering problem. Currently following parameters are valid:

[type] {parameters} ; description

'pw' [theta] [phi] [psi] ; A time-harmonic linearly polarized plane-wave
'focus_rad' [focal] [waist] [na] [normalize] ; Focused radially polarised beam.
'focus_azimut' [focal] [waist] [na] [normalize] ; Focused azimutally polarized beam
'focus_x' [focal] [waist] [na] [normalize] ; x-polarised focused gaussian beam
'focus_y' [focal] [waist] [na] [normalize] ; y-polarised focused gaussian beam
'focus_hg01' [focal] [waist] [na] [normalize] ; x-polarised focused Hermite-gaussian beam, mode HG_01
'dipole' [x] [y] [z] [dx] [dy] [dz] ; An electric dipole source at point (x,y,z) with dipole moment (dx,dy,sz)

The angle theta denotes the elevation angle i.e. angle between the z-axis and direction of propagation k. Angle phi denotes the azimuthal angle measured with respect to x-axis towards the y-axis. Angle psi denotes the polarization angle, measured from k x z. Angles are in degrees. To be more precise:

kx = k0*sin(theta)*cos(phi)
ky = k0*sin(theta)*sin(phi)
kz = k0*cos(theta)

where k0 = 2*pi/lambda

Ex = E0*(cos(psi)*sin(phi) - sin(psi)*cos(theta)*cos(phi))
Ey = E0*(-cos(psi)*cos(phi) - sin(psi)*cos(theta)*sin(phi))
Ez = E0*sin(psi)*sin(theta)

Parameters [focal], [waist] and [na] denote focal length, beam waist and numerical aperture, respectively. If parameter [normalize] is 'true', then the fields are scaled by the maximum electric field at the focus. If parameter [normalize] is 'false', then the fields are normalized to the incident plane-wave incident of the focusing lens.

COMMAND: scan [sz] [d] [npt]  
DESCRIPTION:  
Defines a special set of sources, where a source defined by a single 'ssrc' command is scanned over a square area of side length [d] over plane z = [sz]. The number of scan points is [npt]*[npt]. Example usage:

nsrc 1
ssrc 1 'focus_rad' 1e-3 0.8e-3 0.8 false
scan 0 2e-6 20

COMMAND: nmed [n]  
DESCRIPTION:  
Allocates space for [n] media descriptors.

COMMAND: smed [mindex] [type] [method] {parameters}  
DESCRIPTION:  
Sets the medium of index [mindex]. The medium type [type] may be one of

linear
nlsurf
nlbulk_nonlocal

The paramter [method] describes how a property is set. Valid choices are

value
file

For type 'linear', the {parameters} descrive the complex index of refraction. If method 'value' is used, then {parameters} is the numeric value of the complex index of refraction in Fortran syntax and the value is used for all wavelengths. If method is 'file', then the indices of refraction are read from a file, whose name is given as {parameters}. The file must contain three columns: wavelength in meters, real part of refractive index, imaginary part of refractive index (non-negative).

For type 'nlsurf', only method 'value' is valid. Then {parameters} = [nnn] [ntt] [ttn] where the three numbers are the corresponding second-order surface susceptiblity components for isotropic surfaces.

For type 'nlbulk_nonlocal', only method 'value' is valid. Then {parameters} = [zzz] which is the zzz-component of the dipolar bulk second-order susceptibility (other components are zero).

COMMAND: ndom [n]  
DESCRIPTION:  
Sets the number of domains to [n].

COMMAND: sdom [dindex] [nsurf] [surf_ids] [nvol] [vol_ids] [refind_file] [gf_index]  
DESCRIPTION:  
Sets domain parameters for domain of index [dindex]. Command ndom must be issued before this one. [nsurf] is the number of surfaces whose union constitutes the boundary of the domain in the mesh file. [surf_ids] are integers corresponding to the physical id:s of the surfaces. The integers may be also negative, which is interpreted so that the orientation of the corresponding surface is reversed. [nvol] in the number of volumes whose union constitutes the domain. This may be zero if the volume information is not required. [vol_ids] are integers corresponding to the physical id:s of the volumes (orientation is arbitrary). Parameter [refind_file] is a file name of the refractive index data to be used for the domain. Parameter [gf_index] is an index to a Green's function that has been loaded by command ipgf. A value of -1 corresponds to the nonperiodic Green's function.

It is required that the boundary surface of a domain is oriented so that the normals point into the domain. Thus, if a surface id appears in one sdom expression, it must appear in another one with a minus sign.

It is important to note that the domain, where the incident field is defined, always corresponds to dindex = 1.

COMMAND: fdom  
DESCRIPTION:  
Finalizes the domain setup by dividing the given mesh into proper submeshes while keeping book on the edge connectivity. It also orients the basis function on the submeshes properly. This should be called right before solv.

COMMAND: sbnd [id] [bndname]  
DESCRIPTION:  
Assigns boundary values for edges with physical group [id]. The last argument is a string from the following list

prdx1
prdx2
prdy1
prdy2

which denote the two ends of the unit cell in z = constant plane. E.g. prdx1 corresponds to edge x = -px/2, where px is the period in the direction of the x-axis. Basis coefficients for edges on prdx1 and prdy1 are removed from the system and added to coefficients related to edges prdx2 and prdy2 with proper Bloch phase shifts. These conditions must be set if the mesh touches the unit cell boundaries.
This routine must be called right after the 'mesh' command.

COMMAND: symm [nsubgroups] {names}  
DESCRIPTION:  
Defines the symmetry of the geometry. The symmetry is described by a group, which is generated from the given subgroups. The number of these subgroups is [nsubgroups]. A number of [nsubgroups] names must follow. The list of valid names is


id	identity
mxp	mirror symmetry with respect to the x-plane
myp	mirror symmetry with respect to the y-plane
mzp	mirror symmetry with respect to the z-plane
r[n]	[n]-fold rotation symmetry with respect to the z-axis

It is required that the symmetry group is commutative. Thus, for examples, command

symm 2 mzp r3

is valid but command

symm 2 mxp r3

is not valid. No error checking is done, so this is user's responsibility.

COMMAND: solv  
DESCRIPTION:  
Solves the problem.

COMMAND: nfms [wlindex] [srcindex] [dindex]  
DESCRIPTION:  
Computes the electric and magnetic fields corresponding to source [srcindex] on the surface of domain [dindex] of the particle at wavelength determined by [wlindex]. Produces a msh-file which can be inspected in gmsh.

COMMAND: crst  
DESCRIPTION:  
Computes scattering and absoprtion cross-sections at all available wavelengths and sources.

COMMAND: rcst [wlindex] [srcindex] [ntheta] [nphi]  
DESCRIPTION:  
Computes the bi-static radar cross-section (a.k.a. scattering power per unit solid angle) using the solution denoted by [wlindex] and source of index [srcindex]. Parameters [ntheta] and [nphi] denote how many angular points are evaluated for elevation and azimuthal angles respectively.

COMMAND: rcs2 [wlindex]  
DESCRIPTION:  
Computes the RCS integrated over the solid angle determined by a focused beam numerical aperture (in reflection). The solution at wavelength corresponding to [wlindex] is used.

COMMAND: npgf [n]  
DESCRIPTION:  
Allocate memory for [n] Green's functions. This must be called if periodic problems are to be solved.

COMMAND: ipgf [index] [filename]  
DESCRIPTION:  
Loads the precalculated coefficients for a periodic Green's function from file given by [filename] and assigns the data to positive integer [index]. Command npgf must be issued prior to calling ipgf.

The Green's function data is pre-computed by a MATLAB code made by the author. See the function file gpWlRange.m.

COMMAND: diff [srcindex] [dindex] [orderx] [ordery] [pol] [polangle]  
DESCRIPTION:  
Computes the diffracted power in the given order in a periodic problem. Integer [srcindex] refers to the excitation source. The integer [dindex] denotes the domain, where the calculation is done, i.e., which field expressions are used. Integers [orderx] and [ordery] denote the diffraction orders. [pol] may be 'true' or 'false' to specify whether a linear polarization filter is used when recording the diffracted wave. [polangle] is the polarizer angle in degrees (see diffr.f90 for the actual definition of the angle).

It is assumed that the incident plane-wave propagates from half-space z>0 to half-space z<0.

Example input file
------------------

Say you want to compute the cross-sections of a sphere. Then pass in this input file:

name 'sphere'  
symm 1 id  
mesh 'unitsphere.msh' 200d-9  
wlrg 100 400d-9 1200d-9  
nmed 2  
smed 1 linear value (1.0, 0.0)  
smed 2 linear file aujc.ref  
ndom 2  
sdom 1 1 15 0 1 -1  
sdom 2 1 -15 0 2 -1  
fdom  
nsrc 1  
ssrc 1 'pw' 0 90 0  
solv  
crst  
exit  

Here the physical id of the surface is 15 and is originally oriented so that the normals point into the surrounding space.

A more involved case would be a periodic problem, where a particle stands on a substrate.
The input would look something like this:

name 'myproblem'  
symm 1 id  
wlrg 100 4e-7 10e-7  
mesh 'mymesh.msh' 1e-9  
nsrc 1  
ssrc 1 'pw' 180 0 0  
nmed 3  
smed 1 linear value (1.0, 0.0)  
smed 2 linear file aujc.ref  
smed 3 linear file silica.ref  
npgf 2  
ipgw 1 'vacuum.pgf'  
ipgw 2 'silica.pgf'  
ndom 3  
sdom 1 2 187 188 1 1  
sdom 2 2 -188 -189 2 -1  
sdom 3 2 -187 189 3 2  
fdom  
solv  
diff 1 3 0 0 'false' 0  
exit  

The glass substrate occupies the negative z half-space and the plane-wave is incident in the negative z direction.

File format specifications
--------------------------

Here brackets denote a row or a matrix of data and braces denote an arbitrary sequence of data.

EXTENSION: crs
DESCRIPTION:
Scattering and absorption cross-sections.
CONTENT:
[wavelengths (float)] [scattering cross-sections (float)] [absorption cross-sections (float)]

EXTENSION: dif
DESCRIPTION:
Diffracted intensity in zeroth order normalized to intensity of the incident plane-wave.
CONTENT:
[wavelengths (float)] [intensity (float)]
