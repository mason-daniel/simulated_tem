# Simulated_TEM 

 Compute simulated TEM images in dynamical two beam condition.
 (c) UKAEA April 2024


We recast the Howie-Whelan equations for generating simulated transmission electron microscope
(TEM) images, replacing the dependence on local atomic displacements with atomic positions only.
This allows very rapid computation of simulated TEM images for arbitrarily complex atomistic
configurations of lattice defects and dislocations in the dynamical two beam approximation.

Full details of the physics captured in this model are available at 
https://arxiv.org/abs/2401.14781


## Dependencies

- Fortran compiler and C compiler
- LAPACK
- libpng
- cmake
- openMPI ( recommended )

## Compilation (Quick Start)

This project is designed to use the cmake build environment.
After downloading and unpacking the source, configure an mpi executable using gfortran with
```
cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$PWD/install -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_BUILD_TYPE=Release .
```
or a serial executable with 
```
cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$PWD/install -DCMAKE_SERIAL=True -DCMAKE_BUILD_TYPE=Release . 
```
then compile with 
```
cmake --build build --target install -j
```
and check the build with
```
cd build
ctest
```
The tests should all be very quick, with the final integration test taking a few seconds.

## Running the code

The basic operation of the code is to load in an .xyz format file and produce a 2d array of intensities corresponding to it's simulated TEM image under given diffraction conditions.

Running without any command line options gives a brief guide.
Note that a boolean flag will also take a "-no" form to negate it. ( eg "-png" to produce a .png file and "-nopng" not to produce one. )
```
> mpiexec -n 1 d2bi
usage
^^^^^
    input/output
        -f <char>                     input filename
        -a0 <float_array>             unit cell dimension (A)
        [-lattice <char> ]            lattice type [ default "bcc" ]
        [-o <char> ]                  output filename
        [-dp <char> ]                 location of diffraction pattern file

    imaging conditions
        -g <float_array>               g-vector (reduced recip space direction [hkl] or [hkil])
        -k <float_array>               zone axis (cell space direction, [uvw] or [uvtw])
        [-sg <float_array> ]           value for deviation parameter s_g (leave unset for fixed k-vector)
        [-ng <float_array> ]           value for deviation parameter defined by (g,ng)
        [-dark <bool_array> ]          dark field imaging mode [ default T,T,T,... ]
        [-theta <float> ]              maximum tilt search angle (deg) [ default 5.000 ]

    cell manipulation
        [-xpad ]                       while reading input file, add padding in x-direction to create a surface [ default F ]
        [-ypad ]                       while reading input file, add padding in y-direction to create a surface [ default F ]
        [-zpad ]                       while reading input file, add padding in z-direction to create a surface [ default F ]
        [-surf ]                       check atom extents for foil thickness instead of box size [ default F ]
        [-d <float_array> ]            imaging space extents. Set 0 to use supercell, unset to use atom extents.
        [-addOffset <float_array> ]    add a global atom position offset (in order to stitch together images)
        [-U <float_array> ]            force crystal rotation matrix in column-major order
        [-R <float_array> ]            force foil tilt rotation matrix in column-major order
        [-opxyz ]                      output .xyz file of atom positions after rotation to two-beam condition [ default F ]

    .png image output
        [-png ]                       output .png file [ default T ]
        [-png_min <float> ]           .png intensity 0 level [ default 0.000 ]
        [-png_max <float> ]           .png intensity 1 level (set 0 for automatic) [ default 0.000 ]
        [-png_blur <float> ]          .png gaussian blur radius (A) [ default 2.500 ]

    calculation
        [-xifile <char> ]              location of extinction distances file
        [-xi0 <float> ]                xi_0 (A) [ default 103.9 ]
        [-xig <float> ]                xi_g (A) [ default 207.5 ]
        [-density ]                    use atomic density field in TEM calc [ default F ]
        [-alterg ]                     use deformation gradient/diffraction pattern to adjust reciprocal lattice vectors [ default F ]
        [-M <int_array> ]              voxel grid per unit cell [ default 4,4 ]

    microscope properties
        [-V <float> ]                  accelerating voltage (keV) [ default 200.0 ]
        [-nPrecAngle <int> ]           number of points take for convergent beam precession [ default 1 ]
        [-precAngle <float> ]          angle (radians) for convergent beam precession [ default 0.5000E-02 ]
        [-nTomoAngle <int> ]           tomography tilt steps [ default 1 ]
        [-tomoAngle <float> ]          tomography tilt half angle (deg) [ default 20.00 ]
```

### input/output

#### -f filename 
Required input option. The file containing atom positions. Acceptable formats are lammps write_data() and extended .xyz files containing the second line containing the supercell lattice repeats, eg of the form
```
Lattice="25.3216 0 0 0 25.3216 0 0 0 25.3216" Properties=species:S:1:pos:R:3
```
If extended .xyz is used, the fourth column is optionally a grain index ( the most populous grain is used to define the g-vector direction ). If the fourth column is not present, all atoms are assumed to be in the same grain.
The fifth to thirteenth columns are optionally the deformation gradient local to each atom ( this information can be used to define an average crystal orientation and so define a g-vector direction ).
If these columns are not set, all atoms are assumed to have identity deformation gradient ( ie identity rotation matrix and zero strain tensor ).
If lammps format is used, all atoms are assumed in grain 1 and with identity deformation gradient.
The deformation gradient data can also be input directly with the -dp flag below.

#### -a0 latt_para_a [,latt_param_b,latt_param_c]
Required input option. Defines the lattice constant. This is used to determine the g-vectors. The number of lattice parameters given should be consistent with the symmetry of the lattice (eg hcp will need 3 numbers).

#### \[-lattice lattice_name]
The lattice symmetry, used to determine the g-vectors. Acceptable values currently "bcc","fcc","hcp"

#### \[-o filename]
The root of the output filename. Taken to be the input filename if not given explicitly.
 
#### \[-dp filename]
Location of crystal grain orientation file specifying lattice orientation. Used to determine the orientation of the crystal, in order to determine g-vectors. File in ascii format. First line gives the integer number of grains (typically 1), subsequent lines give the integer grain number (ignored) and then 9 floating point numbers defining the deformation gradient in column major order.
```
number_of_grains
grain_id T_xx T_yx T_zx T_xy T_yy T_zy T_xz T_yz T_zz
...
```

### imaging conditions

Note that multiple imaging conditions can be defined using a comma separated list.
These are computed one at a time. Examples -g 2,0,0,1,1,0,0,-2,0 or -sg 0.00,0.01,0.02.

#### -g h,k,l or -g h,k,i,l

Defines the g-vector direction as a Miller index ( or Miller-Bravais index for hcp crystal ), using 3(4) comma-separated values eg -g 2,0,0. 

#### -k u,v,w or -k u,v,t,w

Defines the zone axis direction as a Miller index ( or Miller-Bravais index for hcp crystal ).

#### \[ -sg deviation_parameter ]

Specify the deviation parameter s<sub>g</sub> in 1/A units. The system will be tilted in the g-vector direction.
Note the requirement s<sub>g</sub> < |g| / (2 π).
If both -sg and -ng are unset, the system will not be tilted ( ie the input k-vector direction is exactly the beam direction ).


#### \[ -ng deviation_parameter ]

Specify the deviation parameter using the (g,n<sub>g</sub>) notation. The system will be tilted in both the g-vector direction and the normal direction to get good two beam imaging conditions.
If n<sub>g</sub> is integer, two beam conditions are sought with just the diffraction spots [000] and n<sub>g</sub> x [hkl] bright. 
If n<sub>g</sub> is non-integer, then two beam conditions are sought for the integer values floor\[n<sub>g</sub>\] and ceiling\[n<sub>g</sub>\], and a spherical interpolation is used to find a tilt between the two.
If both -sg and -ng are unset, the system will not be tilted ( ie the input k-vector direction is exactly the beam direction ).

#### \[ -dark boolean ]

If set true, output the diffracted beam intensity. If false output the transmitted beam intensity.
Note there is no -nodark option, as 

#### \[ -theta angle ]

Maximum angle that the stage can be tilted to get good two-beam conditions. In degrees, default 5 degrees.

### cell manipulation

#### \[ -\{x|y|z\}pad ]

Add empty space to the system by changing the supercell periodic repeat length in the x- (or y- or z-) direction. This is a quick way of introducing a surface, suitable for when the zone axis is along the x- (or y- or z-) direction. If more control over surfaces are needed then change the input atom position file manually.

#### \[ -surf ]

Indicates the input file has a user- introduced surface rather than periodic boundary conditions.

#### \[ -addOffset dx,dy,dz ]

Introduce a specific offset to all atoms ( as a three-vector in A ). Useful for stitching together images.
 
#### \[ -d \<d1 \[d2,\[d3]]>]

Define the dimensions of the imaging space. The image space is a cuboidal box with fixed axes, with the electron beam in the z direction and the sample rotated and periodically repeated as necessary inside this space.
If unset, then the atom extents (after rotation) are used to define image space extents. 
If set to "-d 0" then the supercell (after rotation)is used to define image space extents.
If one parameter is set "-d z" then the z-direction of image space is set to z A.
If two parameter is set "-d x,y" then the x- and y- directions of image space are set to x,y A. The z length is set by the supercell (after rotation).
If three parameters are set "-d x,y,z" then all three dimensions of image space are set.

#### \[ -U u11,u21,u31,u12,u22,u32,u13,u23,u33 ]

Force the rotation of the simulation cell to be a given rotation matrix


#### \[ -R r11,r21,r31,r12,r22,r32,r13,r23,r33 ]

Force the tilt of the simulation cell to be a given rotation matrix

#### \[ -opxyz]

Output the positions of the atoms used in the tem calculation after rotation and periodic copies made. This allows exact matching between (eg) the positions of defects and the image intensity features, or to check if periodic copies are impinging wrongly on the imaging space.

### .png output

The basic code output is an ascii text file, "\<outfile>.\<diffraction_conditions>.dat".
The filename prefix "\<outfile>" is set by the -o flag, or taken to be equal to the input atom position file name if unset.
The diffraction conditions code indicates the g- and k- vectors used to generate the image, and an indication of the deviation parameter s<sub>g</sub> or n<sub>g</sub>.

The format of the .dat file starts with basic computational parameters, then gives the data in column-major order.
```
# comment line (version number)
# k             u v w
# g             h k l   |g|
# M             Mx Mz
# V (keV)       voltage
# xi (A)        xi0 xig
# n_g           ng
# s_g (1/A)     sg
# img space(A)  d1 d2 d3
# scale (A/px)  scale
# rotation      RU11 RU21 RU31 RU12 RU22 RU32 RU13 RU23 RU33
    NX  NY
intensity_data_in_column_major_order
```
where NX,NY are the number of pixels in the image

#### \[-png]

Convert the .dat file into a greyscale 16-bit .png image.

#### \[-png_min min]

Set the black level (rgb 000000) to .dat file intensity level min.

#### \[-png_max max]

Set the white level (rgb ffffff) to .dat file intensity level max.
If max<=min then the greyscale levels are automatically set by the output intensities.

#### \[-png_blur radius]

Adds a small radius gaussian blurring to the output image, smoothing over small discontinuity artefacts in the output which would not be visible in a "real" microscope. Defaults to a0/2, the phase field smoothing.


### calculation

#### \[-xifile filename]

Reads the extinction distances (in A) from the given file.
The extinction distance file is a plain text file with format
```
# comment line (version number)
 element          g        T (K)      V (V)         Re(U_g)         Im(U_g)            xi_g
      CU   0   0   0     300.000  200000.000     22.67198715      0.59415463      190.05484373
      CU   2   0   0     300.000  200000.000     10.10539531      0.53068873      425.95749593
more_values_for_g-vectors
```

#### \[-xi0 xi0 -xig xig]

If the extinction distances file is not availiable, then just the values for the transmitted and diffracted beam can be passed. 

#### \[-density]

Use a simply calculated atomic density field to determine when the electron beam is passing through free space rather than crystal. 
The atomic phase field is computed at each mesh point as the maximum value of 
    $SmoothStep$\[ 2 - $d/\sigma$ ]$
where $d$ is the distance between an atom and the grid mesh point, and $\sigma=a_0/2$ is a characteristic lengthscale.

#### \[-alterg]

USe the orientation of the crystal provided in the .xyz file or in the orientation file (-dp) to determine the direction of the g-vectors. Defaults to true if the orientation file is given.

#### \[-M m \[,mz]]

Defines the spacing of the phase field grid, $x = Exp[ -i g.r ]$.
m gives the aproximate number of divisions per unit cell ( defined by a0 ). If two arguments are given, then the first is the grid divisions per unit cell in the x-y direction and the second is the grid divisions per unit cell in the z direction. Defaults to m=4.

### microscope properties

#### \[-V voltage_in_keV]

Defines the electron beam voltage in keV. Should agree with the extinction distances (-xi0,-xig or -xifile).

#### \[-nPrecAngle n]

Determines that the precession method will be used to smooth the sensitivity of variations in intensity due to small changes in deviation parameter sg. In addition to the standard calculation, a further n-1 g-vector points are taken in a circle around it. The intensity is averaged. Typically converged for n=10.

#### \[-precAngle angle]

Determines the precession angle used in radians. Defaults to 5 mrad.

#### \[nTomoAngle n]

Determines that tomography will be performed by rotating the sample around the g-vector direction.
n steps will be taken.

#### \[tomoAngle angle]

Determines the maximum tilt angle for the tomography sequence (in degrees). The tilt is taken from -angle through to +angle.



