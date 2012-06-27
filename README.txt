SLRA package - a package for mosaic Hankel structured low-rank approximation 
with interfaces to MATLAB/Octave and R
==============================================================================

Copyright (C) 2012 Ivan Markovsky and Konstantin Usevich 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Ivan Markovsky and Konstantin Usevich
School of Electronics and Computer Science
University of Southampton
Southampton SO17 1BJ
United Kingdom
{im,kdu}@ecs.soton.ac.uk

The package is primarily distributed in the form of source code and is hosted 
at <http://github.com/slra/slra/>. 

0) Prerequisites

GSL library should be installed. To use MATLAB/Octave or R interfaces,
the corresponding environments should also be installed.

SLRA package uses LAPACK and BLAS libraries, which are included in MATLAB, 
Octave and R installations, so if you have MATLAB, Octave or R installed, you 
don't need to install LAPACK and BLAS libraries. The package has an optional 
binding to SLICOT library, which can speed up computations in some cases.

The source files of the libraries can be obtained at
- GSL: <http://www.gnu.org/software/gsl/>
- BLAS, LAPACK: <http://www.netlib.org/>
- SLICOT: <http://www.slicot.org/>
GSL, BLAS and LAPACK libraries are also included in repositories for popular 
Linux distributions.

1) Installation instructions

1.1) Get source
  Download from <http://github.com/slra/slra/> and unpack to a directory
1.2) Compile
  Type in the directory:
    make matlab			--- produces MEX binary file for MATLAB
    make octave			--- produces MEX binary file for Octave
    make R			--- produces R package and installs it
1.3) Install
  MATLAB/Octave MEX files (slra.mex<xxx>) can be installed by copying them
  to directory of your *.m scripts or by using "addpath" command.
  
  If you downloaded a precompiled R package, you should run
    install.packages(repos=NULL, pkgs="Rslra_x.x.xxx"); 
  in R console launched in the same directory where the package is.
  
  The R package should be loaded each time before using it by typing
    library(Rslra);
  in R console.    
    
2) Using the package
Package consists of a single "slra" function, which is documented in the
supplied manual (manual.pdf). An experimental MATLAB version described
in the manual is available in "doc" directory.

A standard help for the function is also available by typing:
  help slra			--- in MATLAB/Octave
  ?slra				--- in R
Directories "test_m" and "test_r" contain demo files for MATLAB/Octave and R.

If you use the package in your research, please cite the following reference:

@@TechReport{slra-package,
  author      = {Markovsky, I. and Usevich, K.},
  title       = {Software for weighted structured low-rank approximation},
  institution =  {ECS, Univ. of Southampton},
  number = 	 {339974},
  address =      {\url{http://ecs.soton.ac.uk/339974}},
  year        = {2012}
}

3) Notes for advanced users

The package contains a demo C++ program, which gives an example of using 
C++ interface and tests various SLRA problems (make target "testc"). 
Documentation for C++ interface can be obtained by running Doxygen
in "cpp" subdirectory.

Advanced compilation options can be found in other targets of makefile, but
not all of targets may run on your machine "as is".

If you wish to try the SLICOT library, download it and copy to SLICOT
subdirectory. You can copy only a few needed files (see makefile). 
Use make targets "xxx-slicot-xxx".

Depending on the version of MATLAB, static binding to non-default 
BLAS and LAPACK (for example, ATLAS) can be faster. Use "xxx-static" target as
a base for your compilation instructions.
