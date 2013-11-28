SLRA: package for weighted mosaic Hankel structured low-rank approximation 
==============================================================================
### with interfaces to MATLAB/Octave and R

License
-------

Copyright (C) 2012-13 Ivan Markovsky and Konstantin Usevich 

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

This program contains the Timer class for cross-platform time measurements.
(Copyright (c) 2003 Song Ho Ahn, (song.ahn@gmail.com).)

    Ivan Markovsky and Konstantin Usevich
    Vrije Universiteit Brussel
    Department ELEC
    Pleinlaan 2
    1050 Elsene
    Belgium
    {imarkovs,kusevich}@vub.ac.be

Downloading and installing the package
--------------------------------------
The package is primarily distributed in the form of source code and is hosted 
at <http://github.com/slra/slra/>. For some platforms, precompiled binaries
are also available in the repository.

The package can be downloaded from <http://slra.github.io/software.html>,
or directly from <http://github.com/slra/slra/>.

If the precompiled binaries are available for you platform, you should

  * In MATLAB/Octave: add the whole directory to the MATLAB path
    with the `addpath` command 
  * In R (for Windows): run
  
      `install.packages(repos=NULL, pkgs="Rslra_x.x.xxx");` 
    
	in the R console launched in the same directory where the package is.	  
   
Otherwise, please follow instructions in section *Installing the package from source*

Using the package
-----------------

The package consists of a single `slra` function, which is documented in
the supplied manual `doc/slra.pdf`. 

A standard help for the function is also available by typing:
* `help slra` in MATLAB/Octave
* `?slra` in R

Directories `test_m` and `test_r` contain demo files for MATLAB/Octave and R.

If you use the package in your research, please cite the following reference:

    @Article{slra-software,
		author = {I. Markovsky and K. Usevich},
    	title = {Software for weighted structured low-rank approximation},
    	journal = {J. Comput. Appl. Math.},
		volume = {256},
		pages = {278--292},
		year = {2014},
    }
	
	

Installing the package from source
----------------------------------
### Prerequisites

GSL library should be installed (if not using a precompiled Windows package
for R). To use MATLAB/Octave or R interfaces, the corresponding environments
should also be installed.

SLRA package uses LAPACK and BLAS libraries, which are included in MATLAB, 
Octave and R installations, so if you have MATLAB, Octave or R installed, you 
don't need to install LAPACK and BLAS libraries. The package has an optional 
binding to SLICOT library, which can speed up computations in some cases.

The source files of the libraries can be obtained at
* GSL: <http://www.gnu.org/software/gsl/>
* BLAS, LAPACK: <http://www.netlib.org/>
* SLICOT: <http://www.slicot.org/>
GSL, BLAS and LAPACK libraries are also included in repositories for popular 
Linux distributions.

### Installation from source

1. Get source

   Download from <http://github.com/slra/slra/> and unpack to a directory

2. Compile
	 * type `make matlab` to produce a MEX binary file for MATLAB
	 * type `make matlab-win` to produce a MEX binary file in Windows
	 * type `make octave` to produce a MEX binary file for Octave
	 * type `make R` to produce an R package and install it 

3. Install
   * In MATLAB/Octave you should add the whole directory to the MATLAB path
     with the `addpath` command 
   * The R package should be loaded each time before using it by typing

        `library(Rslra);`

     in the R console.
    
### Special instructions for Windows MEX binary file compilation

The target `matlab-win` uses the MinGW64 compiler. MinGW64 can be installed 
from <http://mingw-w64.sourceforge.net/> or 
from <http://cran.r-project.org/bin/windows/Rtools/>.
MSYS can be installed  from <http://www.mingw.org/wiki/MSYS>. 

Prior to running `make matlab-win` you should
   * compile GSL from source using the MSYS console
   * copy the `gsl` subfolder (with header files) to one of the include paths
     of the MinGW compiler
   * copy the `libgsl.a` and `libgslcblas.a` files to `mex` directory
   * set up the paths to MinGW and MATLAB in `makefile`
   

### Notes for advanced users

The package contains a demo C++ program, which gives an example of using 
C++ interface and tests various SLRA problems (`make testc`). 
Documentation for C++ interface can be obtained by running Doxygen
in `cpp` subdirectory.

Advanced compilation options can be found in other targets of makefile, but
not all of targets may run on your machine "as is".

If you wish to try the SLICOT library, download it and copy to SLICOT
subdirectory. You can copy only a few needed files (see `makefile`). 
Use  `make xxx-slicot-xxx`.

Depending on the version of MATLAB, static binding to non-default 
BLAS and LAPACK (for example, ATLAS) can be faster. Use `xxx-static` target as
a base for your compilation instructions.
