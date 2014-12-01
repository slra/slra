# makefile: SLRA makefile
CC  = gcc  -g -fPIC -static -Wno-write-strings
CCPP  = g++  -g -fPIC -Wno-write-strings
F77 = gcc -g -fPIC -static 
INC_FLAGS = -I./$(SLRA_CPP_DIR) 
OPT_FLAGS = -O2 # -pg 

OCTAVE_MEX = mkoctfile --mex -v -DBUILD_MEX_OCTAVE 
MEX = mex -v -largeArrayDims 
MEXWIN = mex.bat -v -largeArrayDims -DBUILD_MEX_WINDOWS

SLRA_OBJ_FILES=$(shell cat SLRAOBJ.txt)
SLRA_SRC_FILES=$(SLRA_OBJ_FILES:%o=%cpp)

SLRA_CPP_DIR = cpp
MEX_SRC_FILES = mex/slra_mex_obj.cpp mex/SLRAObject.cpp

BUILD_MODE=BUILD_DEFAULT

# Windows-specific defines (for compilation with Mingw64)
MINGWPATH=c:/Rtools/gcc-4.6.3
MATLABROOT=c:/Progra~1/MATLAB/R2013a
WINLIBS= -L$(MATLABROOT)/bin/win64 -L$(MATLABROOT)/extern/lib/win64/microsoft  -lmex -lmx -lmwlapack -lmwblas -leng 
WINGSL=mex/libgsl.a mex/libgslcblas.a
WINMEX=$(MINGWPATH)/bin/g++ -m64 -shared mex/mex.def -DMATLAB_MEX_FILE -DBUILD_MEX_WINDOWS -I$(MATLABROOT)/extern/include -Wno-write-strings  

# Main targets
matlab: clean $(MEX_SRC_FILES) 
	$(MEX) $(INC_FLAGS) $(MEX_SRC_FILES) $(SLRA_SRC_FILES) \
	-lgsl -lgslcblas -lmwlapack -lmwblas -o slra_mex_obj 

matlab-win: $(MEX_SRC_FILES) 
	$(WINMEX) $(INC_FLAGS) $(MEX_SRC_FILES) $(SLRA_SRC_FILES) $(WINGSL) \
	$(WINLIBS) -o slra_mex_obj.mexw64

octave: clean $(MEX_SRC_FILES)
	$(OCTAVE_MEX)  $(INC_FLAGS) $(MEX_SRC_FILES) $(SLRA_SRC_FILES) \
	-lgsl -lgslcblas -o slra_mex_obj.mex

R: BUILD_MODE=BUILD_R_PACKAGE
R: 
	mkdir -p Rslra/src/$(SLRA_CPP_DIR)
	rm -r -f Rslra/*/*.o Rslra/*/*.so
	cp $(SLRA_CPP_DIR)/*.cpp Rslra/src/$(SLRA_CPP_DIR)
	cp $(SLRA_CPP_DIR)/*.h Rslra/src/$(SLRA_CPP_DIR)
	cp SLRAOBJ.txt Rslra/src/SLRAOBJ.txt
	R CMD check Rslra
	R CMD build Rslra
	R CMD INSTALL Rslra

## Targets for advanced users
testc : clean test_c/test.o $(SLRA_OBJ_FILES) 
	$(CCPP)  $(INC_FLAGS) $(OPT_FLAGS) -o test_c/test test_c/test.o \
	$(SLRA_OBJ_FILES) -lgsl -lcblas -llapack -latlas -lblas -lm -lrt

# Targets with SLICOT
SLICOT_SRC_FILES = SLICOT/MA02FD.f  SLICOT/MB02CU.f  \
		SLICOT/MB02CV.f  SLICOT/MB02GD.f 
SLICOT_OBJ_FILES = MA02FD.o  MB02CU.o  MB02CV.o  MB02GD.o 

octave-slicot: SLICOT.a $(MEX_SRC_FILES)
	$(OCTAVE_MEX) $(INC_FLAGS) -DUSE_SLICOT $(MEX_SRC_FILES) \
	$(SLRA_SRC_FILES) -lgsl -lgslcblas -llapack -lblas -o slra.mex

matlab-slicot: BUILD_MODE=MEX_MATLAB
matlab-slicot: $(MEX_SRC_FILES) SLICOT.a
	$(MEX) $(INC_FLAGS) $(MEX_SRC_FILES) $(SLRA_SRC_FILES) SLICOT.a \
	-lgsl -lgslcblas -llapack -lblas -o slra_mex_obj

testc-slicot : BUILD_MODE=USE_SLICOT
testc-slicot : test_c/test.o $(SLRA_OBJ_FILES) SLICOT.a
	$(CCPP) $(INC_FLAGS) $(OPT_FLAGS) -DUSE_SLICOT -o test_c/test test_c/test.o \
	$(SLRA_OBJ_FILES) SLICOT.a -lgsl -lcblas -llapack \
	-latlas -lblas -lm -lgfortran -lrt

# Advanced targets for static compilation
mex-static : SLICOT.a $(MEX_SRC_FILES)
	$(MEX) $(INC_FLAGS) $(MEX_SRC_FILES) $(SLRA_SRC_FILES) \
	/usr/lib/libgsl.a  /usr/lib/atlas-base/libcblas.a  \
	/usr/lib/atlas-base/atlas/liblapack.a \
	/usr/lib/atlas-base/atlas/libblas.a -lgfortran -o slra_mex_obj 

mex-slicot-static : SLICOT.a $(MEX_SRC_FILES)
	$(MEX) $(INC_FLAGS) -DUSE_SLICOT $(MEX_SRC_FILES) $(SLRA_SRC_FILES) \
	SLICOT.a /usr/lib/libgsl.a /usr/lib/atlas-base/libcblas.a \
	/usr/lib/atlas-base/atlas/liblapack.a \
	/usr/lib/atlas-base/atlas/libblas.a -lgfortran -o slra_mex_obj 

testcomp : clean test_comp/test.o $(SLRA_OBJ_FILES) 
	$(CCPP)  $(INC_FLAGS) $(OPT_FLAGS) -o test_comp/test test_comp/test.o \
	$(SLRA_OBJ_FILES) -lgsl -lcblas -llapack -latlas -lblas -lm


# Targets for parts of the package
%.o : %.cpp
	$(CCPP) -D$(BUILD_MODE) $(INC_FLAGS) $(OPT_FLAGS) -c $< -o $@

SLICOT.a : $(SLICOT_SRC_FILES)
	$(F77) $(INC_FLAGS) $(OPT_FLAGS) -c $(SLICOT_SRC_FILES) 
	ar -r SLICOT.a $(SLICOT_OBJ_FILES)

clean : 
	rm -f -r */*.o *.o *.a 

