# makefile: SLRA makefile
CC  = gcc  -g -fPIC -static -Wno-write-strings
CCPP  = g++  -g -fPIC -Wno-write-strings
F77 = gcc -g -fPIC -static 
INC_FLAGS = -I./$(SLRA_INCLUDE_DIR) 
OPT_FLAGS = -O # -pg 

OCTAVE_MEX = mkoctfile --mex -v -DBUILD_MEX_OCTAVE 
MEX = mex -v -compatibleArrayDims -DBUILD_MEX_MATLAB 

SLRA_OBJ_FILES=$(shell cat SLRAOBJ.txt)
SLRA_SRC_FILES=$(SLRA_OBJ_FILES:%o=%cpp)

SLRA_INCLUDE_DIR = slra
MEX_SRC_FILES = mexslra/mex_slra.cpp

BUILD_MODE=BUILD_DEFAULT

# Main targets
mex: $(MEX_SRC_FILES) 
	$(MEX) $(INC_FLAGS) $(MEX_SRC_FILES) $(SLRA_SRC_FILES) \
	-lgsl -lgslcblas -llapack -lblas -o slra

mex-octave: $(MEX_SRC_FILES)
	$(OCTAVE_MEX)  $(INC_FLAGS) $(MEX_SRC_FILES) $(SLRA_SRC_FILES) \
	-lgsl -lgslcblas -llapack -lblas -o slra.mex

R: BUILD_MODE=BUILD_R_PACKAGE
R: 
	mkdir -p Rslra/src/slra
	cp slra/*.cpp Rslra/src/slra
	cp slra/*.h Rslra/src/slra
	cp SLRAOBJ.txt Rslra/src/SLRAOBJ.txt
	R CMD check Rslra
	R CMD build Rslra
	R CMD INSTALL Rslra

testc : test_c/test.o $(SLRA_OBJ_FILES) 
	$(CCPP)  $(INC_FLAGS) $(OPT_FLAGS) -o test_c/test test.o \
	$(SLRA_OBJ_FILES) -lgsl -lcblas -llapack -latlas -lblas -lm

# Targets with SLICOT
SLICOT_SRC_FILES = SLICOT/MA02FD.f  SLICOT/MB02CU.f  \
		SLICOT/MB02CV.f  SLICOT/MB02GD.f 
SLICOT_OBJ_FILES = MA02FD.o  MB02CU.o  MB02CV.o  MB02GD.o 

mex-octave-slicot: SLICOT.a $(MEX_SRC_FILES)
	$(OCTAVE_MEX) $(INC_FLAGS) -DUSE_SLICOT $(MEX_SRC_FILES) \
	$(SLRA_SRC_FILES) -lgsl -lgslcblas -llapack -lblas -o slra.mex

mex-slicot: BUILD_MODE=MEX_MATLAB
mex-slicot: $(MEX_SRC_FILES) SLICOT.a
	$(MEX) $(INC_FLAGS) $(MEX_SRC_FILES) $(SLRA_SRC_FILES) SLICOT.a \
	-lgsl -lgslcblas -llapack -lblas -o slra

testc-slicot : BUILD_MODE=USE_SLICOT
testc-slicot : test_c/test.o $(SLRA_OBJ_FILES) SLICOT.a
	$(CCPP) $(INC_FLAGS) $(OPT_FLAGS) -DUSE_SLICOT -o test_c/test test.o \
	$(SLRA_OBJ_FILES) SLICOT.a -lgsl -lcblas -llapack \
	-latlas -lblas -lm -lgfortran

# Advanced targets for static compilation
mex-static : SLICOT.a $(MEX_SRC_FILES)
	$(MEX) $(INC_FLAGS) $(MEX_SRC_FILES) $(SLRA_SRC_FILES) \
	/usr/lib/libgsl.a  /usr/lib/atlas-base/libcblas.a  \
	/usr/lib/atlas-base/atlas/liblapack.a \
	/usr/lib/atlas-base/atlas/libblas.a -lgfortran -o slra 

mex-slicot-static : SLICOT.a $(MEX_SRC_FILES)
	$(MEX) $(INC_FLAGS) -DUSE_SLICOT $(MEX_SRC_FILES) $(SLRA_SRC_FILES) \
	SLICOT.a /usr/lib/libgsl.a /usr/lib/atlas-base/libcblas.a \
	/usr/lib/atlas-base/atlas/liblapack.a \
	/usr/lib/atlas-base/atlas/libblas.a -lgfortran -o slra 

# Targets for parts of the package
%.o : %.cpp
	$(CCPP) -D$(BUILD_MODE) $(INC_FLAGS) $(OPT_FLAGS) -c $< -o $@

SLICOT.a : $(SLICOT_SRC_FILES)
	$(F77) $(INC_FLAGS) $(OPT_FLAGS) -c $(SLICOT_SRC_FILES) 
	ar -r SLICOT.a $(SLICOT_OBJ_FILES)

clean : 
	rm -r *.o *.mex* *.a

