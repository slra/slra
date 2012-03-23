# makefile: SLRA makefile
CC  = gcc  -g -fPIC -static -Wno-write-strings
CCPP  = g++  -g -fPIC -Wno-write-strings
F77 = gcc -g -fPIC -static 

OCTAVE_MEX = mkoctfile --mex -v -DBUILD_MEX_OCTAVE 
MEX = mex -v -compatibleArrayDims -DBUILD_MEX_MATLAB 

SLICOT_SRC_FILES = SLICOT/MA02FD.f  SLICOT/MB02CU.f  SLICOT/MB02CV.f  SLICOT/MB02GD.f 
SLICOT_OBJ_FILES = MA02FD.o  MB02CU.o  MB02CV.o  MB02GD.o 
SLRA_SRC_FILES = slra/slra.cpp  slra/slra_common.cpp  slra/slra_func.cpp  slra/slra_computation.cpp \
		slra/slra_flex_structure.cpp slra/slra_flex_derivative.cpp slra/slra_flex_cholesky.cpp \
		slra/slra_optimize.cpp  
SLRA_OBJ_FILES = slra.o  slra_common.o slra_func.o slra_computation.o slra_flex_structure.o \
		slra_flex_derivative.o slra_flex_cholesky.o slra_optimize.o 
SLRA_INCLUDE_DIR = slra
SLRA_INCLUDE_FILES = slra/slra.h
MEX_SRC_FILES = mexslra/mex_slra.cpp

INC_FLAGS = -I./$(SLRA_INCLUDE_DIR) 
OPT_FLAGS = -O # -pg 

BUILD_MODE=BUILD_DEFAULT

# Main targets
testc : test.o slra.a 
	$(CCPP)  $(INC_FLAGS) $(OPT_FLAGS) -o test_c/test test.o slra.a   \
	-lgsl -lcblas -llapack -latlas -lblas -lm

mex: $(MEX_SRC_FILES)
	$(MEX) $(INC_FLAGS) $(MEX_SRC_FILES) $(SLRA_SRC_FILES) \
	-lgsl -lgslcblas -llapack -lblas -o slra

mexoct: SLICOT.a $(MEX_SRC_FILES)
	$(OCTAVE_MEX)  $(INC_FLAGS) $(MEX_SRC_FILES) $(SLRA_SRC_FILES) -lgsl -lgslcblas \
	-llapack -lblas  -o slra.mex

# Not yet tested
R: BUILD_MODE=BUILD_R_PACKAGE
R: 
	cp $(SLRA_INCLUDE_FILES) Rslra/src/slra
	cp $(SLICOT_SRC_FILES) Rslra/src/SLICOT
	cp $(SLRA_SRC_FILES) Rslra/src/slra
	R CMD check Rslra
	R CMD build Rslra
	R CMD INSTALL Rslra

# Experimental targets with SLICOT
mexoct-slicot: SLICOT.a $(MEX_SRC_FILES)
	$(OCTAVE_MEX) $(INC_FLAGS) -DUSE_SLICOT $(MEX_SRC_FILES) $(SLRA_SRC_FILES) -lgsl -lgslcblas \
	-llapack -lblas  -o slra.mex

mex-slicot-static : SLICOT.a $(MEX_SRC_FILES)
	$(MEX) $(INC_FLAGS) -DUSE_SLICOT $(MEX_SRC_FILES) $(SLRA_SRC_FILES) SLICOT.a \
	/usr/lib/libgsl.a /usr/lib/atlas-base/libcblas.a \
	/usr/lib/atlas-base/atlas/liblapack.a /usr/lib/atlas-base/atlas/libblas.a -lgortran -o slra 

mex-slicot-dynamic: BUILD_MODE=MEX_MATLAB
mex-slicot-dynamic: slra.a SLICOT.a $(MEX_SRC_FILES)
	$(MEX) $(INC_FLAGS) -DUSE_SLICOT $(MEX_SRC_FILES) slra.a SLICOT.a \
	-lgsl -lcblas -llapack -lblas -lgfortran -o mex_slra
	
testc-slicot : test.o slra.a SLICOT.a
	$(CCPP)  $(INC_FLAGS) $(OPT_FLAGS) -DUSE_SLICOT -o test_c/test test.o slra.a SLICOT.a  \
	-lgsl -lcblas -llapack -latlas -lblas -lm -lgfortran

test.o : test_c/test.cpp $(SLRA_INCLUDE_FILES) 
	$(CCPP) -D$(BUILD_MODE) $(INC_FLAGS) $(OPT_FLAGS) -c test_c/test.cpp

slra.a : $(SLRA_SRC_FILES) $(SLRA_INCLUDE_FILES)
	$(CC) -D$(BUILD_MODE) $(INC_FLAGS) $(OPT_FLAGS) -c $(SLRA_SRC_FILES)
	ar -r slra.a $(SLRA_OBJ_FILES)

SLICOT.a : $(SLICOT_SRC_FILES)
	$(F77) $(INC_FLAGS) $(OPT_FLAGS) -c $(SLICOT_SRC_FILES) 
	ar -r SLICOT.a $(SLICOT_OBJ_FILES)

clean : 
	rm *.o SLICOT.a slra.a *.ps 

