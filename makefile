# makefile: STLS makefile

CC  = gcc  -g -fPIC -static
CCPP  = g++  -g -fPIC 
F77 = gcc -g -fPIC -static 

OCTAVE_MEX = mkoctfile --mex -v -DBUILD_MEX_OCTAVE #mex  -v -compatibleArrayDims
MEX = mex -v -compatibleArrayDims -DBUILD_MEX_MATLAB

SLICOT_SRC_FILES = SLICOT/MA02FD.f  SLICOT/MB02CU.f  SLICOT/MB02CV.f  SLICOT/MB02GD.f SLICOT/MB02MD.f
SLICOT_OBJ_FILES = MA02FD.o  MB02CU.o  MB02CV.o  MB02GD.o MB02MD.o
STLS_SRC_FILES = stls/stls.c  stls/mgsl.c  stls/stls_func.c # stls/stls_func_old.c
STLS_OBJ_FILES = stls.o  mgsl.o stls_func.o  #stls_func_old.o
STLS_INCLUDE_DIR = stls
STLS_INCLUDE_FILES = stls/stls.h
MEX_SRC_FILES = mex/mex_stls.c

INC_FLAGS =  -I./$(STLS_INCLUDE_DIR) # -I/home/kdu/local/include -g
OPT_FLAGS = -O -pg # gprof ./test gmon.out > gmon.txt

BUILD_MODE=BUILD_DEFAULT


# Building octave mex-file
mexoct:   SLICOT.a $(MEX_SRC_FILES)
	$(OCTAVE_MEX)  $(INC_FLAGS) $(MEX_SRC_FILES) $(STLS_SRC_FILES) SLICOT.a -lgsl -lgslcblas \
	-llapack -lblas  -o slra.mex
#	cp -f slra.mex test_m

mex-im-desktop : SLICOT.a $(MEX_SRC_FILES)
	$(MEX) $(INC_FLAGS) $(MEX_SRC_FILES) $(STLS_SRC_FILES) SLICOT.a /usr/lib/libgsl.a /usr/lib/libcblas.a \
	 /usr/lib/atlas-base/atlas/liblapack.a /usr/lib/atlas-base/atlas/libblas.a -lgfortran -o slra
#	cp -f slra.mex* test_m/

mex-im-laptop : BUILD_MODE=MEX_MATLAB
mex-im-laptop : stls.a SLICOT.a $(MEX_SRC_FILES)
	$(MEX) $(INC_FLAGS) $(MEX_SRC_FILES) stls.a SLICOT.a -lgsl -lcblas -llapack -lblas -lgfortran -o slra
#	cp -f slra.mex* test_m/

R: BUILD_MODE=BUILD_R_PACKAGE
R: 
	cp $(STLS_INCLUDE_FILES) rstls/src/stls
	cp $(SLICOT_SRC_FILES) rstls/src/SLICOT
	cp $(STLS_SRC_FILES) rstls/src/stls
	R CMD check rstls
	R CMD build rstls
	R CMD INSTALL rstls

testc : test.o stls.a SLICOT.a
	$(CCPP)  $(INC_FLAGS) $(OPT_FLAGS) -o test_c/test test.o stls.a SLICOT.a \
	/home/kdu/local/lib/liblapack.a \
	/home/kdu/local/lib/libgsl.a \
	/home/kdu/local/lib/libcblas.a \
	/home/kdu/local/lib/libf77blas.a \
	/home/kdu/local/lib/libatlas.a \
	-lgfortran -lm 

#	/home/kdu/src/lapack-3.2.1/lapack_LINUX.a \
#	/home/kdu/src/gsl-1.15/.libs/libgsl.a \
#	/home/kdu/src/CBLAS/lib/cblas_LINUX.a  \
#	/home/kdu/src/lapack-3.2.1/blas_LINUX.a \
#	-lgfortran -lm -lgfortran 
#	/home/kdu/src/gsl-1.15/cblas/.libs/libgslcblas.a \

#testc : test.o stls.a SLICOT.a
#	$(CCPP)  $(INC_FLAGS) $(OPT_FLAGS) -o test_c/test test.o stls.a SLICOT.a \
#	 -lm   -latlas -lgfortran -lf77blas -llapack 	-lcblas -lgslcblas -lgsl


test.o : test_c/test.cpp $(STLS_INCLUDE_FILES) 
	$(CCPP) -D$(BUILD_MODE) $(INC_FLAGS) $(OPT_FLAGS) -c test_c/test.cpp

stls.a : $(STLS_SRC_FILES) $(STLS_INCLUDE_FILES)
	$(CC) -D$(BUILD_MODE) $(INC_FLAGS) $(OPT_FLAGS) -c $(STLS_SRC_FILES)
	ar -r stls.a $(STLS_OBJ_FILES)

SLICOT.a : $(SLICOT_SRC_FILES)
	$(F77) $(INC_FLAGS) $(OPT_FLAGS) -c $(SLICOT_SRC_FILES) 
	ar -r SLICOT.a $(SLICOT_OBJ_FILES)

clean : 
	rm *.o SLICOT.a stls.a

