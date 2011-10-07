# makefile: STLS makefile

CC  = gcc  -g -fPIC -static
CCPP  = g++  -g -fPIC -static 
F77 = gcc -g -fPIC -static
OCTAVE_MEX = mkoctfile --mex -v -DMEX_OCTAVE #mex  -v -compatibleArrayDims

SLICOT_SRC_FILES = SLICOT/MA02FD.f  SLICOT/MB02CU.f  SLICOT/MB02CV.f  SLICOT/MB02GD.f SLICOT/MB02MD.f
SLICOT_OBJ_FILES = MA02FD.o  MB02CU.o  MB02CV.o  MB02GD.o MB02MD.o
STLS_SRC_FILES = stls/stls.c  stls/mgsl.c  stls/stls_func.c # stls/stls_func_old.c
STLS_OBJ_FILES = stls.o  mgsl.o stls_func.o  #stls_func_old.o
STLS_INCLUDE_DIR = stls
STLS_INCLUDE_FILES = stls/stls.h
MEX_SRC_FILES = mex/mex_stls.c


INC_FLAGS =  -I./$(STLS_INCLUDE_DIR) # -I/home/kdu/local/include -g
OPT_FLAGS = -O -pg # gprof ./test gmon.out > gmon.txt

# Building octave mex-file
mexoct :   stls.a SLICOT.a $(MEX_SRC_FILES)
	$(OCTAVE_MEX)  $(INC_FLAGS) $(MEX_SRC_FILES) stls.a SLICOT.a -lgsl -lgslcblas \
	-llapack -lblas  -o stls.mex
	cp -f stls.mex test_m
#	cp -f mex_stls.mex ../stls.mex

R: 
	cp $(STLS_INCLUDE_FILES) rstls/src/stls
	cp $(SLICOT_SRC_FILES) rstls/src/SLICOT
	cp $(STLS_SRC_FILES) rstls/src/stls
	R CMD check rstls
	R CMD build rstls
	R CMD INSTALL stls
	
testc : test.o stls.a SLICOT.a
	$(CCPP)  $(INC_FLAGS) $(OPT_FLAGS) -o test_c/test test.o stls.a SLICOT.a \
	/home/kdu/src/lapack-3.2.1/lapack_LINUX.a \
	/home/kdu/src/gsl-1.15/.libs/libgsl.a \
	/home/kdu/src/gsl-1.15/cblas/.libs/libgslcblas.a \
	/home/kdu/src/CBLAS/lib/cblas_LINUX.a  \
	/home/kdu/src/lapack-3.2.1/blas_LINUX.a \
	-lgfortran -lm -lgfortran 

#testc : test.o stls.a SLICOT.a
#	$(CCPP)  $(INC_FLAGS) $(OPT_FLAGS) -o test_c/test test.o stls.a SLICOT.a \
#	 -lm   -latlas -lgfortran -lf77blas -llapack 	-lcblas -lgslcblas -lgsl


test.o : test_c/test.cpp $(STLS_INCLUDE_FILES) 
	$(CCPP) $(INC_FLAGS) $(OPT_FLAGS) -c test_c/test.cpp

stls.a : $(STLS_SRC_FILES) $(STLS_INCLUDE_FILES)
	$(CC) $(INC_FLAGS) $(OPT_FLAGS) -c $(STLS_SRC_FILES)
	ar -r stls.a $(STLS_OBJ_FILES)

SLICOT.a : $(SLICOT_SRC_FILES)
	$(F77) $(INC_FLAGS) $(OPT_FLAGS) -c $(SLICOT_SRC_FILES) 
	ar -r SLICOT.a $(SLICOT_OBJ_FILES)

clean : 
	rm *.o SLICOT.a stls.a

 
