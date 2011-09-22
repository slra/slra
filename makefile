# makefile: STLS makefile

CC  = gcc -fPIC 
CCPP  = g++ -fPIC 
F77 = gcc -fPIC 
OCTAVE_MEX = mkoctfile --mex -v -DMEX_OCTAVE #mex  -v -compatibleArrayDims

SLICOT_SRC_FILES = SLICOT/MA02FD.f  SLICOT/MB02CU.f  SLICOT/MB02CV.f  SLICOT/MB02GD.f SLICOT/MB02MD.f
SLICOT_OBJ_FILES = MA02FD.o  MB02CU.o  MB02CV.o  MB02GD.o MB02MD.o
STLS_SRC_FILES = stls/stls.c stls/mgsl.c
STLS_OBJ_FILES = stls.o mgsl.o
STLS_INCLUDE_DIR = stls
STLS_INCLUDE_FILES = stls/stls.h
MEX_SRC_FILES = mex/mex_stls.c


INC_FLAGS = -I/home/kdu/local/include -I./$(STLS_INCLUDE_DIR) # -g
OPT_FLAGS = -O #-pg # gprof ./test gmon.out > gmon.txt

# Building octave mex-file
mexoct :   stls.a SLICOT.a $(MEX_SRC_FILES)
	$(OCTAVE_MEX)  $(INC_FLAGS) $(MEX_SRC_FILES) stls.a SLICOT.a -lgsl -lgslcblas \
	-llapack -lblas  -o stls.mex
#	cp -f mex_stls.mex ../stls.mex
	cp -f stls.mex test_m

R: 
	cp $(STLS_INCLUDE_FILES) rstls/src/stls
	cp $(SLICOT_SRC_FILES) rstls/src/SLICOT
	cp $(STLS_SRC_FILES) rstls/src/stls
	R CMD check rstls
	R CMD build rstls
	R CMD INSTALL stls
	
testc : test.o stls.o SLICOT.a
	$(CCPP)  $(INC_FLAGS) $(OPT_FLAGS) -o test_c/test test.o stls.o SLICOT.a \
	-lgfortran -lgsl -lgslcblas -lcblas -lm -lgfortran -llapack -lf77blas -latlas	

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


