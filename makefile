# makefile: SLRA makefile

CC = gcc -g -fPIC -static
CCPP = g++ -g -fPIC 
F77 = gcc -g -fPIC -static 

OCTAVE_MEX = mkoctfile --mex -v -DBUILD_MEX_OCTAVE #mex -v -compatibleArrayDims
MEX = mex -v -compatibleArrayDims -DBUILD_MEX_MATLAB 

SLICOT_SRC_FILES = SLICOT/MA02FD.f SLICOT/MB02CU.f SLICOT/MB02CV.f SLICOT/MB02GD.f SLICOT/MB02MD.f
SLICOT_OBJ_FILES = MA02FD.o MB02CU.o MB02CV.o MB02GD.o MB02MD.o
SLRA_SRC_FILES = slra/slra.c slra/slra_common.c slra/slra_func.c # slra/slra_func_old.c
SLRA_OBJ_FILES = slra.o slra_common.o slra_func.o  #slra_func_old.o
SLRA_INCLUDE_DIR = slra
SLRA_INCLUDE_FILES = slra/slra.h
MEX_SRC_FILES = mex/mex_slra.c

INC_FLAGS = -I./$(SLRA_INCLUDE_DIR) # -I/home/kdu/local/include -g
OPT_FLAGS = -O -pg # gprof ./test gmon.out > gmon.txt

BUILD_MODE=BUILD_DEFAULT

# Building octave mex-file
mexoct: SLICOT.a $(MEX_SRC_FILES)
	$(OCTAVE_MEX)  $(INC_FLAGS) $(MEX_SRC_FILES) $(SLRA_SRC_FILES) SLICOT.a -lgsl -lgslcblas \
	-llapack -lblas  -o slra.mex
#	cp -f slra.mex test_m

mex-im-desktop: SLICOT.a $(MEX_SRC_FILES)
	$(MEX) $(INC_FLAGS) $(MEX_SRC_FILES) $(SLRA_SRC_FILES) SLICOT.a /usr/lib/libgsl.a /usr/lib/libcblas.a \
	 /usr/lib/atlas-base/atlas/liblapack.a /usr/lib/atlas-base/atlas/libblas.a -o slra

mex-im-laptop: BUILD_MODE=MEX_MATLAB
mex-im-laptop: slra.a SLICOT.a $(MEX_SRC_FILES)
	$(MEX) $(INC_FLAGS) $(MEX_SRC_FILES) slra.a SLICOT.a -lgsl -lcblas -llapack -lblas -lgfortran -o slra

R: BUILD_MODE=BUILD_R_PACKAGE
R: 
	cp $(SLRA_INCLUDE_FILES) Rslra/src/slra
	cp $(SLICOT_SRC_FILES) Rslra/src/SLICOT
	cp $(SLRA_SRC_FILES) Rslra/src/slra
	R CMD check Rslra
	R CMD build Rslra
	R CMD INSTALL Rslra

testc: test.o slra.a SLICOT.a
	$(CCPP)  $(INC_FLAGS) $(OPT_FLAGS) -o test_c/test test.o slra.a SLICOT.a \
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

testc-im-desktop : test.o slra.a SLICOT.a
	$(CCPP)  $(INC_FLAGS) $(OPT_FLAGS) -o test_c/test test.o slra.a SLICOT.a \
	 -lm   -latlas -lgfortran -lf77blas -llapack 	-lcblas -lgslcblas -lgsl

test.o : test_c/test.cpp $(SLRA_INCLUDE_FILES) 
	$(CCPP) -D$(BUILD_MODE) $(INC_FLAGS) $(OPT_FLAGS) -c test_c/test.cpp

slra.a : $(SLRA_SRC_FILES) $(SLRA_INCLUDE_FILES)
	$(CC) -D$(BUILD_MODE) $(INC_FLAGS) $(OPT_FLAGS) -c $(SLRA_SRC_FILES)
	ar -r slra.a $(SLRA_OBJ_FILES)

SLICOT.a : $(SLICOT_SRC_FILES)
	$(F77) $(INC_FLAGS) $(OPT_FLAGS) -c $(SLICOT_SRC_FILES) 
	ar -r SLICOT.a $(SLICOT_OBJ_FILES)

slra_.pdf: slra_m.nw
	noweave -delay -index slra_m.nw >  slra_m.tex
	latex  slra_m
	bibtex slra_m
	latex  slra_m
	dvips  slra_m.dvi
	psnup -2 -d slra_m.ps slra_m-2x1.ps
	ps2pdf -sPAPERSIZE=a4  slra_m.ps
#	ps2pdf -sPAPERSIZE=a4  slra_m-2x1.ps

slra_.m:   slra_m.nw
	notangle -R"[[slra.m]]" slra_m.nw > slra_.m

clean : 
	rm *.o SLICOT.a slra.a

