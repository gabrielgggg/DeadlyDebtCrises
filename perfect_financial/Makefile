
tools = intel

nlflib = /users/PAS2282/mihalache2/nlopt-f/build/ifort_EB5B9F16B771E4F5/nlopt-f/
nlfinc = /users/PAS2282/mihalache2/nlopt-f/build/ifort_EB5B9F16B771E4F5/
NLlibs = -L$(nlflib) -L${MKLROOT}/lib -L/users/PAS2282/mihalache2/lib64/
NLcflags = -I$(nlfinc) -I${MKLROOT}/include
# NLlibs = `pkg-config nlopt-f --libs` -L${MKLROOT}/lib
# NLcflags = `pkg-config nlopt-f --cflags` -I${MKLROOT}/include

ifeq ($(tools),intel)
  compiler = ifort
  flags = -O3 -parallel -align -xHost -diag-disable 8291 -std18 -warn all,noexternal $(NLcflags)
  allow_unused = -diag-disable 7712 
  libs = -mkl -lnlopt-f -lnlopt
else ifeq ($(tools),intelLLVM)
  compiler = ifx
  flags = -Ofast -align -xHost -qopenmp -ipo -std18 -g -warn all,noexternal $(NLcflags) 
  allow_unused = -diag-disable 7712 
  libs = -qmkl
else ifeq ($(tools),gnu)
  compiler = gfortran
  flags = -O3  -pipe  -g $(NLcflags) 
  allow_unused = -Wno-unused-parameter
  libs = -llapack -lblas -lm
else ifeq ($(tools),cray)
  compiler = mpifort 
  flags = -homp -O3 -hipa5 -hscalar3 -hvector3 -h wp -h pl=./_pl -h nobounds $(NLcflags) 
  allow_unused = -Wno-unused-parameter
  libs =
else
  $(error Unknown toolchain?)
endif

.PHONY: all clean

all:	bin/perfect

bin/perfect:	src/perfect.f90 perfectMod.o
	$(compiler) $(flags) -o bin/perfect *.o src/perfect.f90 $(NLlibs) ${libs} 

datetime_module.o:	src/datetime_module.f90 
	$(compiler) $(flags) -c src/datetime_module.f90

sim.o:	src/sim.f90 
	$(compiler) $(flags) -c src/sim.f90

NL.o:	src/NL.f90
	$(compiler) $(flags) $(allow_unused) -c src/NL.f90

perfectMod.o:	src/perfectMod.f90 NL.o sim.o datetime_module.o
	$(compiler) $(flags) -c src/perfectMod.f90

clean:
	rm -rf bin/* ./*.mod ./*.o

cleanresults:
	rm -f results/*.tab results/*.bin

cleanmlab:
	rm -f results/*.mat results/*.fig

ccc:	clean cleanresults

run:	all
	DYLD_LIBRARY_PATH=${MKLROOT}/lib ./bin/perfect


