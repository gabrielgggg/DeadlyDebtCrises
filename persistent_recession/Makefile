
tools = intel

ifeq ($(tools),intel)
  compiler = mpifort
  flags = -Ofast -parallel -align -xHost -qopenmp -ipo -std18 -warn all,noexternal
  libs = -mkl
else ifeq ($(tools),cray)
  compiler = mpifort 
  flags = -O3 -h omp,scalar3,vector3,ipa3,contiguous,nobounds,wp,pl=./_pl
  libs =
else ifeq ($(tools),intelLLVM)
  compiler = mpiifort -fc=ifx
  flags = -Ofast -align -xHost -qopenmp -ipo -std18 -g -warn all,noexternal 
  libs = -qmkl
else ifeq ($(tools),gnu)
  compiler = mpifort
  flags = -Ofast -march=native -mtune=native -pipe -fopenmp -flto -Wall -pedantic -g
  libs = -llapack -lblas -lm
else ifeq ($(tools),fujitsu)
  compiler = mpifrt
  flags = -K fast,SVE,openmp,optmsg=2 -N lst=t -flto -fw
  libs = -SSL2 
else
  $(error Unknown toolchain?)
endif

.PHONY: clean cleanresults cleanmlab ccc

bin/pdef:	src/pdef.f90 pdefMod.o VFI.o MPIMod.o datetime_module.o
	$(compiler) $(flags) -o bin/pdef *.o src/pdef.f90 ${libs}

run:	bin/pdef
	mpirun --map-by ppr:1:node ./bin/pdef

Param.o:	src/Param.f90
	$(compiler) $(flags) -c src/Param.f90

sim.o:	src/sim.f90 
	$(compiler) $(flags) -c src/sim.f90

NL.o:	src/NL.f90
	$(compiler) $(flags) -c src/NL.f90

datetime_module.o:	src/datetime_module.f90 
	$(compiler) $(flags) -c src/datetime_module.f90

pdefMod.o:	src/pdefMod.f90 NL.o sim.o Param.o MPIMod.o datetime_module.o
	$(compiler) $(flags) -c src/pdefMod.f90

MPIMod.o:	src/MPIMod.f90 Param.o
	$(compiler) $(flags) -c src/MPIMod.f90

VFI.o:	src/VFI.f90 NL.o sim.o Param.o pdefMod.o MPIMod.o
	$(compiler) $(flags) -c src/VFI.f90
  
clean:
	rm -f bin/* ./*.mod ./*.o ./*.lst
	rm -rf ./_pl

cleanresults:
	rm -f results/*.tab results/*.bin

cleanmlab:
	rm -f results/*.mat results/*.fig

ccc:	clean cleanresults


