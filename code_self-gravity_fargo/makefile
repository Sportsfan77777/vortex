# !! Define Environment variable FARGO_ARCH to reflect your
#architecture for ex.: setenv FARGO_ARCH LINUX (or export
#FARGO_ARCH=LINUX depending on your shell) Possible choices :
#undefined, LINUX, INTEL, MACOSX, EV68, EVMP, SPARC

# Generic IA32 Linux platform for pentium 3 or higher (otherwise edit
#accordingly) FARGO_ARCH must be set to LINUX
CC_LINUX  = gcc
OPT_LINUX = -Wall -ffast-math
OPTSEQ_LINUX = 
PARAOPT_LINUX =
PARACC_LINUX = mpicc

# Intel Linux platform optimized for pentium 4 (otherwise edit
#accordingly) FARGO_ARCH must be set to INTEL (sequential only, MPI
#identical to generic IA32)
CC_INTEL  = icc
OPT_INTEL = 
OPTSEQ_INTEL = -O3 -tpp7 -ipo -unroll -w2 -Wall -wd1418,310,266,810,1419,981,161 
PARAOPT_INTEL = -O3 -march=pentium3 -Wall -ffast-math
PARACC_INTEL = mpicc

# Opteron platform FARGO_ARCH must be set to OPTERON
CC_OPTERON = gcc
OPT_OPTERON = -O3 -m64 -Wall -ffast-math
OPTSEQ_OPTERON =
PARAOPT_OPTERON =
PARACC_OPTERON = mpicc

# Macintosh MAC OS X platform (with gcc) FARGO_ARCH must be set to
#MACOSX
CC_MACOSX  = gcc
OPT_MACOSX = -O3
OPTSEQ_MACOSX = 
PARAOPT_MACOSX =
PARACC_MACOSX = mpicc

# Generic options valid for any platfom which has 'cc' in the path
#These are used if you leave FARGO_ARCH undefined
CC_  = mpiicc
OPT_ = -O3 -xHost
OPTSEQ_ = 
PARAOPT_ =
PARACC_ = mpiicc $(CPPFLAGS)

# Setup for Compaq HP/SC (Digital EV68) with no OpenMP support Set the
#FARGO_ARCH variable to EV68
CC_EV68  = cc
OPT_EV68 = -arch ev6 -fast -O3 -inline all -msg_disable ompdirignored
OPTSEQ_EV68 = 
PARAOPT_EV68 =
PARACC_EV68 = cc

# Setup for Compaq HP/SC (Digital EV68) with OpenMP support Set the
#FARGO_ARCH variable to EVMP
CC_EVMP  = cc
OPT_EVMP = -arch ev6 -fast -O3 -inline all -mp
OPTSEQ_EVMP = 
PARAOPT_EVMP =
PARACC_EVMP = cc

# Generic Sparc (v7) SunStation platform w/ gcc FARGO_ARCH must be set
#to SPARC
CC_SPARC  = gcc
OPT_SPARC = -O3 -Wall
OPTSEQ_SPARC = 
PARAOPT_SPARC =
PARACC_SPARC = mpicc
#
#
#
#--------------------No Changes Needed after this line------------------------
#
#
#
SHELL		=  /bin/sh

MAINOBJ         = LowTasks.o SideEuler.o Output.o Init.o main.o Theo.o\
		  Interpret.o SourceEuler.o TransportEuler.o Stockholm.o\
		  Planet.o RungeKunta.o Viscosity.o Psys.o Force.o var.o\
		  Pframeforce.o split.o merge.o commbound.o fpe.o rebin.o\
		  sgmain.o sginit.o sgdens.o sgkernel.o sgacc.o sgzero.o\
		  sgupdate.o sgvelinit.o sgsysinit.o axilib.o aniso.o\

MPIDUMMY	= mpi_dummy.o
FFTWDUMMY	= fftw_dummy.o

COMP        = $(CC_$(FARGO_ARCH))
PARACOMP    = $(PARACC_$(FARGO_ARCH))
OPT         = $(OPT_$(FARGO_ARCH))
OPTSEQ      = $(OPTSEQ_$(FARGO_ARCH))
PARAOPT     = $(PARAOPT_$(FARGO_ARCH)) -D_PARALLEL -I$(MPI_PREFIX)/include
FFTWOPT	    = -D_FFTW -I$(FFTW_PREFIX)/include
LIBS        = -lm
PARALIBS    = -L$(MPI_PREFIX)/lib -lmpi
FFTWLIBS    = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
FFTWLIBS    = -L$(FFTW_PREFIX)/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw  # The 1st arg is extra
AUTOINCL    = param.h param_noex.h global_ex.h

include	.config

EXENAME        = ../fargo

ARCHIVE		= $(EXENAME:../%=%.tar)
ARCHIVECOMP	= $(EXENAME:../%=%.tar.gz)

SRC = *.c
INCLUDE = *.h

ifeq ($(BUILD),parallelfftw)
COMPILER	= $(PARACOMP)
LIBRARIES	= $(LIBS) $(PARALIBS) $(FFTWLIBS)
OPTIONS		= $(OPT) $(PARAOPT) $(FFTWOPT)
OBJ		= $(MAINOBJ)
endif
ifeq ($(BUILD),parallel)
COMPILER	= $(PARACOMP)
LIBRARIES	= $(LIBS) $(PARALIBS)
OPTIONS		= $(OPT) $(PARAOPT)
OBJ		= $(MAINOBJ) $(FFTWDUMMY)
endif
ifeq ($(BUILD),sequentialfftw)
COMPILER	= $(COMP)
LIBRARIES	= $(LIBS) $(FFTWLIBS)
OPTIONS		= $(OPT) $(OPTSEQ) $(FFTWOPT)
OBJ		= $(MAINOBJ) $(MPIDUMMY)
endif
ifeq ($(BUILD),sequential)
COMPILER	= $(COMP)
LIBRARIES	= $(LIBS)
OPTIONS		= $(OPT) $(OPTSEQ)
OBJ		= $(MAINOBJ) $(MPIDUMMY) $(FFTWDUMMY)
endif

all: conditionalrebuild $(AUTOINCL) $(OBJ) $(EXENAME) archive
	@echo "" 
	@echo ""
	@echo "      NOTE"
	@echo ""
ifeq ($(BUILD),parallelfftw)
	@echo "This built is PARALLEL (MPI) and uses FFTW librairies"
	@echo "If you want to change this,"
	@echo "then you need to issue:"
	@echo "gmake BUILD=parallel"
	@echo "gmake BUILD=sequentialfftw"
	@echo "gmake BUILD=sequential"
endif
ifeq ($(BUILD),parallel)
	@echo "This built is PARALLEL (MPI)"
	@echo "If you want to change this,"
	@echo "then you need to issue:"
	@echo "gmake BUILD=parallelfftw"
	@echo "gmake BUILD=sequentialfftw"
	@echo "gmake BUILD=sequential"
endif
ifeq ($(BUILD),sequentialfftw)
	@echo "This built is SEQUENTIAL and uses FFTW librairies"
	@echo "If you want to change this,"
	@echo "then you need to issue:"
	@echo "gmake BUILD=sequential"
	@echo "gmake BUILD=parallelfftw"
	@echo "gmake BUILD=parallel"
endif
ifeq ($(BUILD),sequential)
	@echo "This built is SEQUENTIAL"
	@echo "If you want to change this,"
	@echo "then you need to issue:"
	@echo "gmake BUILD=sequentialfftw"
	@echo "gmake BUILD=parallelfftw"
	@echo "gmake BUILD=parallel"
endif
	@echo ""

$(EXENAME): $(OBJ)
	$(COMPILER) $(OBJ) $(OPTIONS) -o $(EXENAME) $(LIBRARIES)

.PHONY: conditionalrebuild
ifneq ($(BUILD),$(OLDBUILD))
conditionalrebuild: clean
	@echo "OLDBUILD = $(BUILD)" > .config
	@echo "BUILD = $(BUILD)" >> .config
else
conditionalrebuild:
endif

.oldconfig:
.config:

archive : $(SRC) $(INCL) makefile varparser.pl	
	@echo "Creating ../source.tar.bz2"
	@tar cf ../source.tar *.c
	@tar rf ../source.tar *.h
	@tar rf ../source.tar makefile
	@tar rf ../source.tar varparser.pl
	@bzip2 -9 -f ../source.tar

para:
	@gmake BUILD=parallel
parafftw:
	@gmake BUILD=parallelfftw
seq:
	@gmake BUILD=sequential
seqfftw:
	@gmake BUILD=sequentialfftw

$(AUTOINCL) : var.c global.h makefile varparser.pl
	@./varparser.pl

$(OBJ): mp.h fondam.h param.h param_noex.h types.h makefile

.PHONY: clean mrproper package

mrproper:
	rm -f *.o *~ *.s *.il $(AUTOINCL) $(EXENAME) ../core.*\
	*.tex *.dvi *.pdf *.ps *.log *.aux *.lint $(ARCHIVE)\
	$(ARCHIVECOMP)

clean:
	rm -f *.o *~ *.s *.il

package: $(ARCHIVECOMP)

release:
	@echo "Creating archive fargoadsg.tar.gz for release"
	@cd ../; tar -c\
         --exclude src/adimvalue.c\
         -f fargoadsg.tar src/*.c
	@cd ..; tar rf fargoadsg.tar src/makefile
	@cd ..; tar rf fargoadsg.tar src/*.pl
	@cd ..; tar rf fargoadsg.tar src/*.h
	@cd ..; tar rf fargoadsg.tar in/*.par
	@cd ..; tar rf fargoadsg.tar in/*.cfg
	@cd ..; tar rf fargoadsg.tar idl/*.pro
	@cd ..; tar rf fargoadsg.tar out1/
	@cd ..; gzip -9 fargoadsg.tar

.c.o  :
	$(COMPILER) $*.c -c $(OPTIONS)
