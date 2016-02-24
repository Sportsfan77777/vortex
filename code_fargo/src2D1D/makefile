SHELL         =  /bin/sh

MAINOBJ          =  LowTasks.o SideEuler.o Output.o Init2D1D.o main.o Theo.o\
		 Interpret.o SourceEuler2D1D.o TransportEuler.o\
		 Planet.o RungeKunta.o Viscosity.o Psys.o Force.o\
		 var.o Gframeforce.o split.o merge.o commbound.o\
                 fpe.o rebin.o Output1D.o Pframeforce1D.o\
		 rebin1D.o SideEuler1D.o Theo1D.c\
		 TransportEuler1D.o Viscosity1D.o

DUMMY		= mpi_dummy.o

LIBS            = -lm
#OPT             = -O3 -mcpu=pentiumpro -tpp6 -xK -ipo -unroll -openmp
OPT		= 
#Below: setup for i686/Linux (needs lam-mpi package), with GNU Compiler
#COMP            = icc   --> C'est la VO. Mais y a pas icc sur debussy.
COMP            = cc
#OPT		= -O4 -march=pentium3 -Wall -ffast-math
OPT		=
#OPTSEQ		= -O3 -mcpu=pentiumpro -tpp6 -xK -ipo -unroll -w2 -Wall -wd1418,310,266,810,1419,981
OPTSEQ		= -O3
#-D_TRAP_FPE
PARACOMP	= mpicc
#Below: setup for i686/Linux (needs lam-mpi package), with Intel Compiler, sequential only
#COMP            = icc
#OPT		= -tpp7 -axW -ipo -O3
#OPT		= -g
#Below: setup for Compaq HP/SC (Nickel/Chrome):
#COMP            = cc
#PARACOMP	= cc
#OPT		= -arch ev6 -fast -O3 -inline all -mp
PARALIBS	= -lmpi
PARAOPT		= -D_PARALLEL -O3
AUTOINCL        = param.h param_noex.h global_ex.h

include	.config

EXENAME        = ../fargo2D1D

SRC = *.c
INCLUDE = *.h

ifeq ($(BUILD),parallel)
COMPILER	= $(PARACOMP)
LIBRARIES	= $(LIBS) $(PARALIBS)
OPTIONS		= $(OPT) $(PARAOPT)
OBJ		= $(MAINOBJ)
else
COMPILER	= $(COMP)
LIBRARIES	= $(LIBS)
OPTIONS		= $(OPT) $(OPTSEQ)
OBJ		= $(MAINOBJ) $(DUMMY)
endif

all: conditionalrebuild $(AUTOINCL) $(OBJ) $(EXENAME) archive
	@echo "" 
	@echo ""
	@echo "      NOTE"
	@echo ""
ifeq ($(BUILD),parallel)
	@echo "This built is PARALLEL (MPI)."
else
	@echo "This built is SEQUENTIAL."
endif
	@echo "If you want to change this,"
	@echo "then you need to issue:"
	@echo ""
ifeq ($(BUILD),parallel)
	@echo "gmake BUILD=sequential"
else
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
	@echo "Creating ../.source.tar.bz2"
	@tar cf ../.source.tar *.c
	@tar rf ../.source.tar *.h
	@tar rf ../.source.tar makefile
	@tar rf ../.source.tar varparser.pl
	@gzip -9 -f ../.source.tar

para:
	@gmake BUILD=parallel

seq:
	@gmake BUILD=sequential

$(AUTOINCL) : var.c global.h makefile varparser.pl
	@./varparser.pl

$(OBJ): mp.h fondam.h param.h param_noex.h types.h makefile

.PHONY: clean mrproper
clean:
	rm -f *.o *~ *.s *.il global_ex.h param.h param_noex.h $(EXENAME)
mrproper:
	rm -f *.o *~ *.s *.il

.c.o  :
	$(COMPILER) $*.c -c $(OPTIONS)

