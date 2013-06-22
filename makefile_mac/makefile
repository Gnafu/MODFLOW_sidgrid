# Makefile example for compiling MODFLOW-2005-LGR with VSF extension:
# SID&GRID 2011
# 
#
PROGRAM = mflgr

# Define the Fortran compile flags
 
F90FLAGS= -O2 -Wall -m64
F90= gfortran

#
# Define the C compile flags
# -D_UF defines UNIX naming conventions for mixed language compilation.
#
CFLAGS= -D_UF -O2 -arch x86_64 -Wall
CC= gcc

# Define GMG objects
#
GMG = r_vector.o\
      solvers.o\
      ccfd.o\
      mf2kgmg.o\

# Define the libraries

SYSLIBS= -lc
USRLIB  =

# Define all Fortran object files which make up MODFLOW

SOURCEFILES  = \
        gmg7.f \
        mhc7.f\
        gwf2bfh1.f\
        de47.f \
        pcg7.f \
        sip7.f \
        gmg7.f \
        mhc7.f\
        gwf2uzf1.f \
        gwf2gag7.f \
        gwf2chd7.f \
        gwf2drn7.f \
        gwf2drt7.f \
        gwf2ets7.f \
        gwf2evt7.f \
        gwf2fhb7.f \
        gwf2ghb7.f \
        gwf2hfb7.f \
        gwf2ibs7.f \
        gwf2res7.f \
        gwf2riv7.f \
        gwf2str7.f \
        gwf2sub7.f \
        gwf2swt7.f\
        gwf2wel7.f \
        hufutl7.f \
        obs2bas7.f \
        obs2drn7.f \
        obs2ghb7.f \
        obs2riv7.f \
        obs2chd7.f \
        obs2str7.f \
        parutl7.f \
        gwf2mnw17.f \
        gwf2mnw27.f \
        gwf2mnw2i7.f \
        gwf2hydmod7.f \
        vsf1ref1_soilchar.f \
        vsf2lpf7.f \
        vsf2pnd1.f \
		vsf2rze1.f \
		vsf2sev1.f \
		vsf2spf1.f \
		ats2.f \
		utl7.f \
        lmt7.f \
		mflgr.f \


OBJECTS = \
        gwf2bas7.o \
        gwf2bcf7.o \
        gwf2lpf7.o \
        gwf2huf7.o \
        gwf2lgr1.o \
        gmg7.o \
        mhc7.o\
        gwf2bfh1.o\
        de47.o \
        pcg7.o \
        sip7.o \
        gmg7.o \
        mhc7.o\
        gwfuzfmodule.o \
        gwfsfrmodule.o \
        gwf2lak7.o \
        routing_module.o \
        cascading.o \
        gwf2sfr7.o \
        gwf2uzf1.o \
        gwf2gag7.o \
        gwf2chd7.o \
        gwf2drn7.o \
        gwf2drt7.o \
        gwf2ets7.o \
        gwf2evt7.o \
        gwf2fhb7.o \
        gwf2ghb7.o \
        gwf2hfb7.o \
        gwf2ibs7.o \
        gwf2rch7.o \
        gwf2res7.o \
        gwf2riv7.o \
        gwf2str7.o \
        gwf2sub7.o \
        gwf2swt7.o\
        gwf2wel7.o \
        hufutl7.o \
        obs2bas7.o \
        obs2drn7.o \
        obs2ghb7.o \
        obs2riv7.o \
        obs2chd7.o \
        obs2str7.o \
        parutl7.o \
        gwf2mnw17.o \
        gwf2mnw27.o \
        gwf2mnw2i7.o \
        gwf2hydmod7.o \
        vsf2ref1.o \
        vsf1ref1_soilchar.o \
        vsf2lpf7.o \
        vsf2pnd1.o \
		vsf2rze1.o \
		vsf2sev1.o \
		vsf2spf1.o \
		ats2.o \
		utl7.o \
        lmt7.o \
		mflgr.o \
# Define Task Function

install: mflgr

# Define Task Function

all: mflgr

# Define what mflgr does

preface:
	gfortran -c gwf2bas7.f gwf2bcf7.f gwf2lpf7.f gwf2huf7.f gwfuzfmodule.f gwfsfrmodule.f gwf2lak7.f gwf2lgr1.f  
preface2:
	gfortran -c routing_module.f95 gwf2sfr7.f gwf2rch7.f vsf2ref1.f cascading.f95
preface3: $(SOURCEFILES)
	gfortran -c $(SOURCEFILES)
preface4:
	gfortran -c -D_UF -Ilocation *.c

mflgr: preface preface2 preface3 preface4
	gfortran *.o -o mflgr_vsf

# Object codes

.f.o:
	$(F90) $(F90FLAGS) -c $<	

.c.o:
	$(CC) $(CFLAGS) -c $<

#  end
