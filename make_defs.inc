###############################################################################
#                                                                             #
# make_defs.inc - include file for making libaed.                             #
#                                                                             #
#  Developed by :                                                             #
#      AquaticEcoDynamics (AED) Group                                         #
#      School of Agriculture and Environment                                  #
#      The University of Western Australia                                    #
#                                                                             #
#      http://aquatic.science.uwa.edu.au/                                     #
#                                                                             #
#  Copyright 2013-2026 - The University of Western Australia                  #
#                                                                             #
#   AED is free software: you can redistribute it and/or modify               #
#   it under the terms of the GNU General Public License as published by      #
#   the Free Software Foundation, either version 3 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   AED is distributed in the hope that it will be useful,                    #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

ifeq ($(shell uname),Linux)
  OSTYPE=$(shell uname -s)
else ifeq ($(shell uname),Darwin)
  OSTYPE=$(shell uname -s)
else ifeq ($(shell uname),FreeBSD)
  OSTYPE=$(shell uname -s)
else
  OSTYPE=$(shell uname -o)
endif

ifeq ($(VERSION),)
  VERSION=$(shell grep AED_VERSION ../libaed-water/include/aed.h | head -1 | cut -f2 -d\")
endif
SOVERS=$(shell echo $(VERSION) | cut -f1 -d\.)
VERS=$(shell echo $(VERSION) | cut -f2- -d\.)

objdir=obj
srcdir=src
libdir=lib
moddir=mod
incdir=include

ifeq ($(F90),)
  F90=gfortran
endif
#ifeq ($(MDEBUG),true)
#  DEBUG=true
#endif

#INCLUDES+=-I${incdir}
MDBG_FFLAGS=""

ifeq ($(OSTYPE),Darwin)
  SHARED=-dynamiclib -undefined dynamic_lookup
  so_ext=dylib
else
  SHARED=-shared -Wl,-soname,$(OUTLIB).so.$(SOVERS)
  so_ext=so
endif

ifeq ($(F90),ifort)
  INCLUDES+=-I/opt/intel/include
  DEBUG_FFLAGS=-g -traceback
  OMPFLAG=-qopenmp
  OPT_FFLAGS=-O3
  FFLAGS=-warn all -module ${moddir} -static-intel -mp1 -stand f08 -warn nounused
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-check all -check noarg_temp_created
  endif
  FFLAGS+=-real-size 64
else ifeq ($(F90),ifx)
# INCLUDES+=-I/opt/intel/oneapi/compiler/latest/include
  DEBUG_FFLAGS=-g -traceback -O0
  OMPFLAG=-qopenmp
  OPT_FFLAGS=-O3
  ifeq ($(OSTYPE),Msys)
    FFLAGS=-fpp -warn=all -module=${moddir} -static-intel -mp1 -stand=f23 -warn=nounused
    ifeq ($(WITH_CHECKS),true)
      FFLAGS+=-check=all -check=noarg_temp_created
    endif
    FFLAGS+=-real-size=64
  else
    FFLAGS=-fpp -warn all -module ${moddir} -static-intel -mp1 -stand f23 -warn nounused
    ifeq ($(WITH_CHECKS),true)
      FFLAGS+=-check all -check noarg_temp_created
    endif
    FFLAGS+=-real-size 64
  endif
else ifeq ($(F90),flang)
  ifeq ($(OSTYPE),FreeBSD)
    INCLUDES+=-I../ancillary/freebsd/mod
  endif
  DEBUG_FFLAGS=-g
  OMPFLAG=-fopenmp
  OPT_FFLAGS=-O3
  FFLAGS=-module ${moddir}
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-Mbounds
  endif
  FFLAGS+=-r8
else
  DEBUG_FFLAGS=-g -fbacktrace
  MDBG_FFLAGS=-fsanitize=address
  OMPFLAG=-fopenmp
  OPT_FFLAGS=-O3
  FFLAGS=-Wall -J ${moddir} -ffree-line-length-none -std=f2018
  ifeq ($(OSTYPE),Msys)
    FFLAGS+=-D_WIN32
  endif
  FFLAGS+=-fall-intrinsics -Wno-unused -Wno-unused-dummy-argument
  FFLAGS+=-Wno-missing-include-dirs
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-fcheck=all
  else
    FFLAGS+=-fno-range-check -Wno-integer-division
  endif
  FFLAGS+=-fdefault-real-8 -fdefault-double-8
endif
FFLAGS+=$(DEFINES) $(INCLUDES)

ifeq ($(EXTERNAL_LIBS),shared)
  FFLAGS+=-fPIC
  TARGET = ${libdir}/$(OUTLIB).${so_ext}
else
  FFLAGS+=-fPIE
  TARGET = ${libdir}/$(OUTLIB).a
endif

ifeq ($(DEBUG),true)
  DEBUG_CFLAGS=-g
  OPT_CFLAGS=
  OPT_FFLAGS=
  FFLAGS+=$(DEBUG_FFLAGS)
else
  DEBUG_FFLAGS=
  DEBUG_CFLAGS=
  # OPT_CFLAGS=-O4 -Ofast -frounding-math
  OPT_CFLAGS=-O3
  # OPT_CFLAGS=
  # OPT_FFLAGS=
  FFLAGS+=$(OPT_FFLAGS)
endif

ifeq ($MDEBUG),true)
  FFLAGS+=$(MDBG_FFLAGS)
endif

