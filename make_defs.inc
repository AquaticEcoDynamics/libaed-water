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
## 2026-07: left commented on purpose. GLM/Makefile DOES promote MDEBUG to DEBUG;
## libaed does not. Uncommenting would change what a default debug build emits, so
## it is only noted here. See the ifeq ($(MDEBUG),true) block further down.
#ifeq ($(MDEBUG),true)
#  DEBUG=true
#endif

#INCLUDES+=-I${incdir}
## 2026-07 fix: was MDBG_FFLAGS="" -- make passes the two literal quote characters
## through as an empty shell argument, and ifort/flang then error with
## "no such file or directory: ''". A bare assignment is genuinely empty.
MDBG_FFLAGS=

ifeq ($(OSTYPE),Darwin)
  SHARED=-dynamiclib -undefined dynamic_lookup
  so_ext=dylib
else ifeq ($(OSTYPE),Msys)
  SHARED=-shared -Wl,-soname,$(OUTLIB).dll
  so_ext=dll
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
    FFLAGS=-fpp -warn=all -module=${moddir} -static-intel -stand=f23 -warn=nounused
    ifeq ($(WITH_CHECKS),true)
      FFLAGS+=-check=all -check=noarg_temp_created
    endif
    FFLAGS+=-real-size=64
  else
    FFLAGS=-fpp -warn all -module ${moddir} -static-intel -stand f23 -warn nounused
    ifeq ($(WITH_CHECKS),true)
      FFLAGS+=-check all -check noarg_temp_created
    endif
    FFLAGS+=-real-size 64
  endif
else ifeq ($(F90),flang)
  ifeq ($(OSTYPE),FreeBSD)
    INCLUDES+=-I../ancillary/mod
  endif
  DEBUG_FFLAGS=-g
  OMPFLAG=-fopenmp
  OPT_FFLAGS=-O3
  FFLAGS=-module-dir ${moddir}
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-Mbounds
  endif
  FFLAGS+=-fdefault-real-8 -fdefault-double-8
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
  endif
  # 2026-08-23: -fno-range-check was inside the `else` above, so WITH_CHECKS=true removed it and
  # the tree stopped compiling at aed_core.F90:217, `nan_ = zero_ / zero_` — a deliberate
  # compile-time 0/0 that gfortran rejects without this flag. The two settings are unrelated
  # (verified: `gfortran -fcheck=all -fno-range-check` accepts it), so the flag is unconditional
  # now. Before this, WITH_CHECKS=true could not build at all, which is why no bounds-checked run
  # had ever been made.
  FFLAGS+=-fno-range-check -Wno-integer-division
  FFLAGS+=-fdefault-real-8 -fdefault-double-8
endif
FFLAGS+=$(DEFINES) $(INCLUDES)

FFLAGS+=-fPIC
ifeq ($(EXTERNAL_LIBS),shared)
# FFLAGS+=-fPIC
  TARGET = ${libdir}/$(OUTLIB).${so_ext}
else
# FFLAGS+=-fPIE
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

## 2026-07 fix: this was written "ifeq ($MDEBUG),true)" - "$M" expanded to nothing
## and the parens were unbalanced, so the condition never matched and MDBG_FFLAGS
## has never been applied to ANY AED Fortran source. Corrected to proper syntax.
##
## NOTE: on the MinGW-w64 toolchain used for the Windows build (gcc 16.1.0) there
## are NO sanitizer runtimes shipped - linking with -fsanitize=address fails with
## "cannot find -lasan". So MDEBUG=true is only useful on a Linux build. Now that
## the conditional actually works, setting MDEBUG=true on Windows will break the
## link rather than being silently ignored as before.
##
## See also the commented-out "MDEBUG implies DEBUG" promotion near the top of this
## file (GLM/Makefile does have it). It is left commented deliberately - enabling it
## would change what a default debug build produces.
ifeq ($(MDEBUG),true)
  FFLAGS+=$(MDBG_FFLAGS)
endif

