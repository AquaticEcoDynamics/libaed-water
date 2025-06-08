###############################################################################
#                                                                             #
# Makefile to build libaed-water                                              #
#                                                                             #
#  Developed by :                                                             #
#      AquaticEcoDynamics (AED) Group                                         #
#      School of Agriculture and Environment                                  #
#      The University of Western Australia                                    #
#                                                                             #
#      http://aquatic.science.uwa.edu.au/                                     #
#                                                                             #
#  Copyright 2013 - 2024 - The University of Western Australia                #
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

LIBAEDWTR=aed-water
OUTLIB=lib$(LIBAEDWTR)

INCLUDES=-Iinclude
include make_defs.inc

OBJS=${objdir}/aed_core.o \
     ${objdir}/aed_util.o \
     ${objdir}/aed_bio_utils.o \
     ${objdir}/aed_zoop_utils.o \
     ${objdir}/aed_bio_particles.o \
     ${objdir}/aed_csv_reader.o \
     ${objdir}/aed_carbon.o \
     ${objdir}/aed_dummy.o \
     ${objdir}/aed_gctypes.o \
     ${objdir}/aed_gclib.o \
     ${objdir}/aed_gcsolver.o \
     ${objdir}/aed_geochemistry.o \
     ${objdir}/aed_habitat_water.o \
     ${objdir}/aed_nitrogen.o \
     ${objdir}/aed_noncohesive.o \
     ${objdir}/aed_organic_matter.o \
     ${objdir}/aed_oxygen.o \
     ${objdir}/aed_pathogens.o \
     ${objdir}/aed_pesticides.o \
     ${objdir}/aed_phosphorus.o \
     ${objdir}/aed_phytoplankton.o \
     ${objdir}/aed_phyto_abm.o \
     ${objdir}/aed_sedflux.o \
     ${objdir}/aed_silica.o \
     ${objdir}/aed_totals.o \
     ${objdir}/aed_tracer.o \
     ${objdir}/aed_zooplankton.o \
     ${objdir}/aed_water.o \
     ${objdir}/aed_common.o


include make_rules.inc

${objdir}/aed_external.o: ${srcdir}/aed_external.F90 ${objdir}/aed_core.o ${incdir}/aed.h
	$(F90) $(FFLAGS) -DLIBDEF -c $< -o $@

${objdir}/aed_common.o: ${srcdir}/aed_common.F90 ${srcdir}/aed_core.F90 ${objdir}/aed_external.o ${incdir}/aed.h

