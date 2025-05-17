!###############################################################################
!#                                                                             #
!# aed_water.F90                                                               #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013 - 2025 -  The University of Western Australia               #
!#                                                                             #
!#   AED is free software: you can redistribute it and/or modify               #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   AED is distributed in the hope that it will be useful,                    #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created May 2020                                                            #
!#                                                                             #
!###############################################################################

#include "aed.h"

!###############################################################################
MODULE aed_water
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_sedflux
   USE aed_oxygen
   USE aed_silica
   USE aed_carbon
   USE aed_nitrogen
   USE aed_phosphorus
   USE aed_organic_matter
   USE aed_phytoplankton
   USE aed_zooplankton
   USE aed_tracer
   USE aed_noncohesive
   USE aed_totals
   USE aed_dummy
   USE aed_bio_particles
   USE aed_geochemistry
   USE aed_pathogens
   USE aed_pesticides
   USE aed_habitat_water

   IMPLICIT NONE


   !#---------------------------------------------------------------------------

   PRIVATE   !# By default make everything private

   PUBLIC aed_new_wat_model, aed_print_wat_version

CONTAINS
!===============================================================================


!###############################################################################
FUNCTION aed_new_wat_model(modelname) RESULT(model)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: modelname
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
   CHARACTER(len=4) :: prefix
!
!-------------------------------------------------------------------------------
!BEGIN
   NULLIFY(model)

   SELECT CASE (modelname)
      CASE ('aed_bio_particles');  prefix = 'PTM'; ALLOCATE(aed_bio_particles_data_t::model)
      CASE ('aed_sedflux');        prefix = 'SDF'; ALLOCATE(aed_sedflux_data_t::model)
      CASE ('aed_oxygen');         prefix = 'OXY'; ALLOCATE(aed_oxygen_data_t::model)
      CASE ('aed_silica');         prefix = 'SIL'; ALLOCATE(aed_silica_data_t::model)
      CASE ('aed_carbon');         prefix = 'CAR'; ALLOCATE(aed_carbon_data_t::model)
      CASE ('aed_nitrogen');       prefix = 'NIT'; ALLOCATE(aed_nitrogen_data_t::model)
      CASE ('aed_phosphorus');     prefix = 'PHS'; ALLOCATE(aed_phosphorus_data_t::model)
      CASE ('aed_organic_matter'); prefix = 'OGM'; ALLOCATE(aed_organic_matter_data_t::model)
      CASE ('aed_phytoplankton');  prefix = 'PHY'; ALLOCATE(aed_phytoplankton_data_t::model)
      CASE ('aed_zooplankton');    prefix = 'ZOO'; ALLOCATE(aed_zooplankton_data_t::model)
      CASE ('aed_tracer');         prefix = 'TRC'; ALLOCATE(aed_tracer_data_t::model)
      CASE ('aed_noncohesive');    prefix = 'NCS'; ALLOCATE(aed_noncohesive_data_t::model)
      CASE ('aed_geochemistry');   prefix = 'GEO'; ALLOCATE(aed_geochemistry_data_t::model)
      CASE ('aed_pathogens');      prefix = 'PAT'; ALLOCATE(aed_pathogens_data_t::model)
      CASE ('aed_pesticides');     prefix = 'PST'; ALLOCATE(aed_pesticides_data_t::model)
      CASE ('aed_totals');         prefix = 'TOT'; ALLOCATE(aed_totals_data_t::model)
      CASE ('aed_dummy');          prefix = 'DUM'; ALLOCATE(aed_dummy_data_t::model)
      CASE ('aed_habitat_water');  prefix = 'HAB'; ALLOCATE(aed_habitat_water_data_t::model)
!     CASE DEFAULT;                print *,'*** Unknown module ', TRIM(modelname)
   END SELECT

   IF (ASSOCIATED(model)) THEN
      model%aed_model_name = modelname
      model%aed_model_prefix = prefix
   ENDIF
END FUNCTION aed_new_wat_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_print_wat_version
!-------------------------------------------------------------------------------
!BEGIN
   print*,"    libaed-water version ", TRIM(AED_VERSION)
#ifdef __INTEL_COMPILER
   print*,"    libaed built using intel fortran version ", __INTEL_COMPILER
#else
# ifdef __PGI
   print*,"    libaed built using pgfortran version ", __PGIC__, '.', __PGIC_MINOR__, '.', __PGIC_PATCHLEVEL__
# else
#  ifdef __GNUC__
    print*,"    libaed built using gfortran version ", __GNUC__, '.', __GNUC_MINOR__, '.', __GNUC_PATCHLEVEL__
#  else
#   ifdef __clang__
     print*,"    libaed built using flang version ", __clang_major__, '.', __clang_minor__, '.', __clang_patchlevel__
#   else
     print*,"    libaed built using unknown fortran version "
#   endif
#  endif
# endif
#endif
END SUBROUTINE aed_print_wat_version
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!===============================================================================
END MODULE aed_water
