!###############################################################################
!#                                                                             #
!# aed_common.F90                                                              #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013 - 2023 -  The University of Western Australia               #
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
!# Created Aug 2013                                                            #
!#                                                                             #
!###############################################################################

#include "aed.h"


!###############################################################################
MODULE aed_common
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_water
   USE aed_benthic
   USE aed_demo
   USE aed_riparian
   USE aed_lighting
   USE aed_dev

   IMPLICIT NONE

   !#---------------------------------------------------------------------------

   PRIVATE   !# By default make everything private

   !#---------------------------------------------------------------------------

   PUBLIC aed_define_model, aed_delete, aed_print_version, aed_requested_zones

   !#---------------------------------------------------------------------------

   PUBLIC aed_initialize, aed_initialize_benthic
   PUBLIC aed_calculate, aed_calculate_surface, aed_calculate_benthic
   PUBLIC aed_calculate_riparian, aed_calculate_dry, aed_calculate_column
   PUBLIC aed_light_extinction, aed_light_shading
   PUBLIC aed_equilibrate, aed_mobility, aed_rain_loss
   PUBLIC aed_bio_drag, aed_particle_bgc, aed_inflow_update

   !#---------------------------------------------------------------------------

   !# Re-export these from aed_core.
   PUBLIC aed_model_data_t, aed_variable_t, aed_column_t
   PUBLIC aed_init_core, aed_core_status, aed_get_var
   PUBLIC aed_provide_global, aed_provide_sheet_global

   !#---------------------------------------------------------------------------

   PUBLIC zero_, one_, nan_, secs_per_day, misval_

   !#---------------------------------------------------------------------------

   INTEGER,PARAMETER :: NO_ZONES = 1
   INTEGER,PARAMETER :: DO_ZONES = 2

CONTAINS
!===============================================================================


!###############################################################################
SUBROUTINE aed_build_model(model, namlst, do_prefix)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_model_data_t),POINTER :: model
   INTEGER,INTENT(in)      :: namlst
   LOGICAL,INTENT(in)      :: do_prefix
!
!LOCALS
   CHARACTER(len=4),POINTER :: prefix_p
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL aed_set_current_model(model)
   IF ( do_prefix ) THEN
      prefix_p => model%aed_model_prefix
      CALL aed_set_prefix(prefix_p)
   ENDIF
   CALL model%define(namlst)
   IF ( do_prefix ) THEN
      prefix_p => null()
      CALL aed_set_prefix(prefix_p)
   ENDIF
   CALL aed_set_current_model(null())
END SUBROUTINE aed_build_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION scan_name(modeldef, flags) RESULT(modelname)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in)  :: modeldef
   INTEGER(4), INTENT(out)  :: flags
!
!LOCALS
   INTEGER :: len, i
   CHARACTER(len=64) :: modelname
!
!-------------------------------------------------------------------------------
!BEGIN
   flags = 0
   modelname = ''
   len = LEN_TRIM(modeldef)

   DO i=1,len
      IF (modeldef(i:i) == ':') EXIT
      modelname(i:i) = modeldef(i:i)
   ENDDO

   IF ( i >= len ) RETURN

   DO WHILE ( i <= len )
      IF ( modeldef(i:i+1) == 'nz' ) flags = IOR(flags, NO_ZONES)
      IF ( modeldef(i:i+1) == 'za' ) flags = IOR(flags, DO_ZONES)

      DO WHILE ( i <= len .AND. modeldef(i:i) /= ':' ) ; i = i + 1 ; ENDDO
      IF ( i <= len .AND. modeldef(i:i) == ':' ) i = i + 1
   ENDDO

   modelname = TRIM(modelname)
END FUNCTION scan_name
!===============================================================================


!###############################################################################
SUBROUTINE aed_define_model(modeldef, namlst)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: modeldef
   INTEGER,INTENT(in)      :: namlst
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
   CHARACTER(len=64) :: modelname
   INTEGER :: flags = 0
!
!-------------------------------------------------------------------------------
!BEGIN
   modelname = scan_name(modeldef, flags)
   NULLIFY(model)
   model => aed_new_wat_model(modelname)
   IF (.NOT. ASSOCIATED(model)) model => aed_new_ben_model(modelname)
   IF (.NOT. ASSOCIATED(model)) model => aed_new_dmo_model(modelname)
   IF (.NOT. ASSOCIATED(model)) model => aed_new_rip_model(modelname)
   IF (.NOT. ASSOCIATED(model)) model => aed_new_lgt_model(modelname)
   IF (.NOT. ASSOCIATED(model)) model => aed_new_dev_model(modelname)

   IF ( ASSOCIATED(model) ) THEN
      n_aed_models = n_aed_models + 1
      model%aed_model_id = n_aed_models

      IF ( IAND(flags, DO_ZONES) /= 0 ) &
         model%aed_model_zone_avg = .TRUE.
      IF ( IAND(flags, NO_ZONES) /= 0 ) &
         model%aed_model_zone_avg = .FALSE.

      CALL aed_build_model(model, namlst, .TRUE.)

      ! report the outcome of special tokens/flags that were set
      IF(model%aed_model_zone_avg) THEN
        print *,'          ******************************* '
        print *,'          module ', TRIM(modelname),' : configured for zone_averaging '
        print *,'          ******************************* '
      ENDIF

      IF ( .NOT. ASSOCIATED(model_list) ) model_list => model
      IF ( ASSOCIATED(last_model) ) last_model%next => model
      last_model => model
   ELSE
      print *,'*** Unknown module ', TRIM(modelname)
      STOP
   ENDIF
END SUBROUTINE aed_define_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_print_version
!-------------------------------------------------------------------------------
   CALL aed_print_wat_version
   CALL aed_print_ben_version
   CALL aed_print_dmo_version
   CALL aed_print_rip_version
   CALL aed_print_lgt_version
   CALL aed_print_dev_version
END SUBROUTINE aed_print_version
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION aed_requested_zones(n_aed_vars)
!-------------------------------------------------------------------------------
   INTEGER, INTENT(in) :: n_aed_vars
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
   TYPE  (aed_variable_t)  ,POINTER :: tvar
   INTEGER :: i, j
   LOGICAL :: res = .FALSE., err = .FALSE.
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))

      IF ( model%aed_model_zone_avg ) THEN
         ! This model is set to be a averaged model; check & ensure zavg = true
         res = .TRUE.
         j = 0
         DO i=1, n_aed_vars
            IF ( aed_get_var(i, tvar) ) THEN

               IF (tvar%extern) CYCLE  ! external (environment) not currently able to be updated from zones

               IF ( tvar%model%aed_model_id .EQ. model%aed_model_id ) THEN
                  ! sheet variables in a zone averaged model must have zavg=T
                  IF ( tvar%sheet ) THEN
                    ! non-environent variable sheets can be averaged
                    j = j + 1
                    tvar%zavg = .TRUE.
                    print *,'        zone averaged variable: ',TRIM(tvar%name)
                  ELSE
                    ! averageing requests for env or non-sheets are not possible
                    IF ( tvar%zavg ) err = .TRUE.
                    tvar%zavg = .FALSE.
                  ENDIF

               ELSE

                  ! this model requested zone averaged updates from others
                  IF ( tvar%zavg_req ) THEN
                     IF ( tvar%sheet ) THEN
                        ! non-environent variable sheets can be averaged
                        j = j + 1
                        tvar%zavg = .TRUE.
                        print *,'        zone averaged variable: ', &
                                 TRIM(tvar%name)//'   (linked by',model%aed_model_id,TRIM(model%aed_model_name),')'
                                         !MH glitch here is if order is out then wrong linked model appears
                     ELSE
                        ! averageing requests for env or non-sheets are not possible
                        err = .TRUE.
                        tvar%zavg = .FALSE.
                     ENDIF
                  ENDIF
               ENDIF
               tvar%zavg_req = .FALSE.
            ENDIF
         ENDDO
      ENDIF

      model => model%next
   ENDDO

   aed_requested_zones = res

END FUNCTION aed_requested_zones
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!#                                                                             #
!# These are wrappers for the individual models.                               #
!#                                                                             #
!###############################################################################


!###############################################################################
SUBROUTINE aed_initialize(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%initialize(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed_initialize
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_initialize_benthic(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%initialize_benthic(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed_initialize_benthic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed_calculate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_surface(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate_surface(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed_calculate_surface
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic(column, layer_idx, do_zones)
!-------------------------------------------------------------------------------
! The benthic routine may be grouped in zones by the global do_zone_averaging
! flag, however a new model level flag allows us to not average in zones.
! This routine takes the optional argument "do_zones" which should only be
! passed if do_zone_averaging is on.
! If do_zones is not present we can call every models calculate_benthic
! routine.
! If it IS present we pass true when called from inside the zone calculations,
! but false when called from the normal flux calculation section.
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   LOGICAL,OPTIONAL,INTENT(in) :: do_zones
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   IF ( PRESENT(do_zones) ) THEN
      DO WHILE (ASSOCIATED(model))
         IF ( model%aed_model_zone_avg .EQV. do_zones ) &
            CALL model%calculate_benthic(column, layer_idx)
         model => model%next
      ENDDO
   ELSE
      DO WHILE (ASSOCIATED(model))
         CALL model%calculate_benthic(column, layer_idx)
         model => model%next
      ENDDO
   ENDIF
END SUBROUTINE aed_calculate_benthic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_riparian(column, layer_idx, pc_wet)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: pc_wet
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate_riparian(column, layer_idx, pc_wet)
      model => model%next
   ENDDO
END SUBROUTINE aed_calculate_riparian
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_column(column, layer_map)
   !-------------------------------------------------------------------------------
      TYPE (aed_column_t),INTENT(inout) :: column(:)
      INTEGER,INTENT(in) :: layer_map(:)
   !
   !LOCALS
      CLASS (aed_model_data_t),POINTER :: model
   !-------------------------------------------------------------------------------
      model => model_list
      DO WHILE (ASSOCIATED(model))
         CALL model%calculate_column(column, layer_map)
         model => model%next
      ENDDO
   END SUBROUTINE aed_calculate_column
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !###############################################################################
SUBROUTINE aed_calculate_dry(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%calculate_dry(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed_calculate_dry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_equilibrate(column, layer_idx)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%equilibrate(column, layer_idx)
      model => model%next
   ENDDO
END SUBROUTINE aed_equilibrate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed_validate(column, layer_idx) RESULT(valid)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
   LOGICAL :: valid
!-------------------------------------------------------------------------------
   valid = .TRUE.
   model => model_list
   DO WHILE (ASSOCIATED(model))
      IF (.NOT. model%validate(column, layer_idx)) valid = .FALSE.
      model => model%next
   ENDDO
END FUNCTION aed_validate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_extinction(column, layer_idx, extinction)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   extinction = zero_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%light_extinction(column, layer_idx, extinction)
      model => model%next
   ENDDO
END SUBROUTINE aed_light_extinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_mobility(column, layer_idx, mobility)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   !mobility = zero_ !MH leave this as is in case default settling vals provided
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%mobility(column, layer_idx, mobility)
      model => model%next
   ENDDO
END SUBROUTINE aed_mobility
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_rain_loss(column, layer_idx, infil)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: infil
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   infil = zero_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%rain_loss(column, layer_idx, infil)
      model => model%next
   ENDDO
END SUBROUTINE aed_rain_loss
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_shading(column, layer_idx, shade_frac)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: shade_frac
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
!BEGIN
   shade_frac = one_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%light_shading(column, layer_idx, shade_frac)
      model => model%next
   ENDDO
END SUBROUTINE aed_light_shading
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_bio_drag(column, layer_idx, drag)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: drag
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   drag = zero_
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%bio_drag(column, layer_idx, drag)
      model => model%next
   ENDDO
END SUBROUTINE aed_bio_drag
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_particle_bgc(column, layer_idx, ppid, partcl)
!-------------------------------------------------------------------------------
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   INTEGER,INTENT(inout) :: ppid
   AED_REAL,DIMENSION(:),INTENT(inout) :: partcl
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%particle_bgc(column, layer_idx, ppid, partcl)
      model => model%next
   ENDDO
END SUBROUTINE aed_particle_bgc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_inflow_update(wqinf, temp, salt)
!-------------------------------------------------------------------------------
   AED_REAL,DIMENSION(:),INTENT(inout) :: wqinf
   AED_REAL,             INTENT(inout) :: temp, salt
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%inflow_update(wqinf, temp, salt)
      model => model%next
   ENDDO
END SUBROUTINE aed_inflow_update
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_delete
!-------------------------------------------------------------------------------
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
!-------------------------------------------------------------------------------
   model => model_list
   DO WHILE (ASSOCIATED(model))
      CALL model%delete
      model => model%next
   ENDDO
END SUBROUTINE aed_delete
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!===============================================================================
END MODULE aed_common
