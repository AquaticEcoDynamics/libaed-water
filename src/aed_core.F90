!###############################################################################
!#                                                                             #
!# aed_core.F90                                                                #
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
!# Created June 2013                                                           #
!#                                                                             #
!###############################################################################

#include "aed.h"

MODULE aed_core

   IMPLICIT NONE

   PRIVATE  !# Everything defaults to private

   PUBLIC aed_model_data_t, aed_variable_t, aed_column_t
   PUBLIC aed_init_core, aed_get_var, aed_core_status
   PUBLIC aed_set_current_model, aed_set_prefix
   PUBLIC aed_is_const_var, aed_set_const_var

   PUBLIC aed_define_variable,         aed_define_sheet_variable
   PUBLIC aed_define_diag_variable,    aed_define_sheet_diag_variable
   PUBLIC aed_locate_variable,         aed_locate_sheet_variable

   PUBLIC aed_provide_global,          aed_provide_sheet_global
   PUBLIC aed_locate_global,           aed_locate_sheet_global

   PUBLIC host_has_cell_vel
   PUBLIC zero_, one_, nan_, misval_, secs_per_day
   PUBLIC model_list, last_model
   PUBLIC n_aed_models

   !#---------------------------------------------------------------------------
   TYPE :: aed_prefix_list_t
      CHARACTER(len=4)  :: aed_model_prefix
      CLASS(aed_prefix_list_t),POINTER :: next => null()
   END TYPE aed_prefix_list_t
   !#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !#---------------------------------------------------------------------------
   TYPE :: aed_model_data_t
      INTEGER           :: aed_model_id
      CHARACTER(len=64) :: aed_model_name
      CHARACTER(len=4)  :: aed_model_prefix
      LOGICAL           :: aed_model_zone_avg = .FALSE.
      CLASS(aed_model_data_t),POINTER :: next => null()
      !# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CONTAINS
         procedure :: define             => aed_define
         procedure :: initialize         => aed_initialize
         procedure :: initialize_benthic => aed_initialize_benthic
         procedure :: calculate_surface  => aed_calculate_surface
         procedure :: calculate          => aed_calculate
         procedure :: calculate_benthic  => aed_calculate_benthic
         procedure :: calculate_column   => aed_calculate_column
         procedure :: calculate_riparian => aed_calculate_riparian
         procedure :: calculate_dry      => aed_calculate_dry
         procedure :: equilibrate        => aed_equilibrate
         procedure :: light_extinction   => aed_light_extinction
         procedure :: rain_loss          => aed_rain_loss
         procedure :: light_shading      => aed_light_shading
         procedure :: bio_drag           => aed_bio_drag
         procedure :: particle_bgc       => aed_particle_bgc
         procedure :: mobility           => aed_mobility
         procedure :: validate           => aed_validate
         procedure :: inflow_update      => aed_inflow_update
         procedure :: delete             => aed_delete
      !# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   END TYPE aed_model_data_t
   !#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !#---------------------------------------------------------------------------
   TYPE :: aed_variable_t
      CHARACTER(len=64) :: name
      CLASS(aed_model_data_t), POINTER :: model => null()
      CHARACTER(len=128):: longname
      CHARACTER(len=24) :: units
      AED_REAL          :: initial
      AED_REAL          :: minimum
      AED_REAL          :: maximum
      AED_REAL          :: mobility
      AED_REAL          :: light_extinction
      LOGICAL           :: sheet, diag, extern, found
      LOGICAL           :: top, bot, const
      LOGICAL           :: zavg, zavg_req
      INTEGER           :: particle_link
      CLASS(aed_prefix_list_t),POINTER :: req => null()
   END TYPE aed_variable_t
   !#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !#---------------------------------------------------------------------------
   TYPE :: aed_column_t
      AED_REAL,DIMENSION(:),POINTER :: cell
      AED_REAL,             POINTER :: cell_sheet
      AED_REAL,             POINTER :: flux_atm
      AED_REAL,DIMENSION(:),POINTER :: flux_pel
      AED_REAL,             POINTER :: flux_ben
      AED_REAL,             POINTER :: flux_rip
   END TYPE aed_column_t
   !#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!-------------------------------------------------------------------------------
!MODULE VARIABLES
   INTEGER :: cur_mod_base = 0
   INTEGER :: n_aed_models = 0
   INTEGER :: n_aed_vars = 0, a_vars = 0
   INTEGER :: n_vars = 0, n_sheet_vars = 0
   INTEGER :: n_diags = 0, n_sheet_diags = 0

   TYPE(aed_variable_t),DIMENSION(:),ALLOCATABLE,TARGET :: all_vars

   CHARACTER(len=4),POINTER :: aed_cur_prefix => null()
   CHARACTER(len=80) :: base_aed_directory = '.'
   CLASS(aed_model_data_t), POINTER :: current_model => null()

   LOGICAL :: host_has_cell_vel = .false.

#ifdef f2003
   USE, intrinsic :: iso_fortran_env, ONLY : stdin=>input_unit, &
                                             stdout=>output_unit, &
                                             stderr=>error_unit
#else
#  define stdin  5
#  define stdout 6
#  define stderr 0
#endif
   INTEGER :: log = stderr

   !#---------------------------------------------------------------------------

   CLASS(aed_model_data_t), POINTER :: model_list => null()
   CLASS(aed_model_data_t), POINTER :: last_model => null()

!-------------------------------------------------------------------------------
!CONSTANTS
   AED_REAL,PARAMETER :: zero_ = 0., one_ = 1.
   AED_REAL,PARAMETER :: secs_per_day = 86400.
   AED_REAL,PARAMETER :: misval_ = -9999.

   AED_REAL,PARAMETER :: nan_ = zero_ / zero_

!===============================================================================
CONTAINS


!###############################################################################
INTEGER FUNCTION aed_init_core(dname,have_cell_vel)
!-------------------------------------------------------------------------------
! Initialise the aed model library core routines
!-------------------------------------------------------------------------------
!ARGUMENTS - none
   CHARACTER(len=*) :: dname
   LOGICAL,OPTIONAL :: have_cell_vel
!
!LOCAL
!   AED_REAL :: tmpr = zero_
!
!-------------------------------------------------------------------------------
!BEGIN
   base_aed_directory = dname
   n_aed_vars = 0 ; a_vars = 0
   n_vars = 0;  n_sheet_vars = 0
   n_diags = 0; n_sheet_diags = 0
   IF (PRESENT(have_cell_vel)) host_has_cell_vel = have_cell_vel
   aed_init_core = 0
END FUNCTION aed_init_core
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE display_var(var, idx)
!-------------------------------------------------------------------------------
! Status of the aed model library core routines
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(aed_variable_t),INTENT(in) :: var
   INTEGER :: idx
!
!LOCALS
   CHARACTER(80) :: line
   CHARACTER(80) :: l2
   CLASS(aed_prefix_list_t),POINTER :: req => null()
!
!-------------------------------------------------------------------------------
!BEGIN
! I_o                   FV                    2D                        -                           OGM, PHY, MAG
! Wind                  FV                    2D                        -                           OGM, PHY, MAG
! Temp                  FV                    3D                        -                           OXY,CAR,NIT,OGM, PHY, MAG

!  WRITE(line, *) idx
!  line = TRIM(ADJUSTL(line))
!  line(1:4) = ADJUSTR(line(1:4))

!  line = line(1:4) // "  " // TRIM(var%name) // '                    '
   line = TRIM(var%name) // '                    '
   l2 = line(1:20)
   IF ( var%found .AND. ASSOCIATED(var%model) ) THEN
      line = l2(1:20) // ' ' // var%model%aed_model_name
   ELSE
      line = l2(1:20) // ' ???'
!     print log,'Requested variable ', TRIM(var%name), ' not defined.'
   ENDIF
   line = TRIM(line) // '             '

   l2 = line(1:40)
   IF ( var%sheet ) THEN
      !# IF ( ASSOCIATED(var%model) .and. var%model%aed_model_zone_avg ) THEN
      !# The above will segfault if var%model is null - like pascal, it seems fortran
      !# also evaluates both conditions regardless
      IF ( ASSOCIATED(var%model) ) THEN
         IF ( var%model%aed_model_zone_avg ) THEN
            line = l2(1:40) // ' SZ'
         ELSE
            line = l2(1:40) // ' 2D'
         ENDIF
      ELSE
         line = l2(1:40) // ' 2D'
      ENDIF
   ELSE
      line = l2(1:40) // ' 3D'
   ENDIF

   req => var%req

   IF ( var%extern ) THEN
      line = TRIM(line) // '   ---'
   ELSE IF ( ASSOCIATED(req) ) THEN
      line = TRIM(line) // '   ' // req%aed_model_prefix
      req => req%next
   ENDIF
   IF ( ASSOCIATED(req) ) THEN
      line = TRIM(line) // "   " // req%aed_model_prefix

      DO WHILE ( ASSOCIATED(req%next) )
         req => req%next
         IF (LEN_TRIM(line) >= 75 .AND. ASSOCIATED(req%next)) THEN
             write(log, *) TRIM(line),","
             line(1:51)=' '
             l2 = line(1:51)
             line = l2(1:51) // req%aed_model_prefix
         ELSE
             line = TRIM(line) // ", " // req%aed_model_prefix
         ENDIF
      ENDDO
   ENDIF
   IF ( var%zavg_req ) THEN
      IF ( ASSOCIATED(var%model) .AND. .NOT. var%model%aed_model_zone_avg ) THEN
         line = TRIM(line) // "     (zavg req)"
      ENDIF
   ENDIF

   write(log, *) TRIM(line)
END SUBROUTINE display_var
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION aed_core_status(n_v, n_sv, n_d, n_sd, logit)
!-------------------------------------------------------------------------------
! Status of the aed model library core routines
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(out) :: n_v, n_sv, n_d, n_sd
   LOGICAL,INTENT(in),OPTIONAL :: logit
!
!LOCALS
   INTEGER :: i
!
!-------------------------------------------------------------------------------
!BEGIN
   write(log, *)
   write(log, *) ' ---------------------- AED Variables Summary ----------------------'
   write(log, *) 'Var name           | Module           | Type | ID | Usage (ie who linked to me)'
   write(log, *)
   write(log, *) 'ENVIRONMENT:'
   DO i=1,n_aed_vars
      IF ( all_vars(i)%extern ) CALL display_var(all_vars(i), i)
   ENDDO

   write(log, *)
   write(log, *) 'STATE:'
   n_vars = 0 ; n_sheet_vars = 0
   DO i=1,n_aed_vars
      IF ( .NOT. all_vars(i)%extern .AND. .NOT. all_vars(i)%diag ) THEN
         CALL display_var(all_vars(i), i)
         IF ( all_vars(i)%sheet ) THEN ; n_sheet_vars = n_sheet_vars + 1
         ELSE ; n_vars = n_vars + 1 ; ENDIF
      ENDIF
   ENDDO

   write(log, *)
   write(log, *) 'DIAGNOSTIC:'
   n_diags = 0 ; n_sheet_diags = 0;
   DO i=1,n_aed_vars
      IF ( .NOT. all_vars(i)%extern .AND. all_vars(i)%diag ) THEN
         CALL display_var(all_vars(i), i)
         IF ( all_vars(i)%sheet ) THEN ; n_sheet_diags = n_sheet_diags + 1
         ELSE ; n_diags = n_diags + 1 ; ENDIF
      ENDIF
   ENDDO
   write(log, *)
   write(log, *) ' -------------------------------------------------------------------'
   write(log, *)

   n_v = n_vars;  n_sv = n_sheet_vars
   n_d = n_diags; n_sd = n_sheet_diags
   aed_core_status = n_aed_vars
END FUNCTION aed_core_status
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_set_current_model(model)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS(aed_model_data_t), POINTER :: model
!
!-------------------------------------------------------------------------------
!BEGIN
   current_model => model
END SUBROUTINE aed_set_current_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_set_prefix(prefix)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=4),POINTER, INTENT(in) :: prefix
!
!-------------------------------------------------------------------------------
!BEGIN
   aed_cur_prefix => prefix
END SUBROUTINE aed_set_prefix
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!# These create space for the variables :

!###############################################################################
SUBROUTINE extend_allocated_variables(pcount)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: pcount
!
!LOCALS
   TYPE(aed_variable_t),DIMENSION(:),ALLOCATABLE :: tmp
   INTEGER count, i
!
!-------------------------------------------------------------------------------
!BEGIN
   count = pcount
   IF ( count < 10 ) count = 10
   IF ( ALLOCATED(all_vars) ) THEN
      ALLOCATE(tmp(1:a_vars))
      tmp = all_vars
      DEALLOCATE(all_vars)
      ALLOCATE(all_vars(1:a_vars+count))
      all_vars(1:a_vars) = tmp(1:a_vars)
      DEALLOCATE(tmp);
   ELSE
      ALLOCATE(all_vars(1:count))
   ENDIF

   DO i=1,count
      NULLIFY(all_vars(a_vars+i)%model)
   ENDDO
   all_vars(a_vars+1:a_vars+count)%initial = nan_
   all_vars(a_vars+1:a_vars+count)%minimum = nan_
   all_vars(a_vars+1:a_vars+count)%maximum = nan_
   all_vars(a_vars+1:a_vars+count)%mobility = nan_
   all_vars(a_vars+1:a_vars+count)%light_extinction = nan_

   all_vars(a_vars+1:a_vars+count)%sheet = .false.
   all_vars(a_vars+1:a_vars+count)%diag = .false.
   all_vars(a_vars+1:a_vars+count)%extern = .false.
   all_vars(a_vars+1:a_vars+count)%found = .false.

   all_vars(a_vars+1:a_vars+count)%top = .false.
   all_vars(a_vars+1:a_vars+count)%bot = .false.
   all_vars(a_vars+1:a_vars+count)%zavg = .true.
   all_vars(a_vars+1:a_vars+count)%zavg_req = .false.
   all_vars(a_vars+1:a_vars+count)%const = .false.

   a_vars = a_vars + count

END SUBROUTINE extend_allocated_variables
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed_find_variable(name) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   ret = 0
   IF ( .NOT. ALLOCATED(all_vars) ) RETURN
   DO ret=1,n_aed_vars
      IF ( all_vars(ret)%name == name ) RETURN
      IF ( .NOT. ASSOCIATED(all_vars(ret)%model) ) CYCLE
      IF ( TRIM(all_vars(ret)%model%aed_model_name)//'_'//  &
                                  TRIM(all_vars(ret)%name) == name ) RETURN
!print*,TRIM(name), ' is neither ', TRIM(all_vars(ret)%name), ' nor ',  &
!       TRIM(all_vars(ret)%model%aed_model_name)//'_'//TRIM(all_vars(ret)%name)
   ENDDO
   ret = 0
END FUNCTION aed_find_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE extend_requested(var, prefix)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(aed_variable_t),INTENT(inout) :: var
   CHARACTER(len=4), INTENT(in) :: prefix
!
!LOCALS
   CLASS(aed_prefix_list_t),POINTER :: req => null()
!-------------------------------------------------------------------------------
!BEGIN
   req => var%req

   IF ( .NOT. ASSOCIATED(req) ) THEN
      ALLOCATE(var%req)
      var%req%aed_model_prefix = prefix
      var%req%next => null()
      return
   ENDIF

   !# if we already have it, get out ...
   IF ( req%aed_model_prefix == prefix ) RETURN

   DO WHILE ( ASSOCIATED(req%next) )
      req => req%next
      IF ( req%aed_model_prefix == prefix ) RETURN
   ENDDO
   ALLOCATE(req%next)
   req%next%aed_model_prefix = prefix
   req%next%next => null()
END SUBROUTINE extend_requested
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed_create_variable(name, longname, units, place) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name, longname, units
   LOGICAL, INTENT(in) :: place
!
!LOCALS
   INTEGER :: ret
   CHARACTER(len=64) :: tname
!-------------------------------------------------------------------------------
!BEGIN
   IF ( ASSOCIATED(aed_cur_prefix) .AND. .NOT. place ) THEN
      tname = TRIM(aed_cur_prefix)//"_"//name
   ELSE
      tname = name
   ENDIF

   !# First see if it already exists
   ret = aed_find_variable(tname)
   IF ( ret == 0 ) THEN
      IF (n_aed_vars .GE. a_vars) CALL extend_allocated_variables(n_aed_vars+1-a_vars)
      n_aed_vars = n_aed_vars + 1
      ret = n_aed_vars

      all_vars(ret)%name = tname
!print*,"CREATE Variable '",TRIM(tname),"' long name '",TRIM(longname),"' with units '",TRIM(units),"' at ", ret

!     IF ( ASSOCIATED(current_model) ) THEN
         all_vars(ret)%model => current_model
!     ENDIF
      all_vars(ret)%longname = longname
      all_vars(ret)%units = units

      all_vars(ret)%sheet = .false.
      all_vars(ret)%diag = .false.
      all_vars(ret)%extern = .false.
      all_vars(ret)%found = .false.
      all_vars(ret)%zavg = .false.
      all_vars(ret)%zavg_req = .false.
      all_vars(ret)%top = .false.
      all_vars(ret)%bot = .false.
      all_vars(ret)%particle_link = 0
   ENDIF
   IF ( ASSOCIATED(aed_cur_prefix) ) THEN
      CALL extend_requested(all_vars(ret), aed_cur_prefix)
   ENDIF
END FUNCTION aed_create_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!# These define names and defaults for variables, returning their index:

!###############################################################################
FUNCTION aed_define_variable(name, units, longname, initial, minimum, maximum, mobility) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name, longname, units
   AED_REAL,INTENT(in),OPTIONAL :: initial, minimum, maximum
   AED_REAL,INTENT(in),OPTIONAL :: mobility
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   ret = aed_create_variable(name, longname, units, .false.)

   if ( present(initial) ) all_vars(ret)%initial = initial
   if ( present(minimum) ) all_vars(ret)%minimum = minimum
   if ( present(maximum) ) all_vars(ret)%maximum = maximum
   if ( present(mobility) ) all_vars(ret)%mobility = mobility

!  all_vars(ret)%sheet = .false.
!  all_vars(ret)%diag = .false.
!  all_vars(ret)%extern = .false.
   all_vars(ret)%found = .true.
END FUNCTION aed_define_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed_define_sheet_variable(name, units, longname, initial, minimum, maximum, surf) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name, longname, units
   AED_REAL,INTENT(in),OPTIONAL :: initial, minimum, maximum
   LOGICAL,INTENT(in),OPTIONAL :: surf
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   ret = aed_create_variable(name, longname, units, .false.)

   if ( present(initial) ) all_vars(ret)%initial = initial
   if ( present(minimum) ) all_vars(ret)%minimum = minimum
   if ( present(maximum) ) all_vars(ret)%maximum = maximum

   all_vars(ret)%sheet = .true.
!  all_vars(ret)%diag = .false.
!  all_vars(ret)%extern = .false.
   all_vars(ret)%found = .true.

   IF ( PRESENT(surf) ) THEN
      all_vars(ret)%bot = .NOT. surf
      all_vars(ret)%top = surf
   ELSE
      all_vars(ret)%bot = .TRUE.
   ENDIF
END FUNCTION aed_define_sheet_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed_define_diag_variable(name, units, longname) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name, longname, units
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   ret = aed_create_variable(name, longname, units, .false.)
   n_diags = n_diags + 1

!  all_vars(ret)%sheet = .false.
   all_vars(ret)%diag = .true.
!  all_vars(ret)%extern = .false.
   all_vars(ret)%found = .true.
END FUNCTION aed_define_diag_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed_define_sheet_diag_variable(name, units, longname, surf, zavg) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name, longname, units
   LOGICAL,INTENT(in),OPTIONAL :: surf
   LOGICAL,INTENT(in),OPTIONAL :: zavg
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   ret = aed_create_variable(name, longname, units, .false.)

   n_sheet_diags = n_sheet_diags + 1

   all_vars(ret)%sheet = .true.
   all_vars(ret)%diag = .true.
!  all_vars(ret)%extern = .false.
   all_vars(ret)%found = .true.

   IF ( PRESENT(surf) ) THEN
      all_vars(ret)%bot = .NOT. surf
      all_vars(ret)%top = surf
   ELSE
      all_vars(ret)%bot = .TRUE.
   ENDIF

   IF ( PRESENT(zavg) ) THEN
      all_vars(ret)%zavg = zavg
      all_vars(ret)%zavg_req = .FALSE.
   ELSE
      all_vars(ret)%zavg = .FALSE.
      all_vars(ret)%zavg_req = .FALSE.
   ENDIF

   ret = n_aed_vars
END FUNCTION aed_define_sheet_diag_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed_provide_global(name, longname, units) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name, longname, units
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   ret = aed_create_variable(name, longname, units, .false.)

!  all_vars(ret)%sheet = .false.
!  all_vars(ret)%diag = .false.
   all_vars(ret)%extern = .true.
   all_vars(ret)%found = .true.
END FUNCTION aed_provide_global
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed_provide_sheet_global(name, longname, units, surf) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name, longname, units
   LOGICAL,INTENT(in),OPTIONAL :: surf
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   ret = aed_create_variable(name, longname, units, .false.)

   all_vars(ret)%sheet = .true.
!  all_vars(ret)%diag = .false.
   all_vars(ret)%extern = .true.
   all_vars(ret)%found = .true.

   IF ( PRESENT(surf) ) THEN
      all_vars(ret)%bot = .NOT. surf
      all_vars(ret)%top = surf
   ELSE
      all_vars(ret)%bot = .TRUE.
   ENDIF
END FUNCTION aed_provide_sheet_global
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!# These return the index for the named variable :

!###############################################################################
FUNCTION aed_locate_variable(name) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   IF ( TRIM(name) == '' ) THEN ; ret = 0; RETURN ; ENDIF

   ret = aed_create_variable(name, '', '', .true.)

   IF ( ret /= 0 ) THEN
     IF ( all_vars(ret)%extern .OR. all_vars(ret)%sheet ) ret = 0
   ENDIF
!  IF ( ret /= 0 ) THEN
!    print*,"VARIABLE ",TRIM(name)," LOCATED at ",ret
!  ELSE
!    print*,"VARIABLE ",TRIM(name)," NOT FOUND"
!  ENDIF
END FUNCTION aed_locate_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed_locate_sheet_variable(name,update_from_zone) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name
   LOGICAL, OPTIONAL       :: update_from_zone  ! requesting module using zones
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   IF ( TRIM(name) == '' ) THEN ; ret = 0 ; RETURN ; ENDIF

   ret = aed_create_variable(name, '', '', .true.)

   IF ( ret /= 0 ) THEN
     IF ( all_vars(ret)%extern .OR. .NOT. all_vars(ret)%sheet ) ret = 0
   ENDIF
   IF ( ret /= 0 .AND. PRESENT(update_from_zone)) THEN
     all_vars(ret)%zavg_req = .true.
   ENDIF
!  IF ( ret /= 0 ) THEN
!    print*,"SHEET_VAR ",TRIM(name)," LOCATED at ",ret
!  ELSE
!    print*,"SHEET_VAR ",TRIM(name)," NOT FOUND"
!  ENDIF
END FUNCTION aed_locate_sheet_variable
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed_locate_global(name) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   IF ( TRIM(name) == '' ) THEN ; ret = 0; RETURN ; ENDIF

   IF ( TRIM(name) == "cell_vel" ) THEN
       ret = -1
       IF ( .NOT. host_has_cell_vel ) RETURN
   ENDIF

   ret = aed_create_variable(name, '', '', .true.)
   all_vars(ret)%extern = .true.
END FUNCTION aed_locate_global
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed_locate_sheet_global(name) RESULT(ret)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   IF ( TRIM(name) == '' ) THEN ; ret = 0; RETURN ; ENDIF

   ret = aed_create_variable(name, '', '', .true.)

   all_vars(ret)%sheet = .true.
   all_vars(ret)%extern = .true.
END FUNCTION aed_locate_sheet_global
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION aed_get_var(which, tvar)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: which
   TYPE(aed_variable_t),POINTER,INTENT(out) :: tvar
!
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
   aed_get_var = .TRUE. ; NULLIFY(tvar)
   IF (which > 0 .AND. which <= n_aed_vars) THEN ; tvar => all_vars(which)
   ELSE ; aed_get_var = .FALSE. ; ENDIF
END FUNCTION aed_get_var
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION aed_is_const_var(which)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: which
!
!-------------------------------------------------------------------------------
!BEGIN
   aed_is_const_var = .FALSE.
   IF (which > 0 .AND. which <= n_aed_vars) &
      aed_is_const_var =  all_vars(which)%const
END FUNCTION aed_is_const_var
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_set_const_var(which)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: which
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (which > 0 .AND. which <= n_aed_vars) &
      all_vars(which)%const = .TRUE.
END SUBROUTINE aed_set_const_var
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!#                                                                             #
!# These are dummies for the declarations.                                     #
!#                                                                             #
!###############################################################################


!###############################################################################
SUBROUTINE aed_define(data, namlst)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!-------------------------------------------------------------------------------
!print*,"Default aed_define ", TRIM(data%aed_model_name)
END SUBROUTINE aed_define
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_initialize(data,column, layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------
!print*,"Default aed_initialize ", TRIM(data%aed_model_name)
END SUBROUTINE aed_initialize
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_initialize_benthic(data,column, layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------
!print*,"Default aed_initialize_benthic ", TRIM(data%aed_model_name)
END SUBROUTINE aed_initialize_benthic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate(data,column,layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------
!print*,"Default aed_calculate ", TRIM(data%aed_model_name)
END SUBROUTINE aed_calculate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_surface(data,column,layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------
!print*,"Default aed_calculate_surface ", TRIM(data%aed_model_name)
END SUBROUTINE aed_calculate_surface
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic(data,column,layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------
!print*,"Default aed_calculate_benthic ", TRIM(data%aed_model_name)
END SUBROUTINE aed_calculate_benthic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_column(data,column,layer_map)
   !-------------------------------------------------------------------------------
      CLASS (aed_model_data_t),INTENT(in) :: data
      TYPE (aed_column_t),INTENT(inout) :: column(:)
      INTEGER,INTENT(in) :: layer_map(:)
!-------------------------------------------------------------------------------
!print*,"Default aed_calculate_benthic ", TRIM(data%aed_model_name)
END SUBROUTINE aed_calculate_column
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_riparian(data,column,layer_idx, pc_wet)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: pc_wet
!-------------------------------------------------------------------------------
!print*,"Default aed_calculate_riparian ", TRIM(data%aed_model_name)
END SUBROUTINE aed_calculate_riparian
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_dry(data,column,layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------
!print*,"Default aed_calculate_dry ", TRIM(data%aed_model_name)
END SUBROUTINE aed_calculate_dry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_equilibrate(data,column,layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------
!print*,"Default aed_equilibrate ", TRIM(data%aed_model_name)
END SUBROUTINE aed_equilibrate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION aed_validate(data,column,layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------
!print*,"Default aed_validate ", TRIM(data%aed_model_name)
   aed_validate = .TRUE.
END FUNCTION aed_validate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_mobility(data,column,layer_idx,mobility)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!-------------------------------------------------------------------------------
!print*,"Default aed_mobility ", TRIM(data%aed_model_name)
END SUBROUTINE aed_mobility
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_rain_loss(data, column, layer_idx, infil)
!-------------------------------------------------------------------------------
! Get the soil moisture deficit, so host model can infiltrate rain accordingly
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: infil
!-------------------------------------------------------------------------------
!BEGIN
!print*,"Default aed_rain_loss ", TRIM(data%aed_model_name)
END SUBROUTINE aed_rain_loss
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_shading(data, column, layer_idx, shade_frac)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: shade_frac
!-------------------------------------------------------------------------------
!BEGIN
!print*,"Default aed_light_shading ", TRIM(data%aed_model_name)
END SUBROUTINE aed_light_shading
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_bio_drag(data, column, layer_idx, drag)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: drag
!-------------------------------------------------------------------------------
!BEGIN
!print*,"Default aed_bio_drag ", TRIM(data%aed_model_name)
END SUBROUTINE aed_bio_drag
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_extinction(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!-------------------------------------------------------------------------------
!print*,"Default aed_light_extinction ", TRIM(data%aed_model_name)
END SUBROUTINE aed_light_extinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_particle_bgc(data,column,layer_idx,ppid,partcl)
   CLASS (aed_model_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   INTEGER,INTENT(inout) :: ppid
   AED_REAL,DIMENSION(:),INTENT(inout) :: partcl
!-------------------------------------------------------------------------------
!print*,"Default aed_particle_bgc ", TRIM(data%aed_model_name)
END SUBROUTINE aed_particle_bgc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_inflow_update(data, wqinf, temp, salt)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(in) :: data
   AED_REAL,DIMENSION(:),INTENT(inout) :: wqinf
   AED_REAL,             INTENT(inout) :: temp, salt
!-------------------------------------------------------------------------------
!print*,"Default aed_inflow_update ", TRIM(data%aed_model_name)
END SUBROUTINE aed_inflow_update
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_delete(data)
!-------------------------------------------------------------------------------
   CLASS (aed_model_data_t),INTENT(inout) :: data
!-------------------------------------------------------------------------------
!print*,"Default aed_delete ", TRIM(data%aed_model_name)
END SUBROUTINE aed_delete
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_core
