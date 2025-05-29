!###############################################################################
!#                                                                             #
!# aed_phyto_abm.F90                                                       #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2018 - 2025 -  The University of Western Australia               #
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
!# Created Jul 2018                                                            #
!#                                                                             #
!###############################################################################

#include "aed.h"

!
MODULE aed_phyto_abm
!-------------------------------------------------------------------------------
! aed_phyto_abm --- test particle model
!
! The AED module test contains basic equations that have no dependencies
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util
   USE aed_bio_utils
   USE prob_functions_mod

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_phyto_abm_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_phyto_abm_data_t
      !# Variable identifiers
      INTEGER :: id_ptm_01, id_ptm_02, id_ptm_03, id_ptm_04,   &
                 id_ptm_05, id_ptm_06, id_ptm_07, id_ptm_08,   &
                 id_ptm_09, id_ptm_10, id_ptm_11, id_ptm_12,   &
                 id_ptm_13, id_ptm_14, id_ptm_15, id_ptm_16,   &
                 id_ptm_17, id_ptm_18
      INTEGER :: id_ptm101, id_ptm102, id_ptm103, id_ptm104,   &
                 id_ptm105, id_ptm106, id_ptm107, id_ptm108,   &
                 id_ptm109, id_ptm110, id_ptm111, id_ptm112,   &
                 id_ptm113, id_ptm114, id_ptm115, id_ptm116,   &
                 id_ptm117, id_ptm118
      INTEGER :: id_ptm_00
      INTEGER :: ip_ptm_c
      INTEGER :: id_d_oxy, id_d_dc, id_d_dn, id_d_dp
      INTEGER :: id_oxy,id_amm,id_nit,id_frp,id_doc,id_don,id_dop
      INTEGER :: id_lht, id_larea, id_dep, id_tem, id_par

      AED_REAL :: vvel_new, vvel_old, decay_rate_new, decay_rate_old
      AED_REAL :: X_dwww, X_cdw, X_nc, X_pc, mass_limit

      ! Phytoplankton parameters
      INTEGER  :: num_phytos
      TYPE(phyto_data_t),DIMENSION(:),ALLOCATABLE :: phytos

      CONTAINS
         PROCEDURE :: define             => aed_define_phyto_abm
!        PROCEDURE :: calculate          => aed_calculate_phyto_abm
!        PROCEDURE :: calculate_benthic  => aed_calculate_benthic_phyto_abm
!        PROCEDURE :: calculate_riparian => aed_calculate_riparian_phyto_abm
!        PROCEDURE :: calculate_dry      => aed_calculate_dry_phyto_abm
!        PROCEDURE :: equilibrate        => aed_equilibrate_phyto_abm
!        PROCEDURE :: mobility           => aed_mobility_phyto_abm
!        PROCEDURE :: light_extinction   => aed_light_extinction_phyto_abm
!        PROCEDURE :: delete             => aed_delete_phyto_abm
         PROCEDURE :: particle_bgc       => aed_particle_bgc_phyto_abm
   END TYPE

   INTEGER, PARAMETER :: PTM_MASS   = 15
   INTEGER, PARAMETER :: PTM_VVEL   = 14
   INTEGER, PARAMETER :: PTM_BIRTH  = 17
   INTEGER, PARAMETER :: PTM_AGE    = 18
   INTEGER, PARAMETER :: PTM_STATUS = 19

   ! PTM array reference index locations
   INTEGER, PARAMETER :: STAT = 1  !#define STAT   0 _PTM_STAT_
   INTEGER, PARAMETER :: IDX2 = 2  !#define IDX2   1 _PTM_STAT_
   INTEGER, PARAMETER :: IDX3 = 3  !#define IDX3   2 _PTM_STAT_
   INTEGER, PARAMETER :: LAYR = 4  !#define LAYR   3 _PTM_STAT_
   INTEGER, PARAMETER :: FLAG = 5  !#define FLAG   4 _PTM_STAT_

   INTEGER, PARAMETER :: MASS = 1  !#define MASS   0 _PTM_ENV_
   INTEGER, PARAMETER :: DIAM = 2  !#define DIAM   1 _PTM_ENV_
   INTEGER, PARAMETER :: DENS = 3  !#define DENS   2 _PTM_ENV_
   INTEGER, PARAMETER :: VVEL = 4  !#define VVEL   3 _PTM_ENV_
   INTEGER, PARAMETER :: HGHT = 5  !#define HGHT   4 _PTM_ENV_

   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs


!===============================================================================
CONTAINS


!###############################################################################
INTEGER FUNCTION load_csv(dbase, pd, dbsize)
!-------------------------------------------------------------------------------
   USE aed_csv_reader
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*),INTENT(in) :: dbase
   TYPE(phyto_param_t) :: pd(MAX_PHYTO_TYPES)
   INTEGER,INTENT(out) :: dbsize
!
!LOCALS
   INTEGER :: unit, nccols, ccol, dcol
   CHARACTER(len=32),POINTER,DIMENSION(:) :: csvnames
   CHARACTER(len=32) :: name
   TYPE(AED_SYMBOL),DIMENSION(:),ALLOCATABLE :: values
   INTEGER :: idx_col = 0
   LOGICAL :: meh
   INTEGER :: ret = 0
!
!BEGIN
!-------------------------------------------------------------------------------
   dbsize = 0
   unit = aed_csv_read_header(dbase, csvnames, nccols)
   IF (unit <= 0) THEN
      load_csv = -1
      RETURN !# No file found
   ENDIF

   !# default values for some
   pd%c1           = 0.1240/60.   ! From Chung et al (2014)
   pd%c3           = 0.0230/60.   !  "
   pd%f1           = 0.675        ! Ross and Sharples (2007)
   pd%f2           = 0.750        !  "
   pd%d_phy        = 1e-5

   ALLOCATE(values(nccols))

   DO WHILE ( aed_csv_read_row(unit, values) )
      DO ccol=2,nccols
         dcol = ccol - 1
         pd(dcol)%p_name = csvnames(ccol)

         CALL copy_name(values(1), name)
         SELECT CASE (name)
            CASE ('p_initial')     ; pd(dcol)%p_initial     = extract_double(values(ccol))
            CASE ('p0')            ; pd(dcol)%p0            = extract_double(values(ccol))
            CASE ('w_p')           ; pd(dcol)%w_p           = extract_double(values(ccol))
            CASE ('Xcc')           ; pd(dcol)%Xcc           = extract_double(values(ccol))
            CASE ('R_growth')      ; pd(dcol)%R_growth      = extract_double(values(ccol))
            CASE ('fT_Method')     ; pd(dcol)%fT_Method     = extract_integer(values(ccol))
            CASE ('theta_growth')  ; pd(dcol)%theta_growth  = extract_double(values(ccol))
            CASE ('T_std')         ; pd(dcol)%T_std         = extract_double(values(ccol))
            CASE ('T_opt')         ; pd(dcol)%T_opt         = extract_double(values(ccol))
            CASE ('T_max')         ; pd(dcol)%T_max         = extract_double(values(ccol))
            CASE ('lightModel')    ; pd(dcol)%lightModel    = extract_integer(values(ccol))
            CASE ('I_K')           ; pd(dcol)%I_K           = extract_double(values(ccol))
            CASE ('I_S')           ; pd(dcol)%I_S           = extract_double(values(ccol))
            CASE ('KePHY')         ; pd(dcol)%KePHY         = extract_double(values(ccol))
            CASE ('f_pr')          ; pd(dcol)%f_pr          = extract_double(values(ccol))
            CASE ('R_resp')        ; pd(dcol)%R_resp        = extract_double(values(ccol))
            CASE ('theta_resp')    ; pd(dcol)%theta_resp    = extract_double(values(ccol))
            CASE ('k_fres')        ; pd(dcol)%k_fres        = extract_double(values(ccol))
            CASE ('k_fdom')        ; pd(dcol)%k_fdom        = extract_double(values(ccol))
            CASE ('salTol')        ; pd(dcol)%salTol        = extract_integer(values(ccol))
            CASE ('S_bep')         ; pd(dcol)%S_bep         = extract_double(values(ccol))
            CASE ('S_maxsp')       ; pd(dcol)%S_maxsp       = extract_double(values(ccol))
            CASE ('S_opt')         ; pd(dcol)%S_opt         = extract_double(values(ccol))
            CASE ('simDINUptake')  ; pd(dcol)%simDINUptake  = extract_integer(values(ccol))
            CASE ('simDONUptake')  ; pd(dcol)%simDONUptake  = extract_integer(values(ccol))
            CASE ('simNFixation')  ; pd(dcol)%simNFixation  = extract_integer(values(ccol))
            CASE ('simINDynamics') ; pd(dcol)%simINDynamics = extract_integer(values(ccol))
            CASE ('N_o')           ; pd(dcol)%N_o           = extract_double(values(ccol))
            CASE ('K_N')           ; pd(dcol)%K_N           = extract_double(values(ccol))
            CASE ('X_ncon')        ; pd(dcol)%X_ncon        = extract_double(values(ccol))
            CASE ('X_nmin')        ; pd(dcol)%X_nmin        = extract_double(values(ccol))
            CASE ('X_nmax')        ; pd(dcol)%X_nmax        = extract_double(values(ccol))
            CASE ('R_nuptake')     ; pd(dcol)%R_nuptake     = extract_double(values(ccol))
            CASE ('k_nfix')        ; pd(dcol)%k_nfix        = extract_double(values(ccol))
            CASE ('R_nfix')        ; pd(dcol)%R_nfix        = extract_double(values(ccol))
            CASE ('simDIPUptake')  ; pd(dcol)%simDIPUptake  = extract_integer(values(ccol))
            CASE ('simIPDynamics') ; pd(dcol)%simIPDynamics = extract_integer(values(ccol))
            CASE ('P_0')           ; pd(dcol)%P_0           = extract_double(values(ccol))
            CASE ('K_P')           ; pd(dcol)%K_P           = extract_double(values(ccol))
            CASE ('X_pcon')        ; pd(dcol)%X_pcon        = extract_double(values(ccol))
            CASE ('X_pmin')        ; pd(dcol)%X_pmin        = extract_double(values(ccol))
            CASE ('X_pmax')        ; pd(dcol)%X_pmax        = extract_double(values(ccol))
            CASE ('R_puptake')     ; pd(dcol)%R_puptake     = extract_double(values(ccol))
            CASE ('simSiUptake')   ; pd(dcol)%simSiUptake   = extract_integer(values(ccol))
            CASE ('Si_0')          ; pd(dcol)%Si_0          = extract_double(values(ccol))
            CASE ('K_Si')          ; pd(dcol)%K_Si          = extract_double(values(ccol))
            CASE ('X_sicon')       ; pd(dcol)%X_sicon       = extract_double(values(ccol))

            CASE ('c1')            ; pd(dcol)%c1            = extract_double(values(ccol))
            CASE ('c3')            ; pd(dcol)%c3            = extract_double(values(ccol))
            CASE ('f1')            ; pd(dcol)%f1            = extract_double(values(ccol))
            CASE ('f2')            ; pd(dcol)%f2            = extract_double(values(ccol))
            CASE ('d_phy')         ; pd(dcol)%d_phy         = extract_double(values(ccol))

            CASE DEFAULT ; print *, 'Unknown row "', TRIM(name), '"'
         END SELECT
      ENDDO
   ENDDO

   meh = aed_csv_close(unit)
   !# don't care if close fails

   IF (ASSOCIATED(csvnames)) DEALLOCATE(csvnames)
   IF (ALLOCATED(values))    DEALLOCATE(values)

   dbsize = nccols-1
   load_csv = ret
END FUNCTION load_csv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_phytoplankton_load_params(data, dbase, count, list, settling, resuspension)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_phyto_abm_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(inout)       :: count
   INTEGER,INTENT(in)          :: list(*)
   INTEGER,INTENT(in)          :: settling(*)
   AED_REAL,INTENT(in)         :: resuspension(*)
!
!LOCALS
   INTEGER  :: status
   INTEGER  :: i,tfil, dbsize = 0
   AED_REAL :: minNut

   TYPE(phyto_param_t),ALLOCATABLE :: pd(:)
   NAMELIST /phyto_data/ pd     ! %% type phyto_param_t - see aed_bio_utils
!-------------------------------------------------------------------------------
!BEGIN
    ALLOCATE(pd(MAX_PHYTO_TYPES))
    SELECT CASE (param_file_type(dbase))
       CASE (CSV_TYPE)
           status = load_csv(dbase, pd, dbsize)
       CASE (NML_TYPE)
           print*,"nml format parameter file is deprecated. Please update to CSV format"
           pd%p_name = ''
           open(NEWUNIT=tfil,file=dbase, status='OLD',iostat=status)
           IF (status /= 0) STOP 'Cannot open phyto_data namelist file'
           read(tfil,nml=phyto_data,iostat=status)
           close(tfil)
           DO i=1,MAX_PHYTO_TYPES
              IF (pd(i)%p_name == '') EXIT
              dbsize = dbsize + 1
           ENDDO
       CASE DEFAULT
           print *,'Unknown file type "',TRIM(dbase),'"'; status=1
    END SELECT
    IF (status /= 0) STOP 'Error reading namelist phyto_data'

    data%num_phytos = 0
    ALLOCATE(data%phytos(count))

    DO i=1,count
       IF ( list(i) < 1 .OR. list(i) > dbsize ) EXIT  !# bad index, exit the loop

       data%num_phytos = data%num_phytos + 1
       ! Assign parameters from database to chosen simulated groups
       ! Note:  all rates are provided in values per day,
       !        and are converted in here to values per second.
       data%phytos(i)%p_name       = pd(list(i))%p_name
       data%phytos(i)%p0           = pd(list(i))%p0
       data%phytos(i)%w_p          = pd(list(i))%w_p/secs_per_day
      !data%phytos(i)%settling     = settling(i)
      !data%phytos(i)%resuspension = resuspension(i)
       data%phytos(i)%Xcc          = pd(list(i))%Xcc
       data%phytos(i)%R_growth     = pd(list(i))%R_growth/secs_per_day
       data%phytos(i)%fT_Method    = pd(list(i))%fT_Method
       data%phytos(i)%theta_growth = pd(list(i))%theta_growth
       data%phytos(i)%T_std        = pd(list(i))%T_std
       data%phytos(i)%T_opt        = pd(list(i))%T_opt
       data%phytos(i)%T_max        = pd(list(i))%T_max
       data%phytos(i)%lightModel   = pd(list(i))%lightModel
       data%phytos(i)%I_K          = pd(list(i))%I_K
       data%phytos(i)%I_S          = pd(list(i))%I_S
       data%phytos(i)%KePHY        = pd(list(i))%KePHY
       data%phytos(i)%f_pr         = pd(list(i))%f_pr
       data%phytos(i)%R_resp       = pd(list(i))%R_resp/secs_per_day
       data%phytos(i)%theta_resp   = pd(list(i))%theta_resp
       data%phytos(i)%k_fres       = pd(list(i))%k_fres
       data%phytos(i)%k_fdom       = pd(list(i))%k_fdom
       data%phytos(i)%salTol       = pd(list(i))%salTol
       data%phytos(i)%S_bep        = pd(list(i))%S_bep
       data%phytos(i)%S_maxsp      = pd(list(i))%S_maxsp
       data%phytos(i)%S_opt        = pd(list(i))%S_opt
       data%phytos(i)%simDINUptake = pd(list(i))%simDINUptake
       data%phytos(i)%simDONUptake = pd(list(i))%simDONUptake
       data%phytos(i)%simNFixation = pd(list(i))%simNFixation
       data%phytos(i)%simINDynamics= pd(list(i))%simINDynamics
       data%phytos(i)%N_o          = pd(list(i))%N_o
       data%phytos(i)%K_N          = pd(list(i))%K_N
       data%phytos(i)%X_ncon       = pd(list(i))%X_ncon
       data%phytos(i)%X_nmin       = pd(list(i))%X_nmin
       data%phytos(i)%X_nmax       = pd(list(i))%X_nmax
       data%phytos(i)%R_nuptake    = pd(list(i))%R_nuptake/secs_per_day
       data%phytos(i)%k_nfix       = pd(list(i))%k_nfix  ! should be between 0 and 1.
       data%phytos(i)%R_nfix       = pd(list(i))%R_nfix/secs_per_day
       data%phytos(i)%simDIPUptake = pd(list(i))%simDIPUptake
       data%phytos(i)%simIPDynamics= pd(list(i))%simIPDynamics
       data%phytos(i)%P_0          = pd(list(i))%P_0
       data%phytos(i)%K_P          = pd(list(i))%K_P
       data%phytos(i)%X_pcon       = pd(list(i))%X_pcon
       data%phytos(i)%X_pmin       = pd(list(i))%X_pmin
       data%phytos(i)%X_pmax       = pd(list(i))%X_pmax
       data%phytos(i)%R_puptake    = pd(list(i))%R_puptake/secs_per_day
       data%phytos(i)%simSiUptake  = pd(list(i))%simSiUptake
       data%phytos(i)%Si_0         = pd(list(i))%Si_0
       data%phytos(i)%K_Si         = pd(list(i))%K_Si
       data%phytos(i)%X_sicon      = pd(list(i))%X_sicon

       data%phytos(i)%c1           = pd(list(i))%c1
       data%phytos(i)%c3           = pd(list(i))%c3
       data%phytos(i)%f1           = pd(list(i))%f1
       data%phytos(i)%f2           = pd(list(i))%f2
       data%phytos(i)%d_phy        = pd(list(i))%d_phy

    ENDDO
    DEALLOCATE(pd)
END SUBROUTINE aed_phytoplankton_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_define_phyto_abm(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with AED.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_phyto_abm_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER  :: status

   INTEGER            :: num_phytos = 0
   INTEGER            :: the_phytos(MAX_PHYTO_TYPES) = 0
   INTEGER            :: settling(MAX_PHYTO_TYPES)     = _MOB_CONST_
   AED_REAL           :: resuspension(MAX_PHYTO_TYPES) = 0.
   CHARACTER(len=64)  :: resus_link='NCS_resus'
   CHARACTER(len=64)  :: phyto_particle_link=''            !   For FV API 2.0 (To be implemented)
   CHARACTER(len=64)  :: p_excretion_target_variable='OGM_dop'
   CHARACTER(len=64)  :: p_mortality_target_variable='OGM_pop'
   CHARACTER(len=64)  :: p1_uptake_target_variable='PHS_frp'
   CHARACTER(len=64)  :: p2_uptake_target_variable=''
   CHARACTER(len=64)  :: n_excretion_target_variable='OGM_don'
   CHARACTER(len=64)  :: n_mortality_target_variable='OGM_pon'
   CHARACTER(len=64)  :: n1_uptake_target_variable='NIT_nit'
   CHARACTER(len=64)  :: n2_uptake_target_variable='NIT_amm'
   CHARACTER(len=64)  :: n3_uptake_target_variable=''
   CHARACTER(len=64)  :: n4_uptake_target_variable=''
   CHARACTER(len=64)  :: c_excretion_target_variable='OGM_doc'
   CHARACTER(len=64)  :: c_mortality_target_variable='OGM_poc'
   CHARACTER(len=64)  :: c_uptake_target_variable=''
   CHARACTER(len=64)  :: do_uptake_target_variable='OXY_oxy'
   CHARACTER(len=64)  :: si_excretion_target_variable=''
   CHARACTER(len=64)  :: si_mortality_target_variable=''
   CHARACTER(len=64)  :: si_uptake_target_variable=''
   CHARACTER(len=128) :: dbase='aed_phyto_pars.csv'
   AED_REAL           :: zerolimitfudgefactor = 0.9 * 3600
   AED_REAL           :: R_mpbg = 0.
   AED_REAL           :: R_mpbr = 0.
   AED_REAL           :: R_mpbb = 0.
   AED_REAL           :: I_Kmpb = 100.
   AED_REAL           :: mpb_max = 1000.
   AED_REAL           :: theta_mpb_growth = 1.05
   AED_REAL           :: theta_mpb_resp   = 1.05
   AED_REAL           :: min_rho = 900.
   AED_REAL           :: max_rho = 1200.
   INTEGER            :: do_mpb = 0
   INTEGER            :: n_zones = 0
   AED_REAL           :: active_zones(1000) = 0


!  %% NAMELIST    %%  /aed_phyto_abm/
!  %% Last Checked 20/08/2021
   AED_REAL :: vvel_new = 0.
   AED_REAL :: vvel_old = 0.
   AED_REAL :: decay_rate_new = 0.
   AED_REAL :: decay_rate_old = 0.
   AED_REAL :: mass_limit = 10.
   AED_REAL :: X_cdw = 0.5
   AED_REAL :: X_nc = 0.1
   AED_REAL :: X_pc = 0.01
   AED_REAL :: X_dwww = 1.0

!  From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST    %%  /aed_phyto_abm/



   NAMELIST /aed_phyto_abm/ num_phytos, the_phytos, settling, resuspension,&
                    p_excretion_target_variable,p_mortality_target_variable,   &
                     p1_uptake_target_variable, p2_uptake_target_variable,     &
                    n_excretion_target_variable,n_mortality_target_variable,   &
                     n1_uptake_target_variable,n2_uptake_target_variable,      &
                     n3_uptake_target_variable,n4_uptake_target_variable,      &
                    c_excretion_target_variable,c_mortality_target_variable,   &
                      c_uptake_target_variable, do_uptake_target_variable,     &
                    si_excretion_target_variable,si_mortality_target_variable, &
                      si_uptake_target_variable,                               &
                    dbase, zerolimitfudgefactor,     &
                    do_mpb, R_mpbg, R_mpbr, I_Kmpb, mpb_max, min_rho, max_rho, &
                    resus_link, n_zones, active_zones, diag_level,             &
                    theta_mpb_growth,theta_mpb_resp,                           &
                    phyto_particle_link, R_mpbb,                               &
                    vvel_new, vvel_old, &
                    decay_rate_new, decay_rate_old, &
                    mass_limit, &
                    X_dwww, X_cdw, X_nc, X_pc
!-------------------------------------------------------------------------------
!BEGIN

   ! Initialise

   ! Read the namelist
   read(namlst,nml=aed_phyto_abm,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_phyto_abm'

   print *,"        aed_phyto_abm initialization"

   IF ( num_phytos >1 ) THEN
     print *," STOPPING: aed_phyto_abm only capable of 1 phyto at this time"
     STOP
   ENDIF

   ! Set module parameters
   data%vvel_new       = vvel_new/secs_per_day
   data%vvel_old       = vvel_old/secs_per_day
   data%decay_rate_new = decay_rate_new/secs_per_day
   data%decay_rate_old = decay_rate_old/secs_per_day
   data%mass_limit     = mass_limit
   data%X_dwww         = X_dwww
   data%X_cdw          = X_cdw
   data%X_nc           = X_nc
   data%X_pc           = X_pc


   ! Store group parameter values in the phyto derived type
   CALL aed_phytoplankton_load_params(data,dbase,                                &
                                      num_phytos,the_phytos,settling,resuspension)



   ! Register particle state variables (phytopalnkton agents)
   data%ip_ptm_c = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_c', 'mmol C/particle', 'particle C concentration')

   ! Diagnostic outputs for particle properties
   data%id_ptm_00 = aed_define_diag_variable('total_count', '#', 'particle count')
   data%id_ptm_14 = aed_define_diag_variable('total_vvel', 'm/s', 'sum of particle vvel')
   data%id_ptm_15 = aed_define_diag_variable('total_mass', 'g', 'sum of particle mass')
   data%id_ptm_17 = aed_define_diag_variable('total_birth', 'day', 'sum of birth date')
   data%id_ptm_18 = aed_define_diag_variable('total_age', 'days', 'sum of particle age')
   data%id_ptm114 = aed_define_diag_variable('vvel', 'm/s', 'last particle vvel')
   data%id_ptm115 = aed_define_diag_variable('mass', 'g', 'last particle mass')
   data%id_ptm117 = aed_define_diag_variable('birth', 'day', 'last particle birth time')
   data%id_ptm118 = aed_define_diag_variable('age', 'days', 'last particle age')

   ! Junk properties
   IF( diag_level >= 10 ) THEN
    data%id_ptm_01 = aed_define_diag_variable('bioptm01', '', 'bio_particles 01')
    data%id_ptm_02 = aed_define_diag_variable('bioptm02', '', 'bio_particles 02')
    data%id_ptm_03 = aed_define_diag_variable('bioptm03', '', 'bio_particles 03')
    data%id_ptm_04 = aed_define_diag_variable('bioptm04', '', 'bio_particles 04')
    data%id_ptm_05 = aed_define_diag_variable('bioptm05', '', 'bio_particles 05')
    data%id_ptm_06 = aed_define_diag_variable('bioptm06', '', 'bio_particles 06')
    data%id_ptm_07 = aed_define_diag_variable('bioptm07', '', 'bio_particles 07')
    data%id_ptm_08 = aed_define_diag_variable('bioptm08', '', 'bio_particles 08')
    data%id_ptm_09 = aed_define_diag_variable('bioptm09', '', 'bio_particles 09')
    data%id_ptm_10 = aed_define_diag_variable('bioptm10', '', 'bio_particles 10')
    data%id_ptm_11 = aed_define_diag_variable('bioptm11', '', 'bio_particles 11')
    data%id_ptm_12 = aed_define_diag_variable('bioptm12', '', 'bio_particles 12')
    data%id_ptm_13 = aed_define_diag_variable('bioptm13', '', 'bio_particles 13')
    data%id_ptm_16 = aed_define_diag_variable('bioptm16', '', 'bio_particles 16')

    data%id_ptm101 = aed_define_diag_variable('bioptm101', '', 'bio_particles101')
    data%id_ptm102 = aed_define_diag_variable('bioptm102', '', 'bio_particles102')
    data%id_ptm103 = aed_define_diag_variable('bioptm103', '', 'bio_particles103')
    data%id_ptm104 = aed_define_diag_variable('bioptm104', '', 'bio_particles104')
    data%id_ptm105 = aed_define_diag_variable('bioptm105', '', 'bio_particles105')
    data%id_ptm106 = aed_define_diag_variable('bioptm106', '', 'bio_particles106')
    data%id_ptm107 = aed_define_diag_variable('bioptm107', '', 'bio_particles107')
    data%id_ptm108 = aed_define_diag_variable('bioptm108', '', 'bio_particles108')
    data%id_ptm109 = aed_define_diag_variable('bioptm109', '', 'bio_particles109')
    data%id_ptm110 = aed_define_diag_variable('bioptm110', '', 'bio_particles110')
    data%id_ptm111 = aed_define_diag_variable('bioptm111', '', 'bio_particles111')
    data%id_ptm112 = aed_define_diag_variable('bioptm112', '', 'bio_particles112')
    data%id_ptm113 = aed_define_diag_variable('bioptm113', '', 'bio_particles113')
    data%id_ptm116 = aed_define_diag_variable('bioptm116', '', 'bio_particles116')
   ENDIF

   ! Diagnostic outputs for fluxes into cell
   data%id_d_oxy = aed_define_diag_variable('oxy_flux', 'mmol O2/m3/day','oxygen consumption')
   data%id_d_dc  = aed_define_diag_variable('dc_flux', 'mmol DOC/m3/day','dissolved C release')
   data%id_d_dn  = aed_define_diag_variable('dn_flux', 'mmol N/m3/day','dissolved N release')
   data%id_d_dp  = aed_define_diag_variable('dp_flux', 'mmol P/m3/day','dissolved P release')

   ! Linked state variables
   data%id_oxy = aed_locate_variable('OXY_oxy')
   data%id_amm = aed_locate_variable('NIT_amm')
   data%id_nit = aed_locate_variable('NIT_nit')
   data%id_frp = aed_locate_variable('PHS_frp')
   data%id_doc = aed_locate_variable('OGM_doc')
   data%id_don = aed_locate_variable('OGM_don')
   data%id_dop = aed_locate_variable('OGM_dop')

   ! Environment variables
   data%id_tem = aed_locate_global('temperature')
   data%id_lht = aed_locate_global('layer_ht')
   data%id_par = aed_locate_global('par')
   data%id_larea = aed_locate_sheet_global('layer_area')
   data%id_dep = aed_locate_sheet_global('col_depth')

END SUBROUTINE aed_define_phyto_abm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_particle_bgc_phyto_abm( data,column,layer_idx,ppid,p )
!ARGUMENTS
   CLASS (aed_phyto_abm_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   INTEGER,INTENT(inout) :: ppid
   TYPE (aed_ptm_t),INTENT(inout) :: p
!
!LOCALS
   INTEGER :: n
   AED_REAL :: oxy_flux
   AED_REAL :: decay, area, thickness

   AED_REAL :: Mu_max, WaterTemperature, Light, Depth, Kd, N_Limitation, P_Limitation, D0, D1
   AED_REAL :: f_T, Iz, f_I, Respiration, Mu_net

   AED_REAL, PARAMETER :: buoyancy_age = 86400.
   AED_REAL, PARAMETER :: DT           = 15.*60.
   AED_REAL, PARAMETER :: decay_rate   = 0.15
   AED_REAL, PARAMETER :: fres         = 0.7

!
!-------------------------------------------------------------------------------
!BEGIN

  print *, 'Phyto ABM: ',ppid
  
   ! Check if we are in a new cell, to reset cumulative counters
   IF (ppid == 0) THEN
      _DIAG_VAR_(data%id_ptm_14) = zero_
      _DIAG_VAR_(data%id_ptm_15) = zero_
      _DIAG_VAR_(data%id_ptm_17) = zero_
      _DIAG_VAR_(data%id_ptm_18) = zero_
      _DIAG_VAR_(data%id_d_oxy) = zero_
      _DIAG_VAR_(data%id_d_dc)  = zero_
      _DIAG_VAR_(data%id_d_dn)  = zero_
      _DIAG_VAR_(data%id_d_dp)  = zero_

      IF( diag_level >= 10 ) THEN
       _DIAG_VAR_(data%id_ptm_01) = zero_
       _DIAG_VAR_(data%id_ptm_02) = zero_
       _DIAG_VAR_(data%id_ptm_03) = zero_
       _DIAG_VAR_(data%id_ptm_04) = zero_
       _DIAG_VAR_(data%id_ptm_05) = zero_
       _DIAG_VAR_(data%id_ptm_06) = zero_
       _DIAG_VAR_(data%id_ptm_07) = zero_
       _DIAG_VAR_(data%id_ptm_08) = zero_
       _DIAG_VAR_(data%id_ptm_09) = zero_
       _DIAG_VAR_(data%id_ptm_10) = zero_
       _DIAG_VAR_(data%id_ptm_11) = zero_
       _DIAG_VAR_(data%id_ptm_12) = zero_
       _DIAG_VAR_(data%id_ptm_13) = zero_
       _DIAG_VAR_(data%id_ptm_16) = zero_
      ENDIF
    ENDIF

    ! Increment the particle count for this cell and set to diagnostic
    ppid = ppid + 1
    _DIAG_VAR_(data%id_ptm_00) = ppid !,AED_REAL)   ! total number of particles within a cell

   ! Temporary settings
   Mu_max=1.2  !maximum daily growth rate
   Light=2000.  !surface irradiance
   Kd=1.5  !light extinction coefficient
   N_Limitation=0.8  !limitation by N
   P_Limitation=0.8  !limitation by P
   D0= p%ptm_env(DIAM) !10.  !initial size

   ! Local environmental conditions in this layer
   WaterTemperature= _STATE_VAR_(data%id_tem) !22  !water temperature
   Depth     = _STATE_VAR_S_(data%id_dep) -  _PTM_ENV_(HGHT)  !cyanobacteria depth = water depth-cell height
   thickness = _STATE_VAR_(data%id_lht)
   area      = 1000. !_STATE_VAR_S_(data%id_larea)

   

   print *,'cell depth & temp',Depth, WaterTemperature

   !print *,'p%ptm_env(5),',p%ptm_env(5)

   !print *,'normal_sample(),',normal_sample(data%phytos(1)%T_opt, 3.)
   

   ! Net photosynthesis of cells
   f_T = exp(-((WaterTemperature - 22.) / 5.)**2) ! %temperature limitation term
   Iz = Light * exp(-Kd * Depth)  !Lambert-Beer's law of exponential light extinction
   f_I = (0.219 * Iz) / (0.219 * Iz + 25. + 0.001 * (0.219 * Iz)**2.)  !light limitation term
   Respiration = 0.1 * 1.1**(WaterTemperature - 20.)  !respiration
   Mu_net = Mu_max * f_T * f_I * min(N_Limitation, P_Limitation) - Respiration    !net daily growth rate
   D1 = D0 * 2.**(1. / (log10(2.) / Mu_net * 24.))  !predicted Dolichospermum size

   !print *, ' D_0: ', D0 !disp(['D_0: ', num2str(D0), ' μm'])
   !print *, ' D_1: ', D1 !disp(['D_1: ', num2str(D1), ' μm'])

   _PTM_ENV_(DIAM) = D1                  ! Set particle diameter

   _PTM_VAR_(data%ip_ptm_c) = 0.2       ! Set particle diameter

   print *, ' _PTM_VAR_: ', _PTM_VAR_(data%ip_ptm_c)

   ! Set interactions/fluxes with water properties
   oxy_flux = data%X_dwww * (1e3/12.) * (Mu_net/DT) * data%X_cdw / (area*thickness)  ! mmol C / m3/ s ! CHECK UNITS

   _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) - oxy_flux
   _FLUX_VAR_(data%id_amm) = _FLUX_VAR_(data%id_amm) + oxy_flux * data%X_nc
   _FLUX_VAR_(data%id_frp) = _FLUX_VAR_(data%id_frp) + oxy_flux * data%X_pc
   _FLUX_VAR_(data%id_nit) = _FLUX_VAR_(data%id_nit) + zero_

   ! DOM leakage from particles during decay
  ! _FLUX_VAR_(data%id_doc) = _FLUX_VAR_(data%id_doc) + (1.-fres) * oxy_flux
  ! _FLUX_VAR_(data%id_don) = _FLUX_VAR_(data%id_don) + (1.-fres) * oxy_flux * data%X_nc
  ! _FLUX_VAR_(data%id_dop) = _FLUX_VAR_(data%id_dop) + (1.-fres) * oxy_flux * data%X_pc

   ! Set diagnostics to track cumulative oxygen and nutrient fluxes into a cell
   _DIAG_VAR_(data%id_d_oxy) = _DIAG_VAR_(data%id_d_oxy) - oxy_flux * secs_per_day    ! O2
  ! _DIAG_VAR_(data%id_d_dc)  = _DIAG_VAR_(data%id_d_dc) - (1.-fres)* oxy_flux * secs_per_day ! DOC
  ! _DIAG_VAR_(data%id_d_dn)  = &
  !                          _DIAG_VAR_(data%id_d_dn) - oxy_flux * data%X_nc * secs_per_day    ! DON + NH4 + NO3
  ! _DIAG_VAR_(data%id_d_dp)  = _DIAG_VAR_(data%id_d_dp) - oxy_flux * data%X_pc * secs_per_day ! DOP + FRP

   ! Update particle bouyancy, changing with age
   p%ptm_env(VVEL) = -1.0/86400.

   ! Set general diagnostics, summarising particles in this cell

   ! 1st, Cumulate particle properties (needs to be divided by particle number for average)
   _DIAG_VAR_(data%id_ptm_14) = _DIAG_VAR_(data%id_ptm_14) + p%ptm_env(VVEL)
   _DIAG_VAR_(data%id_ptm_15) = &
                  _DIAG_VAR_(data%id_ptm_15) + p%ptm_env(MASS) ! total particle mass within a cell
   _DIAG_VAR_(data%id_ptm_17) = _DIAG_VAR_(data%id_ptm_17) !+ partcl(PTM_BIRTH)
   _DIAG_VAR_(data%id_ptm_18) = &
                  _DIAG_VAR_(data%id_ptm_18) !+ (partcl(PTM_AGE)-partcl(PTM_BIRTH)) /secs_per_day
   IF( diag_level >= 10 ) THEN
    _DIAG_VAR_(data%id_ptm_01) = _DIAG_VAR_(data%id_ptm_01) + 0. !partcl(1)
    _DIAG_VAR_(data%id_ptm_02) = _DIAG_VAR_(data%id_ptm_02) + 0. !partcl(2)
    _DIAG_VAR_(data%id_ptm_03) = _DIAG_VAR_(data%id_ptm_03) + 0. !partcl(3)
    _DIAG_VAR_(data%id_ptm_04) = _DIAG_VAR_(data%id_ptm_04) + 0. !partcl(4)
    _DIAG_VAR_(data%id_ptm_05) = _DIAG_VAR_(data%id_ptm_05) + 0. !partcl(5)
    _DIAG_VAR_(data%id_ptm_06) = _DIAG_VAR_(data%id_ptm_06) + 0. !partcl(6)
    _DIAG_VAR_(data%id_ptm_07) = _DIAG_VAR_(data%id_ptm_07) + 0. !partcl(7)
    _DIAG_VAR_(data%id_ptm_08) = _DIAG_VAR_(data%id_ptm_08) + 0. !partcl(8)
    _DIAG_VAR_(data%id_ptm_09) = _DIAG_VAR_(data%id_ptm_09) + 0. !partcl(9)
    _DIAG_VAR_(data%id_ptm_10) = _DIAG_VAR_(data%id_ptm_10) + 0. !partcl(10)
    _DIAG_VAR_(data%id_ptm_11) = _DIAG_VAR_(data%id_ptm_11) + 0. !partcl(11)
    _DIAG_VAR_(data%id_ptm_12) = _DIAG_VAR_(data%id_ptm_12) + 0. !partcl(12)
    _DIAG_VAR_(data%id_ptm_13) = _DIAG_VAR_(data%id_ptm_13) + 0. !partcl(13)
    _DIAG_VAR_(data%id_ptm_16) = _DIAG_VAR_(data%id_ptm_16) + 0. !partcl(16)
   ENDIF

   ! 2nd, Set particle property (this will therefore remember last particle only)
   _DIAG_VAR_(data%id_ptm114) = p%ptm_env(VVEL)
   _DIAG_VAR_(data%id_ptm115) = 0. !partcl(PTM_MASS)
   _DIAG_VAR_(data%id_ptm117) = 0. !partcl(PTM_BIRTH)
   _DIAG_VAR_(data%id_ptm118) = 0. !(partcl(PTM_AGE)-partcl(PTM_BIRTH))/secs_per_day
   IF( diag_level >= 10 ) THEN
    _DIAG_VAR_(data%id_ptm101) = 0. !partcl(1)
    _DIAG_VAR_(data%id_ptm102) = 0. !partcl(2)
    _DIAG_VAR_(data%id_ptm103) = 0. !partcl(3)
    _DIAG_VAR_(data%id_ptm104) = 0. !partcl(4)
    _DIAG_VAR_(data%id_ptm105) = 0. !partcl(5)
    _DIAG_VAR_(data%id_ptm106) = 0. !partcl(6)
    _DIAG_VAR_(data%id_ptm107) = 0. !partcl(7)
    _DIAG_VAR_(data%id_ptm108) = 0. !partcl(8)
    _DIAG_VAR_(data%id_ptm109) = 0. !partcl(9)
    _DIAG_VAR_(data%id_ptm110) = 0. !partcl(10)
    _DIAG_VAR_(data%id_ptm111) = 0. !partcl(11)
    _DIAG_VAR_(data%id_ptm112) = 0. !partcl(12)
    _DIAG_VAR_(data%id_ptm113) = 0. !partcl(13)
    _DIAG_VAR_(data%id_ptm116) = 0. !partcl(16)
   ENDIF


END SUBROUTINE aed_particle_bgc_phyto_abm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_phyto_abm
