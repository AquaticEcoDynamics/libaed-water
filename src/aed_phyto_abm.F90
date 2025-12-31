!###############################################################################
!#                                                                             #
!# aed_phyto_abm.F90                                                           #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2018-2025 - The University of Western Australia                  #
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

   USE trait_functions
   USE gmk
   USE mGf90

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
      INTEGER :: id_count
      INTEGER :: id_phyc, id_phyn, id_phyp, id_chl, id_mtopt, id_vtopt, &
                 id_mcdiv, id_vcdiv, id_mlnalpha, id_vlnalpha, &
                 id_cov_ta, id_cov_tl, id_cov_al, id_N_birth,  &
                 id_N_mutate, id_IPAR, id_NPPc, id_cells,      &
                 id_N_death
      INTEGER :: ip_c, ip_n, ip_p, ip_par, ip_tem, ip_no3, ip_frp, ip_chl, ip_num, ip_cdiv, ip_Topt, ip_LnalphaChl, ip_mu_C, ip_fit
      INTEGER :: id_d_oxy, id_d_dc, id_d_dn, id_d_dp, id_d_nit, id_d_pon, id_d_frp, id_d_pop, id_d_poc
      INTEGER :: id_oxy,id_amm,id_nit,id_frp,id_doc,id_don,id_dop,id_poc,id_pon,id_pop
      INTEGER :: id_lht, id_larea, id_dep, id_tem, id_par, id_I0, id_dens, id_yday, id_depth

      AED_REAL :: vvel_new, vvel_old, decay_rate_new, decay_rate_old
      AED_REAL :: X_dwww, X_cdw, X_nc, X_pc, mass_limit

      ! Phytoplankton parameters
      INTEGER  :: num_phytos
      TYPE(phyto_data_t),DIMENSION(:),ALLOCATABLE :: phytos

      CONTAINS
         PROCEDURE :: define              => aed_define_phyto_abm
!        PROCEDURE :: calculate           => aed_calculate_phyto_abm
!        PROCEDURE :: calculate_benthic   => aed_calculate_benthic_phyto_abm
!        PROCEDURE :: calculate_riparian  => aed_calculate_riparian_phyto_abm
!        PROCEDURE :: calculate_dry       => aed_calculate_dry_phyto_abm
!        PROCEDURE :: equilibrate         => aed_equilibrate_phyto_abm
         PROCEDURE :: initialize_particle => aed_particle_initialize_phyto_abm
         PROCEDURE :: particle_bgc        => aed_particle_bgc_phyto_abm
!        PROCEDURE :: mobility            => aed_mobility_phyto_abm
         PROCEDURE :: light_extinction    => aed_light_extinction_phyto_abm
!        PROCEDURE :: delete              => aed_delete_phyto_abm

   END TYPE

   INTEGER, PARAMETER :: PTM_MASS   = 15    ! TBD
   INTEGER, PARAMETER :: PTM_VVEL   = 14    ! TBD
   INTEGER, PARAMETER :: PTM_BIRTH  = 17    ! TBD
   INTEGER, PARAMETER :: PTM_AGE    = 18    ! TBD
   INTEGER, PARAMETER :: PTM_STATUS = 19    ! TBD

   ! PTM array reference index locations
   INTEGER, PARAMETER :: STAT = 1  !#define STAT   0 _PTM_STAT_
   INTEGER, PARAMETER :: IDX2 = 2  !#define IDX2   1 _PTM_STAT_
   INTEGER, PARAMETER :: IDX3 = 3  !#define IDX3   2 _PTM_STAT_
   INTEGER, PARAMETER :: LAYR = 4  !#define LAYR   3 _PTM_STAT_
   INTEGER, PARAMETER :: FLAG = 5  !#define FLAG   4 _PTM_STAT_
   INTEGER, PARAMETER :: PTID = 6  !#define PTID   5 _PTM_STAT_

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

            ! added with addition of PIBM code
            CASE ('X_cinit')       ; pd(dcol)%X_cinit       = extract_double(values(ccol))
            CASE ('X_ninit')       ; pd(dcol)%X_ninit       = extract_double(values(ccol))
            CASE ('X_pinit')       ; pd(dcol)%X_pinit       = extract_double(values(ccol))
            CASE ('X_chlinit')     ; pd(dcol)%X_chlinit     = extract_double(values(ccol))
            CASE ('Cdiv')          ; pd(dcol)%Cdiv          = extract_double(values(ccol))
            CASE ('n0')            ; pd(dcol)%n0            = extract_double(values(ccol))
            CASE ('Lnalphachl')    ; pd(dcol)%Lnalphachl    = extract_double(values(ccol))
            CASE ('RC')            ; pd(dcol)%RC            = extract_double(values(ccol))
            CASE ('RN')            ; pd(dcol)%RN            = extract_double(values(ccol))
            CASE ('RP')            ; pd(dcol)%RP            = extract_double(values(ccol))
            CASE ('RChl')          ; pd(dcol)%RChl          = extract_double(values(ccol))
            CASE ('zeta_N')        ; pd(dcol)%zeta_N        = extract_double(values(ccol))
            CASE ('zeta_P')        ; pd(dcol)%zeta_P        = extract_double(values(ccol))
            CASE ('a1')            ; pd(dcol)%a1            = extract_double(values(ccol))
            CASE ('mort_prob')     ; pd(dcol)%mort_prob     = extract_double(values(ccol))
            CASE ('nx')            ; pd(dcol)%nx            = extract_double(values(ccol))
            CASE ('thetaNmax')     ; pd(dcol)%thetaNmax     = extract_double(values(ccol))
            CASE ('QNmin_a')       ; pd(dcol)%QNmin_a       = extract_double(values(ccol))
            CASE ('QNmin_b')       ; pd(dcol)%QNmin_b       = extract_double(values(ccol))
            CASE ('QNmax_a')       ; pd(dcol)%QNmax_a       = extract_double(values(ccol))
            CASE ('QNmax_b')       ; pd(dcol)%QNmax_b       = extract_double(values(ccol))
            CASE ('QPmin_a')       ; pd(dcol)%QPmin_a       = extract_double(values(ccol))
            CASE ('QPmin_b')       ; pd(dcol)%QPmin_b       = extract_double(values(ccol))
            CASE ('QPmax_a')       ; pd(dcol)%QPmax_a       = extract_double(values(ccol))
            CASE ('QPmax_b')       ; pd(dcol)%QPmax_b       = extract_double(values(ccol))
            CASE ('KN_a')          ; pd(dcol)%KN_a          = extract_double(values(ccol))
            CASE ('KN_b')          ; pd(dcol)%KN_b          = extract_double(values(ccol))
            CASE ('KPho_a')        ; pd(dcol)%KPho_a        = extract_double(values(ccol))
            CASE ('KPho_b')        ; pd(dcol)%KPho_b        = extract_double(values(ccol))


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
       data%phytos(i)%settling     = settling(i)
      !data%phytos(i)%resuspension = resuspension(i)
       data%phytos(i)%Xcc          = pd(list(i))%Xcc
       data%phytos(i)%R_growth     = pd(list(i))%R_growth!ML/secs_per_day not how it's used in PIBM
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

       ! added with addition of PIBM code
       data%phytos(i)%X_cinit       = pd(list(i))%X_cinit
       data%phytos(i)%X_ninit       = pd(list(i))%X_ninit
       data%phytos(i)%X_pinit       = pd(list(i))%X_pinit
       data%phytos(i)%X_chlinit     = pd(list(i))%X_chlinit
       data%phytos(i)%Cdiv          = pd(list(i))%Cdiv
       data%phytos(i)%n0            = pd(list(i))%n0
       data%phytos(i)%Lnalphachl    = pd(list(i))%Lnalphachl
       data%phytos(i)%RC            = pd(list(i))%RC
       data%phytos(i)%RN            = pd(list(i))%RN
       data%phytos(i)%RP            = pd(list(i))%RP
       data%phytos(i)%RChl          = pd(list(i))%RChl
       data%phytos(i)%zeta_N        = pd(list(i))%zeta_N
       data%phytos(i)%zeta_P        = pd(list(i))%zeta_P
       data%phytos(i)%a1            = pd(list(i))%a1
       data%phytos(i)%mort_prob     = pd(list(i))%mort_prob
       data%phytos(i)%nx            = pd(list(i))%nx
       data%phytos(i)%thetaNmax     = pd(list(i))%thetaNmax
       data%phytos(i)%QNmin_a       = pd(list(i))%QNmin_a
       data%phytos(i)%QNmin_b       = pd(list(i))%QNmin_b
       data%phytos(i)%QNmax_a       = pd(list(i))%QNmax_a
       data%phytos(i)%QNmax_b       = pd(list(i))%QNmax_b
       data%phytos(i)%QPmin_a       = pd(list(i))%QPmin_a
       data%phytos(i)%QPmin_b       = pd(list(i))%QPmin_b
       data%phytos(i)%QPmax_a       = pd(list(i))%QPmax_a
       data%phytos(i)%QPmax_b       = pd(list(i))%QPmax_b
       data%phytos(i)%KN_a          = pd(list(i))%KN_a
       data%phytos(i)%KN_b          = pd(list(i))%KN_b
       data%phytos(i)%KPho_a        = pd(list(i))%KPho_a
       data%phytos(i)%KPho_b        = pd(list(i))%KPho_b

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



   ! Register particle state variables (phytoplankton agents)

   !real    :: PAR = 0.01  associated PAR                                       _PTM_VAR_          p(..)%ptm_state(data%ip_par)
   data%ip_par = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_par', 'ummol m2 sec', 'particle layer PAR')

   !real    :: temp= 20.   associated Temperature                               _PTM_VAR_
   data%ip_tem = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_tem', 'degrees C', 'particle layer temperature')

   !real    :: NO3 = 0.1   associated nutrient                                  _PTM_VAR_
   data%ip_no3 = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_no3', 'mmol N', 'particle layer NO3')

   !real    :: FRP = 0.1   associated nutrient                                  _PTM_VAR_
   data%ip_frp = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_frp', 'mmol P', 'particle layer FRP')

   !real    :: C   = 0.02          ! cellular carbon content (pmol; assuming a 1 micron cell)
   data%ip_c = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_c', 'pmol C/cell', 'cell C concentration',initial=data%phytos(1)%X_cinit) ! 0.02
   
   !real    :: N   = 0.02/106.*16. ! cellular nitrogen content (pmol per cell; assuming a 1 micron cell)
   data%ip_n = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_n', 'pmol N/cell', 'cell N concentration',initial=data%phytos(1)%X_ninit) !0.02/106.*16.

   !real    :: P   = 0.02/106.*1. ! cellular phosphorus content (pmol per cell; assuming a 1 micron cell)
   data%ip_p = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_p', 'pmol P/cell', 'cell P concentration',initial=data%phytos(1)%X_pinit) !0.02/106.*1.

   !real    :: Chl = 0.02 * 12/50  ! Cellular Chl content (pg Chl)
   data%ip_chl = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_chl', 'pg Chl', 'cell Chl concentration',initial = data%phytos(1)%X_chlinit) !0.02 * 12/50

   !real    :: num = 5d9           ! Number of cells per superindividual
   data%ip_num = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_num', 'number', 'number of cells/particle',initial = data%phytos(1)%n0) !5d12

   !real :: Cdiv= 0.04d0           !cellular carbon content threshold for division (pmol/cell), can be used as a proxy for size and can be converted to ESD; Phytoplankton half-saturation constant, minimal N:C and maximal N:C ratios are allometric functions of this parameter
                                   !This trait will vary with mutation
   data%ip_Cdiv = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_cdiv', 'pmol/cell', 'cellular carbon content threshold for division',initial = data%phytos(1)%Cdiv) ! 0.04d0

   !real :: Topt = 20.d0           !Optimal temperature
   data%ip_Topt = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_topt', 'degrees C', 'optimal temperature',initial = data%phytos(1)%T_opt)

   !real :: LnalphaChl = -2.3  !log(0.1) !Ln alphaChl (slope of the P-I curve; unit: (W m-2)-1 (gChl molC)-1 d-1 instead of micro mol quanta m-2 s-1)
   data%ip_LnalphaChl = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_lnalphachl', '(W m-2)-1 (gChl molC)-1 d-1', 'slope of the P-I curve',initial = data%phytos(1)%Lnalphachl) !-2.3

   !real :: mu_C = 0.d0              !Carbon-specific growth rate
   data%ip_mu_C = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_mu_c', '??', 'carbon-specific growth rate',initial = 0.d0)

   !real :: fitness = 0.d0     !Fitness represented by net growth rate (growth - mortality)
   data%ip_fit = aed_define_ptm_variable(TRIM(data%phytos(1)%p_name)//'_fit', '??', 'fitness represented by net growth rate (growth - mortality)',initial = 0.d0)



   ! Diagnostic outputs for particle properties
   data%id_count  = aed_define_diag_variable('total_count', '#', 'layer live particle count')             !#MH  N_ in PIBM

   !real      :: PHYC(nlev) = 0d0  !Eulerian concentration of phyto C, mean trait  for each layer
   data%id_phyc = aed_define_diag_variable('mean_C', 'mmol/m3', 'mean Eulerian concentration of phyto particle C', rezero = .false.)

   !real      ::    PHY(nlev) = 0d0
   data%id_phyn = aed_define_diag_variable('mean_N', 'mmol/m3', 'mean Eulerian concentration of phyto particle N')

   !real      ::    PHY(nlev) = 0d0
   data%id_phyp = aed_define_diag_variable('mean_P', 'mmol/m3', 'mean Eulerian concentration of phyto particle P')

   !real      ::    CHL(nlev) = 0d0
   data%id_chl = aed_define_diag_variable('mean_Chl', 'mg/m3', 'mean Eulerian concentration of phyto particle Chl')

   !real      ::   mTopt_(nlev) = 0d0
   data%id_mTopt = aed_define_diag_variable('mean_Topt', 'degrees C/pmol C', 'carbon-weighted mean layer optimal temperature')

   !real      ::   vTopt_(nlev) = 0d0
   data%id_vTopt = aed_define_diag_variable('var_Topt', 'degrees C/pmol C', 'carbon-weighted layer optimal temperature variance')

   !real      ::   mCDiv_(nlev) = 0d0                                                                                                 !#ML check units here
   data%id_mCDiv = aed_define_diag_variable('mean_CDiv', 'log(pmol C/cell)/pmol C', 'carbon-weighted log mean layer cellular carbon content threshold for division')

   !real      ::   vCDiv_(nlev) = 0d0                                                                                                 !#ML check units here
   data%id_vCDiv = aed_define_diag_variable('var_CDiv', 'log(pmol C/cell)/pmol C', 'carbon-weighted log layer cellular carbon content threshold for division variance')

   !real      ::   mlnalpha_(nlev) = 0d0
   data%id_mlnalpha = aed_define_diag_variable('mean_lnalpha', 'log((W m-2)-1 (gChl molC)-1 d-1)/pmol C', 'carbon-weighted mean layer log slope of the P-I curve')

   !real      ::   vlnalpha_(nlev) = 0d0
   data%id_vlnalpha = aed_define_diag_variable('var_lnalpha', 'log((W m-2)-1 (gChl molC)-1 d-1)/pmol C', 'carbon-weighted layer log slope of the P-I curve variance')

   !real      ::   cov_TA(nlev) = 0d0
   data%id_cov_TA = aed_define_diag_variable('cov_TA', 'degrees C/pmol C * log(((W m-2)-1 (gChl molC)-1 d-1))/pmol C', 'covariance between Topt and ln(alpha Chl)')

   !real      ::   cov_TL(nlev) = 0d0
   data%id_cov_TL = aed_define_diag_variable('cov_TL', 'degrees C/pmol C * log(pmol/cell)/pmol C', 'covariance between Topt and Cdiv')

   !real      ::   cov_AL(nlev) = 0d0
   data%id_cov_AL = aed_define_diag_variable('cov_AL', 'log((W m-2)-1 (gChl molC)-1 d-1))/pmol C * log(pmol/cell)/pmol C', 'covariance between ln(alpha Chl) and Cdiv')

   !integer :: N_birth(nlev) = 0  !Number of birth events during one hour                    # ML this may cause problems because is defined as an integer
   data%id_N_birth = aed_define_diag_variable('N_birth', 'number', 'number of births in a layer')

   !integer :: N_mutate(nlev)= 0  !Number of mutation events during one hour at each grid    # ML this may cause problems because is defined as an integer
   data%id_N_mutate = aed_define_diag_variable('N_mutate', 'number', 'number of mutation events in a layer')

   !integer :: N_death(nlev) = 0  !Number of death events during one hour
   data%id_N_death = aed_define_diag_variable('N_death', 'number', 'number of deaths in a layer')

   ! Additional diagnostics from PIBM Geider_Lag.F90
   ! Includes diagnostic outputs for fluxes to and from layers
   data%id_IPAR = aed_define_diag_variable('IPAR', 'umol/m2/s', 'daily integrated PAR at each depth')
   data%id_NPPc = aed_define_diag_variable('NPPc', '(mg C m-3 d-1)', 'C-based phytoplankton production')
   data%id_cells = aed_define_diag_variable('id_cells', 'number', 'cells/m3')
   data%id_d_nit = aed_define_diag_variable('id_d_nit', 'mmol N/m3/day', 'daily flux of NO3 from particles in a layer')
   data%id_d_pon = aed_define_diag_variable('id_d_pon', 'mmol N/m3/day', 'daily flux of PON from particles in a layer')
   data%id_d_frp = aed_define_diag_variable('id_d_frp', 'mmol P/m3/day', 'daily flux of FRP from particles in a layer')
   data%id_d_pop = aed_define_diag_variable('id_d_pop', 'mmol P/m3/day', 'daily flux of POP from particles in a layer')
   data%id_d_poc = aed_define_diag_variable('id_d_poc', 'mmol C/m3/day', 'daily flux of POC from particles in a layer')
   data%id_d_oxy = aed_define_diag_variable('oxy_flux', 'mmol O2/m3/day','oxygen production')
   
   ! Linked state variables
   data%id_oxy = aed_locate_variable('OXY_oxy')
   data%id_amm = aed_locate_variable('NIT_amm')
   data%id_nit = aed_locate_variable('NIT_nit')
   data%id_frp = aed_locate_variable('PHS_frp')
   data%id_doc = aed_locate_variable('OGM_doc')
   data%id_don = aed_locate_variable('OGM_don')
   data%id_dop = aed_locate_variable('OGM_dop')
   data%id_poc = aed_locate_variable('OGM_poc')
   data%id_pon = aed_locate_variable('OGM_pon')
   data%id_pop = aed_locate_variable('OGM_pop')

   ! Environment variables
   data%id_depth = aed_locate_global('depth')
   data%id_tem   = aed_locate_global('temperature')
   data%id_dens  = aed_locate_global('density')
   data%id_lht   = aed_locate_global('layer_ht')
   data%id_par   = aed_locate_global('par')
   data%id_larea = aed_locate_global('layer_area') ! ML currently state var s; use this for now but is only top layer
                                                   ! CAB no! top layer is col_area
   data%id_dep   = aed_locate_sheet_global('col_depth')
   data%id_I0    = aed_locate_sheet_global('par_sf')
   data%id_yday  = aed_locate_sheet_global('yearday')

END SUBROUTINE aed_define_phyto_abm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_particle_initialize_phyto_abm( data,ppid,p )

use params,          only : NInit, sigma_init
use mGf90,           only : srand_mtGaus

!ARGUMENTS
   CLASS (aed_phyto_abm_data_t),INTENT(in) :: data
   TYPE(aed_ptm_t), DIMENSION(:), INTENT(inout) :: p
   INTEGER,INTENT(inout) :: ppid

!

!LOCALS
real      :: oldtt(1) = 0.   !Scratch variable for storing the old trait
real      :: newtt(1) = 0.   !Scratch variable for storing the new trait
real      :: vartt(1,1) = 0.   !Variance of the mutating trait
INTEGER   :: i, v

!

!Number of traits
integer, parameter :: iTopt = 1        !Trait index for Topt
integer, parameter :: iSize = 2        !Trait index for Cdiv: carbon-weighted cellular carbon content threshold for division
integer, parameter :: ialphaChl = 3    !Trait index for optimal light
integer, parameter :: iC = 4           !Trait index for internal C
integer, parameter :: iNO3 = 5         !Trait index for internal NO3
integer, parameter :: iP = 6           !Trait index for internal P
integer, parameter :: iChl = 7         !Trait index for internal chl
integer, parameter :: iNum = 8         !Trait index for number of cells per particle

!
!-------------------------------------------------------------------------------
!BEGIN
DO i = 1, ppid
   If (NInit > 0) Then
      DO v = 1, NInit
            select case(v)
            case(iTopt)
               oldtt(1) = p(i)%ptm_state(data%ip_Topt)
            case(iSize)
               oldtt(1) = log(p(i)%ptm_state(data%ip_Cdiv))
            case(ialphaChl)
               oldtt(1) = p(i)%ptm_state(data%ip_LnalphaChl)
            case(iC)
               oldtt(1) = log(p(i)%ptm_state(data%ip_c))
            case(iNO3)
               oldtt(1) = log(p(i)%ptm_state(data%ip_n))
            case(iP)
               oldtt(1) = log(p(i)%ptm_state(data%ip_p))
            case(iChl)
               oldtt(1) = log(p(i)%ptm_state(data%ip_chl))
            case(iNum)
               oldtt(1) = log(p(i)%ptm_state(data%ip_num))
            case DEFAULT
               stop "Trait index wrong!"
            end select

            vartt(1,1)= sigma_init(v)**2   !Construct the covariance matrix for the selected trait

            !A new Topt is randomly sampled from a Gaussian distribution with mean of previous Topt and SD of sigma
            newtt = srand_mtGaus(1, oldtt, vartt)
            select case(v)
            case(iTopt)
               p(i)%ptm_state(data%ip_Topt) = newtt(1)
            case(iSize)
               p(i)%ptm_state(data%ip_Cdiv) = exp(newtt(1))
            case(ialphaChl)
               p(i)%ptm_state(data%ip_LnalphaChl) = newtt(1)
            case(iC)
               p(i)%ptm_state(data%ip_c) = exp(newtt(1))
            case(iNO3)
               p(i)%ptm_state(data%ip_n) = exp(newtt(1))
            case(iP)
               p(i)%ptm_state(data%ip_p) = exp(newtt(1))
            case(iChl)
               p(i)%ptm_state(data%ip_chl) = exp(newtt(1))
            case(iNum)
               p(i)%ptm_state(data%ip_num) = exp(newtt(1))
            case DEFAULT
               stop "Trait index wrong!"
            end select
         !ENDIF
      ENDDO !End of looping through Traits
   ENDIF !End of if (NTrait > 0)
ENDDO ! End of looping through superindividuals


END SUBROUTINE aed_particle_initialize_phyto_abm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_particle_bgc_phyto_abm( data,column,layer_idx,ppid,p )

!DECLARATIONS FROM PIBM Par2PHY SUBROUTINE
   use params,          only : NTrait, nu, sigma
   use mGf90,           only : srand_mtGaus
   USE Trait_functions, only : TEMPBOL, PHY_C2Vol, palatability

   IMPLICIT NONE

!ARGUMENTS
   CLASS (aed_phyto_abm_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   INTEGER,INTENT(inout) :: ppid
   TYPE(aed_ptm_t), DIMENSION(:), INTENT(inout) :: p

!LOCALS
   INTEGER :: n, status
   AED_REAL :: oxy_flux
   AED_REAL :: decay, area, thickness

   AED_REAL :: Mu_max, WaterTemperature, Light, Depth, Kd, N_Limitation, P_Limitation, D0, D1, par, no3, frp, pw, mu, pw20, mu20
   AED_REAL :: f_T, Iz, f_I, Respiration, Mu_net

   AED_REAL, PARAMETER :: buoyancy_age = 86400.
   AED_REAL, PARAMETER :: DT           = 15.*60.
   AED_REAL, PARAMETER :: decay_rate   = 0.15
   AED_REAL, PARAMETER :: fres         = 0.7

! ADDITIONAL DECLARATIONS FROM PIBM
   real       :: PHYC = 0d0
   real       :: PHYN = 0d0
   real       :: PHYP = 0d0
   real       :: CHL = 0d0
   real      ::  nu_ = 0d0   !Basic Mutation rate
   real      :: cff = 0.d0   !Random number [0,1] for mutation
   real      :: rnm = 0.d0   !Random number [0,1] for mortality
   real      :: oldtt(1) = 0.   !Scratch variable for storing the old trait
   real      :: newtt(1) = 0.   !Scratch variable for storing the new trait
   real      :: vartt(1,1) = 0.   !Variance of the mutating trait
   integer :: k,i,m,ipar
   !Additional variables declared from PIBM variables.F90 file
   integer, parameter :: iTopt = 1        !Trait index for Topt
   integer, parameter :: iSize = 2        !Trait index for Size (ESD)
   integer, parameter :: ialphaChl = 3    !Trait index for optimal light
   !Additional variables declared from PIBM grid.F90 file
   AED_REAL :: Hz != 9962463.581  mean layer area of FCR
   ! More PIBM declarations
   integer, parameter :: nzoo = 0
   INTEGER :: j, zz, v
   INTEGER :: N_ = 0   !Number of particles in each grid
   real    :: Abun_ = 0   !Total abundance in each grid
   real    :: ZOO(NZOO) = 0.
   real    :: DET = 0.
   real    :: tf_z  = 0.
   real    :: Graz  = 0.
   real    :: dC_   = 0.
   real    :: dN_   = 0.
   real    :: dP_   = 0.
   real    :: dChl_ = 0.
   real    :: ESD_ = 0.    ! equivalent spherical diameter (micron)
   real    :: uptake= 0.   !Total NO3 uptake
   real    :: uptake_P = 0.   !Total FRP uptake
   real    :: pp_DZ = 0.
   real    :: pp_DZP = 0.
   real    :: pp_ND = 0.
   real    :: pp_PD = 0.
   real    :: Pmort = 0.
   real    :: Pmort_P = 0.
   real    :: Pmort_C = 0.
   real    :: Cmin  = 0. !Phytoplankton subsistence carbon quota below which the cell will die
   real    :: FZoo(NZOO) = 0.   !The total amount of palatable prey (in Nitrogen)
                             !for each zooplankton size class
   real    :: phyV = 0.   !Phytoplankton cell volume
   real    :: RES     = 0.
   real    :: RES_P   = 0.
   real    :: EGES  = 0.
   real    :: gbar  = 0.
   real    :: INGES(NZOO) = 0.
   real    :: Zmort = 0.
   real    :: Gmatrix(NZOO,NZOO) = 0.d0     !Grazer biomass specific grazing rate matrix
   real,    allocatable :: BN(:)            !The amount of nitrogen in each super-individual
   real,    allocatable :: BP(:)            !The amount of phosphorus in each super-individual
   real,    allocatable :: BC(:)            !The amount of carbon in each super-individual
   real,    allocatable :: Pmatrix(:,:)     !Phytoplankton mortality rates by each zooplankton size class for each superindividual
   real,    parameter   :: eta    = -1.d0*6.6  !Prey refuge parameter for nitrogen
   real,    parameter   :: A_g    = 21.9   !Intercept of the allometric equation of maximal zooplankton grazing rate (Ward et al. 2012)
   real,    parameter   :: B_g    = -0.16  !Slope of the allometric equation of maximal zooplankton grazing rate (Ward et al. 2012)
   real,    parameter   :: mz_g   = 0.d0  !Power of zooplankton mortality following Ward et al. (2013)
   real   :: Nt_min = 0.0 !Minimal amount nitrogen of each super-individual (number of cells * N content per cell)
   real   :: P_min = 0.0 !Minimal amount phosphorus of each super-individual (number of cells * P content per cell)
   integer:: N_PAR = 0
   INTEGER, ALLOCATABLE :: index_(:)    !The indexes of particles in each grid
   INTEGER, ALLOCATABLE :: scratch(:)    !The scratch indexes of particles in each grid
   INTEGER              :: Allocatestatus = 0
   real :: PHY_t = 0d0  !Total phytoplankton N
   ! additional declarations from PIBM time_settings.F90
   real     :: dtdays       = 1.0d0/24 !ML initializing this to 1 hour
   real  :: hour_of_day   != 0
   integer  :: integer_day   != 0
   real  :: real_day
   ! ML new declarations to get things working with PIBM and AED
   integer  :: chk_lyr_new  = 0.
   integer  :: chk_lyr_old  = 0.
   real     :: p_vvel = 0. !vertical velocity calculated in AED (m/s)
   real     :: min_rho = 900.
   real     :: max_rho = 1200.
   integer  :: phy_i = 1.
   real     :: dens_flux = 0.

!End of declaration
!
!-------------------------------------------------------------------------------
!BEGIN

!------------------------------------------------------
! INITIAL CODE FROM AED

   i=1 !ML need to remove this later when can handle more than one phyto group in ABM

   ! Check if we are in a new cell, to reset cumulative counters
   chk_lyr_new = layer_idx
   IF ((chk_lyr_new > chk_lyr_old) .or. chk_lyr_new == 1) THEN
      chk_lyr_old = chk_lyr_new
   ENDIF

   ! Local environmental conditions in this layer
   WaterTemperature= _STATE_VAR_(data%id_tem) !22  !water temperature

   ! PAR IS NOT WORKING! THIS SHOULD BE PULLING id_par AND IT CAN'T BECAUSE THAT IS 0 EVERYWHERE
   ! THIS NEEDS TO BE RESOLVED BEFORE WE CAN USE THE ABM FOR SCIENCING
   par =  _STATE_VAR_S_(data%id_I0) ! _STATE_VAR_(data%id_par) local photosynth. active radiation
   
   no3 = _STATE_VAR_(data%id_nit)        ! local no3
   frp = _STATE_VAR_(data%id_frp)        ! local frp

!------------------------------------------------------
! CODE FROM PIBM BIOLOGY SUBROUTINE

   ! tracking daily integrated PAR in each layer
   _DIAG_VAR_(data%id_IPAR) = _DIAG_VAR_(data%id_IPAR) + par * dtdays   !Unit: W m-2

   ! keep track of what the datetime is
   real_day = _STATE_VAR_S_(data%id_yday)
   integer_day = INT(real_day)
   hour_of_day = NINT((real_day - real(integer_day))*24)

   if (hour_of_day == 0) then ! ML reset NPP and PAR every day; maybe this already happens anyway in GLM-AED?
       _DIAG_VAR_(data%id_NPPc) = zero_  !Reset NPPc_
       _DIAG_VAR_(data%id_IPAR) = zero_  !Reset IPAR
   endif

   ! initialize flux diagnostics at 0
   _DIAG_VAR_(data%id_d_pon) = zero_
   _DIAG_VAR_(data%id_d_nit) = zero_
   _DIAG_VAR_(data%id_d_pop) = zero_
   _DIAG_VAR_(data%id_d_frp) = zero_
   _DIAG_VAR_(data%id_d_poc) = zero_
   _DIAG_VAR_(data%id_d_oxy) = zero_

   !Calculate the number of super-individuals (N_) in this vertical layer and obtain their indices (index_)
   N_PAR = ppid
   N_ = 0

   !Get the indices of particles in this grid (for grazing loss)
   allocate(index_(0), stat=AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating index_***"

   do i = 1, N_PAR 
      if (_PTM_STAT_(i,STAT) == 1) then !Ignore dead super-individuals
         N_ = N_ + 1
         if (N_ == 1) then
             allocate(scratch(1), stat=Allocatestatus)
             IF (AllocateStatus /= 0) STOP "*** Problem in allocating scratch***"
             scratch(1) = _PTM_STAT_(i,STAT)

            !Update index_ and scratch has been deallocated
            call move_alloc(scratch, index_)
         else
             allocate(scratch(size(index_) + 1), stat=Allocatestatus)
             IF (AllocateStatus /= 0) STOP "*** Problem in allocating scratch***"
             scratch(1:size(index_))   = index_
             scratch(size(index_) + 1) = i

            !Update index_ and scratch has been deallocated
            call move_alloc(scratch, index_)
         endif
      endif
   enddo

   ! get volume of layer
   Hz = _STATE_VAR_(data%id_lht)*_STATE_VAR_(data%id_larea)

   !Save number of super-individuals per layer
   _DIAG_VAR_(data%id_count) = N_ 

   !Reset total abundance
   Abun_ = 0d0

   IF (N_ > 0) THEN

      allocate(BN(N_), stat=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating BN***"
      BN(:) = 0d0

      allocate(BC(N_), stat=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating BC***"
      BC(:) = 0d0

      allocate(BP(N_), stat=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating BP***"
      BP(:) = 0d0

      ! calculate the amount of nitrogen and total abundances (cells/m3) in each super-individual
      DO m = 1, N_  

         i = index_(m)

         !The amount of N in the super-individual m                                    ML: looks like pmol to mmol conversion then to m3
         BN(m) = p(i)%ptm_state(data%ip_n) * p(i)%ptm_state(data%ip_num) * 1d-9 / Hz   ! ML BN is only used in zooplankton module

         !The amount of C in the super-individual m
         BC(m) = p(i)%ptm_state(data%ip_c) * p(i)%ptm_state(data%ip_num) * 1d-9 / Hz   ! ML BC is used to calculate fitness

         !The amount of P in the super-individual m                                    ML: looks like pmol to mmol conversion then to m3
         BP(m) = p(i)%ptm_state(data%ip_p) * p(i)%ptm_state(data%ip_num) * 1d-9 / Hz   ! ML BP is not currently used but including here for symmetry in case we add zoops later

         !Count the number of cells in each layer
         Abun_ = Abun_ + p(i)%ptm_state(data%ip_num)

      ENDDO !End of looping through superindividuals
   ENDIF ! if N_ > 0

   ! convert abundance to cells per m-3 for diagnostic
   _DIAG_VAR_(data%id_cells) = Abun_ / Hz !Abundances (cells m-3)

   ! Calculate total phytoplankton nitrogen uptake, mortality, and PP (must after calculation of num(t+dt))
   uptake  = 0d0
   uptake_P  = 0d0
   Pmort   = 0d0
   Pmort_P   = 0d0
   Pmort_C   = 0d0
   oxy_flux = 0d0

   if (N_ > 0) then
      do j = 1, N_
         i = index_(j) 

         ! set particle par, temp, and no3 using environmental conditions from this layer
         p(i)%ptm_state(data%ip_par) = par
         p(i)%ptm_state(data%ip_tem) = WaterTemperature
         p(i)%ptm_state(data%ip_no3) = no3
         p(i)%ptm_state(data%ip_frp) = frp

         ! call the PIBM function that is the primary engine for phytoplankton physiology
         ! this is located in aed_pibm_utils.F90
         call GMK98_Ind_TempSizeLight(p(i)%ptm_state(data%ip_tem), p(i)%ptm_state(data%ip_par), p(i)%ptm_state(data%ip_no3), p(i)%ptm_state(data%ip_frp), p(i)%ptm_state(data%ip_Topt),&
               p(i)%ptm_state(data%ip_c), p(i)%ptm_state(data%ip_n), p(i)%ptm_state(data%ip_p), p(i)%ptm_state(data%ip_chl), p(i)%ptm_state(data%ip_cdiv), exp(p(i)%ptm_state(data%ip_LnalphaChl)),&
               dC_, dN_, dP_, dChl_, ESD_, data%phytos(1)%RC, data%phytos(1)%RN, data%phytos(1)%RP, data%phytos(1)%RChl, data%phytos(1)%zeta_N, data%phytos(1)%zeta_P, data%phytos(1)%a1,&
               data%phytos(1)%R_growth, data%phytos(1)%nx, data%phytos(1)%thetaNmax, data%phytos(1)%QNmin_a, data%phytos(1)%QNmin_b, data%phytos(1)%QNmax_a, data%phytos(1)%QNmax_b, &
               data%phytos(1)%QPmin_a, data%phytos(1)%QPmin_b, data%phytos(1)%QPmax_a, data%phytos(1)%QPmax_b, data%phytos(1)%KN_a, data%phytos(1)%KN_b, data%phytos(1)%KPho_a, data%phytos(1)%KPho_b)

         ! cumulate N and P uptake and oxy flux after running physiology function
         uptake   =   uptake + dN_ * p(i)%ptm_state(data%ip_num) ! Unit: pmol N d-1
         uptake_P =   uptake_P + dP_ * p(i)%ptm_state(data%ip_num) ! Unit: pmol P d-1
         oxy_flux =   oxy_flux + dC_ * p(i)%ptm_state(data%ip_num) ! Unit: pmol C d-1 (assuming 1:1 stoichiometry with O2)
         !NPPc_(k) = NPPc_(k) + dC_ * p_PHY(i)%num *1d-9/Hz(k)*12.d0*dtdays !Unit: mgC m-3 d-1 !ML leaving this PIBM code because of questions below re: NPP diagnostic
         _DIAG_VAR_(data%id_NPPc) = _DIAG_VAR_(data%id_NPPc) + dC_ * p(i)%ptm_state(data%ip_num) * 1d-9 / Hz * 12.d0 * dtdays !Unit: mgC m-3 d-1 ML why multiplied by 12? and also multiplying by dtdays means this is hr-1 (timestep) not d-1

         ! Save carbon-specific growth rate
         p(i)%ptm_state(data%ip_mu_C) = dC_ / p(i)%ptm_state(data%ip_c)

         ! Update cellular C, N, and Chl
         p(i)%ptm_state(data%ip_c)   = p(i)%ptm_state(data%ip_c)   + dC_   * dtdays !Unit: pmol C cell-1 hr-1
         p(i)%ptm_state(data%ip_n)   = p(i)%ptm_state(data%ip_n)   + dN_   * dtdays !Unit: pmol N cell-1 hr-1
         p(i)%ptm_state(data%ip_p)   = p(i)%ptm_state(data%ip_p)   + dP_   * dtdays !Unit: pmol P cell-1 hr-1
         p(i)%ptm_state(data%ip_chl) = p(i)%ptm_state(data%ip_chl) + dChl_ * dtdays !Unit: pg Chl cell-1 hr-1

         ! Update ESD
         p(i)%ptm_env(DIAM) = ESD_*1d-6 !meters

         !# PHYTOPLANKTON CELL DENSITY
         ! density increases during carbohydrate creation (daytime)
         dens_flux = zero_
         IF( par>zero_ ) THEN
            dens_flux = data%phytos(phy_i)%c1 * &
            (one_ - EXP(-par/data%phytos(phy_i)%I_K) ) - data%phytos(phy_i)%c3
         ELSE
            ! darkness
            dens_flux = -data%phytos(phy_i)%c3
         ENDIF

         ! Update density
         p(i)%ptm_env(DENS) = p(i)%ptm_env(DENS) + dens_flux

         ! check maximum/minimum density are not exceeded
         IF( p(i)%ptm_env(DENS)>max_rho ) THEN
            p(i)%ptm_env(DENS)=max_rho
         ENDIF
         IF(p(i)%ptm_env(DENS)<min_rho ) THEN
            p(i)%ptm_env(DENS)=min_rho
         ENDIF

         !# PHYTOPLANKTON SETTLING
         SELECT CASE (data%phytos(1)%settling)

            CASE ( _MOB_OFF_ )
               ! disable settling by setting vertical velocity to 0
               p(i)%ptm_env(VVEL) = zero_

            CASE ( _MOB_CONST_ )
               ! constant settling velocity using user provided value
               p(i)%ptm_env(VVEL) = data%phytos(1)%w_p

            CASE ( _MOB_TEMP_ )
               ! constant settling velocity @20C corrected for density changes
               pw   = _STATE_VAR_(data%id_dens)
               mu   = water_viscosity(WaterTemperature)
               mu20 = 0.001002  ! N s/m2
               pw20 = 998.2000  ! kg/m3 (assuming freshwater)
               p(i)%ptm_env(VVEL) = data%phytos(1)%w_p*mu20*pw / ( mu*pw20 )

            CASE ( _MOB_STOKES_ )
               ! Update p_vvel (Stokes' Law)
               pw   = _STATE_VAR_(data%id_dens)
               mu   = water_viscosity(WaterTemperature)
               p_vvel = -9.807*((ESD_*1d-6)**2.)*( p(i)%ptm_env(DENS)-pw ) / ( 18.*mu ) 
               p(i)%ptm_env(VVEL) = p_vvel

            CASE DEFAULT
               ! unknown settling/migration option selection
               print *, "Settling option not recognized; vertical velocity set to 0"
               p(i)%ptm_env(VVEL) =  zero_

         END SELECT

         !# PHYTOPLANKTON MORTALITY
         ! If celular carbon is lower than the subsistence threshold (Cmin), it dies:
         Cmin = 0.25d0 * p(i)%ptm_state(data%ip_cdiv) ! ML why multiplied by 0.25? this is from PIBM

         if (p(i)%ptm_state(data%ip_c) < Cmin .or. (p(i)%ptm_state(data%ip_num)*p(i)%ptm_state(data%ip_n)) < Nt_min .or. (p(i)%ptm_state(data%ip_num)*p(i)%ptm_state(data%ip_p)) < P_min) then  ! The superindividual Dies
            _DIAG_VAR_(data%id_N_death) = _DIAG_VAR_(data%id_N_death) + 1
            Pmort = Pmort + p(i)%ptm_state(data%ip_n) * p(i)%ptm_state(data%ip_num) !Natural mortality of phytoplankton ==> DET
            Pmort_P = Pmort_P + p(i)%ptm_state(data%ip_p) * p(i)%ptm_state(data%ip_num) !Natural mortality of phytoplankton ==> DET
            Pmort_C = Pmort_C + p(i)%ptm_state(data%ip_c) * p(i)%ptm_state(data%ip_num) !Natural mortality of phytoplankton ==> DET

            p(i)%ptm_state(data%ip_c)   = 0d0
            p(i)%ptm_state(data%ip_n)   = 0d0
            p(i)%ptm_state(data%ip_p)   = 0d0
            p(i)%ptm_state(data%ip_chl) = 0d0
            p(i)%ptm_state(data%ip_num) = 0d0
            _PTM_STAT_(i,STAT) = 0
            _PTM_STAT_(i,FLAG) = 3
         endif

         ! Now probabilistic mortality to compensate because we don't have zoops, disease, or senescence in ABM
         call random_number(rnm)

         IF (rnm < data%phytos(1)%mort_prob) THEN !it dies
            _DIAG_VAR_(data%id_N_death) = _DIAG_VAR_(data%id_N_death) + 1.
            Pmort = Pmort + p(i)%ptm_state(data%ip_n) * p(i)%ptm_state(data%ip_num) !Natural mortality of phytoplankton ==> DET
            Pmort_P = Pmort_P + p(i)%ptm_state(data%ip_p) * p(i)%ptm_state(data%ip_num) !Natural mortality of phytoplankton ==> DET
            Pmort_C = Pmort_C + p(i)%ptm_state(data%ip_c) * p(i)%ptm_state(data%ip_num) !Natural mortality of phytoplankton ==> DET

            p(i)%ptm_state(data%ip_c)   = 0d0
            p(i)%ptm_state(data%ip_n)   = 0d0
            p(i)%ptm_state(data%ip_p)   = 0d0
            p(i)%ptm_state(data%ip_chl) = 0d0
            p(i)%ptm_state(data%ip_num) = 0d0
            _PTM_STAT_(i,STAT) = 0
            _PTM_STAT_(i,FLAG) = 3
         ENDIF

         !Save fitness of each particle (per day) !ML not sure what the purpose of this diagnostic is
         p(i)%ptm_state(data%ip_fit) = (p(i)%ptm_state(data%ip_c)*p(i)%ptm_state(data%ip_num)*1d-9/Hz - BC(j))/BC(j)/dtdays

      enddo !End of looping through particles
  endif !End if of N_ > 0

  ! Convert all updates to be per m3 per second
  uptake = uptake*1d-9/Hz/secs_per_day !Convert uptake to mmol N m-3 s-1
  uptake_P = uptake_P*1d-9/Hz/secs_per_day !Convert uptake_P to mmol P m-3 s-1
  oxy_flux = oxy_flux*1d-9/Hz/secs_per_day !Convert oxy flux to mmol O2 m-3 s-1
  Pmort  =  Pmort*1d-9/Hz/secs_per_day !Convert Pmort to mmol N m-3 s-1
  Pmort_P  =  Pmort_P*1d-9/Hz/secs_per_day !Convert Pmort_P to mmol P m-3 s-1
  Pmort_C  =  Pmort_C*1d-9/Hz/secs_per_day !Convert Pmort_C to mmol C m-3 s-1

  ! Now increment water column properties based on particle fluxes
  _FLUX_VAR_(data%id_nit)   = _FLUX_VAR_(data%id_nit) + (pp_ND + RES - uptake)    ! mmol/m3/s CHECK ML I think this should not be multiplied by dtdays so stays in seconds; pp_ND and RES are zooplankton variables
  _DIAG_VAR_(data%id_d_nit) = (pp_ND + RES - uptake)  * secs_per_day    ! mmol/m3/day
  
  _FLUX_VAR_(data%id_frp)   = _FLUX_VAR_(data%id_frp) + (pp_PD + RES_P - uptake_P)    ! mmol/m3/s CHECK ML I think this should not be multiplied by dtdays so stays in seconds; pp_PD and RES_P are zooplankton variables
  _DIAG_VAR_(data%id_d_frp) = (pp_PD + RES_P - uptake_P)  * secs_per_day    ! mmol/m3/day

  _FLUX_VAR_(data%id_pon) = _FLUX_VAR_(data%id_pon)  + Pmort + (pp_DZ - pp_ND) ! mmol/m3/s CHECK pp_DZ and pp_ND are zooplankton variables
  _DIAG_VAR_(data%id_d_pon) = Pmort + (pp_DZ - pp_ND) * secs_per_day ! mmol/m3/day
  
  _FLUX_VAR_(data%id_pop) = _FLUX_VAR_(data%id_pop)  + Pmort_P + (pp_DZP - pp_PD) ! mmol/m3/s CHECK pp_DZP and pp_PD are zooplankton variables
  _DIAG_VAR_(data%id_d_pop) = Pmort_P + (pp_DZP - pp_PD) * secs_per_day ! mmol/m3/day
  
  _FLUX_VAR_(data%id_poc) = _FLUX_VAR_(data%id_poc)  + Pmort_C  ! mmol/m3/s CHECK would need to add equivalent zoop variables here if added zoops in
  _DIAG_VAR_(data%id_d_poc) = Pmort_C * secs_per_day ! mmol/m3/day would need to add equivalent zoop variables here if added zoops in

  _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + oxy_flux  ! mmol/m3/s CHECK would need to add equivalent zoop variables here if added zoops in
  _DIAG_VAR_(data%id_d_oxy) = oxy_flux * secs_per_day ! mmol/m3/day would need to add equivalent zoop variables here if added zoops in

  ! deallocate things
  if (allocated(index_))  deallocate(index_)
  if (allocated(Pmatrix)) deallocate(Pmatrix)
  if (allocated(BN))      deallocate(BN)
  if (allocated(BP))      deallocate(BP)
  if (allocated(BC))      deallocate(BC)

! END OF CODE FROM PIBM BIOLOGY SUBROUTINE
!------------------------------------------------------


!------------------------------------------------------
! CODE FROM PIBM Par2PHY SUBROUTINE

! initialize diag vars at 0
_DIAG_VAR_(data%id_phyc) = 0.
_DIAG_VAR_(data%id_phyn) = 0.
_DIAG_VAR_(data%id_phyp) = 0.
_DIAG_VAR_(data%id_chl) = 0.
_DIAG_VAR_(data%id_mcdiv) = 0.
_DIAG_VAR_(data%id_mtopt) = 0.
_DIAG_VAR_(data%id_mlnalpha) = 0.

!loop through all super-individuals to divide, mutate, and convert into Eulerian concentrations
DO i = 1, N_  

   ! Cell division: if cellular carbon is above the division threshold, it divides
   IF (p(i)%ptm_state(data%ip_c) >= p(i)%ptm_state(data%ip_cdiv)) THEN  !Divide
      !print *, 'I divided!'
      _DIAG_VAR_(data%id_N_birth) = _DIAG_VAR_(data%id_N_birth) + 1.
      p(i)%ptm_state(data%ip_c)   = p(i)%ptm_state(data%ip_c)/2d0
      p(i)%ptm_state(data%ip_n)   = p(i)%ptm_state(data%ip_n)/2d0
      p(i)%ptm_state(data%ip_p)   = p(i)%ptm_state(data%ip_p)/2d0
      p(i)%ptm_state(data%ip_chl) = p(i)%ptm_state(data%ip_chl)/2d0
      p(i)%ptm_state(data%ip_num) = p(i)%ptm_state(data%ip_num)*2d0

      !Mutation
      If (NTrait > 0) Then
         DO m = 1, NTrait
            nu_ = p(i)%ptm_state(data%ip_num)*nu(m)
            call random_number(cff)

            IF (cff < nu_) THEN !Mutation occurs
               _DIAG_VAR_(data%id_N_mutate) = _DIAG_VAR_(data%id_N_mutate) + 1.

               select case(m)
               case(iTopt)
                   oldtt(1) = p(i)%ptm_state(data%ip_Topt)
               case(iSize)
                   oldtt(1) = log(p(i)%ptm_state(data%ip_Cdiv)) !ML I could understand CDiv being logged to avoid too-big changes in division threshold but why logged below?
               case(ialphaChl)
                   oldtt(1) = p(i)%ptm_state(data%ip_LnalphaChl)
               case DEFAULT
                   stop "Trait index wrong!"
               end select

               vartt(1,1)= sigma(m)**2   !Construct the covariance matrix for the selected trait

               !A new Topt is randomly sampled from a Gaussian distribution with mean of previous Topt and SD of sigma
               newtt = srand_mtGaus(1, oldtt, vartt)
               select case(m)
               case(iTopt)
                   p(i)%ptm_state(data%ip_Topt) = newtt(1)
               case(iSize)
                   p(i)%ptm_state(data%ip_Cdiv) = exp(newtt(1))
               case(ialphaChl)
                   p(i)%ptm_state(data%ip_LnalphaChl) = newtt(1)
               case DEFAULT
                   stop "Trait index wrong!"
               end select
            ENDIF
         ENDDO !End of looping through Traits
      ENDIF !End of if (NTrait > 0)
   ENDIF !End of division

   !Cumulate concentrations of phyto C, N, and Chl, and mean traits for superindividuals in the layer
   PHYC = PHYC + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c)
   PHYN = PHYN + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_n)
   PHYP = PHYP + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_p)
   CHL  = CHL  + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_chl)
   _DIAG_VAR_(data%id_mcdiv)    = _DIAG_VAR_(data%id_mcdiv)    + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * log(p(i)%ptm_state(data%ip_cdiv))
   _DIAG_VAR_(data%id_mtopt)    = _DIAG_VAR_(data%id_mtopt)    + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * p(i)%ptm_state(data%ip_Topt)
   _DIAG_VAR_(data%id_mlnalpha) = _DIAG_VAR_(data%id_mlnalpha) + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * p(i)%ptm_state(data%ip_LnalphaChl)
ENDDO !End of iterating over all super-individuals

!Calculate Eulerian concentrations of phyto C, N, and Chl, and mean traits and covariances for each layer
_DIAG_VAR_(data%id_phyc) = PHYC * 1d-9/Hz   !Convert Unit to mmol/m^3
_DIAG_VAR_(data%id_phyn) = PHYN * 1d-9/Hz   !Convert Unit to mmol/m^3
_DIAG_VAR_(data%id_phyp) = PHYP * 1d-9/Hz   !Convert Unit to mmol/m^3
_DIAG_VAR_(data%id_chl)  = CHL  * 1d-9/Hz   !Convert Unit to mg Chl/m^3
if (PHYC > 0.) then   ! need to make sure we are preventing dividing by 0 (case when there are no particles in layer)
   _DIAG_VAR_(data%id_mcdiv)    = _DIAG_VAR_(data%id_mcdiv)    / PHYC
   _DIAG_VAR_(data%id_mtopt)    = _DIAG_VAR_(data%id_mtopt)    / PHYC
   _DIAG_VAR_(data%id_mlnalpha) = _DIAG_VAR_(data%id_mlnalpha) / PHYC
endif

!Compute trait covariances
_DIAG_VAR_(data%id_vcdiv)    = 0d0
_DIAG_VAR_(data%id_vtopt)    = 0d0
_DIAG_VAR_(data%id_vlnalpha) = 0d0
_DIAG_VAR_(data%id_cov_AL)   = 0d0
_DIAG_VAR_(data%id_cov_TA)   = 0d0
_DIAG_VAR_(data%id_cov_TL)   = 0d0

DO i = 1, N_ !ML removing logs from ip_cdiv here too - pretty sure this is a typo given how it is handled in GMK routine OR should be logged everywhere
   _DIAG_VAR_(data%id_vcdiv)    = _DIAG_VAR_(data%id_vcdiv)    + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * (log(p(i)%ptm_state(data%ip_cdiv)) - _DIAG_VAR_(data%id_mcdiv))**2
   _DIAG_VAR_(data%id_vtopt)    = _DIAG_VAR_(data%id_vtopt)    + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * (p(i)%ptm_state(data%ip_Topt) - _DIAG_VAR_(data%id_mtopt))**2
   _DIAG_VAR_(data%id_vlnalpha) = _DIAG_VAR_(data%id_vlnalpha) + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * (p(i)%ptm_state(data%ip_LnalphaChl) - _DIAG_VAR_(data%id_mlnalpha))**2
   _DIAG_VAR_(data%id_cov_TL)   = _DIAG_VAR_(data%id_cov_TL)   + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * (log(p(i)%ptm_state(data%ip_cdiv)) - _DIAG_VAR_(data%id_mcdiv)) * (p(i)%ptm_state(data%ip_Topt) - _DIAG_VAR_(data%id_mtopt))
   _DIAG_VAR_(data%id_cov_AL)   = _DIAG_VAR_(data%id_cov_AL)   + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * (log(p(i)%ptm_state(data%ip_cdiv)) - _DIAG_VAR_(data%id_mcdiv)) * (p(i)%ptm_state(data%ip_LnalphaChl) - _DIAG_VAR_(data%id_mlnalpha))
   _DIAG_VAR_(data%id_cov_TA)   = _DIAG_VAR_(data%id_cov_TA)   + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * (p(i)%ptm_state(data%ip_LnalphaChl) - _DIAG_VAR_(data%id_mlnalpha)) * (p(i)%ptm_state(data%ip_Topt) - _DIAG_VAR_(data%id_mtopt))
ENDDO

if (PHYC > 0d0) then ! need to make sure we are preventing dividing by 0 (case when there are no particles in layer)
   _DIAG_VAR_(data%id_vcdiv)    = _DIAG_VAR_(data%id_vcdiv)   / PHYC
   _DIAG_VAR_(data%id_vtopt)    = _DIAG_VAR_(data%id_vtopt)   / PHYC
   _DIAG_VAR_(data%id_vlnalpha) = _DIAG_VAR_(data%id_vlnalpha)/ PHYC
   _DIAG_VAR_(data%id_cov_TL)   = _DIAG_VAR_(data%id_cov_TL)  / PHYC
   _DIAG_VAR_(data%id_cov_TA)   = _DIAG_VAR_(data%id_cov_TA)  / PHYC
   _DIAG_VAR_(data%id_cov_AL)   = _DIAG_VAR_(data%id_cov_AL)  / PHYC
  endif

! END PIBM Par2PHY CODE
!------------------------------------------------------

END SUBROUTINE aed_particle_bgc_phyto_abm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_extinction_phyto_abm(data,column,layer_idx,extinction)
!------------------------------------------------------------------------------+
! Get the light extinction coefficient due to biogeochemical variables
!------------------------------------------------------------------------------+
!ARGUMENTS
   CLASS (aed_phyto_abm_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: phy
   INTEGER  :: phy_i
!
!------------------------------------------------------------------------------+
!BEGIN

   DO phy_i=1,data%num_phytos
      ! Retrieve current (local) phytoplankton biomass in this layer
      phy = _DIAG_VAR_(data%id_phyc) !MH will need expanding to multi-groups for phy>1
      
      ! Self-shading with contribution from this phytoplankton concentration.
      extinction = extinction + (data%phytos(phy_i)%KePHY*phy)
   ENDDO
END SUBROUTINE aed_light_extinction_phyto_abm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_phyto_abm
