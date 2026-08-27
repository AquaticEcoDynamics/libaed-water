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
!#  Copyright 2018-2026 : The University of Western Australia                  #
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
      ! DIAGNOSTIC: one entry per simulated phyto group (data%phytos(g)), sized/allocated
      ! and registered (name suffixed by that group's p_name) in aed_define_phyto_abm.
      INTEGER,DIMENSION(:),ALLOCATABLE :: id_count
      INTEGER,DIMENSION(:),ALLOCATABLE :: id_phyc, id_phyn, id_phyp, id_chl, id_mtopt, id_vtopt, &
                 id_mcdiv, id_vcdiv, id_mlnalpha, id_vlnalpha, &
                 id_cov_ta, id_cov_tl, id_cov_al, id_N_birth,  &
                 id_N_mutate, id_NPPc, id_cells,      &
                 id_N_death
      ! id_IPAR stays a single scalar: it is accumulated daily integrated layer PAR
      ! (an environmental forcing), not a per-species phyto response quantity.
      INTEGER :: id_IPAR
      ! DIAGNOSTIC: cumulative (never rezero'd) biological loss/gain counters, since sim
      ! start, so total births/deaths can be attributed by cause without the ~24x
      ! undercount that comes from summing the hourly-reset id_N_birth/id_N_death
      ! diagnostics at daily (or coarser) output frequency.
      INTEGER,DIMENSION(:),ALLOCATABLE :: id_cum_birth, id_cum_death_thresh, id_cum_death_prob
      INTEGER,DIMENSION(:),ALLOCATABLE :: id_cum_cells_decayed
      INTEGER :: ip_c, ip_n, ip_p, ip_par, ip_tem, ip_no3, ip_nh4, ip_frp, ip_chl, ip_num, ip_cdiv, ip_Topt, ip_LnalphaChl
      ! Per-particle carryover of rhoChl (fraction of C production devoted to Chl synthesis)
      ! from the end of the last daylight period, used by GMK98_Ind_TempSizeLight during dark
      ! hours. Registered LAST (after ip_LnalphaChl) so it doesn't shift the fixed VVEL+N offsets
      ! GLM/src/glm_ptm.c uses to read the other particle state variables out of ptm_env.
      INTEGER :: ip_rhochl
      INTEGER,DIMENSION(:),ALLOCATABLE :: id_d_oxy, id_d_nit, id_d_amm, id_d_pon, id_d_frp, id_d_pop, id_d_poc
      INTEGER :: id_d_dc, id_d_dn, id_d_dp
      INTEGER :: id_oxy,id_amm,id_nit,id_frp,id_doc,id_don,id_dop,id_poc,id_pon,id_pop
      INTEGER :: id_lht, id_larea, id_dep, id_tem, id_par, id_I0, id_dens, id_yday, id_depth
      INTEGER :: id_lvol

      ! Phytoplankton parameters
      INTEGER  :: num_phytos
      TYPE(phyto_data_t),DIMENSION(:),ALLOCATABLE :: phytos

      ! FEATURE: minimum active-particle population maintenance. When > 0, and the
      ! domain-wide active particle count (across ALL groups) falls below this value,
      ! the largest active particle IN THE WHOLE DOMAIN (any group) is split (see
      ! aed_split_particle_phyto_abm) to replace losses. This is intentionally a
      ! domain-wide computational/statistical-resolution safeguard, not a per-species
      ! ecological control, so it stays a single scalar (unlike recruit_rate below).
      ! Set via 'target_particle_num' in the &aed_phyto_abm namelist (api.nml).
      INTEGER :: target_particle_num = 0

      ! FEATURE: general biological recruitment. When recruit_rate(g) > 0, new active
      ! particles of group g are introduced at a steady rate (particles/day), independent
      ! of both initial seeding (glm_ptm.c ptm_init_glm) and inflow-sourced particles
      ! (inflow_conc in &particles). Represents an unmodeled external source (e.g.
      ! resuspension/immigration). Each recruit is placed at an independently-random
      ! height in [recruit_height_min(g), recruit_height_max(g)] (absolute height above
      ! the lake-bottom datum, matching ptm_env(HGHT)'s convention - NOT depth-below-
      ! surface like init_depth_min/max) and is assigned trait values drawn the same
      ! way initial-population diversity is generated in aed_particle_initialize_
      ! phyto_abm - a Gaussian perturbation (srand_mtGaus, sigma_init from module
      ! params) around group g's aed_phyto_pars.csv base defaults, rather than around a
      ! stale ptm_state left over from a previous particle life.
      ! Per-group (one value per simulated phyto group, unlike target_particle_num
      ! above): set via 'recruit_rate' (particles/day; 0 = disabled), 'recruit_height_min',
      ! 'recruit_height_max' (m), and 'recruit_max_per_call' in the &aed_phyto_abm
      ! namelist (api.nml). See aed_split_particle_phyto_abm, block (3).
      AED_REAL,DIMENSION(:),ALLOCATABLE :: recruit_rate
      AED_REAL,DIMENSION(:),ALLOCATABLE :: recruit_height_min
      AED_REAL,DIMENSION(:),ALLOCATABLE :: recruit_height_max
      INTEGER,DIMENSION(:),ALLOCATABLE  :: recruit_max_per_call

      ! Bounds on per-particle cell density (kg/m3), applied to ptm_env(DENS) after the
      ! carbohydrate-ballasting update in aed_particle_bgc_phyto_abm. Set via 'min_rho' /
      ! 'max_rho' in the &aed_phyto_abm namelist (api.nml).
      ! Stored here because the namelist variables of the same name were previously read
      ! and then never used: the clamp read two SEPARATE locals declared inside
      ! aed_particle_bgc_phyto_abm and hardcoded to the same 900./1200. defaults, so
      ! setting min_rho/max_rho in the namelist had no effect whatsoever.
      ! APPENDED AT THE END of the component list deliberately - inserting these
      ! mid-structure (directly after 'phytos') produced a segfault in memmove during
      ! model configuration, before this routine's own code was even reached.
      AED_REAL :: min_rho = 900.
      AED_REAL :: max_rho = 1200.

      ! Per-group reseeding floor/target (2026-08): threshold-triggered population
      ! maintenance, complementary to recruit_rate above rather than a replacement for
      ! it. recruit_rate models an ongoing BIOLOGICAL source (continuous rate, its own
      ! height band, participates fully in growth/mortality/biogeochemistry from the
      ! moment it appears). Reseeding is INSTRUMENT MAINTENANCE: if a group's live
      ! particle count falls below particle_reseed_min(g), it is topped back up to
      ! particle_reseed_target(g) in one call. Physiology for the new particles is
      ! seeded via the SAME recruit_particle() Gaussian-draw-around-CSV-defaults
      ! routine recruit_rate already uses - not duplicated CSV literals - since this
      ! model's growth kernel cannot bootstrap a particle from a zeroed C/N/P/Chl/num
      ! state (dC/dt is proportional to C; C=0 is a fixed point, not a starting point).
      ! 0 = disabled (default) for both.
      INTEGER,DIMENSION(:),ALLOCATABLE :: particle_reseed_min
      INTEGER,DIMENSION(:),ALLOCATABLE :: particle_reseed_target

      CONTAINS
         PROCEDURE :: define              => aed_define_phyto_abm
!        PROCEDURE :: calculate           => aed_calculate_phyto_abm
!        PROCEDURE :: calculate_benthic   => aed_calculate_benthic_phyto_abm
!        PROCEDURE :: calculate_riparian  => aed_calculate_riparian_phyto_abm
!        PROCEDURE :: calculate_dry       => aed_calculate_dry_phyto_abm
!        PROCEDURE :: equilibrate         => aed_equilibrate_phyto_abm
         PROCEDURE :: initialize_particle => aed_particle_initialize_phyto_abm
         PROCEDURE :: split_particle      => aed_split_particle_phyto_abm
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
   ! Which data%phytos(:) entry (phyto group/species) this particle belongs to. Set
   ! explicitly whenever a particle is created/reassigned (initial seeding, split,
   ! recruitment) - never inferred from the particle's position in p(:), since free
   ! slots get reused across groups.
   ! IMPORTANT: PTM_STATUS above (=19) is an unrelated, vestigial "TBD" constant and
   ! does NOT size this array - the real ptm_istat array is allocated with exactly
   ! n_ptm_istat elements in aed_ptm.F90 (libaed-api). GRP must not exceed that value;
   ! n_ptm_istat was bumped from 6 to 7 there specifically to make room for GRP.
   INTEGER, PARAMETER :: GRP  = 7

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

   ! DIAGNOSTIC: cumulative count of successful particle splits (the only mechanism in this
   ! module that creates a new active particle / increases the particle count - "births" via
   ! cell division elsewhere in this file only double cell count within an existing particle,
   ! they do not create a new one). aed_split_particle_phyto_abm has no column/layer_idx
   ! context, so this is tracked as a plain module-level counter and logged directly to file
   ! rather than through the usual per-layer _DIAG_VAR_ mechanism.
   INTEGER, SAVE :: ptm_cum_split_count = 0

   ! DIAGNOSTIC: cumulative count of splits performed to maintain the minimum active
   ! particle population (target_particle_num), as distinct from n_split-triggered splits
   ! above. Logged to its own file for the same reason as ptm_cum_split_count.
   INTEGER, SAVE :: ptm_cum_popmaintain_split_count = 0

   ! DIAGNOSTIC: cumulative count of recruitment events (the third, independent
   ! pathway that creates a new active particle - see also (1) and (2) above).
   INTEGER, SAVE  :: ptm_cum_recruit_count = 0

   ! DIAGNOSTIC: cumulative count of reseed events (the fourth, independent pathway -
   ! see block (4) in aed_split_particle_phyto_abm, and particle_reseed_min/target on
   ! the data type).
   INTEGER, SAVE  :: ptm_cum_reseed_count = 0
   ! Fractional-particle accumulator for recruit_rate (particles/day), one per phyto
   ! group since recruit_rate is now per-group. Allocated/sized in aed_define_phyto_abm
   ! once num_phytos is known. See block (3) in aed_split_particle_phyto_abm for why a
   ! fractional accumulator is used instead of a fixed integer-per-call count.
   AED_REAL, SAVE, DIMENSION(:), ALLOCATABLE :: recruit_accumulator


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
            CASE ('thetaPmax')     ; pd(dcol)%thetaPmax     = extract_double(values(ccol))
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
            CASE ('a_c2vol')       ; pd(dcol)%a_c2vol       = extract_double(values(ccol))
            CASE ('b_c2vol')       ; pd(dcol)%b_c2vol       = extract_double(values(ccol))
            CASE ('a_pmax')        ; pd(dcol)%a_pmax        = extract_double(values(ccol))
            CASE ('rho')           ; pd(dcol)%rho           = extract_double(values(ccol))
            CASE ('rho_star')      ; pd(dcol)%rho_star      = extract_double(values(ccol))
            CASE ('b_rho')         ; pd(dcol)%b_rho         = extract_double(values(ccol))
            CASE ('V_s')           ; pd(dcol)%V_s           = extract_double(values(ccol))
            CASE ('Ea0')           ; pd(dcol)%Ea0           = extract_double(values(ccol))
            CASE ('Ed0')           ; pd(dcol)%Ed0           = extract_double(values(ccol))
            CASE ('Ei')            ; pd(dcol)%Ei            = extract_double(values(ccol))
            CASE ('beta')          ; pd(dcol)%beta          = extract_double(values(ccol))
            CASE ('phi')           ; pd(dcol)%phi           = extract_double(values(ccol))
            CASE ('Tau')           ; pd(dcol)%Tau           = extract_double(values(ccol))
            CASE ('Beta_Ainf')     ; pd(dcol)%Beta_Ainf     = extract_double(values(ccol))
            CASE ('Kappa')         ; pd(dcol)%Kappa         = extract_double(values(ccol))
            CASE ('Kd_Ainf')       ; pd(dcol)%Kd_Ainf       = extract_double(values(ccol))
            CASE ('a_')            ; pd(dcol)%a_            = extract_double(values(ccol))
            CASE ('b_')            ; pd(dcol)%b_            = extract_double(values(ccol))
            CASE ('v_')            ; pd(dcol)%v_            = extract_double(values(ccol))
            CASE ('simSplit')      ; pd(dcol)%simSplit      = extract_integer(values(ccol))
            CASE ('n_split')       ; pd(dcol)%n_split       = extract_double(values(ccol))
            CASE ('simMortNumDecay'); pd(dcol)%simMortNumDecay = extract_integer(values(ccol))


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
       data%phytos(i)%thetaPmax     = pd(list(i))%thetaPmax
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
       data%phytos(i)%a_c2vol       = pd(list(i))%a_c2vol
       data%phytos(i)%b_c2vol       = pd(list(i))%b_c2vol
       data%phytos(i)%a_pmax        = pd(list(i))%a_pmax
       data%phytos(i)%rho           = pd(list(i))%rho
       data%phytos(i)%rho_star      = pd(list(i))%rho_star
       data%phytos(i)%b_rho         = pd(list(i))%b_rho
       data%phytos(i)%V_s           = pd(list(i))%V_s
       data%phytos(i)%Ea0           = pd(list(i))%Ea0
       data%phytos(i)%Ed0           = pd(list(i))%Ed0
       data%phytos(i)%Ei            = pd(list(i))%Ei
       data%phytos(i)%beta          = pd(list(i))%beta
       data%phytos(i)%phi           = pd(list(i))%phi
       data%phytos(i)%Tau           = pd(list(i))%Tau
       data%phytos(i)%Beta_Ainf     = pd(list(i))%Beta_Ainf
       data%phytos(i)%Kappa         = pd(list(i))%Kappa
       data%phytos(i)%Kd_Ainf       = pd(list(i))%Kd_Ainf
       data%phytos(i)%a_            = pd(list(i))%a_
       data%phytos(i)%b_            = pd(list(i))%b_
       data%phytos(i)%v_            = pd(list(i))%v_
       data%phytos(i)%simSplit      = pd(list(i))%simSplit
       data%phytos(i)%n_split       = pd(list(i))%n_split
       data%phytos(i)%simMortNumDecay = pd(list(i))%simMortNumDecay

       ! n_split semantics (2026-08, adaptive splitting - see aed_split_particle_phyto_abm):
       !   > 0  fixed threshold - particle splits once it exceeds this many cells
       !        (the only behaviour this module had before; unchanged, and still the
       !        default for any group whose CSV n_split is a normal positive value)
       !   = 0  never split. n_split has no declared default in phyto_param_t, so an
       !        omitted CSV column would otherwise leave this at whatever value the
       !        compiler happens to give an uninitialised AED_REAL - map it to a very
       !        large number instead, so "not set" reliably means "don't split".
       !   < 0  ADAPTIVE (constant particle NUMBER rather than constant particle
       !        WEIGHT): magnitude is not used, only the sign is a sentinel.
       IF (data%phytos(i)%n_split == 0.d0) data%phytos(i)%n_split = 1000000.d0

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
   INTEGER  :: grp_i   ! loop index over simulated phyto groups, for per-group diagnostic registration

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
   INTEGER            :: target_particle_num = 0   ! 0 = disabled; domain-wide, NOT per-group; see FEATURE note on the type definition
   ! Per-group (one value per simulated phyto group; index by CSV column order via the_phytos)
   AED_REAL           :: recruit_rate(MAX_PHYTO_TYPES)         = 0.0d0   ! particles/day; 0 = disabled
   AED_REAL           :: recruit_height_min(MAX_PHYTO_TYPES)   = 0.0d0   ! m, absolute height above lake-bottom datum
   AED_REAL           :: recruit_height_max(MAX_PHYTO_TYPES)   = 0.0d0   ! m, absolute height above lake-bottom datum
   INTEGER            :: recruit_max_per_call(MAX_PHYTO_TYPES) = 1000   ! safety cap, mirrors MAX_POPMAINTAIN_SPLITS_PER_CALL
   ! Per-group reseeding: 0 = disabled (default) for both. target <= 0 (with min > 0)
   ! falls back to target = min, i.e. "top up to exactly the floor". See FEATURE note
   ! on the type definition.
   INTEGER            :: particle_reseed_min(MAX_PHYTO_TYPES)    = 0
   INTEGER            :: particle_reseed_target(MAX_PHYTO_TYPES) = 0

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
                    phyto_particle_link, R_mpbb, target_particle_num,          &
                    recruit_rate, recruit_height_min, recruit_height_max,      &
                    recruit_max_per_call, particle_reseed_min,                 &
                    particle_reseed_target
!-------------------------------------------------------------------------------
!BEGIN

   ! Initialise

   ! Read the namelist
   read(namlst,nml=aed_phyto_abm,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_phyto_abm'

   data%target_particle_num = target_particle_num

   ! Particle density bounds. Defaults (900./1200.) match the values the clamp previously
   ! hardcoded, so an unchanged namelist reproduces earlier results exactly.
   data%min_rho = min_rho
   data%max_rho = max_rho
   IF ( data%max_rho <= data%min_rho ) THEN
      print *,'aed_phyto_abm: ERROR max_rho (',data%max_rho,') must exceed min_rho (', &
              data%min_rho,')'
      STOP 'Invalid min_rho/max_rho in namelist aed_phyto_abm'
   ENDIF

   ! Store group parameter values in the phyto derived type. data%num_phytos may end up
   ! less than the requested num_phytos if the_phytos(:) contains bad/out-of-range CSV
   ! column indices (aed_phytoplankton_load_params EXITs early in that case) - so the
   ! per-group arrays below are sized from data%num_phytos (the actual count), not the
   ! raw namelist num_phytos.
   CALL aed_phytoplankton_load_params(data,dbase,                                &
                                      num_phytos,the_phytos,settling,resuspension)

   ! Per-group recruitment parameters and their fractional-particle accumulator
   ! (target_particle_num stays a single scalar - see FEATURE note on the type definition)
   ALLOCATE(data%recruit_rate(data%num_phytos))
   ALLOCATE(data%recruit_height_min(data%num_phytos))
   ALLOCATE(data%recruit_height_max(data%num_phytos))
   ALLOCATE(data%recruit_max_per_call(data%num_phytos))
   ALLOCATE(recruit_accumulator(data%num_phytos))
   data%recruit_rate(1:data%num_phytos)         = recruit_rate(1:data%num_phytos)
   data%recruit_height_min(1:data%num_phytos)   = recruit_height_min(1:data%num_phytos)
   data%recruit_height_max(1:data%num_phytos)   = recruit_height_max(1:data%num_phytos)
   data%recruit_max_per_call(1:data%num_phytos) = recruit_max_per_call(1:data%num_phytos)
   recruit_accumulator(1:data%num_phytos)       = 0.0d0

   ! Per-group reseeding floor/target (block (4) in aed_split_particle_phyto_abm).
   ! target <= 0 (with min > 0) falls back to target = min - "top up to the floor",
   ! the simplest sensible default when no explicit target is given.
   ALLOCATE(data%particle_reseed_min(data%num_phytos))
   ALLOCATE(data%particle_reseed_target(data%num_phytos))
   data%particle_reseed_min(1:data%num_phytos)    = particle_reseed_min(1:data%num_phytos)
   data%particle_reseed_target(1:data%num_phytos) = particle_reseed_target(1:data%num_phytos)
   DO grp_i = 1, data%num_phytos
      IF (data%particle_reseed_min(grp_i) > 0 .AND. data%particle_reseed_target(grp_i) <= 0) &
         data%particle_reseed_target(grp_i) = data%particle_reseed_min(grp_i)
      IF (data%particle_reseed_min(grp_i) > 0 .AND. data%recruit_rate(grp_i) > 0) &
         print *,'        WARNING aed_phyto_abm: group ',grp_i,' has both recruit_rate > 0 ', &
                 'and particle_reseed_min > 0 - both mechanisms will add particles.'
   ENDDO

   print *,"        aed_phyto_abm initialization"

   ! Register particle state variables (phytoplankton agents). These are single generic
   ! slots shared by every particle regardless of species (unlike the Eulerian
   ! aed_phytoplankton module, each particle already carries its own individual state,
   ! so there's no need for per-group arrays of these) - hence a fixed generic name
   ! rather than one built from any particular group's p_name. The initial= defaults
   ! below (from group 1) are only a transient host-side seed value: every particle's
   ! state is overwritten to match its own group's parameters in
   ! aed_particle_initialize_phyto_abm before any biology runs.

   !real    :: PAR = 0.01  associated PAR                                       _PTM_VAR_          p(..)%ptm_state(data%ip_par)
   data%ip_par = aed_define_ptm_variable('phy_par', 'ummol m2 sec', 'particle layer PAR')

   !real    :: temp= 20.   associated Temperature                               _PTM_VAR_
   data%ip_tem = aed_define_ptm_variable('phy_tem', 'degrees C', 'particle layer temperature')

   !real    :: NO3 = 0.1   associated nutrient                                  _PTM_VAR_
   data%ip_no3 = aed_define_ptm_variable('phy_no3', 'mmol N', 'particle layer NO3')

   !real    :: NH4         associated nutrient                                  _PTM_VAR_
   data%ip_nh4 = aed_define_ptm_variable('phy_nh4', 'mmol N', 'particle layer NH4')

   !real    :: FRP = 0.1   associated nutrient                                  _PTM_VAR_
   data%ip_frp = aed_define_ptm_variable('phy_frp', 'mmol P', 'particle layer FRP')

   !real    :: C   = 0.02          ! cellular carbon content (pmol; assuming a 1 micron cell)
   data%ip_c = aed_define_ptm_variable('phy_c', 'pmol C/cell', 'cell C concentration',initial=data%phytos(1)%X_cinit) ! 0.02

   !real    :: N   = 0.02/106.*16. ! cellular nitrogen content (pmol per cell; assuming a 1 micron cell)
   data%ip_n = aed_define_ptm_variable('phy_n', 'pmol N/cell', 'cell N concentration',initial=data%phytos(1)%X_ninit) !0.02/106.*16.

   !real    :: P   = 0.02/106.*1. ! cellular phosphorus content (pmol per cell; assuming a 1 micron cell)
   data%ip_p = aed_define_ptm_variable('phy_p', 'pmol P/cell', 'cell P concentration',initial=data%phytos(1)%X_pinit) !0.02/106.*1.

   !real    :: Chl = 0.02 * 12/50  ! Cellular Chl content (pg Chl)
   data%ip_chl = aed_define_ptm_variable('phy_chl', 'pg Chl', 'cell Chl concentration',initial = data%phytos(1)%X_chlinit) !0.02 * 12/50

   !real    :: num = 5d9           ! Number of cells per superindividual
   data%ip_num = aed_define_ptm_variable('phy_num', 'number', 'number of cells/particle',initial = data%phytos(1)%n0) !5d12

   !real :: Cdiv= 0.04d0           !cellular carbon content threshold for division (pmol/cell), can be used as a proxy for size and can be converted to ESD; Phytoplankton half-saturation constant, minimal N:C and maximal N:C ratios are allometric functions of this parameter
                                   !This trait will vary with mutation
   data%ip_Cdiv = aed_define_ptm_variable('phy_cdiv', 'pmol/cell', 'cellular carbon content threshold for division',initial = data%phytos(1)%Cdiv) ! 0.04d0

   !real :: Topt = 20.d0           !Optimal temperature
   data%ip_Topt = aed_define_ptm_variable('phy_topt', 'degrees C', 'optimal temperature',initial = data%phytos(1)%T_opt)

   !real :: LnalphaChl = -2.3  !log(0.1) !Ln alphaChl (slope of the P-I curve; unit: (W m-2)-1 (gChl molC)-1 d-1 instead of micro mol quanta m-2 s-1)
   data%ip_LnalphaChl = aed_define_ptm_variable('phy_lnalphachl', '(W m-2)-1 (gChl molC)-1 d-1', 'slope of the P-I curve',initial = data%phytos(1)%Lnalphachl) !-2.3

   ! Per-particle carryover of rhoChl from the end of the last daylight period (used by
   ! GMK98_Ind_TempSizeLight during dark hours). Must be per-particle state, NOT a variable
   ! local to that subroutine - see the long comment in aed_pibm_utils.F90 for why. Registered
   ! last so it doesn't shift the fixed VVEL+N offsets glm_ptm.c uses for the other variables.
   data%ip_rhochl = aed_define_ptm_variable('phy_rhochl', 'gChl/molC', 'carbon production allocated to Chl synthesis, carried over from the last daylight period', initial = 0.0d0)


   ! Diagnostic outputs for particle properties. One entry per simulated phyto group,
   ! with the variable name suffixed by that group's p_name (mirroring the Eulerian
   ! aed_phytoplankton module's id_p(phy_i) pattern) so each species' biomass/traits can
   ! be examined separately.
   ALLOCATE(data%id_count(data%num_phytos))
   ALLOCATE(data%id_phyc(data%num_phytos))
   ALLOCATE(data%id_phyn(data%num_phytos))
   ALLOCATE(data%id_phyp(data%num_phytos))
   ALLOCATE(data%id_chl(data%num_phytos))
   ALLOCATE(data%id_mTopt(data%num_phytos))
   ALLOCATE(data%id_vTopt(data%num_phytos))
   ALLOCATE(data%id_mCDiv(data%num_phytos))
   ALLOCATE(data%id_vCDiv(data%num_phytos))
   ALLOCATE(data%id_mlnalpha(data%num_phytos))
   ALLOCATE(data%id_vlnalpha(data%num_phytos))
   ALLOCATE(data%id_cov_TA(data%num_phytos))
   ALLOCATE(data%id_cov_TL(data%num_phytos))
   ALLOCATE(data%id_cov_AL(data%num_phytos))
   ALLOCATE(data%id_N_birth(data%num_phytos))
   ALLOCATE(data%id_N_mutate(data%num_phytos))
   ALLOCATE(data%id_N_death(data%num_phytos))
   ALLOCATE(data%id_cum_birth(data%num_phytos))
   ALLOCATE(data%id_cum_death_thresh(data%num_phytos))
   ALLOCATE(data%id_cum_death_prob(data%num_phytos))
   ALLOCATE(data%id_cum_cells_decayed(data%num_phytos))
   ALLOCATE(data%id_NPPc(data%num_phytos))
   ALLOCATE(data%id_cells(data%num_phytos))
   ALLOCATE(data%id_d_nit(data%num_phytos))
   ALLOCATE(data%id_d_amm(data%num_phytos))
   ALLOCATE(data%id_d_pon(data%num_phytos))
   ALLOCATE(data%id_d_frp(data%num_phytos))
   ALLOCATE(data%id_d_pop(data%num_phytos))
   ALLOCATE(data%id_d_poc(data%num_phytos))
   ALLOCATE(data%id_d_oxy(data%num_phytos))

   DO grp_i=1,data%num_phytos
      !real      :: N_(nlev) = 0d0  !#MH  N_ in PIBM
      data%id_count(grp_i)  = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_total_count', '#', 'layer live particle count')

      !real      :: PHYC(nlev) = 0d0  !Eulerian concentration of phyto C, mean trait  for each layer
      data%id_phyc(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_C', 'mmol/m3', 'layer concentration of phyto particle C', rezero = .false.)

      !real      ::    PHY(nlev) = 0d0
      data%id_phyn(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_N', 'mmol/m3', 'layer concentration of phyto particle N')

      !real      ::    PHY(nlev) = 0d0
      data%id_phyp(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_P', 'mmol/m3', 'layer concentration of phyto particle P')

      !real      ::    CHL(nlev) = 0d0
      data%id_chl(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_Chl', 'mg/m3', 'layer concentration of phyto particle Chl')

      !real      ::   mTopt_(nlev) = 0d0
      data%id_mTopt(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_mean_Topt', 'degrees C/pmol C', 'carbon-weighted mean layer optimal temperature')

      !real      ::   vTopt_(nlev) = 0d0
      data%id_vTopt(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_var_Topt', 'degrees C/pmol C', 'carbon-weighted layer optimal temperature variance')

      !real      ::   mCDiv_(nlev) = 0d0                                                                                                 !#ML check units here
      data%id_mCDiv(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_mean_CDiv', 'log(pmol C/cell)/pmol C', 'carbon-weighted log mean layer cellular carbon content threshold for division')

      !real      ::   vCDiv_(nlev) = 0d0                                                                                                 !#ML check units here
      data%id_vCDiv(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_var_CDiv', 'log(pmol C/cell)/pmol C', 'carbon-weighted log layer cellular carbon content threshold for division variance')

      !real      ::   mlnalpha_(nlev) = 0d0
      data%id_mlnalpha(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_mean_lnalpha', 'log((W m-2)-1 (gChl molC)-1 d-1)/pmol C', 'carbon-weighted mean layer log slope of the P-I curve')

      !real      ::   vlnalpha_(nlev) = 0d0
      data%id_vlnalpha(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_var_lnalpha', 'log((W m-2)-1 (gChl molC)-1 d-1)/pmol C', 'carbon-weighted layer log slope of the P-I curve variance')

      !real      ::   cov_TA(nlev) = 0d0
      data%id_cov_TA(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_cov_TA', 'degrees C/pmol C * log(((W m-2)-1 (gChl molC)-1 d-1))/pmol C', 'covariance between Topt and ln(alpha Chl)')

      !real      ::   cov_TL(nlev) = 0d0
      data%id_cov_TL(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_cov_TL', 'degrees C/pmol C * log(pmol/cell)/pmol C', 'covariance between Topt and Cdiv')

      !real      ::   cov_AL(nlev) = 0d0
      data%id_cov_AL(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_cov_AL', 'log((W m-2)-1 (gChl molC)-1 d-1))/pmol C * log(pmol/cell)/pmol C', 'covariance between ln(alpha Chl) and Cdiv')

      !integer :: N_birth(nlev) = 0  !Number of birth events during one hour                    # ML this may cause problems because is defined as an integer
      data%id_N_birth(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_N_birth', 'number', 'number of births in a layer')

      !integer :: N_mutate(nlev)= 0  !Number of mutation events during one hour at each grid    # ML this may cause problems because is defined as an integer
      data%id_N_mutate(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_N_mutate', 'number', 'number of mutation events in a layer')

      !integer :: N_death(nlev) = 0  !Number of death events during one hour
      data%id_N_death(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_N_death', 'number', 'number of deaths in a layer')

      ! DIAGNOSTIC: cumulative (rezero=.false.) counters, since simulation start, for full
      ! attribution of particle-count change over time. Unlike N_birth/N_death/N_mutate above
      ! (which are hourly-reset rate diagnostics), these accumulate monotonically and so give
      ! an accurate running total regardless of output-save frequency.
      data%id_cum_birth(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_cum_birth', 'number', &
                                 'cumulative number of births in a layer since sim start', rezero = .false.)
      data%id_cum_death_thresh(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_cum_death_thresh', 'number', &
                                 'cumulative number of Cmin/N/P subsistence-threshold deaths in a layer since sim start', rezero = .false.)
      data%id_cum_death_prob(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_cum_death_prob', 'number', &
                                 'cumulative number of background probabilistic (mort_prob) deaths in a layer since sim start', rezero = .false.)
      ! DIAGNOSTIC: cumulative CELLS (not particle-death events) lost via mort_prob when
      ! simMortNumDecay=1 - see aed_particle_bgc_phyto_abm. Distinct units from
      ! cum_death_prob (event counts), so kept as its own diagnostic.
      data%id_cum_cells_decayed(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_cum_cells_decayed', 'number', &
                                 'cumulative number of cells lost to mort_prob num-decay in a layer since sim start', rezero = .false.)

      ! Additional diagnostics from PIBM Geider_Lag.F90
      ! Includes diagnostic outputs for fluxes to and from layers
      data%id_NPPc(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_NPPc', '(mg C m-3 d-1)', 'C-based phytoplankton production')
      data%id_cells(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_id_cells', 'number', 'cells/m3')
      data%id_d_nit(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_id_d_nit', 'mmol N/m3/day', 'daily flux of NO3 due to particles in a layer')
      data%id_d_amm(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_id_d_amm', 'mmol N/m3/day', 'daily flux of NH4 due to particles in a layer')
      data%id_d_pon(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_id_d_pon', 'mmol N/m3/day', 'daily flux of PON due to particles in a layer')
      data%id_d_frp(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_id_d_frp', 'mmol P/m3/day', 'daily flux of FRP due to particles in a layer')
      data%id_d_pop(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_id_d_pop', 'mmol P/m3/day', 'daily flux of POP due to particles in a layer')
      data%id_d_poc(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_id_d_poc', 'mmol C/m3/day', 'daily flux of POC due to particles in a layer')
      data%id_d_oxy(grp_i) = aed_define_diag_variable(TRIM(data%phytos(grp_i)%p_name)//'_oxy_flux', 'mmol O/m3/day','oxygen production')
   ENDDO

   ! id_IPAR stays a single scalar - accumulated layer PAR forcing, not a per-species quantity
   data%id_IPAR = aed_define_diag_variable('IPAR', 'umol/m2/s', 'daily integrated PAR at each depth')

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
   data%id_lvol  = aed_locate_global('layer_vol')
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
INTEGER   :: g   ! this particle's phyto group

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

! The total particle pool (ppid) is divided among data%num_phytos species by
! round-robin position (NOT contiguous blocks): particle i belongs to group
! ((i-1) mod num_phytos)+1. This is the only place position arithmetic is used to
! determine group membership - every particle's group is then stored explicitly in
! ptm_istat(GRP) and read from there everywhere else (splitting, recruitment,
! biology), so it survives arbitrary later reuse of particle slots.
!
! Round-robin (rather than contiguous blocks) matters because the host only
! activates a contiguous PREFIX of this pool as "initially active" (ptm_addparticles
! in glm_ptm.c); the remainder sits inactive as reserve capacity for later
! splitting/recruitment. Contiguous species blocks can coincide badly with that
! active/inactive split point (in the worst case, one whole species could start
! with zero active particles and never recover, since a split's
! duplicate_and_halve_particle copies the PARENT's group onto whatever free slot it
! finds, so an already-dominant species keeps colonizing another species' reserve
! slots). Round-robin guarantees every contiguous prefix of the pool contains a
! near-even mix of every species, regardless of where that active/inactive cutoff
! falls.
DO i = 1, ppid
   g = MOD(i-1, data%num_phytos) + 1
   p(i)%ptm_istat(GRP) = g

   ! Seed this particle's base state from its OWN group's CSV parameters (not
   ! whatever generic single-group value the host used as the aed_define_ptm_variable
   ! initial= default) before the Gaussian trait perturbation below draws around it.
   p(i)%ptm_state(data%ip_c)          = data%phytos(g)%X_cinit
   p(i)%ptm_state(data%ip_n)          = data%phytos(g)%X_ninit
   p(i)%ptm_state(data%ip_p)          = data%phytos(g)%X_pinit
   p(i)%ptm_state(data%ip_chl)        = data%phytos(g)%X_chlinit
   p(i)%ptm_state(data%ip_num)        = data%phytos(g)%n0
   p(i)%ptm_state(data%ip_Cdiv)       = data%phytos(g)%Cdiv
   p(i)%ptm_state(data%ip_Topt)       = data%phytos(g)%T_opt
   p(i)%ptm_state(data%ip_LnalphaChl) = data%phytos(g)%Lnalphachl

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
SUBROUTINE aed_split_particle_phyto_abm( data,ppid,p )

!ARGUMENTS
   CLASS (aed_phyto_abm_data_t),INTENT(in) :: data
   TYPE(aed_ptm_t), DIMENSION(:), INTENT(inout) :: p
   INTEGER,INTENT(inout) :: ppid

!

!LOCALS
! NOTE: the loop index over phyto groups is deliberately NOT named 'grp'. Fortran
! identifiers are case-insensitive, so a local named 'grp' silently shadows the
! module PARAMETER GRP (= the ptm_istat slot holding a particle's phyto group),
! turning every p(..)%ptm_istat(GRP) in this scope into an indexed read/write at
! whatever the local happens to hold. The same applies to the contained
! recruit_particle's dummy argument, hence 'rgrp' there.
INTEGER   :: i, j, ng
INTEGER   :: n_active, largest_i, n_split_this_call
AED_REAL  :: largest_num
LOGICAL   :: found_slot

! Adaptive (constant-particle-NUMBER) splitting (2026-08). SB's only splitting mode
! was constant-WEIGHT: n_split > 0 is a fixed cells-per-particle threshold, so the
! particle COUNT grows with the population and eventually hits max_particle_num, after
! which splits silently fail (no free slot found) and per-particle statistical weight
! ends up wildly uneven. n_split < 0 (see the sentinel handling in
! aed_phytoplankton_load_params) opts a group INTO this alternative: the threshold
! becomes that group's live-cell total divided by its fair share of particle slots, so
! WEIGHT floats and COUNT is held near capacity instead. n_split > 0 is entirely
! unaffected - this is additive, not a replacement.
INTEGER   :: g, slots_per_grp
AED_REAL, ALLOCATABLE :: grp_cells(:)   ! per-group total live cell count this call
INTEGER,  ALLOCATABLE :: grp_live(:)    ! per-group live PARTICLE count this call
AED_REAL  :: thresh
! A split produces TWO half-size particles, not one full-size one, so a population
! held exactly at the threshold settles to ~0.6x it, not to it - 0.45 (not the more
! obvious 0.9) targets ~75% block occupancy at steady state. The hard occupancy guard
! below is what actually GUARANTEES the count stays bounded, regardless of exactly
! what steady-state ratio this constant achieves.
AED_REAL, PARAMETER :: SPLIT_FILL = 0.45d0
! Safety cap on population-maintenance splits performed in a single call: each split
! requires an O(ppid) scan for the current largest particle plus an O(ppid) scan for a
! free slot, so an uncapped loop could in the worst case cost O(ppid^2) in one call.
! This still allows fast recovery (this routine runs once per host timestep) without
! that pathological cost.
INTEGER, PARAMETER :: MAX_POPMAINTAIN_SPLITS_PER_CALL = 1000

INTEGER   :: n_recruit_this_call
AED_REAL  :: rand01, recruit_height
! Biological timestep in days, taken from the same host-published value that
! aed_particle_bgc_phyto_abm uses for dtdays (aed_core%aed_ptm_dt_days). This
! routine has no column/layer context of its own, but the value is a module
! scalar set by aed_calculate_particles before the models run, so it is
! available here too. Previously hardcoded to 1/24 d, which assumed the host
! calls the particle model exactly once per hour.
AED_REAL :: dtdays_assumed

!

!-------------------------------------------------------------------------------
!BEGIN

! Same host-published biological timestep the BGC routine uses (see dtdays there).
dtdays_assumed = aed_ptm_dt_days


! Adaptive-splitting pre-pass: one O(ppid) sweep to get each group's live cell total
! and live particle count, needed below by any group using n_split < 0. Computed
! unconditionally (cheap; avoids doing it per-group inside the per-particle loop,
! which would cost O(ppid*num_phytos)) and is harmless for groups using the ordinary
! fixed n_split > 0 threshold, which never reads grp_cells/grp_live.
! slots_per_grp is SB's stand-in for Bejoy's n_per_grp: Bejoy derives its per-group
! slot budget from a CONTIGUOUS BLOCK particle layout (ppid/n_grp, one block per
! group). SB deliberately does NOT use contiguous blocks - particles are assigned
! round-robin (see the comment on aed_particle_initialize_phyto_abm) specifically so
! that splitting/recruitment don't have to search within a group's own block. The
! natural equivalent budget under round-robin is simply ppid divided evenly across
! groups; CEILING so a non-divisible ppid doesn't under-budget the group(s) that get
! the remainder.
ALLOCATE(grp_cells(data%num_phytos), grp_live(data%num_phytos))
grp_cells = zero_
grp_live  = 0
slots_per_grp = CEILING(REAL(ppid,KIND(thresh)) / REAL(data%num_phytos,KIND(thresh)))
DO i = 1, ppid
   IF (p(i)%ptm_istat(STAT) > 0) THEN
      g = p(i)%ptm_istat(GRP)
      grp_cells(g) = grp_cells(g) + p(i)%ptm_state(data%ip_num)
      grp_live(g)  = grp_live(g) + 1
   ENDIF
ENDDO

! (1) Per-particle size-triggered splitting: a particle that has grown past n_split
!     divides in two, same as cell division elsewhere, except it also allocates a new
!     active particle rather than doubling cell count in place. Gated by simSplit.
!     Checked per-particle against ITS OWN group's simSplit/n_split (not a single
!     hardcoded group), since different species can have different division thresholds.
DO i = 1, ppid
   g = p(i)%ptm_istat(GRP)
   IF (data%phytos(g)%simSplit > 0) THEN
      IF (data%phytos(g)%n_split < zero_) THEN
         ! ADAPTIVE: threshold is this group's current fair share, so it floats with
         ! the population instead of being fixed.
         thresh = grp_cells(g) / (REAL(slots_per_grp,KIND(thresh)) * SPLIT_FILL)
         ! Hard occupancy guard - this, not the ratio above, is what actually
         ! GUARANTEES the particle count stays bounded: once this group already
         ! fills its notional share of slots, refuse to split it further regardless
         ! of what the threshold computation says.
         IF (REAL(grp_live(g),KIND(thresh)) >= REAL(slots_per_grp,KIND(thresh)) * SPLIT_FILL) &
            thresh = HUGE(thresh)
      ELSE
         thresh = data%phytos(g)%n_split   ! fixed threshold (existing behaviour, unchanged)
      ENDIF
      IF (p(i)%ptm_istat(STAT) > 0 .AND. p(i)%ptm_state(data%ip_num) > thresh) THEN
         j = 1
         found_slot = .FALSE.
         DO
            IF (j > ppid) THEN
               print *, 'aed_split_particle(): WARNING no more available particles; skipping particle split'
               EXIT
            ELSE IF (p(j)%ptm_istat(STAT) == 0) THEN ! SCB change 7/13 - check for inactivity only
               CALL duplicate_and_halve_particle(i, j)
               found_slot = .TRUE.
               EXIT
            END IF
            j = j + 1
         END DO ! end search for available particles loop

         IF (found_slot) THEN
            ! DIAGNOSTIC: record the split event (one of the two pathways that create a
            ! new active particle / increase the particle count - see also (2) below).
            ptm_cum_split_count = ptm_cum_split_count + 1
            BLOCK
               INTEGER :: split_log_unit
               OPEN(NEWUNIT=split_log_unit, FILE='ptm_split_events.csv', &
                    POSITION='append', STATUS='unknown', ACTION='write')
               IF (ptm_cum_split_count == 1) WRITE(split_log_unit,'(A)') 'event_num,parent_particle,new_particle,cells_per_particle_after_split'
               WRITE(split_log_unit,'(I0,A,I0,A,I0,A,ES14.6)') ptm_cum_split_count, ',', i, ',', j, ',', p(i)%ptm_state(data%ip_num)
               CLOSE(split_log_unit)
            END BLOCK
            ! Keep grp_live in sync for the occupancy guard above: this split just
            ! created a new active particle in slot j, in the SAME group as parent i
            ! (duplicate_and_halve_particle copies the parent's group onto it).
            IF (data%phytos(g)%n_split < zero_) grp_live(g) = grp_live(g) + 1
         ELSE
            RETURN ! no free slots anywhere; further scanning this call is pointless
         END IF
      END IF ! end does particle need to split if statement
   END IF ! end this particle's group has simSplit enabled
END DO ! end particle loop
! grp_cells/grp_live are NOT deallocated here - block (4) below (reseeding) also
! reads grp_live. Deallocated once, after block (4), just before CONTAINS.

! (2) Population maintenance: if the domain-wide active particle count has fallen below
!     target_particle_num (namelist parameter, disabled when 0), repeatedly split the
!     largest currently-active particle - by cells per particle, ip_num - into a free
!     slot to replace losses, until the target is met, no free slots remain, or the
!     per-call cap above is reached. Independent of simSplit: this is a population-size
!     trigger, not a per-particle-size trigger.
IF (data%target_particle_num > 0) THEN
   n_active = 0
   DO i = 1, ppid
      IF (p(i)%ptm_istat(STAT) > 0) n_active = n_active + 1
   END DO

   n_split_this_call = 0
   DO WHILE (n_active < data%target_particle_num .AND. n_split_this_call < MAX_POPMAINTAIN_SPLITS_PER_CALL)
      ! find the largest active particle
      largest_i = -1
      largest_num = -1.0d0
      DO i = 1, ppid
         IF (p(i)%ptm_istat(STAT) > 0 .AND. p(i)%ptm_state(data%ip_num) > largest_num) THEN
            largest_num = p(i)%ptm_state(data%ip_num)
            largest_i = i
         END IF
      END DO
      IF (largest_i < 0 .OR. largest_num < 2.0d0) EXIT ! nothing left worth splitting

      ! find a free slot to hold the new particle
      j = 1
      found_slot = .FALSE.
      DO
         IF (j > ppid) THEN
            print *, 'aed_split_particle(): WARNING population below target_particle_num but no free particle slots available'
            EXIT
         ELSE IF (p(j)%ptm_istat(STAT) == 0) THEN
            found_slot = .TRUE.
            EXIT
         END IF
         j = j + 1
      END DO
      IF (.NOT. found_slot) EXIT ! no free slots anywhere; give up until next call

      CALL duplicate_and_halve_particle(largest_i, j)
      n_active = n_active + 1
      n_split_this_call = n_split_this_call + 1

      ! DIAGNOSTIC: record the split event, distinct from the n_split-triggered log above.
      ptm_cum_popmaintain_split_count = ptm_cum_popmaintain_split_count + 1
      BLOCK
         INTEGER :: pm_log_unit
         OPEN(NEWUNIT=pm_log_unit, FILE='ptm_popmaintain_split_events.csv', &
              POSITION='append', STATUS='unknown', ACTION='write')
         IF (ptm_cum_popmaintain_split_count == 1) WRITE(pm_log_unit,'(A)') &
            'event_num,active_count_before_this_split,target_particle_num,parent_particle,new_particle,cells_per_particle_after_split'
         WRITE(pm_log_unit,'(I0,A,I0,A,I0,A,I0,A,I0,A,ES14.6)') &
            ptm_cum_popmaintain_split_count, ',', n_active - 1, ',', data%target_particle_num, ',', &
            largest_i, ',', j, ',', p(largest_i)%ptm_state(data%ip_num)
         CLOSE(pm_log_unit)
      END BLOCK
   END DO
END IF

! (3) General recruitment: an ongoing, rate-controlled source of brand-new active
!     particles, independent of simSplit and target_particle_num above, and
!     independent of initial seeding / inflow-sourced particles (host-side C).
!     Represents an unmodeled external biological source (e.g. resuspension,
!     immigration). Gated by recruit_rate(g) (particles/day, disabled when 0).
!     Per-group (unlike population maintenance in (2) above): each species recruits
!     independently at its own rate, into its own height band.
DO ng = 1, data%num_phytos
   IF (data%recruit_rate(ng) > 0) THEN
      ! NOTE: same fragile-but-precedented dt assumption as elsewhere in this file
      ! (see dtdays in aed_particle_bgc_phyto_abm) - this routine has no column/
      ! layer context and so cannot read the model's actual timestep.
      recruit_accumulator(ng) = recruit_accumulator(ng) + data%recruit_rate(ng) * dtdays_assumed

      n_recruit_this_call = 0
      DO WHILE (recruit_accumulator(ng) >= 1.0d0 .AND. n_recruit_this_call < data%recruit_max_per_call(ng))

         ! find a free slot (any slot regardless of which group previously occupied it -
         ! correctness comes from recruit_particle setting GRP=ng on the new occupant)
         j = 1
         found_slot = .FALSE.
         DO
            IF (j > ppid) THEN
               print *, 'aed_split_particle(): WARNING recruit_rate > 0 but no free particle slots available'
               EXIT
            ELSE IF (p(j)%ptm_istat(STAT) == 0) THEN
               found_slot = .TRUE.
               EXIT
            END IF
            j = j + 1
         END DO
         IF (.NOT. found_slot) EXIT ! no free slots anywhere; keep the fractional
                                     ! accumulator so nothing is lost - try again next call

         ! ---- draw an independently-random height within this group's recruitment band,
         !      matching the ptm_addparticles() random-within-band precedent used
         !      for initial seeding (GLM/src/glm_ptm.c) ----
         CALL RANDOM_NUMBER(rand01)
         recruit_height = data%recruit_height_min(ng) + &
                           rand01 * (data%recruit_height_max(ng) - data%recruit_height_min(ng))

         CALL recruit_particle(j, ng, recruit_height)

         ! DIAGNOSTIC: record the recruitment event, distinct from the split-triggered
         ! logs above (this is the third, independent pathway that creates a new
         ! active particle).
         ptm_cum_recruit_count = ptm_cum_recruit_count + 1
         n_recruit_this_call = n_recruit_this_call + 1
         recruit_accumulator(ng) = recruit_accumulator(ng) - 1.0d0
         BLOCK
            INTEGER :: recruit_log_unit
            OPEN(NEWUNIT=recruit_log_unit, FILE='ptm_recruit_events.csv', &
                 POSITION='append', STATUS='unknown', ACTION='write')
            IF (ptm_cum_recruit_count == 1) WRITE(recruit_log_unit,'(A)') &
               'event_num,new_particle,group,height,cells_per_particle'
            WRITE(recruit_log_unit,'(I0,A,I0,A,I0,A,ES14.6,A,ES14.6)') &
               ptm_cum_recruit_count, ',', j, ',', ng, ',', recruit_height, ',', p(j)%ptm_state(data%ip_num)
            CLOSE(recruit_log_unit)
         END BLOCK
      END DO
   END IF
END DO

! (4) Reseeding: threshold-triggered population maintenance, complementary to (3)
!     above rather than a replacement for it (see the FEATURE note on
!     particle_reseed_min/target in the type definition for the distinction).
!     Gated by particle_reseed_min(g) (disabled when 0). Reuses grp_live from the
!     adaptive-splitting pre-pass at the top of this routine, so no extra O(ppid)
!     scan is needed to know each group's current live count.
DO ng = 1, data%num_phytos
   IF (data%particle_reseed_min(ng) > 0 .AND. grp_live(ng) < data%particle_reseed_min(ng)) THEN
      n_recruit_this_call = 0
      DO WHILE (grp_live(ng) < data%particle_reseed_target(ng) .AND. &
                n_recruit_this_call < data%recruit_max_per_call(ng))

         j = 1
         found_slot = .FALSE.
         DO
            IF (j > ppid) THEN
               print *, 'aed_split_particle(): WARNING particle_reseed_min > 0 but no free particle slots available'
               EXIT
            ELSE IF (p(j)%ptm_istat(STAT) == 0) THEN
               found_slot = .TRUE.
               EXIT
            END IF
            j = j + 1
         END DO
         IF (.NOT. found_slot) EXIT ! no free slots anywhere; try again next call

         ! Reseeded particles use the SAME height band as recruit_rate for this group
         ! (recruit_height_min/max) - there is no separate "initial seeding depth
         ! range" exposed on the Fortran side to draw from instead, and this group's
         ! own recruitment band is a reasonable place to reintroduce it.
         CALL RANDOM_NUMBER(rand01)
         recruit_height = data%recruit_height_min(ng) + &
                           rand01 * (data%recruit_height_max(ng) - data%recruit_height_min(ng))

         CALL recruit_particle(j, ng, recruit_height)

         ! DIAGNOSTIC: record the reseed event, distinct from the recruit-triggered
         ! log above (this is the fourth, independent pathway that creates a new
         ! active particle).
         ptm_cum_reseed_count = ptm_cum_reseed_count + 1
         n_recruit_this_call = n_recruit_this_call + 1
         grp_live(ng) = grp_live(ng) + 1
         BLOCK
            INTEGER :: reseed_log_unit
            OPEN(NEWUNIT=reseed_log_unit, FILE='ptm_reseed_events.csv', &
                 POSITION='append', STATUS='unknown', ACTION='write')
            IF (ptm_cum_reseed_count == 1) WRITE(reseed_log_unit,'(A)') &
               'event_num,new_particle,group,height,cells_per_particle'
            WRITE(reseed_log_unit,'(I0,A,I0,A,I0,A,ES14.6,A,ES14.6)') &
               ptm_cum_reseed_count, ',', j, ',', ng, ',', recruit_height, ',', p(j)%ptm_state(data%ip_num)
            CLOSE(reseed_log_unit)
         END BLOCK
      END DO
   END IF
END DO


DEALLOCATE(grp_cells, grp_live)

CONTAINS

   !###########################################################################
   ! Duplicate particle i's attributes into slot j, halve ip_num on both so the
   ! total cell count represented is conserved, and assign j a new PTID. Shared
   ! by both split pathways above; does not touch the diagnostic counters/logs,
   ! since each pathway logs differently.
   SUBROUTINE duplicate_and_halve_particle(i, j)
      INTEGER, INTENT(in) :: i, j
      AED_REAL :: div
      INTEGER  :: pid

      !first duplicate all attributes except PTID
      p(j)%ptm_istat(STAT) = p(i)%ptm_istat(STAT)
      p(j)%ptm_istat(IDX2) = p(i)%ptm_istat(IDX2)
      p(j)%ptm_istat(IDX3) = p(i)%ptm_istat(IDX3)
      p(j)%ptm_istat(LAYR) = p(i)%ptm_istat(LAYR)
      p(j)%ptm_istat(FLAG) = p(i)%ptm_istat(FLAG)
      p(j)%ptm_istat(GRP)  = p(i)%ptm_istat(GRP)  ! the duplicate/split-off particle is the same species as its parent
      p(j)%ptm_env(MASS) = p(i)%ptm_env(MASS)
      p(j)%ptm_env(DIAM) = p(i)%ptm_env(DIAM)
      p(j)%ptm_env(DENS) = p(i)%ptm_env(DENS)
      p(j)%ptm_env(VVEL) = p(i)%ptm_env(VVEL)
      p(j)%ptm_env(HGHT) = p(i)%ptm_env(HGHT)
      p(j)%ptm_state(data%ip_par) = p(i)%ptm_state(data%ip_par)
      p(j)%ptm_state(data%ip_tem) = p(i)%ptm_state(data%ip_tem)
      p(j)%ptm_state(data%ip_no3) = p(i)%ptm_state(data%ip_no3)
      p(j)%ptm_state(data%ip_nh4) = p(i)%ptm_state(data%ip_nh4)
      p(j)%ptm_state(data%ip_frp) = p(i)%ptm_state(data%ip_frp)
      p(j)%ptm_state(data%ip_c) = p(i)%ptm_state(data%ip_c)
      p(j)%ptm_state(data%ip_n) = p(i)%ptm_state(data%ip_n)
      p(j)%ptm_state(data%ip_p) = p(i)%ptm_state(data%ip_p)
      p(j)%ptm_state(data%ip_chl) = p(i)%ptm_state(data%ip_chl)
      p(j)%ptm_state(data%ip_num) = p(i)%ptm_state(data%ip_num)
      p(j)%ptm_state(data%ip_Cdiv) = p(i)%ptm_state(data%ip_Cdiv)
      p(j)%ptm_state(data%ip_Topt) = p(i)%ptm_state(data%ip_Topt)
      p(j)%ptm_state(data%ip_LnalphaChl) = p(i)%ptm_state(data%ip_LnalphaChl)
      p(j)%ptm_state(data%ip_rhochl) = p(i)%ptm_state(data%ip_rhochl)

      !then split the number of cells for both old and new particles
      p(j)%ptm_state(data%ip_num) = p(i)%ptm_state(data%ip_num) / 2
      p(i)%ptm_state(data%ip_num) = p(i)%ptm_state(data%ip_num) / 2

      !then get a new PTID for new particle
      IF(p(j)%ptm_istat(PTID) < 0) THEN
         p(j)%ptm_istat(PTID) = j;
      ELSE
         div = REAL(p(j)%ptm_istat(PTID) / ppid)
         pid = FLOOR(div)
         p(j)%ptm_istat(PTID) = ppid + pid*ppid + j
      END IF
   END SUBROUTINE duplicate_and_halve_particle

   !###########################################################################
   ! Activate free slot j as a brand-new recruited particle at target_height.
   ! Unlike duplicate_and_halve_particle, there is no parent - PTM host state and
   ! AED biology state are both populated fresh: host-state (LAYR/DIAM/DENS/VVEL)
   ! is copied from whichever other active particle's HGHT is closest to
   ! target_height (best-effort placeholder for LAYR - GLM's ptm_update_layerid(),
   ! called unconditionally every host step, corrects it before it's ever read by
   ! per-layer bgc; falls back to LAYR=1, DIAM/DENS/VVEL=0 if the domain is
   ! entirely empty of other active particles). Biology state (traits + cell C/N/P/
   ! Chl) is drawn the same way aed_particle_initialize_phyto_abm gives the
   ! founding population its diversity - a Gaussian perturbation (srand_mtGaus,
   ! sigma_init from module params) - except centered on the target group's own
   ! aed_phyto_pars.csv defaults (data%phytos(ng)%...) rather than on a stale
   ! ptm_state left behind by whatever previously occupied this slot.
   SUBROUTINE recruit_particle(j, rgrp, target_height)
      use params,  only : NInit, sigma_init
      use mGf90,   only : srand_mtGaus
      INTEGER, INTENT(in)  :: j
      INTEGER, INTENT(in)  :: rgrp   ! which data%phytos(:) group this recruit belongs to
      AED_REAL, INTENT(in) :: target_height
      INTEGER  :: k, best_k, v
      AED_REAL :: best_dist, d
      real     :: oldtt(1), newtt(1), vartt(1,1)
      integer, parameter :: iTopt=1, iSize=2, ialphaChl=3, iC=4, iNO3=5, iP=6, iChl=7, iNum=8

      ! ---- generic PTM host-state fields ----
      p(j)%ptm_istat(STAT) = 1
      p(j)%ptm_istat(IDX2) = 0
      p(j)%ptm_istat(IDX3) = 0
      p(j)%ptm_istat(FLAG) = 0
      p(j)%ptm_istat(GRP)  = rgrp
      p(j)%ptm_env(HGHT)   = target_height
      p(j)%ptm_env(MASS)   = 1.0d0   ! matches ptm_addparticles' hardcoded MASS

      best_k = -1
      best_dist = HUGE(best_dist)
      DO k = 1, ppid
         IF (k == j) CYCLE
         IF (p(k)%ptm_istat(STAT) > 0) THEN
            d = ABS(p(k)%ptm_env(HGHT) - target_height)
            IF (d < best_dist) THEN
               best_dist = d
               best_k = k
            END IF
         END IF
      END DO
      IF (best_k > 0) THEN
         p(j)%ptm_istat(LAYR) = p(best_k)%ptm_istat(LAYR)
         p(j)%ptm_env(DIAM)   = p(best_k)%ptm_env(DIAM)
         p(j)%ptm_env(DENS)   = p(best_k)%ptm_env(DENS)
         p(j)%ptm_env(VVEL)   = p(best_k)%ptm_env(VVEL)
      ELSE
         p(j)%ptm_istat(LAYR) = 1
         p(j)%ptm_env(DIAM)   = 0.0d0
         p(j)%ptm_env(DENS)   = 0.0d0
         p(j)%ptm_env(VVEL)   = 0.0d0
      END IF

      ! ---- PTID: same generation-aware scheme as duplicate_and_halve_particle ----
      IF (p(j)%ptm_istat(PTID) < 0) THEN
         p(j)%ptm_istat(PTID) = j
      ELSE
         p(j)%ptm_istat(PTID) = ppid + FLOOR(REAL(p(j)%ptm_istat(PTID)) / ppid)*ppid + j
      END IF

      ! ---- non-trait environment/state: neutral start, filled in by the next bgc call ----
      p(j)%ptm_state(data%ip_par) = 0.0d0
      p(j)%ptm_state(data%ip_tem) = 0.0d0
      p(j)%ptm_state(data%ip_no3) = 0.0d0
      p(j)%ptm_state(data%ip_nh4) = 0.0d0
      p(j)%ptm_state(data%ip_frp) = 0.0d0
      p(j)%ptm_state(data%ip_rhochl) = 0.0d0

      ! ---- traits + cell C/N/P/Chl/num: base CSV defaults, Gaussian-perturbed ----
      DO v = 1, NInit
         SELECT CASE (v)
         CASE (iTopt);     oldtt(1) = data%phytos(rgrp)%T_opt
         CASE (iSize);     oldtt(1) = log(data%phytos(rgrp)%Cdiv)
         CASE (ialphaChl); oldtt(1) = data%phytos(rgrp)%Lnalphachl
         CASE (iC);        oldtt(1) = log(data%phytos(rgrp)%X_cinit)
         CASE (iNO3);      oldtt(1) = log(data%phytos(rgrp)%X_ninit)
         CASE (iP);        oldtt(1) = log(data%phytos(rgrp)%X_pinit)
         CASE (iChl);      oldtt(1) = log(data%phytos(rgrp)%X_chlinit)
         CASE (iNum);      oldtt(1) = log(data%phytos(rgrp)%n0)
         END SELECT
         vartt(1,1) = sigma_init(v)**2
         newtt = srand_mtGaus(1, oldtt, vartt)
         SELECT CASE (v)
         CASE (iTopt);     p(j)%ptm_state(data%ip_Topt)       = newtt(1)
         CASE (iSize);     p(j)%ptm_state(data%ip_Cdiv)       = exp(newtt(1))
         CASE (ialphaChl); p(j)%ptm_state(data%ip_LnalphaChl) = newtt(1)
         CASE (iC);        p(j)%ptm_state(data%ip_c)          = exp(newtt(1))
         CASE (iNO3);      p(j)%ptm_state(data%ip_n)          = exp(newtt(1))
         CASE (iP);        p(j)%ptm_state(data%ip_p)          = exp(newtt(1))
         CASE (iChl);      p(j)%ptm_state(data%ip_chl)        = exp(newtt(1))
         CASE (iNum);      p(j)%ptm_state(data%ip_num)        = exp(newtt(1))
         END SELECT
      END DO
   END SUBROUTINE recruit_particle

END SUBROUTINE aed_split_particle_phyto_abm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_particle_bgc_phyto_abm( data,column,layer_idx,ppid,p )

!DECLARATIONS FROM PIBM Par2PHY SUBROUTINE
   use params,          only : NTrait, nu, sigma
   use mGf90,           only : srand_mtGaus
   USE Trait_functions, only : PHY_C2Vol

   IMPLICIT NONE

!ARGUMENTS
   CLASS (aed_phyto_abm_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   INTEGER,INTENT(inout) :: ppid
   TYPE(aed_ptm_t), DIMENSION(:), INTENT(inout) :: p

!LOCALS
   INTEGER :: n, status
   AED_REAL, ALLOCATABLE :: oxy_flux(:)   ! per-group
   AED_REAL :: decay, area, thickness

   AED_REAL :: Mu_max, WaterTemperature, Light, Depth, Kd, N_Limitation, P_Limitation, D0, D1, par, no3, nh4, frp, pw, mu, pw20, mu20
   AED_REAL :: f_T, Iz, f_I, Respiration, Mu_net

   AED_REAL, PARAMETER :: buoyancy_age = 86400.
   AED_REAL, PARAMETER :: DT           = 15.*60.
   AED_REAL, PARAMETER :: decay_rate   = 0.15
   AED_REAL, PARAMETER :: fres         = 0.7

! ADDITIONAL DECLARATIONS FROM PIBM
   ! PHYC/PHYN/PHYP/CHL and the uptake/mortality/oxy_flux accumulators below are one
   ! entry per simulated phyto group (data%num_phytos), rather than single layer-wide
   ! scalars, so that biomass/flux diagnostics can be reported per species. Allocated
   ! just below the declaration block, before use.
   real, allocatable :: PHYC(:)
   real, allocatable :: PHYN(:)
   real, allocatable :: PHYP(:)
   real, allocatable :: CHL(:)
   real      ::  nu_ = 0d0   !Basic Mutation rate
   real      :: cff = 0.d0   !Random number [0,1] for mutation
   real      :: rnm = 0.d0   !Random number [0,1] for mortality
   real      :: mort_prob_step = 0.d0   !mort_prob (d-1) converted to a per-call death probability via dtdays
   real      :: num_lost = 0.d0         !cells lost this call under simMortNumDecay=1 (see PHYTOPLANKTON MORTALITY)
   real      :: num_before = 0.d0       !DEBUG: ip_num at top of this particle's per-call update, to detect single-call collapses
   integer, save :: dbg_call_count = 0  !DEBUG: monotonic count of calls to this subroutine, for ordering debug output

   real      :: oldtt(1) = 0.   !Scratch variable for storing the old trait
   real      :: newtt(1) = 0.   !Scratch variable for storing the new trait
   real      :: vartt(1,1) = 0.   !Variance of the mutating trait
   integer :: k,i,m,ipar
   integer :: g, ng   ! g = current particle's phyto group (p(i)%ptm_istat(GRP)); ng = loop index over groups
   integer, allocatable :: N_by_group(:)      !Number of active particles per group in this layer
   real,    allocatable :: Abun_by_group(:)   !Total cell abundance per group in this layer
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
   real    :: ZOO(NZOO) = 0.
   real    :: DET = 0.
   real    :: tf_z  = 0.
   real    :: Graz  = 0.
   real    :: dC_   = 0.
   real    :: dN_   = 0.
   real    :: dNO3_   = 0.
   real    :: dNH4_   = 0.
   real    :: dP_   = 0.
   real    :: dChl_ = 0.
   real    :: ESD_ = 0.    ! equivalent spherical diameter (micron)
   real, allocatable :: uptake(:)   !Per-group DIN uptake
   real, allocatable :: uptake_NO3(:)   !Per-group NO3 uptake
   real, allocatable :: uptake_NH4(:)   !Per-group NH4 uptake
   real, allocatable :: uptake_P(:)   !Per-group FRP uptake
   real    :: pp_DZ = 0.
   real    :: pp_DZP = 0.
   real    :: pp_ND = 0.
   real    :: pp_PD = 0.
   !# Per-particle physiology increments from the GMK98 pass, held so that the
   !# aggregate nutrient demand of every particle in this layer can be compared
   !# against what the layer actually holds, and rationed down if it exceeds it,
   !# BEFORE any of it is committed. See the rationing block below.
   real, allocatable :: sv_dC(:), sv_dN(:), sv_dNO3(:), sv_dNH4(:)
   real, allocatable :: sv_dP(:), sv_dChl(:), sv_ESD(:)
   real    :: dem_NO3, dem_NH4, dem_FRP
   real    :: scl_NO3, scl_NH4, scl_FRP, scl_N
   real, allocatable :: Pmort(:)
   real, allocatable :: Pmort_P(:)
   real, allocatable :: Pmort_C(:)
   real    :: Cmin  = 0. !Phytoplankton subsistence carbon quota below which the cell will die
   real    :: FZoo(NZOO) = 0.   !The total amount of palatable prey (in Nitrogen)
                             !for each zooplankton size class
   real    :: phyV = 0.   !Phytoplankton cell volume
   real    :: RES     = 0. !ML not sure why this is included here because respiration is already accounted for in dN of GMK function, but keeping for posterity, always 0 for now
   real    :: RES_NO3     = 0. !ML not sure why this is included here because respiration is already accounted for in dN of GMK function, but keeping for posterity, always 0 for now
   real    :: RES_NH4     = 0. !ML not sure why this is included here because respiration is already accounted for in dN of GMK function, but keeping for posterity, always 0 for now
   real    :: RES_P   = 0.  !ML not sure why this is included here because respiration is already accounted for in dN of GMK function, but keeping for posterity, always 0 for now
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
   !# Biological timestep in days. Assigned from aed_ptm_dt_days at !BEGIN below
   !# (published by aed_calculate_particles from the value the host passes in),
   !# so this tracks the interval the host is actually advancing over rather
   !# than assuming an hourly call. The initialiser is only a fallback.
   real     :: dtdays       = 1.0d0/24
   real  :: hour_of_day   != 0
   integer  :: integer_day   != 0
   real  :: real_day
   ! ML new declarations to get things working with PIBM and AED
   integer  :: chk_lyr_new  = 0.
   integer  :: chk_lyr_old  = 0.
   real     :: p_vvel = 0. !vertical velocity calculated in AED (m/s)
   ! min_rho/max_rho were declared here as locals hardcoded to 900./1200., which shadowed
   ! the identically-named namelist parameters and made those parameters inert. The clamp
   ! below now reads data%min_rho / data%max_rho directly (set in aed_define_phyto_abm).
   real     :: dens_flux = 0.

!End of declaration
!
!-------------------------------------------------------------------------------
!BEGIN

!------------------------------------------------------
! INITIAL CODE FROM AED

   ! Allocate this call's per-group accumulators/diagnostics-staging arrays (freshly
   ! each call, mirroring the existing index_/BN/BC/BP allocate-per-call pattern below).

   !# Take this call's biological timestep from the host (published via
   !# aed_calculate_particles). Was hardcoded to 1/24 d, which silently assumed
   !# the host calls us exactly once per hour - wrong the moment the particle
   !# model is run once per split_factor sub-step rather than once per host step.
   dtdays = aed_ptm_dt_days

   ALLOCATE(PHYC(data%num_phytos));  PHYC = 0d0
   ALLOCATE(PHYN(data%num_phytos));  PHYN = 0d0
   ALLOCATE(PHYP(data%num_phytos));  PHYP = 0d0
   ALLOCATE(CHL(data%num_phytos));   CHL  = 0d0
   ALLOCATE(uptake(data%num_phytos));      uptake      = 0d0
   ALLOCATE(uptake_NO3(data%num_phytos));  uptake_NO3  = 0d0
   ALLOCATE(uptake_NH4(data%num_phytos));  uptake_NH4  = 0d0
   ALLOCATE(uptake_P(data%num_phytos));    uptake_P    = 0d0
   ALLOCATE(Pmort(data%num_phytos));   Pmort   = 0d0
   ALLOCATE(Pmort_P(data%num_phytos)); Pmort_P = 0d0
   ALLOCATE(Pmort_C(data%num_phytos)); Pmort_C = 0d0
   ALLOCATE(oxy_flux(data%num_phytos)); oxy_flux = 0d0
   ALLOCATE(N_by_group(data%num_phytos));    N_by_group    = 0
   ALLOCATE(Abun_by_group(data%num_phytos)); Abun_by_group = 0d0

   ! Check if we are in a new cell, to reset cumulative counters
   chk_lyr_new = layer_idx
   IF ((chk_lyr_new > chk_lyr_old) .or. chk_lyr_new == 1) THEN
      chk_lyr_old = chk_lyr_new
   ENDIF

   ! Local environmental conditions in this layer
   WaterTemperature = _STATE_VAR_(data%id_tem) !22  !water temperature
   par =  _STATE_VAR_(data%id_par) ! _STATE_VAR_S_(data%id_I0) local photosynth. active radiation
   no3 = _STATE_VAR_(data%id_nit)        ! local no3
   nh4 = _STATE_VAR_(data%id_amm)        ! local nh4
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
       DO ng=1,data%num_phytos
          _DIAG_VAR_(data%id_NPPc(ng)) = zero_  !Reset NPPc_
       ENDDO
       _DIAG_VAR_(data%id_IPAR) = zero_  !Reset IPAR
   endif

   ! initialize flux diagnostics at 0 (one entry per phyto group)
   DO ng=1,data%num_phytos
      _DIAG_VAR_(data%id_d_pon(ng)) = zero_
      _DIAG_VAR_(data%id_d_nit(ng)) = zero_
      _DIAG_VAR_(data%id_d_amm(ng)) = zero_
      _DIAG_VAR_(data%id_d_pop(ng)) = zero_
      _DIAG_VAR_(data%id_d_frp(ng)) = zero_
      _DIAG_VAR_(data%id_d_poc(ng)) = zero_
      _DIAG_VAR_(data%id_d_oxy(ng)) = zero_
   ENDDO

   dbg_call_count = dbg_call_count + 1   !DEBUG

   !Calculate the number of super-individuals (N_) in this vertical layer and obtain their indices (index_)
   N_PAR = ppid
   N_ = 0

   !Get the indices of particles in this grid (for grazing loss)
   allocate(index_(0), stat=AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating index_***"

   do i = 1, N_PAR
      if (_PTM_STAT_(i,STAT) == 1) then !Ignore dead super-individuals ! SCB - ADD water column flag check here?
         N_ = N_ + 1
         N_by_group(p(i)%ptm_istat(GRP)) = N_by_group(p(i)%ptm_istat(GRP)) + 1
         if (N_ == 1) then
             allocate(scratch(1), stat=Allocatestatus)
             IF (AllocateStatus /= 0) STOP "*** Problem in allocating scratch***"
             scratch(1) = i ! _PTM_STAT_(i,STAT) ! SCB - note that Claude identified this as a likely bug - this was the recommended fix

            !Update index_ and scratch has been deallocated
            call move_alloc(scratch, index_)
         else
             allocate(scratch(size(index_) + 1), stat=Allocatestatus)
             IF (AllocateStatus /= 0) STOP "*** Problem in allocating scratch***"
             scratch(1:size(index_)) = index_
             scratch(size(index_) + 1) = i

            !Update index_ and scratch has been deallocated
            call move_alloc(scratch, index_)
         endif
      endif
   enddo
   ! get volume of layer
   !# TRUE layer volume from the host's hypsograph. Was dz*area, which treats the
   !# layer as a prism of constant area and so mis-scales every particle<->water
   !# exchange; Hz appears in the denominator of every flux conversion below.
   Hz = _STATE_VAR_(data%id_lvol)

   !Save number of super-individuals per layer, per group
   DO ng=1,data%num_phytos
      _DIAG_VAR_(data%id_count(ng)) = N_by_group(ng)
   ENDDO

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

         !Count the number of cells in each layer, per group
         Abun_by_group(p(i)%ptm_istat(GRP)) = Abun_by_group(p(i)%ptm_istat(GRP)) + p(i)%ptm_state(data%ip_num)

      ENDDO !End of looping through superindividuals
   ENDIF ! if N_ > 0

   ! convert abundance to cells per m-3 for diagnostic, per group
   DO ng=1,data%num_phytos
      _DIAG_VAR_(data%id_cells(ng)) = Abun_by_group(ng) / Hz !Abundances (cells m-3)
   ENDDO

   ! Calculate total phytoplankton nitrogen uptake, mortality, and PP (must after calculation of num(t+dt))
   ! Start by setting everything to 0 (already allocated/zeroed above, but keep the
   ! reset here too since this section can be re-entered across host sub-steps)
   uptake  = 0d0
   uptake_NO3  = 0d0
   uptake_NH4  = 0d0
   uptake_P  = 0d0
   Pmort   = 0d0
   Pmort_P   = 0d0
   Pmort_C   = 0d0
   oxy_flux = 0d0

   !###########################################################################
   !# PASS 1 of 2: compute each particle's physiology increments WITHOUT
   !# committing them, so the layer's aggregate nutrient demand is known before
   !# any of it is taken. Pass 2 (the loop further below) applies the results,
   !# scaled by the rationing factors computed in between.
   !#
   !# Every particle's uptake is a function of the SAME start-of-step no3/nh4/frp
   !# - no particle can see what the others are taking. Michaelis-Menten throttles
   !# each individual draw but places no bound on the sum, so with thousands of
   !# particles in a layer the aggregate can exceed what is present. The pool then
   !# goes negative and repair_state clamps it back to zero, fabricating mass.
   !#
   !# This mirrors aed_phytoplankton, which has had the same guard for its
   !# Eulerian groups since long before this module existed (see the
   !# "Check uptake values for availability to prevent -ve numbers" block, where
   !# uptake is scaled by available*0.99/upTot when the total exceeds supply).
   !# That module is naturally two-pass - it fills nuptake/puptake arrays, rations,
   !# then applies - whereas this one computed and committed in a single sweep,
   !# which is why the guard was never inherited. Same job, same treatment.
   !###########################################################################
   if (N_ > 0) then
      allocate(sv_dC(N_), sv_dN(N_), sv_dNO3(N_), sv_dNH4(N_), &
               sv_dP(N_), sv_dChl(N_), sv_ESD(N_), stat=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem allocating particle uptake scratch ***"

      dem_NO3 = 0d0 ; dem_NH4 = 0d0 ; dem_FRP = 0d0

      do j = 1, N_
         i = index_(j)
         g = p(i)%ptm_istat(GRP)

         ! set particle par, temp, and no3 using environmental conditions from this layer
         p(i)%ptm_state(data%ip_par) = par
         p(i)%ptm_state(data%ip_tem) = WaterTemperature
         p(i)%ptm_state(data%ip_no3) = no3
         p(i)%ptm_state(data%ip_nh4) = nh4
         p(i)%ptm_state(data%ip_frp) = frp

         call GMK98_Ind_TempSizeLight(p(i)%ptm_state(data%ip_tem),              &
                                      p(i)%ptm_state(data%ip_par),              &
                                      p(i)%ptm_state(data%ip_no3) + p(i)%ptm_state(data%ip_nh4), &
                                      p(i)%ptm_state(data%ip_no3),              &
                                      p(i)%ptm_state(data%ip_nh4),              &
                                      p(i)%ptm_state(data%ip_frp),              &
                                      p(i)%ptm_state(data%ip_c),                &
                                      p(i)%ptm_state(data%ip_n),                &
                                      p(i)%ptm_state(data%ip_p),                &
                                      p(i)%ptm_state(data%ip_chl),              &
                                      p(i)%ptm_state(data%ip_Topt),             &
                                      p(i)%ptm_state(data%ip_cdiv),             &
                                      exp(p(i)%ptm_state(data%ip_LnalphaChl)),  &
                                      dC_, dN_, dNO3_, dNH4_, dP_, dChl_, ESD_, &
                                      data%phytos(g)%RC,                        &
                                      data%phytos(g)%RN,                        &
                                      data%phytos(g)%RP,                        &
                                      data%phytos(g)%RChl,                      &
                                      data%phytos(g)%zeta_N,                    &
                                      data%phytos(g)%zeta_P,                    &
                                      data%phytos(g)%a1,                        &
                                      data%phytos(g)%R_growth,                  &
                                      data%phytos(g)%nx,                        &
                                      data%phytos(g)%thetaNmax,                 &
                                      data%phytos(g)%thetaPmax,                 &
                                      data%phytos(g)%QNmin_a,                   &
                                      data%phytos(g)%QNmin_b,                   &
                                      data%phytos(g)%QNmax_a,                   &
                                      data%phytos(g)%QNmax_b,                   &
                                      data%phytos(g)%QPmin_a,                   &
                                      data%phytos(g)%QPmin_b,                   &
                                      data%phytos(g)%QPmax_a,                   &
                                      data%phytos(g)%QPmax_b,                   &
                                      data%phytos(g)%KN_a,                      &
                                      data%phytos(g)%KN_b,                      &
                                      data%phytos(g)%KPho_a,                    &
                                      data%phytos(g)%KPho_b,                    &
                                      data%phytos(g)%a_c2vol,                   &
                                      data%phytos(g)%b_c2vol,                   &
                                      data%phytos(g)%a_pmax,                    &
                                      data%phytos(g)%rho,                       &
                                      data%phytos(g)%rho_star,                  &
                                      data%phytos(g)%b_rho,                     &
                                      data%phytos(g)%V_s,                       &
                                      data%phytos(g)%Ea0,                       &
                                      data%phytos(g)%Ed0,                       &
                                      data%phytos(g)%Ei,                        &
                                      data%phytos(g)%beta,                      &
                                      data%phytos(g)%phi,                       &
                                      data%phytos(g)%T_std,                     &
                                      data%phytos(g)%Tau,                       &
                                      data%phytos(g)%Beta_Ainf,                 &
                                      data%phytos(g)%Kappa,                     &
                                      data%phytos(g)%Kd_Ainf,                   &
                                      data%phytos(g)%a_,                        &
                                      data%phytos(g)%b_,                        &
                                      data%phytos(g)%v_,                        &
                                      p(i)%ptm_state(data%ip_rhochl))

         sv_dC(j)=dC_ ; sv_dN(j)=dN_ ; sv_dNO3(j)=dNO3_ ; sv_dNH4(j)=dNH4_
         sv_dP(j)=dP_ ; sv_dChl(j)=dChl_ ; sv_ESD(j)=ESD_

         !# Accumulate only NET REMOVAL from the water (positive d* = uptake).
         !# Negative values are a net return to the pool and can never exhaust it,
         !# so they must not offset another particle's demand in the sum.
         if (dNO3_ > 0d0) dem_NO3 = dem_NO3 + dNO3_ * p(i)%ptm_state(data%ip_num)
         if (dNH4_ > 0d0) dem_NH4 = dem_NH4 + dNH4_ * p(i)%ptm_state(data%ip_num)
         if (dP_   > 0d0) dem_FRP = dem_FRP + dP_   * p(i)%ptm_state(data%ip_num)
      enddo

      !# RATION. Convert demand from pmol/day for the layer into the mmol/m3 it
      !# would actually remove over this timestep (the same 1d-9/Hz*dtdays the
      !# commit below applies), then cap it at 99% of what is present - the same
      !# 0.99 headroom aed_phytoplankton uses. Scaling is proportional, so every
      !# particle keeps its share of a pool it cannot collectively exhaust.
      scl_NO3 = 1d0 ; scl_NH4 = 1d0 ; scl_FRP = 1d0
      dem_NO3 = dem_NO3 * 1d-9 / Hz * dtdays
      dem_NH4 = dem_NH4 * 1d-9 / Hz * dtdays
      dem_FRP = dem_FRP * 1d-9 / Hz * dtdays
      if (dem_NO3 > 1d-30 .and. dem_NO3 > max(no3,0d0)) scl_NO3 = max(no3,0d0)*0.99d0/dem_NO3
      if (dem_NH4 > 1d-30 .and. dem_NH4 > max(nh4,0d0)) scl_NH4 = max(nh4,0d0)*0.99d0/dem_NH4
      if (dem_FRP > 1d-30 .and. dem_FRP > max(frp,0d0)) scl_FRP = max(frp,0d0)*0.99d0/dem_FRP

      !# Apply the scaling to the stored increments. Only positive (uptake) terms
      !# are scaled; a negative d* is a return to the water and passes through.
      !# dN must stay equal to dNO3+dNH4 - the cell's nitrogen gain has to match
      !# what was actually removed from the two inorganic pools, or the particle
      !# and the water disagree and mass is created. dC/dChl are scaled by the
      !# same factor as the nitrogen that limits them, keeping C:N:Chl coherent.
      do j = 1, N_
         if (sv_dNO3(j) > 0d0) sv_dNO3(j) = sv_dNO3(j) * scl_NO3
         if (sv_dNH4(j) > 0d0) sv_dNH4(j) = sv_dNH4(j) * scl_NH4
         if (sv_dP(j)   > 0d0) sv_dP(j)   = sv_dP(j)   * scl_FRP
         scl_N = 1d0
         if (sv_dN(j) > 0d0) then
            scl_N = (sv_dNO3(j) + sv_dNH4(j)) / sv_dN(j)
            sv_dN(j) = sv_dNO3(j) + sv_dNH4(j)
         endif
         if (sv_dC(j)   > 0d0) sv_dC(j)   = sv_dC(j)   * min(min(scl_N, scl_FRP), 1d0)
         if (sv_dChl(j) > 0d0) sv_dChl(j) = sv_dChl(j) * min(scl_N, 1d0)
      enddo
   endif

   !###########################################################################
   !# PASS 2 of 2: commit the (possibly rationed) increments.
   !###########################################################################
   if (N_ > 0) then
      do j = 1, N_
         i = index_(j)
         g = p(i)%ptm_istat(GRP)   ! this particle's phyto group - selects which data%phytos(g) parameter set to use
         num_before = p(i)%ptm_state(data%ip_num)   !DEBUG

         dC_ = sv_dC(j) ; dN_ = sv_dN(j) ; dNO3_ = sv_dNO3(j) ; dNH4_ = sv_dNH4(j)
         dP_ = sv_dP(j) ; dChl_ = sv_dChl(j) ; ESD_ = sv_ESD(j)


         ! cumulate N and P uptake and oxy flux after running physiology function (per group)
         uptake(g)   =   uptake(g) + dN_ * p(i)%ptm_state(data%ip_num) ! Unit: pmol N d-1
         uptake_NO3(g)   =   uptake_NO3(g) + dNO3_ * p(i)%ptm_state(data%ip_num) ! Unit: pmol N d-1
         uptake_NH4(g)   =   uptake_NH4(g) + dNH4_ * p(i)%ptm_state(data%ip_num) ! Unit: pmol N d-1
         uptake_P(g) =   uptake_P(g) + dP_ * p(i)%ptm_state(data%ip_num) ! Unit: pmol P d-1
         oxy_flux(g) =   oxy_flux(g) + dC_ * p(i)%ptm_state(data%ip_num) ! Unit: pmol C d-1 (assuming 1:1 stoichiometry with O2)
         !NPPc_(k) = NPPc_(k) + dC_ * p_PHY(i)%num *1d-9/Hz(k)*12.d0*dtdays !Unit: mgC m-3 d-1 !ML leaving this PIBM code because of questions below re: NPP diagnostic
         _DIAG_VAR_(data%id_NPPc(g)) = _DIAG_VAR_(data%id_NPPc(g)) + dC_ * p(i)%ptm_state(data%ip_num) * 1d-9 / Hz * 12.d0 * dtdays !Unit: mgC m-3 d-1 ML why multiplied by 12? and also multiplying by dtdays means this is hr-1 (timestep) not d-1

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
            dens_flux = data%phytos(g)%c1 * &
            (one_ - EXP(-par/data%phytos(g)%I_K) ) - data%phytos(g)%c3
         ELSE
            ! darkness
            dens_flux = -data%phytos(g)%c3
         ENDIF

         ! Update density
         p(i)%ptm_env(DENS) = p(i)%ptm_env(DENS) + dens_flux

         ! check maximum/minimum density are not exceeded
         IF( p(i)%ptm_env(DENS)>data%max_rho ) THEN
            p(i)%ptm_env(DENS)=data%max_rho
         ENDIF
         IF(p(i)%ptm_env(DENS)<data%min_rho ) THEN
            p(i)%ptm_env(DENS)=data%min_rho
         ENDIF

         !# PHYTOPLANKTON VVEL updates based on specified settling mode
         SELECT CASE (data%phytos(g)%settling)

            CASE ( _MOB_OFF_ )
               ! disable settling by setting vertical velocity to 0
               p(i)%ptm_env(VVEL) = zero_

            CASE ( _MOB_CONST_ )
               ! constant settling velocity using user provided value
               p(i)%ptm_env(VVEL) = data%phytos(g)%w_p

            CASE ( _MOB_TEMP_ )
               ! constant settling velocity @20C corrected for density changes
               pw   = _STATE_VAR_(data%id_dens)
               mu   = water_viscosity(WaterTemperature)
               mu20 = 0.001002  ! N s/m2
               pw20 = 998.2000  ! kg/m3 (assuming freshwater)
               p(i)%ptm_env(VVEL) = data%phytos(g)%w_p*mu20*pw / ( mu*pw20 )

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
         ! If cellular carbon is lower than the subsistence threshold (Cmin), it dies:
         Cmin = 0.25d0 * p(i)%ptm_state(data%ip_cdiv) ! ML why multiplied by 0.25? this is from PIBM

         if (p(i)%ptm_state(data%ip_c) < Cmin .or. &
            (p(i)%ptm_state(data%ip_num)*p(i)%ptm_state(data%ip_n)) < Nt_min .or. &
            (p(i)%ptm_state(data%ip_num)*p(i)%ptm_state(data%ip_p)) < P_min .or. &
            p(i)%ptm_state(data%ip_num) < 1.0d0) then  ! The superindividual Dies
         !if (p(i)%ptm_state(data%ip_c) < Cmin .or. (p(i)%ptm_state(data%ip_num)*p(i)%ptm_state(data%ip_n)) < Nt_min .or. (p(i)%ptm_state(data%ip_num)*p(i)%ptm_state(data%ip_p)) < P_min) then  ! The superindividual Dies
            _DIAG_VAR_(data%id_N_death(g)) = _DIAG_VAR_(data%id_N_death(g)) + 1
            _DIAG_VAR_(data%id_cum_death_thresh(g)) = _DIAG_VAR_(data%id_cum_death_thresh(g)) + 1
            Pmort(g) = Pmort(g) + p(i)%ptm_state(data%ip_n) * p(i)%ptm_state(data%ip_num) !Natural mortality of phytoplankton ==> DET
            Pmort_P(g) = Pmort_P(g) + p(i)%ptm_state(data%ip_p) * p(i)%ptm_state(data%ip_num) !Natural mortality of phytoplankton ==> DET
            Pmort_C(g) = Pmort_C(g) + p(i)%ptm_state(data%ip_c) * p(i)%ptm_state(data%ip_num) !Natural mortality of phytoplankton ==> DET

            p(i)%ptm_state(data%ip_c)   = 0d0
            p(i)%ptm_state(data%ip_n)   = 0d0
            p(i)%ptm_state(data%ip_p)   = 0d0
            p(i)%ptm_state(data%ip_chl) = 0d0
            p(i)%ptm_state(data%ip_num) = 0d0
            _PTM_STAT_(i,STAT) = 0
            _PTM_STAT_(i,FLAG) = 3
         endif

         ! Now probabilistic mortality to compensate because we don't have zoops, disease, or senescence in ABM
         ! mort_prob is a per-day rate (d-1), consistent with every other rate parameter in this
         ! module (R_resp, RC, RN, RP, RChl, etc.); convert it to the probability of death over
         ! THIS call (which covers dtdays of simulated time) via the usual exponential-decay
         ! rate-to-probability conversion, rather than comparing the raw daily rate directly.
         mort_prob_step = 1.d0 - exp(-data%phytos(g)%mort_prob * dtdays)

         IF (data%phytos(g)%simMortNumDecay > 0) THEN
            ! ---- num-decay mode: mort_prob shrinks the number of cells the particle
            ! represents (ip_num) by the same exponential fraction each call, instead of
            ! probabilistically killing the whole particle. Per-cell state (ip_c/n/p/chl)
            ! is unchanged - those describe the surviving cells, not the ones lost. The
            ! particle itself only fully dies via the Cmin/Nt_min/P_min subsistence check
            ! above, once its shrinking total N/P/C eventually crosses that threshold.
            IF (_PTM_STAT_(i,STAT) == 1) THEN
               num_lost = p(i)%ptm_state(data%ip_num) * mort_prob_step
               _DIAG_VAR_(data%id_cum_cells_decayed(g)) = _DIAG_VAR_(data%id_cum_cells_decayed(g)) + num_lost
               Pmort(g)   = Pmort(g)   + p(i)%ptm_state(data%ip_n) * num_lost !Natural mortality of phytoplankton ==> DET
               Pmort_P(g) = Pmort_P(g) + p(i)%ptm_state(data%ip_p) * num_lost !Natural mortality of phytoplankton ==> DET
               Pmort_C(g) = Pmort_C(g) + p(i)%ptm_state(data%ip_c) * num_lost !Natural mortality of phytoplankton ==> DET

               p(i)%ptm_state(data%ip_num) = p(i)%ptm_state(data%ip_num) - num_lost
            END IF
         ELSE
            ! ---- default mode: a per-call Bernoulli draw kills the whole particle.
            call random_number(rnm)

            ! Guard against double-counting a particle that already died from the subsistence
            ! threshold check above (STAT was set to 0 there); without this guard, a particle
            ! could be recorded as dying twice in the same call (harmless to mass balance, since
            ! its state was already zeroed, but it inflates the death-cause diagnostics below).
            IF (_PTM_STAT_(i,STAT) == 1 .AND. rnm < mort_prob_step) THEN !it dies
               _DIAG_VAR_(data%id_N_death(g)) = _DIAG_VAR_(data%id_N_death(g)) + 1.
               _DIAG_VAR_(data%id_cum_death_prob(g)) = _DIAG_VAR_(data%id_cum_death_prob(g)) + 1.
               Pmort(g) = Pmort(g) + p(i)%ptm_state(data%ip_n) * p(i)%ptm_state(data%ip_num) !Natural mortality of phytoplankton ==> DET
               Pmort_P(g) = Pmort_P(g) + p(i)%ptm_state(data%ip_p) * p(i)%ptm_state(data%ip_num) !Natural mortality of phytoplankton ==> DET
               Pmort_C(g) = Pmort_C(g) + p(i)%ptm_state(data%ip_c) * p(i)%ptm_state(data%ip_num) !Natural mortality of phytoplankton ==> DET

               p(i)%ptm_state(data%ip_c)   = 0d0
               p(i)%ptm_state(data%ip_n)   = 0d0
               p(i)%ptm_state(data%ip_p)   = 0d0
               p(i)%ptm_state(data%ip_chl) = 0d0
               p(i)%ptm_state(data%ip_num) = 0d0
               _PTM_STAT_(i,STAT) = 0
               _PTM_STAT_(i,FLAG) = 3
            ENDIF
         END IF

         !DEBUG: flag any single-call collapse of ip_num from a meaningful population down to <1 cell
         IF (num_before > 1.0d0 .AND. p(i)%ptm_state(data%ip_num) < 1.0d0) THEN
            print *, 'PTM_NUM_COLLAPSE call=', dbg_call_count, ' yday=', real_day, &
                     ' ptid=', p(i)%ptm_istat(PTID), ' i=', i, ' layer=', layer_idx, &
                     ' num_before=', num_before, ' num_after=', p(i)%ptm_state(data%ip_num), &
                     ' num_lost=', num_lost, ' mort_prob_step=', mort_prob_step, &
                     ' simMortNumDecay=', data%phytos(g)%simMortNumDecay, &
                     ' mort_prob=', data%phytos(g)%mort_prob, ' dtdays=', dtdays, &
                     ' ip_c=', p(i)%ptm_state(data%ip_c), ' ip_n=', p(i)%ptm_state(data%ip_n), &
                     ' ip_p=', p(i)%ptm_state(data%ip_p), ' Cmin=', Cmin, &
                     ' WaterTemp=', WaterTemperature, ' par=', par
         END IF
      enddo !End of looping through particles
  endif !End if of N_ > 0

  ! Convert all updates to be per m3 per second
  uptake = uptake*1d-9/Hz/secs_per_day !Convert uptake to mmol N m-3 s-1 (DIN)
  uptake_NO3 = uptake_NO3*1d-9/Hz/secs_per_day !Convert uptake to mmol N m-3 s-1 (NO3)
  uptake_NH4 = uptake_NH4*1d-9/Hz/secs_per_day !Convert uptake to mmol N m-3 s-1 (NH4)
  uptake_P = uptake_P*1d-9/Hz/secs_per_day !Convert uptake_P to mmol P m-3 s-1
  oxy_flux = oxy_flux*1d-9/Hz/secs_per_day !Convert oxy flux to mmol O2 m-3 s-1
  !# NOTE the different divisor to the uptake terms above: uptake accumulates
  !# dN_*ip_num, and dN_ is already a per-day RATE (pmol N cell-1 d-1), so
  !# /secs_per_day correctly converts it to per-second. Pmort accumulates
  !# ip_n*ip_num - an AMOUNT (pmol N), the whole nitrogen content of the cells
  !# that died during THIS call, which must all be returned to PON now. Dividing
  !# an amount by secs_per_day declares "this mass flows per day", so the ODE
  !# integrating it over the host step (dtdays*secs_per_day seconds) delivered
  !# only a fraction dtdays of it - 1/24 at dt=3600s - and silently destroyed the
  !# other 95.8% of every dead particle's N/P/C. Divide by the actual step length
  !# so the full amount reaches the water column.
  Pmort  =  Pmort*1d-9/Hz/(dtdays*secs_per_day) !Convert Pmort to mmol N m-3 s-1
  Pmort_P  =  Pmort_P*1d-9/Hz/(dtdays*secs_per_day) !Convert Pmort_P to mmol P m-3 s-1
  Pmort_C  =  Pmort_C*1d-9/Hz/(dtdays*secs_per_day) !Convert Pmort_C to mmol C m-3 s-1

  ! Now increment water column properties based on particle fluxes. The bulk Eulerian
  ! nutrient/OM pools are shared across all species, so _FLUX_VAR_ uses the total
  ! across groups (SUM); the _DIAG_VAR_ outputs stay per-group so each species'
  ! contribution can be examined separately.
  _FLUX_VAR_(data%id_nit)   = _FLUX_VAR_(data%id_nit) + (pp_ND + RES_NO3 - SUM(uptake_NO3))    ! mmol/m3/s CHECK ML I think this should not be multiplied by dtdays so stays in seconds; pp_ND and RES are zooplankton variables
  _FLUX_VAR_(data%id_amm)   = _FLUX_VAR_(data%id_amm) + (pp_ND + RES_NH4 - SUM(uptake_NH4))    ! mmol/m3/s CHECK ML I think this should not be multiplied by dtdays so stays in seconds; pp_ND and RES are zooplankton variables
  _FLUX_VAR_(data%id_frp)   = _FLUX_VAR_(data%id_frp) + (pp_PD + RES_P - SUM(uptake_P))    ! mmol/m3/s CHECK ML I think this should not be multiplied by dtdays so stays in seconds; pp_PD and RES_P are zooplankton variables
  _FLUX_VAR_(data%id_pon) = _FLUX_VAR_(data%id_pon)  + SUM(Pmort) + (pp_DZ - pp_ND) ! mmol/m3/s CHECK pp_DZ and pp_ND are zooplankton variables
  _FLUX_VAR_(data%id_pop) = _FLUX_VAR_(data%id_pop)  + SUM(Pmort_P) + (pp_DZP - pp_PD) ! mmol/m3/s CHECK pp_DZP and pp_PD are zooplankton variables
  _FLUX_VAR_(data%id_poc) = _FLUX_VAR_(data%id_poc)  + SUM(Pmort_C)  ! mmol/m3/s CHECK would need to add equivalent zoop variables here if added zoops in
  _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + SUM(oxy_flux)  ! mmol/m3/s CHECK would need to add equivalent zoop variables here if added zoops in

  DO ng=1,data%num_phytos
     _DIAG_VAR_(data%id_d_nit(ng)) = (pp_ND + RES_NO3 - uptake_NO3(ng))  * secs_per_day    ! mmol/m3/day
     _DIAG_VAR_(data%id_d_amm(ng)) = (pp_ND + RES_NH4 - uptake_NH4(ng))  * secs_per_day    ! mmol/m3/day
     _DIAG_VAR_(data%id_d_frp(ng)) = (pp_PD + RES_P - uptake_P(ng))  * secs_per_day    ! mmol/m3/day
     _DIAG_VAR_(data%id_d_pon(ng)) = Pmort(ng) + (pp_DZ - pp_ND) * secs_per_day ! mmol/m3/day
     _DIAG_VAR_(data%id_d_pop(ng)) = Pmort_P(ng) + (pp_DZP - pp_PD) * secs_per_day ! mmol/m3/day
     _DIAG_VAR_(data%id_d_poc(ng)) = Pmort_C(ng) * secs_per_day ! mmol/m3/day would need to add equivalent zoop variables here if added zoops in
     _DIAG_VAR_(data%id_d_oxy(ng)) = oxy_flux(ng) * secs_per_day ! mmol/m3/day would need to add equivalent zoop variables here if added zoops in
  ENDDO

! END OF CODE FROM PIBM BIOLOGY SUBROUTINE
!------------------------------------------------------

!------------------------------------------------------
! CODE FROM PIBM Par2PHY SUBROUTINE

! initialize diag vars at 0 (already allocated above; PHYC/PHYN/PHYP/CHL are per-group arrays)
PHYC = 0
PHYN = 0
PHYP = 0
CHL  = 0
DO ng=1,data%num_phytos
   _DIAG_VAR_(data%id_phyc(ng)) = 0.
   _DIAG_VAR_(data%id_phyn(ng)) = 0.
   _DIAG_VAR_(data%id_phyp(ng)) = 0.
   _DIAG_VAR_(data%id_chl(ng)) = 0.
   _DIAG_VAR_(data%id_mcdiv(ng)) = 0.
   _DIAG_VAR_(data%id_mtopt(ng)) = 0.
   _DIAG_VAR_(data%id_mlnalpha(ng)) = 0.
ENDDO

!loop through all super-individuals to divide, mutate, and convert into Eulerian concentrations
IF (N_ > 0) THEN
   ! DO i = 1, N_ ! SCB - previously, we looped through the first N_ indeces
   
   DO j = 1, N_ ! SCB - new: we loop through all active particles.
   i = index_(j)
   g = p(i)%ptm_istat(GRP)

   ! Cell division: if cellular carbon is above the division threshold, it divides
   IF (p(i)%ptm_state(data%ip_c) >= p(i)%ptm_state(data%ip_cdiv)) THEN  !Divide
      !print *, 'I divided!'
      _DIAG_VAR_(data%id_N_birth(g)) = _DIAG_VAR_(data%id_N_birth(g)) + 1.
      _DIAG_VAR_(data%id_cum_birth(g)) = _DIAG_VAR_(data%id_cum_birth(g)) + 1.
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
               _DIAG_VAR_(data%id_N_mutate(g)) = _DIAG_VAR_(data%id_N_mutate(g)) + 1.

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

   !Cumulate concentrations of phyto C, N, and Chl, and mean traits for superindividuals in the layer (per group)
   PHYC(g) = PHYC(g) + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c)
   PHYN(g) = PHYN(g) + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_n)
   PHYP(g) = PHYP(g) + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_p)
   CHL(g)  = CHL(g) + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_chl)
   _DIAG_VAR_(data%id_mcdiv(g))    = _DIAG_VAR_(data%id_mcdiv(g))    + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * log(p(i)%ptm_state(data%ip_cdiv))
   _DIAG_VAR_(data%id_mtopt(g))    = _DIAG_VAR_(data%id_mtopt(g))    + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * p(i)%ptm_state(data%ip_Topt)
   _DIAG_VAR_(data%id_mlnalpha(g)) = _DIAG_VAR_(data%id_mlnalpha(g)) + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * p(i)%ptm_state(data%ip_LnalphaChl)
   ENDDO
ENDIF !End of iterating over all super-individuals

!Calculate Eulerian concentrations of phyto C, N, and Chl, and mean traits and covariances for each layer, per group
DO ng=1,data%num_phytos
   _DIAG_VAR_(data%id_phyc(ng)) = PHYC(ng) * 1d-9 / Hz   !Convert Unit to mmol/m^3
   _DIAG_VAR_(data%id_phyn(ng)) = PHYN(ng) * 1d-9 / Hz   !Convert Unit to mmol/m^3
   _DIAG_VAR_(data%id_phyp(ng)) = PHYP(ng) * 1d-9 / Hz   !Convert Unit to mmol/m^3
   _DIAG_VAR_(data%id_chl(ng))  = CHL(ng) * 1d-9 / Hz   !Convert Unit to mg Chl/m^3
   if (PHYC(ng) > 0.) then   ! need to make sure we are preventing dividing by 0 (case when there are no particles in layer)
      _DIAG_VAR_(data%id_mcdiv(ng))    = _DIAG_VAR_(data%id_mcdiv(ng))    / PHYC(ng)
      _DIAG_VAR_(data%id_mtopt(ng))    = _DIAG_VAR_(data%id_mtopt(ng))    / PHYC(ng)
      _DIAG_VAR_(data%id_mlnalpha(ng)) = _DIAG_VAR_(data%id_mlnalpha(ng)) / PHYC(ng)
   endif

   !Compute trait covariances
   _DIAG_VAR_(data%id_vcdiv(ng))    = 0d0
   _DIAG_VAR_(data%id_vtopt(ng))    = 0d0
   _DIAG_VAR_(data%id_vlnalpha(ng)) = 0d0
   _DIAG_VAR_(data%id_cov_AL(ng))   = 0d0
   _DIAG_VAR_(data%id_cov_TA(ng))   = 0d0
   _DIAG_VAR_(data%id_cov_TL(ng))   = 0d0
ENDDO

IF (N_ > 0) THEN
   DO j = 1, N_
      i = index_(j)
      g = p(i)%ptm_istat(GRP)
      _DIAG_VAR_(data%id_vcdiv(g)) = _DIAG_VAR_(data%id_vcdiv(g)) + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * (log(p(i)%ptm_state(data%ip_cdiv)) - _DIAG_VAR_(data%id_mcdiv(g)))**2
      _DIAG_VAR_(data%id_vtopt(g)) = _DIAG_VAR_(data%id_vtopt(g)) + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * (p(i)%ptm_state(data%ip_Topt) - _DIAG_VAR_(data%id_mtopt(g)))**2
      _DIAG_VAR_(data%id_vlnalpha(g)) = _DIAG_VAR_(data%id_vlnalpha(g)) + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * (p(i)%ptm_state(data%ip_LnalphaChl) - _DIAG_VAR_(data%id_mlnalpha(g)))**2
      _DIAG_VAR_(data%id_cov_TL(g)) = _DIAG_VAR_(data%id_cov_TL(g)) + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * (log(p(i)%ptm_state(data%ip_cdiv)) - _DIAG_VAR_(data%id_mcdiv(g))) * (p(i)%ptm_state(data%ip_Topt) - _DIAG_VAR_(data%id_mtopt(g)))
      _DIAG_VAR_(data%id_cov_AL(g)) = _DIAG_VAR_(data%id_cov_AL(g)) + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * (log(p(i)%ptm_state(data%ip_cdiv)) - _DIAG_VAR_(data%id_mcdiv(g))) * (p(i)%ptm_state(data%ip_LnalphaChl) - _DIAG_VAR_(data%id_mlnalpha(g)))
      _DIAG_VAR_(data%id_cov_TA(g)) = _DIAG_VAR_(data%id_cov_TA(g)) + p(i)%ptm_state(data%ip_num) * p(i)%ptm_state(data%ip_c) * (p(i)%ptm_state(data%ip_LnalphaChl) - _DIAG_VAR_(data%id_mlnalpha(g))) * (p(i)%ptm_state(data%ip_Topt) - _DIAG_VAR_(data%id_mtopt(g)))
   ENDDO
ENDIF

DO ng=1,data%num_phytos
   IF (PHYC(ng) > 0d0) then ! need to make sure we are preventing dividing by 0 (case when there are no particles in layer)
      _DIAG_VAR_(data%id_vcdiv(ng))    = _DIAG_VAR_(data%id_vcdiv(ng))   / PHYC(ng)
      _DIAG_VAR_(data%id_vtopt(ng))    = _DIAG_VAR_(data%id_vtopt(ng))   / PHYC(ng)
      _DIAG_VAR_(data%id_vlnalpha(ng)) = _DIAG_VAR_(data%id_vlnalpha(ng))/ PHYC(ng)
      _DIAG_VAR_(data%id_cov_TL(ng))   = _DIAG_VAR_(data%id_cov_TL(ng))  / PHYC(ng)
      _DIAG_VAR_(data%id_cov_TA(ng))   = _DIAG_VAR_(data%id_cov_TA(ng))  / PHYC(ng)
      _DIAG_VAR_(data%id_cov_AL(ng))   = _DIAG_VAR_(data%id_cov_AL(ng))  / PHYC(ng)
   ENDIF
ENDDO

! deallocate things
if (allocated(index_))  deallocate(index_)
if (allocated(Pmatrix)) deallocate(Pmatrix)
if (allocated(BN))      deallocate(BN)
if (allocated(BP))      deallocate(BP)
if (allocated(BC))      deallocate(BC)
if (allocated(PHYC))    deallocate(PHYC)
if (allocated(PHYN))    deallocate(PHYN)
if (allocated(PHYP))    deallocate(PHYP)
if (allocated(CHL))     deallocate(CHL)
if (allocated(uptake))      deallocate(uptake)
if (allocated(uptake_NO3))  deallocate(uptake_NO3)
if (allocated(uptake_NH4))  deallocate(uptake_NH4)
if (allocated(uptake_P))    deallocate(uptake_P)
if (allocated(sv_dC))   deallocate(sv_dC)
if (allocated(sv_dN))   deallocate(sv_dN)
if (allocated(sv_dNO3)) deallocate(sv_dNO3)
if (allocated(sv_dNH4)) deallocate(sv_dNH4)
if (allocated(sv_dP))   deallocate(sv_dP)
if (allocated(sv_dChl)) deallocate(sv_dChl)
if (allocated(sv_ESD))  deallocate(sv_ESD)
if (allocated(Pmort))   deallocate(Pmort)
if (allocated(Pmort_P)) deallocate(Pmort_P)
if (allocated(Pmort_C)) deallocate(Pmort_C)
if (allocated(oxy_flux))      deallocate(oxy_flux)
if (allocated(N_by_group))    deallocate(N_by_group)
if (allocated(Abun_by_group)) deallocate(Abun_by_group)

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
      ! Retrieve this group's own (local) phytoplankton biomass in this layer
      phy = _DIAG_VAR_(data%id_phyc(phy_i))

      ! Self-shading with contribution from this phytoplankton group's concentration.
      extinction = extinction + (data%phytos(phy_i)%KePHY*phy)
   ENDDO
END SUBROUTINE aed_light_extinction_phyto_abm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_phyto_abm
