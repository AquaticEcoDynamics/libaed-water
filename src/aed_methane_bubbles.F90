!###############################################################################
!#                                                                             #
!# aed_methane_bubbles.F90                                                     #
!#                                                                             #
!#  This is an empty example source file to start from with new modules.       #
!#  it needs a bit of work to clean out some stuff but leave decent example    #
!#                                                                             #
!#  Further Developed by :                                                     #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2024 - The University of Western Australia                       #
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
!# Created Feb 2024                                                            #
!#                                                                             #
!###############################################################################
!#                                                                             #
! Purpose: this module governs the transportation of bubble (ebullition process) 
!          in the water, which is mainly from "Modeling bubbles and dissolved gases 
!          in the ocean" [Liang et al., 2011]. And some parameters are calculated
!          according to "Bubbles and the air-sea exchange of gases in near-saturation 
!          conditions" [Woolf, D. and Thorpe, S., 1991].
!          The bubble dynamics run in single bubble mode.
!
!---------------------------------------------------------------------------------
!###############################################################################
!
#include "aed.h"
!
!
MODULE aed_methane_bubbles
!-------------------------------------------------------------------------------
! aed_methane_bubbles --- methane bubbles model
!---------------------------------------------------------------------------------

  use iso_fortran_env, only: r8 => REAL64
! use shr_ctrl_mod, e8 => SHR_CTRL_E8, inft => INFINITESIMAL_E8
! use phy_utilities_mod
! use data_buffer_mod
! use shr_param_mod

!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util,ONLY : aed_int2str

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_bubbles_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_bubbles_data_t
      !# Variable identifiers
      INTEGER :: id_pres, id_temp, id_airpres, id_dens, id_dz, id_depth, id_depth_c, id_area
      INTEGER :: id_oxy, id_ch4, id_pco2, id_n2, id_ph
      INTEGER :: id_Fsed_oxy, id_Fsed_dic, id_Ebb_flux_ch4
      INTEGER :: NGAS
      INTEGER :: id_ch4_diss_target
      INTEGER :: id_btmb_ch4_3d 
      INTEGER :: id_dGama, id_sed_load
      INTEGER :: id_ice_b, id_ice_w


      INTEGER,ALLOCATABLE :: id_Vab(:), id_Cbg(:, :), id_bsolu(:), id_bDiff(:), id_schmidt(:), id_bNumb(:), id_bubble_fluxes(:), id_gas_ex(:)

      !# Model parameters
      AED_REAL, ALLOCATABLE :: m_Rb0(:)
      AED_REAL, ALLOCATABLE :: Ebb_flux_Z(:)
      CHARACTER(len=48), allocatable :: gasname(:)

     CONTAINS
         PROCEDURE :: define            => aed_define_bubbles
         PROCEDURE :: calculate         => aed_calculate_bubbles
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_bubbles
         PROCEDURE :: initialize_column => aed_initialize_column_bubbles
         PROCEDURE :: initialize        => aed_initialize_bubbles
         PROCEDURE :: initialize_benthic => aed_initialize_bubbles_benthic
         PROCEDURE :: calculate_column => aed_calculate_column_bubbles
   END TYPE

! MODULE GLOBALS
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

   CLASS (aed_bubbles_data_t),POINTER :: data
   TYPE (aed_column_t),POINTER :: column(:)

   ! Global constants
   AED_REAL :: inft = 1.d-30                                                                                !KK need to make these globals
   AED_REAL :: T0 = 273.15
   AED_REAL :: e8 = 1.d-8    
   AED_REAL :: Gconst = 9.8  !m/s2
   AED_REAL :: inf = 1.d+30  
   AED_REAL :: Pi = 3.14159265d+0
   AED_REAL :: R = 8.314462175d+0  !J/(mole*K)
   AED_REAL :: P0 = 98150.0
   AED_REAL :: Xo2 = 0.21
   AED_REAL :: MasO2 = 32
   !integer, parameter :: Wn2 = 1, Wo2 = 2, Wco2 = 3, Wch4 = 4
   !integer, parameter :: NGAS = 4
   INTEGER :: NRLAYER
   INTEGER  :: n_zones = 4
   AED_REAL :: m_bubbleGasCon, m_iceBubblePool


   !----------------------------------------------------------------------------
   
   TYPE LakeInfo
      integer  :: id                   ! lake id
      integer  :: itype                ! lake type identifier
      AED_REAL :: latitude             ! lake latitude
      AED_REAL :: longitude            ! lake longitude
      AED_REAL :: depth                ! lake depth (m)
      AED_REAL :: Asurf                ! lake surface area (m2)
      AED_REAL :: zalt                 ! lake altitude (km)
      AED_REAL :: kext                 ! chla and CDOM light attenuation (m-1)
      AED_REAL :: excice               ! land excessive ice fraction
      AED_REAL :: hsed                 ! sediment column thickness (m)
      integer  :: thrmkst              ! 1=thermokarst, 0=vice
      integer  :: margin               ! 1=margin, 2=slope, 0=center
      integer  :: hydroconn            ! 1=hydro connect, 0=vice
      integer  :: mixing
      integer  :: roun
   END TYPE

   TYPE(LakeInfo) :: lake_info

   integer, parameter :: Wn2 = 1, Wo2 = 2, Wco2 = 3, Wch4 = 4
   integer, parameter :: NGAS = 4

   !integer, parameter :: NLAKTYPE = 5

   !AED_REAL, parameter :: G = 9.8
   !AED_REAL, parameter :: INFINITESIMAL_E8 = 1.d-30
   !AED_REAL, parameter :: inft = 1.d-30
   !AED_REAL, parameter :: INFINITE_E8 = 1.d+30
   !AED_REAL, parameter :: inf = 1.d+30
   !AED_REAL, parameter :: Pi = 3.14159265d+0
   !AED_REAL, parameter :: R = 8.314462175d+0
   AED_REAL, parameter :: Roul = 1000.0, Roui = 931.0
   !AED_REAL, parameter :: T0 = 273.15
   !AED_REAL, parameter :: P0 = 98150.0
   !AED_REAL, parameter :: Xo2 = 0.21
   !AED_REAL, parameter :: MasO2 = 32
   !AED_REAL, parameter :: LKpH(NLAKTYPE) = (/6.8, 5.0, 7.0, 5.0, 6.5/)
   AED_REAL, parameter :: SOLO2(11) = (/14.6, 12.8, 11.3, 10.1, 9.1, 8.3, &
                           7.6, 7.0, 6.5, 6.0, 5.6/)

   !AED_REAL, allocatable :: m_waterTemp(:), m_sedTemp(:)
   !AED_REAL, allocatable :: m_dVsc(:)
   !AED_REAL, allocatable :: m_Zw(:)
   !AED_REAL, allocatable :: m_waterSubCon(:,:)

   !AED_REAL :: m_bubbleGasCon, m_iceBubblePool

   !integer :: NWLAYER, NSLAYER, NRLAYER
   !integer :: WATER_LAYER
   
   ! time step (s)
   AED_REAL :: tdelta

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed_define_bubbles(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with AED
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_bubbles_data_t),INTENT(inout) :: data

   INTEGER :: gas
   
!
!LOCALS
   CHARACTER(len=64) :: rnum
   INTEGER  :: r, g
   
   INTEGER           :: status

   INTEGER           :: nradius        = 1

   CHARACTER(len=64)  :: ch4_diss_target_variable=''
   CHARACTER(len=64)  :: ch4_btmb_target_variable=''
   
   NAMELIST /aed_methane_bubbles/nradius,ch4_diss_target_variable,ch4_btmb_target_variable

!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_bubbles configuration"

   !# Read the namelist
   read(namlst,nml=aed_methane_bubbles,iostat=status)
   IF (status /= 0) THEN
      print *,'Error reading namelist for &aed_methane_bubbles'
      STOP
   ENDIF


   NRLAYER = nradius                        !number of radius groups
   
   

       
   !# Register diagnostic variables
   ALLOCATE(data%id_Cbg(NGAS, NRLAYER+1))   ! bubble gas amount [mol/m3] 
   ALLOCATE(data%id_bsolu(NGAS))            ! gas solubility [mol/m3/Pa]
   ALLOCATE(data%id_bDiff(NGAS))               ! diffusivity [m2/s]
   ALLOCATE(data%id_schmidt(NGAS))             ! schmidt number [m2/s]
   ALLOCATE(data%id_bNumb(NRLAYER+1))          ! bubble amount 
   ALLOCATE(data%id_Vab(NRLAYER+1))         ! bouyant velocity  [m/s]
   ALLOCATE(data%id_bubble_fluxes(NGAS))    ! total flux of each gas across bubble sizes
   ALLOCATE(data%id_gas_ex(NGAS))           ! total gas exchange to a water layer for each gas

   
   
   ALLOCATE(data%gasname(NGAS))
   data%gasname = (/'Wn2 ', 'Wo2 ', 'Wco2', 'Wch4'/)

   data%id_dGama = aed_define_diag_variable('dGama_', 'N/m', 'surface tension')
   data%id_sed_load = aed_define_diag_variable('ch4_sed', 'mmol/d', "sediment methane")
   
   DO r=1,NRLAYER+1
      CALL aed_int2str(r, rnum)
      
      ! Bulk bubble variables, for each bubble size class
      data%id_Vab(r) = aed_define_diag_variable('Vab_'//rnum, 'units', 'name ')
      data%id_bNumb(r) = aed_define_diag_variable('bNumb_'//rnum, '', 'bubble amount')


      ! Gas specific allocations, for each bubble size class
      DO g=1,NGAS
         data%id_Cbg(g,r) = aed_define_diag_variable('Cbg_'//TRIM(data%gasname(g))//'_'//TRIM(rnum), 'units', 'name ')
      ENDDO

   ENDDO
   DO g=1,NGAS
      data%id_bsolu(g) = aed_define_diag_variable('bsolu_'//TRIM(data%gasname(g)), 'mol/m3/Pa', 'gas solubility')
      data%id_bDiff(g) = aed_define_diag_variable('bDiff_'//TRIM(data%gasname(g)), 'm2/s', 'gas diffusivity')
      data%id_schmidt(g) = aed_define_diag_variable('schmidt_'//TRIM(data%gasname(g)), 'm2/s', 'schmidt number')
      data%id_bubble_fluxes(g) = aed_define_sheet_diag_variable('bub_fluxes_'//TRIM(data%gasname(g)), 'units/m2', 'name ') ! same as bubble_fluxes(NGAS) update name later based on gas names KK: done bu had to shorten bubble_fluxes_ to bub_fluxes
      data%id_gas_ex(g) = aed_define_diag_variable('gas_exchange_'//TRIM(data%gasname(g)), 'unit', 'total gas exchange from bubbles') ! same as m_gasExchange(NGAS,WATER_LAYER+1)
   ENDDO


   !# Register external state variable dependencies
   data%id_oxy = aed_locate_variable('OXY_oxy')
   data%id_ch4 = aed_locate_variable('CH4_ch4') ! KK: This can be retired
   data%id_pco2 = aed_locate_variable('CAR_pCO2') ! KK change this to pCO2
   data%id_ph = aed_locate_variable('CAR_pH')

   data%id_btmb_ch4_3d = aed_locate_variable('CH4_ch4_ebb_dsfv')   ! mmol /m2 /s
   data%id_ch4_diss_target = aed_locate_variable(ch4_diss_target_variable)


   !# Register environmental dependencies
   data%id_pres = aed_locate_global('pressure')
   data%id_airpres = aed_locate_sheet_global('air_pres')
   data%id_temp = aed_locate_global('temperature')
   data%id_depth_c = aed_locate_sheet_global('col_depth')
   data%id_dens = aed_locate_global('density')
   data%id_dz   = aed_locate_global('layer_ht')
   data%id_area = aed_locate_global('layer_area')
   data%id_depth= aed_locate_global('depth')
   data%id_ice_b = aed_locate_sheet_global('delzBlueIce')
   data%id_ice_w = aed_locate_sheet_global('delzWhiteIce')
   

   ! BUBBLE JOBS
   allocate(data%m_Rb0(NRLAYER+1))
   !allocate(data%Ebb_flux_Z(4))           !allocate temporary Vbg
   call InitializeBubbleStateVariables()  !Setting m_Rb0, tdelta

CONTAINS

   subroutine InitializeBubbleStateVariables()
      implicit none
      AED_REAL :: minR, maxR, rratio, dr
      integer :: rr


      minR = 1.0_r8 !2.5_r8     ! mm  !MH > goto nml                            !KK bubbles measured at the surface of the studied lakes had diameters within 5-20 mm
      maxR = 5.0_r8 !10.0_r8    ! mm  !MH > goto nml
      dr = (maxR - minR) / dble(NRLAYER)                                !increment in radius between consecutive layers of bubbles. 
      data%m_Rb0(:) = (/(minR+dr*(rr-1), rr = 1, NRLAYER+1)/)           
      data%m_Rb0(:) = 1.0e-3_r8 * data%m_Rb0(:)  ! mm => m              !Convert mm to meters (1e-3 mm/mm)
      rratio = minR / 0.1_r8
      tdelta = 2.0_r8 / (0.1418*rratio**2 + 0.05579*rratio + 0.7794)
      !m_iceBubblePool = 0.0_r8
      !m_topbflux = 0.0_r8
      !m_upbflux = 0.0_r8
   end subroutine


END SUBROUTINE aed_define_bubbles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE aed_initialize_bubbles(data, column, layer_idx) 
   !-------------------------------------------------------------------------------
   !ARGUMENTS
      CLASS (aed_bubbles_data_t),INTENT(in) :: data
      TYPE (aed_column_t),INTENT(inout) :: column(:)
      INTEGER,INTENT(in) :: layer_idx

END SUBROUTINE aed_initialize_bubbles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_initialize_column_bubbles(data, column, layer_map) !MH Change layer_idx layer_map
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_bubbles_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_map(:)
!

END SUBROUTINE aed_initialize_column_bubbles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_initialize_bubbles_benthic(data,column,layer_idx)
!-------------------------------------------------------------------------------
   CLASS (aed_bubbles_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!-------------------------------------------------------------------------------



END SUBROUTINE aed_initialize_bubbles_benthic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_bubbles(data,column,layer_idx)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_bubbles_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL :: ch4_sed

!
!-------------------------------------------------------------------------------
!BEGIN

   _FLUX_VAR_(data%id_ch4) = _FLUX_VAR_(data%id_ch4) + _DIAG_VAR_(data%id_gas_ex(4))/secs_per_day
   

   

END SUBROUTINE aed_calculate_bubbles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE aed_calculate_column_bubbles(data,column,layer_map)
   !-------------------------------------------------------------------------------
   !ARGUMENTS
      CLASS (aed_bubbles_data_t),INTENT(in) :: data
      TYPE (aed_column_t),INTENT(inout) :: column(:)
      INTEGER,INTENT(in) :: layer_map(:)
   !
   !-------------------------------------------------------------------------------
   !LOCALS

      integer  :: rIndx, locIndx, top, bot, layer_idx, g, ii, k
      integer  :: icol, nscol
      AED_REAL :: temp, vsc, depth, dz, dens, airpres, ph
      AED_REAL :: pressure, gas_state_var, Hn2
      CHARACTER(len=64) :: gas, gasname
      AED_REAL :: rr, wb
      AED_REAL :: gama, pos
      AED_REAL :: ice_water_interface
      AED_REAL :: ratio(NGAS), ndm(NGAS), Vbg(NGAS), Cbg1(NRLAYER+1, size(layer_map)), Vbg0(NGAS, NRLAYER+1, size(layer_map))
      AED_REAL :: flux(NGAS), Cbg_sub(NRLAYER+1, size(layer_map)), gas_ex_up(size(layer_map)), area(size(layer_map)), bNumb(NRLAYER+1, size(layer_map))
      !new arrays
      AED_REAL :: Cbg(NGAS, NRLAYER+1, size(layer_map)), m_gasExchange(NGAS, size(layer_map)), m_Az(size(layer_map))
      AED_REAL :: m_btmbflux(NGAS, size(layer_map)), Vab(NRLAYER+1,size(layer_map))
      AED_REAL :: tmp1, da_sum, topAz, m_dZw, sumAz, m_topbflux
      AED_REAL :: vec, gas_exchange
      AED_REAL :: dr_tmp1, dr_tmp2, pr_tmp
      AED_REAL :: reynold(NGAS), k_ndm(NGAS), peclet(NGAS), nusselt(NGAS)
      AED_REAL :: exlayer
  !-----------------------------------------------------------------------------------
      gas_state_var = 0.

   !External variables
      airpres = _STATE_VAR_S_(data%id_airpres) * 100.0 !KK convert from hPa to Pa

   !Layer map
      bot = layer_map(SIZE(layer_map))
      top = layer_map(1) 
      ice_water_interface = 0.
      nscol = size(layer_map)

      m_btmbflux(:, :) = 0.0                          !KK set the bub flux to 0
      
      
   !print*, "This replicates GetBubbleReleaseMagnitude": Producing bubbles based on sediment fluxes                       !con(NGAS,NRLAYER+1)==Vbg(NGAS, NRLAYER+1), unit: umol/(m3*mm)
      
   ! LAYER LOOP NEEDED HERE
   do icol = 1, nscol, 1
      layer_idx = icol
      temp = _STATE_VAR_(data%id_temp) + 273.15                                        
      vsc = CalcDynamicViscosity(temp)                                                          
      gama = CalcSurfaceTension(temp)
      dens = _STATE_VAR_(data%id_dens)
      depth = _STATE_VAR_(data%id_depth)
      sumAz = _STATE_VAR_(data%id_area)
      pressure = airpres + dens*Gconst*depth                                                    ! KK pressure = m_surfData%pressure + Roul*G*(btmZw(indx)- m_surfData%dzsurf)                                       !pressure = m_surfData%pressure + Roul*G*lake_info%depth
      tmp1 = 0.0_r8
      ! Set CH4 bottom bubble flux for layer, using linked flux
      m_btmbflux(4, icol) = _DIAG_VAR_(data%id_btmb_ch4_3d)*1e-3     ! KK mmol to mol; set CH4 bub flux to diag
      m_btmbflux(1, icol) = 0.15 * m_btmbflux(4, icol)                                          !4.96*1e-8                                
      !print*, "diag1", m_btmbflux(4, icol)
      _DIAG_VAR_(data%id_sed_load) = m_btmbflux(4, icol) * secs_per_day * sumAz * 1000
      !print*, "sed_load", _DIAG_VAR_(data%id_sed_load)
      do rIndx = 1, NRLAYER+1, 1                                                                
         rr = data%m_Rb0(rIndx)
         wb = CalcBuoyantVelocity(rr, vsc/dens)
         Vab(rIndx,icol) = wb
         tmp1 = tmp1 + (pressure*rr + 2.0*gama) * wb                                            !KK [N/m]*[m/s]
         do g = 1, NGAS, 1
            Vbg0(g, rIndx, icol) = m_btmbflux(g, icol) * (pressure*rr + 2.0*gama)                  !KK [mol/m2/s]*[N/m]
         end do
         !print*, "Vbg0_init", Vbg0(4, rIndx, locIndx)
      end do
      Vbg0(:, :, icol) = Vbg0(:, :, icol) / (tmp1 + inft)
      !print*, "Atmp", Atmp
   end do

   do locIndx = bot, top, 1
      layer_idx = locIndx
      do rIndx = 1, NRLAYER+1, 1
         _DIAG_VAR_(data%id_Vab(rIndx)) = Vab(rIndx,locIndx)
      end do
   end do
   
   do icol = 1, nscol, 1
      layer_idx = icol
      temp = _STATE_VAR_(data%id_temp) + 273.15                                        
      vsc = CalcDynamicViscosity(temp)                                                          
      gama = CalcSurfaceTension(temp)
      dens = _STATE_VAR_(data%id_dens)
      depth = _STATE_VAR_(data%id_depth)
      do rIndx = 1, NRLAYER+1, 1
         rr = data%m_Rb0(rIndx)
         pressure = airpres + dens*Gconst*depth + 2.0*gama/rr
         do g = 1, NGAS, 1
            bNumb(rIndx,icol) = 0.75*R*temp/(Pi*pressure*rr**3.0)* &                     ! KK why is this not just 0.75, 0.75d-3
               sum(Vbg0(:,rIndx,icol))                                                      !KK this sums the gases for rIndx and icol
         end do
      end do  
   end do
   
   do locIndx = bot, top, 1
      layer_idx = locIndx
      do rIndx = 1, NRLAYER+1, 1
         _DIAG_VAR_(data%id_bNumb(rIndx)) = bNumb(rIndx,locIndx)
      end do
   end do
   
   

   !print*, "This replicates Bubble Module Setup"
   da_sum = 0.0
   topAz = 0.0
   !print*, "Here"
   do locIndx = bot, top, 1                                                                                                                                  
      layer_idx = locIndx
      temp = MAX(_STATE_VAR_(data%id_temp) + 273.15, T0) 
      pH = _STATE_VAR_(data%id_ph)
      sumAz = _STATE_VAR_(data%id_area) 
      da_sum = da_sum + sumAz
      m_Az(locIndx) = da_sum
      do g = 1, NGAS, 1
         _DIAG_VAR_(data%id_bsolu(g)) = CalcHenrySolubility(g,temp,pH)
         _DIAG_VAR_(data%id_bDiff(g)) = CalcGasDiffusivityInWater(g,temp)
         _DIAG_VAR_(data%id_schmidt(g)) = CalcSchmidtNumber(g,temp)
      end do
   end do
   
   topAz = da_sum
   
   !print*, "topAz", topAz

!print*, "This replicates bubble dynamics"   
   
   
   m_gasExchange(:, :) = 0.0                                                                
   m_topbflux = 0.0
   da_sum = 0.0
   sumAz = 0.0
   !print*, "here2"
   do icol = 1, nscol, 1  
      layer_idx = icol  
      !print*, "layer_idx", layer_idx
      sumAz = _STATE_VAR_(data%id_area)   
      m_dZw = _STATE_VAR_(data%id_dz)                                                      
      if (m_btmbflux(4, icol) < 1.d-12) then 
          m_gasExchange(4,icol) = m_gasExchange(4,icol) + m_btmbflux(4, icol)*sumAz/m_Az(icol)/m_dZw
          cycle
      end if                         
      m_topbflux = m_topbflux + m_btmbflux(4, icol)*sumAz/topAz                       
      Cbg(:,:,:) = 0.0 
      if (icol > 1) then                                            
         Cbg(:,:,icol-1) = Vbg0(:,:,icol)  
      end if  
                                                      
                                                          
      do rIndx = 1, NRLAYER+1, 1  
         !print*, "rIndx", rIndx                 
         locIndx = icol   
         layer_idx = locIndx 
         !print*, "layer_idx1", layer_idx  
         depth     = - _STATE_VAR_(data%id_depth) 
         m_dZw = _STATE_VAR_(data%id_dz)                                                           ! KK start from the bottom                                                                                          
         pos = depth - m_dZw 
         !print*, "pos", pos
         rr =  data%m_Rb0(rIndx)
         Vbg = Vbg0(:,rIndx,icol)
         !print*, "Vbg", Vbg
         do while (pos<ice_water_interface .and. rr>e8) 
            if ( pos > depth) then
               !print*, "in the loop"
               !print*, "depth", depth
               !print*, "locIndx", locIndx
               Cbg(:, rIndx, locIndx) = Vbg
               !print*, "Cbg", Cbg(:, rIndx, locIndx)
               locIndx = locIndx + 1
               layer_idx = locIndx
               depth     = - _STATE_VAR_(data%id_depth)
               !print*, "depth1", depth
            end if
            !print*, "Here I am"
            temp = MAX(_STATE_VAR_(data%id_temp) + 273.15, T0)
            dens = _STATE_VAR_(data%id_dens)
            vsc = CalcDynamicViscosity(temp)
            wb = CalcBuoyantVelocity(rr, vsc/dens)
            !print*, "wb", wb
            gama = CalcSurfaceTension(temp)                                            
            pressure = airpres - dens*Gconst*pos + 2.0*gama/rr                                     !KK with the pressure their position is references from the top so it's reveresed
            ! gas exchange rate  
            do g=1, NGAS, 1
               peclet(g) = rr * wb / _DIAG_VAR_(data%id_bDiff(g))  
               reynold(g) = peclet(g) / _DIAG_VAR_(data%id_schmidt(g))
               if (reynold(g)<=1.0) then
                  nusselt(g) = sqrt(2.0*Pi*peclet(g)/3.0)
               else
                  nusselt(g) = 2.0*sqrt(peclet(g)/Pi)                                                       !0.45*(reynold(g)**(1.0/6.0))*(peclet(g)**(1.0/3.0))
               endif
               k_ndm(g) = -4.0 * Pi * rr * _DIAG_VAR_(data%id_bDiff(g)) * nusselt(g)       
               call CalcGasRatioInSingleBubble(Vbg(:), ratio(:))   
               gas = TRIM(data%gasname(g))
               SELECT CASE (gas)
                  CASE ('Wch4')
                     gas_state_var = _STATE_VAR_(data%id_ch4)/1000.0                                  !unit conversion mmol/m3 to mol/m3                              
                     !print*, "ch4", gas_state_var
                  CASE ('Wo2')
                     gas_state_var = _STATE_VAR_(data%id_oxy)/1000.0
                     !print*, "oxy", gas_state_var
                  CASE ('Wco2')
                     gas_state_var = (_STATE_VAR_(data%id_pco2)*101325.0)/(R*temp)                    !converting pCO2 (atm) to mol/m3
                     !print*, "co2", gas_state_var
                  CASE ('Wn2')
                     Hn2 = 6.4e-6 * exp(1300.0*((1/temp)-(1/298.0)))                                  !mole/m3/Pa
                     !print*, "Hn2", Hn2
                     gas_state_var = Hn2 * 0.78 * P0
                     !gas_state_var = Hn2 * 40530.0                                                    !partial pressure of N2 is 40% of atmospheric pressure (Peltola et al. 2017)
                     !print*, "n2_conc", gas_state_var
               END SELECT 
               ndm(g) = k_ndm(g) * (_DIAG_VAR_(data%id_bsolu(g))*pressure*ratio(g) - &
                  gas_state_var)  
               Vbg(g) = Vbg(g) + bNumb(rIndx, icol) * ndm(g) * tdelta  
            end do     
            where (Vbg<0._r8) Vbg = 0._r8   
            pr_tmp = 3.0 * airpres - 3.0 * dens*Gconst*pos + 4.0*gama/rr   
            dr_tmp1 = 0.75*R*temp/(Pi*rr**2.0)/pr_tmp                                              !KK why is this not just 0.75   , 0.75d-3
            dr_tmp2 = rr*dens*Gconst*wb/pr_tmp     
            rr = rr + (dr_tmp1*sum(ndm)+dr_tmp2)*tdelta
            ! new position
            pos = pos + wb*tdelta
            !print*, "pos1", pos
         end do
         if (rr>e8 .and. locIndx==top) then                                                           ! KK handling the top layer: I am assigning last Vbg to Cbg in the top layer
            layer_idx = locIndx
            Cbg(:, rIndx, locIndx) = Vbg
            !print*, "Cbgtop", Cbg(:, rIndx, locIndx) 
         endif
      end do


      do ii = icol, top, 1
         layer_idx = ii
         !print*, "layer_idx5", layer_idx
         !print*, "icol", icol
         m_dZw = _STATE_VAR_(data%id_dz)
         ! size-integrated bubble gas and gas exchange
         do rIndx = 1, NRLAYER+1, 1
            !print*, "rIndx", rIndx
            !print*, "ii", ii
            !print*, "icol", icol
            if (ii == 1) then
               exlayer = (Vbg0(4, rIndx, ii) - Cbg(4,rIndx,ii))*Vab(rIndx,icol)
            else
               exlayer = (Cbg(4,rIndx,ii-1) - Cbg(4,rIndx,ii))*Vab(rIndx,icol)
            end if
            if (exlayer < 0.0) then
               !print*, "exlayer", exlayer
            end if
            if (ii==icol) then
               sumAz = _STATE_VAR_(data%id_area)
               !print*, "PathA"
               !print*, "sumAz", sumAz
               m_gasExchange(4,ii) = m_gasExchange(4,ii) + exlayer/m_dZw
            else
               m_gasExchange(4,ii) = m_gasExchange(4,ii) + exlayer*sumAz/ &
                  m_Az(ii)/m_dZw
                  !print*, "sumAz", sumAz
                  !print*, "m_Az", m_Az(ii)
                  !print*, "m_dZw", m_dZw
                  !print*, "PathB", sumAz/m_Az(ii)/m_dZw
                  !print*, "sumAzB", sumAz
            end if
            !print*, "gasex", m_gasExchange(4,ii)
            _DIAG_VAR_(data%id_Cbg(4,rIndx)) = Cbg(4,rIndx,ii)
         end do
         !print*, "gas_exchangemid", m_gasExchange(4,ii)
      end do


   end do

   ! update bubble flux
   !print*, "mtopbfluxpre", m_topbflux
   !print*, "releasenogasex", m_topbflux * 1000 * secs_per_day * topAz
   do locIndx = bot, top, 1
      layer_idx = locIndx
      !print*, "layer_idx", layer_idx
      m_dZw = _STATE_VAR_(data%id_dz)
      !print*, "m_dZw", m_dZw
      m_topbflux = m_topbflux - m_gasExchange(4,locIndx)*m_dZw* &
         m_Az(locIndx)/topAz
      !print*, "m_topbflux", m_topbflux
      !print*, "topAz", topAz
      _DIAG_VAR_(data%id_gas_ex(4)) = m_gasExchange(4,locIndx) * 1000 * secs_per_day
   !print*, "gas_ex_final", m_gasExchange(4,locIndx)
   !print*, "gas_ex_diag", _DIAG_VAR_(data%id_gas_ex(4))
   end do

   if (_STATE_VAR_S_(data%id_ice_b) + _STATE_VAR_S_(data%id_ice_w) > 1.0e-4) then
      ! Redirect blocked bubble flux to top water layer as dissolved CH4
      layer_idx = top
      m_dZw = _STATE_VAR_(data%id_dz)
      _DIAG_VAR_(data%id_gas_ex(4)) = _DIAG_VAR_(data%id_gas_ex(4)) + &
                                       m_topbflux * 1000 * secs_per_day / m_dZw
      m_topbflux = 0.0_r8
   end if


   _DIAG_VAR_S_(data%id_bubble_fluxes(4)) = m_topbflux * 1000 * secs_per_day
   !print*, "top_ch4", _DIAG_VAR_S_(data%id_bubble_fluxes(4))
   
   CONTAINS      
      real function CalcDynamicViscosity(temp)  !KK Does this need be calced for each water layer in initialize?
      implicit none
      AED_REAL, intent(in) :: temp    ! units: K
      AED_REAL :: dVsc  ! units: kg/s/m
  
      !if (temp<T0-e8) THEN
      !  CalcDynamicViscosity = inf 
      !else
        CalcDynamicViscosity = 2.414d-5 * (10.0**(247.8/(temp-140.0)))
      !endif
      end function

      


     




   
 END SUBROUTINE aed_calculate_column_bubbles
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_bubbles(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED test.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_bubbles_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL :: ebb_sed, sed_flux(16)
   INTEGER :: i
!
!LOCALS
   ! Temporary variables

!-------------------------------------------------------------------------------
!BEGIN


END SUBROUTINE aed_calculate_benthic_bubbles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!Physics functions
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate water surface tension from Wikipedia "Surface_tension"
   !  Surface tension is expressed as a surface energy per unit area or the work
   !  required to increase unit surface area at constant temperature.
   !
   !------------------------------------------------------------------------------
   function CalcSurfaceTension(temp)
      implicit none
      AED_REAL, intent(in) :: temp
      AED_REAL :: CalcSurfaceTension               ! units: N/m
      AED_REAL, parameter :: r0 = 75.64d-3, r25 = 71.97d-3, r50 = 67.91d-3
      AED_REAL :: par
   
      if (temp<T0) then
         !CalcSurfaceTension = 0.0_r8
         CalcSurfaceTension = r0
      else if (temp>=T0 .and. temp<T0+25) then
         par = (temp-T0)/25.0
         CalcSurfaceTension = r0*(1-par) + r25*par
      else if (temp>=T0+25 .and. temp<T0+50) then
         par = (temp-T0-25.0)/25.0
         CalcSurfaceTension = r25*(1-par) + r50*par
      else
         CalcSurfaceTension = r50
      end if
      return

      end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate buoyant rising speed of bubble from "Bubbles and the
   !          air-sea exchange of gases in near-saturation conditions" [Woolf,
   !          D. and Thorpe, S., 1991].
   !
   !------------------------------------------------------------------------------
   function CalcBuoyantVelocity(radius, vsc)
      implicit none
      AED_REAL, intent(in) :: radius      ! units: m
      AED_REAL, intent(in) :: vsc         ! kinematic viscosity (m2/s)
      AED_REAL :: CalcBuoyantVelocity     ! units: m/s               
      AED_REAL :: rr, xx, yy

      if (vsc>1.0d10) then
         CalcBuoyantVelocity = 1.e-30_r8
      else 
         xx = Gconst * (radius**3.0_r8) / (vsc**2.0_r8)
         yy = 10.82_r8 / xx
         CalcBuoyantVelocity = (2.0_r8*(radius**2.0_r8)*Gconst/9.0_r8/vsc) * &
            ((yy**2.0_r8+2.0_r8*yy)**0.5_r8-yy)
      end if
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate gas solubility (henry law constant), which are from
   !          "Compilation of Henry law constants for inorganic and organic
   !          species of potential importance in environmental chemistry (version 3)"
   !          [Sander, 1999]
   !   For Bunsen coefficient, C[gas] = Bunsen * (P[gas] / (R * T)) * water_content
   !         (if it is pure water, water_content = 1)
   !
   !------------------------------------------------------------------------------
   function CalcHenrySolubility(gas, temp, pH)
      implicit none
      integer, intent(in) :: gas
      AED_REAL, intent(in) :: temp        ! units: K
      AED_REAL, intent(in) :: pH          ! units: n/a
      AED_REAL :: CalcHenrySolubility     ! units: umol/(m3*Pa) (M = mole/L)
      AED_REAL :: hi, kc1, kc2, par
      integer :: indx
      integer :: Wn2 = 1, Wo2 = 2, Wco2 = 3, Wch4 = 4


      if (gas==Wn2) then
         CalcHenrySolubility = 6.1d-6*exp(1300*(1/temp-1/298.0))
      else if (gas==Wo2) then
         if (temp>=T0 .and. temp<=T0+50) then
            indx = int((temp-T0)/5) + 1
            indx = min(indx, 10)
            par = (temp - T0 - 5*indx + 5) / 5.0
            CalcHenrySolubility = (SOLO2(indx+1)*par + SOLO2(indx)*(1-par)) &
                                    / MasO2 / P0 / Xo2 
         else
            CalcHenrySolubility = 1.3d-5*exp(1500*(1/temp-1/298.0))
         end if
      else if (gas==Wco2) then
         CalcHenrySolubility = 3.4d-4*exp(2400*(1/temp-1/298.0))
         hi = 10**(-pH)    ! Concentration of hydrogen ion
         ! rate constant of dissolved CO2 for first and second dissolution
         kc1 = 4.3d-7*exp(-921.4*(1/temp-1/298.0))
         kc2 = 4.7d-11*exp(-1787.4*(1/temp-1/298.0))
         CalcHenrySolubility = CalcHenrySolubility*(1.0+kc1/hi+kc1*kc2/hi**2)
      else if (gas==Wch4) then
         CalcHenrySolubility = 1.3d-5*exp(1700*(1/temp-1/298.0))
      end if
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate gas diffusivity in water and air, from Appendix A of J.
   !          Tang's thesis. In D. Woolf's "Bubbles and the air-sea exchange of
   !          gases in near saturation conditions", the diffusivity have different
   !          values for those gases (Table 1):
   !          Dn2 = 1.8d-9
   !          Do2 = 1.7d-9
   !          Dco2 = 1.3d-9
   !          Dch4 = 0.5d-9
   !          For ice: "Change in CO2 concentration and O2/N2 ratio in ice
   !          cores due to molecular diffusion" (Bereiter et al., 2009)
   !
   !------------------------------------------------------------------------------
   function CalcGasDiffusivityInWater(gas, temp)
      implicit none
      integer, intent(in) :: gas
      AED_REAL, intent(in) :: temp                 ! units: K
      AED_REAL :: CalcGasDiffusivityInWater        ! units: m2/s

      if (gas==Wn2) then
         CalcGasDiffusivityInWater = 2.57d-9*(temp/273.0)
      else if (gas==Wo2) then
         CalcGasDiffusivityInWater = 2.4d-9*(temp/298.0)
      else if (gas==Wco2) then
         CalcGasDiffusivityInWater = 1.81d-6*exp(-2032.6/temp)
      else if (gas==Wch4) then
         CalcGasDiffusivityInWater = 1.5d-9*(temp/298.0)
      end if
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate Schmidt number for N2, O2, CO2 and CH4
   !          Reference: "Relationship between wind-speed and gas-exchange over
   !          the ocean", R. Wanninkhof, 1992 (N2 and O2) and "Measurement of the
   !          diffusion coefficients of sparingly soluble gases in water", B.
   !          Jahne, et al., 1987 (CO2 and CH4).
   !
   !------------------------------------------------------------------------------
   function CalcSchmidtNumber(gas, temp)
      implicit none
      integer, intent(in) :: gas
      AED_REAL, intent(in) :: temp     ! units: K
      AED_REAL :: CalcSchmidtNumber
      AED_REAL :: T

      T = min(temp-T0,30.0)     ! to celcius
      if (T<0) then
         CalcSchmidtNumber = inf
      else
         if (gas==Wn2) then
            CalcSchmidtNumber = 1970.7-131.45*T+4.139*T**2-0.052106*T**3
         else if (gas==Wo2) then
            CalcSchmidtNumber = 1800.6-120.1*T+3.7818*T**2-0.047608*T**3
         else if (gas==Wco2) then
            CalcSchmidtNumber = 1911-113.7*T+2.967*T**2-0.02943*T**3
         else if (gas==Wch4) then
            CalcSchmidtNumber = 1898-110.1*T+2.834*T**2-0.02791*T**3
         end if
      end if
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate Nusselt number defined as the ratio between the total gas
   !          flux and the molecular diffusivity flux across the bubble surface,
   !          which is from "Bubbles and the air-sea exchange of gases in near-
   !          saturation conditions" [Woolf D. and Thorpe S., 1991].
   !
   !------------------------------------------------------------------------------
   function CalcPecletNumber(gas, radius, temp, vsc)
      implicit none
      integer, intent(in) :: gas
      real(r8), intent(in) :: radius         ! units: m
      real(r8), intent(in) :: temp           ! units: K
      real(r8), intent(in) :: vsc            ! units: m2/s
      real(r8) :: CalcPecletNumber
      real(r8) :: wb, diffusivity

      wb = CalcBuoyantVelocity(radius,vsc)
      diffusivity = CalcGasDiffusivityInWater(gas,temp)
      CalcPecletNumber = radius*wb/(diffusivity+inft)
      return
   end function

   function CalcReynoldNumber(gas, radius, temp, vsc)
      implicit none
      integer, intent(in) :: gas
      real(r8), intent(in) :: radius         ! units: m
      real(r8), intent(in) :: temp           ! units: K
      real(r8), intent(in) :: vsc            ! units: m2/s
      real(r8) :: CalcReynoldNumber
      real(r8) :: Pe, Schmidt                ! Pe = Re*Sc

      Pe = CalcPecletNumber(gas,radius,temp,vsc)
      Schmidt = CalcSchmidtNumber(gas,temp)
      CalcReynoldNumber = Pe/Schmidt
      return
   end function
   

   function CalcNusseltNumber(gas, radius, temp, vsc)
      implicit none
      integer, intent(in) :: gas
      real(r8), intent(in) :: radius         ! units: m
      real(r8), intent(in) :: temp           ! units: K
      real(r8), intent(in) :: vsc            ! units: ms/s
      real(r8) :: CalcNusseltNumber
      real(r8) :: PeN, ReN

      ReN = CalcReynoldNumber(gas,radius,temp,vsc)
      PeN = CalcPecletNumber(gas,radius,temp,vsc)
      if (ReN<=1.0) then
         CalcNusseltNumber = sqrt(2*Pi*PeN/3)
      else
         CalcNusseltNumber = 2.0*sqrt(PeN/Pi)
         !CalcNusseltNumber = 0.45*(ReN**(1.0/6.0))*(PeN**(1.0/3.0))
      end if
      return
   end function

   subroutine CalcGasRatioInSingleBubble(bubbleGasCon, ratio)
      implicit none
      AED_REAL, intent(in) :: bubbleGasCon(NGAS)
      AED_REAL, intent(out) :: ratio(NGAS)
      AED_REAL :: mgas_tot

      mgas_tot = abs(sum(bubbleGasCon))
      if (mgas_tot<1.d-12) then
         ratio = (/0.0, 0.0, 0.0, 1.0/)
      else
         ratio = bubbleGasCon / mgas_tot
      end if
   end subroutine

   function GetBubbleAmount(bubbleGasCon, radius, temp, pressure)
      implicit none
      AED_REAL, intent(in) :: bubbleGasCon(NGAS)
      AED_REAL, intent(in) :: radius
      AED_REAL, intent(in) :: temp
      AED_REAL, intent(in) :: pressure
      AED_REAL :: GetBubbleAmount   ! unit: bubble/(m3*mm)
      AED_REAL :: gases

      gases = sum(bubbleGasCon) 
      GetBubbleAmount = 1.0e-6_r8*gases*0.75*R*temp/ &
                  (Pi*pressure*(radius**3))
      return
   end function

   !function CalcDynamicViscosity(temp)  !KK Does this need be calced for each water layer in initialize?
      !implicit none
      !AED_REAL, intent(in) :: temp    ! units: K
      !AED_REAL :: CalcDynamicViscosity  ! units: kg/s/m
      !AED_REAL :: e8, T0, inf
      
      
      !e8 = 1.d-8
      !T0 = 273.15
      !inf = 1.d+30
 
      !if (temp<T0-e8) THEN
        !CalcDynamicViscosity = inf 
      !else
        !CalcDynamicViscosity = 2.414d-5 * (10.0**(247.8/(temp-140.0)))
      !endif
    !end function

END MODULE aed_methane_bubbles
