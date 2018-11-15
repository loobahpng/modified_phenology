!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_soil

  ! 
  ! Description: 
  !   Variables that hold soil parameters for JSBACH
  ! 
  ! Current Code Owner: jsbach_admin
  ! 
  ! History: 
  !  
  ! Version   Date        Comment 
  ! -------   ----        ------- 
  ! 0.1       2001/07/01  Original code. Reiner Schnur
  ! 
  ! Code Description: 
  !   Language:           Fortran 90. 
  !   Software Standards: "European Standards for Writing and  
  !     Documenting Exchangeable Fortran 90 Code". 
  !
  ! ******* Vs. 1.1 - June 2007 - Stefan Hagemann
  ! *** Implementation of 5 layer soil scheme
  ! ***   Note that the former soil moisture bucket 'moisture' is now defined as
  ! ***   root zone soil moisture. For simplicity of the implementation, the
  ! ***   array definitions of 'moisture' were kept as they are although a
  ! ***   multi-layer definition makes no sense anymore. Later nsoil=1 dimension may be
  ! ***   eliminated and ntsoil (=5 currently) may be set to nsoil.
  ! ***   Currently there is an unused array called root_moisture which should be removed later.
  ! *** Some unused arrays are eliminated.
  !
  ! ******* Vs. 1.2 - December 2009 - Stefan Hagemann
  ! *** New arrays implemented to allow changes in the percolation formulation of the 
  ! *** Final version of 5 layer scheme
  !
  ! ******* Vs. 1.3 - January 2010 - Stefan Hagemann
  ! *** Diagnostic subroutine for global water balance at each grid box.
  ! *** This is currently incomplete only for glaciated gridboxes. 
  !
  ! ******* Vs. 1.4 - March 2010 - Stefan Hagemann
  ! *** Adding a few output variables, such as Bare Soil and Skin reservoir evap.
  ! 
  ! ******* Vs. 1.5 - June 2012 - Stefan Hagemann
  ! *** Introducing a separate bare soil water storage bare_soil_moisture
  ! *** that is used to calculate Bare Soil evaporation
  ! 
  ! ******* Vs. 1.6 - May 2015 - Stefan Hagemann
  ! *** Bugfix for water balance diagnostic call for permafrost version
  ! 
  ! Modules used: 
  ! 
  USE mo_netCDF,           ONLY: FILE_INFO, BELOWSUR, SOILLEV, SNOWLEV
  USE mo_jsbach,           ONLY: debug
  USE mo_jsbach_grid,      ONLY: kstart, kend, nidx
  USE mo_linked_list,      ONLY: t_stream
  USE mo_mpi,              ONLY: p_parallel, p_parallel_io, p_bcast, p_io
  USE mo_io_units,         ONLY: nout
  USE mo_kind,             ONLY: dp 
  USE mo_exception,        ONLY: finish, message, message_text, em_warn
  USE mo_soil_temperature, ONLY: update_soil_temperature, update_soil_permafrost, hcond_none, hcond_johansen, get_liquid_max
  USE mo_soil_hydrology,   ONLY: update_surface_hydrology, update_soil_hydrology

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: soil_type, soil_param_type, init_soil, update_soil, soil_diagnostics, get_soil_diag

  ! Global (i.e. public) Declarations: 
  ! Global Type Definitions: 
  !! Structure to hold global soil parameters for each land cell and soil layer.
  TYPE soil_param_type
   !
   ! Variables describing soil characteristics and being read from initialization files
   !
   ! (nland)
   REAL(dp), POINTER, DIMENSION(:) :: &
        Ds,                           &           !! Fraction of maximum subsurface flow rate
        Dsmax,                        &           !! Maximum subsurface flow rate [mm/day]
        Ws,                           &           !! Fraction of maximum soil moisture where non-linear baseflow occurs
        B_infilt,                     &           !! Infiltration parameter
        AvgTemp,                      &           !! Average soil temperature used as the bottom
                                                  !!    boundary for soil heat flux solutions
        Roughness ,                   &           !! Surface roughness of bare soil
        VolHeatCapacity,              &           !! Volumetric heat capacity of soil [J/m**3/K]
        ThermalDiffusivity                        !! Thermal diffusivity of soil [m**2/s]
   REAL(dp), POINTER    :: InitMoisture(:)        !! Initial root zone soil moisture [m]
   REAL(dp), POINTER    :: MaxMoisture(:)         !! Maximum root zone moisture content (field capacity) [m]

   ! (nland), multi-layer soil scheme
   REAL(dp), POINTER    :: hyd_cond_sat(:)        !! Saturated hydraulic conductivity [m/s]
   REAL(dp), POINTER    :: PotMoisture(:)         !! Potential Moisture [m]
   REAL(dp), POINTER    :: bclapp(:)              !! exponent b in Clapp and Hornberger 
   REAL(dp), POINTER    :: RootDepth(:)           !! Rooting Depth in [m] 
   REAL(dp), POINTER    :: SoilDepth(:)           !! Soil Depth above bedrock [m]
   REAL(dp), POINTER    :: SoilPorosity(:)        !! Volumetric soil porosity [m/m]   
   REAL(dp), POINTER    :: FieldCapacity(:)       !! Volumetric soil field capacity [m/m]   
   REAL(dp), POINTER    :: WiltingPoint(:)        !! Volumetric soil wilting point [m/m]   
   REAL(dp), POINTER    :: PoreSizeIndex(:)       !! Soil Pore size distribution index    
   REAL(dp), POINTER    :: HeatConductivity(:)    !! Heat conductivity of dry soil [J/(m*s*K)]
   ! (nland,nsoil) - multi-layer soil scheme
   REAL(dp), POINTER    :: InitLayerMoisture(:,:) !! Initial soil moisture of layers [m]

  END TYPE soil_param_type
  !
  !----------------------------------------------------------------------------------------------------------
  ! Soil state variables
  TYPE soil_type

     INTEGER          :: nsoil                 !! Number of soil layers
     INTEGER          :: ntsoil                !! Number of soil layers for soil temperature calculation
     INTEGER          :: ntiles                !! Number of soil tiles
     INTEGER          :: nsnow                 !! Number of snow layers
     REAL(dp), ALLOCATABLE  :: cdel(:)         !! soil layer depth [m]
     REAL(dp), ALLOCATABLE  :: clev(:)         !! depth of soil layer lower boundaries [m]
     REAL(dp), ALLOCATABLE  :: cmid(:)         !! depth of soil layer mid points [m]
     REAL(dp), POINTER, DIMENSION(:)   ::   &  !! (nland)
       csat,                                &  !!
       cair,                                &  !!
       relative_humidity_air,               &  !! Relative humidity of the air in the lowest layer of the atmosphere
       z0m,                                 &  !! roughness length momentum (needed for offline simulations)
       z0h,                                 &  !! roughness length heat and moisture (needed for offline simulations)
       zril_old,                            &  !! Richardson number of previous time step (needed for offline simulations)
       csat_transpiration,                  &  !!
       snow_age,                            &  !! non-dimensional age of snow
       reduced_evap_acc,                    &  !! Potential for reducing evaporation (accumulated) - multi-layer scheme only
       bare_soil_evap_acc,                  &  !! Bare soil evaporation (accumulated) - diagnostic
       snow_evap_acc,                       &  !! Snow evaporation (accumulated) - diagnostic
       skin_reservoir_evap_acc,             &  !! Skin reservoir evaporation (accumulated) - diagnostic
       evap_deficit_acc,                    &  !! Evaporation deficit flux due to inconsistent treatment of snow and skin evaporation
       precip_acc,                          &  !! Precipitation received by JSBACH (accumulated) - diagnostic
       water_balance,                       &  !! Soil water balance = Fluxes - changes in surface water storage (diagnostic)
       storage_pre,                         &  !! Water storage of previous timestep (needed for water balance)
       bare_soil_moisture                      !! Soil Moisture in bare soil part of gridbox

     REAL(dp), POINTER, DIMENSION(:,:) ::   &  !! (nland, ntiles)
       skin_reservoir,                      &  !! Water content of skin reservoir
       albedo,                              &  !! background albedo (ECHAM albedo scheme) or albedo_soil_vis (JSBACH albedo scheme)
       albedo_vegetation_vis,               &  !! albedo of vegetation in the visible range
       albedo_vegetation_nir,               &  !! albedo of vegetation in the NIR range
       albedo_soil_vis,                     &  !! albedo of soil in the visible range
       albedo_soil_nir,                     &  !! albedo of soil in the NIR range
       snow,                                &  !! Snow on ground [m water equivalent]
       snow_fract,                          &  !! Fraction of snow covered ground
       snow_accum_acc,                      &  !! Snow accumulation (budget) at non-glacier points [kg/(m^2 s)] (accumulated)
       snow_melt,                           &  !! Snow melt [m]
       snow_melt_acc,                       &  !! Snow melt (accumulated) [kg/(m^2 s)]
       runoff_acc,                          &  !! Total runoff (surface runoff + drainage) at non-glacier points (accumulated)
       drainage_acc,                        &  !! Drainage at non-glacier points (accumulated) [kg/(m^2 s)]
       glacier_depth,                       &  !! Glacier depth (including snow) [m water equivalent]
       glacier_precip_minus_evap_acc,       &  !! Precipitation minus evaporation for glaciers (accumulated)
       glacier_runoff_acc,                  &  !! Glacier runoff (rain+snow/ice melt) (accumulated) [kg/(m^2 s)]
       surface_temperature,                 &  !! Temperature of the surface
       surface_temperature_acc,             &  !! Accumulated temperature of the surface
       surface_temperature_old,             &  !! Temperature of the surface at time step t-dt
       surface_temperature_unfiltered,      &  !! Temperature of the surface at time step t+dt (unfiltered)
       radiative_temperature,               &  !! Temp for radiation derived from dry_satic_energy_new
       sat_surface_specific_humidity,       &  !! Saturated surface specific humidity
       evapotranspiration,                  &  !! Evaporation and transpiration (for time step)
       evapotranspiration_acc,              &  !! Evaporation and transpiration (accumulated)
       evaporation_pot,                     &  !! Potential evaporation (for time step)
       transpiration,                       &  !! Transpiration (for time step)
       transpiration_acc,                   &  !! Transpiration (accumulated)
       sensible_heat_flux,                  &  !! Sensible heat flux (for time step) [W/m^2]
       sensible_heat_acc,                   &  !!                  (accumulated)
       latent_heat_flux,                    &  !! Latent heat flux (for time step) [W/m^2]
       latent_heat_acc,                     &  !!                  (accumulated)
       ground_heat_flux,                    &  !! Ground heat flux (for time step) [W/m^2]
       ground_heat_flux_acc,                &  !!                  (accumulated)
       heat_capacity,                       &  !! Surface heat capacity [J/m**2/K]
       dry_static_energy_new,               &  !! ...
       qair,                                &  !! surface specific humidity (for timestep) [g/g]           
       qair_acc,                            &  !! surface specific humidity (accumulated) [g/g]           
       wetskin_fract,                       &  !! Mean Wet Skin Fraction
       relative_moisture,                   &  !! relative water content of the soil
       moisture,                            &  !! Moisture content of the unfrozen soil within the root zone [m]
       ! Variables for permafrost model
       nsnow_layer,                         &  !! previous number of snow layers
       snow_density,                        &  !! flexible snow density
       psnow_old,                           &  !! previous snow water equivalent
       k_snow,                              &  !! snow heat conductivity
       c_snow,                              &  !! snow volumetric heat capacity
       thaw_depth
 
     REAL(dp), POINTER, DIMENSION(:,:)   :: &  !! (nland, nsoil) - multi-layer soil scheme
       layer_moisture                           !! layer_moisture(i) = Soil Moisture content in i. soil layer [m]

     REAL(dp), POINTER, DIMENSION(:,:,:) :: &  !! (nland, ntsoil, ntiles)
       soil_temperature,                    &  !! Soil temperature
       c_soil_temperature,                  &  !! Coefficient to compute soil_temperature
       d_soil_temperature,                  &  !!                 "
       ! Variables for permafrost model
       ice_content,                         &  !! ice content of the soil layer
       heatcap,                             &  !! soil volumetric heat capacity
       heatcond,                            &  !! soil heat conductivity
       snow_temperature,                    &  !! Snow layer temperature
       c_snow_temperature,                  &  !! Coefficient to compute snow temperature
       d_snow_temperature                      !!                 "

  END TYPE soil_type

!----------------------------------------------------------------------------------------------------------
  TYPE soil_diag_type
    REAL(dp), POINTER, DIMENSION(:) ::   &
      surface_temperature,               &
      surface_temperature_acc,           &
      surface_radiative_temp,            &
      sat_surface_specific_humidity,     &
      ground_heat_flux_acc,              &
      ground_heat_flux_accw,             &
      heat_capacity,                     &
      evapotranspiration_acc,            &
      transpiration_acc,                 &
      sensible_heat_acc,                 &
      latent_heat_acc,                   &
      albedo,                            &
      skin_reservoir,                    &
      snow,                              &
      snow_accum_acc,                    &
      snow_accum_accw,                   &
      snow_fract,                        &
      snow_melt,                         & ! for dry deposition modules
      snow_melt_acc,                     &
      snow_melt_accw,                    &
      glacier_melt_accw,                 &
      runoff_acc,                        &
      drainage_acc,                      &
      runoff_accw,                       &
      drainage_accw,                     &
      glacier_depth,                     &
      glacier_precip_minus_evap_acc,     &
      glacier_precip_minus_evap_accw,    &
      glacier_runoff_acc,                &
      glacier_runoff_accw,               &
      net_radiation,                     &  !! net radiation at surface [Wm-2]
      qair,                              &  !! realised surface humidity (needed for offline model forcing) [-]
      qair_acc,                          &  !! realised surface humidity accumulated (needed for offline model forcing) [-]
      tair,                              &  !! lowest level temperature [K]  (needed for offline forcing with mvstreams)
      tair_max,                          &  !! lowest level max temperature [K] (needed for offline model forcing)
      tair_min,                          &  !! lowest level min temperature [K] (needed for offline model forcing)
      wind_lal,                          &  !! lowest level wind speed  [m/s] (needed for offline forcing with mvstreams)
      wind_lal_acc,                      &  !! lowest level wind speed  [m/s] (needed for offline model forcing)
      reduced_evap_acc,                  &  !! Potential for reducing evaporation (accumulated) - multi-layer soil scheme
      ! Soil water balance diagnostic variables = Fluxes - changes in surface water storage
      bare_soil_evap_acc,                &  !! Bare soil evaporation (accumulated) - diagnostic
      snow_evap_acc,                     &  !! Snow evaporation (accumulated) - diagnostic
      skin_reservoir_evap_acc,           &  !! Skin reservoir evaporation (accumulated) - diagnostic
      evap_deficit_acc,                  &  !! Evaporation deficit flux due to inconsistent treatment of snow and skin evaporation
      precip_acc,                        &  !! Precipitation received by JSBACH (accumulated) - diagnostic
      water_balance,                     &
      storage_pre,                       &
      csat,                              &  !!
      cair,                              &  !!
      csat_transpiration,                &  !!
      wetskin_fract,                     &  !! Mean Wet Skin Fraction - diagnostic
      moisture,                          &
      ! Variables for permafrost model
      nsnow_layer,                       &
      snow_density,                      &
      snow_conductivity,                 &
      snow_capacity,                     &
      thaw_depth

    REAL(dp), POINTER, DIMENSION(:,:) :: &
      soil_temperature,                  &
      ice_content,                       &
      heatcap,                           &
      heatcond,                          &
      snow_temperature

  END TYPE soil_diag_type
  !
  !----------------------------------------------------------------------------------------------------------
  TYPE soil_options_type
    INTEGER  :: nsoil                    !! Number of (water) soil layers
    LOGICAL  :: withPermafrost           !! Use permafrost model?
    LOGICAL  :: ldiag                    !! Output soil water balance diagnostics
    REAL(dp) :: SkinReservoirMax         !! Maximum content of skin reservoir over bare soil (set to 0. to disable)
    REAL(dp) :: MoistureFractCritical    !! Fractional soil moisture content at the critical point (fraction of maximum moisture)
    REAL(dp) :: MoistureFractWilting     !! Fractional soil moisture content at the wilting point (fraction of maximum moisture)
    REAL(dp) :: MaxMoistureLimit         !! Upper limit for maximum soil moisture content
    REAL(dp) :: CriticalSnowDepth        !! Critical snow depth for correction of new surface temperature [m water equivalent]
    LOGICAL  :: lbsoil                   !! Separate handling of bare soil moisture for bare soil evaporation
    INTEGER  :: nsnow                    !! Maximum number of snow layers
    LOGICAL  :: lsnow                    !! Consider snow in soil thermal calculations
    LOGICAL  :: ldynsnow                 !! Calculate snow parameters dynamically
    LOGICAL  :: lorganic                 !! Consider organic layers in soil thermal calculations
    LOGICAL  :: ldynorg                  !! Calculate organic layer parameters dynamically
    LOGICAL  :: lfreeze                  !! Consider freezing and thawing in thermal soil calculations
    LOGICAL  :: lsupercool               !! Allow for supercooled soil water
    LOGICAL  :: lheatcap                 !! Calculate soil heat capacity dynamically
    LOGICAL  :: lheatcond                !! Calculate soil heat conductivity dynamically
    INTEGER  :: HcondScheme              !! Scheme for dry soil heat conductivity
    LOGICAL  :: HeatCapMap               !! Read dry soil heat capacity from initial file (T) or calculate from soil type (F)
  END TYPE soil_options_type
  !
  !----------------------------------------------------------------------------------------------------------
  INTEGER,  SAVE :: nsoil  = -1    !! Number of soil layers (for water balance)
  INTEGER,  SAVE :: ntsoil = -1    !! Number of soil layers (for energy balance)
  INTEGER,  SAVE :: nsnow  = -1    !! Number of snow layers (for energy balance)
  LOGICAL,  SAVE :: ldiag_soil = .true. ! Add extra diagnostics for soil to land stream output
  TYPE(t_stream),  POINTER :: IO_soil      !! Memory stream for soil model state
  TYPE(t_stream),  POINTER :: IO_diag      !! Memory stream for soil diagnostic output
  TYPE(t_stream),  POINTER :: IO_accw      !! Memory stream for weighted accumulated soil diagnostic output
  TYPE(FILE_INFO), SAVE    :: soil_file    !! Input file for soil parameters

  TYPE(soil_diag_type),    SAVE :: soil_diag
  TYPE(soil_options_type), SAVE :: soil_options


CONTAINS
  !
  !=================================================================================================
  !
  SUBROUTINE config_soil()

    USE mo_jsbach,           ONLY: nml_unit
    USE mo_namelist,         ONLY: position_nml, POSITIONED
    USE mo_netcdf,           ONLY: nf_max_name
    USE mo_util_string,      ONLY: tolower

    !! locals
    INTEGER :: read_status, f_unit

    !! Namelist Parameters
    INTEGER  ::  nsoil             !! Number of soil (water) layers. Intented to be the same as the temperature layers (ntsoil)
    LOGICAL  ::  with_permafrost   !! Use permafrost model?
    LOGICAL  ::  ldiag             !! Output soil water balance diagnostics
    REAL(dp) ::  skin_res_max      !! Maximum content of skin reservoir of bare soil [m]
    REAL(dp) ::  moist_crit_fract  !! Critical value of soil moisutre
    REAL(dp) ::  moist_wilt_fract  !! Soil moisture content at permanent wilting point
    REAL(dp) ::  moist_max_limit   !! Upper limit for maximum soil moisture content
    REAL(dp) ::  crit_snow_depth   !! Critical snow depth for correction of surface temperature for melt 
    LOGICAL  ::  lbsoil            !! Separate handling of bare soil moisture for bare soil evaporation
    INTEGER  ::  nsnow             !! Maximum number of snow layers
    LOGICAL  ::  lsnow             !! Consider snow in soil thermal calculations
    LOGICAL  ::  ldynsnow          !! Calculate snow parameters dynamically
    LOGICAL  ::  lorganic          !! Consider organic layers in soil thermal calculations
    LOGICAL  ::  ldynorg           !! Calculate organic layer parameters dynamically
    LOGICAL  ::  lfreeze           !! Consider freezing and thawing in thermal soil calculations
    LOGICAL  ::  lsupercool        !! Allow for supercooled soil water
    LOGICAL  ::  lheatcap          !! Calculate soil heat capacity dynamically
    LOGICAL  ::  lheatcond         !! Calculate soil heat conductivity dynamically
    CHARACTER(nf_max_name) :: hcond_scheme !! scheme for dry soil heat conductivity ('none'/'johansen')
    LOGICAL  ::  heat_cap_map     !! read heat capacity from initial file (T) or calculate it from fao soil types (F)

    INCLUDE 'soil_ctl.inc'

    !! Read namelist soil_ctl

    IF (p_parallel_io) THEN

      ! define default values
      nsoil            = 1
      ldiag            = .TRUE.
      skin_res_max     = 2.E-4_dp
      moist_crit_fract = 0.75_dp
      moist_wilt_fract = 0.35_dp
      moist_max_limit  = -1.0_dp
      crit_snow_depth  = 5.85036E-3_dp
      lbsoil           = .FALSE.
      with_permafrost  = .FALSE.
      nsnow            = 5
      lsnow            = .TRUE. 
      ldynsnow         = .FALSE.
      lorganic         = .FALSE.
      ldynorg          = .TRUE.
      lfreeze          = .TRUE.
      lsupercool       = .FALSE.
      lheatcap         = .TRUE.
      lheatcond        = .TRUE.
      hcond_scheme     = 'JOHANSEN'
      heat_cap_map     = .FALSE.

      f_unit = position_nml ('SOIL_CTL', nml_unit, status=read_status)
      SELECT CASE (read_status)
      CASE (POSITIONED)
        READ (f_unit, soil_ctl)
        CALL message('config_soil', 'Namelist SOIL_CTL: ')
        WRITE(nout, soil_ctl)
      END SELECT

      soil_options%nsoil = nsoil
      WRITE(message_text,*) 'Soil layers: ',soil_options%nsoil
      CALL message('config_soil',message_text)

      soil_options%withPermafrost = with_permafrost
      WRITE (message_text,*) 'With permafrost model: ', soil_options%withPermafrost
      CALL message('config_coil', message_text)
      IF (soil_options%withPermafrost .AND. nsoil == 1) &
         CALL finish ('config_soil', 'Permafrost calculation is only possible with the multi-layer soil scheme')
      IF (.NOT. with_permafrost .AND. nsnow > 1) THEN
         WRITE(message_text,*) 'Model is run without permafrost and snow is simulated by a single layer (NSNOW = 1)'
         CALL message('config_coil', message_text)
      END IF

      soil_options%ldiag = ldiag
      WRITE(message_text,*) 'Soil diagnostics: ',soil_options%ldiag
      CALL message('config_soil',message_text)

      soil_options%SkinReservoirMax = skin_res_max
      WRITE(message_text,*) 'Maximum content of (bare soil) skin reservoir: ', soil_options%SkinReservoirMax
      CALL message('config_soil', message_text)

      soil_options%MoistureFractCritical = moist_crit_fract
      WRITE(message_text,*) 'Fractional soil moisture at critical point: ', soil_options%MoistureFractCritical
      CALL message('config_soil', message_text)
     
      soil_options%MoistureFractWilting = moist_wilt_fract
      WRITE(message_text,*) 'Fractional soil moisture at permanent wilting point: ', soil_options%MoistureFractWilting
      CALL message('config_soil', message_text)

      soil_options%MaxMoistureLimit = moist_max_limit
      WRITE(message_text,*) 'Upper limit for maximum soil moisture content: ', soil_options%MaxMoistureLimit
      CALL message('config_soil', message_text)

      soil_options%CriticalSnowDepth = crit_snow_depth
      WRITE(message_text,*) 'Critical snow depth for correction of surface temperature for melt: ', soil_options%CriticalSnowDepth
      CALL message('config_soil', message_text)

      soil_options%lbsoil = lbsoil
      WRITE(message_text,*) 'Bare Soil Moisture Separation: ',soil_options%lbsoil
      CALL message('config_soil',message_text)
      IF (lbsoil .AND. nsoil == 1) THEN
         WRITE(message_text,*) 'WARNING: Option lbsoil is only tested with the multi-layer model'
         CALL message ('config_soil',message_text)
      END IF

      soil_options%nsnow = nsnow
      IF (with_permafrost .AND. nsnow <= 2) THEN
         WRITE(message_text,*) 'ERROR: You need at least three snow layers.'
         CALL finish ('config_soil',message_text)
      END IF
      soil_options%lsnow = with_permafrost .AND. lsnow
      soil_options%ldynsnow = with_permafrost .AND. ldynsnow
      soil_options%lorganic = with_permafrost .AND. lorganic
      soil_options%ldynorg = with_permafrost .AND. ldynorg
      soil_options%lfreeze = with_permafrost .AND. lfreeze
      soil_options%lsupercool = with_permafrost .AND. lsupercool
      soil_options%lheatcap = with_permafrost .AND. lheatcap
      soil_options%lheatcond = with_permafrost .AND. lheatcond
      SELECT CASE (TRIM(tolower(hcond_scheme)))
      CASE ('none')
        soil_options%HcondScheme = hcond_none
      CASE('johansen')
        ! the scheme only produces plausible results in combination with the dynamic calculation of the heat conductivity
        IF (lheatcond) THEN
          soil_options%HcondScheme = hcond_johansen
        ELSE
          WRITE(message_text,*) ' Johanson scheme only works in case of dynamic heat conductivities'
          CALL finish('config_soil', message_text)
        END IF
      CASE default
        WRITE(message_text,*) 'Option "', TRIM(hcond_scheme), '" not suppoerted as hcond_scheme.'
        CALL finish('config_soil', message_text)
      END SELECT
      soil_options%HeatCapMap = heat_cap_map
      WRITE(message_text,*) 'Map with volumetric heat capacity of dry soil read from file: ', soil_options%HeatCapMap
      CALL message('config_soil', message_text)
    END IF

    IF (p_parallel) THEN
      CALL p_bcast(soil_options%nsoil,p_io)
      CALL p_bcast(soil_options%WithPermafrost, p_io)
      CALL p_bcast(soil_options%ldiag,p_io)
      CALL p_bcast(soil_options%SkinReservoirMax,p_io)
      CALL p_bcast(soil_options%MoistureFractCritical,p_io)
      CALL p_bcast(soil_options%MoistureFractWilting,p_io)
      CALL p_bcast(soil_options%MaxMoistureLimit,p_io)
      CALL p_bcast(soil_options%CriticalSnowDepth,p_io)
      CALL p_bcast(soil_options%lbsoil,p_io)
      CALL p_bcast(soil_options%nsnow,p_io)
      CALL p_bcast(soil_options%lsnow,p_io)
      CALL p_bcast(soil_options%ldynsnow,p_io)
      CALL p_bcast(soil_options%lorganic,p_io)
      CALL p_bcast(soil_options%ldynorg,p_io)
      CALL p_bcast(soil_options%lfreeze,p_io)
      CALL p_bcast(soil_options%lsupercool,p_io)
      CALL p_bcast(soil_options%lheatcap,p_io)
      CALL p_bcast(soil_options%lheatcond,p_io)
      CALL p_bcast(soil_options%HcondScheme,p_io)
      CALL p_bcast(soil_options%HeatCapMap,p_io)
    END IF

  END SUBROUTINE config_soil

!=================================================================================================
!
  SUBROUTINE soil_init_io(grid, soil, IO_file_name)

    USE mo_netCDF, ONLY : add_dim
    USE mo_io,  ONLY: IO_open, IO_READ
    USE mo_jsbach_grid, ONLY: grid_type
    USE mo_netcdf,      ONLY: IO_inq_dimid, IO_inq_dimlen, nf_max_name, IO_inq_varid, IO_get_var_double

    TYPE(grid_type),         INTENT(in)    :: grid
    TYPE(soil_type),         INTENT(inout) :: soil
    CHARACTER(NF_MAX_NAME), INTENT(in) :: IO_file_name
    INTEGER :: IO_file_id, IO_var_id, IO_dim_id
    INTEGER               :: i, znlon, znlat
    REAL(dp)              :: snow_layer(soil%nsnow)
    REAL(dp), ALLOCATABLE :: cmid_in_cm(:)
    LOGICAL :: lcdelFound

    IF (debug) CALL message('soil_init_io','')

    IF (p_parallel_io) THEN

       ! Open ini file
       WRITE(message_text,*) 'Reading soil parameters from ', TRIM(IO_file_name)
       CALL message('init_soil_io',message_text)
       soil_file%opened = .FALSE.
       CALL IO_open(TRIM(IO_file_name), soil_file, IO_READ)
       IO_file_id = soil_file%file_id

       ! Check resolution
       CALL IO_inq_dimid  (IO_file_id, 'lat', IO_dim_id)
       CALL IO_inq_dimlen (IO_file_id, IO_dim_id, znlat)
       CALL IO_inq_dimid  (IO_file_id, 'lon', IO_dim_id)
       CALL IO_inq_dimlen (IO_file_id, IO_dim_id, znlon)

       IF (znlon /= grid%nlon .OR. znlat /= grid%nlat) THEN
          WRITE(message_text,*) 'Unexpected resolution:', znlon, znlat
          CALL finish('init_soil_io', message_text)
       ENDIF

       ! check number of soil layers
       IF (soil%nsoil > 1) THEN
          CALL IO_inq_dimid  (IO_file_id, 'soillev', IO_dim_id)
          CALL IO_inq_dimlen (IO_file_id, IO_dim_id, nsoil)
          IF (soil%nsoil /= nsoil) THEN
             WRITE(message_text,*) &
                  'Number of soil levels in ', TRIM(soil_file%file_name),' does not match namelist: ', nsoil, soil%nsoil
             CALL finish('init_soil_io', message_text)
          ENDIF
          ! thermal soil layers need to match hydrological soil layers
          soil%ntsoil = nsoil
       ELSE
          ! use five soil layers for thermal calculations with bucket scheme
          soil%ntsoil = 5
       END IF
    END IF
    IF (p_parallel) CALL p_bcast(soil%ntsoil,p_io)

    nsoil  = soil%nsoil
    ntsoil = soil%ntsoil
    nsnow  = soil%nsnow

    WRITE(message_text,*) 'Number of soil layers (water) :', nsoil
    CALL message('soil_init_io', message_text)
    WRITE(message_text,*) 'Number of soil layers (energy):', ntsoil
    CALL message('soil_init_io', message_text)


    ! Read layer thickness from initial file

    ALLOCATE (soil%cdel(ntsoil))
    ALLOCATE (soil%clev(ntsoil))
    ALLOCATE (soil%cmid(ntsoil))

    lcdelfound = .FALSE.
    IF (nsoil > 1) THEN
       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'soil_layer_depth', IO_var_id,lcdelfound)
          IF (lcdelfound) &
             CALL IO_get_var_double(IO_file_id, IO_var_id, soil%cdel)
       ENDIF
       IF (p_parallel) THEN
         CALL p_bcast(lcdelfound,p_io)
         CALL p_bcast(soil%cdel,p_io)
       ENDIF
    ENDIF
    IF (.NOT. lcdelfound) THEN
       ! for compatibility with jsbach initial files without soil layer dimension
       !   thickness of the five thermal soil layers 
       soil%cdel = (/0.065_dp, 0.254_dp, 0.913_dp, 2.902_dp, 5.700_dp/)
    END IF

    ! calculate lower boundary soil levels
    soil%clev(1) = soil%cdel(1)
    DO i=2,ntsoil
       soil%clev(i) = soil%clev(i-1) + soil%cdel(i)
    END DO

    ! calculate soil level mid points
    soil%cmid(1) = 0.5_dp * soil%cdel(1)
    DO i=2,ntsoil
       soil%cmid(i) = soil%clev(i-1) + 0.5_dp * soil%cdel(i)
    END DO

    ! Add dimensions
    IF (debug) CALL message('soil_init_io','Adding dimensions')

    ! levels need to be defined as integer numbers; digits after the decimal point are lost.
    ALLOCATE(cmid_in_cm(ntsoil))
    cmid_in_cm=REAL(NINT(soil%cmid*100),dp)

    IF (nsoil > 1) THEN
       CALL add_dim ("soil_layer",    nsoil,  longname='soil levels (water)', levtyp = 71, &
            units = 'cm', value = cmid_in_cm(:), indx = SOILLEV)
    ELSE
       CALL add_dim ("soil_layer",    nsoil,  longname='soil levels (water)', levtyp = 71, &
            units = 'cm', value = (/ 1._dp /), indx = SOILLEV)
    END IF
    IF (nsnow > 1) THEN
       DO i = 1, nsnow
         snow_layer(i) = REAL(i,dp)
       END DO
       CALL add_dim ("snow_layer",    nsnow,  longname='snow layers', levtyp = 74, &
            units = '1', value = snow_layer(:), indx = SNOWLEV)
    ELSE
       CALL add_dim ("snow_layer",    nsnow,  longname='snow layers', levtyp = 74, &
            units = '1', value = (/ 1._dp /), indx = SNOWLEV)
    END IF

    CALL add_dim ("belowsurface", ntsoil, longname='soil levels (energy)', levtyp = 111, &
         units = 'cm', value = cmid_in_cm(:), indx = BELOWSUR)

    WRITE(message_text,*) 'soil levels in output [cm]: ', FLOOR(cmid_in_cm)
    CALL message('soil_init_io', message_text)
    DEALLOCATE(cmid_in_cm)

  END SUBROUTINE soil_init_io
  !
  !=================================================================================================
  SUBROUTINE soil_init_memory(g_nland, l_nland, ntiles, useDynveg, soil, isStandalone, &
                              fileformat, fileztype, diag_stream, accw_stream, stream)

    USE mo_jsbach,      ONLY : missing_value, lpost_echam
    USE mo_linked_list, ONLY : LAND, TILES
    USE mo_memory_base, ONLY : new_stream, default_stream_setting, &
                               add =>add_stream_element
    USE mo_netCDF,      ONLY : max_dim_name
    USE mo_output,      ONLY : land_table

    INTEGER, INTENT(in)               :: g_nland, l_nland, ntiles
    LOGICAL, INTENT(in)               :: isStandalone
    INTEGER, INTENT(in)               :: fileformat
    INTEGER, INTENT(in)               :: fileztype
    LOGICAL, INTENT(in)               :: useDynveg
    TYPE(soil_type), INTENT(inout)    :: soil
    TYPE(t_stream), POINTER           :: diag_stream, accw_stream
    TYPE(t_stream), POINTER, OPTIONAL :: stream

    INTEGER                     :: dim1p(1), dim1(1)
    INTEGER                     :: dim2p(2), dim2(2)
    INTEGER                     :: dim3p(2), dim3(2)
    INTEGER                     :: dim4p(3), dim4(3)
    INTEGER                     :: dim5p(3), dim5(3)
    INTEGER                     :: dim6p(2), dim6(2)
    INTEGER                     :: dim7p(2), dim7(2)
    INTEGER                     :: dim8p(3), dim8(3)
    CHARACTER(LEN=max_dim_name) :: dim1n(1), dim2n(2), dim3n(2), dim4n(3), dim5n(3), dim6n(2), dim7n(2), dim8n(3)

    IF (ntiles < 1) CALL finish('soil_init_memory', 'tiles not initialized')
    IF (nsoil  < 1) CALL finish('soil_init_memory', 'soil layers not initialized (water balance)')
    IF (ntsoil < 1) CALL finish('soil_init_memory', 'soil layers not initialized (energy balance)')

    IF (ASSOCIATED(diag_stream)) THEN
       IO_diag => diag_stream
    ELSE
       CALL finish('soil_init_memory', 'Diagnostic stream not present')
    END IF

    IF (ASSOCIATED(accw_stream)) THEN
       IO_accw => accw_stream
    ELSE
       CALL finish('soil_init_memory', 'Aggregated diagnostic stream not present')
    END IF

    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'soil', filetype=fileformat, ztype=fileztype)
          ! Set default stream options
          CALL default_stream_setting(stream, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE.)
       ENDIF
       IO_soil => stream
    ELSE
       ! Add new stream
       CALL new_stream(IO_soil, 'soil', filetype=fileformat, ztype=fileztype)
       ! Set default stream options
       CALL default_stream_setting(IO_soil, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF


    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    dim2p = (/ l_nland, nsoil /)
    dim2  = (/ g_nland, nsoil /)
    dim2n(1) = 'landpoint'
    dim2n(2) = 'soil_layer'

    dim3p = (/ l_nland, ntiles /)
    dim3  = (/ g_nland, ntiles /)
    dim3n(1) = 'landpoint'
    dim3n(2) = 'tiles'

    dim4p = (/ l_nland, nsoil, ntiles /)
    dim4  = (/ g_nland, nsoil, ntiles /)
    dim4n(1) = 'landpoint'
    dim4n(2) = 'soil_layer'
    dim4n(3) = 'tiles'

    dim5p = (/ l_nland, ntsoil, ntiles /)
    dim5  = (/ g_nland, ntsoil, ntiles /)
    dim5n(1) = 'landpoint'
    dim5n(2) = 'belowsurface'
    dim5n(3) = 'tiles'

    dim6p = (/ l_nland, ntsoil /)
    dim6  = (/ g_nland, ntsoil /)
    dim6n(1) = 'landpoint'
    dim6n(2) = 'belowsurface'

    dim7p = (/ l_nland, nsnow /)
    dim7  = (/ g_nland, nsnow /)
    dim7n(1) = 'landpoint'
    dim7n(2) = 'snow_layer'

    dim8p = (/ l_nland, nsnow, ntiles /)
    dim8  = (/ g_nland, nsnow, ntiles /)
    dim8n(1) = 'landpoint'
    dim8n(2) = 'snow_layer'
    dim8n(3) = 'tiles'

    ! code numbers of this routine are 20-89, using the GRIB jsbach_table
    ! code numbers of this routine are 35-35, 39-40, 42, 44, 47, 49-50, 55, 57, 59-68, 72-74, 76-83, 85-86, 91-105, 
    !                                  using the GRIB land_table
    ! code numbers of this routine are 160-161, 206, 215, 218, 221-222, 228, using the GRIB accw_table
    CALL add(IO_soil, 'surface_temperature',     soil%surface_temperature,   longname='Land Surface Temperature', &
             units='K',                    ldims=dim3p, gdims=dim3, dimnames=dim3n, code=35, lpost=.FALSE.)
    CALL add(IO_soil, 'surface_temperature_acc', soil%surface_temperature_acc,   longname='Land Surface Temperature accumulated', &
             units='K',                    ldims=dim3p, gdims=dim3, dimnames=dim3n, code=139, laccu=.TRUE., lpost=.FALSE.,        &
             contnorest=.TRUE.)
    CALL add(IO_soil, 'surface_temperature_old', soil%surface_temperature_old, &
                                           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=36, lpost=.FALSE.)
    CALL add(IO_soil, 'surface_temp_unfiltered', soil%surface_temperature_unfiltered, &
                                           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=37, lpost=.FALSE.)
    CALL add(IO_soil, 'surface_radiative_temp',  soil%radiative_temperature, longname='Surface Radiative Temperature', &
                                           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=38, lpost=.FALSE.)
    CALL add(IO_soil, 'surface_qsat',    soil%sat_surface_specific_humidity, longname='Soil Specific Humidity at Saturation', &
             units='',                     ldims=dim3p, gdims=dim3, dimnames=dim3n, code=39, lpost=.FALSE.)
    CALL add(IO_soil, 'ground_heat_flux_inst',   soil%ground_heat_flux,      longname='Ground Heat Flux (inst.)', &
             units='W m-2',                ldims=dim3p, gdims=dim3, dimnames=dim3n, code=40, lpost=.FALSE.)
    CALL add(IO_soil, 'ground_heat_flux',        soil%ground_heat_flux_acc,  longname='Ground Heat Flux (ave.)', &
             units='W m-2',                ldims=dim3p, gdims=dim3, dimnames=dim3n, code=41, laccu=.TRUE., lpost=.FALSE.)
    CALL add(IO_soil, 'heat_capacity',           soil%heat_capacity,         longname='Heat Capacity of the Uppermost Soil Layer', &
             units='J m-2 K-1',            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=42, lpost=.FALSE.)
    CALL add(IO_soil, 'evapotranspiration_inst', soil%evapotranspiration,    longname='Evapotranspiration (inst.)', &
             units='kg m-2 s-1',           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=43, lpost=.FALSE.)
    CALL add(IO_soil, 'evapotranspiration',      soil%evapotranspiration_acc,longname='Evapotranspiration (avg)', &
             units='kg m-2 s-1',           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=44, laccu=.TRUE., lpost=.FALSE.)
    CALL add(IO_soil, 'evaporation_pot_inst',    soil%evaporation_pot,       longname='Potential Evaporation', &
             units='kg m-2 s-1',           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=45, lpost=.FALSE.)
    CALL add(IO_soil, 'transpiration_inst',      soil%transpiration,         longname='Transpiration (inst.)', &
             units='kg m-2 s-1',           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=75, lpost=.FALSE.)
    CALL add(IO_soil, 'transpiration',           soil%transpiration_acc,     longname='Transpiration (avg)', &
             units='kg m-2 s-1',           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=76, laccu=.TRUE., lpost=.FALSE.)
    CALL add(IO_soil, 'sensible_heat_flux_inst', soil%sensible_heat_flux,    longname='Sensible Heat Flux (inst.)', &
             units='W m-2',                ldims=dim3p, gdims=dim3, dimnames=dim3n, code=46, lpost=.FALSE.)
    CALL add(IO_soil, 'sensible_heat_flux',      soil%sensible_heat_acc,     longname='Sensible Heat Flux (avg)', &
             units='W m-2',                ldims=dim3p, gdims=dim3, dimnames=dim3n, code=47, laccu=.TRUE., lpost=.FALSE.)
    CALL add(IO_soil, 'latent_heat_flux_inst',   soil%latent_heat_flux,      longname='Latent Heat Flux (inst.)', &
             units='W m-2',                ldims=dim3p, gdims=dim3, dimnames=dim3n, code=48, lpost=.FALSE.)
    CALL add(IO_soil, 'latent_heat_flux',        soil%latent_heat_acc,       longname='Latent Heat Flux (avg)', &
             units='W m-2',                ldims=dim3p, gdims=dim3, dimnames=dim3n, code=49, laccu=.TRUE., lpost=.FALSE.)
    CALL add(IO_soil, 'soil_albedo',             soil%albedo,                longname='Soil Albedo', &
             units='',                     ldims=dim3p, gdims=dim3, dimnames=dim3n, code=50, lpost=.FALSE., lrerun=.FALSE.)
    CALL add(IO_soil, 'albedo_veg_vis',          soil%albedo_vegetation_vis, &
                                           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=51, lpost=.FALSE., lrerun=.FALSE.)
    CALL add(IO_soil, 'albedo_veg_nir',          soil%albedo_vegetation_nir, &
                                           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=52, lpost=.FALSE., lrerun=.FALSE.)
    CALL add(IO_soil, 'albedo_soil_vis',         soil%albedo_soil_vis, &
                                           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=53, lpost=.FALSE., lrerun=.FALSE.)
    CALL add(IO_soil, 'albedo_soil_nir',         soil%albedo_soil_nir, &
                                           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=54, lpost=.FALSE., lrerun=.FALSE.)
    CALL add(IO_soil, 'skin_reservoir',          soil%skin_reservoir,        longname='Water Content of the Skin Reservoir', &
             units='m',                    ldims=dim3p, gdims=dim3, dimnames=dim3n, code=55, lpost=.TRUE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_soil, 'soil_moisture',           soil%moisture,              longname='Water Content of the Soil', &
             units='m',                    ldims=dim3p, gdims=dim3, dimnames=dim3n, code=57, lpost=.FALSE.)
    CALL add(IO_soil, 'rel_soil_moisture',       soil%relative_moisture,     longname='Relative Water Content of the Soil', &
             units='',                     ldims=dim3p, gdims=dim3, dimnames=dim3n, code=58, lpost=.FALSE.)
    CALL add(IO_soil, 'snow',                    soil%snow,                  longname='Snow Depth in Water Equivalent', &
             units='m',                    ldims=dim3p, gdims=dim3, dimnames=dim3n, code=59, lpost=.FALSE.)
    CALL add(IO_soil, 'snow_fract',              soil%snow_fract,            longname='Snow Cover Fraction', &
             units='',                     ldims=dim3p, gdims=dim3, dimnames=dim3n, code=60, lpost=.FALSE.)
    CALL add(IO_soil, 'snow_accumulation', soil%snow_accum_acc,              longname='Snow accumulation (acc.)', &
             units='kg m-2 s-1',           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=61, laccu=.TRUE., lpost=.FALSE.)
    CALL add(IO_soil, 'snow_melt_inst',          soil%snow_melt,             longname='Snow Melt Flux', &
             units='m',                    ldims=dim3p, gdims=dim3, dimnames=dim3n, code=1, laccu=.FALSE., &
             contnorest=.TRUE., lpost=.FALSE.)
    CALL add(IO_soil, 'snow_melt',               soil%snow_melt_acc,         longname='Snow Melt Flux', &
             units='kg m-2 s-1',           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=62, laccu=.TRUE., lpost=.FALSE.)
    CALL add(IO_soil, 'snow_age',                soil%snow_age,              longname='Non-dimensional Snow Age', &
             units='',                     ldims=dim1p, gdims=dim1, dimnames=dim1n, code=77, laccu=.FALSE., &
             contnorest=.TRUE., lpost=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_soil, 'runoff',                  soil%runoff_acc,            longname='Surface Runoff and Drainage', &
             units='kg m-2 s-1',           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=63, laccu=.TRUE., lpost=.FALSE.)
    CALL add(IO_soil, 'drainage',                soil%drainage_acc,          longname='Drainage', &
             units='kg m-2 s-1',           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=64, laccu=.TRUE., lpost=.FALSE.)
    CALL add(IO_soil, 'glacier_depth',           soil%glacier_depth,         longname='Glacier Depth', &
             units='m',                    ldims=dim3p, gdims=dim3, dimnames=dim3n, code=65, lpost=.FALSE.)
    CALL add(IO_soil, 'glacier_precip_minus_evap', soil%glacier_precip_minus_evap_acc, &
             longname='Precipitation minus Evaporation over Glaciers', &
             units='m',                    ldims=dim3p, gdims=dim3, dimnames=dim3n, code=66, laccu=.TRUE., lpost=.FALSE.)
    CALL add(IO_soil, 'glacier_runoff',          soil%glacier_runoff_acc,    longname='Glacier Runoff',  &
             units='kg m-2 s-1',           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=67, laccu=.TRUE., lpost=.TRUE.,  &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_soil, 'soil_temperature',        soil%soil_temperature,      longname='Soil Temperature', &
             units='K',                    ldims=dim5p, gdims=dim5, dimnames=dim5n, code=68, lpost=.FALSE., leveltype=BELOWSUR)
    CALL add(IO_soil, 'c_soil_temp_coef',        soil%c_soil_temperature, &
                                           ldims=dim5p, gdims=dim5, dimnames=dim5n, code=69, lpost=.FALSE., leveltype=BELOWSUR)
    CALL add(IO_soil, 'd_soil_temp_coef',        soil%d_soil_temperature, &
                                           ldims=dim5p, gdims=dim5, dimnames=dim5n, code=70, lpost=.FALSE., leveltype=BELOWSUR)
    CALL add(IO_soil, 'dry_static_energy_new',   soil%dry_static_energy_new, &
                                           ldims=dim3p, gdims=dim3, dimnames=dim3n, code=71, lpost=.FALSE.)

    CALL add(IO_soil, 'csat',                    soil%csat, dim1p, dim1, dimnames=dim1n, code=72, lpost=.FALSE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_soil, 'cair',                    soil%cair, dim1p, dim1, dimnames=dim1n, code=73, lpost=.FALSE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_soil, 'qair',   soil%qair, &
             units='kg/kg',                ldims=dim3p, gdims=dim3, dimnames=dim3n, code=87, lpost=.FALSE., lrerun=.FALSE.)
    CALL add(IO_soil, 'qair_acc',   soil%qair_acc, &
             units='kg/kg',                ldims=dim3p, gdims=dim3, dimnames=dim3n, code=88, laccu=.TRUE., lpost=.FALSE., &
             contnorest=.TRUE.)
    CALL add(IO_soil, 'z0m',   soil%z0m, &
             units='m',                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=93, laccu=.FALSE., lpost=.FALSE., &
             contnorest=.TRUE.)
    CALL add(IO_soil, 'z0h',   soil%z0h, &
             units='m',                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=94, laccu=.FALSE., lpost=.FALSE., &
             contnorest=.TRUE.)
    CALL add(IO_soil, 'zril_old',   soil%zril_old, &
             units='',                 ldims=dim1p, gdims=dim1, dimnames=dim1n, code=109,laccu=.FALSE., lpost=.FALSE., &
             lrerun=.TRUE., contnorest=.TRUE.)
    IF (useDynveg) THEN
       CALL add(IO_soil, 'relative_humidity_air_inst',soil%relative_humidity_air,longname='Relative Humidity', &
            units='',                      ldims=dim1p, gdims=dim1, dimnames=dim1n, code=56, lpost=.TRUE.,     &
            lmiss=.TRUE., missval=missing_value)
    ENDIF
    CALL add(IO_soil, 'csat_transpiration',soil%csat_transpiration, dim1p, dim1, dimnames=dim1n, code=74, lpost=.FALSE., &
            lmiss=.TRUE., missval=missing_value)

    IF (ldiag_soil) THEN
      ! all these variables are douplicated in the IO_diag (i.e. land) stream, therefore lpost=F here.
      CALL add(IO_soil, 'water_balance',         soil%water_balance, longname='Soil Water Balance',            &
           units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=78,  &
           lpost=.FALSE., lmiss=.TRUE., missval=missing_value, contnorest=.TRUE.)
      CALL add(IO_soil, 'water_storage_pre',     soil%storage_pre, longname='Water storage of previous time step', &
           units='m',                        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=77,      &
           lpost=.FALSE., lmiss=.TRUE., missval=missing_value, contnorest=.TRUE.)
      CALL add(IO_soil, 'bare_soil_evaporation', soil%bare_soil_evap_acc, longname='Bare Soil Evaporation (avg)', &
           units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=79,     &
           lpost=.FALSE., contnorest=.TRUE.)
      CALL add(IO_soil, 'snow_evaporation',      soil%snow_evap_acc, longname='Snow Evaporation (avg)',        &
           units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=80,  &
           lpost=.FALSE., contnorest=.TRUE.)
      CALL add(IO_soil, 'skin_evaporation',      soil%skin_reservoir_evap_acc, longname='Skin Reservoir Evaporation (avg)', &
           units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=81,               &
           lpost=.FALSE., contnorest=.TRUE.)
      CALL add(IO_soil, 'evaporation_deficit',   soil%evap_deficit_acc, longname='Evaporation deficit flux (avg)', &
           units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=102,     &
           lpost=.FALSE., contnorest=.TRUE.)
      CALL add(IO_soil, 'precip_acc',            soil%precip_acc, longname='Precipitation received by JSBACH', &
           units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=83,  &
           lpost=.FALSE., lmiss=.TRUE., missval=missing_value, contnorest=.TRUE.)
      IF (nsoil > 1) & 
        CALL add(IO_soil, 'reduced_evaporation', soil%reduced_evap_acc, longname='Potential for reducing Evaporation (avg)', &
             units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=82,              &
             lpost=.FALSE., contnorest=.TRUE.)
    ENDIF
    CALL add(IO_soil, 'wetskin_fract',              soil%wetskin_fract,            longname='Wet Skin Fraction', &
             units='',                     ldims=dim3p, gdims=dim3, dimnames=dim3n, code=60, lpost=.FALSE., contnorest=.true.)

    IF (nsoil > 1) THEN
      CALL add(IO_soil, 'layer_moisture', soil%layer_moisture,longname='Soil Moisture Content', &
               units='m',        ldims=dim2p, gdims=dim2, dimnames=dim2n, code=84, lpost=.TRUE., leveltype=SOILLEV, &
               contnorest=.TRUE., lmiss=.TRUE., missval=missing_value)
    ENDIF

    IF (soil_options%lbsoil) THEN
       CALL add(IO_soil, 'bare_soil_moisture',         soil%bare_soil_moisture, longname='Soil Moisture for bare soil fraction', &
             units='m',                   ldims=dim1p, gdims=dim1, dimnames=dim1n, code=89, laccu=.FALSE., &
             contnorest=.TRUE., lpost=.TRUE., lmiss=.TRUE., missval=missing_value)
    ENDIF

    CALL add(IO_soil, 'soil_ice_content',   soil%ice_content, longname='Ice content of the soil layers',                     &
             units='m',             ldims=dim5p, gdims=dim5, dimnames=dim5n, code=20, lpost=.FALSE., leveltype=BELOWSUR,     &
             contnorest=.TRUE.)
    IF (soil_options%withPermafrost) THEN
      CALL add(IO_soil, 'soil_heatcap',       soil%heatcap,     longname='Heat capacity of the soil layers (modified)',      &
               units='J m-3 K-1',     ldims=dim5p, gdims=dim5, dimnames=dim5n, code=21, lpost=.FALSE., leveltype=BELOWSUR)
      CALL add(IO_soil, 'soil_heatcond',      soil%heatcond,    longname='Heat conductivity  of the soil layers (modified)', &
               units='J m-1 s-1 K-1', ldims=dim5p, gdims=dim5, dimnames=dim5n, code=22, lpost=.FALSE., leveltype=BELOWSUR)
      CALL add(IO_soil, 'c_snow_temperature', soil%c_snow_temperature,               &
                                      ldims=dim8p, gdims=dim8, dimnames=dim8n, code=23, lpost=.FALSE.)
      CALL add(IO_soil, 'd_snow_temperature', soil%d_snow_temperature,               &
                                      ldims=dim8p, gdims=dim8, dimnames=dim8n, code=24, lpost=.FALSE.)
      CALL add(IO_soil, 'snow_temperature',   soil%snow_temperature,               &
                                      ldims=dim8p, gdims=dim8, dimnames=dim8n, code=25, lpost=.FALSE.)
      CALL add(IO_soil, 'thaw_depth',         soil%thaw_depth,           &
               units='m',             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=29, lpost=.FALSE.)
      CALL add(IO_soil, 'nsnow_layer',        soil%nsnow_layer,              &
               units='',              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=30, lpost=.FALSE.)
      CALL add(IO_soil, 'snow_density',       soil%snow_density,          &
               units='kg/m^3',        ldims=dim3p, gdims=dim3, dimnames=dim3n, code=31, lpost=.FALSE.)
      CALL add(IO_soil, 'psnow_old',          soil%psnow_old,              &
               units='m',             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=32, lpost=.FALSE.)
      CALL add(IO_soil, 'k_snow',             soil%k_snow,               &
               units='W/K/m',         ldims=dim3p, gdims=dim3, dimnames=dim3n, code=33, lpost=.FALSE.)
      CALL add(IO_soil, 'c_snow',             soil%c_snow,               &
               units='J/K/m3',        ldims=dim3p, gdims=dim3, dimnames=dim3n, code=34, lpost=.FALSE.)
    ENDIF

    CALL add(IO_diag, 'surface_temperature', soil_diag%surface_temperature,   longname='Land Surface Temperature',             &
         units='K',                        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=35, lpost=lpost_echam, &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'surface_temperature_acc', soil_diag%surface_temperature_acc,                                         &
         longname='Land Surface Temperature accumulated', units='K', ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE., &
         code=139, lpost=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'surface_radiative_temp', soil_diag%surface_radiative_temp,longname='Surface Radiative Temperature', &
         units='K',                        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=36, &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'surface_qsat',     soil_diag%sat_surface_specific_humidity,longname='Soil Specific Humidity at Saturation',&
         units='kg/kg',                    ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=39, lpost=.FALSE.,        &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'ground_heat_flux',        soil_diag%ground_heat_flux_acc,  longname='Ground Heat Flux (avg)', &
         units='W/m^2',                    ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=40, lpost=.TRUE.,         &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'heat_capacity',           soil_diag%heat_capacity,    longname='Heat Capacity of the Uppermost Soil Layer',&
         units='J m-2 K-1',                ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=42, lpost=.FALSE.,        &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'evapotranspiration',      soil_diag%evapotranspiration_acc,longname='Evapotranspiration (avg)',            &
         units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=44,                       &
         lmiss=.TRUE., missval=missing_value, bits=24)
    CALL add(IO_diag, 'transpiration',           soil_diag%transpiration_acc,     longname='Transpiration (avg)',                 &
         units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=76,  lpost=.TRUE.,        &
         lmiss=.TRUE., missval=missing_value, bits=24)
    CALL add(IO_diag, 'sensible_heat_flx',       soil_diag%sensible_heat_acc,     longname='Sensible Heat Flux (avg)',            &
         units='W/m^2',                    ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=47,  lpost=lpost_echam,   &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'latent_heat_flx',         soil_diag%latent_heat_acc,       longname='Latent Heat Flux (avg)',              &
         units='W/m^2',                    ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=49,  lpost=lpost_echam,   &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'soil_albedo',             soil_diag%albedo,                longname='Soil Albedo',                         &
         units='',                         ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=50,  lpost=.FALSE.,       &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'soil_moisture',           soil_diag%moisture,              longname='Water Content of the Soil',           &
         units='m',                        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=57, lpost=lpost_echam,    &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'skin_reservoir',          soil_diag%skin_reservoir,        longname='Water Content of the Skin Reservoir', &
         units='m',                        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=55, lpost=lpost_echam,    &
         lmiss=.TRUE., missval=missing_value, bits=24)
    CALL add(IO_diag, 'snow',                    soil_diag%snow,                  longname='Snow Depth in Water Equivalent',   &
         units='m',                        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=59, lpost=lpost_echam, &
         lmiss=.TRUE., missval=missing_value, bits=24)
    CALL add(IO_diag, 'snow_accumulation',       soil_diag%snow_accum_acc,        longname='Snow Accumulation (acc.)',         &
         units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=61, lpost=lpost_echam, &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'snow_fract',              soil_diag%snow_fract,            longname='Snow Cover Fraction', &
         units='',                         ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=60,       &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'snow_melt_inst',          soil_diag%snow_melt,             longname='Snow Melt',                        &
         units='m',                        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=1, lpost=lpost_echam,  &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'snow_melt',               soil_diag%snow_melt_acc,         longname='Snow Melt',                        &
         units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=62, lpost=lpost_echam, &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'runoff',                  soil_diag%runoff_acc,            longname='Surface Runoff and Drainage',      &
         units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=63, lpost=lpost_echam, &
         lmiss=.TRUE., missval=missing_value, bits=24)
    CALL add(IO_diag, 'drainage',                soil_diag%drainage_acc,          longname='Drainage',                         &
         units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=64, lpost=lpost_echam, &
         lmiss=.TRUE., missval=missing_value, bits=24)
    CALL add(IO_diag, 'glacier_depth',           soil_diag%glacier_depth,         longname='Glacier Depth',                    &
         units='m',                        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=65, lpost=lpost_echam, &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'glacier_precip_minus_evap', soil_diag%glacier_precip_minus_evap_acc,                                    &
         longname='Precipitation minus Evaporation over Glaciers',                                                             &
         units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=66, lpost=lpost_echam, &
         lmiss=.TRUE., missval=missing_value, bits=24)
    CALL add(IO_diag, 'glacier_runoff',          soil_diag%glacier_runoff_acc,    longname='Glacier Runoff',   &
         units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=67,    &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'soil_temperature',        soil_diag%soil_temperature,      longname='Soil Temperature', &
         units='K',                        ldims=dim6p, gdims=dim6, dimnames=dim6n, laccu=.FALSE., code=68,    &
         lmiss=.TRUE., missval=missing_value)
         !sz added for sitelevel evaluation
    CALL add(IO_diag, 'net_radiation',          soil_diag%net_radiation,          longname='Net Radiation at Surface', &
         units='W/m^2',                    ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=85,  lpost=lpost_echam, &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'qair_lal', soil_diag%qair, longname='Specific humidity lowest atmospheric level',     &
         units='m/s', ldims=dim1p, gdims=dim1, dimnames=dim1n, code=108, laccu=.FALSE., lpost=.FALSE., contnorest=.TRUE., &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'surface_qair',          soil_diag%qair_acc,          longname='Surface specific humidity, accu.', &
         units='kg/kg',                    ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=91,  lpost=.TRUE., &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'tair', soil_diag%tair, longname='Temperature of lowest atmospheric level',                       &
         units='K', ldims=dim1p, gdims=dim1, dimnames=dim1n, code=103, laccu=.FALSE., lpost=.FALSE., contnorest=.TRUE., &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'tair_max', soil_diag%tair_max, longname='Maximum temperature lowest atmospheric level',         &
         units='K', ldims=dim1p, gdims=dim1, dimnames=dim1n, code=104, laccu=.FALSE., lpost=.TRUE., contnorest=.TRUE., &
         lmiss=.TRUE., reset=-99._dp, missval=missing_value)
    CALL add(IO_diag, 'tair_min', soil_diag%tair_min, longname='Minimum temperature lowest atmospheric level',         &
         units='K', ldims=dim1p, gdims=dim1, dimnames=dim1n, code=105, laccu=.FALSE., lpost=.TRUE., contnorest=.TRUE., &
         lmiss=.TRUE., reset=999._dp, missval=missing_value)
    CALL add(IO_diag, 'wind_lal', soil_diag%wind_lal, longname='Windspeed lowest atmospheric level',                      &
         units='m/s', ldims=dim1p, gdims=dim1, dimnames=dim1n, code=111, laccu=.FALSE., lpost=.FALSE., contnorest=.TRUE., &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'wind_lal_acc', soil_diag%wind_lal_acc, longname='Windspeed lowest atmospheric level',              &
         units='m/s', ldims=dim1p, gdims=dim1, dimnames=dim1n, code=106, laccu=.TRUE., lpost=.TRUE., contnorest=.TRUE.,   &
         lmiss=.TRUE., missval=missing_value)



    IF (ldiag_soil) THEN
      CALL add(IO_diag, 'water_balance', soil_diag%water_balance,    longname='Soil Water Balance', &
           units='kg/m^2s', ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=78,        &
           lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
      CALL add(IO_diag, 'water_storage_pre', soil_diag%storage_pre,  longname='Water storage of previous time step', &
           units='m',    ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=77, lpost=.TRUE.)
      CALL add(IO_diag, 'bare_soil_evaporation', soil_diag%bare_soil_evap_acc, longname='Bare Soil Evaporation (avg)', &
           units='kg/m^2s', ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=79,                           &
           lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
      CALL add(IO_diag, 'snow_evaporation', soil_diag%snow_evap_acc, longname='Snow Evaporation (avg)', &
           units='kg/m^2s', ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=80,            &
           lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
      CALL add(IO_diag, 'skin_evaporation',soil_diag%skin_reservoir_evap_acc, longname='Skin Reservoir Evaporation (avg)', &
           units='kg/m^2s', ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=81,                               &
           lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
      CALL add(IO_diag, 'evaporation_deficit',soil_diag%evap_deficit_acc,     longname='Evaporation deficit flux (avg)', &
           units='kg/m^2s', ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=102,                            &
           lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
      CALL add(IO_diag, 'precip_acc', soil_diag%precip_acc, longname='Precipitation received by JSBACH', &
           units='kg/m^2s', ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=83,             &
           lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
      IF (nsoil > 1) &
          CALL add(IO_diag, 'reduced_evaporation',soil_diag%reduced_evap_acc, longname='Potential for reducing Evaporation (avg)',&
               units='kg/m^2s', ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=82,                                  &
               lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
      CALL add(IO_diag, 'wetskin_fract',              soil_diag%wetskin_fract,            longname='Mean Wet Skin Fraction', &
           units='',                         ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=86,       &
           lmiss=.TRUE., missval=missing_value)
      CALL add(IO_diag, 'csat',                       soil_diag%csat,        longname='csat exchange coefficient',  &      
           units='',                         ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=72,       &
           lmiss=.TRUE., missval=missing_value)
      CALL add(IO_diag, 'cair',                       soil_diag%cair,        longname='cair exchange coefficient',  &      
           units='',                         ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=73,       &
           lmiss=.TRUE., missval=missing_value)
      CALL add(IO_diag, 'csat_transpiration',         soil_diag%csat_transpiration,  longname='csat transpiration coefficient', &
           units='',                         ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=74,       &
           lmiss=.TRUE., missval=missing_value)
      IF (soil_options%withPermafrost) THEN
        CALL add(IO_diag, 'soil_ice_content', soil_diag%ice_content,      longname='Soil Ice content of layer i [m]',          &
             units='m',      ldims=dim6p, gdims=dim6, dimnames=dim6n, laccu=.FALSE., code=92,   &
             lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
        CALL add(IO_diag, 'heatcap',          soil_diag%heatcap,          longname='Soil heat capacity of layer i [J/m^3K]',   &
             units='J/m^3K', ldims=dim6p, gdims=dim6, dimnames=dim6n, laccu=.FALSE., code=93,   &
             lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
        CALL add(IO_diag, 'heatcond',         soil_diag%heatcond,         longname='Soil heat conductivity of layer i [J/ms]', &
             units='J/ms',   ldims=dim6p, gdims=dim6, dimnames=dim6n, laccu=.FALSE., code=94,   &
             lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
        CALL add(IO_diag, 'nsnow_layer',      soil_diag%nsnow_layer,      longname='Number of snow layers',                    &
             units='',       ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE.,  code=96,  &
             lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
        CALL add(IO_diag, 'snow_density',     soil_diag%snow_density,     longname='Snow density',                             &
             units='kg/m3',  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE.,  code=97,  &
             lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
        CALL add(IO_diag, 'thaw_depth',       soil_diag%thaw_depth,       longname='Active layer depth [cm]',                  &
             units='cm',     ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE.,  code=98,  &
             lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
        CALL add(IO_diag, 'snow_temperature', soil_diag%snow_temperature, longname='Snow Temperature',                         &
             units='K',      ldims=dim7p, gdims=dim7, dimnames=dim7n, laccu=.FALSE., code=99,   &
             lmiss=.TRUE.,               missval=missing_value) 
        CALL add(IO_diag, 'snow_conductivity',soil_diag%snow_conductivity,longname='Snow heat conductivity [W/K/m]',           &
             units='W/K/m',  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE.,  code=100, &
             lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
        CALL add(IO_diag, 'snow_capacity',    soil_diag%snow_capacity,    longname='Snow heat capacity [J/K/m3]',              &
             units='J/K/m3', ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE.,  code=101, &
             lpost=.TRUE., lmiss=.TRUE., missval=missing_value, bits=24)
      ENDIF
    ENDIF

    CALL add(IO_accw, 'ground_heat_flux',        soil_diag%ground_heat_flux_accw,  longname='Ground Heat Flux (acc+weighted)', &
         units='W/m^2', ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=206, lpost=.NOT.isStandalone, &
         lmiss=.TRUE., missval=missing_value, contnorest=.TRUE.)
    CALL add(IO_accw, 'snow_accumulation',       soil_diag%snow_accum_accw,        longname='Snow Accumulation (acc+weighted)', &
         units='kg/m^2s',ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=222, lpost=.NOT.isStandalone, &
         lmiss=.TRUE., missval=missing_value, contnorest=.TRUE.)
    CALL add(IO_accw, 'snow_melt',               soil_diag%snow_melt_accw,        longname='Land Snow Melt (acc+weighted)',  &
         units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=218, &
         lpost=.NOT.isStandalone, lmiss=.TRUE., missval=missing_value, contnorest=.TRUE.)
    CALL add(IO_accw, 'glacier_melt',            soil_diag%glacier_melt_accw,     longname='Glacier Ice Melt (acc+weighted)', &
         units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=228, &
         lpost=.NOT.isStandalone, lmiss=.TRUE., missval=missing_value, contnorest=.TRUE.)
    CALL add(IO_accw, 'glacier_runoff',          soil_diag%glacier_runoff_accw,   longname='Glacier Runoff (acc+weighted)',   &
         units='kg/m^2s',                  ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=215,    &
         lpost=.NOT.isStandalone, lmiss=.TRUE., missval=missing_value, contnorest=.TRUE.)
    CALL add(IO_accw, 'runoff',                  soil_diag%runoff_accw,     longname='Surface Runoff and Drainage (acc+weighted)', &
         units='kg/m^2s', ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=160, &
         lpost=.NOT.isStandalone, lmiss=.TRUE., missval=missing_value, bits=24, contnorest=.TRUE.)
    CALL add(IO_accw, 'drainage',                soil_diag%drainage_accw,    longname='Drainage (acc+weighted)',     &
         units='kg/m^2s', ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=161, &
         lpost=.NOT.isStandalone, lmiss=.TRUE., missval=missing_value, bits=24, contnorest=.TRUE.)
    CALL add(IO_accw, 'glacier_precip_minus_evap', soil_diag%glacier_precip_minus_evap_accw,                                  &
         longname='P-E over land ice (acc+weighted)',                                                                         &
         units='kg/m^2s', ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.TRUE.,  code=221, lpost=.NOT.isStandalone, &
         lmiss=.TRUE., missval=missing_value, bits=24, contnorest=.TRUE.)

  END SUBROUTINE soil_init_memory
  !
  !=================================================================================================
  !
  SUBROUTINE init_soil(grid, domain, soil_param, soil, IO_file_name, isStandalone, isRestart, useAlbedo, &
       useDynveg, fileformat, fileztype, IO_diag_stream, IO_accw_stream, IO_stream)

    ! Initialize soil

    ! Called from jsbach_init

    USE mo_jsbach_grid, ONLY: grid_type, domain_type
    USE mo_netcdf,      ONLY: NF_NOERR, nf_inq_varid, nf_get_var_double, &
                              IO_inq_varid, IO_get_var_double, nf_max_name
    USE mo_tr_scatter,  ONLY: scatter_field
    USE mo_temp,        ONLY: zreal2d, zreal3d, zreal2d_ptr, zzreal2d

    TYPE(grid_type),         INTENT(in)    :: grid
    TYPE(domain_type),       INTENT(in)    :: domain
    TYPE(soil_param_type),   INTENT(inout) :: soil_param
    TYPE(soil_type),         INTENT(inout) :: soil
    CHARACTER(nf_max_name),  INTENT(in)    :: IO_file_name
    LOGICAL,                 INTENT(in)    :: isStandalone
    LOGICAL,                 INTENT(in)    :: isRestart
    LOGICAL,                 INTENT(in)    :: UseAlbedo
    LOGICAL,                 INTENT(in)    :: UseDynveg
    INTEGER,                 INTENT(in)    :: fileformat       ! output file format
    INTEGER,                 INTENT(in)    :: fileztype
    TYPE(t_stream), POINTER            :: IO_diag_stream, IO_accw_stream
    TYPE(t_stream), POINTER, OPTIONAL  :: IO_stream

    REAL(dp), PARAMETER :: VolHeatCap (0:6) = (/2.25e+6_dp,1.93e+6_dp,2.10e+6_dp,2.25e+6_dp,2.36e+6_dp,2.48e+6_dp,-1._dp/)
    REAL(dp), PARAMETER :: ThermalDiff(0:6) = (/7.4e-7_dp ,8.7e-7_dp ,8.0e-7_dp ,7.4e-7_dp ,7.1e-7_dp ,6.7e-7_dp ,-1._dp/)

    INTEGER               :: ntiles
    INTEGER               :: IO_file_id, IO_var_id
    INTEGER               :: i, status
    REAL(dp), POINTER     :: fao(:)            !! FAO soil type flag [0...5]
    REAL(dp), ALLOCATABLE :: surface_temperature_clim(:,:)  !! Climatological values of surface temperature


    IF (debug) CALL message('init_soil','')

    CALL config_soil()

    ldiag_soil  = soil_options%ldiag
    ntiles      = soil%ntiles
    soil%nsoil  = soil_options%nsoil
    soil%nsnow  = soil_options%nsnow

    CALL soil_init_io (grid, soil, IO_file_name)
    IO_file_id = soil_file%file_id

    CALL soil_init_memory(grid%nland, domain%nland, ntiles, useDynveg, soil, isStandalone, fileformat, fileztype, &
         IO_diag_stream, IO_accw_stream, stream=IO_stream)

    IF (.NOT. ASSOCIATED(IO_soil)) &
         CALL finish('init_soil', 'No memory stream for soil')

    ! Read global soil parameters

    ! Temporary storage for local domain fields
    ALLOCATE(zreal2d(domain%ndim,domain%nblocks))

    ! Temporary storage for global fields
    IF (p_parallel_io) ALLOCATE(zzreal2d(grid%nlon,grid%nlat))
    IF (p_parallel_io) ALLOCATE( zreal3d(grid%nlon,grid%nlat,nsoil))

    ! Field capacity
    ALLOCATE(soil_param%MaxMoisture(domain%nland))
    IF (p_parallel_io) THEN
       CALL IO_inq_varid(IO_file_id, 'maxmoist', IO_var_id)
       CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
    END IF
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_field(zreal2d_ptr, zreal2d)
    soil_param%MaxMoisture(:) = PACK(zreal2d, MASK=domain%mask)

    IF (soil_options%MaxMoistureLimit > 0.0_dp) THEN
       soil_param%MaxMoisture(:) = MIN(soil_param%MaxMoisture(:), soil_options%MaxMoistureLimit)
    END IF

    ! Initial soil moisture
    ALLOCATE(soil_param%InitMoisture(domain%nland))
    IF (p_parallel_io) THEN
       CALL IO_inq_varid(IO_file_id, 'init_moist', IO_var_id)
       CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
    ENDIF
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_field(zreal2d_ptr, zreal2d)
    soil_param%InitMoisture(:) = PACK(zreal2d, MASK=domain%mask)


    IF (nsoil > 1) THEN

      ! Initialization of soil moisture layers
      ALLOCATE(soil_param%InitLayerMoisture(domain%nland,nsoil))
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'layer_moist', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
      ENDIF
      NULLIFY(zreal2d_ptr)
      DO i=1,nsoil
         IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
         CALL scatter_field(zreal2d_ptr, zreal2d)
         soil_param%InitLayerMoisture(:,i) = PACK(zreal2d, MASK=domain%mask)
      ENDDO

      !(nland) variables for multi-layer scheme

      ! Saturated hydraulic conductivity [m/s]
      ALLOCATE(soil_param%hyd_cond_sat(domain%nland))
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'hyd_cond_sat', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      ENDIF
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      CALL scatter_field(zreal2d_ptr, zreal2d)
      soil_param%hyd_cond_sat = PACK(zreal2d, MASK=domain%mask)

      ! Matrix Potential [m]
      ALLOCATE(soil_param%PotMoisture(domain%nland))
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'moisture_pot', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      ENDIF
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      CALL scatter_field(zreal2d_ptr, zreal2d)
      soil_param%PotMoisture = PACK(zreal2d, MASK=domain%mask)

      ! exponent b in Clapp and Hornberger
      ALLOCATE(soil_param%bclapp(domain%nland))
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'bclapp', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      ENDIF
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      CALL scatter_field(zreal2d_ptr, zreal2d)
      soil_param%bclapp = PACK(zreal2d, MASK=domain%mask)

      ! Volumetric soil porosity [m/m]
      ALLOCATE(soil_param%SoilPorosity(domain%nland))
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'soil_porosity', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      ENDIF
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      CALL scatter_field(zreal2d_ptr, zreal2d)
      soil_param%SoilPorosity = PACK(zreal2d, MASK=domain%mask)

      ! Volumetric soil field capacity [m/m]
      ALLOCATE(soil_param%FieldCapacity(domain%nland))
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'soil_field_cap', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      ENDIF
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      CALL scatter_field(zreal2d_ptr, zreal2d)
      soil_param%FieldCapacity = PACK(zreal2d, MASK=domain%mask)

      ! Volumetric soil permanent wilting point [m/m]
      ALLOCATE(soil_param%WiltingPoint(domain%nland))
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'wilting_point', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      ENDIF
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      CALL scatter_field(zreal2d_ptr, zreal2d)
      soil_param%WiltingPoint = PACK(zreal2d, MASK=domain%mask)

      ! Soil Pore size distribution index
      ALLOCATE(soil_param%PoreSizeIndex(domain%nland))
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'pore_size_index', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      ENDIF
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      CALL scatter_field(zreal2d_ptr, zreal2d)
      soil_param%PoreSizeIndex = PACK(zreal2d, MASK=domain%mask)

      ! Rooting Depth in [m] 
      ALLOCATE(soil_param%RootDepth(domain%nland))
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'root_depth', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      ENDIF
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      CALL scatter_field(zreal2d_ptr, zreal2d)
      soil_param%RootDepth = PACK(zreal2d, MASK=domain%mask)

      ! Soil Depth in [m]
      ALLOCATE(soil_param%SoilDepth(domain%nland))
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'soil_depth', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      ENDIF
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      CALL scatter_field(zreal2d_ptr, zreal2d)
      soil_param%SoilDepth = PACK(zreal2d, MASK=domain%mask)

      soil_param%RootDepth = MIN(soil_param%RootDepth,soil_param%SoilDepth)
    ENDIF

    IF (p_parallel_io) DEALLOCATE(zreal3d)

    ! (nland) variables

    IF (soil_options%withPermafrost) THEN
       ! Heat conductivity of dry soil [J/(m*s*K)]
       ALLOCATE(soil_param%HeatConductivity(domain%nland))
       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'heat_conductivity', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
       CALL scatter_field(zreal2d_ptr, zreal2d)
       soil_param%HeatConductivity = PACK(zreal2d, MASK=domain%mask)
    END IF

    ! Derive heat capacity and conductivity from FAO soil types
    ! Get FAO soil types
    ALLOCATE(fao(domain%nland))
    IF (p_parallel_io) THEN
       CALL io_inq_varid(IO_file_id, 'fao', IO_var_id)
       CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
    END IF
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_field(zreal2d_ptr, zreal2d)
    fao = PACK(zreal2d, MASK=domain%mask)

    ALLOCATE(soil_param%ThermalDiffusivity(domain%nland))
    fao(:) = MAX(REAL(LBOUND(VolHeatCap,1),dp),MIN(REAL(UBOUND(VolHeatCap,1),dp),fao(:)))
    soil_param%ThermalDiffusivity(:) = ThermalDiff(NINT(fao(:)))

    ! Volumetric heat capacity of dry soil [J/(m^3*K)]
    ALLOCATE(soil_param%VolHeatCapacity(domain%nland))
    IF (soil_options%HeatCapMap) THEN
       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'heat_capacity', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
       CALL scatter_field(zreal2d_ptr, zreal2d)
       soil_param%VolHeatCapacity(:) = PACK(zreal2d, MASK=domain%mask)
    ELSE
       soil_param%VolHeatCapacity(:) = VolHeatCap(NINT(fao(:)))
       IF (ANY(soil_param%VolHeatCapacity < 0.0_dp)) &
            CALL finish('init_soil','FAO soil types must be between 0 and 5')
    END IF
    DEALLOCATE(fao)

    ! Surface roughness length
    ALLOCATE(soil_param%Roughness(domain%nland))
    IF (p_parallel_io) THEN
       CALL IO_inq_varid(IO_file_id, 'roughness_length_oro', IO_var_id)
       CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
    ENDIF
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_field(zreal2d_ptr, zreal2d)
    soil_param%Roughness = PACK(zreal2d, MASK=domain%mask)

    ! Fraction of maximum subsurface flow rate
    ALLOCATE(soil_param%Ds(domain%nland))

    ! Average temperature at damping depth
    ALLOCATE(soil_param%AvgTemp(domain%nland))

    DEALLOCATE(zreal2d)

    soil_diag%surface_temperature              = 0._dp
    soil_diag%surface_temperature_acc          = 0._dp
    soil_diag%surface_radiative_temp           = 0._dp
    soil_diag%sat_surface_specific_humidity    = 0._dp
    soil_diag%ground_heat_flux_acc             = 0._dp
    soil_diag%ground_heat_flux_accw            = 0._dp
    soil_diag%heat_capacity                    = 0._dp
    soil_diag%evapotranspiration_acc           = 0._dp
    soil_diag%transpiration_acc                = 0._dp
    soil_diag%sensible_heat_acc                = 0._dp
    soil_diag%latent_heat_acc                  = 0._dp
    soil_diag%albedo                           = 0._dp
    soil_diag%moisture                         = 0._dp
    soil_diag%skin_reservoir                   = 0._dp
    soil_diag%snow                             = 0._dp
    soil_diag%snow_accum_acc                   = 0._dp
    soil_diag%snow_accum_accw                  = 0._dp
    soil_diag%snow_fract                       = 0._dp
    soil_diag%snow_melt                        = 0._dp
    soil_diag%snow_melt_acc                    = 0._dp
    soil_diag%snow_melt_accw                   = 0._dp
    soil_diag%glacier_melt_accw                = 0._dp
    soil_diag%runoff_acc                       = 0._dp
    soil_diag%drainage_acc                     = 0._dp
    soil_diag%runoff_accw                      = 0._dp
    soil_diag%drainage_accw                    = 0._dp
    soil_diag%glacier_depth                    = 0._dp
    soil_diag%glacier_precip_minus_evap_acc    = 0._dp
    soil_diag%glacier_precip_minus_evap_accw   = 0._dp
    soil_diag%glacier_runoff_acc               = 0._dp
    soil_diag%glacier_runoff_accw              = 0._dp
    soil_diag%soil_temperature                 = 0._dp
    soil_diag%net_radiation                    = 0._dp
    soil_diag%qair                             = 0._dp
    soil_diag%qair_acc                         = 0._dp
    soil_diag%wind_lal_acc                     = 0._dp

    IF (ldiag_soil) THEN
      soil_diag%water_balance                  = 0._dp
      soil_diag%storage_pre                    = 0._dp
      soil_diag%bare_soil_evap_acc             = 0._dp
      soil_diag%snow_evap_acc                  = 0._dp
      soil_diag%skin_reservoir_evap_acc        = 0._dp
      soil_diag%evap_deficit_acc               = 0._dp
      soil_diag%precip_acc                     = 0._dp
      IF (nsoil > 1) soil_diag%reduced_evap_acc = 0._dp
      soil_diag%wetskin_fract                  = 0._dp
      soil_diag%csat                           = 0._dp
      soil_diag%cair                           = 0._dp
      soil_diag%csat_transpiration             = 0._dp
    ENDIF
    IF (soil_options%withPermafrost) THEN
      soil_diag%snow_temperature                 = 0._dp
      soil_diag%ice_content                      = 0._dp
      soil_diag%heatcap                          = 0._dp
      soil_diag%heatcond                         = 0._dp
      soil_diag%nsnow_layer                      = 0._dp
      soil_diag%snow_conductivity                = 0.1_dp
      soil_diag%snow_capacity                    = 526500._dp
    ENDIF

    ALLOCATE(zreal2d(domain%ndim,domain%nblocks))

    IF (.NOT. UseAlbedo) THEN
    ! background albedo
      IF (p_parallel_io) THEN
         status = nf_inq_varid(IO_file_id,'albedo', IO_var_id)
         IF (status /= NF_NOERR) THEN
            CALL message('init_soil', 'Initializing background albedo')
            zzreal2d(:,:) = 0.0_dp
         ELSE
            status = nf_get_var_double(IO_file_id, IO_var_id, zzreal2d)
         ENDIF
      ENDIF
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      CALL scatter_field(zreal2d_ptr, zreal2d)
      DO i=1,ntiles
         soil%albedo(:,i) = PACK(zreal2d, MASK=domain%mask)   
      END DO
    ELSE
    ! albedo of vegetation and soil for NIR and visible range
      IF (p_parallel_io) THEN
         status = nf_inq_varid(IO_file_id,'albedo_soil_vis', IO_var_id)
         IF (status /= NF_NOERR) THEN
            CALL message('init_soil', 'Initializing albedo soil visible')
            zzreal2d(:,:) = 0.0_dp
         ELSE
            status = nf_get_var_double(IO_file_id, IO_var_id, zzreal2d)
         ENDIF
      ENDIF
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      CALL scatter_field(zreal2d_ptr, zreal2d)
      DO i=1,ntiles
         soil%albedo_soil_vis(:,i) = PACK(zreal2d, MASK=domain%mask)
         soil%albedo(:,i) = soil%albedo_soil_vis(:,i)
      END DO

      IF (p_parallel_io) THEN
         status = nf_inq_varid(IO_file_id,'albedo_soil_nir', IO_var_id)
         IF (status /= NF_NOERR) THEN
            CALL message('init_soil', 'Initializing albedo soil NIR')
            zzreal2d(:,:) = 0.0_dp
         ELSE
            status = nf_get_var_double(IO_file_id, IO_var_id, zzreal2d)
         ENDIF
      ENDIF
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      CALL scatter_field(zreal2d_ptr, zreal2d)
      DO i=1,ntiles
         soil%albedo_soil_nir(:,i) = PACK(zreal2d, MASK=domain%mask)
      END DO

      IF (p_parallel_io) THEN
         status = nf_inq_varid(IO_file_id,'albedo_veg_vis', IO_var_id)
         IF (status /= NF_NOERR) THEN
            CALL message('init_soil', 'Initializing albedo vegetation visible')
            zzreal2d(:,:) = 0.0_dp
         ELSE
            status = nf_get_var_double(IO_file_id, IO_var_id, zzreal2d)
         ENDIF
      ENDIF
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      CALL scatter_field(zreal2d_ptr, zreal2d)
      DO i=1,ntiles
         soil%albedo_vegetation_vis(:,i) = PACK(zreal2d, MASK=domain%mask)   
      END DO

      IF (p_parallel_io) THEN
         status = nf_inq_varid(IO_file_id,'albedo_veg_nir', IO_var_id)
         IF (status /= NF_NOERR) THEN
            CALL message('init_soil', 'Initializing albedo vegetation NIR')
            zzreal2d(:,:) = 0.0_dp
         ELSE
            status = nf_get_var_double(IO_file_id, IO_var_id, zzreal2d)
         ENDIF
      ENDIF
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      CALL scatter_field(zreal2d_ptr, zreal2d)
      DO i=1,ntiles
         soil%albedo_vegetation_nir(:,i) = PACK(zreal2d, MASK=domain%mask)   
      END DO
    END IF

    ! If this is a restart run we can exit now since the model state variables are read from restart file
    IF (isRestart) THEN
       DEALLOCATE(zreal2d)
       IF (p_parallel_io) DEALLOCATE(zzreal2d)
       RETURN
    ENDIF

    ! If this is not a restart run initialize soil model state from values read from initialization file

    ! Set initial value for root zone soil moisture
    DO i=1,ntiles
       soil%moisture(:,i) = soil_param%InitMoisture(:)
    END DO

    IF (nsoil > 1) THEN
    ! Set initial value for the soil moisture layers in the multi-layer scheme
      soil%layer_moisture(:,:) = soil_param%InitLayerMoisture(:,:)
      WRITE(message_text,*) nsoil,' soil layers layer_moisture are initialized from initial file'
      CALL message('init_soil', message_text)
    ENDIF

    IF (soil_options%lbsoil) THEN
      ! Set initial value for bare soil soil moisture
      IF (nsoil==1) THEN
        soil%bare_soil_moisture(:) =  MAX(0.0_dp, soil_param%InitMoisture(:) -        &
           (soil_param%MaxMoisture(:) - MIN(0.1_dp, soil_param%MaxMoisture(:))) )
      ELSE ! Multilayer soil: Assume bare soil moisture corresponds to that of the upper layer
        soil%bare_soil_moisture(:) = soil%layer_moisture(:,1)
      ENDIF    
    ENDIF

    ! Initial snow depth
    IF (p_parallel_io) THEN
       status = nf_inq_varid(IO_file_id, 'snow', IO_var_id)
       IF (status /= NF_NOERR) THEN
          CALL message('init_soil', 'Initializing snow to zero')
          zzreal2d(:,:) = 0.0_dp
       ELSE
          status = nf_get_var_double(IO_file_id, IO_var_id, zzreal2d)
       ENDIF
    ENDIF
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_field(zreal2d_ptr, zreal2d)
    DO i=1,ntiles
       soil%snow(:,i) = PACK(zreal2d, MASK=domain%mask)
    END DO

    ! Climatological values of surface temperature
    ALLOCATE (surface_temperature_clim(domain%nland,12))

    IF (p_parallel_io) THEN
       ALLOCATE(zreal3d(grid%nlon,grid%nlat,12))
       CALL IO_inq_varid(IO_file_id, 'surf_temp', IO_var_id)
       CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
    ENDIF
    NULLIFY(zreal2d_ptr)
    DO i=1,12
       IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
       CALL scatter_field(zreal2d_ptr, zreal2d)
       surface_temperature_clim(:,i) = PACK(zreal2d, MASK=domain%mask)
    ENDDO
    IF (p_parallel_io) DEALLOCATE(zreal3d)

    ! Initialize soil and surface temperature
    soil%skin_reservoir = 0._dp

    CALL ini_soil_temp(domain%nland, ntsoil, &
         soil%cmid(:),                       &
         surface_temperature_clim(:,:),      &
         soil%soil_temperature(:,:,1),       &
         soil%surface_temperature(:,1))

    DEALLOCATE (surface_temperature_clim)

    IF (ntiles > 1) THEN
      DO i=2,ntiles
        soil%soil_temperature(:,:,i)  = soil%soil_temperature(:,:,1)
        soil%surface_temperature(:,i) = soil%surface_temperature(:,1)
      END DO
    ENDIF

    soil%surface_temperature_old          = soil%surface_temperature
    soil%surface_temperature_unfiltered   = soil%surface_temperature

    soil%latent_heat_flux                 = 0._dp
    soil%latent_heat_acc                  = 0._dp
    soil%sensible_heat_flux               = 0._dp
    soil%sensible_heat_acc                = 0._dp

    soil%heat_capacity                    = 0._dp
    soil%ground_heat_flux                 = 0._dp
    soil%ground_heat_flux_acc             = 0._dp

    soil%snow_accum_acc                   = 0._dp
    soil%snow_melt                        = 0._dp
    soil%snow_melt_acc                    = 0._dp
    soil%snow_age                         = 0._dp
    soil%glacier_depth                    = 0._dp
    soil%glacier_precip_minus_evap_acc    = 0._dp
    soil%runoff_acc                       = 0._dp
    soil%drainage_acc                     = 0._dp
    soil%glacier_runoff_acc               = 0._dp
    soil%qair                             = 0._dp
    soil%qair_acc                         = 0._dp
    soil%z0m                              = 0._dp
    soil%z0h                              = 0._dp
    soil%zril_old                         = 0._dp
    IF (ldiag_soil) THEN
      soil%water_balance                  = 0._dp
      soil%storage_pre                    = 0._dp
      soil%bare_soil_evap_acc             = 0._dp
      soil%snow_evap_acc                  = 0._dp
      soil%skin_reservoir_evap_acc        = 0._dp
      soil%evap_deficit_acc               = 0._dp
      soil%precip_acc                     = 0._dp
      IF (nsoil > 1) soil%reduced_evap_acc = 0._dp
    ENDIF

    soil%ice_content         (:,:,:)      = 0._dp
    IF (soil_options%withPermafrost) THEN
      soil%heatcap           (:,:,:)      = 0._dp
      soil%heatcond          (:,:,:)      = 0._dp
      soil%nsnow_layer       (:,:)        = 0._dp
      soil%k_snow            (:,:)        = 0.1_dp
      soil%c_snow            (:,:)        = 526500._dp
      soil%snow_density      (:,:)        = 100._dp
      soil%psnow_old         (:,:)        = 0._dp
      soil%c_snow_temperature(:,:,:)      = 0._dp
      soil%d_snow_temperature(:,:,:)      = 0._dp
      soil%snow_temperature  (:,:,:)      = SPREAD(soil%surface_temperature(:,:),DIM=2, NCOPIES=soil_options%nsnow)
    ENDIF

    IF (useDynveg) soil%relative_humidity_air = 0._dp

    IF (p_parallel_io) DEALLOCATE(zzreal2d)
    DEALLOCATE(zreal2d)

  END SUBROUTINE init_soil

  SUBROUTINE update_soil(nidx, surface, soil, soil_param, useDynveg, isStandalone, &
       canopy_conductance_max, cover_fract_old,                      &
       lai, canopy_snow, canopy_snow_fract,                          &
       p_echam_cdrag, t_Acoef, t_Bcoef, q_Acoef, q_Bcoef, air_temperature, air_moisture,        &
       surface_pressure, windspeed, wind10, rad_longwave_down, rad_shortwave_net, rain, snow, &
       cair, csat, p_echam_zchl, tte_corr_avg, canopy_conductance_limited,        &
       glac_runoff_evap, surf_runoff_hd, drainage_hd, &
       HeightWind, HeightHumidity, coef_ril_tm1, coef_ril_t, coef_ril_tp1, z0m, z0h)

    USE mo_land_surface, ONLY: land_surface_type
    USE mo_jsbach_constants, ONLY : Gravity, RhoH2O, SpecificHeatDryAirConstPressure, &
         SpecificHeatVaporConstPressure, Emissivity, StefanBoltzmann,                 &
         LatentHeatVaporization, LatentHeatSublimation, GasConstantDryAir,            &
         GasConstantWaterVapor, zsigfac, zepsec, zqsncr, tmelt, sn_depth_min                                   
    USE mo_atmosphere,   ONLY: sat_specific_humidity
    USE mo_utils,        ONLY: average_tiles

    USE mo_physc2,         ONLY: cvdifts       ! Factor for time step weighting in ECHAM5
    USE mo_time_control,   ONLY: lstart, delta_time, time_step_len
    USE mo_semi_impl,      ONLY: eps
    USE mo_param_switches, ONLY: lsurf
    USE mo_exception,      ONLY: message
    USE mo_diag_global_wb, ONLY: diag_global_wb

    INTEGER,                  INTENT(IN)    :: nidx
    TYPE(land_surface_type),  INTENT(IN)    :: surface
    TYPE(soil_type),          INTENT(INOUT) :: soil
    TYPE(soil_param_type),    INTENT(IN)    :: soil_param
    LOGICAL,                  INTENT(IN)    :: useDynveg
    LOGICAL,                  INTENT(IN)    :: isStandalone
    REAL(dp), DIMENSION(:,:), INTENT(IN)    :: &  ! Dimension (nidx,ntiles)
         lai
    REAL(dp), DIMENSION(:,:), INTENT(IN)    :: &  ! Dimension (nidx,ntiles)
         canopy_conductance_max,               &  ! Unstressed canopy resistance
         cover_fract_old                          ! cover fraction of previous timestep
    REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: &  ! Dimension (nidx,ntiles)
         canopy_snow,                          &  ! Snow depth in canopy
         canopy_snow_fract
    REAL(dp), DIMENSION(:,:), INTENT(OUT)   :: &  ! Dimension (nidx,ntiles)
         canopy_conductance_limited               ! resistance limited by water avaiability (water stress)
    REAL(dp), DIMENSION(:),   INTENT(IN)    :: &  ! Dimension nidx
         p_echam_cdrag,                        &
         t_Acoef, t_Bcoef,                     &
         q_Acoef, q_Bcoef,                     &
         air_temperature,                      &  ! Temperature at lowest atmospheric level
         air_moisture,                         &
         surface_pressure,                     &  ! Surface pressure
         windspeed,                            &  ! Wind speed
         wind10,                               &  ! 10m wind speed
         rad_longwave_down,                    &  ! Longwave radiation down
         rad_shortwave_net,                    &  ! Net shortwave radiation
         rain, snow                               ! Rain and snow fall [kg/(m**2s)]


    REAL(dp), DIMENSION(:),   INTENT(IN)    :: p_echam_zchl

    REAL(dp), DIMENSION(:),   INTENT(OUT)   :: &  ! Dimension nidx
         cair, csat, tte_corr_avg, &
         glac_runoff_evap, surf_runoff_hd, drainage_hd ! INPUT for HD-Model in coupled ocean model

    REAL(dp), OPTIONAL,                 INTENT(IN)   :: HeightWind, HeightHumidity, & ! INPUT for offline model
                                                        coef_ril_tm1, coef_ril_t, coef_ril_tp1
    REAL(dp), OPTIONAL, DIMENSION(:),   INTENT(IN)   :: z0m, z0h                  ! INPUT for offline model

! Local variables
    REAL(dp), PARAMETER :: forest_fract_min_org_layer = 0.33_dp   ! Minimum forest frection required organic layer to be present


    REAL(dp), DIMENSION(nidx) ::               &
         zchl,                                 & ! zchl
         cdrag,                                & ! cdrag
         surface_temperature_old,              & ! Averaged surface temperature to calculate Richardson number  
         surface_temperature_upd,              & ! Averaged surface temperature to calculate Richardson number  
         net_radiation,                        & ! Net radiation at surface
         dry_static_energy,                    & ! Surface dry static energy (= C_p * T )
         dry_static_energy_new,                & ! New dry static energy
         surface_qsat_new,                     & ! New surface saturated specific humidity
         air_qsat,                             & ! Saturated specific humidity at lowest atmospheric level
         zdqsl,                                & ! Sensitivity of saturated surface specific humidity wrt temperature
         zcpq,                                 & ! Conversion factor for humidity from dry static energy
         snow_avg,                             & ! Snow [m water equivalent] averaged over all tiles
         nsnow_layer_avg,                      & ! average number of snow layers on the different tiles
         evapotranspiration_soil_avg,          & ! Evapotranspiration without that from snow and the skin reservoir (box avg)
         csat_transpiration                      ! fraction of grid box that contributes to transpiration
                                                 ! (considered to be completely wet)
    REAL(dp), DIMENSION(nidx) ::                                  &
         sat_surf_specific_hum_avg,                               &
         ground_heat_flux_avg, heat_capacity_avg,                 &
         surface_temperature_avg, surface_temperature_unfilt_avg, &
         soil_temperature_avg, soil_moisture_avg,                 &
         snow_melt_avg, glacier_melt_avg,                         &
         snow_accum_avg, glacier_runoff_avg, ice_content_avg

    REAL(dp), DIMENSION(nidx,soil%ntiles) ::   &
         skin_reservoir_max,                   &
         skinres_soil_max,                     &
         skinres_canopy_max,                   &
         wet_skin_fract,                       &
         glacier_precip_minus_evap,            & ! P-E for glaciers [m]
         surface_runoff, drainage,             & ! Surface runoff and drainage for HD model [m]
         melt_water_excess,                    & ! water from snow melting which exceeds skin reservoir and infiltrates soil
         snow_melt, glacier_melt,              & ! Snow/ice melt on land points, resp. glacier points
         snow_accum, glacier_runoff,           &
         qsat_fact, qair_fact,                 & ! Factors for implicit coupling
         air_dry_static_energy_new,            &
         air_moisture_new,                     &
         evaporation_pot,                      & ! Potential evaporation
         evapotranspiration_no_snow,           & ! Evapotranspiration without that from snow
         evapotranspiration_soil,              & ! Evapotranspiration without that from snow and the skin reservoir
         canopy_resistance,                    & ! Water-limited canopy resistance
         soil_moisture_root,                   & ! Soil moisture in root zone
         soil_moisture_root_max,               & ! Field capacity in root zone
         water_stress_factor,                  & ! Water stress factor (=1: no stress, =0: infinite stress)
         relative_humidity,                    & ! Relative humidity (Eq. 3.3.2.9 ECHAM3 Manual)
         tte_corr,                             &
         sat_surface_specific_hum_old,         &
         qsat_veg, qair_veg,                   &
         csat_tiles, cair_tiles,               &
         qsat_transpiration, csat_transpiration_tiles
#ifdef __PGI
    REAL(dp), DIMENSION(nidx,soil%ntiles) :: moisture_max
#endif

    LOGICAL, DIMENSION(nidx,soil%ntiles) :: &
         soil_mask                          ! True if not glacier and not lake

    REAL(dp) :: zcons30, ztpfac2, ztpfac3, vtmpc2
    REAL(dp) :: hlp1
    INTEGER  :: ntiles, kidx0, kidx1, itile, isoil, i, j
    ! multi-layer scheme variables and variables for extra water balance diagnostics
    REAL(dp) :: field_capacity(nidx, soil%nsoil, soil%ntiles)  ! Field capacity of soil layer i [m]
    INTEGER  :: jllog                                 ! Switch/Index of water balance test grid box
!    INTEGER  :: jlpoint                              ! Absolute kpoints Index of water balance test grid box
    REAL(dp) :: skin_res_evap_tile(nidx,soil%ntiles)  ! Skin reservoir evaporation from update_surf_up per tile [m]
    REAL(dp) :: skin_reservoir_evap(nidx)             ! Total-Skin reservoir evaporation [m]
    REAL(dp) :: snow_evap(nidx)                       ! Total snow evaporation [kg/m**2/s]
    REAL(dp) :: snow_evap_tile(nidx,soil%ntiles)      ! Evaporation over snow per tile [kg/m**2/s]
    REAL(dp) :: layer_moisture_tile(nidx, soil%nsoil, soil%ntiles)
                                                      ! Temporary field to distribute layer_moisture to different tiles [m]
    REAL(dp) :: reduced_evap(nidx,soil%ntiles)        ! Diagnostic evaporation reduction obtained from multi-layer scheme
                                                      ! --> is currently distributed between the soil layers [kg/m**2/s]
    REAL(dp) :: evap_deficit_tile(nidx,soil%ntiles)   ! Evaporation deficit flux per tile [m]
    REAL(dp) :: evap_deficit_avg(nidx)                ! Evaporation deficit flux average [m]
    REAL(dp) :: dum(nidx)                             ! Array for various temporary data
!
!   The atmosphere sees only the average skin reservoir regarding evapotranspiration. 
!   Thus, tiling makes no sense and leads to errors in Back-calculations of bare soil evap.
    REAL(dp) :: SkinReservoir(nidx)                   ! Mean skin reservoir [m]
    REAL(dp) :: SkinReservoir_tile(nidx,soil%ntiles)  ! Mean skin reservoir distributed on tiles [m]
    REAL(dp) :: WetSkin_frac(nidx)                    ! Mean wet skin fraction
    REAL(dp) :: MaxSkinResCap(nidx)                   ! Mean maximum skin reservoir capacity[m]
!   Fields for calculating changes of the Bare soil moisture
    REAL(dp) :: bare_soil_moisture_pre(nidx)          ! Bare soil moisture previous time step [m]
    REAL(dp) :: soil_moisture_upper_pre(nidx)         ! Upper layer soil moisture previous time step [m]
    REAL(dp) :: bare_soil_evap(nidx)                  ! Bare soil evaporation in current time step [m]
    REAL(dp) :: zdeltaws(nidx)                        ! Array for temporary differences 
    REAL(dp) :: zupper_bucket(nidx)                   ! Temporary for upper 10 cm bucket moisture
    REAL(dp) :: liquid_max(nidx)                      ! Amount of supercooled water 

    LOGICAL  :: org_layer(nidx)                       ! organic layer mask

    EXTERNAL update_surfacetemp

    ntiles = soil%ntiles
    kidx0  = kstart
    kidx1  = kend

    zcons30 = 1._dp / (cvdifts*Gravity*time_step_len)
    ztpfac2 = 1._dp / cvdifts
    ztpfac3 = 1._dp - ztpfac2
    vtmpc2  = SpecificHeatVaporConstPressure / SpecificHeatDryAirConstPressure - 1._dp

    zchl(1:nidx) = p_echam_zchl(1:nidx)
    cdrag(1:nidx) = p_echam_cdrag(1:nidx)
  
    soil_mask = .NOT. surface%is_glacier(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) &
         .AND. surface%is_present(kidx0:kidx1,:)


    ! update the forest fraction based on the cover fraction of tiles that are considered forest
    surface%forest_fract(kidx0:kidx1) = 0._dp
    DO itile = 1, ntiles
      WHERE (surface%is_forest(kidx0:kidx1,itile))
        surface%forest_fract(kidx0:kidx1) = surface%forest_fract(kidx0:kidx1)  &
                                          + surface%cover_fract(kidx0:kidx1,itile) * surface%veg_ratio_max(kidx0:kidx1)
      END WHERE
    END DO
    org_layer(:)     = .FALSE.
    IF (soil_options%lorganic) THEN
      WHERE (surface%forest_fract(kidx0:kidx1) > forest_fract_min_org_layer)
        org_layer(:) = .TRUE.
      END WHERE
    END IF


    ! Storing forcing variables on the lowest atmospheric model level
    soil_diag%qair(kidx0:kidx1)        = air_moisture(1:nidx) 
    soil_diag%tair(kidx0:kidx1)        = air_temperature(1:nidx)
    soil_diag%tair_max(kidx0:kidx1)    = MAX(soil_diag%tair_max(kidx0:kidx1), air_temperature(1:nidx))  
    soil_diag%tair_min(kidx0:kidx1)    = MIN(soil_diag%tair_min(kidx0:kidx1), air_temperature(1:nidx))
    soil_diag%wind_lal(kidx0:kidx1)    = windspeed(1:nidx)
    soil_diag%wind_lal_acc(kidx0:kidx1)= soil_diag%wind_lal_acc(kidx0:kidx1) + windspeed(1:nidx) * delta_time


    !Update the water and ice reservoirs after cover fractions have changed 
    !----------------------------------------------------------------------------------------

    IF (ANY(ABS(cover_fract_old(1:nidx,:) - surface%cover_fract(kidx0:kidx1,:)) > 0.0_dp)) THEN
      CALL update_h2o_reservoirs(nidx, ntiles,                                           &
                                 cover_fract_old, surface%cover_fract(kidx0:kidx1,:),    &
                                 canopy_snow(1:nidx,:), soil%skin_reservoir(kidx0:kidx1,:))
    END IF
    !----------------------------------------------------------------------------------------
    IF (lstart) THEN
       csat(1:nidx) = .5_dp
       cair(1:nidx) = .5_dp
       csat_transpiration(1:nidx) = .5_dp
    ELSE
       csat(1:nidx) = soil%csat(kidx0:kidx1)
       cair(1:nidx) = soil%cair(kidx0:kidx1)
       csat_transpiration(1:nidx) = soil%csat_transpiration(kidx0:kidx1)
    END IF
    !----------------------------------------------------------------------------------------
    ! Soil moisture in root zone
    DO itile=1,ntiles
       soil_moisture_root(:,itile) = soil%moisture(kidx0:kidx1,itile)
       soil_moisture_root_max(:,itile) = soil_param%MaxMoisture(kidx0:kidx1)
    END DO
    soil_moisture_root     = MERGE(soil_moisture_root,     0._dp, surface%is_vegetation(kidx0:kidx1,:))
    soil_moisture_root_max = MERGE(soil_moisture_root_max, 0._dp, surface%is_vegetation(kidx0:kidx1,:))

    !----------------------------------------------------------------------------------------
    ! Upper layer and bare soil moisture of previous time step
    IF (soil_options%lbsoil) THEN
      bare_soil_moisture_pre(1:nidx) = soil%bare_soil_moisture(kidx0:kidx1)
      IF (nsoil > 1) THEN
        soil_moisture_upper_pre(1:nidx) = soil%layer_moisture(kidx0:kidx1,1)
      ELSE
        CALL average_tiles(soil%moisture(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:)     &
              .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:),  &
              soil_moisture_upper_pre(1:nidx) )
      ENDIF
    ENDIF
    !------------------------------------------------------------------------------------------

    ! Surface dry static energy
    !    ... (loop 331 in vdiff)
    DO itile=1,ntiles
       soil%sat_surface_specific_humidity(kidx0:kidx1,itile) = &
            sat_specific_humidity(soil%surface_temperature(kidx0:kidx1,itile), surface_pressure(1:nidx))
       IF (ANY(soil%sat_surface_specific_humidity(kidx0:kidx1,itile) > 0.99_dp*HUGE(1._dp))) THEN
          CALL message('sat_specific_humidity', 'lookup table overflow', all_print=.TRUE., level=em_warn)
       ENDIF
       sat_surface_specific_hum_old(1:nidx,itile) = soil%sat_surface_specific_humidity(kidx0:kidx1,itile)
    END DO

    CALL average_tiles(soil%surface_temperature(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
                       .AND. .NOT. surface%is_lake(kidx0:kidx1,:),                               &
                       surface%cover_fract(kidx0:kidx1,:), surface_temperature_avg(:))
    CALL average_tiles(soil%sat_surface_specific_humidity(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. surface%is_lake(kidx0:kidx1,:),                                                       &
         surface%cover_fract(kidx0:kidx1,:),sat_surf_specific_hum_avg(:))
    !
    ! Modify minimum canopy resistance (no water stress) according to water limitation in the soil root zone 
    canopy_resistance          = 1.e20_dp
    water_stress_factor        = 0._dp
    canopy_conductance_limited = 1._dp/1.e20_dp

    DO itile=1,ntiles
      DO i=1,nidx
        j=kidx0+i-1
        IF (surface%is_vegetation(j,itile)) THEN
           water_stress_factor(i,itile) = calc_water_stress_factor(soil_moisture_root(i,itile), soil_moisture_root_max(i,itile), &
                                                 soil_options%MoistureFractCritical, soil_options%MoistureFractWilting)
           IF (water_stress_factor(i,itile) > EPSILON(1._dp) .AND. canopy_conductance_max(i,itile) > EPSILON(1._dp) .AND. &
                air_moisture(i) .LE. sat_surf_specific_hum_avg(i)) THEN
              canopy_resistance(i,itile) = 1._dp / (canopy_conductance_max(i,itile) * water_stress_factor(i,itile) + 1.e-20_dp)
           ELSE
              canopy_resistance(i,itile) = 1.e20_dp
           END IF
           canopy_conductance_limited(i,itile) = 1._dp/MAX(canopy_resistance(i,itile),1.e-20_dp)
        END IF
      END DO
    END DO

    ! Sensitivity of saturated surface specific humidity to temperature
    zdqsl = (sat_specific_humidity(surface_temperature_avg(:) + 0.001_dp, surface_pressure(1:nidx)) - &
             sat_surf_specific_hum_avg(:)) * 1000._dp

    ! Compute skin reservoirs
    skinres_soil_max(:,:)   =  MERGE(soil_options%SkinReservoirMax, 0._dp, .NOT. surface%is_glacier(kidx0:kidx1,:))
    skinres_canopy_max(:,:) =  MERGE(soil_options%SkinReservoirMax * lai(:,:)                    &
                              * SPREAD(surface%veg_ratio_max(kidx0:kidx1),DIM=2,ncopies=ntiles), &
                              0._dp, .NOT. surface%is_glacier(kidx0:kidx1,:))
    skin_reservoir_max(:,:) = skinres_soil_max(:,:) + skinres_canopy_max(:,:)

    evaporation_pot = 0._dp
    do itile=1,ntiles
      do i=1,nidx
        j=kidx0+i-1
        IF (skin_reservoir_max(i,itile) > EPSILON(1._dp)) THEN
           canopy_snow_fract(i,itile) = MIN(1._dp, canopy_snow(i,itile) / skin_reservoir_max(i,itile))
        ELSE
           canopy_snow_fract(i,itile) = 0._dp
        END IF
        IF (soil_mask(i,itile)) THEN
           wet_skin_fract(i,itile) = MERGE(MIN(1._dp, soil%skin_reservoir(j,itile) / skin_reservoir_max(i,itile)), &
                0.0_dp, soil_options%SkinReservoirMax > EPSILON(1._dp))
          
           ! Roesch et al, 2002, Climate Dynamics
           soil%snow_fract(j,itile) = zqsncr * TANH(soil%snow(j,itile) * 100._dp) *            &
                SQRT(soil%snow(j,itile) * 1000._dp / (soil%snow(j,itile) * 1000._dp + zepsec + &
                zsigfac * surface%oro_std_dev(j)))
           IF (soil%snow_fract(j,itile) < EPSILON(1._dp) .AND. canopy_snow_fract(i,itile) > EPSILON(1._dp)) &
                soil%snow_fract(j,itile) = canopy_snow_fract(i,itile)
           ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than 
           ! equivalent snow water content from soil and canopy; same for skin reservoir
           ! Potential evaporation using old values of air and surface humidity
           evaporation_pot(i,itile) = zcons30 * cdrag(i) * &
                (soil%sat_surface_specific_humidity(j,itile) - air_moisture(i))
           IF (soil%snow_fract(j,itile) > 0._dp) THEN
              soil%snow_fract(j,itile) = soil%snow_fract(j,itile) /                              &
                   MAX(1._dp, soil%snow_fract(j,itile) * evaporation_pot(i,itile) * delta_time / &
                   (RhoH2O*(soil%snow(j,itile) + canopy_snow(i,itile))))
           END IF
           IF (wet_skin_fract(i,itile) > 0._dp) THEN
              wet_skin_fract(i,itile) = wet_skin_fract(i,itile) /                                          &
                   MAX(1._dp, (1._dp - soil%snow_fract(j,itile)) * evaporation_pot(i,itile) * delta_time / &
                   (RhoH2O * MAX(EPSILON(1._dp),soil%skin_reservoir(j,itile))))
           END IF
           
        ELSE IF (surface%is_glacier(j,itile)) THEN
           wet_skin_fract(i,itile)  = 0.0_dp
           soil%snow_fract(j,itile) = 1.0_dp
        ELSE
           wet_skin_fract(i,itile)  = 0.0_dp
           soil%snow_fract(j,itile) = 0.0_dp
        END IF
      END DO
    END DO
    !----------------------------------------------------------------------------------------------------------------------

    ! Note: at the moment, surface_temperature and sat_surface_specific_humidity are the same for all tiles in a grid box
    dry_static_energy(:) = surface_temperature_avg(:) * SpecificHeatDryAirConstPressure * &
         (1._dp+ vtmpc2 * ( csat(:) * sat_surf_specific_hum_avg(:) + (1._dp - cair(:)) * air_moisture(:)))
    
    ! Conversion factor for humidity from dry static energy
    zcpq(:) = dry_static_energy(:) / surface_temperature_avg(:)

    net_radiation(:) = rad_shortwave_net(:) + &
         rad_longwave_down(:) - Emissivity * StefanBoltzmann * surface_temperature_avg(:)**4

    ! Accumulate net radiation for JSBACH offline output
    soil_diag%net_radiation(kidx0:kidx1) = soil_diag%net_radiation(kidx0:kidx1) + net_radiation(:) * delta_time

    ! Compute new surface temperature and moisture
    dry_static_energy_new = 0._dp
    surface_qsat_new      = 0._dp

    CALL average_tiles(soil%ground_heat_flux(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. surface%is_lake(kidx0:kidx1,:),                                          &
         surface%cover_fract(kidx0:kidx1,:),ground_heat_flux_avg(:))
    CALL average_tiles(soil%heat_capacity(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. surface%is_lake(kidx0:kidx1,:) ,                                      &
         surface%cover_fract(kidx0:kidx1,:),heat_capacity_avg(:))

    IF (.NOT. lstart .AND. isStandalone) THEN
       CALL update_surfacetemp(nidx, zcpq(:),                                                &
            t_Acoef(:), t_Bcoef(:), q_Acoef(:), q_Bcoef(:),                                  &
            dry_static_energy(:), sat_surf_specific_hum_avg(:), zdqsl(:),                    &
            net_radiation(:), ground_heat_flux_avg(:),                                       &
            zcons30*cdrag(:), cair(:), csat(:),                                              &
            SUM(surface%cover_fract(kidx0:kidx1,:) * soil%snow_fract(kidx0:kidx1,:), DIM=2), &
            heat_capacity_avg(:),                                                            &
            dry_static_energy_new(:))

!
!      *** Calculate Richardson number nfiltered with updated surface temperature
!!!       surface_temperature_temp(:) = (surface_temperature_unfilt_avg(:) + surface_temperature_avg(:) ) * 0.5_dp
       CALL average_tiles(soil%surface_temperature_old(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
                       .AND. .NOT. surface%is_lake(kidx0:kidx1,:),                               &
                       surface%cover_fract(kidx0:kidx1,:), surface_temperature_old(:))

       ! Unfiltered land surface temperature at t+dt used to update drag coefficient
       surface_temperature_unfilt_avg(:) = (ztpfac2 * dry_static_energy_new(:) + ztpfac3 * dry_static_energy(:)) / zcpq(:)

       surface_temperature_upd(:) = surface_temperature_old(:) + &
                  eps * (surface_temperature_avg(:)                                &
                  - 2._dp * surface_temperature_old(:)                           &
                  + surface_temperature_unfilt_avg(:) ) 

       CALL update_cdrag( &                         
            nidx, &                                           
            HeightWind, HeightHumidity, &
            coef_ril_tm1, coef_ril_t, coef_ril_tp1, &
            air_temperature(:), air_moisture(:), windspeed(:), &
            surface_temperature_avg(:), surface_temperature_upd(:), &
            surface_pressure(:), &
            csat(:), cair(:), &
            z0m(:), z0h(:), &
            soil%zril_old(kidx0:kidx1), &
            cdrag(:), zchl(:))

    ENDIF


    IF (.NOT. lstart) THEN
       CALL update_surfacetemp(nidx, zcpq(:),                                                &
            t_Acoef(:), t_Bcoef(:), q_Acoef(:), q_Bcoef(:),                                  &
            dry_static_energy(:), sat_surf_specific_hum_avg(:), zdqsl(:),                    &
            net_radiation(:), ground_heat_flux_avg(:),                                       &
            zcons30*cdrag(:), cair(:), csat(:),                                              &
            SUM(surface%cover_fract(kidx0:kidx1,:) * soil%snow_fract(kidx0:kidx1,:), DIM=2), &
            heat_capacity_avg(:),                                                            &
            dry_static_energy_new(:))
    ELSE
       dry_static_energy_new(:) = dry_static_energy(:)
    ENDIF

    ! New unfiltered land surface temperature at t+dt (\tilde(X)^(t+1))
    soil%surface_temperature_unfiltered(kidx0:kidx1,:) = &
         SPREAD((ztpfac2 * dry_static_energy_new(:) + ztpfac3 * dry_static_energy(:)) / zcpq(:), NCOPIES=ntiles, DIM=2)

    ! Correction for snowmelt
    CALL average_tiles(soil%snow(kidx0:kidx1,1:ntiles),                                 &
         surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:),  &
         surface%cover_fract(kidx0:kidx1,:), snow_avg(:))
    WHERE (snow_avg(:) > soil_options%CriticalSnowDepth .OR. ANY(surface%is_glacier(kidx0:kidx1,:), DIM=2)) &
         dry_static_energy_new(:) = MIN(dry_static_energy_new(:)/zcpq(:), tmelt) * zcpq(:)

    ! Correction for melting or freezing of ice in the uppermost soil layer, if there is no snow layer above
    IF (soil_options%withPermafrost .AND. soil_options%lfreeze) THEN

       ! Correction for melting of soil ice
       CALL average_tiles(soil%nsnow_layer(kidx0:kidx1,1:ntiles),                           &
            surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) ,  &
            surface%cover_fract(kidx0:kidx1,:), nsnow_layer_avg(:))
       CALL average_tiles(soil%ice_content(kidx0:kidx1,1,1:ntiles),                         &
            surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) ,  &
            surface%cover_fract(kidx0:kidx1,:), ice_content_avg(:))
       WHERE (nsnow_layer_avg(:) == 0._dp .AND. ice_content_avg(:) > soil_options%CriticalSnowDepth)  &
          dry_static_energy_new(:) = MIN(dry_static_energy_new(:)/zcpq(:), tmelt) * zcpq(:)

       ! Correction for freezing of soil water
       liquid_max(:) = get_liquid_max(nidx, soil_options%lsupercool, dry_static_energy_new(:)/zcpq(:),             &
            org_layer(:), soil%cdel(1), soil_param%SoilPorosity(kidx0:kidx1), soil_param%PotMoisture(kidx0:kidx1), &
            soil_param%BCLAPP(kidx0:kidx1))
       WHERE (nsnow_layer_avg(:) == 0._dp &
            .AND. soil%layer_moisture(kidx0:kidx1,1) > MAX(soil_options%CriticalSnowDepth, liquid_max(:)))
         dry_static_energy_new(:) = MAX(dry_static_energy_new(:)/zcpq(:), tmelt) * zcpq(:)
       END WHERE
    END IF

    ! Compute new surface qsat (after correction of dry static energy for snow melt!)
    ! Was previously in update_surfacetemp: pqsnew(:) = pqsold(:) + zicp(:) * pdqsold(:) * (psnew(:) - psold(:))

    IF (.NOT. lstart) THEN
       surface_qsat_new(:) = sat_surf_specific_hum_avg(:) + zdqsl(:) * (dry_static_energy_new(:) - dry_static_energy(:)) / zcpq(:)
    ELSE
       surface_qsat_new(:) = sat_surf_specific_hum_avg(:)
    END IF
    DO itile=1,ntiles
       IF (.NOT. lstart) THEN
           soil%sat_surface_specific_humidity(kidx0:kidx1,itile) = surface_qsat_new(:)
       END IF
       soil%dry_static_energy_new(kidx0:kidx1,itile) = dry_static_energy_new(:)
    END DO

    DO itile=1,ntiles
       ! Compute temperature and moisture at lowest atmospheric level by back-substitution
       air_dry_static_energy_new(:,itile) = t_Acoef(:) * dry_static_energy_new(:) + t_Bcoef(:)
       air_moisture_new(:,itile) = q_Acoef(:) * soil%sat_surface_specific_humidity(kidx0:kidx1,itile) + q_Bcoef(:)
    END DO
    soil%qair(kidx0:kidx1,1:ntiles) = air_moisture_new(1:nidx,1:ntiles)

    DO itile=1,ntiles
       soil%radiative_temperature(kidx0:kidx1,itile) = dry_static_energy_new(:) / zcpq(:)
    END DO

    ! Compute sensible heat flux
    DO itile=1,ntiles
       soil%sensible_heat_flux(kidx0:kidx1,itile) = zcons30 * cdrag(:) * &
                                                    (air_dry_static_energy_new(:,itile) - dry_static_energy_new(:))
    END DO

    ! Compute evaporation and latent heat flux
    DO itile=1,ntiles
       ! Evapotranspiration
       soil%evapotranspiration(kidx0:kidx1,itile) = zcons30 * cdrag(:)  &
       * (cair(:) * air_moisture_new(:,itile) - csat(:)  * soil%sat_surface_specific_humidity(kidx0:kidx1,itile))
       ! Transpiration
       soil%transpiration(kidx0:kidx1,itile) = zcons30 * cdrag(:)  &
       * csat_transpiration(:) * (air_moisture_new(:,itile) - soil%sat_surface_specific_humidity(kidx0:kidx1,itile))
       ! Potential evaporation
       soil%evaporation_pot(kidx0:kidx1,itile) = zcons30 * cdrag(:)  &
       * (air_moisture_new(:,itile) - soil%sat_surface_specific_humidity(kidx0:kidx1,itile))
       ! Evaporation over snow
       IF ((nsoil > 1) .OR. (ldiag_soil)) &
         snow_evap_tile(1:nidx,itile) = soil%snow_fract(kidx0:kidx1,itile) * soil%evaporation_pot(kidx0:kidx1,itile)
    END DO

    ! Latent heat flux
    soil%latent_heat_flux(kidx0:kidx1,:) = LatentHeatVaporization  * soil%evapotranspiration(kidx0:kidx1,:) &
    + (LatentHeatSublimation - LatentHeatVaporization) * soil%snow_fract(kidx0:kidx1,:) * soil%evaporation_pot(kidx0:kidx1,:)
!
!   Initial values for multi-layer scheme
    IF (nsoil > 1 .OR. ldiag_soil) THEN
!      jllog=0
!      DO i=kidx0,kidx1
!        SELECT CASE (kpoints(i))
!           CASE(2314)
!              jllog=i - kidx0 + 1
!        END SELECT
!      ENDDO
!     endif
!   *** Note that the atmosphere does not see tiles, thus the evaporative fluxes are
!          valid for the gridbox average. Therefore, separating the water fluxes per tile is 
!          leading to errors in the water balance, especially for the multi-layer hydrology scheme.
      CALL average_tiles(wet_skin_fract(1:nidx, :),                                       &
           surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
           surface%cover_fract(kidx0:kidx1,:), WetSkin_frac(1:nidx) )
      CALL average_tiles(skin_reservoir_max(1:nidx, :),                                   &
           surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
           surface%cover_fract(kidx0:kidx1,:), MaxSkinResCap(1:nidx) )
    ENDIF

    IF ((nsoil > 1) .OR. (ldiag_soil)) THEN
      CALL average_tiles(soil%skin_reservoir(kidx0:kidx1,:),                              &
           surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
           surface%cover_fract(kidx0:kidx1,:), SkinReservoir(1:nidx) )
!
!   *** Save sum of water storages from previous time step
      IF (nsoil > 1) THEN
        DO itile=1,ntiles
           SkinReservoir_tile(1:nidx, itile) = SkinReservoir(1:nidx)
        ENDDO
      ENDIF
    ENDIF

    !
    !--------------------------------------------------------------------------------------------------------
    ! 
    tte_corr                  = 0._dp
    snow_accum                = 0._dp
    snow_melt                 = 0._dp
    glacier_melt              = 0._dp
    glacier_runoff            = 0._dp
    glacier_precip_minus_evap = 0._dp
    melt_water_excess         = 0._dp
    surface_runoff            = 0._dp
    drainage                  = 0._dp

    DO itile=1,ntiles

       IF (nsoil > 1) SkinReservoir(1:nidx) = soil%skin_reservoir(kidx0:kidx1,itile)

       CALL update_surface_hydrology( nidx,                             &  ! INTENT(in)
            surface%is_present(kidx0:kidx1,itile),                      &
            surface%is_glacier(kidx0:kidx1,itile),                      &
            air_moisture(1:nidx),  wind10(1:nidx),                      &
            air_temperature(1:nidx),                                    &
            skinres_canopy_max(1:nidx,itile),                           &
            skin_reservoir_max(1:nidx,itile),                           &
            soil%snow_fract(kidx0:kidx1,itile),                         &
            soil%heat_capacity(kidx0:kidx1,itile),                      &
            soil%evapotranspiration(kidx0:kidx1,itile),                 &
            soil%evaporation_pot(kidx0:kidx1,itile),                    &
            rain, snow,                                                 &
            soil%surface_temperature_unfiltered(kidx0:kidx1,itile),     &  ! INTENT(inout)
            soil%skin_reservoir(kidx0:kidx1,itile),                     &
            soil%snow(kidx0:kidx1,itile), canopy_snow(1:nidx,itile),    &
            soil%glacier_depth(kidx0:kidx1,itile),                      &
            tte_corr(1:nidx,itile),                                     &
            snow_accum(1:nidx,itile), snow_melt(1:nidx,itile),          &  ! INTENT(out)
            glacier_melt(1:nidx,itile), glacier_runoff(1:nidx,itile),   &
            glacier_precip_minus_evap(:,itile),                         &
            evapotranspiration_no_snow(1:nidx,itile),                   &
            melt_water_excess(1:nidx,itile),                            & 
            evap_deficit_tile(1:nidx, itile) )

       ! update_soil_hydrology is called for each tile, multi-layer hydrology routine is without tiles
       !        as this is incompatible with REMO
       !        --> Dummy Separation into tiles before update_surf_up and tile averaging
       !            afterwards using new temporary field zlayer_moisture_tile

       IF (nsoil > 1) THEN
         layer_moisture_tile(1:nidx, :, itile) = soil%layer_moisture(kidx0:kidx1, :)

         SkinReservoir_tile(1:nidx,    itile) = SkinReservoir_tile(1:nidx, itile) + &
                                                 soil%skin_reservoir(kidx0:kidx1,itile) - SkinReservoir(1:nidx)

         IF (soil_options%withPermafrost) THEN
             CALL update_soil_permafrost(nidx, soil%ntsoil                            &   ! INTENT(in)
                             , soil_options%nsnow, soil_options%lsnow                 &
                             , soil_options%ldynsnow                                  &
                             , soil_options%ldynorg                                   &
                             , soil_options%lfreeze, soil_options%lsupercool          &
                             , soil_options%lheatcap, soil_options%lheatcond          &
                             , soil_options%HcondScheme                               &
                             , surface%is_present(kidx0:kidx1,itile)                  &
                             , surface%is_glacier(kidx0:kidx1,itile)                  &
                             , soil%cdel, soil%cmid                                   &
                             , soil%snow(kidx0:kidx1,itile)                           &
                             , soil_param%ThermalDiffusivity(kidx0:kidx1)             &
                             , soil_param%VolHeatCapacity(kidx0:kidx1)                &
                             , soil_param%HeatConductivity(kidx0:kidx1)               &
                             , soil_param%SoilDepth(kidx0:kidx1)                      &
                             , soil_param%PotMoisture(kidx0:kidx1)                    &
                             , soil_param%BCLAPP(kidx0:kidx1)                         &
                             , soil_param%SoilPorosity(kidx0:kidx1)                   &
                             , soil_param%FieldCapacity(kidx0:kidx1)                  &
                             , org_layer(1:nidx)                                      &
                             , soil%surface_temperature_unfiltered(kidx0:kidx1,itile) &   ! INTENT(inout)
                             , soil%c_soil_temperature(kidx0:kidx1,:,itile)           &
                             , soil%d_soil_temperature(kidx0:kidx1,:,itile)           &
                             , soil%soil_temperature(kidx0:kidx1,:,itile)             &
                             , soil%c_snow_temperature(kidx0:kidx1,:,itile)           &
                             , soil%d_snow_temperature(kidx0:kidx1,:,itile)           &
                             , soil%snow_temperature(kidx0:kidx1,:,itile)             &
                             , soil%ice_content(kidx0:kidx1,:,itile)                  &
                             , layer_moisture_tile(1:nidx, :, itile)                  &
                             , soil%snow_density(kidx0:kidx1,itile)                   &
                             , soil%psnow_old(kidx0:kidx1,itile)                      &
                             , soil%nsnow_layer(kidx0:kidx1,itile)                    &
                             , soil%heatcap(kidx0:kidx1,:,itile)                      &   ! INTENT(out)
                             , soil%heatcond(kidx0:kidx1,:,itile)                     &
                             , soil%c_snow(kidx0:kidx1,itile)                         &
                             , soil%k_snow(kidx0:kidx1,itile)                         &
                             , soil%thaw_depth(kidx0:kidx1,itile)                     &
                             , soil%heat_capacity(kidx0:kidx1,itile)                  &
                             , soil%ground_heat_flux(kidx0:kidx1,itile) )
         ELSE
           CALL update_soil_temperature(nidx, soil%ntsoil                                   &   ! INTENT(in)
                              , surface%is_present(kidx0:kidx1,itile)                       &
                              , surface%is_glacier(kidx0:kidx1,itile)                       &
                              , soil%cdel(:), soil%cmid(:)                                  &
                              , soil%surface_temperature_unfiltered(kidx0:kidx1,itile)      &
                              , soil%snow(kidx0:kidx1,itile)                                &
                              , soil_param%ThermalDiffusivity(kidx0:kidx1)                  &
                              , soil_param%VolHeatCapacity(kidx0:kidx1)                     &
                              , soil%c_soil_temperature(kidx0:kidx1,:,itile)                &   ! INTENT(inout)
                              , soil%d_soil_temperature(kidx0:kidx1,:,itile)                &
                              , soil%soil_temperature(kidx0:kidx1,:,itile)                  &
                              , soil%heat_capacity(kidx0:kidx1,itile)                       &   ! INTENT(out)
                              , soil%ground_heat_flux(kidx0:kidx1,itile) )
         ENDIF

         CALL update_soil_hydrology(nidx, nsoil,            &   ! INTENT(in)
            surface%is_present(kidx0:kidx1,itile),          &
            surface%is_glacier(kidx0:kidx1,itile),          &
            soil%snow_fract(kidx0:kidx1,itile),             &
            WetSkin_frac(1:nidx),                           &
            MaxSkinResCap(1:nidx),                          &
            soil_param%MaxMoisture(kidx0:kidx1),            &
            surface%oro_std_dev(kidx0:kidx1),               &
            soil%evaporation_pot(kidx0:kidx1,itile),        &
            soil%soil_temperature(kidx0:kidx1,1,itile),     &
            rain(1:nidx),                                   &
            melt_water_excess(1:nidx,itile),                &
            evapotranspiration_no_snow(1:nidx,itile),       &
            SkinReservoir_tile(1:nidx, itile),             &   ! INTENT(inout)
            soil%moisture(kidx0:kidx1,itile),               &
            evap_deficit_tile(1:nidx, itile),               &
            evapotranspiration_soil(1:nidx,itile),          &   ! INTENT(out)
            surface_runoff(1:nidx,itile),                   &
            drainage(1:nidx,itile),                         &
            skin_res_evap_tile(1:nidx,itile),               &
            ! multi-layer hydrology fields (optional)
            soil%cdel,                                      &   ! INTENT(in)
            soil_param%RootDepth(kidx0:kidx1),              &
            soil_param%SoilDepth(kidx0:kidx1),              &
            soil_param%hyd_cond_sat(kidx0:kidx1),           &
            soil_param%PotMoisture(kidx0:kidx1),            &
            soil_param%bclapp(kidx0:kidx1),                 &
            soil_param%SoilPorosity(kidx0:kidx1),           &
            soil%transpiration(kidx0:kidx1,itile),          &
            soil_param%FieldCapacity(kidx0:kidx1),          &
            soil_param%WiltingPoint(kidx0:kidx1),           &
            soil_param%PoreSizeIndex(kidx0:kidx1),          &
            org_layer(1:nidx),                              &
            layer_moisture_tile(1:nidx, :, itile),          &   ! INTENT(inout)
            field_capacity(1:nidx, :, itile),               &
            jllog,                                          &
            reduced_evap(1:nidx,itile),                     &   ! INTENT(out)
            ! permafrost fields (optional)
            soil%ice_content(kidx0:kidx1,:,itile))              ! INTENT(inout)

       ELSE ! 1 layer soil model below

         CALL update_soil_temperature(nidx, soil%ntsoil,               &   ! INTENT(in)
              surface%is_present(kidx0:kidx1,itile),                   &
              surface%is_glacier(kidx0:kidx1,itile),                   &
              soil%cdel(:), soil%cmid(:),                              &
              soil%surface_temperature_unfiltered(kidx0:kidx1,itile),  &
              soil%snow(kidx0:kidx1,itile),                            &
              soil_param%ThermalDiffusivity(kidx0:kidx1),              &
              soil_param%VolHeatCapacity(kidx0:kidx1),                 &
              soil%c_soil_temperature(kidx0:kidx1,:,itile),            &   ! INTENT(inout)
              soil%d_soil_temperature(kidx0:kidx1,:,itile),            &
              soil%soil_temperature(kidx0:kidx1,:,itile),              &
              soil%heat_capacity(kidx0:kidx1,itile),                   &   ! INTENT(out)
              soil%ground_heat_flux(kidx0:kidx1,itile) )

         CALL update_soil_hydrology(nidx, nsoil,              &   ! INTENT(in)
              surface%is_present(kidx0:kidx1,itile),          &
              surface%is_glacier(kidx0:kidx1,itile),          &
              soil%snow_fract(kidx0:kidx1,itile),             &
              wet_skin_fract(:,itile),                        &
              skin_reservoir_max(:,itile),                    &
              soil_param%MaxMoisture(kidx0:kidx1),            &
              surface%oro_std_dev(kidx0:kidx1),               &
              soil%evaporation_pot(kidx0:kidx1,itile),        &
              soil%soil_temperature(kidx0:kidx1,1,itile),     &
              rain(1:nidx),                                   &
              melt_water_excess(1:nidx,itile),                &
              evapotranspiration_no_snow(1:nidx,itile),       &
              soil%skin_reservoir(kidx0:kidx1,itile),         &   ! INTENT(inout)
              soil%moisture(kidx0:kidx1,itile),               &
              evap_deficit_tile(1:nidx, itile),               &
              evapotranspiration_soil(1:nidx,itile),          &   ! INTENT(out)
              surface_runoff(1:nidx,itile),                   &
              drainage(1:nidx,itile),                         &
              skin_res_evap_tile(1:nidx,itile))
       ENDIF
    END DO                  ! End of loop over tiles

    ! soil%snow_melt is the combined snow/ice melt on land and glacier points, here in m water equivalent.
    soil%snow_melt(kidx0:kidx1,:) = MERGE(snow_melt(1:nidx,:), glacier_melt(1:nidx,:), .NOT. surface%is_glacier(kidx0:kidx1,:))

    CALL average_tiles(tte_corr, surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
         surface%cover_fract(kidx0:kidx1,:), tte_corr_avg)
    CALL average_tiles(soil%moisture(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:)              &
                      .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                      soil_moisture_avg)
    CALL average_tiles(soil%snow(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:)                    &
                      .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                      snow_avg)
    CALL average_tiles(soil%surface_temperature_unfiltered(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:) &
                      .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:),        &
                      surface_temperature_unfilt_avg)
    IF ((nsoil > 1) .OR. (ldiag_soil)) THEN
      CALL average_tiles(skin_res_evap_tile(1:nidx,:), surface%is_present(kidx0:kidx1,:)                &
                        .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                        skin_reservoir_evap)
      CALL average_tiles(evapotranspiration_soil(1:nidx,:), surface%is_present(kidx0:kidx1,:)   &
                        .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                        evapotranspiration_soil_avg)
      CALL average_tiles(snow_evap_tile(1:nidx,:), surface%is_present(kidx0:kidx1,:)                    &
                        .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                        snow_evap)
    ENDIF
    IF (nsoil > 1) THEN
      DO i=1, soil%nsoil
        CALL average_tiles(layer_moisture_tile(1:nidx,i,:), surface%is_present(kidx0:kidx1,:)           &
                        .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                        soil%layer_moisture(kidx0:kidx1,i) )
      END DO

      ! update the skin_reservoir per tile in a consistent way
      CALL average_tiles(SkinReservoir_tile(1:nidx, :),                                  &
           surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
           surface%cover_fract(kidx0:kidx1,:), SkinReservoir(1:nidx) )
      DO itile=1,ntiles
        DO i=1, nidx
          j=kidx0+i-1
          IF (MaxSkinResCap(i) > EPSILON(1._dp)) THEN
            soil%skin_reservoir(j, itile) = MIN(SkinReservoir(i)/MaxSkinResCap(i), 1._dp) * skin_reservoir_max(i, itile)
          ELSE    
            soil%skin_reservoir(j, itile) = 0._dp
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    ! Test that soil moisture has still a positive value
    IF (MINVAL(soil_moisture_avg(:)) < 0._dp) THEN
      CALL message('update_soil', 'negative soil moisture occurred!', all_print=.TRUE., level=em_warn)
      soil_moisture_avg(:) = MAX(0._dp, soil_moisture_avg(:))
    END IF

    CALL average_tiles(soil%soil_temperature(kidx0:kidx1,1,:), surface%is_present(kidx0:kidx1,:) &
             .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:),     &
             soil_temperature_avg)
    DO itile = 1, ntiles
      soil%moisture(kidx0:kidx1,itile) = soil_moisture_avg(:)
      soil%snow(kidx0:kidx1,itile) = snow_avg(:)
      soil%surface_temperature_unfiltered(kidx0:kidx1,itile) = surface_temperature_unfilt_avg(:)
      ! uppermost soil layer is updated with surface temperature, if no snow is present
      WHERE (snow_avg(:) > sn_depth_min)
        soil%soil_temperature(kidx0:kidx1,1,itile) = soil_temperature_avg
      ELSEWHERE
        soil%soil_temperature(kidx0:kidx1,1,itile) = surface_temperature_unfilt_avg(:)
      END WHERE
    END DO
    DO isoil = 2, soil%ntsoil
      CALL average_tiles(soil%soil_temperature(kidx0:kidx1,isoil,:), surface%is_present (kidx0:kidx1,:) &
                         .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:),&
                         soil_temperature_avg)
      DO itile = 1, ntiles
        soil%soil_temperature(kidx0:kidx1,isoil,itile) = soil_temperature_avg
      END DO
    END DO

    DO isoil=1,soil%ntsoil
       CALL average_tiles(soil%ice_content(kidx0:kidx1,isoil,:), surface%is_present(kidx0:kidx1,:)  &
                    .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                    ice_content_avg)
       DO itile=1,ntiles
          soil%ice_content(kidx0:kidx1,isoil,itile) = ice_content_avg
       END DO
    END DO

    ! Time filter for surface temperature
    IF (lsurf) THEN
       IF (.NOT. lstart) THEN
!          WHERE (surface%is_present(kidx0:kidx1,:))
             soil%surface_temperature(kidx0:kidx1,:) = soil%surface_temperature_old(kidx0:kidx1,:) + &
                  eps * (soil%surface_temperature(kidx0:kidx1,:)                                     &
                  - 2._dp * soil%surface_temperature_old(kidx0:kidx1,:)                              &
                  + soil%surface_temperature_unfiltered(kidx0:kidx1,:))             
             soil%surface_temperature_old(kidx0:kidx1,:) = soil%surface_temperature_unfiltered(kidx0:kidx1,:)
!          END WHERE
       ELSE
          soil%surface_temperature(kidx0:kidx1,:) = soil%surface_temperature_old(kidx0:kidx1,:)
       END IF
    ELSE
       soil%surface_temperature_old(kidx0:kidx1,:) = soil%surface_temperature(kidx0:kidx1,:)
       soil%surface_temperature_unfiltered(kidx0:kidx1,:) = soil%surface_temperature(kidx0:kidx1,:)
    END IF

! HAG {
!       *** Water balance diagnostics for jllog != 0

!    IF (nsoil == 5) THEN
!      jllog = 0
!
!    *** kpoints start with: 183,184,?,265
!!!     DO i=kidx0,kidx1
!!!        WRITE(message_text,*) 'Limits: ', kidx0, kidx1, ' und ', i, kpoints(i), p_pe
!!!        WRITE(message_text,*) 'Limits: ', kidx0, kidx1, ' und ', nidx, ' = ', nidx
!!!        CALL message('Update_soil', message_text)
!!!     ENDDO

!
!    *** kpoints start with: 183,184,?,265
!      DO i=kidx0,kidx1
!        SELECT CASE (grid%kpoints(i))
!          CASE(584)
!            jllog=i
!            jlpoint=584
!!!            WRITE(message_text,*) jllog, jlpoint, p_pe
!!!            CALL message('Update_soil: ', message_text)
!            EXIT
!          CASE(2314)
!            jllog=i
!            jlpoint=2314
!!!            WRITE(*,*) 'Update_soil: ', jllog, jlpoint, p_pe
!            EXIT
!        END SELECT
!      ENDDO
!
!     IF (jllog .NE. 0) THEN
!
!       jllog = jllog - kidx0 + 1
!       CALL DIAG_WATER_BALANCE(nidx, soil%ntsoil, ntiles, &
!          delta_time, surface%is_present(kidx0:kidx1,:),  &
!          surface%is_glacier(kidx0:kidx1,:),              &
!          surface%cover_fract(kidx0:kidx1,:),             &
!          jllog, jlpoint,                                 &
!          rain(1:nidx),                                   &
!          snow(1:nidx),                                   &
!          surface_runoff(1:nidx,:),                       &
!          drainage(1:nidx,:),                             &  
!          soil%evapotranspiration(kidx0:kidx1,:),         &  
!          snow_avg(1:nidx),                               &
!          soil_moisture_avg(1:nidx),                      &
!          soil%skin_reservoir(kidx0:kidx1,:),             &
!          canopy_snow(1:nidx,:),                          &
!          soil%layer_moisture(kidx0:kidx1, :),            &            
!          soil%transpiration(kidx0:kidx1, :),             &  ! [kg/m**2/s]
!          evapotranspiration_soil_avg(1:nidx),            &  ! [m]
!          soil%evaporation_pot(kidx0:kidx1,:),            &
!          skin_reservoire_evap(1:nidx),                   &  ! [m]
!          snow_evap(1:nidx),                              &  ! [kg/m**2/s = mm/s]
!          soil%surface_temperature(kidx0:kidx1,:),        &
!          reduced_evap(1:nidx, :)                         &
!          )
!     ENDIF
!   endif


    ! Calculate Soil water balance for each grid box. The routine diagnoses errors in the water balance.

    IF (ldiag_soil) THEN
      IF (nsoil > 1) THEN

        CALL diag_global_wb(nidx, ntiles, nsoil,         &
          lstart, delta_time,                            &
          surface%is_present(kidx0:kidx1,:),             &
          surface%cover_fract(kidx0:kidx1,:),            &
          rain(1:nidx),                                  &
          snow(1:nidx),                                  &
          surface_runoff(1:nidx,:) * rhoh2o/delta_time,  &   ! [m] -> [kg m-2 s-1]
          drainage(1:nidx,:)       * rhoh2o/delta_time,  &  
          soil%evapotranspiration(kidx0:kidx1,:),        &  
          glacier_runoff(1:nidx,:) * rhoh2o/delta_time,  &
          snow_avg(1:nidx),                              &
          soil%layer_moisture(kidx0:kidx1,1:nsoil),      &
          soil%ice_content(kidx0:kidx1,1:nsoil,:),       &
          soil%skin_reservoir(kidx0:kidx1,:),            &
          canopy_snow(1:nidx,:),                         &
          soil%glacier_depth(kidx0:kidx1,:),             &
          soil%storage_pre(kidx0:kidx1),                 &
          soil%water_balance(kidx0:kidx1)                &
        )

        CALL average_tiles(reduced_evap(1:nidx, :), surface%is_present(kidx0:kidx1,:)       &
            .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
             dum(1:nidx) )
        soil%reduced_evap_acc(kidx0:kidx1) = soil%reduced_evap_acc(kidx0:kidx1) + dum(1:nidx) * delta_time

      ELSE

        CALL diag_global_wb(nidx, ntiles, nsoil,         &
          lstart, delta_time,                            &
          surface%is_present(kidx0:kidx1,:),             &
          surface%cover_fract(kidx0:kidx1,:),            &
          rain(1:nidx),                                  &
          snow(1:nidx),                                  &
          surface_runoff(1:nidx,:) * rhoh2o/delta_time,  &
          drainage(1:nidx,:)       * rhoh2o/delta_time,  &  
          soil%evapotranspiration(kidx0:kidx1,:),        &  
          glacier_runoff(1:nidx,:) * rhoh2o/delta_time,  &
          snow_avg(1:nidx),                              &
          soil_moisture_avg(1:nidx),                     &
          soil%ice_content(kidx0:kidx1,1:nsoil,:),       &
          soil%skin_reservoir(kidx0:kidx1,:),            &
          canopy_snow(1:nidx,:),                         &
          soil%glacier_depth(kidx0:kidx1,:),             &
          soil%storage_pre(kidx0:kidx1),                 &
          soil%water_balance(kidx0:kidx1)                &
       )

      ENDIF ! nsoil==1

      ! Accumulate soil water balance and evaporation fluxes
      CALL average_tiles(evap_deficit_tile(1:nidx, :), surface%is_present(kidx0:kidx1,:) &
          .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:),      &
          evap_deficit_avg(1:nidx) )
      soil%evap_deficit_acc(kidx0:kidx1) = soil%evap_deficit_acc(kidx0:kidx1) + &
                                           evap_deficit_avg(1:nidx) * RhoH2O

      CALL average_tiles(soil%transpiration(kidx0:kidx1, :), surface%is_present(kidx0:kidx1,:) &
          .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:),      &
         dum(1:nidx) )
      soil%bare_soil_evap_acc(kidx0:kidx1) = soil%bare_soil_evap_acc(kidx0:kidx1)   &    
          + (evapotranspiration_soil_avg(1:nidx)-evap_deficit_avg(1:nidx)) * RhoH2O - dum(1:nidx) * delta_time

      soil%snow_evap_acc(kidx0:kidx1) = soil%snow_evap_acc(kidx0:kidx1) + snow_evap(1:nidx) * delta_time
      soil%skin_reservoir_evap_acc(kidx0:kidx1) = soil%skin_reservoir_evap_acc(kidx0:kidx1) + &
                                                   skin_reservoir_evap(1:nidx) * RhoH2O

      ! Store precipitation received by JSBACH
      soil%precip_acc(kidx0:kidx1) = soil%precip_acc(kidx0:kidx1)   &  
                                  + ((rain(1:nidx)+snow(1:nidx)) * delta_time)
    ENDIF ! ldiag_soil

    !--------------------------------------------------------------------------------------------------------
    ! Re-calculate snow and wet skin fractions with updated moisture parameters due to time mismatch with old echam
    !--------------------------------------------------------------------------------------------------------

    DO itile=1,ntiles
      DO i=1,nidx
        j=kidx0+i-1
        IF (soil_mask(i,itile)) THEN
           wet_skin_fract(i,itile) = MERGE(MIN(1._dp, soil%skin_reservoir(j,itile) / skin_reservoir_max(i,itile)), &
                0.0_dp, soil_options%SkinReservoirMax > EPSILON(1._dp))

           ! Roesch et al, 2002, Climate Dynamics
           soil%snow_fract(j,itile) = zqsncr * TANH(soil%snow(j,itile) * 100._dp) *     &
                SQRT(soil%snow(j,itile)*1000. / (soil%snow(j,itile)*1000._dp + zepsec + &
                zsigfac * surface%oro_std_dev(j)))

           IF (soil%snow_fract(j,itile) .LT. EPSILON(1._dp) .AND. canopy_snow_fract(i,itile) .GE. EPSILON(1._dp)) &
                soil%snow_fract(j,itile) = canopy_snow_fract(i,itile)
           ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than
           ! equivalent snow water content from soil and canopy; same for skin reservoir
           ! Potential evaporation using old values of air and surface humidity
           evaporation_pot(i,itile) = zcons30 * cdrag(i) * &
                (sat_surface_specific_hum_old(i,itile) - air_moisture(i))
           IF (soil%snow_fract(j,itile) > 0._dp .AND. soil%snow(j,itile) + canopy_snow(i,itile) > 0._dp) THEN
              soil%snow_fract(j,itile) = soil%snow_fract(j,itile) / &
                   MAX(1._dp, soil%snow_fract(j,itile) * evaporation_pot(i,itile) * delta_time / &
                   (RhoH2O*(soil%snow(j,itile) + canopy_snow(i,itile))))
           END IF
          
           IF (wet_skin_fract(i,itile) > 0._dp) THEN
              wet_skin_fract(i,itile) = wet_skin_fract(i,itile) / &
                   MAX(1._dp, (1._dp - soil%snow_fract(j,itile)) * evaporation_pot(i,itile) * delta_time / &
                   (RhoH2O * MAX(EPSILON(1._dp),soil%skin_reservoir(j,itile))))
           END IF
        ELSE IF (surface%is_glacier(j,itile)) THEN
           wet_skin_fract(i,itile)  = 0.0_dp
           soil%snow_fract(j,itile) = 1.0_dp
        ELSE
           wet_skin_fract(i,itile)  = 0.0_dp
           soil%snow_fract(j,itile) = 0.0_dp
        END IF
      END DO
    END DO

    soil%wetskin_fract(kidx0:kidx1,:) = wet_skin_fract(1:nidx,:)
!
!   *** Calculate changes for separate bare soil moisture storage
    CALL average_tiles(soil%transpiration(kidx0:kidx1, :), surface%is_present(kidx0:kidx1,:) &
        .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:),      &
        dum(1:nidx) )
!
!   *** [evapotranspiration_soil_avg] = m, [soil%transpiration] = [dum] = kg/m^2/s = mm/s
    bare_soil_evap(1:nidx) = evapotranspiration_soil_avg(1:nidx) - dum(1:nidx) * delta_time / RhoH2O

    IF (soil_options%lbsoil) THEN
      IF (nsoil > 1) THEN
        zdeltaws(1:nidx) = soil%layer_moisture(kidx0:kidx1,1) - soil_moisture_upper_pre(1:nidx)
      ELSE   
        zdeltaws(1:nidx) = soil_moisture_avg(1:nidx) - soil_moisture_upper_pre(1:nidx)
      ENDIF
      zdeltaws(1:nidx) = zdeltaws(1:nidx) - bare_soil_evap(1:nidx)
      CALL average_tiles(Surface%veg_ratio(kidx0:kidx1, :), surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:),      &
         dum(1:nidx) )

      WHERE(dum < 1._dp)
        soil%bare_soil_moisture(kidx0:kidx1) = bare_soil_moisture_pre(1:nidx) + zdeltaws(1:nidx) +    &
           bare_soil_evap(1:nidx) / ( 1.0_dp - dum(1:nidx)) 
      ELSEWHERE
!       *** This is not necessarily correct for bucket model, 
!       *** but limitation to available moisture accounts for this
        soil%bare_soil_moisture(kidx0:kidx1) = bare_soil_moisture_pre(1:nidx) + zdeltaws(1:nidx)
      END WHERE
      soil%bare_soil_moisture(kidx0:kidx1) = MAX(soil%bare_soil_moisture(kidx0:kidx1), 0._dp)
      IF (nsoil > 1) THEN
        soil%bare_soil_moisture(kidx0:kidx1) = MIN(soil%bare_soil_moisture(kidx0:kidx1), soil%layer_moisture(kidx0:kidx1,1))
      ELSE
        zupper_bucket(1:nidx) =  MAX(0.0_dp, soil_moisture_avg(1:nidx) -                  &
          (soil_param%MaxMoisture(kidx0:kidx1) - MIN(0.1_dp, soil_param%MaxMoisture(kidx0:kidx1))) )
        soil%bare_soil_moisture(kidx0:kidx1) = MIN(soil%bare_soil_moisture(kidx0:kidx1),   &
           zupper_bucket(1:nidx) )
      ENDIF
!
!     *** Diffusion of water between vegetated and bare soil part of gridbox
      IF (nsoil == 1) then
!
!       *** In the standard bucket scheme, soil type dependent parameters
!       *** are not known --> Diffusion not possible       
!!!       CALL baresoildiff(1, nidx, nidx, delta_time,  &
!!!         zupper_bucket(nidx),                   &
!!!         soil%bare_soil_moisture(kidx0:kidx1),  &
!!!         dum(1:nidx),                           &
!!!         cdel(1),                               &
!!!         soil_param%hyd_cond_sat(kidx0:kidx1),  &            
!!!         soil_param%SoilPorosity(kidx0:kidx1),  &            
!!!         soil_param%bclapp(kidx0:kidx1),        &            
!!!         soil_param%PotMoisture(kidx0:kidx1)   )
      ELSE 
        CALL baresoildiff(1, nidx, nidx, delta_time,  &
         soil%layer_moisture(kidx0:kidx1, 1),   &
         soil%bare_soil_moisture(kidx0:kidx1),  &
         dum(1:nidx),                           &
         soil%cdel(1),                          &
         soil_param%hyd_cond_sat(kidx0:kidx1),  &            
         soil_param%SoilPorosity(kidx0:kidx1),  &            
         soil_param%bclapp(kidx0:kidx1),        &            
         soil_param%PotMoisture(kidx0:kidx1)   )
      ENDIF
    ENDIF

    ! Calculate near surface relative humidity 
    IF (nsoil == 1) then
      ! Calculate relative humidity from water content in upper 10 cm of soil bucket
      IF (soil_options%lbsoil) THEN
        ! Calculate relative humidity from bare soil part water content in upper 10 cm of soil bucket
        relative_humidity = calc_relative_humidity_bsoil(soil%bare_soil_moisture(kidx0:kidx1),  &
                     soil_param%MaxMoisture(kidx0:kidx1), soil_options%MoistureFractWilting, ntiles)

      ELSE
        ! Calculate relative humidity from water content in upper 10 cm of soil bucket
#ifdef __PGI
        moisture_max = SPREAD(soil_param%MaxMoisture(kidx0:kidx1), NCOPIES=ntiles, DIM=2)
        relative_humidity = calc_relative_humidity(soil%moisture(kidx0:kidx1,:),moisture_max, & 
                                                 soil_options%MoistureFractWilting)
#else
        relative_humidity = calc_relative_humidity(soil%moisture(kidx0:kidx1,:), &
           SPREAD(soil_param%MaxMoisture(kidx0:kidx1), NCOPIES=ntiles, DIM=2),   &
           soil_options%MoistureFractWilting)
#endif
      ENDIF
    ELSE
      IF (soil_options%lbsoil) THEN
        ! Calculate relative humidity from water content in bare soil part of uppermost soil layer
        relative_humidity = calc_relative_humidity_upper(soil%bare_soil_moisture(kidx0:kidx1),  &
                                                       field_capacity(1:nidx,1,1:ntiles), ntiles)
      ELSE
        ! Calculate relative humidity from water content in uppermost soil layer
        relative_humidity = calc_relative_humidity_upper(soil%layer_moisture(kidx0:kidx1,1),  &
                                                       field_capacity(1:nidx,1,1:ntiles), ntiles)
      ENDIF
    ENDIF

    !------------------------------------------------------------------------------------------------------

    qsat_fact               = 0._dp
    qair_fact               = 0._dp
    qsat_veg(:,:)           = 0._dp
    qair_veg(:,:)           = 0._dp
    qsat_transpiration(:,:) = 0._dp

    DO itile=1,ntiles
      DO i=1,nidx
        j=kidx0+i-1 
        IF (surface%is_present(j,itile)) THEN
           IF (soil%moisture(j,itile) > & 
                (soil_options%MoistureFractWilting * soil_param%MaxMoisture(j))) THEN
              qsat_veg(i,itile) = soil%snow_fract(j,itile) + (1._dp - soil%snow_fract(j,itile)) *   &
                                  (wet_skin_fract(i,itile) + (1._dp - wet_skin_fract(i,itile)) /    &
                                  (1._dp + zchl(i) * canopy_resistance(i,itile) *                   &
                                  MAX(1.0_dp,windspeed(i))))
              qsat_transpiration(i,itile) = (1._dp - soil%snow_fract(j,itile)) *                    &
                                            (1._dp - wet_skin_fract(i,itile)) /                     &
                                            (1._dp + zchl(i) * canopy_resistance(i,itile) *         &
                                            MAX(1.0_dp,windspeed(i)))
           ELSE
              qsat_veg(i,itile) = soil%snow_fract(j,itile) + (1._dp - soil%snow_fract(j,itile)) * &
                                  wet_skin_fract(i,itile)
              qsat_transpiration(i,itile) = 0._dp
           END IF
           qair_veg(i,itile) = qsat_veg(i,itile)
        END IF
      END DO
    END DO
        
    WHERE (relative_humidity(1:nidx,:) > SPREAD(air_moisture(1:nidx), NCOPIES=ntiles, DIM=2) / &
         SPREAD(sat_surf_specific_hum_avg(1:nidx), NCOPIES=ntiles, DIM=2) .AND.                &
         relative_humidity(1:nidx,:) > 1.e-10_dp )

       qsat_fact = soil%snow_fract(kidx0:kidx1,:) +                                 & 
                   (1._dp - soil%snow_fract(kidx0:kidx1,:)) *                       &
                   (wet_skin_fract(1:nidx,:) + (1._dp - wet_skin_fract(1:nidx,:)) * &
                   relative_humidity(1:nidx,:))
       qair_fact = 1._dp

    ELSEWHERE (surface%is_present(kidx0:kidx1,:))

       qsat_fact = soil%snow_fract(kidx0:kidx1,:) + &
                   (1._dp - soil%snow_fract(kidx0:kidx1,:)) * wet_skin_fract(1:nidx,:)
       qair_fact = qsat_fact

    END WHERE

    WHERE(SPREAD(air_moisture(1:nidx), NCOPIES=ntiles, DIM=2) > &
         SPREAD(sat_surf_specific_hum_avg(1:nidx), NCOPIES=ntiles, DIM=2))
       qsat_fact = 1._dp
       qair_fact = 1._dp
    END WHERE
  
    csat_tiles(1:nidx,:) = Surface%veg_ratio(kidx0:kidx1,:) * qsat_veg(1:nidx,:) +            &
                           (1._dp - Surface%veg_ratio(kidx0:kidx1,:)) * qsat_fact(1:nidx,:)
    cair_tiles(1:nidx,:) = Surface%veg_ratio(kidx0:kidx1,:) * qair_veg(1:nidx,:) +            &
                           (1._dp - Surface%veg_ratio(kidx0:kidx1,:)) * qair_fact(1:nidx,:)
    csat_transpiration_tiles(1:nidx,:) = Surface%veg_ratio(kidx0:kidx1,:) * qsat_transpiration(1:nidx,:)

    ! ECHAM5 compatibility: one surface temperature for whole grid box

    CALL average_tiles(csat_tiles(1:nidx,:), surface%is_present(kidx0:kidx1,:) .AND. .NOT.                &
                       surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), csat(1:nidx))
    CALL average_tiles(cair_tiles(1:nidx,:), surface%is_present(kidx0:kidx1,:) .AND. .NOT.                &
                       surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), cair(1:nidx))
    CALL average_tiles(csat_transpiration_tiles(1:nidx,:), surface%is_present(kidx0:kidx1,:) .AND. .NOT.  &
                       surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), csat_transpiration(1:nidx))

    soil%csat(kidx0:kidx1) = csat(1:nidx)
    soil%cair(kidx0:kidx1) = cair(1:nidx)
    soil%csat_transpiration(kidx0:kidx1) = csat_transpiration(1:nidx)

    !------------------------------------------------------------------------------------------------------------
    ! END RECALC QSAT and QAir
    !------------------------------------------------------------------------------------------------------------

    ! Average and/or accumulate fluxes

    ! INPUT for HD-Model in coupled case:----------
    glac_runoff_evap = 0._dp
    surf_runoff_hd   = 0._dp
    drainage_hd      = 0._dp

    soil%surface_temperature_acc(kidx0:kidx1,:) =  soil%surface_temperature_acc(kidx0:kidx1,:)          &
                                                 + soil%surface_temperature(kidx0:kidx1,:) * delta_time

    CALL average_tiles(surface_runoff(1:nidx,1:ntiles), surface%is_present(kidx0:kidx1,1:ntiles), &
                       surface%cover_fract(kidx0:kidx1,1:ntiles), surf_runoff_hd(1:nidx))
    CALL average_tiles(drainage(1:nidx,1:ntiles), surface%is_present(kidx0:kidx1,1:ntiles), &
                       surface%cover_fract(kidx0:kidx1,1:ntiles), drainage_hd(1:nidx))
    CALL average_tiles(glacier_precip_minus_evap(1:nidx,1:ntiles), surface%is_present(kidx0:kidx1,1:ntiles), & 
                       surface%cover_fract(kidx0:kidx1,1:ntiles), glac_runoff_evap(1:nidx))

    ! Note: water fluxes are in m water equivalent; to convert to kg/m^2s multiply by RhoH2O/delta_time,
    !       and for accumulation this is again multiplied by delta_time, i.e. just multiply by RhoH2O
    WHERE (surface%is_present(kidx0:kidx1,:))
      soil%ground_heat_flux_acc(kidx0:kidx1,:) = soil%ground_heat_flux_acc(kidx0:kidx1,:) + &
                                                   soil%ground_heat_flux(kidx0:kidx1,:) * delta_time
      soil%runoff_acc        (kidx0:kidx1,:) = soil%runoff_acc(kidx0:kidx1,:) + &
                                                 (surface_runoff(1:nidx,:) + drainage(1:nidx,:)) * RhoH2O
      soil%drainage_acc      (kidx0:kidx1,:) = soil%drainage_acc(kidx0:kidx1,:) + drainage(1:nidx,:) * RhoH2O
      soil%glacier_precip_minus_evap_acc(kidx0:kidx1,:) = soil%glacier_precip_minus_evap_acc(kidx0:kidx1,:) + &
                                                           glacier_precip_minus_evap(1:nidx,:) * RhoH2O
      soil%snow_accum_acc     (kidx0:kidx1,:) = soil%snow_accum_acc(kidx0:kidx1,:) + snow_accum(1:nidx,:) * RhoH2O
      soil%snow_melt_acc      (kidx0:kidx1,:) = soil%snow_melt_acc(kidx0:kidx1,:) + soil%snow_melt(kidx0:kidx1,:) * RhoH2O
      soil%glacier_runoff_acc (kidx0:kidx1,:) = soil%glacier_runoff_acc(kidx0:kidx1,:) + glacier_runoff(1:nidx,:) * RhoH2O
    END WHERE


    ! calculate diagnostic variables for accw stream (time average over the output interval and weighted with land fraction)
    !   and other diagnostic variables that are time and tile averages, and had not been calculated above

    CALL average_tiles(soil%ground_heat_flux(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. surface%is_lake(kidx0:kidx1,:),                                          &
         surface%cover_fract(kidx0:kidx1,:), ground_heat_flux_avg(:))
    CALL average_tiles(soil%ground_heat_flux_acc(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. surface%is_lake(kidx0:kidx1,:),                                          &
         surface%cover_fract(kidx0:kidx1,:), soil_diag%ground_heat_flux_acc(kidx0:kidx1))
    soil_diag%ground_heat_flux_accw(kidx0:kidx1) = soil_diag%ground_heat_flux_accw(kidx0:kidx1) + &
                                               surface%land_fract(kidx0:kidx1) * ground_heat_flux_avg(1:nidx) * delta_time

    CALL average_tiles(soil%runoff_acc(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                       soil_diag%runoff_acc(kidx0:kidx1))
    soil_diag%runoff_accw  (kidx0:kidx1)   = soil_diag%runoff_accw(kidx0:kidx1) + surface%land_fract(kidx0:kidx1) * &
                                               (surf_runoff_hd(1:nidx) + drainage_hd(1:nidx)) * RhoH2O

    CALL average_tiles(soil%drainage_acc(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                       soil_diag%drainage_acc(kidx0:kidx1))
    soil_diag%drainage_accw(kidx0:kidx1)   = soil_diag%drainage_accw(kidx0:kidx1) + surface%land_fract(kidx0:kidx1) * &
                                               drainage_hd(1:nidx) * RhoH2O

    CALL average_tiles(soil%glacier_precip_minus_evap_acc(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:), &
                       surface%cover_fract(kidx0:kidx1,:), soil_diag%glacier_precip_minus_evap_acc(kidx0:kidx1))
    soil_diag%glacier_precip_minus_evap_accw(kidx0:kidx1) = soil_diag%glacier_precip_minus_evap_accw(kidx0:kidx1) + &
                                               surface%land_fract(kidx0:kidx1) * glac_runoff_evap(1:nidx) * RhoH2O

    CALL average_tiles(soil%snow_accum_acc(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                       soil_diag%snow_accum_acc(kidx0:kidx1))
    CALL average_tiles(snow_accum(1:nidx,:), surface%is_present(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                       snow_accum_avg(1:nidx))
    soil_diag%snow_accum_accw(kidx0:kidx1)   = soil_diag%snow_accum_accw(kidx0:kidx1) + surface%land_fract(kidx0:kidx1) * &
                                                 snow_accum_avg(1:nidx) * RhoH2O

    CALL average_tiles(soil%snow_melt(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                       soil_diag%snow_melt(kidx0:kidx1))
    CALL average_tiles(soil%snow_melt_acc(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                       soil_diag%snow_melt_acc(kidx0:kidx1))
    snow_melt_avg = 0._dp
    glacier_melt_avg = 0._dp
    CALL average_tiles(snow_melt(1:nidx,:), surface%is_present(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                       snow_melt_avg(1:nidx))
    CALL average_tiles(glacier_melt(1:nidx,:), surface%is_present(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                       glacier_melt_avg(1:nidx))
    soil_diag%snow_melt_accw(kidx0:kidx1)   = soil_diag%snow_melt_accw(kidx0:kidx1) + surface%land_fract(kidx0:kidx1) * &
                                                 snow_melt_avg(1:nidx) * RhoH2O
    soil_diag%glacier_melt_accw(kidx0:kidx1)= soil_diag%glacier_melt_accw(kidx0:kidx1) + surface%land_fract(kidx0:kidx1) * &
                                                 glacier_melt_avg(1:nidx) * RhoH2O

    CALL average_tiles(glacier_runoff(1:nidx,:), surface%is_present(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                       glacier_runoff_avg(1:nidx))
    CALL average_tiles(soil%glacier_runoff_acc(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:), &
                       surface%cover_fract(kidx0:kidx1,:), soil_diag%glacier_runoff_acc(kidx0:kidx1))
    soil_diag%glacier_runoff_accw(kidx0:kidx1) = soil_diag%glacier_runoff_accw(kidx0:kidx1) + surface%land_fract(kidx0:kidx1) * &
                                                   glacier_runoff_avg(1:nidx) * RhoH2O

    ! Calculation of relative humidity (of the air)
    IF (useDynveg) THEN
       air_qsat(1:nidx) = sat_specific_humidity(surface_temperature_unfilt_avg(1:nidx), surface_pressure(1:nidx))
       IF (ANY(air_qsat(1:nidx) > 0.99_dp*HUGE(1._dp))) &
          CALL message('sat_specific_humidity', 'lookup table overflow', all_print=.TRUE., level=em_warn)
       soil%relative_humidity_air(kidx0:kidx1) = 100._dp * air_moisture(1:nidx)     &
         * ((1._dp - GasConstantDryAir / GasConstantWaterVapor)                     &
         * air_qsat(1:nidx) + (GasConstantDryAir / GasConstantWaterVapor))          &
         / (air_qsat(1:nidx) * ((1._dp - GasConstantDryAir / GasConstantWaterVapor) &
         * air_moisture(1:nidx) + (GasConstantDryAir / GasConstantWaterVapor)))
    END IF

    ! Snow age
    DO i=1,nidx
       j=kidx0+i-1
       IF (MAXVAL(soil%snow_fract(j,:)) > EPSILON(1._dp)) THEN
          hlp1 = MIN(1._dp,EXP(5000._dp * ((1._dp / tmelt) - (1._dp / surface_temperature_unfilt_avg(i)))))
          soil%snow_age(j) = MAX(0.0_dp,(soil%snow_age(j) +                               &
                             ((hlp1**10._dp + hlp1 + 0.3_dp) * delta_time * 1.e-06_dp)) * &
                             (1._dp - (snow(i) * delta_time * 0.5_dp))) 
       ELSE
          soil%snow_age(j) = 0._dp
       END IF
    END DO

  END SUBROUTINE update_soil

  ELEMENTAL FUNCTION calc_water_stress_factor(moisture, moisture_max, fract_critical, fract_wilting) RESULT(stress)

    REAL(dp), INTENT(in)  :: moisture
    REAL(dp), INTENT(in)  :: moisture_max
    REAL(dp), INTENT(in)  :: fract_critical, fract_wilting
    REAL(dp)              :: stress

    REAL(dp) :: moisture_critical, moisture_wilting

    moisture_critical = fract_critical * moisture_max
    moisture_wilting  = fract_wilting  * moisture_max

    stress = MAX(0._dp, MIN(1._dp, (moisture - moisture_wilting) / (moisture_critical - moisture_wilting) ))

  END FUNCTION calc_water_stress_factor

  FUNCTION calc_relative_humidity(moisture, moisture_max, fract_wilting) RESULT(rel_humidity)

    USE mo_jsbach_constants,  ONLY: pi

    REAL(dp), INTENT(in)  :: moisture(:,:)
    REAL(dp), INTENT(in)  :: moisture_max(:,:)
    REAL(dp), INTENT(in)  :: fract_wilting
    REAL(dp)              :: rel_humidity(SIZE(moisture,1),SIZE(moisture,2))

    REAL(dp) moisture_limit(SIZE(moisture,1),SIZE(moisture,2)), evap_stop(SIZE(moisture,1),SIZE(moisture,2))

    evap_stop      = MIN(0.1_dp, moisture_max)
    moisture_limit = moisture_max - evap_stop

    WHERE (moisture > moisture_limit .AND. moisture > fract_wilting*moisture_max)
       rel_humidity = 0.5_dp * (1._dp - COS((moisture-moisture_limit) * pi / evap_stop))
    ELSEWHERE
       rel_humidity = 0._dp
    END WHERE

  END FUNCTION calc_relative_humidity

  FUNCTION calc_relative_humidity_bsoil(moisture, moisture_max, fract_wilting, ntiles) RESULT(rel_humidity)

    USE mo_jsbach_constants,  ONLY: pi

    REAL(dp), INTENT(in)  :: moisture(:)
    REAL(dp), INTENT(in)  :: moisture_max(:)
    INTEGER, INTENT(in)   :: ntiles              !! Number of soil tiles
    REAL(dp), INTENT(in)  :: fract_wilting
    REAL(dp)              :: rel_humidity(SIZE(moisture), ntiles)

    INTEGER               :: i
    REAL(dp)              :: evap_stop(SIZE(moisture))
    REAL(dp)              :: rhum(SIZE(moisture))

    evap_stop      = MIN(0.1_dp, moisture_max)

    WHERE (moisture > 0._dp .AND. moisture > fract_wilting*evap_stop)
       rhum = 0.5_dp * (1._dp - COS(moisture * pi / evap_stop))
    ELSEWHERE
       rhum = 0._dp
    END WHERE
    DO i = 1, ntiles
      rel_humidity(:,i) = rhum(:)
    ENDDO

  END FUNCTION calc_relative_humidity_bsoil

  FUNCTION calc_relative_humidity_upper(moisture, moisture_max, ntiles) RESULT(rel_humidity)

    USE mo_jsbach_constants,  ONLY: pi

    REAL(dp), INTENT(in)  :: moisture(:)
    REAL(dp), INTENT(in)  :: moisture_max(:,:)
    INTEGER, INTENT(in)   :: ntiles              !! Number of soil tiles
    REAL(dp)              :: rel_humidity(SIZE(moisture), ntiles)

    INTEGER               :: i
    REAL(dp)              :: moisture_upper(SIZE(moisture))
    REAL(dp)              :: rhum(SIZE(moisture))

    DO i = 1, ntiles
       moisture_upper(:) = DMIN1(moisture(:), moisture_max(:,i))
       WHERE (moisture_upper(:) > 0._dp)
          rhum(:) = 0.5_dp * (1._dp - COS((moisture_upper(:)) * pi / moisture_max(:,i)))
       ELSEWHERE
          rhum(:) = 0._dp
       END WHERE
       rel_humidity(:,i) = rhum(:)
    ENDDO
  
  END FUNCTION calc_relative_humidity_upper

  SUBROUTINE soil_diagnostics(surface, soil)

    USE mo_land_surface, ONLY: land_surface_type
    USE mo_utils,        ONLY: average_tiles
    USE mo_time_control, ONLY: delta_time

    TYPE(land_surface_type), INTENT(in) :: surface
    TYPE(soil_type),         INTENT(inout) :: soil

    !Local variables
    LOGICAL  :: mask(nidx,soil%ntiles)
    REAL(dp) :: fract(nidx,soil%ntiles)
    INTEGER  :: ntiles, i
    INTEGER  :: kidx0, kidx1

    ntiles  = soil%ntiles
    kidx0   = kstart
    kidx1   = kend

    ! Accumulate variables
    WHERE (surface%is_present(kidx0:kidx1,:))
       ! Evapotranspiration, transpiration and latent/sensible heat fluxes
       soil%evapotranspiration_acc(kidx0:kidx1,:) = soil%evapotranspiration_acc(kidx0:kidx1,:) + &
            soil%evapotranspiration(kidx0:kidx1,:) * delta_time
       soil%transpiration_acc(kidx0:kidx1,:) = soil%transpiration_acc(kidx0:kidx1,:) + &
            soil%transpiration(kidx0:kidx1,:) * delta_time
       soil%latent_heat_acc(kidx0:kidx1,:) = soil%latent_heat_acc(kidx0:kidx1,:) + &
            soil%latent_heat_flux(kidx0:kidx1,:) * delta_time
       soil%sensible_heat_acc(kidx0:kidx1,:) = soil%sensible_heat_acc(kidx0:kidx1,:) + &
            soil%sensible_heat_flux(kidx0:kidx1,:) * delta_time
       soil%qair_acc(kidx0:kidx1,:) = soil%qair_acc(kidx0:kidx1,:) + &
            soil%qair(kidx0:kidx1,:) * delta_time
    END WHERE

    ! Compute grid box averages
    mask  = surface%is_present(kidx0:kidx1,1:ntiles)
    fract = surface%cover_fract(kidx0:kidx1,1:ntiles)

    CALL average_tiles(soil%surface_temperature(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%surface_temperature(kidx0:kidx1))

    CALL average_tiles(soil%surface_temperature_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%surface_temperature_acc(kidx0:kidx1))

    CALL average_tiles(soil%radiative_temperature(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%surface_radiative_temp(kidx0:kidx1))
    
    CALL average_tiles(soil%sat_surface_specific_humidity(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%sat_surface_specific_humidity(kidx0:kidx1))

    CALL average_tiles(soil%heat_capacity(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%heat_capacity(kidx0:kidx1))

    CALL average_tiles(soil%evapotranspiration_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%evapotranspiration_acc(kidx0:kidx1))

    CALL average_tiles(soil%transpiration_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%transpiration_acc(kidx0:kidx1))

    CALL average_tiles(soil%sensible_heat_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%sensible_heat_acc(kidx0:kidx1))

    CALL average_tiles(soil%latent_heat_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%latent_heat_acc(kidx0:kidx1))

    CALL average_tiles(soil%albedo(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%albedo(kidx0:kidx1))

    CALL average_tiles(soil%skin_reservoir(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%skin_reservoir(kidx0:kidx1))

    CALL average_tiles(soil%snow(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%snow(kidx0:kidx1))

    CALL average_tiles(soil%snow_fract(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%snow_fract(kidx0:kidx1))

    CALL average_tiles(soil%glacier_depth(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%glacier_depth(kidx0:kidx1))

    CALL average_tiles(soil%moisture(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%moisture(kidx0:kidx1))

    DO i=1,ntsoil
       CALL average_tiles(soil%soil_temperature(kidx0:kidx1,i,1:ntiles), mask, fract, &
            soil_diag%soil_temperature(kidx0:kidx1,i))
    END DO
    IF (soil_options%withPermafrost) THEN
       DO i=1,soil_options%nsnow
          CALL average_tiles(soil%snow_temperature(kidx0:kidx1,i,1:ntiles), mask, fract, &
               soil_diag%snow_temperature(kidx0:kidx1,i))
       END DO
       DO i=1,soil_options%nsoil
          CALL average_tiles(soil%ice_content(kidx0:kidx1,i,1:ntiles), mask, fract, &
               soil_diag%ice_content(kidx0:kidx1,i))
          CALL average_tiles(soil%heatcap(kidx0:kidx1,i,1:ntiles), mask, fract, &
               soil_diag%heatcap(kidx0:kidx1,i))
          CALL average_tiles(soil%heatcond(kidx0:kidx1,i,1:ntiles), mask, fract, &
               soil_diag%heatcond(kidx0:kidx1,i))
       END DO
       CALL average_tiles(soil%nsnow_layer(kidx0:kidx1,1:ntiles), mask, fract, &
               soil_diag%nsnow_layer(kidx0:kidx1))
       CALL average_tiles(soil%snow_density(kidx0:kidx1,1:ntiles), mask, fract, &
               soil_diag%snow_density(kidx0:kidx1))
       CALL average_tiles(soil%thaw_depth(kidx0:kidx1,1:ntiles), mask, fract, &
               soil_diag%thaw_depth(kidx0:kidx1))
       CALL average_tiles(soil%k_snow(kidx0:kidx1,1:ntiles), mask, fract, &
               soil_diag%snow_conductivity(kidx0:kidx1))
       CALL average_tiles(soil%c_snow(kidx0:kidx1,1:ntiles), mask, fract, &
               soil_diag%snow_capacity(kidx0:kidx1))
    ENDIF

    CALL average_tiles(soil%qair_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%qair_acc(kidx0:kidx1))
         
    IF (ldiag_soil) THEN
      IF (nsoil > 1) THEN
        soil_diag%reduced_evap_acc(kidx0:kidx1) = soil%reduced_evap_acc(kidx0:kidx1)
      ENDIF
      soil_diag%water_balance(kidx0:kidx1)   = soil_diag%water_balance(kidx0:kidx1) + soil%water_balance(kidx0:kidx1) * delta_time
      soil_diag%storage_pre(kidx0:kidx1)     = soil_diag%storage_pre(kidx0:kidx1) + soil%storage_pre(kidx0:kidx1) * delta_time
      soil_diag%bare_soil_evap_acc(kidx0:kidx1)       = soil%bare_soil_evap_acc(kidx0:kidx1)
      soil_diag%snow_evap_acc(kidx0:kidx1)            = soil%snow_evap_acc(kidx0:kidx1) 
      soil_diag%skin_reservoir_evap_acc(kidx0:kidx1) = soil%skin_reservoir_evap_acc(kidx0:kidx1) 
      soil_diag%evap_deficit_acc(kidx0:kidx1)         = soil%evap_deficit_acc(kidx0:kidx1) 
      soil_diag%precip_acc(kidx0:kidx1)               = soil%precip_acc(kidx0:kidx1) 
      CALL average_tiles(soil%wetskin_fract(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%wetskin_fract(kidx0:kidx1))
      soil_diag%csat(kidx0:kidx1) = soil%csat(kidx0:kidx1)
      soil_diag%cair(kidx0:kidx1) = soil%cair(kidx0:kidx1)
      soil_diag%csat_transpiration(kidx0:kidx1) = soil%csat_transpiration(kidx0:kidx1)
    ENDIF

  END SUBROUTINE soil_diagnostics

  FUNCTION get_soil_diag(kdim, what) RESULT(outfield)

    ! only 1-d at the moment

    INTEGER, INTENT(in) :: kdim
    CHARACTER(len=*), INTENT(in) :: what
    REAL(dp) :: outfield(kdim)

    INTEGER :: k0, k1

    k0 = kstart
    k1 = kend

    IF (kdim /= (k1-k0+1)) CALL finish('get_soil_diag','Wrong dimension')

    SELECT CASE(TRIM(what))
    CASE('surface_temperature')
       outfield(1:kdim) = soil_diag%surface_temperature(k0:k1)
    CASE('surface_radiative_temperature')
       outfield(1:kdim) = soil_diag%surface_radiative_temp(k0:k1)
    CASE('sat_surface_specific_humidity')
       outfield(1:kdim) = soil_diag%sat_surface_specific_humidity(k0:k1)
    CASE('runoff_acc')
       outfield(1:kdim) = soil_diag%runoff_acc(k0:k1)
    CASE('glacier_depth')
       outfield(1:kdim) = soil_diag%glacier_depth(k0:k1)
    CASE('snow_melt_inst')
       outfield(1:kdim) = soil_diag%snow_melt(k0:k1)
    CASE('glacier_precip_minus_evap_acc')
       outfield(1:kdim) = soil_diag%glacier_precip_minus_evap_acc(k0:k1)
    CASE('water_balance')
       IF (.not. ldiag_soil) CALL finish('get_soil_diag','Parameter only available for soil diagnostics')
       outfield(1:kdim) = soil_diag%water_balance(k0:k1)
    CASE default
       CALL finish('get_soil_diag','Unknown request')
    END SELECT

  END FUNCTION get_soil_diag


  SUBROUTINE update_h2o_reservoirs(nidx, ntiles,                   &
                                 cover_fract_old, cover_fract_new, &
                                 snow_can, skin_res)

    INTEGER, INTENT(IN) :: nidx, ntiles
    REAL(dp), DIMENSION(nidx, ntiles)       , INTENT(IN)        :: & 
            cover_fract_old, cover_fract_new

    REAL(dp), DIMENSION(nidx, ntiles)       , INTENT(INOUT)     :: & 
          snow_can, skin_res

    INTEGER :: itile 
    REAL(dp), DIMENSION(nidx)                                   :: & 
            fract_redis,                                           &
            snow_can_redis, skin_res_redis
    REAL(dp), DIMENSION(nidx, ntiles)                           :: & 
            fract_change, fract_redis_tile

    fract_change(:,:)         = cover_fract_new(:,:) - cover_fract_old(:,:)
    fract_redis_tile(:,:)     = MAX(0.0_dp, fract_change(:,:) * (-1._dp))
    fract_redis(:)            = SUM(fract_redis_tile(:,:),DIM=2)
    snow_can_redis(:)         = 0.0_dp
    skin_res_redis(:)         = 0.0_dp

    DO itile=1,ntiles
      WHERE(fract_change(:,itile) < 0.0_dp)
          snow_can_redis(:)   = snow_can_redis(:) &
                              + snow_can(:,itile) * fract_redis_tile(:,itile) / fract_redis(:)
          skin_res_redis(:)   = skin_res_redis(:) &
                              + skin_res(:,itile) * fract_redis_tile(:,itile) / fract_redis(:)
      END WHERE
    ENDDO
    DO itile=1,ntiles
      WHERE(fract_change(:,itile) > 0.0_dp)
        snow_can(:,itile) = (snow_can(:,itile) * cover_fract_old(:,itile) & 
                          +  snow_can_redis(:) * fract_change(:,itile)) &
                          /  cover_fract_new(:,itile)
        skin_res(:,itile) = (skin_res(:,itile) * cover_fract_old(:,itile) & 
                          +  skin_res_redis(:) * fract_change(:,itile)) &
                          /  cover_fract_new(:,itile)
      END WHERE
    ENDDO

  END SUBROUTINE update_h2o_reservoirs



  SUBROUTINE ini_soil_temp(klon, ntsoil, cmid, tslclim, tsoil, tsl)
    !
    ! Initialize soil and surface temperature 
    !
    ! Method:
    !
    ! Starting from the tslclim field temperatures are set in relation
    ! to the depth of the soil layer and position of the initial
    ! day in the annual cycle.
    !
    ! tsl is at 0.07 m
    ! thickness of layers 0.065, 0.254, 0.913, 2.902, 5.700 m (s. soiltemp)
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, April 2000, original source
    ! L. Dumenil, MPI, June 1989, original source of *initemp*
    !             which is now part of this routine
    ! I. Kirchner, MPI, November 2000, date/time control
    ! U. Schulzweida, MPI, May 2002, blocking (nproma)
    ! E. Roeckner, MPI, Sep 2002, initialization of mixed layer ocean
    ! T. Raddatz, MPI, Mai 2005, adaption to JSBACH
    ! 
    USE mo_jsbach_constants,     ONLY: pi
    USE mo_radiation_parameters, ONLY: nmonth
    USE mo_time_control,         ONLY: get_date_components, ndaylen, start_date, get_month_len, get_year_len

    INTEGER,  INTENT(IN)     ::  klon                 !! number of land grid points
    INTEGER,  INTENT(IN)     ::  ntsoil               !! number of thermal soil layers
    REAL(dp), INTENT(IN)     ::  cmid(ntsoil)         !! depth of soil layer mid points
    REAL(dp), INTENT(IN)     ::  tslclim(klon,12)     !! climatological surface temperature for each kalendar month
    REAL(dp), INTENT(OUT)    ::  tsoil(klon,ntsoil)   !! temperature of the soil
    REAL(dp), INTENT(OUT)    ::  tsl(klon)            !! surface temperature at time step t

    !  Local
    REAL(dp) :: zdmax(klon), zkap, zsqrt, zday, yearl
    REAL(dp) :: znmea(klon), zrange(klon)
    INTEGER  ::  jg, jl, im
    INTEGER  ::  iyday
    INTEGER  ::  jmax(klon), nmomid(12), nmonthl(12)
    INTEGER  ::  yr, mo, dy, hr, mn, se

    REAL(dp), PARAMETER ::  cmid_offset = 0.07_dp

    ! Computational constants

    zkap = 7.5E-7_dp

    !
    !    Date handling
    !

    ! get year, month, day of start date

    CALL get_date_components(start_date, yr, mo, dy, hr, mn, se)

    ! get length of the year

    yearl = get_year_len(yr)

    ! get length of month

    DO im = 1,12
       nmonthl(im) = REAL(get_month_len(yr,im),dp)
       nmomid(im)  = nmonthl(im) * 0.5_dp
    END DO

    ! calendar day of the start date

    IF (nmonth == 0) THEN
       ! annual cycle run
       iyday = dy
       DO im = 1, mo-1
          iyday = iyday + nmonthl(im)
       END DO
    ELSE
       ! perpetual months
       iyday = nmomid(nmonth)
       DO im = 1, nmonth-1
          iyday = iyday + nmonthl(im)
       END DO
    END IF
    zday = REAL(iyday,dp)
    zsqrt = SQRT(zkap * yearl * REAL(ndaylen,dp) / pi)

    ! calendar month of maximum surface temperature

    jmax(:) = MAXLOC(tslclim(:,:),DIM=2)

    ! difference between months of maximum and minimum surface temperature

    zrange(:) = MAXVAL(tslclim(:,:),DIM=2) - MINVAL(tslclim(:,:),DIM=2)

    ! calendar day of maximum surface temperature 

    zdmax(:) = 0._dp
    DO jl= 1,klon
       DO im = 1,jmax(jl)
          zdmax(jl) = zdmax(jl) + REAL(nmonthl(im),dp)
       END DO
       zdmax(jl) = zdmax(jl) - REAL(nmomid(jmax(jl)),dp)
    END DO

    !
    !    Initialize soil temperatures

    ! calculate annual mean surface temperature

    znmea(:) = SUM(tslclim(:,:),DIM=2) / 12._dp

    ! set soil temperature of the five layers

    DO jg = 1,ntsoil
       tsoil(:,jg) = znmea(:) + 0.5_dp * zrange(:) * EXP(-(cmid(jg)-cmid_offset) / zsqrt)    &
            * COS(2._dp * pi * (zday - zdmax(:)) / yearl - (cmid(jg)-cmid_offset) / zsqrt)
    END DO

    !  set all time levels of surface temperature to uppermost soil temp.

    tsl(:)   = tsoil(:,1)

  END SUBROUTINE ini_soil_temp

  SUBROUTINE update_cdrag( &                         
       kland, &                                           
       HeightWind, HeightHumidity, &
       coef_ril_tm1, coef_ril_t, coef_ril_tp1, &
       temp_air, qair, windspeed, &
       surface_temp, surface_temp_upd, &
       pressure, &
       csat, cair, &
       z0m, z0h, &
       zril_old, cdrag, zchl)
    
    USE mo_time_control,     ONLY: time_step_len
    USE mo_physc2,           ONLY: cb, cc, cvdifts, ckap    
    USE mo_atmosphere,       ONLY: sat_specific_humidity
    USE mo_jsbach_constants, ONLY: SpecificHeatDryAirConstPressure, GasConstantDryAir, Gravity, vtmpc1


    INTEGER,               INTENT(in)    :: kland                     !! Number of land points in vectors
    REAL(dp),              INTENT(in)    :: HeightWind                !! Defines lowest layer height-> where measurements are taken
    REAL(dp),              INTENT(in)    :: HeightHumidity            !! Defines lowest layer height-> where measurements are taken
    REAL(dp),              INTENT(in)    :: coef_ril_tm1              !! Weighting factor for richardson numbers at different steps,
    REAL(dp),              INTENT(in)    :: coef_ril_t                !!   that are used to calculate a drag coef. that approximates
    REAL(dp),              INTENT(in)    :: coef_ril_tp1              !!   the drag coef. at time t, but helps to maintain stability.
    REAL(dp),              INTENT(in)    :: temp_air(kland)           !! Lowest level air temperature [Kelvin]
    REAL(dp),              INTENT(in)    :: qair(kland)               !! Lowest level specific humidity
    REAL(dp),              INTENT(in)    :: windspeed(kland)          !! Lowest level windspeed
    REAL(dp),              INTENT(in)    :: surface_temp(kland)       !! Surface temperature
    REAL(dp),              INTENT(in)    :: surface_temp_upd(kland)   !! Surface temperature
    REAL(dp),              INTENT(in)    :: pressure(kland)           !! Surface pressure
    REAL(dp),              INTENT(in)    :: csat(kland)
    REAL(dp),              INTENT(in)    :: cair(kland)
    REAL(dp),              INTENT(in)    :: z0h(kland)
    REAL(dp),              INTENT(in)    :: z0m(kland)
    REAL(dp),              INTENT(inout) :: zril_old(kland)           !! Richardson Number at previous timestep
    REAL(dp),              INTENT(out)   :: cdrag(kland)              !! drag coefficient updated
    REAL(dp),              INTENT(out)   :: zchl(kland)               !! drag coefficient updated

    REAL(dp)                   :: zkappa, zsigma, zcons9, zcons11, zcons12
    REAL(dp), DIMENSION(kland) :: ztvir,  ztvd, ztvs, ztvs_upd, air_pressure,      &
                                  geopot_surf, geopot_surf_hum, windshear_squared, &
                                  qsat_surf, qsat_surf_upd, zril, zril_initial, zril_upd, &
                                  zchnl, zcfnchl, zcfhl, zucfhl, zscfl, zcons

    !-------------------------------------------------------------------------------------------------
    ! Some constants from echam/mo_surface_land.f90
    !-------------------------------------------------------------------------------------------------
    zkappa        = GasConstantDryAir / SpecificHeatDryAirConstPressure

    !-------------------------------------------------------------------------------------------------
    ! Some constants from echam/mo_surface_land.f90
    !-------------------------------------------------------------------------------------------------
    zcons9        = 3._dp * cb
    zcons11       = 3._dp * cb * cc
    zcons12       = cvdifts * time_step_len * Gravity / GasConstantDryAir
   
    ! virtual potential air temperature (see mo_surface_boundary.f90)
    ! according to Saucier, WJ Principles of Meteoroligical Analyses
    ! tv = t * (1 + 0.61 * q * t)    ! virtual temperature
    ! td = t * ( 1000 / p_mb ) ^ R/cdp  ! potential temperature
    ! tvd = tair * (100000/p_pa)^zkappa * 1 + vtmpc1 * qair) ! virtual potential temperature
    zsigma    = 0.99615_dp

!   *** Calculate squared wind shear as minimum wind speed square from echam/mo_surface_land.f90:precalc_land
    windshear_squared(:) = MAX(windspeed(:)**2, 1._dp)
    
    air_pressure(:) = zsigma * pressure(:)

    ! virtual potential air temperature (see mo_surface_boundary.f90)
    ! according to Saucier, WJ Principles of Meteoroligical Analyses
    ! tv = t * (1 + 0.61 * q * t)    ! virtual temperature
    ztvir(:)  = temp_air(:) * ( 1._dp + vtmpc1 * qair(:) )

    ! geopotential of the surface layer (see echam's auxhybc.f90 & geopot.f90)
    ! If HeightWind is set, then the measurement height + an offset from the average vegetation height
    ! is used. Otherwise, the code defaults to the half-level of ECHAM's lowest atmospheric layer
    IF(HeightWind > 0._dp) THEN
       geopot_surf(:)  = HeightWind * Gravity
       geopot_surf_hum(:) = HeightHumidity * Gravity ! HeightHumidity and HeightTemperature have to be equal
    ELSE
      geopot_surf(:)  = ztvir(:) *  GasConstantDryAir * LOG(1._dp / zsigma)
      geopot_surf_hum(:) = geopot_surf(:)
    ENDIF

    
    ztvd(:)   = temp_air(:) * ( 100000._dp/air_pressure(:))**zkappa * ( 1._dp + vtmpc1 * qair(:) )
    
    ! virtual potential surface temperature for timestep t and t + 1
    qsat_surf(:) = sat_specific_humidity(surface_temp(:),pressure(:))
    ztvs = surface_temp(:) * ( 100000._dp/pressure(:))**zkappa * &
           ( 1._dp + vtmpc1 * (csat(:) * qsat_surf(:) + ( 1._dp - cair(:) ) * qair(:)))

    qsat_surf_upd(:) = sat_specific_humidity(surface_temp_upd(:),pressure(:))
    ztvs_upd = surface_temp_upd(:) * ( 100000._dp/pressure(:))**zkappa * &
           ( 1._dp + vtmpc1 * (csat(:) * qsat_surf_upd(:) + ( 1._dp - cair(:) ) * qair(:)))

    ! Richardson number (dry, Brinkop & Roeckner 1995, Tellus)
    ! ztvd, ztvs are now virtual potential temperatures, changed by Thomas Raddatz 07.2014
    ! for timestep t and t + 1 
    zril_initial(:) = geopot_surf(:) * ( ztvd(:) - ztvs(:) )  / ( windshear_squared(:) * (ztvd(:)+ztvs(:))/2 )
    zril_upd(:) = geopot_surf(:) * ( ztvd(:) - ztvs_upd(:) )  &
               / ( windshear_squared(:) * (ztvd(:)+ztvs_upd(:))/2 )

    ! Richrdson number needs to be filtered for offline simulations in order to maintain stability 
    zril(:) = coef_ril_tm1 * zril_old(:)     &
            + coef_ril_t   * zril_initial(:) &
            + coef_ril_tp1 * zril_upd(:)
    zril_old(:) = zril(:)    

    !------------------------------------------------------------------------------------
    ! Approximation of cdrag
    !------------------------------------------------------------------------------------

    ! Neutral drag coefficient for momentum and heat
    zchnl(:) = ckap**2 / (LOG(1._dp + geopot_surf(:) / (Gravity * z0m(:) )) &
    			* LOG( ( Gravity * z0m(:) + geopot_surf_hum(:) ) / (Gravity * z0h(:) )) )

    ! account for stable/unstable case: helper variables
    zscfl(:) = SQRT (  1._dp + 5._dp * ABS(zril(:)))
    zucfhl(:) = 1._dp / (1._dp + zcons11 * zchnl(:) * SQRT(ABS(zril(:)) * (1._dp  &
            + geopot_surf_hum(:) / (Gravity * z0h(:)))))

    ! ignoring cloud water correction (see mo_surface_land.f90)
    zcons(:) = zcons12 * pressure(:) / ( temp_air(:) * (1._dp + vtmpc1 * qair(:)))
    zcfnchl(:)  = zcons(:) * SQRT(windshear_squared(:)) * zchnl(:)

    ! Stable / Unstable case
    WHERE ( zril(:) .GT. 0._dp )
       zcfhl(:) = zcfnchl(:) / (1._dp + zcons9 * zril(:) * zscfl(:))
    ELSEWHERE
       zcfhl(:) = zcfnchl(:) * (1._dp - zcons9 * zril(:) * zucfhl(:))
    ENDWHERE

    cdrag(1:kland) = zcfhl(:)
    zchl(1:kland)  = zcfhl(:) / zcfnchl(:) * zchnl(:)

  END SUBROUTINE update_cdrag


END MODULE mo_soil
