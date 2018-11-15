!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_soil_hydrology

  USE mo_time_control,     ONLY: delta_time, lstart
  USE mo_kind,             ONLY: dp
  USE mo_exception,        ONLY: finish, message, message_text, em_warn

  IMPLICIT NONE

  PRIVATE

  PUBLIC update_surface_hydrology, update_soil_hydrology

  REAL(dp), PARAMETER ::  zwdmin = 0.05_dp
  REAL(dp), PARAMETER ::  zwdcri = 0.9_dp
  REAL(dp), PARAMETER ::  zdrmin = 0.001_dp / (3600._dp*1000._dp)  ! factor for minimum drainage (~ 10 hours)
  REAL(dp), PARAMETER ::  zdrmax = 0.1_dp   / (3600._dp*1000._dp)  ! factor for maximum drainage (~ 4 days)
  REAL(dp), PARAMETER ::  zdrexp = 1.5_dp
  REAL(dp), PARAMETER ::  zfcmin = 1.e-10_dp              ! minimum field capacity (or rather epsilon?)
  INTEGER :: max2d(2)                                     ! temporarily needed for for maximum locations

CONTAINS

  SUBROUTINE update_surface_hydrology ( nidx,                    &
       is_land, is_glacier,  q_air, wind_10m, t_air,             &
       skinres_canopy_max, skinres_max, sfract_srf, hcap_grnd,   &
       evapotrans, evapopot, rain, snow,                         &
       t_srf_unfilt, skinres, snow_soil, canopy_snow,            &
       glacier_depth, tte_corr,                                  &
       snow_accum, snow_melt_soil, glacier_melt, glacier_runoff, &
       pme_glacier, evapotrans_no_snow, melt_to_soil, evap_deficit)

    !----------------------------------------------------------------------------------------------
    ! Routine based on the former JSBACH3 routine update_surf_down
    !                                             ----------------
    !     rewritten in Mai 2016 by V. Gayler
    !
    !----------------------------------------------------------------------------------------------

    USE mo_param_switches,   ONLY: lsurf
    USE mo_jsbach_constants, ONLY: cvinter, tmelt, cpd, gravity, rhoh2o, vtmpc2, alf


    INTEGER,  INTENT(in) :: nidx                       ! vector length / chunk dimension
    LOGICAL,  INTENT(in) :: is_land(nidx)              ! land mask
    LOGICAL,  INTENT(in) :: is_glacier(nidx)           ! glacier mask
    REAL(dp), INTENT(in) :: q_air(nidx)                ! specific humidity at lowest atmospheric level
    REAL(dp), INTENT(in) :: wind_10m(nidx)             ! wind speed at 10m height [m/s]
    REAL(dp), INTENT(in) :: t_air(nidx)                ! lowest layer atmosphere temperature [K]
    REAL(dp), INTENT(in) :: skinres_canopy_max(nidx)   ! capacity of the canopy skin reservoir [m]
    REAL(dp), INTENT(in) :: skinres_max(nidx)
    REAL(dp), INTENT(in) :: sfract_srf(nidx)           ! surface snow fraction []
    REAL(dp), INTENT(in) :: hcap_grnd(nidx)            ! heat capacity of the uppermost soil layer [J m-2 K-1]
    REAL(dp), INTENT(in) :: evapotrans(nidx)           ! evapotranspiration (including sublimation) [kg m-2 s-1]
    REAL(dp), INTENT(in) :: evapopot(nidx)             ! potential evaporation [kg m-2 s-1]
    REAL(dp), INTENT(in) :: rain(nidx)                 ! liquid precipitation [kg m-2 s-1]
    REAL(dp), INTENT(in) :: snow(nidx)                 ! solid precipitation [kg m-2 s-1]

    REAL(dp), INTENT(inout) :: t_srf_unfilt(nidx)      ! surface tempemperature [K] (unfiltered)
    REAL(dp), INTENT(inout) :: skinres(nidx)           ! water content of the skin reservoir (vegetation and bare soil) [m]
    REAL(dp), INTENT(inout) :: snow_soil(nidx)         ! snow depth at the ground [m water equivalent]
    REAL(dp), INTENT(inout) :: canopy_snow(nidx)       ! snow depth on the canopy [m water equivalent]
    REAL(dp), INTENT(inout) :: glacier_depth(nidx)     ! glacier depth (snow and ice) [m water equivalent]
    REAL(dp), INTENT(inout) :: tte_corr(nidx)          ! temperature tendency correction due to snow melt

    REAL(dp), INTENT(out)   :: snow_accum(nidx)        ! snow budget at non-glacier points [m water equivalent]
    REAL(dp), INTENT(out)   :: snow_melt_soil(nidx)    ! snow/ice melt at land points (excluding canopy) [m water equivalent]
    REAL(dp), INTENT(out)   :: glacier_melt(nidx)      ! Snow/ice melt at glacier points [m water equivalent]
    REAL(dp), INTENT(out)   :: glacier_runoff(nidx)    ! glacier runoff (rain+snow/ice melt, no calving) [m water equivalent]
    REAL(dp), INTENT(out)   :: pme_glacier(nidx)       ! precipitation minus sublimation on glacier [m water equivalent]
    REAL(dp), INTENT(out)   :: evapotrans_no_snow(nidx)  ! evapotranspiration without snow evaporation [m]
    REAL(dp), INTENT(out)   :: melt_to_soil(nidx)      ! melt water available for infiltration into the soil [m water equivalent]
    REAL(dp), INTENT(out)   :: evap_deficit(nidx)      ! Evaporation deficit flux due to inconsistent treatment of snow evap. [m]

    !
    !  local variables
    !
    REAL(dp) :: rain_in_m (nidx)         ! amount of rainfall within time step [m] 
    REAL(dp) :: new_snow(nidx)           ! amount of snowfall within time step [m]
    REAL(dp) :: evapotrans_in_m(nidx)    ! amount of evapotranspiration within time step [m]
    REAL(dp) :: evap_snow_in_m(nidx)     ! amount of sublimation from snow within time step [m]
    REAL(dp) :: snow_melt_can(nidx)      ! snow melt from canopy [m]
    REAL(dp) :: snow_melt_pot(nidx)      ! potential amount of snow melt [m]
    REAL(dp) :: new_snow_canopy(nidx)    ! amount of snowfall on canopy within time step [m]
    REAL(dp) :: new_snow_soil(nidx)      ! amount of snowfall to soil within time step [m]
    REAL(dp) :: exp_t(nidx), exp_w(nidx) ! exponents needed for unloading of snow due to temperature / wind
    REAL(dp) :: snow_blown(nidx)         ! amount of snow blown from canopy to the ground
    REAL(dp) :: factor(nidx)             ! factor to calculate cooling due to snow/ice melt
    REAL(dp) :: canopy_pre_snow(nidx)    ! snow depth on the canopy prior to its update
    REAL(dp) :: evap_snow_pot_in_m(nidx) ! potential snow evaporation

    !
    !  Parameters  - compare chapter 2 "Land Physics" of the JSBACH3 documentation
    !
    REAL(dp), PARAMETER :: zc1 = tmelt - 3._dp ! [K]
    REAL(dp), PARAMETER :: zc2=1.87E5_dp       ! [Ks]
    REAL(dp), PARAMETER :: zc3=1.56E5_dp       ! [m]

    !----------------------------------------------------------------------------------------------

    ! Initializations
    !-----------------

    pme_glacier(:)      = 0._dp
    glacier_runoff(:)   = 0._dp
    melt_to_soil(:)     = 0._dp
    snow_accum(:)       = 0._dp       ! new_snow
    snow_melt_soil(:)   = 0._dp
    glacier_melt(:)     = 0._dp
    snow_melt_can(:)    = 0._dp


    ! Convert water fluxes to m water equivalent within the time step
    !-----------------------------------------------------------------

    rain_in_m(:)          = rain(:)       * delta_time/rhoh2o
    new_snow(:)           = snow(:)       * delta_time/rhoh2o   !in_m       ! snow_fall
    evapotrans_in_m(:)    = evapotrans(:) * delta_time/rhoh2o
    evap_snow_pot_in_m(:) = sfract_srf(:) * evapopot(:) * delta_time/rhoh2o
    evapotrans_no_snow(:) = evapotrans_in_m(:) - evap_snow_pot_in_m(:)


    IF (lsurf) THEN
     
      !  Budgets of snow (canopy, ground, glaciers)
      !---------------------------------------------

      WHERE (is_land(:) .AND. .NOT. is_glacier(:))

        ! amount of snow remaining on the canopy
        new_snow_canopy(:) = MIN(new_snow(:) * cvinter, skinres_canopy_max(:) - canopy_snow(:))

        ! remaining snow falls on the ground 
        new_snow_soil(:) = new_snow(:) - new_snow_canopy(:)

        ! update snow on the canopy
        ! Note: evaporation happens from canopy, as long as there is enough snow
        canopy_pre_snow(:) = canopy_snow(:)
        canopy_snow(:) = MIN(MAX(0._dp, canopy_pre_snow(:) + new_snow_canopy(:) + evap_snow_pot_in_m(:)), skinres_canopy_max(:))
        evap_snow_in_m(:) = evap_snow_pot_in_m(:) + canopy_pre_snow(:) + new_snow_canopy(:) - canopy_snow(:)

        ! unloading of snow from the canopy due to melting
        exp_t(:) = MAX(0._dp, t_air(:)-zc1)/zc2 * delta_time
        snow_melt_can(:) = canopy_snow(:) * (1._dp-EXP(-exp_t(:)))
        canopy_snow(:) = canopy_snow(:) - snow_melt_can(:)

        ! unloading of snow from the canopy due to wind
        exp_w(:) = wind_10m(:)/zc3 * delta_time
        snow_blown(:) = canopy_snow(:) * (1._dp-EXP(-exp_w(:)))
        canopy_snow(:)   = canopy_snow(:)   - snow_blown(:)
        new_snow_soil(:) = new_snow_soil(:) + snow_blown(:)

        ! temperature correction term due to snow melt
        factor(:) = alf / (cpd * (1._dp + vtmpc2*MAX(0.0_dp,q_air(:))))
        tte_corr(:) = snow_melt_can(:) * rhoh2o * gravity/delta_time * factor(:)

      ELSEWHERE
        canopy_snow(:) = 0._dp
      END WHERE


      !  Snowfall and sublimation on non-glacier land
      !-----------------------------------------------

      WHERE (is_land(:) .AND. .NOT. is_glacier(:))
        snow_soil(:) = snow_soil(:) + new_snow_soil(:) + evap_snow_in_m(:)
      ELSEWHERE
        snow_soil(:) = 0._dp
      ENDWHERE

      ! Correction if there was too much snow evaporation
      WHERE (snow_soil(:) < 0._dp)
        evapotrans_no_snow(:) = evapotrans_no_snow(:) + snow_soil(:)
        evap_deficit(:) = snow_soil(:)
        snow_soil(:) = 0._dp
      ELSEWHERE
        evap_deficit(:) = 0._dp
      END WHERE


      !  Snowfall and sublimation on glaciers
      !---------------------------------------

      WHERE (is_glacier(:))
        glacier_depth(:)  = glacier_depth(:) + new_snow(:) + evap_snow_pot_in_m(:)   ! glacier depth [m]
        pme_glacier(:)    = rain_in_m(:)     + new_snow(:) + evapotrans_in_m(:)      ! P-E on glaciers [m]
        glacier_runoff(:) = rain_in_m(:)                                             ! there is no infiltration on glacies
      END WHERE


      !  Snow and glacier melt
      !------------------------

      IF (.NOT. lstart) THEN      ! hcap_ground = 0. at model start
        WHERE (is_land(:))
          WHERE (t_srf_unfilt(:) > tmelt)
            snow_melt_pot(:) = hcap_grnd(:) * (t_srf_unfilt(:)-tmelt)/(alf*rhoh2o) !  potential snow and ice melt [m]
            WHERE (is_glacier(:))
              glacier_melt(:)   = snow_melt_pot(:)                          ! there is an unlimited amount of snow
              glacier_depth(:)  = glacier_depth(:) - glacier_melt(:)        ! reduce glacier depth according to melting
              glacier_runoff(:) = glacier_runoff(:) + glacier_melt(:)       ! add melt water to the runoff
              t_srf_unfilt(:) = tmelt                                       ! surface temp. on glaciers is limited to 0 degC
            ELSEWHERE (snow_soil(:) > 0._dp)
              snow_melt_soil(:) = MIN(snow_melt_pot(:), snow_soil(:))            ! snow melt limited by actual snow depth [m]
              snow_soil(:) = snow_soil(:) - snow_melt_soil(:)                    ! reduce snow depth according to melting
              t_srf_unfilt(:) = t_srf_unfilt(:) - snow_melt_soil(:)*alf*rhoh2o/hcap_grnd(:)  ! re-calculation of surface temp.
            END WHERE
          END WHERE
        END WHERE
      END IF


      !  Snow budget and meltwater (glacier-free land only)
      !-----------------------------------------------------

      WHERE (is_land(:) .AND. .NOT. is_glacier(:))
        skinres(:) = skinres(:) + snow_melt_can(:)                                    ! Add melt water from canopy to skin reservoir
        melt_to_soil(:) = snow_melt_soil(:) + MAX(0._dp, skinres(:)-skinres_max(:))   ! Excess water enters the soil, and
        skinres(:) = MIN(skinres_max(:), skinres(:))                                  ! skin reservoir is limited to the maximum
        snow_accum(:) = new_snow(:) + evap_snow_pot_in_m(:) - snow_melt_soil(:) - snow_melt_can(:)
      END WHERE

    END IF

    RETURN

  END SUBROUTINE update_surface_hydrology
  !------------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------------
  SUBROUTINE update_soil_hydrology (nidx, nsoil,                   &
       is_land, is_glacier, sfract_srf, wfract_srf,                &
       skinres_max, ws_max_mineral, oro_stddev, evapopot, t_soil,  &
       rain, melt_to_soil, evapotrans_no_snow,                     &
       skinres, ws, evap_deficit,                                  &
       evapotrans_soil, runoff, drainage, evap_skin,               &
       ! multi-layer hydrology fields (optional)
       dz, root_depth, soil_depth, hyd_cond_sat, matrix_pot,       &
       bclapp, vol_porosity, trans, vol_field_cap, vol_p_wilt,     &
       pore_size_index, org_layer, ws_l, field_cap_l, jllog,       &
       reduce_evap,                                                &
       ! permafrost fields (optional)
       soil_ice)

    !------------------------------------------------------------------------------------------------

    USE mo_param_switches,            ONLY: lsurf
    USE mo_jsbach_constants,          ONLY: cvinter, rhoh2o, tmelt
#ifdef STANDALONE
    USE mo_jsbach_comm_to_echam5mods, ONLY: nlat
#else
    USE mo_control,                   ONLY: ngl
#endif
 
    INTEGER,  INTENT(in)    :: nidx                  ! vector length
    INTEGER,  INTENT(in)    :: nsoil                 ! number of soil layers
    LOGICAL,  INTENT(in)    :: is_land(nidx)         ! land mask
    LOGICAL,  INTENT(in)    :: is_glacier(nidx)      ! glacier mask
    REAL(dp), INTENT(in)    :: sfract_srf(nidx)      ! snow cover fraction []
    REAL(dp), INTENT(in)    :: wfract_srf(nidx)      ! wet surface fraction []
    REAL(dp), INTENT(in)    :: skinres_max(nidx)     ! maximum size of the skin reservoir [m]
    REAL(dp), INTENT(in)    :: ws_max_mineral(nidx)  ! water holding capacity of the non organic soil [m]
    REAL(dp), INTENT(in)    :: oro_stddev(nidx)      ! subgrid standard deviation of the orographie [m]
    REAL(dp), INTENT(in)    :: evapopot(nidx)        ! potential evaporation/sublimation [kg m-2 s-1]
    REAL(dp), INTENT(in)    :: t_soil(nidx)          ! temperature of the uppermost soil layer [K]
    REAL(dp), INTENT(in)    :: rain(nidx)            ! rainfall [kg m-2 s-1]
    REAL(dp), INTENT(in)    :: melt_to_soil(nidx)    ! meltwater for infiltration in the soil [m]
    REAL(dp), INTENT(in)    :: evapotrans_no_snow(nidx) ! evapotranspiration without snow evaporation [m]

    REAL(dp), INTENT(inout) :: skinres(nidx)         ! water in the skin reservoir [m]
    REAL(dp), INTENT(inout) :: ws(nidx)              ! soil water content (root zone) [m] 
    REAL(dp), INTENT(inout) :: evap_deficit(nidx)    ! evaporation deficit due to inconsistent treatment of snow & skin evap. [m]
    REAL(dp), INTENT(out)   :: evapotrans_soil(nidx) ! evapotranspiration from the soil [m]
    REAL(dp), INTENT(out)   :: runoff(nidx)          ! runoff (does NOT include drainage) [m]
    REAL(dp), INTENT(out)   :: drainage(nidx)        ! drainage [m]

    REAL(dp), INTENT(out)   :: evap_skin(nidx)       ! evaporation from the skin reservoir [m]

    ! variables needed with multi-layer soil hydrology
    REAL(dp), INTENT(in),    OPTIONAL :: dz(nsoil)               ! thickness of the soil layers [m]
    REAL(dp), INTENT(in),    OPTIONAL :: root_depth(nidx)        ! rooting depth [m]
    REAL(dp), INTENT(in),    OPTIONAL :: soil_depth(nidx)        ! soil depth until bedrock [m]
    REAL(dp), INTENT(in),    OPTIONAL :: hyd_cond_sat(nidx)      ! hydraulic conductivity of saturated soil [m/s]
    REAL(dp), INTENT(in),    OPTIONAL :: matrix_pot(nidx)        ! soil matrix potential [m]
    REAL(dp), INTENT(in),    OPTIONAL :: bclapp(nidx)            ! exponent B in Clapp and Hornberge
    REAL(dp), INTENT(in),    OPTIONAL :: vol_porosity(nidx)      ! volumetric soil porosity [m/m]   following Cosby et al.
    REAL(dp), INTENT(in),    OPTIONAL :: trans(nidx)
    REAL(dp), INTENT(in),    OPTIONAL :: vol_field_cap(nidx)     ! volumetric soil field capacity [m/m]
    REAL(dp), INTENT(in),    OPTIONAL :: vol_p_wilt(nidx)        ! volumetric soil wilting point [m/m]
    REAL(dp), INTENT(in),    OPTIONAL :: pore_size_index(nidx)   ! soil pore size distribution index
    LOGICAL,  INTENT(in),    OPTIONAL :: org_layer(nidx)         ! mask for organic layer presence []
    REAL(dp), INTENT(inout), OPTIONAL :: ws_l(nidx,nsoil)        ! soil moisture content of each layer [m]
    REAL(dp), INTENT(inout), OPTIONAL :: field_cap_l(nidx,nsoil) ! field capacity of each soil layer [m]
    INTEGER,  INTENT(inout), OPTIONAL :: jllog                   ! gridbox number for logoutput in routinr soilhyd

    REAL(dp), INTENT(out),   OPTIONAL :: reduce_evap(nidx)       ! diagnostic desired reduction of evapotranspiration
                                                                 !   due to limited soil water content
    ! variable needed with permafrost sheme
    REAL(dp), INTENT(inout), OPTIONAL :: soil_ice(nidx,nsoil)    ! soil ice content in each soil layer [m]
 
    !
    !  local variables
    !
    INTEGER  :: jl, jk                 ! looping indces
    REAL(dp) :: rain_in_m(nidx)        ! amount of rain falling within time step [m]
    REAL(dp) :: rain_to_soil(nidx)     ! amount of rain going to the soil [m]
    REAL(dp) :: infilt(nidx)           ! infiltration [m]
    REAL(dp) :: rain_to_skinres(nidx)  ! amount of rain going into the skin reservoir
    REAL(dp) :: ws_crit(nidx)          ! soil moisture above which drainage strongly increases
    REAL(dp) :: ws_min_drain(nidx)     ! minimum soil moisture for drainage
    REAL(dp) :: ws_incl_ice(nidx)      ! soil moisture including ice [m]
    REAL(dp) :: ws_rel                 ! relative soil moisture (0: wilting point, 1: field capacity)
    REAL(dp) :: ice_rootzone(nidx)     ! amount of ice in the rootzone [m]
    REAL(dp) :: steepness              ! sub-grid scale steepness of the orographie
    REAL(dp) :: sigma_0, sigma_max     ! parameters for infiltration calculation
    REAL(dp) :: zbws, zb1, zbm, zconw1 !   ''
    REAL(dp) :: zvol                   !
    REAL(dp) :: water_to_soil          ! amount of water going into the soil [m]
    REAL(dp) :: droot(nidx,nsoil)      ! rooted depth within the soil layers (until root_depth) [m]
    REAL(dp) :: dsoil(nidx,nsoil)      ! soil depth until bedrock within the soil layers [m]
    REAL(dp) :: wsat_l(nidx,nsoil)     ! maximum water storage in each soil layer (porosity) [m]
    REAL(dp) :: p_wilt_l(nidx,nsoil)   ! wilting point of the soil layer
    REAL(dp) :: evap_soil(nidx)        ! evaporation from the soil (without transpiration) [m]
    REAL(dp) :: trans_in_m(nidx)       ! transpiration within time step [m]
    REAL(dp) :: evap_skin_pot(nidx)    ! potential transpiration from the skin reservoir [m] 
    REAL(dp) :: skinres_pre_rain(nidx) ! amount of water in the skin reservoir prior to the update [m]
    REAL(dp) :: ws_max(nidx)           ! water holding capacity of the soil [m]

    INTEGER,  PARAMETER :: ilog  = 0               ! Switch for debugging output

    !----------------------------------------------------------------------------------------------

    ! Initializations
    !-----------------

    runoff(:)          = 0._dp
    drainage(:)        = 0._dp
    reduce_evap(:)     = 0._dp
    evap_skin(:)       = 0._dp
    evapotrans_soil(:) = 0._dp

    !
    !      Parameters
    !
    sigma_0 = 100._dp
#ifndef STANDALONE
    sigma_max = 1000._dp*64._dp/ngl
#else
    sigma_max = 1000._dp*64._dp/nlat
    IF (nlat == 2) sigma_max = 1000._dp*64._dp/96   ! for site level jsbach: T63 value
#endif

    ! Convert water fluxes to m water equivalent within the time step
    !-----------------------------------------------------------------

    rain_in_m(:) = rain(:) * delta_time / rhoh2o
    trans_in_m(1:nidx) = trans(1:nidx)  * delta_time / rhoh2o

    WHERE (is_land(:))
      evap_skin_pot(:) = (1._dp-sfract_srf(:))*wfract_srf(:)*evapopot(:) * delta_time/rhoh2o
    END WHERE


    ! if organic layer is used, ws_max will be modified below
    ws_max(:) = ws_max_mineral(:)

    ! Initialization of working arrays for multi-layer soil scheme
    !--------------------------------------------------------------
    IF (nsoil > 1) THEN
      CALL soildef(nidx, nsoil, is_land, dz, root_depth, soil_depth,            &
           vol_porosity, vol_field_cap, vol_p_wilt, org_layer, soil_ice, ws_l,  &
           ws_max, droot, dsoil, wsat_l, field_cap_l, p_wilt_l, ice_rootzone, ws)
    ENDIF

    IF (lsurf) THEN

      !--------------------------------------------------------------------------------------------
      !  water budgets
      !--------------------------------------------------------------------------------------------

      !  skin reservoir (vegetation and bare soil)
      ! ----------------

      WHERE (is_land(:) .AND. .NOT. is_glacier(:))

        ! interception of rain
        rain_to_skinres(:) = MIN(rain_in_m(:) * cvinter, skinres_max(:)-skinres(:))
        rain_to_soil(:) = rain_in_m(:) - rain_to_skinres(:)
        skinres_pre_rain(:) = skinres(:)
        skinres(:) = MIN(MAX(0._dp, skinres_pre_rain(:) + rain_to_skinres + evap_skin_pot(:)), skinres_max(:))

        ! note: at this point ignoring the amount of dew that exceeds the skinreservoir
        !       as dew is calculated later based on evap_soil (f(evapotrans_soil))
        evap_skin(:) = skinres(:) - (skinres_pre_rain(:) + rain_to_skinres(:))
        evapotrans_soil(:) = evapotrans_no_snow(:) - evap_skin(:)
        WHERE(evap_skin_pot(:) < 0._dp)
          evap_deficit(:) = evap_deficit(:) + evap_skin_pot(:) - evap_skin(:)
        END WHERE
      ELSEWHERE
        skinres(:) = 0._dp     ! no water in skin reservoir on glaciers
      END WHERE

      ! calculation of bare soil evaporation
      IF (nsoil > 1) evap_soil(:) = evapotrans_soil(:) - trans_in_m(:)

      !  soil water reservoir
      ! ----------------------

      ! infiltration
      !
      ! f(w) = 1 - (1 - w/w_max)**b
      !
      ! b: shape parameter, defining the sub-gid scale steepness of the orographie:
      !
      ! b = (sigma_z - sigma_0) / (sigma_z + sigma_max)
      !
      ! sigma_z: standard deviation of topographic height
      ! sigma_0: minimum value (100 m); below b = 0.01
      ! sigma_max: resolution dependent maximum

      DO jl = 1, nidx

        infilt(jl) = 0._dp

        IF (is_land(jl) .AND. .NOT. is_glacier(jl)) THEN
          ws_min_drain(jl) = zwdmin * ws_max(jl)            ! minimum amount of soil water for drainage

          steepness = MAX(0._dp, oro_stddev(jl)-sigma_0) / (oro_stddev(jl)+sigma_max)

          zbws   = MAX(MIN(steepness,0.5_dp),0.01_dp)   ! wie passt das zur Doku???
          zb1    = 1._dp+zbws
          zbm    = 1._dp/zb1
          zconw1 = ws_max(jl)*zb1

          !  surface runoff, infiltration and evaporation from soil
          ! --------------------------------------------------------

          IF (nsoil == 1) THEN          ! bucket scheme

            ! positive evapotranspiration (dew) is added to water_to_soil,
            !   negative evapotranspiration (real evapotranspiration) is taken from the soil water
            IF (evapotrans_soil(jl) >= 0.0_dp) THEN
              water_to_soil = melt_to_soil(jl) + rain_to_soil(jl) + evapotrans_soil(jl)
            ELSE
              ws(jl)        = ws(jl)    + evapotrans_soil(jl)
              water_to_soil = melt_to_soil(jl) + rain_to_soil(jl)
            END IF

          ELSE                          ! multi-layer soil hydrology

            ! positive values of evaporation and transpiration (dew) are added to water_to_soil
            ! Note: negative fluxes change soil moisture in routine soilchange
            water_to_soil = melt_to_soil(jl) + rain_to_soil(jl)
            IF (evap_soil(jl) >= 0.0_dp)  water_to_soil = water_to_soil + evap_soil(jl)
            IF (trans_in_m(jl) >= 0.0_dp) water_to_soil = water_to_soil + trans_in_m(jl)

          ENDIF

          IF (t_soil(jl) < tmelt) THEN
            runoff(jl) = water_to_soil       !  no infiltration as uppermost soil layer is frozen -> runoff

          ELSE

            ws_incl_ice(jl) = ws(jl) + ice_rootzone(jl)                    ! soil water including ice
            IF (water_to_soil > 0._dp .AND. ws_incl_ice(jl) > ws_min_drain(jl)) THEN  ! soil water content above limit for drainage

              ws_rel = MIN(1._dp, ws_incl_ice(jl) / ws_max(jl))            !  relative soil moisture
              zvol   = (1._dp-ws_rel)**zbm - water_to_soil/zconw1          ! ???

              runoff(jl) = water_to_soil - (ws_max(jl)-ws_incl_ice(jl))    ! if >0: amount of water exceeding water holding capacity
              IF (zvol > 0._dp) runoff(jl) = runoff(jl) + ws_max(jl)*zvol**zb1
              runoff(jl) = MAX(MIN(runoff(jl), water_to_soil), 0._dp)      ! no negative runoff, runoff limited by available water
              infilt(jl) = water_to_soil - runoff(jl)                      ! water not going into runoff is infiltrated
            ELSE
              runoff(jl) = 0._dp                                           ! no runoff as soil is too dry
              infilt(jl) = water_to_soil                                   ! all water infiltrated
            END IF

            ! In the bucket scheme infiltration is added to soil moisture.
            ! In multi-layer scheme, soil moisture change due to infiltration is considered in routine soilchange.
            IF (nsoil == 1) ws(jl) = ws(jl) + infilt(jl)

          END IF
        ENDIF
      ENDDO


      ! drainage and runoff
      !---------------------

      IF (nsoil == 1) THEN                 ! bucket scheme

        WHERE (is_land(:) .AND. .NOT. is_glacier(:))

          ws_crit(:) = zwdcri * ws_max(:)          ! soil moisture above which drainage strongly increases
          ws_min_drain(:) = zwdmin * ws_max(:)     ! minimum soil moisture for drainage 

          WHERE (ws(:) <= ws_min_drain(:))
            drainage(:) = 0._dp
          ELSEWHERE
            WHERE (t_soil(:) > tmelt)
              drainage(:) = zdrmin * MAX(ws(:)/ws_max(:), 0.7_dp)    ! until 70% drainage stays the same
              WHERE (ws(:) > ws_crit(:))
                drainage(:) = drainage(:) + (zdrmax-zdrmin)*((ws(:)-ws_crit(:))/(ws_max(:)-ws_crit(:)))**zdrexp
              END WHERE
              drainage(:) = drainage(:) * delta_time                 ! drainage within time step [m]
              drainage(:) = MIN(drainage(:), ws(:)-ws_min_drain(:))  ! drainage limited by minimum soil moisture
              ws(:) = ws(:) - drainage(:)
            ELSEWHERE
              drainage(:) = 0._dp                  ! no drainage in frozen soil
            END WHERE
          END WHERE
          !!!vg: drainage is not added to runoff in multi-layer scheme !!!! 
          runoff(:) = runoff(:) + drainage(:) + MAX(ws(:)-ws_max(:), 0._dp)
          ws(:)  = MIN(ws_max(:), ws(:))
        ELSEWHERE
          ws(:) = 0._dp
        END WHERE

      ELSE         ! multi-layer scheme

        ! initilalization of drainage 
        drainage(:) = 0._dp

        ! Changes by evaporation fluxes
        CALL digest_evapotrans(nidx, nsoil, is_land, root_depth, droot, dsoil, &
             field_cap_l, trans_in_m, evap_soil, ws_l, soil_ice, reduce_evap)

        ! For numerical reasons, half of the water is infiltrated 
        ! before soilhyd and the other half of it afterwards
        infilt(:) = infilt(:) / 2._dp

        ! update soil moisture according to infiltration
        !    - first half of the water
        CALL digest_infilt(nidx, nsoil, is_land, field_cap_l, infilt, &
             ws_l, drainage)

        ! Drainage calculation with multi-layer soil scheme
        CALL soilhyd(nidx, nsoil, delta_time, ilog, jllog, ws_l, dsoil, &
             wsat_l, field_cap_l, hyd_cond_sat, vol_porosity, bclapp,   &
             matrix_pot, drainage, pore_size_index, p_wilt_l, dz,       &
             soil_ice, org_layer)

        ! update soil moisture according to infiltration
        !    - second half of the water
        CALL digest_infilt(nidx, nsoil, is_land, field_cap_l, infilt, &
             ws_l, drainage)


        ! recalculate soil moisture from the soil moisture of the individual levels
        ! Note: ws only includes the soil moisture within the root zone
        ws(:) = 0.0_dp
        DO jk=1, nsoil
          WHERE (is_land(:))
            WHERE (droot(:,jk) >= dz(jk))
              ws(:) = ws(:) + ws_l(:,jk)
            ELSEWHERE (droot(:,jk) > 0._dp)
              ws(:) = ws(:) + ws_l(:,jk) * droot(:,jk)/dsoil(:,jk)
            END WHERE
          END WHERE
        END DO

      END IF ! multi-layer scheme

    END IF ! lsurf

  END SUBROUTINE update_soil_hydrology
  !------------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------------
  SUBROUTINE soilhyd(nidx, nsoil, delta_time, ilog, jllog, ws_l, dsoil,      &
       wsat_l, field_cap_l, hyd_cond_sat, vol_porosity, bclapp, matrix_pot,  &
       drainage, pore_size_index, p_wilt_l, dz, soil_ice, org_layer)



    !----------------------------------------------------------------------------------------------
    !
    !  SOIL HYDROLOGY - Routine to be build into REMO to calculate soil hydrology for 5 soil layers
    !
    !  Calculation of percolation and vertical diffusion of water
    !
    !  Programming and development: Stefan Hagemann (MPI-M)
    !----------------------------------------------------------------------------------------------
    USE mo_jsbach_constants,          ONLY: vol_vpor_org, hyd_cond_sat_org, bclapp_org, &
                                            matrix_pot_org, pore_size_index_org 

    INTEGER,  INTENT(IN)            :: nidx                    ! vector length
    INTEGER,  INTENT(IN)            :: nsoil                   ! number of soil layers
    REAL(dp), INTENT(IN)            :: delta_time              ! model time step legth [s]
    INTEGER,  INTENT(IN)            :: ilog                    !
    INTEGER,  INTENT(IN)            :: jllog                   !
    REAL(dp), INTENT(IN)            :: dsoil(nidx, nsoil)      ! soil depth until bedrock within each layer [m]
    REAL(dp), INTENT(IN)            :: wsat_l(nidx,nsoil)      ! Maximum storage in each soil layer (porosity) [m]
    REAL(dp), INTENT(IN)            :: field_cap_l(nidx,nsoil) ! field capacity of the soil layer [m]
    REAL(dp), INTENT(IN)            :: hyd_cond_sat(nidx)      ! hydraulic conductivity of saturated soil [m/s]
    REAL(dp), INTENT(IN)            :: vol_porosity(nidx)      ! volumetric soil porosity [m/m]
    REAL(dp), INTENT(IN)            :: bclapp(nidx)            ! exponent B in Clapp and Hornberger
    REAL(dp), INTENT(IN)            :: matrix_pot(nidx)        ! matrix potential [m]
    REAL(dp), INTENT(IN)            :: pore_size_index(nidx)   ! soil pore size distribution index used in Van Genuchten
    REAL(dp), INTENT(IN)            :: p_wilt_l(nidx,nsoil )   ! wilting point of the soil layer
    REAL(dp), INTENT(IN)            :: dz(nsoil)               ! vertical depth of each soil layer
    REAL(dp), INTENT(INOUT)         :: ws_l(nidx, nsoil)       ! soil moisture of each layer [m]
    REAL(dp), INTENT(INOUT)         :: drainage(nidx)          ! amount of drainage within timestep [m]
    REAL(dp), INTENT(IN), OPTIONAL  :: soil_ice(nidx, nsoil)   ! amount of ice in the soil [m]
    LOGICAL,  INTENT(IN), OPTIONAL  :: org_layer(nidx)         ! mask for organic layer []

    INTEGER,  PARAMETER :: iredu  = 0         ! reduction of conductivity with depth yes/no (1/0)
    REAL(dp), PARAMETER :: credu  = 1.0_dp

    INTEGER  :: i                             ! looping index
    REAL(dp) :: percolation(nidx,0:nsoil)     ! downward flux from one soil layer to the layer below [m/s]
                                              !     and later amount of perculation per time step [m]
    REAL(dp) :: diffus_l(nidx,nsoil)          ! diffusivity at mid-depth of the layer
    REAL(dp) :: ws_min_drain(nidx)            ! minimum soil moisture for drainage
    REAL(dp) :: ws_rel(nidx)                  ! normalized (~relative) soil mositure
    REAL(dp) :: ws_crit(nidx)                 ! soil moisture above which drainage strongly increases [m]
    REAL(dp) :: drain_crit(nidx)              ! drainage above a soil moisture of ws_crit [m/s]
    REAL(dp) :: ws_tmp(nidx)                  ! temporary array for soil moisture
    REAL(dp) :: drain_l(nidx)                 ! drainage from the current layer
    LOGICAL  :: soil_above(nidx)              ! is there a soil layer above? 
    LOGICAL  :: soil_below(nidx)              ! does the layer below have a soil fraction?
    REAL(dp) :: zda(nidx,nsoil)               ! diffusion coefficients
    REAL(dp) :: zdb(nidx,nsoil)               !     ''
    REAL(dp) :: zdc(nidx,nsoil)               !     ''
    REAL(dp) :: ztri(nidx,nsoil)              !     ''
    REAL(dp) :: diffusivity(nidx,nsoil-1)     ! diffusivity between two layers, i.e. at the bottom of each layer
    REAL(dp) :: zkredu(nsoil)                 ! reduction factor of hydraulic conductivity with depth for percolation
    REAL(dp) :: zdiflog(10), wslog(nsoil)     ! needed for debugging output

    REAL(dp) :: hyd_cond_sat_l(nidx,nsoil)    ! hydraulic conductivity of layer 1 [m/s]
    REAL(dp) :: matrix_pot_l(nidx,nsoil)      ! matrix potential of layer 1 [[m]
    REAL(dp) :: porosity_l(nidx,nsoil)        ! soil porosity of layer 1 [m]
    REAL(dp) :: bclapp_l(nidx,nsoil)          ! exponent B in Clapp and Hornberger of layer 1
    REAL(dp) :: pore_size_index_l(nidx,nsoil) ! soil pore size distribution index used in Van Genuchten of layer 1
    REAL(dp) :: zmvg_l(nidx,nsoil)            ! exponent M in Van Genuchten scheme = pore_size_index / (pore_size_index + 1)

    ! define reduction factor of hydraulic conductivity with depth (for percolation)
    IF (iredu == 1) THEN
      zkredu(1) = (1._dp + EXP(-credu*dz(1)) )/2._dp
      DO i = 2, nsoil
        zkredu(i) = (EXP(-credu*SUM(dz(1:i-1))) + EXP(-credu*SUM(dz(1:i)))) / 2._dp
      END DO
    ELSE
      ! conductivity is constant with depth
      zkredu(1:nsoil) = 1._dp
    END IF

    hyd_cond_sat_l(:,:)   = SPREAD(hyd_cond_sat(:),   DIM=2,ncopies=nsoil)
    porosity_l(:,:)       = SPREAD(vol_porosity(:),   DIM=2,ncopies=nsoil)
    bclapp_l(:,:)         = SPREAD(bclapp(:),         DIM=2,ncopies=nsoil)
    matrix_pot_l(:,:)     = SPREAD(matrix_pot(:),     DIM=2,ncopies=nsoil)
    pore_size_index_l(:,:)= SPREAD(pore_size_index(:),DIM=2,ncopies=nsoil)
    zmvg_l(:,:)           = pore_size_index_l(:,:) / (pore_size_index_l(:,:)+1._dp)
    WHERE(org_layer(:))
        hyd_cond_sat_l(:,1)   = hyd_cond_sat_org
        porosity_l(:,1)       = vol_vpor_org
        bclapp_l(:,1)         = bclapp_org
        matrix_pot_l(:,1)     = matrix_pot_org
        pore_size_index_l(:,1)= pore_size_index_org
        zmvg_l(:,1)           = pore_size_index_l(:,1) / (pore_size_index_l(:,1)+1._dp) 
    END WHERE

    !----------------------------------------------------------------------------------------------
    ! Percolation = Gravitational drainage
    !
    ! de Vrese: percolation scheme using Van Genuchten-formulation and Euler upstream
    ! ----------------------------------------------------------------------------------------------
    ! Van Genuchten-formulation
    !
    ! percolation
    !     = relative hydraulic conductivity (K_r) * hydraulic conductivity at staturation (hydr_cond_sat)
    !
    !    K_r = SQRT(THETA) * (1-(1-THETA**1/m)**m)**2
    !
    ! with
    !  THETA = (theta - theta_r)/(theta_s - theta_r)
    !
    !  THETA:   normalised soil moisture     -> ws_rel
    !  theta:   current soil moisture        -> ws_l
    !  theta_s: soil moisture at stauration  -> field_cap
    !  theta_r: residual soil moisture       -> p_wilt

    percolation(:,nsoil) = 0._dp
    DO i = 1, nsoil-1
      WHERE (field_cap_l(:,i) > zfcmin .AND. field_cap_l(:,i) > p_wilt_l(:,i))
        ws_rel(:) = MAX(ws_l(:,i) - p_wilt_l(:,i), 0._dp) / (field_cap_l(:,i) - p_wilt_l(:,i))
        ws_rel(:) = MIN(ws_rel(:), 1._dp)
        percolation(:,i) =  hyd_cond_sat_l(:,i)*zkredu(i) * SQRT(ws_rel(:)) &
                           * (1._dp - (1._dp - ws_rel(:)**(1._dp/zmvg_l(:,i)))**zmvg_l(:,i))**2._dp &
                           * delta_time     ! [m/s] -> [m]
      ELSEWHERE
        percolation(:,i) = 0._dp
      END WHERE
    END DO


    ! drainage of the lowest layer above bedrock replaces the percolation calculated above

    DO i = 1, nsoil

      ! is there a layer above or a layer with soil below?
      soil_above(:) = (i > 1)
      IF (i == nsoil) THEN
        soil_below(:) = .false.
      ELSE
        soil_below(:) = dsoil(:,i+1) > 0._dp
      END IF

      ws_crit(:) = field_cap_l(:,i) * zwdcri     ! soil moisture above which drainage strongly increases

      ! calculate drainage of the lowest layer above bedrock
      WHERE (dsoil(:,i) > 0._dp .AND. .NOT. soil_below(:))          ! lowest layer above bedrock

        WHERE (ws_l(:,i) <= p_wilt_l(:,i))
          drain_l(:) = 0._dp              ! no drainage below the wilting point
        ELSEWHERE
          WHERE (field_cap_l(:,i) > p_wilt_l(:,i))
            ! drain_min = minimum factor * normalized soil moisture
            drain_l(:) = zdrmin * (ws_l(:,i) - p_wilt_l(:,i)) / (field_cap_l(:,i) - p_wilt_l(:,i))
          ELSEWHERE
            drain_l(:) = 0._dp
          END WHERE
          WHERE (ws_l(:,i) > ws_crit(:))
            drain_crit(:) = (zdrmax-zdrmin) * ((ws_l(:,i)-ws_crit(:)) / (wsat_l(:,i) - ws_crit(:)))**zdrexp
            drain_l(:) = drain_l(:) + MIN(drain_crit(:), (ws_l(:,i)-ws_crit(:))/delta_time)
          END WHERE
          ! drainage is limited by the available water of the layer
          drain_l(:) = MIN(drain_l(:), (ws_l(:,i)-p_wilt_l(:,i))/delta_time)
        END WHERE
        percolation(:,i) = drain_l(:) * delta_time  ! [m/s] -> [m]

      ELSEWHERE (dsoil(:,i) > 0._dp)

        ! Drainage is calculated with percolation removed from WS_L
        !   de Vrese --- WS_L corrected for incomming percolation

        ! Remove percolating water from the corresponding soil layer and add percolation from the above layer
        WHERE (.NOT. soil_above(:))       ! surface layer
          ws_tmp(:) = ws_l(:,i) - percolation(:,i)
        ELSEWHERE
          ws_tmp(:) = MIN(field_cap_l(:,i), ws_l(:,i) - percolation(:,i) + percolation(:,i-1))
        END WHERE

        WHERE (ws_tmp(:) > p_wilt_l(:,i))          ! drainage occurs only above wilting point
          WHERE (field_cap_l(:,i) > p_wilt_l(:,i))
          ! drain_min = minimum factor * normalized soil moisture
            drain_l(:) = zdrmin * (ws_tmp(:) - p_wilt_l(:,i)) / (field_cap_l(:,i)-p_wilt_l(:,i))
          ELSEWHERE
            drain_l(:) = 0._dp
          END WHERE
          WHERE (ws_tmp(:) > ws_crit(:))
            drain_crit(:) = (zdrmax-zdrmin) * ((ws_tmp(:) - ws_crit(:)) / (wsat_l(:,i) - ws_crit(:)))**zdrexp
            drain_l(:) = drain_l(:) + MIN(drain_crit(:), (ws_tmp(:) - ws_crit(:))/delta_time)
          END WHERE
          ! drainage is limited by the soil moisture within the layer
          drain_l(:) = MIN(drain_l(:) * delta_time, (ws_tmp(:) - p_wilt_l(:,i)), (field_cap_l(:,i) * (1._dp - zwdmin)))
          ws_l(:,i) = ws_l(:,i) - drain_l(:)       ! remove drainage of the layer from soil moisture
          drainage(:) = drainage(:) + drain_l(:)   ! and add it to the total drainage
        END WHERE
      END WHERE

    END DO  ! soil level

    !---------------------------------------------------------------------------------------------- 
    !  water diffusion
    ! -----------------

    ! numerical recipes/Richtmyer & -Morton diffusion

    ! calculating the diffusivity at mid depth of layer i
    DO i = 1, nsoil

      ws_rel(:) = 0._dp
        WHERE (ws_l(:,i) > 0._dp)
          ws_rel(:) =  ws_l(:,i) / dsoil(:,i)  &
                      / (porosity_l(:,i) - (soil_ice(:,i) / dsoil(:,i)))
        END WHERE

      diffus_l(:,i) = 0._dp
      WHERE (ws_rel(:) > zfcmin)
        diffus_l(:,i) =  bclapp_l(:,i) * hyd_cond_sat_l(:,i) * matrix_pot_l(:,i)  &
                        / (ws_l(:,i) / dsoil(:,i)) *  ws_rel(:) ** (bclapp_l(:,i)+3._dp)
      END WHERE
    END DO

    ! calculating the diffusivity at the bottom of layer i
    DO i = 1, nsoil-1
      WHERE (diffus_l(:,i+1) == 0._dp .OR. diffus_l(:,i) == 0._dp)
        diffusivity(:,i) = 0._dp
      ELSEWHERE
        diffusivity(:,i) = diffus_l(:,i+1) * diffus_l(:,i) / &
                        ((diffus_l(:,i+1)*dsoil(:,i) + diffus_l(:,i)*dsoil(:,i+1))/(dsoil(:,i+1)+dsoil(:,i)))
      END WHERE
    END DO

    ! calculation of diffusion coefficients

    ! deepest layer
    zda(:,nsoil) = 0._dp
    WHERE (dsoil(:,nsoil) > 0._dp)
      zdb(:,nsoil) = diffusivity(:,nsoil-1) * delta_time / dsoil(:,nsoil) / (dsoil(:,nsoil)+dsoil(:,nsoil-1)) * 2._dp
    ELSEWHERE
      zdb(:,nsoil) = 0._dp
    END WHERE

    ! middle layers
    DO i = nsoil-1, 2, -1
      WHERE (dsoil(:,i) > 0._dp)
        WHERE (dsoil(:,i+1) > 0._dp)
          zda(:,i) = diffusivity(:,i) * delta_time / dsoil(:,i) / (dsoil(:,i)+dsoil(:,i+1)) * 2._dp
        ELSEWHERE
          zda(:,i) = 0._dp
        END WHERE
        zdb(:,i) = diffusivity(:,i-1) * delta_time / dsoil(:,i) / (dsoil(:,i)+dsoil(:,i-1)) * 2._dp
      ELSEWHERE
        zda(:,i) = 0._dp
        zdb(:,i) = 0._dp
      END WHERE
    END DO

    ! uppermost layer
    zdb(:, 1) = 0._dp
    WHERE (dsoil(:,1) > 0)
      WHERE (dsoil(:,2) > 0)
        zda(:,1) = diffusivity(:,1) * delta_time / dsoil(:,1) / (dsoil(:,1)+dsoil(:,2)) * 2._dp
      ELSEWHERE
        zda(:,1) = 0._dp
      END WHERE
    ELSEWHERE
      zda(:,1) = 0._dp
    END WHERE

    ! temporary storage for diagnostic printing
    IF (ilog == 1 .OR. ilog == 3) wslog(1:nsoil) = ws_l(jllog, 1:nsoil)


    ! ---------------------
    !  update soil wetness  - with respect to diffusion
    ! ---------------------
    ! routine TRIDIAG from Numerical Recipes, p. 43
    !   -ZDA = CI, -ZDB=AI, ZDA+ZDB+1=BI, WS_L(T)=RI, WS_L(T+1)=UI

    WHERE (dsoil(:,1) > 0._dp)
      zdc(:,1) = zda(:,1) + zdb(:,1) + 1._dp
      ws_l(:,1) = ws_l(:,1) / zdc(:,1)
    END WHERE
 
    ! decomposition and forward substitution
    DO i = 2, nsoil
      WHERE (dsoil(:,i) > 0._dp)
        ztri(:,i) = -zda(:,i-1) / zdc(:,i-1)
        zdc(:,i)  = zda(:,i) + zdb(:,i) + 1._dp + zdb(:,i)*ztri(:,i)
        ws_l(:,i) = (ws_l(:,i) / dsoil(:,i) + zdb(:,i) * ws_l(:,i-1)/dsoil(:,i-1)) / zdc(:,i) * dsoil(:,i)
      END WHERE
    END DO

    ! backsubstitution
    DO i = nsoil-1, 1, -1
      WHERE (dsoil(:,i+1) > 0._dp)
        ws_l(:,i) = ws_l(:,i) - ztri(:,i+1) * ws_l(:,i+1) * (dsoil(:,i)/dsoil(:,i+1))
      END WHERE
    END DO

    ! temporary storage for diagnostic printing
    IF (ilog == 1 .OR. ilog == 3) THEN
      zdiflog(1:nsoil) = wslog(1:nsoil) - ws_l(jllog,1:nsoil)
    END IF


    ! treatment of layer overflow due to diffusion
    !    limit soil moisture to field capacity and add overflow to drainage
    DO i = 1, nsoil
      WHERE (dsoil(:,i) > 0._dp .AND. ws_l(:,i) > field_cap_l(:,i))
        drainage(:) = drainage(:) + ws_l(:,i)-field_cap_l(:,i)
        ws_l(:,i) = field_cap_l(:,i)
      END WHERE
    END DO

    !----------------------------------------------------------------------------------------------
    ! debugging output for grid cell jllog
    IF (ilog == 2 .OR. ilog == 3) THEN
      write(*,'(A10,1X, 11(E16.9E2, 1X) )')  &
        "SOILHYD1: ", hyd_cond_sat(jllog)*1000._dp, ws_l(jllog,1:nsoil)*1000._dp,  &
        field_cap_l(jllog,1:nsoil)*1000._dp
    END IF
    !----------------------------------------------------------------------------------------------


    ! ---------------------
    !  update soil wetness  - with respect to percolation
    ! ---------------------

    ! percolation per time step must be smaller than the usable soil layer water content
    DO i = 1, nsoil
      WHERE (dsoil(:,i) > 0._dp .AND. percolation(:,i) > 0._dp)
        ws_min_drain(:) = zwdmin * field_cap_l(:,i)              ! minimum soil moisture for drainage
        WHERE (ws_l(:,i)-percolation(:,i) < ws_min_drain(:))     ! more percolation than water available 
          WHERE (ws_l(:,i) > ws_min_drain(:))                    !   there is water for perculation
            percolation(:,i) = ws_l(:,i) - ws_min_drain(:)       !   all available water perculates
            ws_l(:,i) = ws_min_drain(:)                          !   soil moisture is limited to minimum
          ELSEWHERE
            percolation(:,i) = 0._dp
          END WHERE
        ELSEWHERE                                                ! there is enough soil water
          ws_l(:,i) = ws_l(:,i) - percolation(:,i)
        END WHERE
      END WHERE
    END DO

    ! Percolation is added to downward layer or to drainage, if lowest layer or bedrock is reached.
    DO i = 1, nsoil-1
      WHERE (dsoil(:,i+1) > 0._dp)
        ws_l(:,i+1) = ws_l(:,i+1) + percolation(:,i)
      ELSEWHERE             ! bedrock below
        drainage(:) = drainage(:) + percolation(:,i)
        percolation(:,i) = 0._dp
      END WHERE
    END DO

    ! Percolation overflow is not allowed from layers below the surface layer.
    ! Instead it is assumed that water piles upwards from saturated layers.
    DO i = nsoil, 2, -1
      WHERE (dsoil(:,i) > 0._dp .AND. ws_l(:,i) > field_cap_l(:,i))
        ! soil moisture above field capacity:
        !   reduce percolation from the above layer by the amount that does not fit into the below layer
        percolation(:,i-1) = MAX(percolation(:,i-1) - ws_l(:,i) + field_cap_l(:,i), 0._dp)
        ws_l(:,i-1) = ws_l(:,i-1) + ws_l(:,i) - field_cap_l(:,i)
        ws_l(:,i) = field_cap_l(:,i)
      END WHERE
    END DO

    drainage(:) = drainage(:) + percolation(:,nsoil)

    ! Overflow from the surface layer enters drainage. This should not happen,
    !   this is just for security.
    WHERE (ws_l(:,1) > field_cap_l(:,1))
      drainage(:) = drainage(:) + ws_l(:,1) - field_cap_l(:,1)
      ws_l(:,1) = field_cap_l(:,1)
    END WHERE


    !----------------------------------------------------------------------------------------------
    ! debugging output for grid cell jllog
    IF (ilog == 1 .OR. ilog == 3) THEN
      DO i=1, nsoil
        zdiflog(nsoil+i) = percolation(jllog,i)
      END DO
      write(*,'(A10,1X,11(E15.6E3,1X))') "SOILFLUX",zdiflog(1:10)*1000._dp, drainage(jllog)*1000._dp
      write(*,'(A10,1X,10(E16.9E2,1X))') "SOILHYD2: ", ws_l(jllog,1:nsoil)*1000._dp, dsoil(jllog,1:nsoil)*1000._dp
    END IF

  END SUBROUTINE soilhyd
  !------------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------------
  SUBROUTINE soildef (nidx, nsoil, is_land, dz, root_depth, soil_depth,                &
       vol_porosity, vol_field_cap, vol_p_wilt, org_layer, soil_ice, soil_water,       &
       ws_max, droot, dsoil, wsat_l, field_cap_l, p_wilt_l, ice_rootzone, water_rootzone)
    !----------------------------------------------------------------------------------------------
    !
    !  Routine to calculate soil layer arrays from arrays comprising the entire soil, and the other
    !  way round.
    !
    !     ******** Version 1.0 - Juli 2006
    !              Programming and development: Stefan Hagemann
    !     ******** Version 1.1 - Juni 2007
    !              Optimization : Ralf Podzun
    !     ******** Version 1.2 - September 2007
    !              A loop optimized : Stefan Hagemann
    !     ******** Version 1.3 - Oktober 2007
    !              Minor technical correction without effect: Stefan Hagemann
    !     ******** Version 1.4 - January 2008
    !              Minor technical correction without effect: Stefan Hagemann
    !     ******** Version 1.5 - February 2008 - Stefan Hagemann
    !              Minor technical correction that may effect only 1-2 gridboxes
    !              where ROOT_DEPTH=0. (near lakes)
    !     ******** Version 1.6 - March 2009 - Stefan Hagemann
    !              Volumetric soil field capacity vol_field_cap and wilting point vol_p_wilt are
    !              given via input list and do not need to be computed in soildef.
    !              new output: P_WILT_L, WSMX eliminated.
    !
    !----------------------------------------------------------------------------------------------

    USE mo_jsbach_constants,  ONLY: vol_vpor_org, vol_field_cap_org, vol_p_wilt_org

    INTEGER,  INTENT(IN)   :: nidx                    ! vector length
    INTEGER,  INTENT(IN)   :: nsoil                   ! number of soil layers
    LOGICAL,  INTENT(IN)   :: is_land(nidx)           ! land mask
    REAL(dp), INTENT(IN)   :: dz(nsoil)               ! soil layers thickness [m]
    REAL(dp), INTENT(IN)   :: root_depth(nidx)        ! rooting depth in [m]
    REAL(dp), INTENT(IN)   :: soil_depth(nidx)        ! soil depth until bedrock [m]
    REAL(dp), INTENT(IN)   :: vol_porosity(nidx)      ! volumetric porosity [m water / m soil]
    REAL(dp), INTENT(IN)   :: vol_field_cap(nidx)     ! volumetric field capacity [m water / m soil]
    REAL(dp), INTENT(IN)   :: vol_p_wilt(nidx)        ! volumetric permanent wilting point [m water / m soil]
    LOGICAL,  INTENT(IN)   :: org_layer(nidx)         ! mask for organic layer [/]
    REAL(dp), INTENT(IN)   :: soil_ice(nidx,nsoil)    ! with permafrost: amount of soil ice [m]
    REAL(dp), INTENT(IN)   :: soil_water(nidx,nsoil)  ! with permafrost: amount of soil water [m]    
    REAL(dp), INTENT(INOUT):: ws_max(nidx)            ! root zone water holding capacity [m]
    REAL(dp), INTENT(OUT)  :: droot(nidx,nsoil)       ! rooted depth per soil layer (until root_depth) [m]
    REAL(dp), INTENT(OUT)  :: dsoil(nidx,nsoil)       ! soil depth per soil layer (until bedrock) [m]
    REAL(dp), INTENT(OUT)  :: wsat_l(nidx,nsoil)      ! maximum storage capacity per soil layer (porosity) [m]
    REAL(dp), INTENT(OUT)  :: field_cap_l(nidx,nsoil) ! field capacity per soil layer [m]
    REAL(dp), INTENT(OUT)  :: p_wilt_l(nidx,nsoil)    ! texture depending wilting point per soil layer [m]
    REAL(dp), INTENT(OUT)  :: ice_rootzone(nidx)      ! with permafrost: amount of soil ice in the root zone [m]
    REAL(dp), INTENT(OUT)  :: water_rootzone(nidx)    ! with permafrost: amount of soil water in the root zone [m]

    INTEGER  :: i                          ! looping index 
    REAL(dp) :: ldepth(0:nsoil)            ! depth of the lower boundary of each soil layer (from soil surface)
    REAL(dp) :: vol_porosity_top(nidx)     ! top layer volumetric porosity [m water / m soil]
    REAL(dp) :: vol_field_cap_top(nidx)    ! top layer volumetric field capacity [m water / m soil]
    REAL(dp) :: vol_p_wilt_top(nidx)       ! top layer volumetric permanent wilting point [m water / m soil]

    !----------------------------------------------------------------------------------------------
    !  initializations

    dsoil(:,:) = 0._dp
    droot(:,:) = 0._dp

    !----------------------------------------------------------------------------------------------

    ! layer depth (from the soil surface to the lower boundary of each layer)
    !-------------
    ldepth(0) = 0._dp
    DO i = 1, nsoil
      ldepth(i) = ldepth(i-1) + dz(i)
    END DO

    ! rooting depth within each layer
    !---------------
    DO i = 1, nsoil
      WHERE (is_land(:))
        WHERE (root_depth(:) >= ldepth(i))
          droot(:,i) = dz(i)
        ELSEWHERE (root_depth(:) > ldepth(i-1))
          droot(:,i) = root_depth(:) - ldepth(i-1)
        ELSEWHERE
          droot(:,i) = 0._dp
        END WHERE
      END WHERE
    END DO

    ! soil depth until bedrock within each layer
    !------------
    DO i = 1, nsoil
      WHERE (is_land(:))
        WHERE (soil_depth(:) >= ldepth(i))
          dsoil(:,i) = dz(i)
        ELSEWHERE (soil_depth(:) > ldepth(i-1))
          dsoil(:,i) = soil_depth(:) - ldepth(i-1)
        ELSEWHERE
          dsoil(:,i) = 0._dp
        END WHERE
      END WHERE
    END DO

    ! water storage capacity of each layer
    !------------------------
      ! consideration of soil ice

      WHERE (org_layer)
        vol_porosity_top  (:)= vol_vpor_org
        vol_field_cap_top (:)= vol_field_cap_org
        vol_p_wilt_top    (:)= vol_p_wilt_org
      ELSEWHERE
        vol_porosity_top  (:)= vol_porosity (:)
        vol_field_cap_top (:)= vol_field_cap(:)
        vol_p_wilt_top    (:)= vol_p_wilt   (:)
      END WHERE

      wsat_l(:,1) = MAX(vol_porosity_top(:) * dsoil(:,1) - soil_ice(:,1), 0._dp)
      field_cap_l(:,1) = MAX(vol_field_cap_top(:) * dsoil(:,1) - soil_ice(:,1), 0._dp)
      WHERE (vol_porosity_top(:) > 0._dp)
        p_wilt_l(:,1) = wsat_l(:,1) * vol_p_wilt_top(:) / vol_porosity_top(:)
      ELSEWHERE
        p_wilt_l(:,1) = 0._dp
      END WHERE

      DO i = 2, nsoil
        wsat_l(:,i) = MAX(vol_porosity(:) * dsoil(:,i) - soil_ice(:,i), 0._dp)
        field_cap_l(:,i) = MAX(vol_field_cap(:) * dsoil(:,i) - soil_ice(:,i), 0._dp)
        WHERE (vol_porosity(:) > 0._dp)
          p_wilt_l(:,i) = wsat_l(:,i) * vol_p_wilt(:) / vol_porosity(:)
        ELSEWHERE
          p_wilt_l(:,i) = 0._dp
        END WHERE
      END DO

      ! root zone ice content
      ice_rootzone(:)   = 0._dp
      water_rootzone(:) = 0._dp
      DO i = 1, nsoil
        WHERE (is_land(:))
          WHERE (droot(:,i) >= dz(i))
            ice_rootzone(:)   = ice_rootzone(:)   + soil_ice(:,i)
            water_rootzone(:) = water_rootzone(:) + soil_water(:,i)
          ELSEWHERE (droot(:,i) > 0._dp)
            ice_rootzone(:)   = ice_rootzone(:)   + soil_ice(:,i)   * droot(:,i)/dsoil(:,i)
            water_rootzone(:) = water_rootzone(:) + soil_water(:,i) * droot(:,i)/dsoil(:,i)
          END WHERE
        END WHERE
      END DO

     ! root zone water holding capacity where organic layer
      WHERE (is_land(:).AND. org_layer)
        ws_max(:) = 0._dp
      END WHERE
      DO i = 1, nsoil
        WHERE (is_land(:) .AND. org_layer)
          WHERE (droot(:,i) >= dz(i))
            ws_max(:) = ws_max(:) + field_cap_l(:,i)
          ELSEWHERE (droot(:,i) > 0._dp)
            ws_max(:) = ws_max(:) + field_cap_l(:,i) * droot(:,i)/dsoil(:,i)
          END WHERE
        END WHERE
      END DO

  END SUBROUTINE soildef
  !------------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------------
  SUBROUTINE digest_evapotrans(nidx, nsoil, is_land, root_depth, droot, dsoil, &
       field_cap_l, trans_in_m, evap_soil, ws_l, soil_ice, reduce_evap)
    !----------------------------------------------------------------------------------------------
    !  The routine is based on routine soilchange by Stefan Hagemann. It replaces the part for
    !  evapotranspiration changes (isch=1)
    !
    !  Update soil moisture due to evaporation and transpiration fluxes. The fluxes had been
    !  calculated with the old bucket scheme, and are not necessarily consistent to soil moisture
    !  of the multi-layer scheme.
    !  Therefore evaporation might need to be reduced (reduce_evap). However this term is only
    !  exported for diagnostics and the water deficit is distributed between the soil layers.
    !----------------------------------------------------------------------------------------------

    ! arguments

    INTEGER,  INTENT(in)     :: nidx                    ! vector length
    INTEGER,  INTENT(in)     :: nsoil                   ! number of soil layers
    LOGICAL,  INTENT(in)     :: is_land(nidx)           ! land mask
    REAL(dp), INTENT(in)     :: root_depth(nidx)        ! rooting depth [m]
    REAL(dp), INTENT(in)     :: droot(nidx, nsoil)      ! root depth within layer [m]
    REAL(dp), INTENT(in)     :: dsoil(nidx, nsoil)      ! soil depth until bedrock within layer [m]
    REAL(dp), INTENT(in)     :: field_cap_l(nidx,nsoil) ! field capacity of the layer
    REAL(dp), INTENT(in)     :: trans_in_m(nidx)        ! transpiration [m]
    REAL(dp), INTENT(in)     :: evap_soil(nidx)         ! bare soil evaporation [m]

    REAL(dp), INTENT(inout)  :: ws_l(nidx,nsoil)        ! water content of the soil layer [m]
    REAL(dp), INTENT(inout)  :: soil_ice(nidx,nsoil)    ! ice content of the soil layer [m]
    REAL(dp), INTENT(inout)  :: reduce_evap(nidx)

    ! local variables

    INTEGER  :: jk
    REAL(dp) :: deficit_evap(nidx)   ! water deficit due to evaporation
    REAL(dp) :: deficit_trans(nidx)  ! water deficit due to tranpiration
    REAL(dp) :: deficit_l(nidx)      ! water deficit within a soil layer
    REAL(dp) :: remaining(nidx)      ! water remaining in the soil layer
    REAL(dp) :: rootfract(nidx)      ! fraction of the soil layer within the root zone
    REAL(dp) :: fixed(nidx)          ! amount of water below the wilting point, not available for plant
    REAL(dp) :: trans_layer(nidx)    ! transpiration from the layer
    REAL(dp) :: ws(nidx)             ! water within all soil layers

    ! parameters

    REAL(dp), PARAMETER :: zwilt  = 0.35_dp 

    !----------------------------------------------------------------------------------------------
    !  update soil moisture due to evaporation and transpiration fluxes
    !----------------------------------------------------------------------------------------------

    !  bare soil evaporation
    ! -----------------------

    deficit_evap(:) = 0._dp   ! water missing for evaporation due to limited soil water
    WHERE (is_land(:))
      WHERE (evap_soil(:) < 0._dp)            ! pos. evaporation flux (not dew)
        ws_l(:,1) = ws_l(:,1) + evap_soil(:)  ! take water from uppermost soil layer 
        WHERE (ws_l(:,1) < 0._dp)             ! if there is not enough soil water
          deficit_evap(:) = ws_l(:,1)         !    this is the amount of water missing
          ws_l(:,1) = 0._dp                   !    soil moisture cannot be negative
        END WHERE
      END WHERE
    END WHERE

    !  transpiration
    ! ---------------

    deficit_trans(:) = 0._dp
    deficit_l(:) = 0._dp

    ! sealed grid cells (rooting_depth zero) and positive transpiration
    !    no water reachable for transpiration, add flux to deficit term
    WHERE (root_depth(:) <= 0._dp .AND. trans_in_m(:) < 0._dp)
      deficit_trans(:) = deficit_trans(:) + trans_in_m(:)
    END WHERE

    DO jk =1, nsoil
      WHERE (is_land(:))

        ! grid cell not sealed and positive transpiration (not dew)
        !    take water equally from all over the root zone
        WHERE (root_depth(:) > 0._dp .AND. trans_in_m(:) < 0._dp)
            
          ! water needed for transpiration from the current layer 
          trans_layer(:) = trans_in_m(:) * droot(:,jk)/root_depth(:)

          WHERE (droot(:,jk) > 0._dp)

            rootfract(:) = droot(:,jk) / dsoil(:,jk)             ! fraction of the soil layer that is within the root zone
            fixed(:) = field_cap_l(:,jk) * zwilt * rootfract(:)  ! amount of water in the root zone, not available for plants

            WHERE (ws_l(:,jk) * rootfract(:) >= fixed(:))        ! there is water in the root zone available for plants
              remaining(:) = ws_l(:,jk) * rootfract(:) + trans_layer(:)  ! water that would remain after transpiration
              WHERE (remaining(:) < fixed(:))                                ! transpiration exceeds available water
                deficit_l(:) = remaining(:) - fixed(:)
                ws_l(:,jk) = fixed(:)     + ws_l(:,jk) * (1._dp-rootfract(:))  
              ELSEWHERE                                                        ! enough soil water
                deficit_l(:) = 0._dp
                ws_l(:,jk) = remaining(:) + ws_l(:,jk) * (1._dp-rootfract(:))
              END WHERE
            ELSEWHERE                                            ! too dry: no water available for transpiration
              deficit_l(:) = trans_layer(:)
            END WHERE

            ! sum up deficits from the different layers
            deficit_trans(:) = deficit_trans(:) + deficit_l(:)

          END WHERE
        END WHERE

      END WHERE
    END DO

    ! -----------------------------------------
    !  reduction needed for evapotranspiration
    ! -----------------------------------------
    WHERE (is_land(:))
      reduce_evap(:) = reduce_evap(:) + deficit_evap(:) + deficit_trans(:)
    END WHERE

    ! --------------------------
    !  soil moisture correction
    ! --------------------------
    ! reduce_evap is only used for diagnostics. If the therm is not zero, soil moisture must be
    ! changed to close the land surface water balance.

    ws(:) = SUM(ws_l(:,:), DIM=2) + SUM(soil_ice(:,:), DIM=2)  ! water and ice  within the column

    DO jk=1, nsoil

      WHERE (reduce_evap(:) < 0._dp .AND. ws(:) > 0._dp)
        ! reduce soil moisture and ice relative to the water content of each layer
        ws_l(:,jk)     = MAX(0._dp, ws_l(:,jk)     + reduce_evap(:) * ws_l(:,jk)/ws(:))
        soil_ice(:,jk) = MAX(0._dp, soil_ice(:,jk) + reduce_evap(:) * soil_ice(:,jk)/ws(:))
      END WHERE
    END DO

    ! make sure there is no negative soil ice
    IF (ANY(soil_ice(:,:) < 0._dp)) THEN
      WRITE (message_text,*) 'negative soil ice: ', MINVAL(soil_ice), ' at ', MINLOC(soil_ice)
      CALL finish ('digest_evapotrans', message_text)
!vg       CALL message ('digest_evapotrans', message_text, all_print=.TRUE., level=em_warn)
!vg       soil_ice(:,:) = MAX(soil_ice(:,:), 0._dp)
    END IF

    ! make sure there is no negative soil moisture and soil moisture does not exceed field capacity
    IF (ANY(ws_l(:,:) < 0._dp)) THEN
      WRITE (message_text,*) 'negative soil moisture: ', MINVAL(ws_l), ' at ', MINLOC(ws_l)
      CALL message ('digest_evapotrans', message_text, all_print=.TRUE., level=em_warn)
      CALL finish ('digest_evapotrans', message_text)
!vg       ws_l(:,:) = MAX(ws_l(:,:), 0._dp)
    END IF
    IF (ANY(ws_l(:,:) > field_cap_l(:,:) + zfcmin)) THEN
      max2d(:) = MAXLOC(ws_l(:,:) - field_cap_l(:,:))
      WRITE (message_text,*) 'liquid water content above field capacity: at ',max2d(:) ,': ws_l: ', ws_l(max2d(1),max2d(2)), &
           ' field_cap_l: ', field_cap_l(max2d(1),max2d(2))
      CALL finish ('digest_evapotrans', message_text)
    END IF 

  END SUBROUTINE digest_evapotrans
  !------------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------------
  SUBROUTINE digest_infilt (nidx, nsoil, is_land, field_cap_l, infilt, ws_l, drainage)
  !------------------------------------------------------------------------------------------------
  !  The routine is based on routine soilchange by Stefan Hagemann. It replaces the part for 
  !  infiltration changes (isch=2)   
  !------------------------------------------------------------------------------------------------

    ! arguments

    INTEGER,  INTENT(in)     :: nidx                    ! vector length
    INTEGER,  INTENT(in)     :: nsoil                   ! number of soil layers
    LOGICAL,  INTENT(in)     :: is_land(nidx)           ! land mask
    REAL(dp), INTENT(in)     :: field_cap_l(nidx,nsoil) ! field capacity of the layer
    REAL(dp), INTENT(in)     :: infilt(nidx)            ! infiltration in [m/s]

    REAL(dp), INTENT(inout)  :: ws_l(nidx,nsoil)        ! water content of the soil layer [m]
    REAL(dp), INTENT(inout)  :: drainage(nidx)          ! drainage [m]

    ! local variables

    INTEGER  :: jk
    REAL(dp) :: excess(nidx)         ! excess water after infiltration

    !------------------------------------------------------------------------------------------------
    !  update soil moisture according to infiltration
    !------------------------------------------------------------------------------------------------

    excess(:) = infilt(:)                   ! amount of water for infiltration at the surface [m]
    DO  jk = 1, nsoil
      WHERE (is_land(:))
        ws_l(:,jk) = ws_l(:,jk) + excess(:)
        WHERE (ws_l(:,jk) > field_cap_l(:,jk))             ! soil moisture exceeds field capacity
          excess(:) = ws_l(:,jk) - field_cap_l(:,jk)       !   excess water -> goes into level below
          ws_l(:,jk) = field_cap_l(:,jk)                   !   limit soil moisture to field capacity
        ELSEWHERE
          excess(:) = 0._dp
        END WHERE
      END WHERE
    END DO
  
    ! if all soil layers are at field capacity and there is still excess water, it is added to drainage
    WHERE (is_land(:) .AND. excess(:) > 0._dp)
      drainage(:) = drainage(:) + excess(:)
    END WHERE

    ! make sure there is no negative soil moisture
    IF (ANY(ws_l(:,:) < 0._dp)) THEN
      WRITE (message_text,*) 'negative soil moisture: ', MINVAL(ws_l), ' at ', MINLOC(ws_l)
      CALL finish ('digest_infilt', message_text)
    END IF

  END SUBROUTINE digest_infilt
  !------------------------------------------------------------------------------------------------

END MODULE
