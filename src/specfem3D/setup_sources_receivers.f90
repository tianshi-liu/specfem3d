!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================


  subroutine setup_sources_receivers()

  use specfem_par

  implicit none

  ! checks if anything to do
  if (.not. IO_compute_task) return

  ! setup for point search
  call setup_point_search_arrays()

  ! sets number of timestep parameters
  call setup_timesteps()

  ! locates sources and determines simulation start time t0
  call setup_sources()

  ! reads in stations file and locates receivers
  call setup_receivers()

  ! pre-compute source arrays
  call setup_sources_precompute_arrays()

  ! pre-compute receiver interpolation factors
  call setup_receivers_precompute_intp()

  ! write source and receiver VTK files for Paraview
  call setup_sources_receivers_VTKfile()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Total number of samples for seismograms = ',NSTEP
    write(IMAIN,*)
    if (NSOURCES > 1) then
      write(IMAIN,*) 'Using ',NSOURCES,' point sources'
      write(IMAIN,*)
    endif
    call flush_IMAIN()
  endif

  ! frees memory
  if (.not. INVERSE_FWI_FULL_PROBLEM) then
    deallocate(xyz_midpoints)
    deallocate(anchor_iax,anchor_iay,anchor_iaz)
    ! no need to keep mesh adjacency after point searches
    deallocate(neighbors_xadj,neighbors_adjncy)
    if (.not. DO_BRUTE_FORCE_POINT_SEARCH) call setup_free_kdtree()
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine setup_sources_receivers

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_point_search_arrays()

  use constants
  use specfem_par

  implicit none

  ! local parameters
  integer :: ispec,iglob,ier

  ! prepares midpoints coordinates
  allocate(xyz_midpoints(NDIM,NSPEC_AB),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating array xyz_midpoints')
  xyz_midpoints(:,:) = 0.d0

  ! store x/y/z coordinates of center point
  do ispec = 1,NSPEC_AB
    iglob = ibool(MIDX,MIDY,MIDZ,ispec)
    xyz_midpoints(1,ispec) =  dble(xstore(iglob))
    xyz_midpoints(2,ispec) =  dble(ystore(iglob))
    xyz_midpoints(3,ispec) =  dble(zstore(iglob))
  enddo

  ! define (i,j,k) indices of the control/anchor points
  allocate(anchor_iax(NGNOD),anchor_iay(NGNOD),anchor_iaz(NGNOD),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating array anchor_i**')
  anchor_iax(:) = 0; anchor_iay(:) = 0; anchor_iaz(:) = 0

  call hex_nodes_anchor_ijk_NGLL(NGNOD,anchor_iax,anchor_iay,anchor_iaz,NGLLX,NGLLY,NGLLZ)

  ! builds search tree
  if (.not. DO_BRUTE_FORCE_POINT_SEARCH) call setup_search_kdtree()

  end subroutine setup_point_search_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_timesteps()

  use specfem_par, only: myrank,NSTEP,SIMULATION_TYPE, &
    NTSTEP_BETWEEN_OUTPUT_SEISMOS,NTSTEP_BETWEEN_READ_ADJSRC, &
    SAVE_ALL_SEISMOS_IN_ONE_FILE, &
    ASDF_FORMAT,HDF5_FORMAT

  implicit none

  ! subsets used to save seismograms must not be larger than the whole time series
  if (NTSTEP_BETWEEN_OUTPUT_SEISMOS > NSTEP) NTSTEP_BETWEEN_OUTPUT_SEISMOS = NSTEP

  ! subsets used to save adjoint sources must not be larger than the whole time series,
  ! otherwise we waste memory
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! limits length
    if (NTSTEP_BETWEEN_READ_ADJSRC > NSTEP) NTSTEP_BETWEEN_READ_ADJSRC = NSTEP
  endif

  ! check
  if (SAVE_ALL_SEISMOS_IN_ONE_FILE) then
    ! requires to have full length of seismograms
    if (NTSTEP_BETWEEN_OUTPUT_SEISMOS < NSTEP) then
      print *, 'Error: Setting SAVE_ALL_SEISMOS_IN_ONE_FILE to .true. requires NTSTEP_BETWEEN_OUTPUT_SEISMOS >= NSTEP'
      call exit_MPI(myrank, &
                    'Setting SAVE_ALL_SEISMOS_IN_ONE_FILE not supported for current NTSTEP_BETWEEN_OUTPUT_SEISMOS in Par_file')
    endif
  endif
  if (ASDF_FORMAT) then
    ! ASDF storage requires to have full length of seismograms
    if (NTSTEP_BETWEEN_OUTPUT_SEISMOS < NSTEP) then
      print *, 'Error: Setting ASDF_FORMAT to .true. requires NTSTEP_BETWEEN_OUTPUT_SEISMOS >= NSTEP'
      call exit_MPI(myrank,'Error: Setting ASDF_FORMAT to .true. requires NTSTEP_BETWEEN_OUTPUT_SEISMOS >= NSTEP')
    endif
  endif
  if (HDF5_FORMAT) then
    ! HDF5 storage requires to have full length of seismograms
    if (NTSTEP_BETWEEN_OUTPUT_SEISMOS < NSTEP) then
      print *, 'Error: Setting HDF5_FORMAT to .true. requires NTSTEP_BETWEEN_OUTPUT_SEISMOS >= NSTEP'
      call exit_MPI(myrank,'Error: Setting HDF5_FORMAT to .true. requires NTSTEP_BETWEEN_OUTPUT_SEISMOS >= NSTEP')
    endif
  endif

  end subroutine setup_timesteps

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources()

  use specfem_par

  implicit none

  ! local parameters
  integer :: isource,ier

  ! initializes onset time (depends on hdur and source time function)
  t0 = 0.d0

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'sources: ',NSOURCES
    ! finite fault simulations ignore CMT and force sources
    if (HAS_FINITE_FAULT_SOURCE) &
      write(IMAIN,*) '         finite fault simulation'
    call flush_IMAIN()
  endif

  ! allocate arrays for source
  allocate(islice_selected_source(NSOURCES), &
           ispec_selected_source(NSOURCES), &
           Mxx(NSOURCES), &
           Myy(NSOURCES), &
           Mzz(NSOURCES), &
           Mxy(NSOURCES), &
           Mxz(NSOURCES), &
           Myz(NSOURCES), &
           xi_source(NSOURCES), &
           eta_source(NSOURCES), &
           gamma_source(NSOURCES), &
           tshift_src(NSOURCES), &
           hdur(NSOURCES), &
           hdur_Gaussian(NSOURCES), &
           utm_x_source(NSOURCES), &
           utm_y_source(NSOURCES), &
           nu_source(NDIM,NDIM,NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2049')
  if (ier /= 0) stop 'error allocating arrays for sources'

  ! initializes arrays
  islice_selected_source(:) = -1
  ispec_selected_source(:) = 0
  Mxx(:) = 0.d0; Myy(:) = 0.d0; Mzz(:) = 0.d0
  Mxy(:) = 0.d0; Mxz(:) = 0.d0; Myz(:) = 0.d0
  xi_source(:) = 0.d0; eta_source(:) = 0.d0; gamma_source(:) = 0.d0
  tshift_src(:) = 0.d0; hdur(:) = 0.d0; hdur_Gaussian(:) = 0.d0
  utm_x_source(:) = 0.d0; utm_y_source(:) = 0.d0
  nu_source(:,:,:) = 0.d0

  if (USE_FORCE_POINT_SOURCE) then
    allocate(force_stf(NSOURCES), &
             factor_force_source(NSOURCES), &
             comp_dir_vect_source_E(NSOURCES), &
             comp_dir_vect_source_N(NSOURCES), &
             comp_dir_vect_source_Z_UP(NSOURCES),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2053')
  else
    allocate(force_stf(1), &
             factor_force_source(1), &
             comp_dir_vect_source_E(1), &
             comp_dir_vect_source_N(1), &
             comp_dir_vect_source_Z_UP(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2057')
  endif
  if (ier /= 0) stop 'error allocating arrays for force point sources'
  force_stf(:) = 0
  factor_force_source(:) = 0.d0
  comp_dir_vect_source_E(:) = 0.d0
  comp_dir_vect_source_N(:) = 0.d0
  comp_dir_vect_source_Z_UP(:) = 0.d0

  ! sets the size of user_source_time_function array
  if (USE_EXTERNAL_SOURCE_FILE) then
    NSTEP_STF = NSTEP
    NSOURCES_STF = NSOURCES
  else
    ! We don't need the array user_source_time_function : use a small dummy array
    NSTEP_STF = 1
    NSOURCES_STF = 1
  endif
  ! allocate array that contains the user defined source time function
  allocate(user_source_time_function(NSTEP_STF, NSOURCES_STF),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2058')
  if (ier /= 0) stop 'error allocating arrays for user sources time function'
  user_source_time_function(:,:) = 0.0_CUSTOM_REAL

  ! fused wavefield simulations
  call get_run_number_of_the_source()

  ! checks if anything left to do
  if (INVERSE_FWI_FULL_PROBLEM) then
    ! fwi will determine acquisition receivers later
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'running inverse FWI full problem, will determine sources later...'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
    ! all done
    return
  endif

  ! locate sources in the mesh
  call locate_source()

  ! count number of sources located in this slice
  nsources_local = 0
  do isource = 1, NSOURCES
    if (myrank == islice_selected_source(isource)) nsources_local = nsources_local + 1
  enddo

  ! checks if source is in an acoustic element and exactly on the free surface because pressure is zero there
  call setup_sources_check_acoustic()

  ! determines onset time
  call setup_stf_constants()

  ! prints source time functions to output files
  if (PRINT_SOURCE_TIME_FUNCTION) call print_stf_file()

  end subroutine setup_sources
!
!-------------------------------------------------------------------------------------------------
!
! When NB_RUNS_ACOUSTIC_GPU > 1, the source file CMTSOLUTION (or FORCESOLUTION) actually contains the sources for all the runs.
! This routine is intended to get the array that contains the run number of each source described in the source file.
! The line i of the file run_number_of_the_source contains the run number \in [ 0;NB_RUNS_ACOUSTIC_GPU-1] of the source i

  subroutine get_run_number_of_the_source()

  use constants
  use specfem_par, only: run_number_of_the_source,NSOURCES

  implicit none
  character(len=MAX_STRING_LEN) :: filename,string
  integer :: ier,isource,icounter,id_run

  allocate(run_number_of_the_source(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2059')
  run_number_of_the_source(:) = 0

  if (NB_RUNS_ACOUSTIC_GPU > 1) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'fused wavefield simulation:'
      write(IMAIN,*) '  NB_RUNS_ACOUSTIC_GPU = ',NB_RUNS_ACOUSTIC_GPU
      write(IMAIN,*) '  NSOURCES             = ',NSOURCES
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! reads file DATA/run_number_of_the_source
    filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'run_number_of_the_source'
    open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error opening run_number_of_the_source file'
    endif

    ! Checks if the number of lines is correct
    icounter = 0
    do while (ier == 0)
      read(IIN,"(a)",iostat=ier) string
      if (ier == 0) icounter = icounter + 1
    enddo
    close(IIN)
    if (icounter /= NSOURCES) stop 'Error total number of lines in run_number_of_the_source file is not equal to NSOURCES'

    ! Fills the array run_number_of_the_source
    open(unit=IIN,file=trim(filename),status='old',action='read')
    ! reads run number for each source
    do isource = 1,NSOURCES
      read(IIN,"(a)") string
      ! reads run id: for each source entry in CMTSOLUTION/FORCESOLUTION a corresponding line in run_number_of_the_source
      !               with a run id must be given. the run id must be between 0 and NB_RUNS_ACOUSTIC_GPU-1
      read(string,*) id_run
      ! checks that id between 0 and NB_RUNS_ACOUSTIC_GPU - 1
      if (id_run < 0 .or. id_run >= NB_RUNS_ACOUSTIC_GPU) then
        print *,'Error: run id entry in run_number_of_the_source file must be between 0 and NB_RUNS_ACOUSTIC_GPU-1'
        print *,'    setting: NB_RUNS_ACOUSTIC_GPU = ',NB_RUNS_ACOUSTIC_GPU
        print *,'             file line entry ',isource,' has invalid run id ',id_run
        stop 'Error invalid run id for source line in run_number_of_the_source file'
      endif
      ! sets source id for run
      run_number_of_the_source(isource) = id_run
    enddo

    ! user output
    if (myrank == 0) then
      do isource = 1,NSOURCES
        write(IMAIN,*) '  source ',isource,' assigned to wavefield run number ',run_number_of_the_source(isource)
      enddo
      write(IMAIN,*)
      ! warning
      if (NSOURCES /= NB_RUNS_ACOUSTIC_GPU) then
        write(IMAIN,*) '  ***'
        write(IMAIN,*) '  *** WARNING: number of sources ',NSOURCES, &
                       ' does not match number of runs ',NB_RUNS_ACOUSTIC_GPU,' ***'
        write(IMAIN,*) '  ***'
        write(IMAIN,*)
      endif
      call flush_IMAIN()
    endif
  endif ! NB_RUNS_ACOUSTIC_GPU

  end subroutine get_run_number_of_the_source

!
!-------------------------------------------------------------------------------------------------
!
  subroutine setup_stf_constants()

  use specfem_par
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_movie

  implicit none

  ! local parameters
  double precision :: t0_acoustic
  integer :: isource,ispec

  ! initializes simulation start time t0
  t0 = 0.d0

  ! checks if anything to do, finite fault simulations ignore CMT and force sources
  if (HAS_FINITE_FAULT_SOURCE) return

  if (abs(minval(tshift_src)) > TINYVAL) then
    ! this should be a warning; users can decide to do this purposely
    ! (in particular for applications outside of seismology, e.g. in imaging or in non-destructive testing)
    write(IMAIN,*)
    write(IMAIN,*) 'INFORMATION: no tshift_src is equal to zero, thus the origin time is not t0, it is changed by tshift_src'
    write(IMAIN,*)
  endif

  ! filter source time function by Gaussian with hdur = HDUR_MOVIE when outputing movies or shakemaps
  if (MOVIE_SURFACE .or. MOVIE_VOLUME .or. CREATE_SHAKEMAP) then
    hdur = sqrt(hdur**2 + HDUR_MOVIE**2)
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Each source is being convolved with HDUR_MOVIE = ',HDUR_MOVIE
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! one cannot use a Heaviside source for the movies
    ! checks if CMT source time function is a Heaviside
    ! (high-frequency oscillations don't look good in movies)
    if (.not. USE_FORCE_POINT_SOURCE .and. &
        .not. USE_EXTERNAL_SOURCE_FILE) then
      if (minval(hdur(:)) < TINYVAL) &
        stop 'Error hdur too small for movie creation, movies do not make sense for Heaviside source'
    endif
  endif

  ! convert the half duration for triangle STF to the one for Gaussian STF
  hdur_Gaussian(:) = hdur(:) / SOURCE_DECAY_MIMIC_TRIANGLE

  ! define t0 as the earliest start time
  if (USE_FORCE_POINT_SOURCE) then
    ! point force sources
    ! (might start depending on the frequency given by hdur)
    ! note: point force sources will give the dominant frequency in hdur, thus the main period is 1/hdur.
    !       also, these sources might use a Ricker source time function instead of a Gaussian.
    !       For a Ricker source time function, a start time ~1.2 * main_period is a good choice.
    t0 = 0.d0
    do isource = 1,NSOURCES
      select case(force_stf(isource))
      case (0)
        ! Gaussian source time function
        t0 = min(t0,1.5d0 * (tshift_src(isource) - hdur(isource)))
      case (1)
        ! Ricker source time function
        t0 = min(t0,1.2d0 * (tshift_src(isource) - 1.0d0/hdur(isource)))
      case (2)
        ! Heaviside
        t0 = min(t0,1.5d0 * (tshift_src(isource) - hdur(isource)))
      case (3)
        ! Monochromatic
        t0 = 0.d0
      case (4)
        ! Gaussian source time function by Meschede et al. (2011)
        t0 = min(t0,1.5d0 * (tshift_src(isource) - hdur(isource)))
      case (5)
        ! Brune
        ! This needs to be CHECKED!!!
        t0 = min(t0,1.5d0 * (tshift_src(isource) - hdur(isource)))
      case (6)
        ! Smoothed Brune
        ! This needs to be CHECKED!!!
        t0 = min(t0,1.5d0 * (tshift_src(isource) - hdur(isource)))
      case default
        stop 'unsupported force_stf value!'
      end select
    enddo
    ! start time defined as positive value, will be subtracted
    t0 = - t0
  else
    ! moment tensors
    if (USE_RICKER_TIME_FUNCTION) then
      ! note: sources will give the dominant frequency in hdur,
      !       thus the main period is 1/hdur.
      !       for a Ricker source time function, a start time ~1.2 * dominant_period is a good choice
      t0 = - 1.2d0 * minval(tshift_src(:) - 1.0d0/hdur(:))
    else
      ! (based on Heaviside functions)
      ! note: an earlier start time also reduces numerical noise due to a
      !          non-zero offset at the beginning of the source time function
      t0 = - 2.0d0 * minval(tshift_src(:) - hdur(:))   ! - 1.5d0 * minval(tshift_src-hdur)
    endif
  endif

  ! uses an earlier start time if source is acoustic with a Gaussian source time function
  t0_acoustic = 0.0d0
  do isource = 1,NSOURCES
    if (myrank == islice_selected_source(isource)) then
      ispec = ispec_selected_source(isource)
      if (ispec_is_acoustic(ispec)) then
        ! uses an earlier start time
        if (USE_FORCE_POINT_SOURCE) then
          if (force_stf(isource) == 0) then
            ! Gaussian
            t0_acoustic = - 3.0d0 * ( tshift_src(isource) - hdur(isource) )
          endif
        else
          if (.not. USE_RICKER_TIME_FUNCTION) then
            ! Gaussian STF by default
            t0_acoustic = - 3.0d0 * ( tshift_src(isource) - hdur(isource) )
          endif
        endif
        if (t0_acoustic > t0 ) t0 = t0_acoustic
      endif
    endif
  enddo
  ! passes maximum value to all processes
  ! note: t0 is defined positive and will be subtracted from simulation time (it-1)*DT
  t0_acoustic = t0
  call max_all_all_dp(t0_acoustic,t0)

  !! VM VM for external source the time will begin with simulation
  if (USE_EXTERNAL_SOURCE_FILE) then
    t0 = 0.d0
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'External STF:'
      write(IMAIN,*) '  simulation start time set to zero: ', t0
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! checks if user set USER_T0 to fix simulation start time
  ! note: USER_T0 has to be positive
  if (USER_T0 > 0.d0) then
    ! user cares about origin time and time shifts of the CMTSOLUTION
    ! and wants to fix simulation start time to a constant start time
    ! time 0 on time axis will correspond to given origin time

    ! notifies user
    if (myrank == 0) then
      write(IMAIN,*) 'USER_T0: ',USER_T0
      write(IMAIN,*) 't0: ',t0,'min_tshift_src_original: ',min_tshift_src_original
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! checks if automatically set t0 is too small
    ! note: min_tshift_src_original can be a positive or negative time shift (minimum from all tshift)
    if (t0 <= USER_T0 + min_tshift_src_original) then
      ! by default, tshift_src(:) holds relative time shifts with a minimum time shift set to zero
      ! re-adds (minimum) original time shift such that sources will kick in
      ! according to their absolute time shift
      tshift_src(:) = tshift_src(:) + min_tshift_src_original

      ! sets new simulation start time such that
      ! simulation starts at t = - t0 = - USER_T0
      t0 = USER_T0

      ! notifies user
      if (myrank == 0) then
        write(IMAIN,*) '  set new simulation start time: ', - t0
        write(IMAIN,*)
        call flush_IMAIN()
      endif
    else
      ! start time needs to be at least t0 for numerical stability
      ! notifies user
      if (myrank == 0) then
        write(IMAIN,*) 'error: USER_T0 is too small'
        write(IMAIN,*) '       must make one of three adjustements:'
        write(IMAIN,*) '       - increase USER_T0 to be at least: ',t0-min_tshift_src_original
        write(IMAIN,*) '       - decrease time shift in CMTSOLUTION file'
        write(IMAIN,*) '       - decrease hdur in CMTSOLUTION file'
        call flush_IMAIN()
      endif
      call exit_mpi(myrank,'error USER_T0 is set but too small')
    endif
  else if (USER_T0 < 0.d0) then
    if (myrank == 0) then
      write(IMAIN,*) 'error: USER_T0 is negative, must be set zero or positive!'
      call flush_IMAIN()
    endif
    call exit_mpi(myrank,'error negative USER_T0 parameter in constants.h')
  endif

  end subroutine setup_stf_constants

!
!-------------------------------------------------------------------------------------------------
!


  subroutine setup_sources_check_acoustic()

! checks if source is in an acoustic element and exactly on the free surface because pressure is zero there

  use specfem_par
  use specfem_par_acoustic
  implicit none

  integer :: isource,ixmin,ixmax,iymin,iymax,izmin,izmax,iface,ispec
  logical :: is_on,is_on_all

  ! outputs a warning in case of an acoustic source lying on the free surface
  do isource = 1,NSOURCES
    ! checks if source is close to face
    is_on = .false.

    ! only receivers in this process
    if (myrank == islice_selected_source(isource)) then

      ispec = ispec_selected_source(isource)
      ! only if receiver is in an acoustic element
      if (ispec_is_acoustic(ispec)) then

        ! checks with free surface face
        do iface = 1,num_free_surface_faces

          if (ispec == free_surface_ispec(iface)) then

            ! determine face
            ixmin = minval( free_surface_ijk(1,:,iface) )
            ixmax = maxval( free_surface_ijk(1,:,iface) )

            iymin = minval( free_surface_ijk(2,:,iface) )
            iymax = maxval( free_surface_ijk(2,:,iface) )

            izmin = minval( free_surface_ijk(3,:,iface) )
            izmax = maxval( free_surface_ijk(3,:,iface) )

            ! xmin face
            if (ixmin == 1 .and. ixmax == 1) then
              if (xi_source(isource) < -0.99d0) is_on = .true.
            ! xmax face
            else if (ixmin == NGLLX .and. ixmax == NGLLX) then
              if (xi_source(isource) > 0.99d0) is_on = .true.
            ! ymin face
            else if (iymin == 1 .and. iymax == 1) then
              if (eta_source(isource) < -0.99d0) is_on = .true.
            ! ymax face
            else if (iymin == NGLLY .and. iymax == NGLLY) then
              if (eta_source(isource) > 0.99d0) is_on = .true.
            ! zmin face
            else if (izmin == 1 .and. izmax == 1) then
              if (gamma_source(isource) < -0.99d0) is_on = .true.
            ! zmax face
            else if (izmin == NGLLZ .and. izmax == NGLLZ) then
              if (gamma_source(isource) > 0.99d0) is_on = .true.
            endif

          endif ! free_surface_ispec
        enddo ! iface
      endif ! ispec_is_acoustic
    endif ! islice_selected_rec

    ! user output
    call any_all_l( is_on, is_on_all )
    if (myrank == 0 .and. is_on_all) then
      write(IMAIN,*) '********************************************************************'
      write(IMAIN,*) '*** source: ',isource,'                                          ***'
      write(IMAIN,*) '*** Warning: acoustic source located exactly on the free surface ***'
      write(IMAIN,*) '*** will be zeroed                                               ***'
      write(IMAIN,*) '********************************************************************'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  enddo ! num_free_surface_faces

  end subroutine setup_sources_check_acoustic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_receivers()

  use specfem_par
  use specfem_par_acoustic

  implicit none

  ! local parameters
  integer :: nrec_simulation,nrec_filtered
  integer :: nrec_tot_found
  integer :: irec,isource,ier
  character(len=MAX_STRING_LEN) :: rec_filename,filtered_rec_filename
  character(len=MAX_STRING_LEN) :: path_to_add

  ! adjoint sources/receivers
  integer :: icomp,itime,nadj_files_found,nadj_files_found_tot
  real(kind=CUSTOM_REAL) :: junk
  character(len=3),dimension(NDIM) :: comp
  character(len=MAX_STRING_LEN) :: filename
  character(len=MAX_STRING_LEN) :: adj_source_file

  ! initializes
  nrec = 0
  nrec_filtered = 0
  nrec_local = 0

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'receivers:'
    call flush_IMAIN()
  endif

  ! gets number of stations
  if (INVERSE_FWI_FULL_PROBLEM) then
    ! stations will be determined later...
    nrec = 0
    nrec_local = 0
  else
    ! reads in station file
    if (SIMULATION_TYPE == 1) then
      rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS'
      filtered_rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_FILTERED'
    else
      rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_ADJOINT'
      filtered_rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_ADJOINT_FILTERED'
    endif

    ! see if we are running several independent runs in parallel
    ! if so, add the right directory for that run
    ! (group numbers start at zero, but directory names start at run0001, thus we add one)
    ! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      rec_filename = path_to_add(1:len_trim(path_to_add))//rec_filename(1:len_trim(rec_filename))
      filtered_rec_filename = path_to_add(1:len_trim(path_to_add))//filtered_rec_filename(1:len_trim(filtered_rec_filename))
    endif

    call station_filter(rec_filename,filtered_rec_filename,nrec_filtered)

    ! sets actual number of stations
    nrec = nrec_filtered

    if (nrec < 1) call exit_MPI(myrank,'need at least one receiver')
  endif

  call synchronize_all()

  ! user output
  if (myrank == 0) then
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      write(IMAIN,*) '  Total number of receivers = ', nrec
    else
      write(IMAIN,*) '  Total number of adjoint sources = ', nrec
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! allocate memory for receiver arrays, i.e. stations given in STATIONS file
  allocate(islice_selected_rec(nrec), &
           ispec_selected_rec(nrec), &
           xi_receiver(nrec), &
           eta_receiver(nrec), &
           gamma_receiver(nrec), &
           nu_rec(NDIM,NDIM,nrec), stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2067')
  if (ier /= 0) stop 'error allocating arrays for receivers'
  ! initializes arrays
  islice_selected_rec(:) = -1
  ispec_selected_rec(:) = 0
  xi_receiver(:) = 0.d0; eta_receiver(:) = 0.d0; gamma_receiver(:) = 0.d0
  nu_rec(:,:,:) = 0.0d0

  allocate(station_name(nrec), &
           network_name(nrec), &
           stlat(nrec), &
           stlon(nrec), &
           stele(nrec), &
           stbur(nrec),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating receiver arrays')
  stlat(:) = 0.d0; stlon(:) = 0.d0; stele(:) = 0.d0; stbur(:) = 0.d0
  station_name(:) = ""; network_name(:) = ""

  ! checks if anything left to do
  if (INVERSE_FWI_FULL_PROBLEM) then
    ! fwi will determine acquisition receivers later
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'running inverse FWI full problem, will determine receivers later...'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
    ! all done
    return
  endif

  ! locate receivers in the mesh
  call locate_receivers(filtered_rec_filename,utm_x_source(1),utm_y_source(1))

  ! count number of receivers located in this slice
  nrec_local = 0
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    ! number of receivers are given by stations
    ! in STATIONS (forward runs) or STATIONS_ADJOINT (kernel runs) file
    nrec_simulation = nrec
    do irec = 1,nrec
      if (myrank == islice_selected_rec(irec)) nrec_local = nrec_local + 1
    enddo
  else
    ! adjoint simulation:
    ! station locations (in STATIONS_ADJOINT file) become adjoint sources
    ! and source locations (in CMTSOLUTION file) become adjoint "receivers"
    nrec_simulation = NSOURCES
    do isource = 1, NSOURCES
      if (myrank == islice_selected_source(isource)) nrec_local = nrec_local + 1
    enddo
  endif

  ! check that the sum of the number of receivers in each slice is nrec
  call sum_all_i(nrec_local,nrec_tot_found)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'found a total of ',nrec_tot_found,' receivers in all the slices'
    write(IMAIN,*)
  endif

  ! checks
  if (myrank == 0) then
    if (nrec_tot_found /= nrec_simulation) then
      call exit_MPI(myrank,'problem when dispatching the receivers')
    endif
  endif
  call synchronize_all()

  ! checks if acoustic receiver is exactly on the free surface because pressure is zero there
  call setup_receivers_check_acoustic()

  ! counter for adjoint receiver stations in local slice, used to allocate adjoint source arrays
  nadj_rec_local = 0

  ! counts receivers for adjoint simulations
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! ADJOINT simulations
    ! gets number of local adjoint sources, i.e. located in this slice (nadj_rec_local)
    if (.not. SU_FORMAT) then
      ! prepares channel names
      do icomp = 1,NDIM
        call write_channel_name(icomp,comp(icomp))
      enddo

      ! temporary counter to check if any files are found at all
      nadj_files_found = 0

      do irec = 1,nrec
        ! adjoint receiver station
        if (myrank == islice_selected_rec(irec)) then
          ! checks that the source slice number is okay
          if (islice_selected_rec(irec) < 0 .or. islice_selected_rec(irec) > NPROC-1) &
            call exit_MPI(myrank,'something is wrong with the source slice number in adjoint simulation')

          ! updates counter
          nadj_rec_local = nadj_rec_local + 1

          ! checks **net**.**sta**.**MX**.adj files for correct number of time steps
          if (READ_ADJSRC_ASDF) then
            ! ASDF format
            call check_adjoint_sources_asdf(irec,nadj_files_found)
          else
            ! ASCII format
            adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))

            ! loops over file components E/N/Z
            do icomp = 1,NDIM
              ! trace name
              ! format: **net**.**sta**.**comp**.adj
              !         example: DB.X01.HXX.adj
              filename = trim(adj_source_file) // '.' // comp(icomp) // '.adj'

              ! adds path to SEM/ folder relative to OUTPUT_FILES/
              filename = OUTPUT_FILES(1:len_trim(OUTPUT_FILES))// '/../SEM/' // trim(filename)

              ! reads in adjoint trace
              open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)
              if (ier == 0) then
                ! checks length of file
                itime = 0
                do while (ier == 0)
                  read(IIN,*,iostat=ier) junk,junk
                  if (ier == 0) itime = itime + 1
                enddo

                ! checks length
                if (itime /= NSTEP) then
                  print *,'adjoint source error: ',trim(filename),' has length',itime,' but should be',NSTEP
                  call exit_MPI(myrank, &
                    'file '//trim(filename)//' has wrong length, please check your adjoint sources and your simulation duration')
                endif

                ! updates counter for found files
                nadj_files_found = nadj_files_found + 1
              endif
              ! closes file
              close(IIN)
            enddo
          endif
        endif
      enddo

    else
      ! SU_FORMAT file
      ! adjoint sources given in single SU_FORMAT file;
      ! skip counting, because only one file per component per proc in SU_FORMAT
      nadj_rec_local = nrec_local
      nadj_files_found = nrec_local
    endif !if (.not. SU_FORMAT)

    ! checks if any adjoint source files found at all
    call sum_all_i(nadj_files_found,nadj_files_found_tot)
    if (myrank == 0) then
      ! user output
      write(IMAIN,*) '    ',nadj_files_found_tot,' adjoint component trace files found in all slices'
      write(IMAIN,*)
      call flush_IMAIN()

      ! main process checks if any adjoint files found
      if (nadj_files_found_tot == 0) then
        print *,'Error no adjoint traces found: ',nadj_files_found_tot
        print *,'in directory : ',OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//'/../SEM/'
        if (.not. SU_FORMAT .and. .not. READ_ADJSRC_ASDF) then
          print *,'with endings : ', '**.'//comp(1)//'.adj',' ','**.'//comp(2)//'.adj',' ','**.'//comp(3)//'.adj'
        endif
        print *
        call exit_MPI(myrank,'no adjoint traces found, please check adjoint sources in directory SEM/')
      endif
    endif

    ! initializes adjoint sources
    allocate(source_adjoint(NDIM,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2073')
    if (ier /= 0) stop 'error allocating array source_adjoint'
    source_adjoint(:,:,:) = 0.0_CUSTOM_REAL

    ! note:
    ! computes adjoint sources in chunks/blocks during time iterations.
    ! we moved it to compute_add_sources_viscoelastic.f90 & compute_add_sources_acoustic.f90,
    ! because we may need to read in adjoint sources block by block
  endif

  end subroutine setup_receivers


!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_receivers_check_acoustic()

! checks if acoustic receiver is exactly on the free surface because pressure is zero there

  use specfem_par
  use specfem_par_acoustic
  implicit none

  integer :: irec,ixmin,ixmax,iymin,iymax,izmin,izmax,iface,ispec
  logical :: is_on,is_on_all

  ! outputs a warning in case the receiver is lying on the free surface
  do irec = 1,nrec

    ! checks if receiver is close to face
    is_on = .false.

    ! only receivers in this process
    if (myrank == islice_selected_rec(irec)) then

      ispec = ispec_selected_rec(irec)
      ! only if receiver is in an acoustic element
      if (ispec_is_acoustic(ispec)) then

        ! checks with free surface face
        do iface = 1,num_free_surface_faces

          if (ispec == free_surface_ispec(iface)) then

            ! determine face
            ixmin = minval( free_surface_ijk(1,:,iface) )
            ixmax = maxval( free_surface_ijk(1,:,iface) )

            iymin = minval( free_surface_ijk(2,:,iface) )
            iymax = maxval( free_surface_ijk(2,:,iface) )

            izmin = minval( free_surface_ijk(3,:,iface) )
            izmax = maxval( free_surface_ijk(3,:,iface) )

            ! xmin face
            if (ixmin == 1 .and. ixmax == 1) then
              if (xi_receiver(irec) < -0.99d0) is_on = .true.
            ! xmax face
            else if (ixmin == NGLLX .and. ixmax == NGLLX) then
              if (xi_receiver(irec) > 0.99d0) is_on = .true.
            ! ymin face
            else if (iymin == 1 .and. iymax == 1) then
              if (eta_receiver(irec) < -0.99d0) is_on = .true.
            ! ymax face
            else if (iymin == NGLLY .and. iymax == NGLLY) then
              if (eta_receiver(irec) > 0.99d0) is_on = .true.
            ! zmin face
            else if (izmin == 1 .and. izmax == 1) then
              if (gamma_receiver(irec) < -0.99d0) is_on = .true.
            ! zmax face
            else if (izmin == NGLLZ .and. izmax == NGLLZ) then
              if (gamma_receiver(irec) > 0.99d0) is_on = .true.
            endif

          endif ! free_surface_ispec
        enddo ! iface
      endif ! ispec_is_acoustic
    endif ! islice_selected_rec

    ! user output
    call any_all_l( is_on, is_on_all )
    if (myrank == 0 .and. is_on_all) then
      ! limits user output if too many receivers
      if (nrec < 1000 .and. (.not. SU_FORMAT )) then
        write(IMAIN,*) '**********************************************************************'
        write(IMAIN,*) '*** station:',irec,'                                               ***'
        write(IMAIN,*) '*** Warning: acoustic receiver located exactly on the free surface ***'
        write(IMAIN,*) '*** Warning: tangential component will be zero there               ***'
        write(IMAIN,*) '**********************************************************************'
        write(IMAIN,*)
        call flush_IMAIN()
      endif
    endif

  enddo ! num_free_surface_faces

  end subroutine setup_receivers_check_acoustic


!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_precompute_arrays()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_acoustic

  implicit none

  integer :: isource,ispec,ier

  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas
  double precision :: norm,comp_x,comp_y,comp_z
  double precision :: factor_source

  double precision :: sizeval
  logical :: does_source_encoding

  ! user output info
  ! sources
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    ! user output
    if (myrank == 0) then
      ! note: all process allocate the full sourcearrays array
      ! sourcearrays(NDIM,NGLLX,NGLLY,NGLLZ,NSOURCES)
      sizeval = dble(NSOURCES) * dble(NDIM * NGLLX * NGLLY * NGLLZ * CUSTOM_REAL / 1024. / 1024. )
      ! outputs info
      write(IMAIN,*) 'source arrays:'
      write(IMAIN,*) '  number of sources is ',NSOURCES
      write(IMAIN,*) '  size of source array                 = ', sngl(sizeval),'MB'
      write(IMAIN,*) '                                       = ', sngl(sizeval/1024.d0),'GB'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! for source encoding (acoustic sources only so far)
  if (USE_SOURCE_ENCODING) then
    allocate(pm1_source_encoding(NSOURCES),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2068')
  else
    allocate(pm1_source_encoding(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2069')
  endif
  if (ier /= 0) stop 'error allocating arrays for sources'
  pm1_source_encoding(:) = 1._CUSTOM_REAL
  does_source_encoding = .false.

  ! forward simulations
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    allocate(sourcearray(NDIM,NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2070')
    allocate(sourcearrays(NSOURCES,NDIM,NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2071')
    if (ier /= 0) stop 'error allocating array sourcearray'

    ! compute source arrays
    do isource = 1,NSOURCES

      ! initializes
      sourcearray(:,:,:,:) = ZERO

      ! check that the source slice number is okay
      if (islice_selected_source(isource) < 0 .or. islice_selected_source(isource) > NPROC-1) &
        call exit_MPI(myrank,'something is wrong with the source slice number')

      ! compute source arrays in source slice
      if (myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        ! compute Lagrange polynomials at the source location
        call lagrange_any(xi_source(isource),NGLLX,xigll,hxis,hpxis)
        call lagrange_any(eta_source(isource),NGLLY,yigll,hetas,hpetas)
        call lagrange_any(gamma_source(isource),NGLLZ,zigll,hgammas,hpgammas)

        if (USE_FORCE_POINT_SOURCE) then
          ! use of FORCESOLUTION files
          !
          ! note: for use_force_point_source xi/eta/gamma are also in the range [-1,1], for exact positioning
          factor_source = factor_force_source(isource)

          ! elastic sources
          ! note: point forces in poroelastic element act on both, solid and fluid parts (see compute_add_sources_poroelastic)
          !       we allow thus the force to have a given direction.
          if (ispec_is_elastic(ispec) .or. ispec_is_poroelastic(ispec)) then
            ! length of component vector
            norm = dsqrt( comp_dir_vect_source_E(isource)**2 &
                        + comp_dir_vect_source_N(isource)**2 &
                        + comp_dir_vect_source_Z_UP(isource)**2 )
            ! checks norm of component vector
            if (norm < TINYVAL) then
              call exit_MPI(myrank,'error force point source: component vector has (almost) zero norm')
            endif
            ! normalizes given component vector
            comp_x = comp_dir_vect_source_E(isource)/norm
            comp_y = comp_dir_vect_source_N(isource)/norm
            comp_z = comp_dir_vect_source_Z_UP(isource)/norm

            call compute_arrays_source_forcesolution(sourcearray,hxis,hetas,hgammas,factor_source, &
                                                     comp_x,comp_y,comp_z,nu_source(:,:,isource))
          endif

          ! acoustic sources
          if (ispec_is_acoustic(ispec)) then
            ! dipole
            if (DIPOLE_SOURCE_IN_FLUID) then
              ! length of component vector
              norm = dsqrt( comp_dir_vect_source_E(isource)**2 &
                          + comp_dir_vect_source_N(isource)**2 &
                          + comp_dir_vect_source_Z_UP(isource)**2 )
              comp_x = comp_dir_vect_source_E(isource)/norm
              comp_y = comp_dir_vect_source_N(isource)/norm
              comp_z = comp_dir_vect_source_Z_UP(isource)/norm
              write(*,*) " DIPOLE FLUID ", comp_x, comp_y, comp_z

              call compute_arrays_source_forcesolution_fluid(ispec, sourcearray, hxis,hetas,hgammas,hpxis,hpetas,hpgammas, &
                                                             factor_source,comp_x,comp_y,comp_z,nu_source(:,:,isource))
            else
              ! point force
              ! identical source array components in x,y,z-direction
              ! source array interpolated on all element GLL points
              ! we use an tilted force defined by its magnitude and the projections
              ! of an arbitrary (non-unitary) direction vector on the E/N/Z_UP basis
              comp_x = 1.0d0
              comp_y = 1.0d0
              comp_z = 1.0d0

              call compute_arrays_source_forcesolution(sourcearray,hxis,hetas,hgammas,factor_source, &
                                                       comp_x,comp_y,comp_z,nu_source(:,:,isource))
            endif
          endif ! acoustic

        else
          ! use of CMTSOLUTION files

          ! elastic or poroelastic moment tensor source
          if (ispec_is_elastic(ispec) .or. ispec_is_poroelastic(ispec)) then
            call compute_arrays_source_cmt(ispec,sourcearray, &
                                           hxis,hetas,hgammas,hpxis,hpetas,hpgammas, &
                                           Mxx(isource),Myy(isource),Mzz(isource), &
                                           Mxy(isource),Mxz(isource),Myz(isource))
          endif

          ! acoustic case
          if (ispec_is_acoustic(ispec)) then
            ! scalar moment of moment tensor values read in from CMTSOLUTION
            ! note: M0 by Dahlen and Tromp, eq. 5.91
            factor_source = 1.d0/sqrt(2.d0) * sqrt( Mxx(isource)**2 + Myy(isource)**2 + Mzz(isource)**2 &
                                                  + 2.d0 * ( Myz(isource)**2 + Mxz(isource)**2 + Mxy(isource)**2) )

            ! scales source such that it would be equivalent to explosion source moment tensor,
            ! where Mxx=Myy=Mzz, others Mxy,.. = zero, in equivalent elastic media
            ! (and getting rid of 1/sqrt(2) factor from scalar moment tensor definition above)
            factor_source = factor_source * sqrt(2.d0) / sqrt(3.d0)

            ! source encoding
            ! determines factor +/-1 depending on sign of moment tensor
            ! (see e.g. Krebs et al., 2009. Fast full-wavefield seismic inversion using encoded sources,
            !   Geophysics, 74 (6), WCC177-WCC188.)
            if (USE_SOURCE_ENCODING) then
              pm1_source_encoding(isource) = sign(1.0d0,Mxx(isource))
              does_source_encoding = .true.
            endif

            comp_x = 1.0d0
            comp_y = 1.0d0
            comp_z = 1.0d0

            ! source array interpolated on all element GLL points
            call compute_arrays_source_forcesolution(sourcearray,hxis,hetas,hgammas,factor_source, &
                                                     comp_x,comp_y,comp_z,nu_source(:,:,isource))
          endif

        endif

        ! stores source excitations
        sourcearrays(isource,:,:,:,:) = sourcearray(:,:,:,:)

      endif
    enddo
  else
    ! SIMULATION_TYPE == 2
    ! allocate dummy array (needed for subroutine calls)
    allocate(sourcearrays(1,1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2072')
    if (ier /= 0) stop 'error allocating dummy sourcearrays'
  endif

  ! user info
  if (USE_SOURCE_ENCODING) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'using source encoding:'
      if (does_source_encoding) then
        write(IMAIN,*) '  sources have been encoded'
      else
        write(IMAIN,*) '  source encoding has no effect (only supported for acoustic sources)'
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! frees dynamically allocated memory
  if (allocated(factor_force_source)) &
    deallocate(factor_force_source)
  if (allocated(comp_dir_vect_source_E)) &
    deallocate(comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP)

  end subroutine setup_sources_precompute_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_receivers_precompute_intp()

  use specfem_par
  use specfem_par_elastic, only: ispec_is_elastic

  implicit none

  ! local parameters
  integer :: irec,irec_local,isource,ier
  integer,dimension(0:NPROC-1) :: tmp_rec_local_all
  integer :: maxrec,maxproc(1)
  double precision :: sizeval,size_s
  ! interpolants
  double precision,dimension(NGLLX) :: hxir
  double precision,dimension(NGLLY) :: hetar
  double precision,dimension(NGLLZ) :: hgammar
  double precision :: hpxir(NGLLX), hpetar(NGLLY),hpgammar(NGLLZ)

  integer :: i,j,k,ispec
  double precision :: hxi,heta,hgamma
  logical :: has_receiver_in_elastic_domain

  ! note: for adjoint simulations (SIMULATION_TYPE == 2),
  !         nrec_local     - is set to the number of sources (CMTSOLUTIONs), which act as "receiver" locations
  !                          for storing receiver (hxir/..) arrays, seismograms or strains
  !
  !         nadj_rec_local - determines the number of adjoint sources, i.e., number of station locations (STATIONS_ADJOINT), which
  !                          act as sources to drive the adjoint wavefield

  ! check if we need to save seismos
  if (SIMULATION_TYPE == 3 .and. (.not. SAVE_SEISMOGRAMS_IN_ADJOINT_RUN)) then
    do_save_seismograms = .false.
  else
    do_save_seismograms = .true.
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'seismograms:'
    if (do_save_seismograms) then
      if (WRITE_SEISMOGRAMS_BY_MAIN) then
        write(IMAIN,*) '  seismograms written by main process only'
      else
        write(IMAIN,*) '  seismograms written by all processes'
      endif
      write(IMAIN,*)
      write(IMAIN,*) '  Total number of simulation steps (NSTEP)                       = ',NSTEP
      write(IMAIN,*) '  writing out seismograms at every NTSTEP_BETWEEN_OUTPUT_SEISMOS = ',NTSTEP_BETWEEN_OUTPUT_SEISMOS
      write(IMAIN,*) '  number of subsampling steps for seismograms                    = ',NTSTEP_BETWEEN_OUTPUT_SAMPLE
      write(IMAIN,*) '  Total number of samples for seismograms                        = ',NSTEP/NTSTEP_BETWEEN_OUTPUT_SAMPLE
      write(IMAIN,*)
    else
      write(IMAIN,*) '  seismograms will not be saved'
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! resets flags if no seismograms are needed
  if (.not. do_save_seismograms) then
    SAVE_SEISMOGRAMS_PRESSURE = .false.
    SAVE_SEISMOGRAMS_DISPLACEMENT = .false.
    SAVE_SEISMOGRAMS_VELOCITY = .false.
    SAVE_SEISMOGRAMS_ACCELERATION = .false.
    SAVE_SEISMOGRAMS_STRAIN = .false.
  endif

  ! checks SU_FORMAT output length
  if (SU_FORMAT .and. (NSTEP/NTSTEP_BETWEEN_OUTPUT_SAMPLE > 32768)) then
    print *
    print *,"!!! BEWARE !!! Two many samples for SU format ! The .su file created will be unusable"
    print *
    call exit_MPI(myrank,'Error: Two many samples for SU format ! The .su file created will be unusable')
  endif

  ! safety check for strain seismograms
  if (SAVE_SEISMOGRAMS_STRAIN .and. WRITE_SEISMOGRAMS_BY_MAIN) then
    call exit_MPI(myrank,'SAVE_SEISMOGRAMS_STRAIN works only correctly when WRITE_SEISMOGRAMS_BY_MAIN = .false.')
  endif

  ! seismogram array length (to write out time portions of the full seismograms)
  nlength_seismogram = NTSTEP_BETWEEN_OUTPUT_SEISMOS / NTSTEP_BETWEEN_OUTPUT_SAMPLE

  ! statistics about allocation memory for seismograms & source_adjoint seismograms
  ! gather from secondarys on main
  call gather_all_singlei(nrec_local,tmp_rec_local_all,NPROC)

  ! user output
  if (myrank == 0) then
    ! determines maximum number of local receivers and corresponding rank
    maxrec = maxval(tmp_rec_local_all(:))
    ! note: MAXLOC will determine the lower bound index as '1'.
    maxproc = maxloc(tmp_rec_local_all(:)) - 1

    ! seismograms array size in MB
    ! seismograms need seismograms(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS/NTSTEP_BETWEEN_OUTPUT_SAMPLE)
    size_s = dble(maxrec) * dble(NDIM) * dble(nlength_seismogram * CUSTOM_REAL / 1024. / 1024. )
    sizeval = 0.d0
    if (SAVE_SEISMOGRAMS_DISPLACEMENT) sizeval = sizeval + size_s
    if (SAVE_SEISMOGRAMS_VELOCITY) sizeval = sizeval + size_s
    if (SAVE_SEISMOGRAMS_ACCELERATION) sizeval = sizeval + size_s
    if (SAVE_SEISMOGRAMS_PRESSURE) sizeval = sizeval + dble(NB_RUNS_ACOUSTIC_GPU) * size_s
    ! strain seismogram needs seismograms_eps(NDIM,NDIM,nrec_local,NSTEP/NTSTEP_BETWEEN_OUTPUT_SAMPLE)
    if (SAVE_SEISMOGRAMS_STRAIN) sizeval = sizeval + dble(NDIM) * size_s

    ! outputs info
    write(IMAIN,*) '  maximum number of local receivers is ',maxrec,' in slice ',maxproc(1)
    if (do_save_seismograms) then
      write(IMAIN,*) '  size of maximum seismogram array       = ', sngl(sizeval),'MB'
      write(IMAIN,*) '                                         = ', sngl(sizeval/1024.d0),'GB'
    endif
    write(IMAIN,*)
    call flush_IMAIN()

    ! checks if receivers are available
    if (maxrec == 0 .and. .not. INVERSE_FWI_FULL_PROBLEM) then
      print *,'Error: no receivers in slices found'
      print *,'Please check if stations in DATA/STATIONS are placed correctly.'
      call exit_MPI(myrank,'Error no receivers in slices found')
    endif
  endif

  ! adjoint sources
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! note: nadj_rec_local is the number of "adjoint sources".
    !       that is, local receiver locations which will act as adjoint source locations.

    ! gather from secondarys on main
    call gather_all_singlei(nadj_rec_local,tmp_rec_local_all,NPROC)

    ! user output
    if (myrank == 0) then
      ! determines maximum number of local receivers and corresponding rank
      maxrec = maxval(tmp_rec_local_all(:))
      ! note: MAXLOC will determine the lower bound index as '1'.
      maxproc = maxloc(tmp_rec_local_all(:)) - 1

      ! source_adjoint size in MB
      !source_adjoint(NDIM,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC)
      sizeval = dble(maxrec) * dble(NDIM * NTSTEP_BETWEEN_READ_ADJSRC * CUSTOM_REAL / 1024. / 1024. )

      ! outputs info
      write(IMAIN,*) 'adjoint source arrays:'
      write(IMAIN,*) '  reading adjoint sources at every NTSTEP_BETWEEN_READ_ADJSRC = ',NTSTEP_BETWEEN_READ_ADJSRC
      write(IMAIN,*) '  maximum number of local adjoint sources is ',maxrec,' in slice ',maxproc(1)
      write(IMAIN,*) '  size of maximum adjoint source array = ', sngl(sizeval),'MB'
      write(IMAIN,*) '                                       = ', sngl(sizeval/1024.d0),'GB'
      write(IMAIN,*)
      call flush_IMAIN()

      ! check if adjoint stations are available
      if (maxrec == 0 .and. .not. INVERSE_FWI_FULL_PROBLEM) then
        print *,'Error: no adjoint sources in slices found'
        print *,'Please check if adjoint stations in DATA/STATIONS_ADJOINT are placed correctly.'
        call exit_MPI(myrank,'Error no adjoint sources in slices found')
      endif
    endif
  endif

  ! receivers: for forward/kernel simulations, receivers are determined by the STATIONS files,
  !            for pure adjoint simulations, CMT source locations become "receiver" locations.
  !
  ! following arrays are used to handle stations at "receiver" locations,
  ! needed for forward and kernel simulations
  !
  ! needs to be allocate for subroutine calls (even if nrec_local == 0)
  allocate(number_receiver_global(nrec_local),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2074')
  if (ier /= 0) stop 'error allocating array number_receiver_global'
  number_receiver_global(:) = 0

  ! stores local receivers interpolation factors
  if (nrec_local > 0) then
    ! allocate Lagrange interpolators for receivers
    allocate(hxir_store(NGLLX,nrec_local),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2075')
    allocate(hetar_store(NGLLY,nrec_local),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2076')
    allocate(hgammar_store(NGLLZ,nrec_local),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2077')
    if (ier /= 0) stop 'error allocating array hxir_store etc.'
    hxir_store(:,:) = 0.0; hetar_store(:,:) = 0.0; hgammar_store(:,:) = 0.0

    ! allocates derivatives
    if (SIMULATION_TYPE == 2) then
      allocate(hpxir_store(NGLLX,nrec_local),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2078')
      allocate(hpetar_store(NGLLY,nrec_local),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2079')
      allocate(hpgammar_store(NGLLZ,nrec_local),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2080')
      if (ier /= 0) stop 'error allocating array hpxir_store'
      hpxir_store(:,:) = 0.0; hpetar_store(:,:) = 0.0; hpgammar_store(:,:) = 0.0
    endif

    ! define local to global receiver numbering mapping
    number_receiver_global(:) = 0
    irec_local = 0
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      do irec = 1,nrec
        if (myrank == islice_selected_rec(irec)) then
          irec_local = irec_local + 1
          number_receiver_global(irec_local) = irec
        endif
      enddo
    else
      ! note: pure adjoint runs (SIMULATION_TYPE == 2):
      !       uses CMT source locations to store "receiver" seismograms
      do isource = 1,NSOURCES
        if (myrank == islice_selected_source(isource)) then
          irec_local = irec_local + 1
          number_receiver_global(irec_local) = isource
        endif
      enddo
    endif

    ! define and store Lagrange interpolators at all the receivers
    do irec_local = 1,nrec_local
      irec = number_receiver_global(irec_local)

      if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
        ! receiver positions
        call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
        call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
        call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
      else
        ! source positions
        ! note: pure adjoint runs (SIMULATION_TYPE == 2):
        !       uses CMT source locations to store "receiver" seismograms
        call lagrange_any(xi_source(irec),NGLLX,xigll,hxir,hpxir)
        call lagrange_any(eta_source(irec),NGLLY,yigll,hetar,hpetar)
        call lagrange_any(gamma_source(irec),NGLLZ,zigll,hgammar,hpgammar)
      endif

      ! stores interpolators
      hxir_store(:,irec_local) = hxir(:)
      hetar_store(:,irec_local) = hetar(:)
      hgammar_store(:,irec_local) = hgammar(:)

      ! stores derivatives
      if (SIMULATION_TYPE == 2) then
        hpxir_store(:,irec_local) = hpxir(:)
        hpetar_store(:,irec_local) = hpetar(:)
        hpgammar_store(:,irec_local) = hpgammar(:)
      endif
    enddo
  else
    ! dummy arrays
    ! VM VM need to allocate Lagrange interpolators for receivers with 0 because it is used
    ! in calling subroutines parmeters. (otherwise it can be crash at runtime).
    allocate(hxir_store(1,1), &
             hetar_store(1,1), &
             hgammar_store(1,1))
    ! allocates derivatives
    if (SIMULATION_TYPE == 2) then
      allocate(hpxir_store(1,1), &
               hpetar_store(1,1), &
               hpgammar_store(1,1))
    endif
  endif ! nrec_local > 0

  ! seismograms
  if (nrec_local > 0 .and. do_save_seismograms) then
    ! note: nrec_local is the number of local "receivers" for which we will store seismograms
    !
    ! allocate seismogram array
    if (SAVE_SEISMOGRAMS_DISPLACEMENT) then
      allocate(seismograms_d(NDIM,nrec_local,nlength_seismogram),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2087')
    else
      allocate(seismograms_d(1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2088')
    endif
    if (ier /= 0) stop 'error allocating array seismograms_d'

    if (SAVE_SEISMOGRAMS_VELOCITY) then
      allocate(seismograms_v(NDIM,nrec_local,nlength_seismogram),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2089')
    else
      allocate(seismograms_v(1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2090')
    endif
    if (ier /= 0) stop 'error allocating array seismograms_v'

    if (SAVE_SEISMOGRAMS_ACCELERATION) then
      allocate(seismograms_a(NDIM,nrec_local,nlength_seismogram),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2091')
    else
      allocate(seismograms_a(1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2092')
    endif
    if (ier /= 0) stop 'error allocating array seismograms_a'

    if (SAVE_SEISMOGRAMS_PRESSURE) then
      !NB_RUNS_ACOUSTIC_GPU is set to 1 by default in constants.h
      !daniel todo: use single component only, i.e., NDIM == 1, for pressure is enough
      allocate(seismograms_p(NDIM,nrec_local*NB_RUNS_ACOUSTIC_GPU,nlength_seismogram),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2093')
    else
      allocate(seismograms_p(1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2094')
    endif
    if (ier /= 0) stop 'error allocating array seismograms_p'

    if (SAVE_SEISMOGRAMS_STRAIN) then
      allocate(seismograms_eps(NDIM,NDIM,nrec_local,nlength_seismogram),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2087')
    else
      allocate(seismograms_eps(1,1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2088')
    endif
    if (ier /= 0) stop 'error allocating array seismograms_d'

    ! initialize seismograms
    seismograms_d(:,:,:) = 0._CUSTOM_REAL
    seismograms_v(:,:,:) = 0._CUSTOM_REAL
    seismograms_a(:,:,:) = 0._CUSTOM_REAL
    seismograms_p(:,:,:) = 0._CUSTOM_REAL
    seismograms_eps(:,:,:,:) = 0._CUSTOM_REAL
  else
    ! dummy allocations
    allocate(seismograms_d(1,1,1), &
             seismograms_v(1,1,1), &
             seismograms_a(1,1,1), &
             seismograms_p(1,1,1), &
             seismograms_eps(1,1,1,1))
  endif

  ! adjoint sources
  ! optimizing arrays for adjoint sources
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! local adjoint sources arrays
    if (nadj_rec_local > 0) then
      ! determines adjoint sources arrays
      if (SIMULATION_TYPE == 2) then
        ! pure adjoint simulations
        allocate(number_adjsources_global(nadj_rec_local),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating number_adjsources_global array')

        ! addressing from local to global receiver index
        number_adjsources_global(:) = 0
        irec_local = 0
        do irec = 1,nrec
          ! add the source (only if this proc carries the source)
          if (myrank == islice_selected_rec(irec)) then
            irec_local = irec_local + 1
            number_adjsources_global(irec_local) = irec
          endif
        enddo
        if (irec_local /= nadj_rec_local) stop 'Error invalid number of nadj_rec_local found'

        ! allocate Lagrange interpolators for adjoint sources
        !
        ! note: adjoint sources for SIMULATION_TYPE == 2 and 3 are located at the receivers,
        !       however, the interpolator arrays hxir_store are used for "receiver" locations which are different
        !       for pure adjoint or kernel simulations
        !
        !       we will thus allocate interpolator arrays especially for adjoint source locations. for kernel simulations,
        !       these would be the same as hxir_store, but not for pure adjoint simulations.
        !
        ! storing these arrays is cheaper than storing a full (i,j,k) array for each element
        allocate(hxir_adjstore(NGLLX,nadj_rec_local),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating hxir_adjstore array')
        allocate(hetar_adjstore(NGLLY,nadj_rec_local),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating hetar_adjstore array')
        allocate(hgammar_adjstore(NGLLZ,nadj_rec_local),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating hgammar_adjstore array')

        ! define and store Lagrange interpolators at all the adjoint source locations
        do irec_local = 1,nadj_rec_local
          irec = number_adjsources_global(irec_local)

          ! receiver positions (become adjoint source locations)
          call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
          call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
          call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)

          ! stores interpolators
          hxir_adjstore(:,irec_local) = hxir(:)
          hetar_adjstore(:,irec_local) = hetar(:)
          hgammar_adjstore(:,irec_local) = hgammar(:)
        enddo
      else
        ! kernel simulations (SIMULATION_TYPE == 3)
        ! adjoint source arrays and receiver arrays are the same, no need to allocate new arrays, just point to the existing ones
        number_adjsources_global => number_receiver_global
        hxir_adjstore => hxir_store
        hetar_adjstore => hetar_store
        hgammar_adjstore => hgammar_store
      endif
    else
      ! dummy arrays
      if (SIMULATION_TYPE == 2) then
        allocate(number_adjsources_global(1), &
                 hxir_adjstore(1,1), &
                 hetar_adjstore(1,1), &
                 hgammar_adjstore(1,1))
      else
        number_adjsources_global => number_receiver_global
        hxir_adjstore => hxir_store
        hetar_adjstore => hetar_store
        hgammar_adjstore => hgammar_store
      endif
    endif ! nadj_rec_local
  endif

  ! safety check
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! adjoint source in this partitions
    if (nadj_rec_local > 0) then
      ! added this because Lei Zhang may have detected a problem
      do irec_local = 1, nadj_rec_local
        irec = number_adjsources_global(irec_local)
        if (irec <= 0) stop 'Error invalid irec for local adjoint source'
        ! adds source array
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              hxi = hxir_adjstore(i,irec_local)
              heta = hetar_adjstore(j,irec_local)
              hgamma = hgammar_adjstore(k,irec_local)
              ! debug
              !print *,'hxi/heta/hgamma = ',hxi,heta,hgamma,irec_local,i,j,k
              ! checks if array values valid
              ! Lagrange interpolators shoud be between [-0.2,1.0]
              if (abs(hxi) > 1.1d0 .or. abs(heta) > 1.1d0 .or. abs(hgamma) > 1.1d0) then
                print *,'ERROR: trying to use arrays hxir_adjstore/hetar_adjstore/hgammar_adjstore with irec_local = ', &
                        irec_local,' but these arrays are invalid!'
                print *,'  invalid hxi/heta/hgamma = ',hxi,heta,hgamma
                call exit_MPI_without_rank('ERROR: trying to use arrays hxir_adjstore/hetar_adjstore/hgammar_adjstore &
                                           &but these arrays are invalid!')
              endif
            enddo
          enddo
        enddo
      enddo
    endif
  endif

  ! seismogram check
  if (SAVE_SEISMOGRAMS_PRESSURE .and. ELASTIC_SIMULATION .and. ATTENUATION) then
    ! pressure output for receivers in elastic domains is only calculated based on non-attenuation stresses
    has_receiver_in_elastic_domain = .false.
    do irec_local = 1,nrec_local
      ! gets global number of that receiver
      irec = number_receiver_global(irec_local)
      ! spectral element in which the receiver is located
      if (SIMULATION_TYPE == 2) then
        ! adjoint "receivers" are located at CMT source positions
        ! note: we take here xi_source,.. when FASTER_RECEIVERS_POINTS_ONLY is set
        ispec = ispec_selected_source(irec)
      else
        ! receiver located at station positions
        ispec = ispec_selected_rec(irec)
      endif
      if (ispec_is_elastic(ispec)) then
        has_receiver_in_elastic_domain = .true.
        exit
      endif
    enddo
    ! warning
    if (has_receiver_in_elastic_domain) then
      write(IMAIN,*) '  ***'
      write(IMAIN,*) '  *** Warning: Pressure seismogram output for receivers in elastic domains&
                                     &is valid only for non-attenuation case ***'
      write(IMAIN,*) '  ***'
      call flush_IMAIN()
    endif
  endif

  end subroutine setup_receivers_precompute_intp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_receivers_VTKfile()

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,OUTPUT_FILES,IOUT_VTK,myrank

  use specfem_par, only: SIMULATION_TYPE,NSOURCES,nrec,NGNOD,NSPEC_AB,NGLOB_AB, &
    xstore,ystore,zstore,ibool, &
    ispec_selected_source,islice_selected_source,xi_source,eta_source,gamma_source, &
    ispec_selected_rec,islice_selected_rec,xi_receiver,eta_receiver,gamma_receiver

  implicit none

  ! local parameters
  double precision :: shape3D(NGNOD)
  double precision :: xil,etal,gammal
  double precision :: xmesh,ymesh,zmesh
  real(kind=CUSTOM_REAL),dimension(NGNOD) :: xelm,yelm,zelm
  integer :: ia,ispec,isource,irec,ier,totalpoints

  character(len=MAX_STRING_LEN) :: filename,filename_new,system_command,system_command1,system_command2

  ! determines number of points for vtk file
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    totalpoints = NSOURCES + nrec
  else
    ! pure adjoint simulation only needs receivers
    totalpoints = nrec
  endif

  if (myrank == 0) then
    ! vtk file
    open(IOUT_VTK,file=trim(OUTPUT_FILES)//'/sr.vtk',status='unknown',iostat=ier)
    if (ier /= 0) stop 'error opening sr.vtk file'
    ! vtk header
    write(IOUT_VTK,'(a)') '# vtk DataFile Version 2.0'
    write(IOUT_VTK,'(a)') 'Source and Receiver VTK file'
    write(IOUT_VTK,'(a)') 'ASCII'
    write(IOUT_VTK,'(a)') 'DATASET POLYDATA'
    write(IOUT_VTK, '(a,i6,a)') 'POINTS ', totalpoints, ' float'
  endif

  ! sources
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    do isource = 1,NSOURCES
      ! spectral element id
      ispec = ispec_selected_source(isource)

      ! gets element anchor nodes
      if (myrank == islice_selected_source(isource)) then
        ! find the coordinates of the anchor nodes of the element
        call eval_shape3D_element_anchors(xelm,yelm,zelm,ispec,ibool,xstore,ystore,zstore,NSPEC_AB,NGLOB_AB)
      endif

      ! main collects corner locations
      if (islice_selected_source(isource) /= 0) then
        if (myrank == 0) then
          call recvv_cr(xelm,NGNOD,islice_selected_source(isource),0)
          call recvv_cr(yelm,NGNOD,islice_selected_source(isource),0)
          call recvv_cr(zelm,NGNOD,islice_selected_source(isource),0)
        else if (myrank == islice_selected_source(isource)) then
          call sendv_cr(xelm,NGNOD,0,0)
          call sendv_cr(yelm,NGNOD,0,0)
          call sendv_cr(zelm,NGNOD,0,0)
        endif
      endif

      if (myrank == 0) then
        ! get the 3-D shape functions
        xil = xi_source(isource)
        etal = eta_source(isource)
        gammal = gamma_source(isource)
        call eval_shape3D_single(shape3D,xil,etal,gammal,NGNOD)

        ! interpolates source locations
        xmesh = 0.d0
        ymesh = 0.d0
        zmesh = 0.d0
        do ia = 1,NGNOD
          xmesh = xmesh + shape3D(ia)*xelm(ia)
          ymesh = ymesh + shape3D(ia)*yelm(ia)
          zmesh = zmesh + shape3D(ia)*zelm(ia)
        enddo

        ! writes out to VTK file
        write(IOUT_VTK,'(3e18.6)') xmesh,ymesh,zmesh
      endif
    enddo ! NSOURCES
  endif

  ! receivers
  do irec = 1,nrec
    ispec = ispec_selected_rec(irec)

    ! find the coordinates of the anchor (eight corners for NGNOD=8) nodes of the element
    if (myrank == islice_selected_rec(irec)) then
      call eval_shape3D_element_anchors(xelm,yelm,zelm,ispec,ibool,xstore,ystore,zstore,NSPEC_AB,NGLOB_AB)
    endif

    ! main collects corner locations
    if (islice_selected_rec(irec) /= 0) then
      if (myrank == 0) then
        call recvv_cr(xelm,NGNOD,islice_selected_rec(irec),0)
        call recvv_cr(yelm,NGNOD,islice_selected_rec(irec),0)
        call recvv_cr(zelm,NGNOD,islice_selected_rec(irec),0)
      else if (myrank == islice_selected_rec(irec)) then
        call sendv_cr(xelm,NGNOD,0,0)
        call sendv_cr(yelm,NGNOD,0,0)
        call sendv_cr(zelm,NGNOD,0,0)
      endif
    endif

    if (myrank == 0) then
      ! get the 3-D shape functions
      xil = xi_receiver(irec)
      etal = eta_receiver(irec)
      gammal = gamma_receiver(irec)
      call eval_shape3D_single(shape3D,xil,etal,gammal,NGNOD)

      ! interpolates receiver locations
      xmesh = 0.d0
      ymesh = 0.d0
      zmesh = 0.d0
      do ia=1,NGNOD
        xmesh = xmesh + shape3D(ia)*xelm(ia)
        ymesh = ymesh + shape3D(ia)*yelm(ia)
        zmesh = zmesh + shape3D(ia)*zelm(ia)
      enddo

      ! writes out to VTK file
      write(IOUT_VTK,'(3e18.6)') xmesh,ymesh,zmesh
    endif
  enddo

  ! closes vtk file
  if (myrank == 0) then
    write(IOUT_VTK,*)
    close(IOUT_VTK)

    ! creates additional receiver and source files
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      ! extracts receiver locations
      filename = trim(OUTPUT_FILES)//'/sr.vtk'
      filename_new = trim(OUTPUT_FILES)//'/receiver.vtk'

      ! vtk file for receivers only
      write(system_command, &
  "('awk ',a1,'{if (NR < 5) print $0;if (NR == 5)print ',a1,'POINTS',i6,' float',a1,';if &
      &(NR > 5+',i6,')print $0}',a1,'<',a,'>',a)")&
      "'",'"',nrec,'"',NSOURCES,"'",trim(filename),trim(filename_new)

      ! extracts source locations
      filename_new = trim(OUTPUT_FILES)//'/source.vtk'

      write(system_command1, &
  "('awk ',a1,'{if (NR < 5) print $0;if (NR == 5)print ',a1,'POINTS',i6,' float',a1,';')") &
        "'",'"',NSOURCES,'"'

      write(system_command2, &
  "('if (NR > 5 && NR < 6+',i6,')print $0}END{print ',a,'}',a1,'<',a,'>',a)") &
        NSOURCES,'" "',"'",trim(filename),trim(filename_new)

      system_command = trim(system_command1)//trim(system_command2)

    endif
  endif

  end subroutine setup_sources_receivers_VTKfile

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_search_kdtree()

  use specfem_par

  use kdtree_search, only: kdtree_setup,kdtree_set_verbose, &
    kdtree_num_nodes,kdtree_nodes_location,kdtree_nodes_index

  implicit none
  integer :: ier,ispec,iglob
  ! determines tree points
  logical :: use_midpoints_only
  integer :: inodes,i,j,k

  ! checks if anything to do
  if (DO_BRUTE_FORCE_POINT_SEARCH) return

  ! kdtree search

  ! determines tree size
  if (NSPEC_AB > 1000000) then
    ! high-resolution mesh
    ! only midpoints for search, should be sufficient to get accurate location
    use_midpoints_only = .true.
  else
    ! low-resolution mesh
    ! uses element's inner points
    use_midpoints_only = .false.
  endif

  ! sets total number of tree points
  if (use_midpoints_only) then
    ! small tree size
    kdtree_num_nodes = NSPEC_AB
  else
    ! uses all internal GLL points for search tree
    ! internal GLL points ( 2 to NGLLX-1 )
    kdtree_num_nodes = NSPEC_AB * (NGLLX-2)*(NGLLY-2)*(NGLLZ-2)
  endif

  ! allocates tree arrays
  allocate(kdtree_nodes_location(NDIM,kdtree_num_nodes),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2096')
  if (ier /= 0) stop 'Error allocating kdtree_nodes_location arrays'

  allocate(kdtree_nodes_index(kdtree_num_nodes),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2097')
  if (ier /= 0) stop 'Error allocating kdtree_nodes_index arrays'

  ! tree verbosity
  if (myrank == 0) call kdtree_set_verbose()

  ! prepares search arrays, each element takes its internal GLL points for tree search
  kdtree_nodes_index(:) = 0
  kdtree_nodes_location(:,:) = 0.0

  ! adds tree nodes
  inodes = 0
  if (use_midpoints_only) then
    ! fills kd-tree arrays
    do ispec = 1,NSPEC_AB
      ! counts nodes
      inodes = inodes + 1
      if (inodes > kdtree_num_nodes) stop 'Error index inodes bigger than kdtree_num_nodes'

      ! adds node index ( index points to same ispec for all internal GLL points)
      kdtree_nodes_index(ispec) = ispec

      ! adds node location (of midpoint)
      iglob = ibool(MIDX,MIDY,MIDZ,ispec)
      kdtree_nodes_location(1,ispec) = xstore(iglob)
      kdtree_nodes_location(2,ispec) = ystore(iglob)
      kdtree_nodes_location(3,ispec) = zstore(iglob)
    enddo
  else
    ! all internal GLL points
    do ispec = 1,NSPEC_AB
      do k = 2,NGLLZ-1
        do j = 2,NGLLY-1
          do i = 2,NGLLX-1
            iglob = ibool(i,j,k,ispec)

            ! counts nodes
            inodes = inodes + 1
            if (inodes > kdtree_num_nodes) stop 'Error index inodes bigger than kdtree_num_nodes'

            ! adds node index (index points to same ispec for all internal GLL points)
            kdtree_nodes_index(inodes) = ispec

            ! adds node location
            kdtree_nodes_location(1,inodes) = xstore(iglob)
            kdtree_nodes_location(2,inodes) = ystore(iglob)
            kdtree_nodes_location(3,inodes) = zstore(iglob)
          enddo
        enddo
      enddo
    enddo
  endif
  if (inodes /= kdtree_num_nodes) stop 'Error index inodes does not match kdtree_num_nodes'

  ! creates kd-tree for searching
  ! serial way
  !do i = 0,NPROCTOT_VAL-1
  !  if (myrank == i) then
  !    print *,'kd-tree setup for process: ',myrank
  !    call kdtree_setup()
  !  endif
  !  call synchronize_all()
  !enddo
  ! parallel way
  call kdtree_setup()

  ! synchronizes all mpi-processes
  call synchronize_all()

  end subroutine setup_search_kdtree

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_free_kdtree()

  use specfem_par
  use kdtree_search, only: kdtree_delete,kdtree_nodes_location,kdtree_nodes_index

  implicit none

  ! checks if anything to do
  if (DO_BRUTE_FORCE_POINT_SEARCH) return
  if (INVERSE_FWI_FULL_PROBLEM) return

  ! deletes tree arrays
  deallocate(kdtree_nodes_location)
  deallocate(kdtree_nodes_index)

  ! deletes search tree nodes
  call kdtree_delete()

  end subroutine setup_free_kdtree
