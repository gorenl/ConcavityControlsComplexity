program DivideAndCapture

  ! Beta version of new cascade algorithm that incorporates river capture
  ! but not the computation of precise divide locations
  ! Jean BRAUN 20 july 2009

  use definitions

  implicit none

  ! for a definition of the derived types see the file module_definitions.f90

  type (geom) geometry
  type (del) delaunay
  !type (del) rddelaunay
  type (netw) network
  type (stck) stack
  type (timr) timer
  type (parm) params
  character*11  title
  integer old_ncapture, i, j, k, m, l, n
  double precision li,lj
  double precision ymin,ymax
  integer node1,node2,node3,cat1,cat2,cat3
  ! sets a common for timing of various subroutines
  ! the results are printed to the screen at the end of the run

  common /global_timer/ timer

  call start_timing
  call initialize_parameters (params) ! initializes the model parameters
  geometry%surface_share=0.d0
  print*, 'going to reinit_geometry',params%read_restart
  if (params%read_restart) then 
     call reinit_geometry(geometry,params,network,delaunay)
  else
     call initialize_geometry (geometry,params) ! initializes the model geometry and history
  endif

  if (.not.params%read_restart) then
     delaunay%ntriangles=0
     call calculate_delaunay (geometry,delaunay) ! calculates Delaunay and neighbour (triangle) information
  endif

  if (params%read_restart) then
     call find_network (geometry,delaunay,network,1) ! find neighbour list
  else
     network%nnode=0
     call find_network (geometry,delaunay,network,0) ! find neighbour list and donor/receiver list
  endif

  if (.not.params%read_restart) then
     call find_precipitation (geometry,params) ! orography (not fully working yet)
  endif

  !if (.not.params%read_restart) then
  stack%nnode=0
  call find_order (geometry,network,stack,delaunay,params) ! find stack order to perform the erosion and cascade algorithms
  !endif
  call find_strahler(geometry,network,stack)
  call VTK (geometry,delaunay,params,network,stack)
  if (params%ascii) call write_ascii(geometry,params,network)
  call output_geometry(geometry,params,network,delaunay)


  !more frequent global data on the system is output to these files
  open (80,file='ASCII/capture_data',status='unknown')
  write(80,'(e15.7,i10)') params%time, geometry%ncapture
  old_ncapture = geometry%ncapture
  do while (params%time.lt.params%tfinal) ! start of time loop


     !print*, '1'
     params%deltat=min(params%deltat,params%tfinal-params%time)
     params%time=params%time+params%deltat
     params%istep=params%istep+1


     call uplift_and_advect (geometry,params) ! move and uplift the grid according to the prescribed velocity field
     !print*, '2'
     if (params%move_points) then !with horizontal advection
        call calculate_delaunay (geometry,delaunay)     !redo delaunay triangulation
        !print*, '3'
        call find_network (geometry,delaunay,network,0)  ! find neighbour list and donor/receiver list
        !print*, '4'
     endif

     if (params%add_nodes) then 
        call add_remove_nodes(geometry,network,params,delaunay)
        !print*, '5'
        call calculate_delaunay (geometry,delaunay)     !redo delaunay triangulation
        !print*, '6'
        call find_network (geometry,delaunay,network,0) ! find neighbour list and donor/receiver list
        !print*, '7'      
     endif

     if (params%capture.or.params%divide.or.params%small_divide) then
        if (params%transient_divide) call update_erosion_rate_history(geometry,params)
        call captures_and_divides (geometry,network,params,delaunay) ! find potential captures and adjusts network accordingly  
     endif
     !print*, '8'
     call find_precipitation (geometry,params) ! orography (not fully working yet)
     call find_surface (geometry,network,params) ! find surface area attached to each point
     !call find_polygon_surface (geometry,network,params,delaunay) ! find surface area attached to each point
     if (params%add_nodes) deallocate(stack%order)
     call find_order (geometry,network,stack,delaunay,params) ! find stack order to perform the erosion and cascade algorithms
     call find_catchment (geometry,network,stack,delaunay,params) ! find catchments
     call find_discharge (geometry,network,stack) ! compute discharge (improved cascade algorithm)
     call erode (geometry,network,stack,delaunay,params) ! erode


     if ((params%istep/params%freq)*params%freq.eq.params%istep) then
        !  call calculate_river_divide_delaunay (geometry,rddelaunay) !calc delaunay of grid and divides 	
        call find_strahler(geometry,network,stack) 
        call VTK (geometry,delaunay,params,network,stack)
        if (params%ascii) call write_ascii(geometry,params,network)
        call output_geometry(geometry,params,network,delaunay)
        print*, 'max number of neighbors is:', maxval(geometry%nb)
     endif

     !frequent data is output here
     if (old_ncapture.ne.geometry%ncapture.or.MOD(params%istep,100).eq.0) then
        write(80,'(e15.7,i10)') params%time, geometry%ncapture
        old_ncapture = geometry%ncapture
     endif

  enddo ! end of time loop
  
  ymin=minval(geometry%y)
  ymax=maxval(geometry%y)
  open (81,file='ASCII/main_divide',status='unknown')
  do i=1,delaunay%ntriangles
     node1=delaunay%icon(1,i)
     node2=delaunay%icon(2,i)
     node3=delaunay%icon(3,i)
     cat1=geometry%catchment(node1)
     cat2=geometry%catchment(node2)
     cat3=geometry%catchment(node3)
     if (network%ndon(node1).eq.0.and.network%ndon(node2).eq.0) then
        if (geometry%y(cat1).eq.ymin.and.geometry%y(cat2).eq.ymax &
             .or.geometry%y(cat1).eq.ymax.and.geometry%y(cat2).eq.ymin) then
           write(81,*) geometry%x(node1), geometry%y(node1)
           write(81,*) geometry%x(node2), geometry%y(node2)
        endif
     endif
     if (network%ndon(node3).eq.0.and.network%ndon(node2).eq.0) then
        if (geometry%y(cat2).eq.ymin.and.geometry%y(cat3).eq.ymax &
             .or.geometry%y(cat2).eq.ymax.and.geometry%y(cat3).eq.ymin) then
           write(81,*) geometry%x(node3), geometry%y(node3)
           write(81,*) geometry%x(node2), geometry%y(node2)
        endif
     endif
     if (network%ndon(node1).eq.0.and.network%ndon(node3).eq.0) then
        if (geometry%y(cat1).eq.ymin.and.geometry%y(cat3).eq.ymax &
             .or.geometry%y(cat1).eq.ymax.and.geometry%y(cat3).eq.ymin) then
           write(81,*) geometry%x(node1), geometry%y(node1)
           write(81,*) geometry%x(node3), geometry%y(node3)
        endif
     endif
  enddo
  close (81)

  write(80,'(e15.7,i10)') params%time, geometry%ncapture
  close (80)

  !call show (geometry,delaunay,network,stack,params)
  call find_strahler(geometry,network,stack)
  call VTK (geometry,delaunay,params,network,stack)
  if (params%ascii) call write_ascii(geometry,params,network)
  call output_z_tau(params,geometry,network,stack)
  call output_geometry(geometry,params,network,delaunay)
  call find_hack(geometry,network,stack)
  print*,params%istep,'Steps'
  call show_timing

end program DivideAndCapture
