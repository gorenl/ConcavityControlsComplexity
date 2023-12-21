subroutine initialize_parameters (params)

! This routine is used to initialize the model parameters
! i.e. those that do not have a spatial variation (I;E. not attached to a node)
! All of these parameters are stored in a (parm) derived type params

use definitions

implicit none

type (parm) params

double precision u,p,slope,k,sloped,area, sloped_total
integer nsteps, nout

call time_in ('initialize_parameters')
   
!params%deltat=.1d3                               ! time step length (in years)
params%deltat=.1d4 
params%time=0.d0                                 ! time at start (should always be 0)
!params%tfinal=.1d4                               ! final time (in years)
params%tfinal= 10.d7                               ! final time (in years)
!params%tfinal= 700
params%freq=2001                           ! frequency of plots in number of time steps 
params%istep=0                                   ! time step counter (should always be set to 0)
params%n=1.d0                                       ! slope exponent in fluvial incision law
params%m=0.1d0                                    ! area (or discharge) exponent in fluvial incision law
params%h=2.d0                                    ! Hack's exponent
params%hmn=1.d0-params%h*params%m/params%n       ! exponent combination parameter (do not change)
params%ka=2.d0/3.d0                              ! Hack's constant (no unit)
!params%ka=1082.d0
params%xc=1600.d0                                ! hillslope length (in the vicinity of divides) (in m)
params%tanthetac=.3839d0                            ! tangent of hillslope slope
params%rainfall_available=30.d0                  ! first orographic parameter (not fully implemented)
params%rainfall_minimum=.5d0                     ! second orographic parameter (not fully implemented)
params%rainfall_height=1.d0                      ! this is used for constant rainfall parameter
params%diffusivity=.3d0                           ! diffusivity for hilltops
params%diffusivity=params%diffusivity*2.d0          !  double diffusivity as it always appears with a 2*
params%min_erosion_rate=.1d-4                     ! minimum allowable erosion rate for diffusion channel head
!params%min_tan_head_slope = 3.1d-2                   ! minimum slope of channel head (used as a threshold to diffusion)
params%min_tan_head_slope = 1.125d-2
params%lmax=2.25d3                                 ! maximum distance between two points before adding a node
params%amin=.172d7                                ! minimum catchment head area before node remove_d
!params%amin=.1d5
params%ldivmax=3.d3                              ! maximum divide length for adding nodes
params%uplift_scalar1=5.0d-4                      ! uplift rate as constant in interior
params%uplift_scalar2=5.0d-4                      ! uplift rate on boundary (used only for channel head calc)
params%max_adv = 1.0d-3                          ! max advection velocity 
params%k_scalar1=9.2828d-4                            ! erodibility constant

params%plot_triangles=.FALSE.                      ! flag to plot triangles
params%plot_receiver=.TRUE.                        ! flag to plot node to receiver connections
params%plot_donors=.FALSE.                         ! flag to plot node to donors connections
params%plot_no_receiver=.TRUE.                     ! flag to plot nodes with no receiver
params%plot_no_donors=.FALSE.                      ! flag to plot nodes with no donor
params%write_discharge=.FALSE.                     ! flag to write discharge
params%write_stack_order=.FALSE.                   ! flag to write stack order
params%plot_catchment=.FALSE.                      ! flag to plot catchments
params%plot_height=.TRUE.                          ! flag to plot height contours
params%plot_precipitation=.FALSE.                  ! flag to plot precipitation contours

params%move_points=.FALSE.                         ! flag to allow the horizontal motion of nodes
params%capture=.FALSE.                              ! flag to allow for captures
params%divide=.FALSE.                               ! flag to allow for divide calculations
params%small_divide=.FALSE.                         ! flag to allow for small divide calculations (case where divide is very close to zj)
params%diffusion=.FALSE.                            ! flag to allow diffusive hilltops
params%transient_divide=.TRUE.                     ! flag to calculate transient elevation of fluvial part of divides 
params%num_bins=11
params%sample_per_bin=1000
params%add_nodes=.TRUE.                            ! flag to allow adding nodes

params%read_restart =.FALSE.                          !flag to read GeoFile to restart run from the end of previosu run
params%num_restart = 0

params%ascii = .TRUE.

!
!  echo parameters to screen
!

print*, ''
print*, 'DAC in CASCADE mode'
print*, 'hard code change precipitation'
print*, 'Timestepping parameters:'
print*, 'delta t =', params%deltat
print*, 'endtime = ',params%tfinal
print*,  'frequency of output =', params%freq
nsteps=params%tfinal/params%deltat
nout=nsteps/params%freq+1
print*, nsteps, ' timesteps with ',nout,' outputs'

print*, ''
print*, 'Fluvial erosion parameters: '
print*, 'n = ',params%n
print*, 'm = ',params%m
print*, 'Hacks law h = ',params%h
print*, 'Hacks Law ka = ',params%ka
print*, 'K = ',params%k_scalar1
print*, 'Precip = ',params%rainfall_height

print*, ''
print*, 'Channel head erosion parameters: '
print*, 'Xc = ',params%xc
print*, 'tan theta = ',params%tanthetac
print*, 'Diffusivity = ',params%diffusivity/2.
print*, 'Min erosion rate for diffusional channel head = ',params%min_erosion_rate

print*, ''
print*, 'Uplift Parameters: '
print*, 'Interior Uplift rate = ',params%uplift_scalar1
print*, 'Boundary Uplift rate = ',params%uplift_scalar2
print*, ''
print*, 'Adding and Removing Node parameters: '
print*, 'Maximum allowable channel length = ',params%lmax
print*, 'Minimum allowable channel head area = ',params%amin
print*, 'Max Divide length = ',params%ldivmax
print*, ''
print*, ''
print*, ''
! calculate at what length the channel slope will become steeper than the imposed channel head
! this is a function of incision rate which is specified in subroutine uplift_advect and precip given in find_precip
! and k given in initialize_geometry. 
! Give them again here as a temp variables just for this calc. This is done in SS here.

u=params%uplift_scalar1
p=params%rainfall_height
k=params%k_scalar1
slope=(u/(k*params%xc**(params%h*params%m)*p**params%m*params%ka**params%m))**(1./params%n)
slope=180*atan(slope)/3.141592
sloped=u*params%xc/(params%diffusivity/2.)
sloped_total=u*params%xc/(params%diffusivity)
sloped=180*atan(sloped)/3.141592
sloped_total=180*atan(sloped_total)/3.141592
if(slope.gt.(180*atan(params%tanthetac)/3.141592))then
print*, '!! warning channel head slope is less steep than channel'
endif
if(sloped_total.lt.slope)then
print*, 'Warning!!! diffusive channel head slope is less steep than channel'
endif
print*, 'max slope of channel: ', slope
print*, 'slope of channel head: ', (180*atan(params%tanthetac)/3.141592)
print*, 'max diffusive slope of channel head: ', sloped_total

! check the channel head area

area=params%ka*params%xc**params%h
if(params%amin.lt.area)print*,'Warning!!!! minimum channel head area too small'
print*, 'Channel Head Area = ',area/1.d6
print*, 'Remesh minimum area used = ',params%amin/1.d6

!check linearity when using transient divide solution
if (params%transient_divide.and.params%n.ne.1) then
   print*, 'Cant solve transient divide with n=',params%n
   print*, 'Use slope exponent n=1'
endif


call time_out ('initialize_parameters')

return

end subroutine initialize_parameters
