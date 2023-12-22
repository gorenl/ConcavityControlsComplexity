subroutine uplift_and_advect (geometry,params)

  ! Update the node geometry (for both horizontal and vertical tectonic advection)
  ! Note that horizontal advection may be turned off by setting the params%move_points flag to .FALSE.
  ! in the initialize_parameters routine

  use definitions

  implicit none

  type (geom) geometry
  type (parm) params
  integer i
  double precision y1,y2,y3,y4,x1,x2
 

  call time_in ('uplift_and_advect')

  y1 = minval(geometry%y)
  y2 = maxval(geometry%y)

  !constant rate of uplift
  geometry%w=params%uplift_scalar1
  geometry%u = 0.
  geometry%v = 0.

  ! use a different value for boundary nodes
  ! this controls captures aroundn boundary
  ! a value of zero causes many captures to boundary, but creates oscilations by removing length dependence of divide
  ! a larger value gives fewer




  do i=1,geometry%nnode
     if(geometry%fix(i).ne.1) geometry%z(i)=geometry%z(i)+params%deltat*geometry%w(i)
     if(params%move_points) then 
        geometry%x(i)=geometry%x(i)+params%deltat*geometry%u(i) 
        geometry%y(i)=geometry%y(i)+params%deltat*geometry%v(i)
     endif
  enddo

  call time_out ('uplift_and_advect')

  return

end subroutine uplift_and_advect
