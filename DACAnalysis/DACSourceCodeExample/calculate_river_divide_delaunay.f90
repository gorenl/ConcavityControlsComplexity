subroutine calculate_river_divide_delaunay (geometry,rddelaunay)

! interface routine to the delaun routine written by Malcolm Sambridge
! calculates the Delaunay triangulation around a set of points
! the number of triangles is stored in delaunay%ntriangles
! the triangulation is stored in delaunay%icon and the 
! neighboouring triangle information in delaunay%neighbours

! This is a "utility" routine and should not be changed

use definitions

implicit none

type (geom) geometry
type (del) rddelaunay

real*8,dimension(:,:),allocatable::points
integer,dimension(:,:),allocatable::e,v
integer,dimension(:),allocatable::vis_tlist,vis_elist,add_tlist
integer num,numtri_max,numtri,nv_max,mode,nfirst,itstart
logical,dimension(:),allocatable::inactive,subset
real*8 eps
integer i

call time_in ('calculate_river_divide_delaunay')

num=geometry%nnode + geometry%ndivide
numtri_max=num*3
nv_max=num
mode=0
nfirst=1
itstart=1
eps=1.d-10

allocate (points(2,num))
allocate (e(3,numtri_max),v(3,numtri_max))
allocate (vis_tlist(nv_max),vis_elist(nv_max),add_tlist(nv_max))
allocate (inactive(num),subset(num))

do i=1,geometry%nnode
points(1,i)=geometry%x(i)
points(2,i)=geometry%y(i)
enddo
do i=1,geometry%ndivide
points(1,geometry%nnode+i)=geometry%xdiv(i)
points(2,geometry%nnode+i)=geometry%ydiv(i)
enddo
inactive=.FALSE.
subset=.TRUE.

call delaun (points,num,e,v,numtri,numtri_max, &
             vis_tlist,vis_elist,add_tlist,eps,nv_max, &
             mode,inactive,nfirst,itstart,subset)

if (rddelaunay%ntriangles.ne.0) deallocate (rddelaunay%icon,rddelaunay%neighbours)

rddelaunay%ntriangles=numtri
print*, 'number of rddelaunay triangles is ',rddelaunay%ntriangles
allocate (rddelaunay%icon(3,rddelaunay%ntriangles),rddelaunay%neighbours(3,rddelaunay%ntriangles))
rddelaunay%icon=v(:,1:numtri)
rddelaunay%neighbours=e(:,1:numtri)

deallocate (points,e,v,vis_tlist,vis_elist,add_tlist,inactive,subset)

call time_out ('calculate_river_divide_delaunay')

return

end subroutine calculate_river_divide_delaunay
