subroutine show (geometry,delaunay,network,stack,params)

! subroutine to display results of computations

use definitions

implicit none

type (del) delaunay
type (geom) geometry
type (netw) network
type (stck) stack
type (parm) params
integer window(4),red(230),green(230),blue(230)
integer ichoice,xchoice,it,k,i,ic,icp,flag,ncontours
external xchoice
real x,y,xl,yl,rad,fact
character d*3
double precision val,v1,v2,rat

ncontours=10
if (ncontours.gt.20) stop 'not enough color slots'

window(1)=100
window(2)=100
window(3)=600
window(4)=600
red=0
green=0
blue=0

red(2)=255
green(3)=255
blue(4)=255
red(5)=255
blue(5)=255

  do i=101,100+ncontours/2
  fact=float(i-101)/(ncontours/2)
  red(i)=int(255*fact)
  green(i)=255
  blue(i)=0
  enddo
  do i=101+ncontours/2,100+ncontours
  fact=float(i-(101+ncontours/2))/(ncontours/2)
  red(i)=255
  green(i)=255-int(255*fact)
  blue(i)=0
  enddo

  do i=101+ncontours,100+ncontours/2+ncontours
  fact=float(i-101+ncontours)/(ncontours/2)
  red(i)=255
  green(i)=0
  blue(i)=int(255*fact)
  enddo
  do i=101+ncontours+ncontours/2,100+ncontours+ncontours
  fact=float(i-(101+ncontours+ncontours/2))/(ncontours/2)
  red(i)=255-int(255*fact)
  green(i)=0
  blue(i)=255
  enddo

xl=maxval(geometry%x)
yl=maxval(geometry%y)

call xopen (window,char(0))
call xcmap (red,green,blue)
call xscale (0,window(3),0,window(4),-xl*0.05,xl*1.05,yl*1.05,-yl*0.05,1)

if (params%plot_height) then
v1=minval(geometry%z)
v2=maxval(geometry%z)
print*,'Topo min-max',v1,v2
  do it=1,delaunay%ntriangles
    do i=1,ncontours
    val=v1+v2*float(i)/(ncontours+1)
    flag=0
      do k=1,3
      ic=delaunay%icon(k,it)
      icp=delaunay%icon(mod(k,3)+1,it)
        if ((geometry%z(ic)-val)*(geometry%z(icp)-val).lt.0.d0) then
        rat=geometry%z(icp)-geometry%z(ic)
          if (rat.ne.0.d0) then
          rat=(val-geometry%z(ic))/rat
          x=geometry%x(ic)+rat*(geometry%x(icp)-geometry%x(ic))
          y=geometry%y(ic)+rat*(geometry%y(icp)-geometry%y(ic))
            if (flag.eq.0) then
            call xplot (x,y,3)
            flag=1
            else
            call xpen (100+i)
            call xplot (x,y,2)
            flag=0
            endif
          endif
        endif
      enddo
    enddo
  enddo
endif

if (params%plot_precipitation) then
v1=minval(geometry%precipitation)
v2=maxval(geometry%precipitation)
print*,'Precip min-max',v1,v2
  do it=1,delaunay%ntriangles
    do i=1,ncontours
    val=v1+v2*float(i)/(ncontours+1)
    flag=0
      do k=1,3
      ic=delaunay%icon(k,it)
      icp=delaunay%icon(mod(k,3)+1,it)
        if ((geometry%precipitation(ic)-val)*(geometry%precipitation(icp)-val).lt.0.d0) then
        rat=geometry%precipitation(icp)-geometry%precipitation(ic)
          if (rat.ne.0.d0) then
          rat=(val-geometry%precipitation(ic))/rat
          x=geometry%x(ic)+rat*(geometry%x(icp)-geometry%x(ic))
          y=geometry%y(ic)+rat*(geometry%y(icp)-geometry%y(ic))
            if (flag.eq.0) then
            call xplot (x,y,3)
            flag=1
            else
            call xpen (100+ncontours+i)
            call xplot (x,y,2)
            flag=0
            endif
          endif
        endif
      enddo
    enddo
  enddo
endif

if (params%plot_triangles) then
call xpen (1)
  do it=1,delaunay%ntriangles
  x=geometry%x(delaunay%icon(3,it))
  y=geometry%y(delaunay%icon(3,it))
  call xplot (x,y,3)
    do k=1,3
    x=geometry%x(delaunay%icon(k,it))
    y=geometry%y(delaunay%icon(k,it))
    call xplot (x,y,2)
    enddo
  enddo
endif

if (params%plot_receiver) then
call xpen (4)
  do i=1,network%nnode
    if (network%receiver(i).ne.0) then
    x=geometry%x(i)
    y=geometry%y(i)
    call xplot (x,y,3)
    x=geometry%x(network%receiver(i))
    y=geometry%y(network%receiver(i))
    call xplot (x,y,2)
    endif
  enddo
endif

if (params%plot_donors) then
call xpen (1)
  do i=1,network%nnode
    do k=1,network%ndon(i)
    x=geometry%x(i)
    y=geometry%y(i)
    call xplot (x,y,3)
    x=geometry%x(network%donors(k,i))
    y=geometry%y(network%donors(k,i))
    call xplot (x,y,2)
    enddo
  enddo
endif

if (params%plot_no_receiver) then
call xpen (3)
rad=xl/100.
  do i=1,network%nnode
    if (network%receiver(i).eq.0) then
    x=geometry%x(i)
    y=geometry%y(i)
    call xcirclef (x,y,rad)
    endif
  enddo
endif

if (params%plot_no_donors) then
call xpen (4)
rad=xl/100.
  do i=1,network%nnode
    if (network%ndon(i).eq.0) then
    x=geometry%x(i)
    y=geometry%y(i)
    call xcirclef (x,y,rad)
    endif
  enddo
endif

if (params%write_discharge) then
call xpen (1)
  do i=1,geometry%nnode
  x=geometry%x(i)
  y=geometry%y(i)
  write (d,'(i3)') int(geometry%discharge(i))
  call xsymbol (x,y,1,d,3)
  enddo
endif

if (params%write_stack_order) then
call xpen (4)
  do i=1,stack%nnode
  x=geometry%x(stack%order(i))
  y=geometry%y(stack%order(i))
  write (d,'(i3)') i
  call xsymbol (x,y,1,d,3)
  enddo
endif

if (params%plot_catchment) then
call xpen (3)
  do it=1,delaunay%ntriangles
    do k=1,3
    ic=delaunay%icon(k,it)
    icp=delaunay%icon(mod(k,3)+1,it)
      if (geometry%catchment(ic).ne.geometry%catchment(icp)) then
      x=(geometry%x(ic)+geometry%x(icp))/2.d0
      y=(geometry%y(ic)+geometry%y(icp))/2.d0
      call xplot (x,y,3)
      x=sum(geometry%x(delaunay%icon(:,it)))/3.d0
      y=sum(geometry%y(delaunay%icon(:,it)))/3.d0
      call xplot (x,y,2)
      endif
    enddo
  enddo
endif

call xsave ()
ichoice=xchoice('Continue            ',1)
call xclose ()

return

end subroutine show
