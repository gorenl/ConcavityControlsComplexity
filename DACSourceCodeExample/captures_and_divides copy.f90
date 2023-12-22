subroutine captures_and_divides (geometry,network,params,delaunay,stack)

! Compute the capture of nodes using Sean Willett's steady-state analytical solution

! I have also coded the algorithm for the computation of divide position (and height)
! between node i and node j (i below j)
! From that position, one computes an array called surface_share that contains the 
! the relative position of the divide along the line connecting the two points
! with the convention:
!   - divide in the middle : surface_share=0
!   - divide at node j : surface_share=1
!   - divide at node i : surface_share=-1
! this array is then used in the subroutine find_surface to adjust the surfaces
! according to the divide positions

! Note that three flags control the behaviour of this routine (they are set in initialize_parameters)
!   params%capture determines whether captures are computed
!   params%divide determines whether divide position are computed
!   params%small_divide determines whether small divide (ie near one of the 2 points) positions are computed

use definitions

implicit none

type (geom) geometry
type (netw) network
type (parm) params
type (del) delaunay
type (stck) stack

integer i,j,k,ncapture,iter,ndivide,nsmall_divide,capnode(geometry%nnode),maxpass,iii,icap
double precision l,fact1,fact2,fact3,fact4,zt,nx,ny,xdi,xdj,fx,dfdx,zdi,xx,capelev(geometry%nnode)
integer temp,ii,ishuffle(2*geometry%nnode),ncapglobe

call time_in ('captures_and_divides')
maxpass=40
ncapture=0
ncapglobe=0
ndivide=0
nsmall_divide=0
geometry%surface_share=0.d0

! first checks whether the river is not flowing up-hill

do i=1,geometry%nnode
  if (network%receiver(i).ne.0) then
  j=network%receiver(i)
   if (geometry%z(j).gt.geometry%z(i)) then
    print*,'flowing uphill',i,j
!    print*,geometry%x(i),geometry%y(i),geometry%z(i),geometry%erosion_rate(i)
!    print*,geometry%x(j),geometry%y(j),geometry%z(j),geometry%erosion_rate(j)
    network%receiver(i)=0
    ncapture=ncapture+1
    ncapglobe=ncapglobe+1
   endif
 endif
enddo

! check whether captures happen
! first shuffle node order
do i=1,geometry%nnode
  call random_number(xx)
  j=1+int(xx*(geometry%nnode-1))
  if (j.gt.0) then
  temp=i
  ishuffle(i)=j
  ishuffle(j)=temp
  endif
enddo

if (params%hmn.eq.0.d0) then ! first in case where hm/n=1
! loop for passes until no captures occur
do iii=1, maxpass
ncapture=0
  do ii=1,geometry%nnode ! loop over the nodes
   i=ishuffle(ii)
   if (geometry%fix(i).ne.1) then
    do k=geometry%nb(i),1,-1 ! loop over all neighbours
    j=geometry%nn(k,i)
      if (network%receiver(i).ne.j .and. network%receiver(j).ne.i) then ! only for connections that are not part of the drainage network
        if (geometry%z(i).lt.geometry%z(j)) then ! check that i is lower than j
        l=sqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length of segment betwenn i and j
        fact1=(geometry%erosion_rate(i)/geometry%k(i)/ &
              (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
        fact2=log(params%xc/l)
        zt=geometry%z(i)+params%xc*params%tanthetac-fact1*fact2 ! zt is test elevation
          if (zt.lt.geometry%z(j)) then ! first case : capture
            if (params%capture) then
            ncapture=ncapture+1
            ncapglobe=ncapglobe+1
            network%receiver(j)=i
!            geometry%z(j)=zt
            capnode(ncapture)=j
            capelev(ncapture)=zt
            endif
          elseif (zt.lt.geometry%z(j)+params%xc*params%tanthetac) then ! second case : small divide
            if(iii.eq.1)then
               if (params%small_divide) then
               nsmall_divide=nsmall_divide+1 ! keeps track how many small divide are computed
               nx=(geometry%x(j)-geometry%x(i))/l ! nx,ny is the direction of the segment (normalized)
               ny=(geometry%y(j)-geometry%y(i))/l
               xdi=l-(zt-geometry%z(j))/params%tanthetac ! xdi is distance to divide
               geometry%surface_share(k,i)=xdi/l*2.d0-1.d0
               zdi=geometry%z(i)+params%xc*params%tanthetac & ! zdi is divide height
                  -fact1*log(params%xc/xdi)
               endif
            endif
          else ! third case : normal divide
            if(iii.eq.1)then
              if (params%divide) then
              ndivide=ndivide+1 ! keeps track of how many divides are computed
              nx=(geometry%x(j)-geometry%x(i))/l ! nx,ny is the direction of the segment (normalized)
              ny=(geometry%y(j)-geometry%y(i))/l
              xdi=params%xc ! first guess for xdi
              iter=0
	          fx=l
	          dfdx=1.d0
                do while (abs(fx/dfdx).gt.l/1.d6) ! Newton-Raphson iterative algorithm to solve simultaneous nonlinear equations
                iter=iter+1
                  if (iter.gt.100) then ! in case of poor convergence output to the screen
                  print*,'l',l
                  print*,'z',geometry%z(i),geometry%z(j)
                  print*,'erosionrate',geometry%erosion_rate(i),geometry%erosion_rate(j)
                  print*,'precipitation',geometry%precipitation(i),geometry%precipitation(j)
                  print*,'k',geometry%k(i),geometry%k(j)
                  print*,'xdi,fx,dfdx',xdi,fx,dfdx
                  stop 'no convergence in finding divide'
                  endif
                xdj=l-xdi
                fact2=log(params%xc/xdi)
                fact3=(geometry%erosion_rate(j)/geometry%k(j)/ &
                      (geometry%precipitation(j)*params%ka)**params%m)**(1.d0/params%n)
                fact4=log(params%xc/xdj)
                fx=geometry%z(j)-geometry%z(i)+(fact3*fact4-fact1*fact2)
                dfdx=-(fact1/xdi+fact3/xdj)
                xdi=xdi-fx/dfdx
                enddo
              geometry%surface_share(k,i)=xdi/l*2.d0-1.d0
              zdi=geometry%z(i)+params%xc*params%tanthetac-fact1*log(params%xc/xdi) ! zdi is divide height
            endif
          endif
        endif
      endif
      endif
    enddo
   endif
  enddo
if(ncapture.eq.0)go to 1492
do icap=1,ncapture
geometry%z(capnode(icap))=capelev(icap)
enddo
!print*, iii,ncapture, ' captures'
enddo
1492 continue
else ! first in case where hm/n is not =1
! this section has been used/tested (FH)
  do i=1,geometry%nnode ! loop over the nodes
    do k=1,geometry%nb(i) ! loop over all neighbours
    j=geometry%nn(k,i)
      if (network%receiver(i).ne.j .and. network%receiver(j).ne.i) then ! only for connections that are not part of the drainage network
        if (geometry%z(i).lt.geometry%z(j)) then ! check that i is lower than j
        l=sqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length of segment betwenn i and j
        fact1=(geometry%erosion_rate(i)/geometry%k(i)/ &
              (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
        fact2=l**params%hmn*(1.d0-(params%xc/l)**params%hmn)/params%hmn
        zt=geometry%z(i)+params%xc*params%tanthetac+fact1*fact2 ! zt is test elevation
          if (zt.lt.geometry%z(j)) then ! first case : capture
            if (params%capture) then
            ncapture=ncapture+1
            network%receiver(j)=i
            geometry%z(j)=zt
            endif
          elseif (zt.lt.geometry%z(j)+params%xc*params%tanthetac) then ! second case : small divide
            if (params%small_divide) then
            nsmall_divide=nsmall_divide+1 ! keeps track how many small divide are computed
            nx=(geometry%x(j)-geometry%x(i))/l ! nx,ny is the direction of the segment (normalized)
            ny=(geometry%y(j)-geometry%y(i))/l
            xdi=l-(zt-geometry%z(j))/params%tanthetac ! xdi is distance to divide
            geometry%surface_share(k,i)=xdi/l*2.d0-1.d0
            zdi=geometry%z(i)+params%xc*params%tanthetac & ! zdi is divide height
               +fact1*xdi**params%hmn*(1.d0-(params%xc/xdi)**params%hmn)/params%hmn
            endif
          else ! third case : normal divide
            if (params%divide) then
            ndivide=ndivide+1 ! keeps track how many divide are computed
            nx=(geometry%x(j)-geometry%x(i))/l ! nx,ny is the direction of the segment (normalized)
            ny=(geometry%y(j)-geometry%y(i))/l
            xdi=params%xc ! first guess for xdi
            iter=0
            fx=l
            dfdx=1.d0
              do while (abs(fx/dfdx).gt.l/1.d6) ! Newton-Raphson iterative algorithm to solve simultaneous nonlinear equations
              iter=iter+1
                if (iter.gt.100) then ! in case of poor convergence output to the screen
                print*,'Nodes',i,j
                print*,'l, xdi',l,xdi
                print*,'z',geometry%z(i),geometry%z(j)
                print*,'erosionrate',geometry%erosion_rate(i),geometry%erosion_rate(j)
                print*,'precipitation',geometry%precipitation(i),geometry%precipitation(j)
                print*,'k',geometry%k(i),geometry%k(j)
                print*,fact1,fact2,fact3,fact4
                print*,'xdi,fx,dfdx',xdi,fx,dfdx
                stop 'no convergence in finding divide'
                endif
              xdj=l-xdi
              fact2=xdi**params%hmn*(1.d0-(params%xc/xdi)**params%hmn)/params%hmn
              fact3=(geometry%erosion_rate(j)/geometry%k(j)/ &
                    (geometry%precipitation(j)*params%ka)**params%m)**(1.d0/params%n)
              fact4=xdj**params%hmn*(1.d0-(params%xc/xdj)**params%hmn)/params%hmn
              fx=geometry%z(j)-geometry%z(i)+(fact3*fact4-fact1*fact2)
              dfdx=-(fact1*xdi**(params%hmn-1.d0)+fact3*xdj**(params%hmn-1.d0))
              xdi=xdi-fx/dfdx
              enddo
            geometry%surface_share(k,i)=xdi/l*2.d0-1.d0
            zdi=geometry%z(i)+params%xc*params%tanthetac & ! zdi is divide height
               +fact1*xdi**params%hmn*(1.d0-(params%xc/xdi)**params%hmn)/params%hmn
            endif
          endif
        endif
      endif
    enddo
  enddo
endif

! print result to screen

if (ndivide.ne.0) then
!print*,ndivide,'divides'
!print*,nsmall_divide,'small divides'
geometry%ndivide=ndivide
!print*,'surface_share',minval(geometry%surface_share),maxval(geometry%surface_share)
endif

! print result to screen and adapt network to reflec stream captures

if (params%capture) then
!print*,ncapglobe,'captures'
geometry%ncapture=ncapglobe
network%donors=0
network%ndon=0
  do i=1,network%nnode
  k=network%receiver(i)
    if (k.ne.0) then
    network%ndon(k)=network%ndon(k)+1
    network%donors(network%ndon(k),k)=i
    endif
  enddo
!call show (geometry,delaunay,network,stack,params)
endif

call time_out ('captures_and_divides')

return

end subroutine captures_and_divides
