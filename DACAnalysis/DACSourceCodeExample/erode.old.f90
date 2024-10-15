subroutine erode (geometry,network,stack,delaunay,params)

  ! Subroutine to compute fluvial erosion using a fully implicit method
  ! for n=1 or n=2
  ! Other cases not considered yet

  ! see attached document for the exact algorithms used

  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (stck) stack
  type (parm) params
  type (del) delaunay
  integer i,is,j,lakeflag,icount, jcount,dep_nnode, k, l, m, next_node_ind,loc_don,loc_rec,n,flag_not_in_list
  double precision fact,deltal,factp,erosion,sedheight,zi,zj,areai,areaj,length,depgrad,delz, erosion2, slopek
  double precision rtsafe, sedV, next_node_z, sumS
  integer, dimension(:),allocatable::dep_index
  double precision, dimension(:),allocatable::ztest
  double precision, dimension(:,:),allocatable::A
  double precision, dimension(:),allocatable::b

  call time_in ('erode')

  ! set a depositional gradient for basins

  depgrad=.1d-9
  !depgrad=.1d-5

  ! check for rivers flowing uphill
  do i = 1, geometry%nnode
     if(network%receiver(i).ne.0.and.geometry%fix(i).ne.1)then 
        if(geometry%z(i).lt.geometry%z(network%receiver(i)))then
           print*, 'Warning River flowing uphill in erode.f90',geometry%z(i),geometry%z(network%receiver(i))
        endif
     endif
  enddo



  !NEW ALGORITHM FOR DEPOSITION IN LAKES. L.G.


  do i=1,geometry%nnode
     if(network%receiver(i).eq.0.and.geometry%fix(i).ne.1)then   !find lake nodes
        sedV=params%deltat*geometry%sediment_flux(i) !volume of sediments to be deposited in that lake
        if (sedV.eq.0.) cycle
        !first try to deposit within the node itself
        allocate(dep_index(geometry%nnode),ztest(geometry%nnode))
        dep_nnode=1
        dep_index(dep_nnode)=i
        ztest(dep_nnode)=geometry%z(i) + sedV/geometry%surface(i)    

152     if (dep_nnode.gt.500) pause 'too many depositing nodes'
        do j=1,dep_nnode !loop over the nodes that are being sedimenting
           k=dep_index(j) 
           if (network%ndon(k).ne.0) then
              do l=1,network%ndon(k) !loop over their donors
                 m=network%donors(l,k)
                 flag_not_in_list=1
                 do n=j,dep_nnode
                    if (dep_index(n).eq.m) flag_not_in_list=0
                 enddo                    
                 if (flag_not_in_list.eq.1.and.geometry%z(m).lt.ztest(j)) then !if a donor is lower than a test elevation of its receiver, need to redistribute the sediments.
                    !print*, 'the donor',m,'with elev',geometry%z(m),'is lower than receiver',dep_index(j),'with test elev',ztest(j)
                    !print*, dep_index(j) ,'should be equal to',network%receiver(m), 'and to',k
                    go to 151 !add a new node to the sedimenting node and test new elevations
                 endif
              enddo
           endif
        enddo
        go to 153 !if the loop was exited normally it means that there are no reverse gradients

151     next_node_z=maxval(geometry%z)+1.
        next_node_ind=-1
        do j=1,dep_nnode !add to the arrays the next node according to elevation order
           k=dep_index(j)
           if (network%ndon(k).ne.0) then
              do l=1,network%ndon(k) !loop over donors
                 m=network%donors(l,k)
                 if (geometry%z(m).lt.next_node_z)then
                    !check that it is not already in the list
                    flag_not_in_list=1
                    do n=j,dep_nnode
                       if (dep_index(n).eq.m) flag_not_in_list=0
                    enddo
                    if (flag_not_in_list.eq.1) then
                       next_node_z = geometry%z(m)
                       next_node_ind = m
                    endif
                 endif
              enddo
           endif
        enddo
        !after the loop we have found the lowest donor. Now we can add it to the list
        dep_nnode=dep_nnode+1
        dep_index(dep_nnode)=next_node_ind
        !print*, 'adding node',next_node_ind,'to the list'
        !build new Ax=b problem
        if (dep_nnode-1.ne.1) deallocate(A,b)
        allocate(A(dep_nnode,dep_nnode),b(dep_nnode))
        A=0.
        b=0.
        !build first row
        do j=1,dep_nnode 
           A(1,j)=geometry%surface(dep_index(j))
        enddo
        b(1)=sedV
        !built next rows
        do j=2,dep_nnode
           loc_don=dep_index(j)
           loc_rec=network%receiver(loc_don)
           do k=1,j !find the reciever in the sedimenting nodes array
              if (dep_index(k).eq.loc_rec) then
                 A(j,j)=1.
                 A(j,k)=-1.
                 length=dsqrt((geometry%x(loc_don)-geometry%x(loc_rec))**2+(geometry%y(loc_don)-geometry%y(loc_rec))**2)
                 b(j)=depgrad*length - geometry%z(loc_don) + geometry%z(loc_rec)
                 exit
              endif
           enddo
        enddo
        !at this stage A and b should be full

        !print them for debug purpose
        !print*, 'A is a matrix of size',dep_nnode,dep_nnode
!!$        print*, 'A=['
!!$        do j=1,dep_nnode
!!$           write(*,100)(A(j,k), k=1,dep_nnode)
!!$        enddo
!!$        print*,']'
!!$        print*,' b=['
!!$        write(*,100) (b(j),j=1,dep_nnode)
!!$        print*,']'
!!$100     format(F4.2)

        call gaussj(A,dep_nnode,dep_nnode,b,1,1)
        !now b contains the result for the delz that should be sedimented in each node
        do j=1,dep_nnode
           ztest(j)=geometry%z(dep_index(j))+b(j)
        enddo
        go to 152 !check again that the gradients have note reversed

153     sumS=0. !for debuging to make sure that the correct volume is being sedimented
        do  j=1,dep_nnode
           sumS=sumS+(ztest(j)-geometry%z(dep_index(j)))*geometry%surface(dep_index(j)) 
        enddo
        !if (sumS.ne.sedV) then
           !print*, 'Available sediments volume:', sedV
           !print*, 'trying to deposit volume of:',sumS
        !endif
        if (dep_nnode.gt.100) then
           print*, 'depositing', sedV, 'sediments'
           print*, 'over',dep_nnode,'nodes'
        endif
        do j=1,dep_nnode !here the new elevation is being updated
           geometry%z(dep_index(j)) = ztest(j)
        enddo
        deallocate(dep_index,ztest)
        if (dep_nnode.gt.1) deallocate(A,b)
     endif
  enddo











  !
  ! Fill lakes algorithm -sdw-15.06.2010
  !    first deposit all sediment at low point of drainage basin
  !
!!$  do i=1,geometry%nnode
!!$     if(network%receiver(i).eq.0.and.geometry%fix(i).ne.1)then   !find lake nodes
!!$        sedheight=params%deltat*geometry%sediment_flux(i)/geometry%surface(i)  !maximum sed available
!!$        geometry%z(i)=geometry%z(i)+sedheight
!!$     endif
!!$  enddo

  !    then pass through stack looking for negative gradients and diffusing them until water flows down original path

!!$  lakeflag=1
!!$  icount=0
!!$  do while (lakeflag.eq.1)
!!$     icount=icount+1
!!$     lakeflag=0
!!$     do is=1,stack%nnode
!!$        i=stack%order(is)
!!$        j=network%receiver(i)
!!$        if(j.ne.0.and.geometry%fix(i).ne.1)then   
!!$           if((geometry%z(j).gt.geometry%z(i)))then    !water now flowing uphill
!!$              lakeflag=1
!!$
!!$              !  correct height of both points so that sediment surface takes a specified gradient downhill
!!$
!!$              length=dsqrt((geometry%x(j)-geometry%x(i))**2+(geometry%y(j)-geometry%y(i))**2)
!!$              delz=length*depgrad/2.
!!$              zj=geometry%z(j)
!!$              zi=geometry%z(i)
!!$              areaj=geometry%surface(j)
!!$              areai=geometry%surface(i)
!!$              geometry%z(i)=(areaj*(zj+delz)+areai*(zi-delz))/(areaj+areai)+delz
!!$              geometry%z(j)=geometry%z(i)-2*delz
!!$           endif
!!$        endif
!!$     enddo
!!$     if(icount.gt.5000)then
!!$        print*, 'Warning: icount in lake fill loop exceeded'
!!$        go to 10
!!$     endif
!!$  enddo
!!$10 continue



  if (params%n.eq.1) then
     do is=1,stack%nnode
        i=stack%order(is)
        if (geometry%fix(i).ne.1 .and. network%receiver(i).ne.0) then
           j=network%receiver(i)
           deltal=sqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2)
           if ((geometry%z(i)-geometry%z(j))/deltal.le.depgrad)then !no erosion its a lake
              erosion=geometry%z(i)
           else
              fact=params%deltat*geometry%k(i)/deltal*geometry%discharge(i)**params%m
              if (fact.gt.1.d0) then
                 !             print*, 'erosion factor =',fact,'>1'
              endif
!!$           if (fact.gt.1.d0) then
!!$              print*,'Time step too large'
!!$              print*,'maximum value accepted is ',params%deltat/fact
!!$              print*,'Parameter values are :'
!!$              print*,'k:',geometry%k(i)
!!$              print*,'discharge:',geometry%discharge(i)
!!$              print*,'n:',params%n
!!$              print*,'m:',params%m
!!$              print*,'deltat:',params%deltat
!!$              print*,'Node geometry (',i,')'
!!$              print*,'xi:',geometry%x(i)
!!$              print*,'yi:',geometry%y(i)   
!!$              print*,'zi:',geometry%z(i)    
!!$              print*,'fixi',geometry%fix(i)
!!$              print*,'Node geometry (',j,')'
!!$              print*,'xj:',geometry%x(j)
!!$              print*,'yj:',geometry%y(j)
!!$              print*,'zj:',geometry%z(j)
!!$              print*,'fixj',geometry%fix(j)
!!$              print*,'slope',(geometry%z(i)-geometry%z(j))/deltal
!!$              stop
!!$           endif
           ! "erosion" is the new height of node i
              erosion=(geometry%z(i)+fact*geometry%z(j))/(1.d0+fact)
           ! the next 6 or so lines are patches to catch errors in gradient or overextrapolation of rate..
           endif
           if (geometry%z(i).lt.geometry%z(j)) then
              !    geometry%z(i)=geometry%z(j)
              print*, 'error in erode.f90 -reciever higher than donor',i,j
           endif
           !    if(erosion.gt.geometry%z(i))then
           !       print*, 'erosion greater than previous elevation'
           erosion=min(erosion,geometry%z(i))
           !    endif
           if (erosion.le.geometry%z(j)) then
              !      print*, 'erosion less than zj'
              erosion=geometry%z(j)
           endif
           !  end checks
           geometry%erosion_rate(i)=-(erosion-geometry%z(i))/params%deltat

           if(geometry%erosion_rate(i).lt.0)then
              geometry%erosion_rate(i)=0.
           else
              geometry%z(i)=erosion
           endif
        endif
     enddo
  elseif (params%n.eq.2) then
     do is=1,stack%nnode
        i=stack%order(is)
        if (geometry%fix(i).ne.1 .and. network%receiver(i).ne.0) then
           j=network%receiver(i)
           deltal=sqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2)
           if ((geometry%z(i)-geometry%z(j))/deltal.le.depgrad)then !no erosion its a lake
              erosion=geometry%z(i)
           else
              fact=params%deltat*geometry%k(i)/deltal**params%n*geometry%discharge(i)**params%m
              factp=(1.d0-2.d0*fact*geometry%z(j))**2-4.d0*fact*(fact*geometry%z(j)**2-geometry%z(i))           
              if (factp.le.0.d0) stop 'factp negative in erode.f90'
              geometry%erosion_rate(i)=0.d0
              if (fact.ne.0.d0) then
                 erosion = (2*geometry%z(j)*fact - 1. + sqrt(factp))/2.d0/fact !first root
                 !erosion2 = (2*geometry%z(j)*fact - 1. - sqrt(factp))/2.d0/fact !second root  
              endif
           endif
           if (erosion .gt. geometry%z(i)) then ! I think that due to tiny numerical errors
              geometry%erosion_rate(i) = 0.0d0
              print*, 'erosion causes the node to rise'
           else
              geometry%erosion_rate(i) = (geometry%z(i) - erosion)/params%deltat
              geometry%z(i)=erosion
           endif
        endif
     enddo
  else !n > 2
     do is=1,stack%nnode
        i=stack%order(is)
        if (geometry%fix(i).ne.1 .and. network%receiver(i).ne.0) then
           j=network%receiver(i)           
           deltal=sqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2)
           if ((geometry%z(i)-geometry%z(j))/deltal.le.depgrad)then !no erosion its a lake
              erosion=geometry%z(i)
           else
              fact = params%deltat*geometry%k(i)*geometry%discharge(i)**params%m 
              erosion = rtsafe(geometry%z(i),geometry%z(j),1.0d-10,deltal,fact,params%n)

!!!! checking that all is normal
!!$           if ((params%istep/params%freq)*params%freq.eq.params%istep) then
!!$              print*, i,geometry%x(i),geometry%y(i),geometry%z(i)
!!$              print*, j,geometry%x(j),geometry%y(j),geometry%z(j)
!!$              print*, deltal,geometry%k(i),geometry%discharge(i),params%m,params%n
!!$              print*, params%deltat,erosion,(geometry%z(i) - erosion)/params%deltat
!!$              print*, '^^^^^^^^'
!!$           endif
           endif
           if (erosion.lt.geometry%z(j))then
              print*, 'source is now lower than reciever'
           endif

           if (erosion .gt. geometry%z(i)) then ! I think that due to tiny numerical errors
              geometry%erosion_rate(i) = 0.0d0
              print*, 'erosion causes the node to rise'
           else
              geometry%erosion_rate(i) = (geometry%z(i) - erosion)/params%deltat
              geometry%z(i)=erosion
           endif
        endif
     enddo
  endif

  ! set erosion rate to zero for nodes with no receiver (lakes)
  ! set the erosion rate for boundary nodes to specified value
  !   this is redundent - done above already
  do i=1, geometry%nnode
     if(network%receiver(i).eq.0) geometry%erosion_rate(i)=0.0
     if(geometry%fix(i).eq.1) geometry%erosion_rate(i)=geometry%w(i)



  enddo



  call time_out ('erode')

  return

end subroutine erode


!!!!!!!!!!!!!!!!!! subroutine rtsafe Newton-Raphson Method with bisection !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! Based on Numerical Recipes in Fortran (page 359)     !!!!!!!!!!!!!!!!!!

function rtsafe(x1,x2,xacc,deltal,fact,n)

integer MAXIT
double precision rtsafe, x1, x2, xacc, deltal, fact
integer n
!EXTERNAL funcd
PARAMETER (MAXIT = 100)
integer k
double precision df,dx,dxold,f,fh,fl,temp,xh,xl


call funcd(x1,fl,df,deltal,x1,x2,fact,n)
call funcd(x2,fh,df,deltal,x1,x2,fact,n)
if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.))&
stop 'root must be bracketed in rtsafe'
if (fl.eq.0.)then
   rtsafe = x1
   return
else if (fh.eq.0.)then
   rtsafe = x2
   return
else if (fl.lt.0.) then
   xl = x1
   xh = x2
else
   xh = x1
   xl = x2
endif
rtsafe = .5*(x1+x2)
dxold = abs(x2-x1)
dx=dxold
call funcd(rtsafe,f,df,deltal,x1,x2,fact,n)
do k= 1,MAXIT
   if (((rtsafe-xh)*df - f)*((rtsafe - xl)*df - f).gt.0.&
        .or. abs(2.*f).gt.abs(dxold*df)) then
      dxold = dx
      dx = 0.5*(xh-xl)
      rtsafe = xl+dx
      if(xl.eq.rtsafe) return
   else
      dxold=dx
      dx=f/df
      temp=rtsafe
      rtsafe = rtsafe-dx
      if (temp.eq.rtsafe) return
   endif
   if (abs(dx).lt.xacc) return
   call funcd(rtsafe,f,df,deltal,x1,x2,fact,n)
   if(f.lt.0.) then
      xl = rtsafe
   else
      xh = rtsafe
   endif
!   print*, 'at time',k,'guess is', rtsafe, dx, xacc
enddo
stop 'rtsafe exceeding maximum iterations'
return
end FUNCTION rtsafe


!!!!!!!!!!!! subroutine funcd !!!!!!!!!!!!!!!!!
!!! returns the erosion function and its derivative at point x

subroutine funcd (x,f,df,deltal,zi,zj,fact,n)
 
  use definitions
  implicit none
  double precision x,f,df,deltal,zi,zj,fact
  integer n
  f = x - zi + fact*((x-zj)/deltal)**n
  df = 1 + fact*n*(1/deltal)*((x-zj)/deltal)**(n-1)
  return 
end subroutine funcd
