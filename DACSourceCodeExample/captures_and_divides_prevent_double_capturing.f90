subroutine captures_and_divides (geometry,network,params,delaunay)
  
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
  

  ! capture related variables
  ! ncapture - counts the number of captures in a single swip over the nodes.
  ! ncapglobe - counts the number of captures in a single time step.
  ! geometry%ncapture - counts the comulative number of capture during the simulation. 

  use definitions
  
  implicit none
  
  type (geom) geometry
  type (netw) network
  type (parm) params
  type (del) delaunay
  
  integer i,j,k,index,ncapture,iter,ndivide,nsmall_divide,capnode(geometry%nnode),maxpass,iii,icap,icap2,jcap
  double precision l,fact1,fact2,fact3,fact4,zt,xdi,xdj,fx,dfdx,zdi,xx,capelev(geometry%nnode),ztprime,ltest,ztd,xd,x1,x2,ztoler,xc,zc,zdj
  double precision zi1,zi2,zj1,zj2,erodrate,avedivide
  integer temp,ii,ishuffle(2*geometry%nnode),ncapglobe,idon,recnode(geometry%nnode),icount,icountdif,icountslope,dflag
  double precision local_slope
  integer it, ncount, ntri
  logical nfound
  double precision alpha, theta, A, B, D, DD, delz, zd
  DOUBLE PRECISION, PARAMETER :: pi = 3.1415927

  integer,dimension(:),allocatable::orig_receiver
  double precision, dimension(:),allocatable::orig_z
  integer counter

  call time_in ('captures_and_divides')

  maxpass=500
  ncapture=0
  ncapglobe=0
  ndivide=0
  nsmall_divide=0
  geometry%surface_share=0.d0
  avedivide=0.0

  allocate (orig_receiver(geometry%nnode))
  allocate (orig_z(geometry%nnode))
  do i = 1,geometry%nnode
     orig_receiver(i) = network%receiver(i)
     orig_z(i) = geometry%z(i)
  enddo

  ! first checks whether the river is flowing up-hill
  ! if so, cut and make a lake, call this a capture
  
  do i=1,geometry%nnode
     if (network%receiver(i).ne.0.and.geometry%fix(i).ne.1) then
        j=network%receiver(i)
        if (geometry%z(j).gt.geometry%z(i)) then
           print*,'flowing uphill - cut river  make lake ',i,j
           !    print*,geometry%x(i),geometry%y(i),geometry%z(i),geometry%erosion_rate(i)
           !    print*,geometry%x(j),geometry%y(j),geometry%z(j),geometry%erosion_rate(j)
           !  cut river, make new lake
           ! set node i reciever to 0 and erosion rate to 0
           network%receiver(i)=0
           geometry%erosion_rate(i)=0.0 
           ncapglobe=ncapglobe+1
        endif
     endif
  enddo
  !
  !
  ! first shuffle node order
  do i=1, geometry%nnode
     ishuffle(i)=i
  enddo
  do i=1,geometry%nnode
     call random_number(xx)
     j=1+int(xx*(geometry%nnode-1))
     if(j.le.0)j=1
     if(j.gt.geometry%nnode)j=geometry%nnode
     call random_number(xx)
     k=1+int(xx*(geometry%nnode-1))
     if(k.le.0)k=1
     if(k.gt.geometry%nnode)k=geometry%nnode
     temp=ishuffle(j)
     ishuffle(j)=ishuffle(k)
     ishuffle(k)=temp
  enddo
  
  ! check that shuffle is complete
  
  !do i=1, geometry%nnode
  !  do j=1,geometry%nnode
  !   if(i.eq.ishuffle(j))go to 4444
  !  enddo
  !print*, 'missing node in shuffle', i
  !4444 continue
  !enddo
  
  ! check for captures
  
  counter = 0
  ncapture = 1
  if (params%hmn.eq.0.d0) then ! first in case where hm/n=1
     ! loop for passes until no captures occur
     do iii=1, maxpass
        ncapture=0
        do ii=1,geometry%nnode ! loop over the nodes
           i=ishuffle(ii)
           do k=geometry%nb(i),1,-1 ! loop over all neighbours
              j=geometry%nn(k,i)
              if(network%receiver(j).ne.i.and.network%receiver(i).ne.j) then ! only for connections that are not part of the drainage network
                 !if (geometry%z(i).lt.geometry%z(j)) then! check that i is lower than j
                 if (geometry%z(i).lt.geometry%z(j).and.orig_receiver(j).ne.i) then ! check that i is lower than j and than j didnt use to flow to i in this time step
                    l=sqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length betwenn i and j
                    if(l.gt.params%xc)then
                       fact1=(geometry%erosion_rate(i)/geometry%k(i)/ &
                       (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                       fact2=log(params%xc/l)
                       zt=params%xc*params%tanthetac ! zt is channel head elevation
                       if(params%diffusion)then
                          erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)
                          ztd=erodrate*params%xc**2/params%diffusivity
                          if(ztd.lt.zt)then
                             zt=ztd
                          endif
                       endif
                       zt=zt+geometry%z(i)-fact1*fact2 ! zt is test elevation
                    else
                       zt=l*params%tanthetac
                       if(params%diffusion)then
                          erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)
                          ztd=erodrate*l**2/params%diffusivity
                          if(ztd.lt.zt)then
                             zt=ztd
                          endif
                       endif
                       zt=zt+geometry%z(i)
                    endif
                    if(zt.lt.geometry%z(i))then
                       print*, 'ERROR - zt capturing a lower point'
                       print*,'l', l
                       print*,'erosion rate', geometry%erosion_rate(i)
                       print*,'fact1, fact2', fact1, fact2
                       print*,'zt', zt
                    endif
                    if (zt.lt.geometry%z(j)) then ! first case : capture
                       if (params%capture) then
                          if (orig_receiver(j).eq.i) then
                             print*, i, 'is capturing its former donor', j, 'FORBIDDEN'
                          endif
                          ncapture=ncapture+1
                          ncapglobe=ncapglobe+1
                          recnode(ncapture)=i
                          capnode(ncapture)=j
                          capelev(ncapture)=zt
                       endif
                    endif
                 endif
              endif
           enddo !end of loop over neighbors
        enddo ! end of loop over nodes
        if(ncapture.eq.0)go to 1492
        ! check for two nodes capturing the same node and select the lower capture elevation
        if(ncapture.gt.1) then
           do icap=1,ncapture
              do jcap=icap+1,ncapture
                 if(capnode(icap).eq.capnode(jcap))then
                    if(capelev(jcap).lt.capelev(icap))then 
                       capelev(icap)=geometry%z(capnode(icap))
                    else
                       recnode(jcap)=recnode(icap)
                       capnode(jcap)=capnode(icap)
                       capelev(jcap)=capelev(icap)
                       capelev(icap)=geometry%z(capnode(icap))
                    endif
                 endif
              enddo
           enddo
        endif
        do icap=1,ncapture
           geometry%erosion_rate(capnode(icap))=(geometry%z(capnode(icap))-capelev(icap))/params%deltat
           !      if(geometry%erosion_rate(capnode(icap)).lt.0)print*, 'negative erosion rate'
           geometry%z(capnode(icap))=capelev(icap)
           network%receiver(capnode(icap))=recnode(icap)
        enddo
        
        !if(ncapture.gt.00)print*, iii,ncapture, ' captures'
        !if (ncapture .ne. 0) print*, 'ncapture is', ncapture
     enddo ! end of loop of maxpasses
     print*, 'warning max passes exceeded'
1492 continue
     
     !if(ncapglobe.gt.5)print*, iii, ' passes with ', ncapglobe, ' total captures in timestep'
     !  end loop for captures
     !
     !         ******************  start loop for divide heights      **********************************
     !
     ! set a tolerance for the ridge height convergence
     ztoler=.1
     icountdif=0
     icountslope=0
     do i=1,geometry%nnode ! loop over the nodes
        do k=geometry%nb(i),1,-1 ! loop over all neighbours
           j=geometry%nn(k,i)
           if(network%receiver(j).ne.i.and.network%receiver(i).ne.j) then ! only for connections that are not part of the drainage network
              if (geometry%z(i).lt.geometry%z(j)) then ! check that i is lower than j
                 l=sqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length of segment betwenn i and j
                 ndivide=ndivide+1
                 avedivide=avedivide+l
                 ! start here with new diffusion algorithm
                 


                 if (orig_receiver(j).ne.i) then !if didn't use to be a channel in this time step
                    ! give initial guess for divide and endpoints
                    x1=0.
                    x2=l
                    zdi=0
                    zdj=2*ztoler
                    icount=0
                    ! find ridge heights at opposing node positions
                    !  find ridge height from node i
                    xc=min(params%xc,l)
                    erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)
                    if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/erodrate))then 
                       zc=erodrate*xc**2/params%diffusivity
                    else
                       zc=xc*params%tanthetac
                    endif
                    if(l.le.params%xc)then
                       zi2=geometry%z(i)+zc
                    else
                       fact2=log(params%xc/((l)))
                       fact3=(geometry%erosion_rate(i)/geometry%k(i)/ &
                       (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                       zi2=geometry%z(i)+zc-fact2*fact3
                    endif
                    zi1=geometry%z(i)
                    !  find ridge height from node j
                    xc=min(params%xc,(l))
                    erodrate=max(geometry%erosion_rate(j),params%min_erosion_rate)
                    if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/erodrate))then 
                       zc=erodrate*xc**2/params%diffusivity
                    else
                       zc=xc*params%tanthetac
                    endif
                    if((l).le.params%xc)then
                       zj1=geometry%z(j)+zc
                    else
                       fact2=log(params%xc/((l)))
                       fact3=(geometry%erosion_rate(j)/geometry%k(j)/ &
                       (geometry%precipitation(j)*params%ka)**params%m)**(1.d0/params%n)
                       zj1=geometry%z(j)+zc-fact2*fact3
                    endif
                    zj2=geometry%z(j)
                    if((zi1-zj1).ge.0.or.(zi2-zj2).le.0)then
                       print*, 'no root in initial conditions'
                       print*, icount, zi1,zj1,zi2,zj2
                    endif
                    ! iterate using false point method until divide height is found within specified tolerance
                    do while (dabs(zdi-zdj).gt.ztoler)
                       icount=icount+1
                       dflag=0
                       !  find root of function zdi-zdj
                       !  find new guess at divide position from false point method
                       !  unless ftc after 100 iterations, then revert to bisector method
                       if(icount.lt.100)then
                          xd=(x1*(zi2-zj2)-x2*(zi1-zj1))/(zi2-zj2-zi1+zj1)
                       else
                          xd=x1+.5*(x2-x1)
                       endif
                       if((zi1-zj1).ge.0.or.(zi2-zj2).le.0)then
                          print*, 'error no root found for divide problem'
                          print*, icount, zi1,zj1,zi2,zj2
                          xd=l/2.
                          go to 333
                       endif
                       
                       !  find ridge height from node i
                       xc=min(params%xc,xd)
                       erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)
                       if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/erodrate))then 
                          zc=erodrate*xc**2/params%diffusivity
                          dflag=1
                       else
                          zc=xc*params%tanthetac
                       endif
                       if(xd.le.params%xc)then
                          zdi=geometry%z(i)+zc
                       else
                          fact2=log(params%xc/((xd)))
                          fact3=(geometry%erosion_rate(i)/geometry%k(i)/ &
                          (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                          zdi=geometry%z(i)+zc-fact2*fact3
                       endif
                       !  find ridge height from node j
                       xc=min(params%xc,(l-xd))
                       erodrate=max(geometry%erosion_rate(j),params%min_erosion_rate)
                       if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/erodrate))then 
                          zc=erodrate*xc**2/params%diffusivity
                          dflag=1
                       else
                          zc=xc*params%tanthetac
                       endif
                       if((l-xd).le.params%xc)then
                          zdj=geometry%z(j)+zc
                       else
                          fact2=log(params%xc/((l-xd)))
                          fact3=(geometry%erosion_rate(j)/geometry%k(j)/ &
                          (geometry%precipitation(j)*params%ka)**params%m)**(1.d0/params%n)
                          zdj=geometry%z(j)+zc-fact2*fact3
                       endif
                       !  reset bounds of interval
                       if(zdj.gt.zdi)then
                          x1=xd
                          zi1=zdi
                          zj1=zdj
                       else
                          x2=xd
                          zi2=zdi
                          zj2=zdj
                       endif
                       !                 if(x1.eq.x2)print*, icount, x1,x2,geometry%z(i), geometry%z(j),zdi,zdj,l
                       !                if(x1.eq.x2)stop
                       !                 if(icount.gt.500)then
                       !                    !    print*, 'failed to find root at ', icount, zdi,zdj,xd,l
                       !                        write(*,222) icount, zdi,zdj,x1,x2, xd,l
                       !          222  format ('failed to find root at ',i5, 6e15.6)
                       !                        print*, 'no root at ', icount, x1,x2,geometry%z(i), geometry%z(j), zdi,zdj,l,xc, & 
                       !                        geometry%erosion_rate(i), geometry%erosion_rate(j)
                       !                        if(xd.le.0.0.or.xd.ge.l)xd=l/2.
                       !                        go to 333
                       !                  endif
                       
                    enddo
333                 continue
                    if(dflag.eq.1)then
                       icountdif=icountdif+1
                    else
                       icountslope=icountslope+1
                    endif
                    !            print*, 'icount ',icount
                    geometry%surface_share(k,i)=xd/l*2.d0-1.d0
                    ! Assign the two neigboring nodes to each divide
                    geometry%nndivnode(ndivide,1) = i
                    geometry%nndivnode(ndivide,2) = j
                    ! calculate the location of the divide
                    local_slope = (geometry%y(j)-geometry%y(i))/(geometry%x(j)-geometry%x(i))
                    if (geometry%x(i).lt.geometry%x(j))then
                       geometry%xdiv(ndivide) = geometry%x(i) + xd/sqrt(1.0+local_slope**2)
                    else
                       geometry%xdiv(ndivide) = geometry%x(i) - xd/sqrt(1.0+local_slope**2) 
                    endif
                    geometry%ydiv(ndivide) = local_slope*(geometry%xdiv(ndivide)-geometry%x(i))+ geometry%y(i)
                    geometry%zdiv(ndivide) = max(zdi,zdj) 
                    if (geometry%zdiv(ndivide).lt.geometry%z(i).or.&
                         &geometry%zdiv(ndivide).lt.geometry%z(j)) then
                       print*, 'divide', ndivide, 'is lower than neighboring verteces'
                       print*, geometry%zdiv(ndivide), geometry%z(i), geometry%z(j)
                    endif

                 else !if used to be a channel in this time step by the upper node was captured
                    print*, 'former channel divide with'
                    print*, 'i', geometry%z(i), orig_z(i)
                    print*, 'j', geometry%z(j), orig_z(j)
                    alpha = atan((orig_z(j) - orig_z(i))/l)
                    theta = atan(params%tanthetac)
                    A = orig_z(j) - geometry%z(j)
                    B = A*l/(orig_z(j) - orig_z(i))
                    D = B*sin(alpha)/sin(pi-alpha-theta)
                    DD = sqrt(A**2+B**2)
                    xd = l-DD+D
                    delz = D*sin(alpha)
                    zd = geometry%z(j) + delz

                    if(dflag.eq.1)then
                       icountdif=icountdif+1
                    else
                       icountslope=icountslope+1
                    endif
                    !            print*, 'icount ',icount
                    geometry%surface_share(k,i)=xd/l*2.d0-1.d0
                    ! Assign the two neigboring nodes to each divide
                    geometry%nndivnode(ndivide,1) = i
                    geometry%nndivnode(ndivide,2) = j
                    ! calculate the location of the divide
                    local_slope = (geometry%y(j)-geometry%y(i))/(geometry%x(j)-geometry%x(i))
                    if (geometry%x(i).lt.geometry%x(j))then
                       geometry%xdiv(ndivide) = geometry%x(i) + xd/sqrt(1.0+local_slope**2)
                    else
                       geometry%xdiv(ndivide) = geometry%x(i) - xd/sqrt(1.0+local_slope**2) 
                    endif
                    geometry%ydiv(ndivide) = local_slope*(geometry%xdiv(ndivide)-geometry%x(i))+ geometry%y(i)
                    geometry%zdiv(ndivide) = zd 
                    if (geometry%zdiv(ndivide).lt.geometry%z(i).or.&
                         &geometry%zdiv(ndivide).lt.geometry%z(j)) then
                       print*, 'divide', ndivide, 'is lower than neighboring verteces'
                       print*, geometry%zdiv(ndivide), geometry%z(i), geometry%z(j)
                    endif
                    print*, 'artificial divide at: x, y, z:'
                    print*, geometry%xdiv(ndivide),geometry%ydiv(ndivide),geometry%zdiv(ndivide)
                    print*, 'i coordinates x y z:', i
                    print*, geometry%x(i),geometry%y(i), geometry%z(i) 
                    print*, 'j coordinates x y z:', j
                    print*, geometry%x(j),geometry%y(j), geometry%z(j) 

                 endif
                 if (geometry%zdiv(ndivide).gt. 10000) then
                    print*, '**************',geometry%zdiv(ndivide)
                 endif
                    
                 ! define the height of the centres of Delaunay circumcircles to be 
                 ! the average height of the divides surrounding it	    
                 !           first find the two triangles that share this divide	
                 it = 0
                 nfound = .TRUE.
                 do while (nfound)
                    it = it+1 	 
                    if (delaunay%icon(1,it).eq.i.or.&
                         &delaunay%icon(2,it).eq.i.or.&
                         &delaunay%icon(3,it).eq.i) then
                       if (delaunay%icon(1,it).eq.j.or.&
                            &delaunay%icon(2,it).eq.j.or.&
                            &delaunay%icon(3,it).eq.j) then
                          geometry%nndivtri(ndivide,1) = it
                          delaunay%centers(3,it) = delaunay%centers(3,it) + geometry%zdiv(ndivide)
                          delaunay%numdivides(1,it) = delaunay%numdivides(1,it)+1
                          delaunay%numdivides(delaunay%numdivides(1,it) +1,it) = ndivide
                          ncount = 0
                          ! next loop we be performed max 3 times 
                          do while (nfound)
                             ncount = ncount+1
                             if (ncount.gt.3) then
                                print*, 'cant find neighboring triangle'
                                print*, 'domain dimensions:', minval(geometry%x),maxval(geometry%x),minval(geometry%y),maxval(geometry%y)
                                print*,'Node geometry (',i,')'
                                print*,'xi:',geometry%x(i)
                                print*,'yi:',geometry%y(i)   
                                print*,'zi:',geometry%z(i)    
                                print*,'fixi',geometry%fix(i)
                                print*,'Node geometry (',j,')'
                                print*,'xj:',geometry%x(j)
                                print*,'yj:',geometry%y(j)
                                print*,'zj:',geometry%z(j)
                                print*,'fixj',geometry%fix(j)
                                print*,'Divide',ndivide,'geometry'
                                print*,'xj:',geometry%xdiv(ndivide)
                                print*,'yj:',geometry%ydiv(ndivide)
                                print*,'zj:',geometry%zdiv(ndivide)
                                print*,'Triangle',it,'geometry'
                                print*,delaunay%icon(1,it),delaunay%icon(2,it),delaunay%icon(3,it)
                                print*, 'neighbors geometry'
                                print*, 'First neighbor is:', delaunay%neighbours(1,it)
                                print*,delaunay%icon(1,delaunay%neighbours(1,it)),delaunay%icon(2,delaunay%neighbours(1,it)),delaunay%icon(3,delaunay%neighbours(1,it))
                                print*, 'Second neighbor is:', delaunay%neighbours(2,it)
                                print*,delaunay%icon(1,delaunay%neighbours(2,it)),delaunay%icon(2,delaunay%neighbours(2,it)),delaunay%icon(3,delaunay%neighbours(2,it))
                                print*, 'Third neighbor is:', delaunay%neighbours(2,it)
                                print*,delaunay%icon(1,delaunay%neighbours(3,it)),delaunay%icon(2,delaunay%neighbours(3,it)),delaunay%icon(3,delaunay%neighbours(3,it))                              
                                nfound = .FALSE.
                             endif
                             ntri = 	delaunay%neighbours(ncount,it)
                             if (ntri.gt.0 .and. ntri .lt. delaunay%ntriangles+1)then
                                if (delaunay%icon(1,ntri).eq.i.or.&
                                     &delaunay%icon(2,ntri).eq.i.or.&
                                     &delaunay%icon(3,ntri).eq.i) then
                                   if (delaunay%icon(1,ntri).eq.j.or.&
                                        &delaunay%icon(2,ntri).eq.j.or.&
                                        &delaunay%icon(3,ntri).eq.j) then
                                      geometry%nndivtri(ndivide,2) = ntri
                                      delaunay%centers(3,ntri) = &
                                           &delaunay%centers(3,ntri) + geometry%zdiv(ndivide)
                                      delaunay%numdivides(1,ntri) = delaunay%numdivides(1,ntri)+1
                                      delaunay%numdivides(delaunay%numdivides(1,ntri) +1,ntri) = ndivide
                                      nfound = .FALSE.
                                   endif
                                endif
                             endif
                          enddo
                       endif
                    endif
                 enddo

              endif
           endif
        enddo
     enddo
     
     if ((params%istep/params%freq)*params%freq.eq.params%istep) then
        !    print*, ' number of diffusive divides:', icountdif
        !    print*, ' number of slope divides:', icountslope
        !    print*, ' Average Distance between nodes: ', avedivide/dble(ndivide)
     endif
     
     
     
     !!            ltest=(geometry%z(i)-geometry%z(j)+2*params%xc*params%tanthetac)/params%tanthetac
     !!            if(l.lt.ltest)then  ! criterion for supershort divide 
!!!  supershort divide
     !!              xdi=(geometry%z(j)-geometry%z(i)+l*params%tanthetac)/(2.*params%tanthetac)
     !!              geometry%surface_share(k,i)=xdi/l*2.d0-1.d0
     !!              if(geometry%surface_share(k,i).lt.-1.0.or.geometry%surface_share(k,i).gt.1.0)print*, 'warning supershort divide error',geometry%z(i), &
     !!                    geometry%z(j),xdi,l
     !!            else
!!! criterion for short divide
     !!             fact1=(geometry%erosion_rate(i)/geometry%k(i)/ &
     !!              (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
     !!             fact2=log(params%xc/(l-params%xc*params%tanthetac))
     !!             if(l.gt.2.0*params%xc)then
     !!              ztprime=geometry%z(i)+params%xc*params%tanthetac-fact1*fact2 ! ztprime is test elevation for short divide calc
     !!             else
     !!               ztprime=(l-params%xc)*params%tanthetac
     !!             endif
     !!             if (ztprime.lt.(geometry%z(j)+params%xc*params%tanthetac)) then ! second case : small divide
     !!               if (params%small_divide) then
     !!                nsmall_divide=nsmall_divide+1 ! keeps track how many small divide are computed
     !!                xdi=l-params%xc ! first guess for xdi
     !!                iter=0
     !!                fx=l
     !!                dfdx=1.d0
     !!                do while (abs(fx/dfdx).gt.l/1.d9) ! Newton-Raphson iterative algorithm to solve simultaneous nonlinear equations
     !!                  iter=iter+1
     !!                  if (iter.gt.100) then ! in case of poor convergence output to the screen
     !!                   print*,'l',l
     !!                   print*,'z',geometry%z(i),geometry%z(j)
     !!                   print*,'erosionrate',geometry%erosion_rate(i),geometry%erosion_rate(j)
     !!                   print*,'precipitation',geometry%precipitation(i),geometry%precipitation(j)
     !!                   print*,'k',geometry%k(i),geometry%k(j)
     !!                   print*,'xdi,fx,dfdx',xdi,fx,dfdx
     !!                   stop 'no convergence in finding small divide'
     !!                  endif
     !!                 xdj=l-xdi
     !!                 fact2=log((params%xc/xdi))
     !!                 fact3=params%tanthetac*(l-params%xc-xdi)
     !!                 fx=geometry%z(j)-geometry%z(i)+fact3+fact1*fact2
     !!                 dfdx=-fact1/xdi-params%tanthetac
     !!                 xdi=xdi-fx/dfdx
     !!                enddo
!!!               if(xdi.lt.0)xdi=0.
     !!               if(xdi.gt.l)xdi=l
     !!                geometry%surface_share(k,i)=xdi/l*2.d0-1.d0
     !!                zdi=geometry%z(i)+params%xc*params%tanthetac & ! zdi is divide height
     !!                   -fact1*log(params%xc/xdi)
     !!                 if(geometry%surface_share(k,i).lt.-1.0.or.geometry%surface_share(k,i).gt.1.0)then
     !!                  print*, 'warning short divide error', geometry%z(i), &
     !!                    geometry%z(j),xdi,l,zdi
     !!                    fact2=log(params%xc/l)
     !!                    zt=geometry%z(i)+params%xc*params%tanthetac-fact1*fact2 !
     !!                    print*, 'zt = ',zt
     !!                 endif
     !!               endif
     !!             else   ! normal divide
     !!               if (params%divide) then
     !!                ndivide=ndivide+1 ! keeps track of how many divides are computed
     !!                xdi=l/2. ! first guess for xdi
     !!                iter=0
     !!                fx=l
     !!                dfdx=1.d0
     !!                do while (abs(fx/dfdx).gt.l/1.d9) ! Newton-Raphson iterative algorithm to solve simultaneous nonlinear equations
     !!                  iter=iter+1
     !!                  if (iter.gt.100) then ! in case of poor convergence output to the screen
     !!                   print*,'l',l
     !!                   print*,'z',geometry%z(i),geometry%z(j)
     !!                   print*,'erosionrate',geometry%erosion_rate(i),geometry%erosion_rate(j)
     !!                   print*,'precipitation',geometry%precipitation(i),geometry%precipitation(j)
     !!                   print*,'k',geometry%k(i),geometry%k(j)
     !!                   print*,'xdi,fx,dfdx',xdi,fx,dfdx
     !!                   stop 'no convergence in finding normal divide'
     !!                  endif
     !!                  xdj=l-xdi
     !!                  fact2=log(params%xc/((xdi)))
     !!                  fact3=(geometry%erosion_rate(j)/geometry%k(j)/ &
     !!                      (geometry%precipitation(j)*params%ka)**params%m)**(1.d0/params%n)
     !!                  fact4=log(params%xc/((xdj)))
     !!                  fx=geometry%z(j)-geometry%z(i)+(fact3*fact4-fact1*fact2)
     !!                  dfdx=-(fact1/xdi+fact3/xdj)
     !!                  xdi=xdi-fx/dfdx
     !!                enddo
     !!                geometry%surface_share(k,i)=xdi/l*2.d0-1.d0
     !!                if(geometry%surface_share(k,i).lt.-1.0.or.geometry%surface_share(k,i).gt.1.0)print*, 'warning normal divide error', geometry%z(i), &
     !!                    geometry%z(j), xdi, l
     !!                zdi=geometry%z(i)+params%xc*params%tanthetac-fact1*log(params%xc/xdi) ! zdi is divide height
     !!               endif
     !!             endif !small divide
     !!            endif !supershort
     !!
     !!
     !!          endif
     !!        endif
     !!      enddo
!!!    endif
     !!  enddo
     !  replace above for diffusion change
     
     
  else ! if mh/n not equal 1
     do iii=1, maxpass
        ncapture=0
        do ii=1,geometry%nnode ! loop over the nodes
           i=ishuffle(ii)
           do k=geometry%nb(i),1,-1 ! loop over all neighbours
              j=geometry%nn(k,i)
              if(network%receiver(j).ne.i.and.network%receiver(i).ne.j) then ! only for connections that are not part of the drainage network
                 !if (geometry%z(i).lt.geometry%z(j)) then ! check that i is lower than j
                 if (geometry%z(i).lt.geometry%z(j).and.orig_receiver(j).ne.i) then ! check that i is lower than j and than j didnt use to flow to i in this time step
                    l=sqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length betwenn i and j
                    if(l.gt.params%xc) then
                       fact1=(geometry%erosion_rate(i)/geometry%k(i)/ &
                            (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                       fact2=(1./params%hmn)*((l)**(params%hmn)-(params%xc)**(params%hmn))
                       zt=params%xc*params%tanthetac ! zt is channel head elevation
                       if(params%diffusion)then
                          !   erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)
                          !   ztd=erodrate*params%xc**2/params%diffusivity
                          !   if(ztd.lt.zt)then
                          !      zt=ztd
                          !   endif
                       endif
                       zt=zt+geometry%z(i)+fact1*fact2 ! zt is test elevation
                    else
                       zt=l*params%tanthetac
                       if(params%diffusion)then
                          !  erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)
                          !  ztd=erodrate*l**2/params%diffusivity
                          !  if(ztd.lt.zt)then
                          !     zt=ztd
                          !  endif
                       endif
                       zt=zt+geometry%z(i)
                    endif
                    if(zt.lt.geometry%z(i))then
                       print*, 'ERROR - zt capturing a lower point'
                       print*,'l', l
                       print*,'erosion rate', geometry%erosion_rate(i)
                       print*,'fact1, fact2', fact1, fact2
                       print*,'zt', zt
                    endif
                    if (zt.lt.geometry%z(j)) then ! first case : capture
                       if (params%capture) then
                          if (orig_receiver(j).eq.i) then
                             print*, i, 'is capturing its former donor', j, 'FORBIDDEN'
                          endif
                          ncapture=ncapture+1
                          ncapglobe=ncapglobe+1
                          recnode(ncapture)=i
                          capnode(ncapture)=j
                          capelev(ncapture)=zt
                       endif
                    endif
                 endif
              endif
           enddo
        enddo
        if(ncapture.eq.0)go to 1493
        ! check for two nodes capturing the same node and select the lower capture elevation
        if(ncapture.gt.1) then
           do icap=1,ncapture
              do jcap=icap+1,ncapture
                 if(capnode(icap).eq.capnode(jcap))then
                    if(capelev(jcap).lt.capelev(icap))then 
                       capelev(icap)=geometry%z(capnode(icap))
                    else
                       recnode(jcap)=recnode(icap)
                       capnode(jcap)=capnode(icap)
                       capelev(jcap)=capelev(icap)
                       capelev(icap)=geometry%z(capnode(icap))
                    endif
                 endif
              enddo
           enddo
        endif
        do icap=1,ncapture
           geometry%erosion_rate(capnode(icap))=(geometry%z(capnode(icap))-capelev(icap))/params%deltat
           !      if(geometry%erosion_rate(capnode(icap)).lt.0)print*, 'negative erosion rate'
           geometry%z(capnode(icap))=capelev(icap)
           network%receiver(capnode(icap))=recnode(icap)
        enddo
     enddo
     print*, 'warning max passes exceeded'
1493 continue
     !         ******************  start loop for divide heights      **********************************
     ! May 2011 - modification: When the two neighbors that now share a divide used to be a channel until this 
     !            time step, the divide finding routine is wrong because:
     !            1) the contributing areas are not corrected 
     !            2) We prevent capture by the original reciever (double capture) and for such cases there is no solution 
     !               to the divide finding problem.
     !             In such a case we artificialy set the divide in the intersection of the previous channel with a line
     !             sloping theta leaving the node that was captured in this time step.
     !
     !
     ! set a tolerance for the ridge height convergence
     ztoler=.1
     icountdif=0
     icountslope=0
     do i=1,geometry%nnode ! loop over the nodes
        do k=geometry%nb(i),1,-1 ! loop over all neighbours
           j=geometry%nn(k,i)
           if(network%receiver(j).ne.i.and.network%receiver(i).ne.j) then ! only for connections that are not part of the drainage network
              if (geometry%z(i).lt.geometry%z(j)) then ! check that i is lower than j
                 l=sqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length of segment betwenn i and j
                 ndivide=ndivide+1
                 avedivide=avedivide+l

                 if (orig_receiver(j).ne.i) then !if didn't use to be a channel in this time step
                    ! give initial guess for divide and endpoints
                    x1=0.
                    x2=l
                    zdi=0
                    zdj=2*ztoler
                    icount=0
                    ! find ridge heights at opposing node positions
                    !  find ridge height from node i
                    xc=min(params%xc,l)
                    erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)
                    if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/erodrate))then 
                       zc=erodrate*xc**2/params%diffusivity
                    else
                       zc=xc*params%tanthetac
                    endif
                    if(l.le.params%xc)then
                       zi2=geometry%z(i)+zc
                    else
                       fact2=(1./params%hmn)*(l**(params%hmn)-(params%xc)**(params%hmn))
                       fact3=(geometry%erosion_rate(i)/geometry%k(i)/ &
                            (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                       zi2=geometry%z(i)+zc+fact2*fact3
                    endif
                    zi1=geometry%z(i)
                    !  find ridge height from node j
                    xc=min(params%xc,(l))
                    erodrate=max(geometry%erosion_rate(j),params%min_erosion_rate)
                    if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/erodrate))then 
                       zc=erodrate*xc**2/params%diffusivity
                    else
                       zc=xc*params%tanthetac
                    endif
                    if((l).le.params%xc)then
                       zj1=geometry%z(j)+zc
                    else
                       fact2=(1./params%hmn)*(l**(params%hmn)-(params%xc)**(params%hmn))
                       fact3=(geometry%erosion_rate(j)/geometry%k(j)/ &
                            (geometry%precipitation(j)*params%ka)**params%m)**(1.d0/params%n)
                       zj1=geometry%z(j)+zc+fact2*fact3
                    endif
                    zj2=geometry%z(j)
                    if((zi1-zj1).ge.0.or.(zi2-zj2).le.0)then
                       print*, 'no root in initial conditions'
                       print*, icount, zi1,zj1,zi2,zj2
                    endif
                    ! iterate using false point method until divide height is found within specified tolerance
                    do while (dabs(zdi-zdj).gt.ztoler)

                       icount=icount+1
                       dflag=0
                       !  find root of function zdi-zdj
                       !  find new guess at divide position from false point method
                       !  unless ftc after 100 iterations, then revert to bisector method
                       if(icount.lt.100)then
                          xd=(x1*(zi2-zj2)-x2*(zi1-zj1))/(zi2-zj2-zi1+zj1)
                       else
                          xd=x1+.5*(x2-x1)
                       endif
                       if((zi1-zj1).ge.0.or.(zi2-zj2).le.0)then
                          print*, 'error no root found for divide problem'
                          print*, icount, zi1,zj1,zi2,zj2
                          xd=l/2.
                          go to 334
                       endif

                       !  find ridge height from node i
                       xc=min(params%xc,xd)
                       erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)
                       if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/erodrate))then 
                          zc=erodrate*xc**2/params%diffusivity
                          dflag=1
                       else
                          zc=xc*params%tanthetac
                       endif
                       if(xd.le.params%xc)then
                          zdi=geometry%z(i)+zc
                       else
                          !fact2=log(params%xc/((xd)))
                          fact2 = (1./params%hmn)*(xd**(params%hmn)-(params%xc)**(params%hmn))
                          fact3=(geometry%erosion_rate(i)/geometry%k(i)/ &
                               (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                          zdi=geometry%z(i)+zc+fact2*fact3
                       endif
                       !  find ridge height from node j
                       xc=min(params%xc,(l-xd))
                       erodrate=max(geometry%erosion_rate(j),params%min_erosion_rate)
                       if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/erodrate))then 
                          zc=erodrate*xc**2/params%diffusivity
                          dflag=1
                       else
                          zc=xc*params%tanthetac
                       endif
                       if((l-xd).le.params%xc)then
                          zdj=geometry%z(j)+zc
                       else
                          !fact2=log(params%xc/((l-xd)))
                          fact2 = (1./params%hmn)*((l-xd)**(params%hmn)-(params%xc)**(params%hmn))
                          fact3=(geometry%erosion_rate(j)/geometry%k(j)/ &
                               (geometry%precipitation(j)*params%ka)**params%m)**(1.d0/params%n)
                          zdj=geometry%z(j)+zc+fact2*fact3
                       endif
                       !  reset bounds of interval
                       if(zdj.gt.zdi)then
                          x1=xd
                          zi1=zdi
                          zj1=zdj
                       else
                          x2=xd
                          zi2=zdi
                          zj2=zdj
                       endif
                    enddo

                    if (max(zdi,zdj).gt.20000) then
                       print*, 'node i in height',geometry%z(i) 
                       print*, 'node j in height',geometry%z(j)
                       print*, 'Erosion number is',fact3 
                       print*, '1-hm/n is',params%hmn
                       print*, 'l is', l, 'xd is', xd 
                       print*, 'xc is', params%xc, 'tantheta is', params%tanthetac
                       print*, 'divide from i is', zdi, 'divide from j is', zdj
                       print*, 'erosion rate', geometry%erosion_rate(j)
                       print*, 'k', geometry%k(j)
                       print*, 'precipitation', geometry%precipitation(j)
                    endif

334                 continue
                    if(dflag.eq.1)then
                       icountdif=icountdif+1
                    else
                       icountslope=icountslope+1
                    endif
                    geometry%surface_share(k,i)=xd/l*2.d0-1.d0
                    ! Assign the two neigboring nodes to each divide
                    geometry%nndivnode(ndivide,1) = i
                    geometry%nndivnode(ndivide,2) = j
                    ! calculate the location of the divide
                    if (geometry%x(j).eq.geometry%x(i)) then
                       geometry%xdiv(ndivide) = geometry%x(i)
                       geometry%ydiv(ndivide) = min(geometry%y(i),geometry%y(j)) + xd
                    else
                       local_slope = (geometry%y(j)-geometry%y(i))/(geometry%x(j)-geometry%x(i))
                       if (geometry%x(i).lt.geometry%x(j))then
                          geometry%xdiv(ndivide) = geometry%x(i) + xd/sqrt(1.0+local_slope**2)
                       else
                          geometry%xdiv(ndivide) = geometry%x(i) - xd/sqrt(1.0+local_slope**2) 
                       endif
                       geometry%ydiv(ndivide) = local_slope*(geometry%xdiv(ndivide)-geometry%x(i))+ geometry%y(i)
                    endif
                    geometry%zdiv(ndivide) = max(zdi,zdj) 
                    if (geometry%zdiv(ndivide).lt.geometry%z(i).or.&
                         &geometry%zdiv(ndivide).lt.geometry%z(j)) then
                       print*, 'divide', ndivide, 'is lower than neighboring verteces'
                       print*, geometry%zdiv(ndivide), geometry%z(i), geometry%z(j)
                    endif

                 else !if used to be a channel in this time step (but upper node was captured)


                 endif





                 ! define the height of the centres of Delaunay circumcircles to be 
                 ! the average height of the divides surrounding it	    
                 !           first find the two triangles that share this divide	
                 it = 0
                 nfound = .TRUE.
                 do while (nfound)
                    it = it+1 	 
                    if (delaunay%icon(1,it).eq.i.or.&
                         &delaunay%icon(2,it).eq.i.or.&
                         &delaunay%icon(3,it).eq.i) then
                       if (delaunay%icon(1,it).eq.j.or.&
                            &delaunay%icon(2,it).eq.j.or.&
                            &delaunay%icon(3,it).eq.j) then
                          geometry%nndivtri(ndivide,1) = it
                          delaunay%centers(3,it) = delaunay%centers(3,it) + geometry%zdiv(ndivide)
                          delaunay%numdivides(1,it) = delaunay%numdivides(1,it)+1
                          delaunay%numdivides(delaunay%numdivides(1,it) +1,it) = ndivide
                          ncount = 0
                          ! next loop we be performed max 3 times 
                          do while (nfound)
                             ncount = ncount+1
                             ntri = delaunay%neighbours(ncount,it)
                             if (ntri.gt.0 .and. ntri .lt. delaunay%ntriangles+1)then
                                if (delaunay%icon(1,ntri).eq.i.or.&
                                     &delaunay%icon(2,ntri).eq.i.or.&
                                     &delaunay%icon(3,ntri).eq.i) then
                                   if (delaunay%icon(1,ntri).eq.j.or.&
                                        &delaunay%icon(2,ntri).eq.j.or.&
                                        &delaunay%icon(3,ntri).eq.j) then
                                      geometry%nndivtri(ndivide,2) = ntri
                                      delaunay%centers(3,ntri) = &
                                           &delaunay%centers(3,ntri) + geometry%zdiv(ndivide)
                                      index = delaunay%numdivides(1,ntri) +1
                                      delaunay%numdivides(1,ntri) = delaunay%numdivides(1,ntri)+1
                                      delaunay%numdivides(delaunay%numdivides(1,ntri) +1,ntri) = ndivide
                                      nfound = .FALSE.
                                   endif
                                endif
                             endif
                          enddo
                       endif
                    endif
                 enddo
              endif
           endif
        enddo
     enddo




     ! case where hm/n is not =1
     ! this section has been used/tested (FH)
     ! needs to be restructured like above (sdw)
!!!!  do i=1,geometry%nnode ! loop over the nodes
!!!!    do k=1,geometry%nb(i) ! loop over all neighbours
!!!!    j=geometry%nn(k,i)
!!!!      if (network%receiver(i).ne.j .and. network%receiver(j).ne.i) then ! only for connections that are not part of the drainage network
!!!!        if (geometry%z(i).lt.geometry%z(j)) then ! check that i is lower than j
!!!!        l=sqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length of segment betwenn i and j
!!!!        fact1=(geometry%erosion_rate(i)/geometry%k(i)/ &
!!!!              (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
!!!!        fact2=l**params%hmn*(1.d0-(params%xc/l)**params%hmn)/params%hmn
!!!!        zt=geometry%z(i)+params%xc*params%tanthetac+fact1*fact2 ! zt is test elevation
!!!!          if (zt.lt.geometry%z(j)) then ! first case : capture
!!!!            if (params%capture) then
!!!!            ncapture=ncapture+1
!!!!            network%receiver(j)=i
!!!!            geometry%z(j)=zt
!!!!            endif
!!!!          else
!!!!            if (params%small_divide) then
!!!!            if (zt.lt.geometry%z(j)+params%xc*params%tanthetac) then ! second case : small divide
!!!!            nsmall_divide=nsmall_divide+1 ! keeps track how many small divide are computed
!!!!            xdi=l-(zt-geometry%z(j))/params%tanthetac ! xdi is distance to divide
!!!!            geometry%surface_share(k,i)=xdi/l*2.d0-1.d0
!!!!            zdi=geometry%z(i)+params%xc*params%tanthetac & ! zdi is divide height
!!!!               +fact1*xdi**params%hmn*(1.d0-(params%xc/xdi)**params%hmn)/params%hmn
!!!!            endif
!!!!          else ! third case : normal divide
!!!!            if (params%divide) then
!!!!            ndivide=ndivide+1 ! keeps track how many divide are computed
!!!!            xdi=params%xc ! first guess for xdi
!!!!            iter=0
!!!!            fx=l
!!!!            dfdx=1.d0
!!!!              do while (abs(fx/dfdx).gt.l/1.d6) ! Newton-Raphson iterative algorithm to solve simultaneous nonlinear equations
!!!!              iter=iter+1
!!!!                if (iter.gt.100) then ! in case of poor convergence output to the screen
!!!!                print*,'Nodes',i,j
!!!!                print*,'l, xdi',l,xdi
!!!!                print*,'z',geometry%z(i),geometry%z(j)
!!!!                print*,'erosionrate',geometry%erosion_rate(i),geometry%erosion_rate(j)
!!!!                print*,'precipitation',geometry%precipitation(i),geometry%precipitation(j)
!!!!                print*,'k',geometry%k(i),geometry%k(j)
!!!!                print*,fact1,fact2,fact3,fact4
!!!!                print*,'xdi,fx,dfdx',xdi,fx,dfdx
!!!!                stop 'no convergence in finding divide'
!!!!                endif
!!!!              xdj=l-xdi
!!!!              fact2=xdi**params%hmn*(1.d0-(params%xc/xdi)**params%hmn)/params%hmn
!!!!              fact3=(geometry%erosion_rate(j)/geometry%k(j)/ &
!!!!                    (geometry%precipitation(j)*params%ka)**params%m)**(1.d0/params%n)
!!!!              fact4=xdj**params%hmn*(1.d0-(params%xc/xdj)**params%hmn)/params%hmn
!!!!              fx=geometry%z(j)-geometry%z(i)+(fact3*fact4-fact1*fact2)
!!!!              dfdx=-(fact1*xdi**(params%hmn-1.d0)+fact3*xdj**(params%hmn-1.d0))
!!!!              xdi=xdi-fx/dfdx
!!!!              enddo
!!!!            geometry%surface_share(k,i)=xdi/l*2.d0-1.d0
!!!!            zdi=geometry%z(i)+params%xc*params%tanthetac & ! zdi is divide height
!!!!               +fact1*xdi**params%hmn*(1.d0-(params%xc/xdi)**params%hmn)/params%hmn
!!!!            endif
!!!!          endif
!!!!        endif
!!!!      endif
!!!!      endif
!!!!    enddo
!!!!  enddo
  endif
  
  ! new piece - sdw 
  ! set contributing area down channel to .2 of distance towards receiver
  !  do i=1,geometry%nnode ! loop over the nodes
  !    do k=geometry%nb(i),1,-1 ! loop over all neighbours
  !    j=geometry%nn(k,i)
  !          if(network%receiver(i).eq.j) geometry%surface_share(k,i)=.2d0*2.d0-1.d0
  !    enddo
  !  enddo
  ! print result to screen

  do i = 1,delaunay%ntriangles
     if (delaunay%numdivides(1,i).ne.0) then
        delaunay%centers(3,i) = delaunay%centers(3,i)/delaunay%numdivides(1,i)
     endif
  enddo

  
  if (ndivide.ne.0) then
     !print*,ndivide,'divides'
     !print*,nsmall_divide,'small divides'
     geometry%ndivide=ndivide
     !print*,'surface_share',minval(geometry%surface_share),maxval(geometry%surface_share)
  endif
  
  ! print result to screen and adapt network to reflec stream captures
  
  if (params%capture.and.ncapglobe.gt.0) then
     print*,ncapglobe,'captures'
     geometry%ncapture=geometry%ncapture+ncapglobe
     network%donors=0
     network%ndon=0
     do i=1,network%nnode
        k=network%receiver(i)
        if (k.ne.0) then
           network%ndon(k)=network%ndon(k)+1
           network%donors(network%ndon(k),k)=i
        endif
     enddo
  endif
  
  !!check for uphill rivers
  do i=1,geometry%nnode
     if (network%receiver(i).ne.0.and.geometry%fix(i).ne.1) then
        j=network%receiver(i)
        if (geometry%z(j).gt.geometry%z(i)) then
           print*,'flowing uphill at end of captures',i,j
        endif
     endif
  enddo

  deallocate(orig_receiver)
  deallocate(orig_z)
  

  call time_out ('captures_and_divides')
  
  return
  
end subroutine captures_and_divides
