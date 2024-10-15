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

  integer,dimension(:),allocatable::orig_receiver
  double precision, dimension(:),allocatable::orig_z
  integer counter

  integer cpn, num_samples, num_full_bin,sample_in_last_bin,h 
  double precision local_erosioni, local_erosionj, tau,ave_erosion 

  call time_in ('captures_and_divides')

  maxpass=500
  ncapture=0
  ncapglobe=0
  ndivide=0
  nsmall_divide=0
  geometry%surface_share=0.d0
  avedivide=0.0d0

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
           geometry%erosion_rate(i)=0.0d0 
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
                 if (geometry%z(i).lt.geometry%z(j)) then! check that i is lower than j 
                    l=dsqrt((geometry%x(i)-geometry%x(j))**2.d0+(geometry%y(i)-geometry%y(j))**2.d0) ! l is length betwenn i and j
                    if(l.gt.params%xc)then
                       if (params%transient_divide) then
                          if (params%h*params%m.eq.1) then
                             tau=log(l/params%xc)/ &
                                  (geometry%k(i)*(geometry%precipitation(i)*params%ka)**(params%m))
                          else
                             tau= (l**(1.d0-params%h*params%m) - params%xc**(1.d0-params%h*params%m)) /&
                                  (geometry%k(i)*(geometry%precipitation(i)*params%ka)**(params%m))/(1.d0-params%h*params%m)
                          endif
                          num_samples=INT(min(tau,params%time)/params%deltat)
                          if (num_samples.eq.0) then
                             fact1=(geometry%erosion_rate(i)/geometry%k(i)/ &
                                  (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                          else
                             if ((params%num_bins-1)*params%sample_per_bin.lt.num_samples) then
                                num_full_bin=params%num_bins
                                sample_in_last_bin=num_samples-(params%num_bins-1)*params%sample_per_bin
                             else
                                num_full_bin = INT(num_samples/params%sample_per_bin)+1
                                sample_in_last_bin=num_samples-(num_full_bin-1)*params%sample_per_bin
                             endif
                             ave_erosion=0.d0
                             if (num_full_bin.eq.1) then
                                ave_erosion=geometry%erosion_rate_history(i,1)
                             else
                                ave_erosion=0.d0
                                do h=1,num_full_bin-1
                                   ave_erosion=ave_erosion +geometry%erosion_rate_history(i,h)*params%sample_per_bin
                                enddo
                                ave_erosion=ave_erosion +geometry%erosion_rate_history(i,num_full_bin)*sample_in_last_bin
                                ave_erosion=ave_erosion/((num_full_bin-1)*params%sample_per_bin + sample_in_last_bin)
                             endif
                             fact1=(ave_erosion/geometry%k(i)/ &
                                  (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                          endif
                       else
                          fact1=(geometry%erosion_rate(i)/geometry%k(i)/ &
                          (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                       endif
                       fact2=log(params%xc/l)
                       zt=params%xc*params%tanthetac ! zt is channel head elevation
                       if(params%diffusion)then
                          !erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)
                          if (params%transient_divide) then
                             ztd=ave_erosion*params%xc**2.d0/params%diffusivity
                          else
                             ztd=geometry%erosion_rate(i)*params%xc**2.d0/params%diffusivity
                          endif
                          if(ztd.lt.zt)then
                             zt=max(ztd,params%xc*params%min_tan_head_slope)
                          endif
                       endif
                       zt=zt+geometry%z(i)-fact1*fact2 ! zt is test elevation
                    else
                       zt=l*params%tanthetac
                       if(params%diffusion)then
                          !erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)
                          ztd=geometry%erosion_rate(i)*l**2.d0/params%diffusivity
                          if(ztd.lt.zt)then
                             zt=max(ztd,l*params%min_tan_head_slope)
                          endif
                       endif
                       zt=zt+geometry%z(i)
                    endif
                    if(zt.lt.geometry%z(i))then
                       print*, 'ERROR - zt capturing a lower point'
                       print*,'l', l
                       if (params%transient_divide.and.l.gt.params%xc) then
                          print*,'ave erosion rate', ave_erosion
                       else
                          print*,'erosion rate', geometry%erosion_rate(i)
                       endif
                       print*,'fact1, fact2', fact1, fact2
                       print*,'zt', zt
                    endif
                    if (zt.lt.geometry%z(j)) then ! first case : capture
                       if (params%capture) then
                          if (orig_receiver(j).eq.i) then
                             !print*, i, 'is caturing its former donor', j
                          endif
                          ncapture=ncapture+1
                          ncapglobe=ncapglobe+1
                          recnode(ncapture)=i
                          capnode(ncapture)=j
                          capelev(ncapture)=zt
                          !print*, 'zt=',zt
                          !print*, geometry%z(j)-geometry%z(i)
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
           !*** Old calculations of erosion rates during capture - L.G.***
           !geometry%erosion_rate(capnode(icap))=(geometry%z(capnode(icap))-capelev(icap))/params%deltat
           !geometry%erosion_rate(capnode(icap))=(orig_z(capnode(icap))-capelev(icap))/params%deltat

           !*** New calculations of erosion rates during capture - L.G.***
           ! In the calculation of the erosion rate add the elevation drop fue to capture to previous 
           ! elevation drops both due to othetr captures in this loop and due to fluvial erosion in the former time step


           if (params%transient_divide) then
              cpn=capnode(icap) !cpn is the captured node
              ! need to update the first bin of the captured node with the correct erosion rate
              ! due to capture. First remove the last added erosion rate
              num_full_bin = INT(params%istep/params%sample_per_bin)+1
              if (num_full_bin.eq.1) then
                 sample_in_last_bin=params%istep
                 geometry%erosion_rate_history(cpn,1)= &
                      geometry%erosion_rate_history(cpn,1) - geometry%erosion_rate(cpn)/sample_in_last_bin
              else
                 geometry%erosion_rate_history(cpn,1)= &
                      geometry%erosion_rate_history(cpn,1) - geometry%erosion_rate(cpn)/(params%sample_per_bin)
              endif
           endif


        !print*, 'erosion rate before update', geometry%erosion_rate(capnode(icap))
        geometry%erosion_rate(capnode(icap)) = (geometry%erosion_rate(capnode(icap))*params%deltat +&
             &geometry%z(capnode(icap))-capelev(icap))/params%deltat
        
        !add to the first bin the updated erosion rate
        if (params%transient_divide) then
              cpn=capnode(icap) !cpn is the captured node
              ! need to update the first bin of the captured node with the correct erosion rate
              ! due to capture. First remove the last added erosion rate
              num_full_bin = INT(params%istep/params%sample_per_bin)+1
              if (num_full_bin.eq.1) then
                 sample_in_last_bin=params%istep
                 geometry%erosion_rate_history(cpn,1)= &
                      geometry%erosion_rate_history(cpn,1) + geometry%erosion_rate(cpn)/sample_in_last_bin
              else
                 geometry%erosion_rate_history(cpn,1)= &
                      geometry%erosion_rate_history(cpn,1) + geometry%erosion_rate(cpn)/(params%sample_per_bin)
              endif
           endif


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
  !print*, 'finished with capture'
  ztoler=.1d0
  icountdif=0
  icountslope=0
  do i=1,geometry%nnode ! loop over the nodes
     do k=geometry%nb(i),1,-1 ! loop over all neighbours
        j=geometry%nn(k,i)
        if(network%receiver(j).ne.i.and.network%receiver(i).ne.j) then ! only for connections that are not part of the drainage network
           if (geometry%z(i).lt.geometry%z(j)) then ! check that i is lower than j                
              l=dsqrt((geometry%x(i)-geometry%x(j))**2.d0+(geometry%y(i)-geometry%y(j))**2.d0) ! l is length of segment betwenn i and j
              ndivide=ndivide+1
              avedivide=avedivide+l
              ! start here with new diffusion algorithm
              
              ! give initial guess for divide and endpoints
              x1=0.d0
              x2=l
              zdi=0.d0
              zdj=2.d0*ztoler
              icount=0
              ! find ridge heights at opposing node positions
              !  find ridge height from node i
              xc=min(params%xc,l)
              !erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)
              if (params%transient_divide.and.l.gt.params%xc) then
                 if (params%h*params%m.eq.1) then
                    tau=log(l/params%xc)/ &
                         (geometry%k(i)*(geometry%precipitation(i)*params%ka)**(params%m))                   
                 else
                    tau= (l**(1.d0-params%h*params%m) - params%xc**(1.d0-params%h*params%m)) /&
                         (geometry%k(i)*(geometry%precipitation(i)*params%ka)**(params%m))/(1.d0-params%h*params%m)                  
                 endif
                 !if (tau.lt.params%time) print*, 'tau is smaller'
                 num_samples=INT(min(tau,params%time)/params%deltat)
                 if (num_samples.eq.0) then
                    local_erosioni=geometry%erosion_rate(i)
                 else
                    if ((params%num_bins-1)*params%sample_per_bin.lt.num_samples) then
                       num_full_bin=params%num_bins
                       sample_in_last_bin=num_samples-(params%num_bins-1)*params%sample_per_bin
                    else
                       num_full_bin = INT(num_samples/params%sample_per_bin)+1
                       sample_in_last_bin=num_samples-(num_full_bin-1)*params%sample_per_bin
                    endif
                    ave_erosion=0.d0
                    if (num_full_bin.eq.1) then
                       ave_erosion=geometry%erosion_rate_history(i,1)
                    else
                       ave_erosion=0.d0
                       do h=1,num_full_bin-1
                          ave_erosion=ave_erosion +geometry%erosion_rate_history(i,h)*params%sample_per_bin
                       enddo
                       ave_erosion=ave_erosion +geometry%erosion_rate_history(i,num_full_bin)*sample_in_last_bin
                       ave_erosion=ave_erosion/((num_full_bin-1)*params%sample_per_bin + sample_in_last_bin)
                    endif
                    local_erosioni=ave_erosion
                 endif
                 if (params%h*params%m.eq.1) then
                    tau=log(l/params%xc)/ &
                         (geometry%k(j)*(geometry%precipitation(j)*params%ka)**(params%m))                   
                 else
                    tau= (l**(1.d0-params%h*params%m) - params%xc**(1.d0-params%h*params%m)) /&
                         (geometry%k(j)*(geometry%precipitation(j)*params%ka)**(params%m))/(1.d0-params%h*params%m)                  
                 endif
                 !if (tau.lt.params%time) print*, 'tau is smaller'
                 num_samples=INT(min(tau,params%time)/params%deltat)
                 if (num_samples.eq.0) then
                    local_erosionj=geometry%erosion_rate(j)
                 else
                    if ((params%num_bins-1)*params%sample_per_bin.lt.num_samples) then
                       num_full_bin=params%num_bins
                       sample_in_last_bin=num_samples-(params%num_bins-1)*params%sample_per_bin
                    else
                       num_full_bin = INT(num_samples/params%sample_per_bin)+1
                       sample_in_last_bin=num_samples-(num_full_bin-1)*params%sample_per_bin
                    endif
                    ave_erosion=0.d0
                    if (num_full_bin.eq.1) then
                       ave_erosion=geometry%erosion_rate_history(j,1)
                    else
                       ave_erosion=0.d0
                       do h=1,num_full_bin-1
                          ave_erosion=ave_erosion +geometry%erosion_rate_history(j,h)*params%sample_per_bin
                       enddo
                       ave_erosion=ave_erosion +geometry%erosion_rate_history(j,num_full_bin)*sample_in_last_bin
                       ave_erosion=ave_erosion/((num_full_bin-1)*params%sample_per_bin + sample_in_last_bin)
                    endif
                    local_erosionj=ave_erosion
                 endif
              else
                 local_erosioni=geometry%erosion_rate(i)
                 local_erosionj=geometry%erosion_rate(j)
              endif
              if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/local_erosioni))then 
                 zc=max(local_erosioni*xc**2.d0/params%diffusivity,xc*params%min_tan_head_slope)   
              else
                 zc=xc*params%tanthetac
              endif
              if(l.le.params%xc)then
                 zi2=geometry%z(i)+zc
              else
                 fact2=log(params%xc/((l)))
                 fact3=(local_erosioni/geometry%k(i)/ &
                 (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                 zi2=geometry%z(i)+zc-fact2*fact3
              endif
              zi1=geometry%z(i)
              if (zi2.eq. geometry%z(j)) then
                 print*, 'divide at height of node'
                 print*, i,j
                 print*, geometry%fix(i), geometry%fix(j)
                 print*, network%receiver(i),network%receiver(j)
                 print*, geometry%z(i), geometry%z(j)
                 print*, geometry%erosion_rate(i), geometry%erosion_rate(j)
                 print*, local_erosioni,local_erosionj
                 print*, l,params%xc
                 print*, zc, fact2, fact3
                 print*, params%diffusivity*params%tanthetac/local_erosioni,local_erosioni*xc**2.d0/params%diffusivity,xc*params%min_tan_head_slope
              endif
              !  find ridge height from node j
              xc=min(params%xc,(l))
              !erodrate=max(geometry%erosion_rate(j),params%min_erosion_rate)
              if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/local_erosionj))then 
                 zc=max(local_erosionj*xc**2.d0/params%diffusivity,xc*params%min_tan_head_slope)
              else
                 zc=xc*params%tanthetac
              endif
              if((l).le.params%xc)then
                 zj1=geometry%z(j)+zc
              else
                 fact2=log(params%xc/((l)))
                 fact3=(local_erosionj/geometry%k(j)/ &
                 (geometry%precipitation(j)*params%ka)**params%m)**(1.d0/params%n)
                 zj1=geometry%z(j)+zc-fact2*fact3
              endif
              zj2=geometry%z(j)
              if((zi1-zj1).ge.0.d0.or.(zi2-zj2).le.0.d0)then
                 print*, 'no root in initial conditions'
                 print*, i,j, geometry%z(i), geometry%z(j)
                 print*, icount, zi1,zj1,zi2,zj2
                 print*, local_erosioni,local_erosionj
                 print*, l,zc,fact2*fact3
              endif
              ! iterate using false point method until divide height is found within specified tolerance
              do while (dabs(zdi-zdj).gt.ztoler)
                 if (icount.gt.10000) stop 'icount > 10000'

                 icount=icount+1
                 dflag=0
                 !  find root of function zdi-zdj
                 !  find new guess at divide position from false point method
                 !  unless ftc after 100 iterations, then revert to bisector method
                 if(icount.lt.100)then
                    xd=(x1*(zi2-zj2)-x2*(zi1-zj1))/(zi2-zj2-zi1+zj1)
                 else
                    xd=x1+.5d0*(x2-x1)
                 endif
                 if((zi1-zj1).ge.0.d0.or.(zi2-zj2).le.0.d0)then
                    print*, 'error no root found for divide problem'
                    print*, icount, zi1,zj1,zi2,zj2
                    xd=l/2.d0
                    go to 333
                 endif
                 
                 !  find ridge height from node i
                 xc=min(params%xc,xd)
                 !erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)
                 if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/local_erosioni))then 
                    zc=max(local_erosioni*xc**2.d0/params%diffusivity,xc*params%min_tan_head_slope)
                    dflag=1
                 else
                    zc=xc*params%tanthetac
                 endif
                 if(xd.le.params%xc)then
                    zdi=geometry%z(i)+zc
                 else
                    fact2=log(params%xc/((xd)))
                    fact3=(local_erosioni/geometry%k(i)/ &
                    (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                    zdi=geometry%z(i)+zc-fact2*fact3
                 endif
                 !  find ridge height from node j
                 xc=min(params%xc,(l-xd))
                 !erodrate=max(geometry%erosion_rate(j),params%min_erosion_rate)
                 if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/local_erosionj))then 
                    zc=max(local_erosionj*xc**2.d0/params%diffusivity,xc*params%min_tan_head_slope)
                    dflag=1
                 else
                    zc=xc*params%tanthetac
                 endif
                 if((l-xd).le.params%xc)then
                    zdj=geometry%z(j)+zc
                 else
                    fact2=log(params%xc/((l-xd)))
                    fact3=(local_erosionj/geometry%k(j)/ &
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
333           continue
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
                 geometry%xdiv(ndivide) = geometry%x(i) + xd/dsqrt(1.0d0+local_slope**2.d0)
              else
                 geometry%xdiv(ndivide) = geometry%x(i) - xd/dsqrt(1.0d0+local_slope**2.d0) 
              endif
              geometry%ydiv(ndivide) = local_slope*(geometry%xdiv(ndivide)-geometry%x(i))+ geometry%y(i)
              geometry%zdiv(ndivide) = max(zdi,zdj) 

              if (geometry%zdiv(ndivide).gt.10000.0d0) then
                 print*, 'divide elevation',geometry%zdiv(ndivide)
                 print*, 'neighboring nodes',i,j
                 print*, 'elevation of neighbors', geometry%z(i),geometry%z(j)
                 print*, 'distance between neighbors', l
                 print*, 'erosion rates', geometry%erosion_rate(i),geometry%erosion_rate(j) 
                 print*, 'using for erosion rate',local_erosioni,local_erosionj
                 print*, 'topobuiling', zc,fact2,fact3
              endif



              if (geometry%zdiv(ndivide).lt.geometry%z(i).or.&
                   &geometry%zdiv(ndivide).lt.geometry%z(j)) then
                 print*, 'divide', ndivide, 'is lower than neighboring verteces'
                 print*, geometry%zdiv(ndivide), geometry%z(i), geometry%z(j)

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
     print*, ' number of diffusive divides:', icountdif
     print*, ' number of slope divides:', icountslope
     !    print*, ' Average Distance between nodes: ', avedivide/dble(ndivide)
  endif
!print*, 'finished with divides'

else ! if mh/n not equal 1
  do iii=1, maxpass
     ncapture=0
     do ii=1,geometry%nnode ! loop over the nodes
        i=ishuffle(ii)
        do k=geometry%nb(i),1,-1 ! loop over all neighbours
           j=geometry%nn(k,i)
           if(network%receiver(j).ne.i.and.network%receiver(i).ne.j) then ! only for connections that are not part of the drainage network
              if (geometry%z(i).lt.geometry%z(j)) then ! check that i is lower than j
                 l=dsqrt((geometry%x(i)-geometry%x(j))**2.d0+(geometry%y(i)-geometry%y(j))**2.d0) ! l is length betwenn i and j
                 if(l.gt.params%xc) then
                    if (params%transient_divide) then
                       if (params%h*params%m.eq.1) then
                          tau=log(l/params%xc)/ &
                               (geometry%k(i)*(geometry%precipitation(i)*params%ka)**(params%m))
                       else
                          tau= (l**(1.d0-params%h*params%m) - params%xc**(1.d0-params%h*params%m)) /&
                               (geometry%k(i)*(geometry%precipitation(i)*params%ka)**(params%m))/(1.d0-params%h*params%m)
                       endif
                       num_samples=INT(min(tau,params%time)/params%deltat)
                       if (num_samples.eq.0) then
                          fact1=(geometry%erosion_rate(i)/geometry%k(i)/ &
                                  (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                       else
                          if ((params%num_bins-1)*params%sample_per_bin.lt.num_samples) then
                             num_full_bin=params%num_bins
                             sample_in_last_bin=num_samples-(params%num_bins-1)*params%sample_per_bin
                          else
                             num_full_bin = INT(num_samples/params%sample_per_bin)+1
                             sample_in_last_bin=num_samples-(num_full_bin-1)*params%sample_per_bin
                          endif
                          ave_erosion=0.d0
                          if (num_full_bin.eq.1) then
                             ave_erosion=geometry%erosion_rate_history(i,1)
                          else
                             ave_erosion=0.d0
                             do h=1,num_full_bin-1
                                ave_erosion=ave_erosion +geometry%erosion_rate_history(i,h)*params%sample_per_bin
                             enddo
                             ave_erosion=ave_erosion +geometry%erosion_rate_history(i,num_full_bin)*sample_in_last_bin
                             ave_erosion=ave_erosion/((num_full_bin-1)*params%sample_per_bin + sample_in_last_bin)
                          endif
                          fact1=(ave_erosion/geometry%k(i)/ &
                               (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                       endif
                    else
                       fact1=(geometry%erosion_rate(i)/geometry%k(i)/ &
                            (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)  
                    endif
                    fact2=(1.d0/params%hmn)*((l)**(params%hmn)-(params%xc)**(params%hmn))
                    zt=params%xc*params%tanthetac ! zt is channel head elevation
                    if(params%diffusion)then
                       !erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)  
                       if (params%transient_divide) then
                          ztd=ave_erosion*params%xc**2.d0/params%diffusivity
                       else
                          ztd=geometry%erosion_rate(i)*params%xc**2.d0/params%diffusivity   
                       endif
                       if(ztd.lt.zt)then
                          zt=max(ztd,params%xc*params%min_tan_head_slope)
                       endif
                    endif
                    zt=zt+geometry%z(i)+fact1*fact2 ! zt is test elevation
                 else
                    zt=l*params%tanthetac
                    if(params%diffusion)then
                       !erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)
                       ztd=geometry%erosion_rate(i)*l**2.d0/params%diffusivity
                       if(ztd.lt.zt)then
                          zt=max(ztd,l*params%min_tan_head_slope)
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
                    if (orig_receiver(j).eq.i) then
                       !print*, i, 'is caturing its former donor', j
                    endif
                    if (params%capture) then
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
        !*** Old calculations of erosion rates during capture - L.G.***
        !geometry%erosion_rate(capnode(icap))=(geometry%z(capnode(icap))-capelev(icap))/params%deltat
        !geometry%erosion_rate(capnode(icap))=(orig_z(capnode(icap))-capelev(icap))/params%deltat

        !*** New calculations of erosion rates during capture - L.G.***
        ! In the calculation of the erosion rate add the elevation drop fue to capture to previous 
        ! elevation drops both due to othetr captures in this loop and due to fluvial erosion in the former time step
       
        if (params%transient_divide) then
           cpn=capnode(icap) !cpn is the captured node
           ! need to update the first bin of the captured node with the correct erosion rate
           ! due to capture. First remove the last added erosion rate
           num_full_bin = INT(params%istep/params%sample_per_bin)+1
           if (num_full_bin.eq.1) then
              sample_in_last_bin=params%istep
              geometry%erosion_rate_history(cpn,1)= &
                      geometry%erosion_rate_history(cpn,1) - geometry%erosion_rate(cpn)/sample_in_last_bin
           else
              geometry%erosion_rate_history(cpn,1)= &
                   geometry%erosion_rate_history(cpn,1) - geometry%erosion_rate(cpn)/(params%sample_per_bin)
           endif
        endif

        geometry%erosion_rate(capnode(icap)) = (geometry%erosion_rate(capnode(icap))*params%deltat +&
             &geometry%z(capnode(icap))-capelev(icap))/params%deltat

        !add to the first bin the updated erosion rate
        if (params%transient_divide) then
           cpn=capnode(icap) !cpn is the captured node
           ! need to update the first bin of the captured node with the correct erosion rate
           ! due to capture. First remove the last added erosion rate
           num_full_bin = INT(params%istep/params%sample_per_bin)+1
           if (num_full_bin.eq.1) then
              sample_in_last_bin=params%istep
              geometry%erosion_rate_history(cpn,1)= &
                   geometry%erosion_rate_history(cpn,1) + geometry%erosion_rate(cpn)/sample_in_last_bin
           else
              geometry%erosion_rate_history(cpn,1)= &
                   geometry%erosion_rate_history(cpn,1) + geometry%erosion_rate(cpn)/(params%sample_per_bin)
           endif
        endif


        geometry%z(capnode(icap))=capelev(icap)
        network%receiver(capnode(icap))=recnode(icap)
     enddo
  enddo
  print*, 'warning max passes exceeded'
1493 continue
  !         ******************  start loop for divide heights      **********************************
  !
  ! set a tolerance for the ridge height convergence
  ztoler=.1d0
  icountdif=0
  icountslope=0
  do i=1,geometry%nnode ! loop over the nodes
     do k=geometry%nb(i),1,-1 ! loop over all neighbours
        j=geometry%nn(k,i)
        if(network%receiver(j).ne.i.and.network%receiver(i).ne.j) then ! only for connections that are not part of the drainage network
           if (geometry%z(i).lt.geometry%z(j)) then ! check that i is lower than j
              l=dsqrt((geometry%x(i)-geometry%x(j))**2.d0+(geometry%y(i)-geometry%y(j))**2.d0) ! l is length of segment betwenn i and j
              ndivide=ndivide+1
              avedivide=avedivide+l
              ! start here with new diffusion algorithm

              ! give initial guess for divide and endpoints
              x1=0.d0
              x2=l
              zdi=0.d0
              zdj=2.d0*ztoler
              icount=0
              ! find ridge heights at opposing node positions
              !  find ridge height from node i
              xc=min(params%xc,l)
              !erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)    
               if (params%transient_divide.and.l.gt.params%xc) then
                 if (params%h*params%m.eq.1) then
                    tau=log(l/params%xc)/ &
                         (geometry%k(i)*(geometry%precipitation(i)*params%ka)**(params%m))                   
                 else
                    tau= (l**(1.d0-params%h*params%m) - params%xc**(1.d0-params%h*params%m)) /&
                         (geometry%k(i)*(geometry%precipitation(i)*params%ka)**(params%m))/(1.d0-params%h*params%m)                  
                 endif
                 !if (tau.lt.params%time) print*, 'tau is smaller'
                 num_samples=INT(min(tau,params%time)/params%deltat)
                 if (num_samples.eq.0) then
                    local_erosioni=geometry%erosion_rate(i)
                 else
                    if ((params%num_bins-1)*params%sample_per_bin.lt.num_samples) then
                       num_full_bin=params%num_bins
                       sample_in_last_bin=num_samples-(params%num_bins-1)*params%sample_per_bin
                    else
                       num_full_bin = INT(num_samples/params%sample_per_bin)+1
                       sample_in_last_bin=num_samples-(num_full_bin-1)*params%sample_per_bin
                    endif
                    ave_erosion=0.d0
                    if (num_full_bin.eq.1) then
                       ave_erosion=geometry%erosion_rate_history(i,1)
                    else
                       ave_erosion=0.d0
                       do h=1,num_full_bin-1
                          ave_erosion=ave_erosion +geometry%erosion_rate_history(i,h)*params%sample_per_bin
                       enddo
                       ave_erosion=ave_erosion +geometry%erosion_rate_history(i,num_full_bin)*sample_in_last_bin
                       ave_erosion=ave_erosion/((num_full_bin-1)*params%sample_per_bin + sample_in_last_bin)
                    endif
                    local_erosioni=ave_erosion
                 endif
                 if (params%h*params%m.eq.1) then
                    tau=log(l/params%xc)/ &
                         (geometry%k(j)*(geometry%precipitation(j)*params%ka)**(params%m))                   
                 else
                    tau= (l**(1.d0-params%h*params%m) - params%xc**(1.d0-params%h*params%m)) /&
                         (geometry%k(j)*(geometry%precipitation(j)*params%ka)**(params%m))/(1.d0-params%h*params%m)                  
                 endif
                 !if (tau.lt.params%time) print*, 'tau is smaller'
                 num_samples=INT(min(tau,params%time)/params%deltat)
                 if (num_samples.eq.0) then
                    local_erosionj=geometry%erosion_rate(j)
                 else
                    if ((params%num_bins-1)*params%sample_per_bin.lt.num_samples) then
                       num_full_bin=params%num_bins
                       sample_in_last_bin=num_samples-(params%num_bins-1)*params%sample_per_bin
                    else
                       num_full_bin = INT(num_samples/params%sample_per_bin)+1
                       sample_in_last_bin=num_samples-(num_full_bin-1)*params%sample_per_bin
                    endif
                    ave_erosion=0.d0
                    if (num_full_bin.eq.1) then
                       ave_erosion=geometry%erosion_rate_history(j,1)
                    else
                       ave_erosion=0.d0
                       do h=1,num_full_bin-1
                          ave_erosion=ave_erosion +geometry%erosion_rate_history(j,h)*params%sample_per_bin
                       enddo
                       ave_erosion=ave_erosion +geometry%erosion_rate_history(j,num_full_bin)*sample_in_last_bin
                       ave_erosion=ave_erosion/((num_full_bin-1)*params%sample_per_bin + sample_in_last_bin)
                    endif
                    local_erosionj=ave_erosion
                 endif
              else
                 local_erosioni=geometry%erosion_rate(i)
                 local_erosionj=geometry%erosion_rate(j)  
              endif
              if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/local_erosioni))then 
                 zc=max(local_erosioni*xc**2.d0/params%diffusivity,xc*params%min_tan_head_slope)
              else
                 zc=xc*params%tanthetac
              endif
              if(l.le.params%xc)then
                 zi2=geometry%z(i)+zc
              else
                 fact2=(1.d0/params%hmn)*(l**(params%hmn)-(params%xc)**(params%hmn))
                 fact3=(local_erosioni/geometry%k(i)/ &
                      (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                 zi2=geometry%z(i)+zc+fact2*fact3
              endif
              zi1=geometry%z(i)
              !  find ridge height from node j
              xc=min(params%xc,(l))
              !erodrate=max(geometry%erosion_rate(j),params%min_erosion_rate)
              if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/local_erosionj))then 
                 zc=max(local_erosionj*xc**2.d0/params%diffusivity,xc*params%min_tan_head_slope)
              else
                 zc=xc*params%tanthetac
              endif
              if((l).le.params%xc)then
                 zj1=geometry%z(j)+zc
              else
                 fact2=(1.d0/params%hmn)*(l**(params%hmn)-(params%xc)**(params%hmn))
                 fact3=(local_erosionj/geometry%k(j)/ &
                      (geometry%precipitation(j)*params%ka)**params%m)**(1.d0/params%n)
                 zj1=geometry%z(j)+zc+fact2*fact3
              endif
              zj2=geometry%z(j)
              if((zi1-zj1).ge.0.d0.or.(zi2-zj2).le.0.d0)then
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
                    xd=l/2.d0
                    go to 334
                 endif

                 !  find ridge height from node i
                 xc=min(params%xc,xd)
                 !erodrate=max(geometry%erosion_rate(i),params%min_erosion_rate)
                 if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/local_erosioni))then 
                    zc=max(local_erosioni*xc**2.d0/params%diffusivity,xc*params%min_tan_head_slope)
                    dflag=1
                 else
                    zc=xc*params%tanthetac
                 endif
                 if(xd.le.params%xc)then
                    zdi=geometry%z(i)+zc
                 else
                    !fact2=log(params%xc/((xd)))
                    fact2 = (1.d0/params%hmn)*(xd**(params%hmn)-(params%xc)**(params%hmn))
                    fact3=(local_erosioni/geometry%k(i)/ &
                         (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                    zdi=geometry%z(i)+zc+fact2*fact3
                 endif
                 !  find ridge height from node j
                 xc=min(params%xc,(l-xd))
                 !erodrate=max(geometry%erosion_rate(j),params%min_erosion_rate)
                 if(params%diffusion.and.xc.lt.(params%diffusivity*params%tanthetac/local_erosionj))then 
                    zc=max(local_erosionj*xc**2.d0/params%diffusivity,xc*params%min_tan_head_slope)
                    dflag=1
                 else
                    zc=xc*params%tanthetac
                 endif
                 if((l-xd).le.params%xc)then
                    zdj=geometry%z(j)+zc
                 else
                    !fact2=log(params%xc/((l-xd)))
                    fact2 = (1.d0/params%hmn)*((l-xd)**(params%hmn)-(params%xc)**(params%hmn))
                    fact3=(local_erosionj/geometry%k(j)/ &
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

334           continue
              if(dflag.eq.1)then
                 icountdif=icountdif+1
              else
                 icountslope=icountslope+1
              endif


              !! checking that ridge calculation is ok
!!$                 if ((params%istep/params%freq)*params%freq.eq.params%istep) then
!!$                    print*, i,geometry%x(i),geometry%y(i),geometry%z(i)
!!$                    print*, j,geometry%x(j),geometry%y(j),geometry%z(j)
!!$                    print*, xc*params%tanthetac, fact3, params%hmn, l 
!!$                    print*, xd, zdi,zdj
!!$                    print*, '^^^^^^^^'
!!$                 endif

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
                    geometry%xdiv(ndivide) = geometry%x(i) + xd/dsqrt(1.0d0+local_slope**2.d0)
                 else
                    geometry%xdiv(ndivide) = geometry%x(i) - xd/dsqrt(1.0d0+local_slope**2.d0) 
                 endif
                 geometry%ydiv(ndivide) = local_slope*(geometry%xdiv(ndivide)-geometry%x(i))+ geometry%y(i)
              endif
              geometry%zdiv(ndivide) = max(zdi,zdj) 
              if (geometry%zdiv(ndivide).lt.geometry%z(i).or.&
                   &geometry%zdiv(ndivide).lt.geometry%z(j)) then
                 print*, 'divide', ndivide, 'is lower than neighboring verteces'
                 print*, geometry%zdiv(ndivide), geometry%z(i), geometry%z(j)
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
   !print*,ncapglobe,'captures'
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


!print*, maxval(geometry%zdiv)

call time_out ('captures_and_divides')

return

end subroutine captures_and_divides
