subroutine add_remove_nodes(geometry,network,params,delaunay)

  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (parm) params
  type (del) delaunay

  double precision l,xx,xd,fact3,fact2
  integer i,j,k,h, nrem,nrem1,nnold,nadd1,nadd2,nadd3,nadd4,nnodecheck,icount,jj,flag
  double precision xmin,xmax,ymin,ymax,min_dist,curr_dist
  integer,dimension(:),allocatable::bneighbor

  call time_in ('add_remove_nodes')



 

  if (params%move_points) then 
     xmin = minval(geometry%x)
     xmax = maxval(geometry%x)
     ymin = minval(geometry%y)
     ymax = maxval(geometry%y)
     allocate (bneighbor(geometry%nnode_max))
  endif

  ! This section is to add nodes on channels where the channel connect two nodes that are not neighbors.
  nadd3=0
  nnold=geometry%nnode
  do i=1,nnold ! loop over the nodes
     if(geometry%fix(i).ne.1.and.network%receiver(i).ne.0)then
        flag=0
        do k=geometry%nb(i),1,-1 ! loop over the neighbors
           if(network%receiver(i).eq.geometry%nn(k,i))flag=1  ! receiver is a neighbor
        enddo
        if(flag.eq.0)then
           j=network%receiver(i)
           nadd3=nadd3+1
           geometry%nnode=geometry%nnode+1
           if (geometry%nnode.gt.geometry%nnode_max) then
              print*,'problem while adding nodes'
              STOP 'too many nodes added. Increase nnode_max'
           endif
           geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.
           geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.
           geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.
           print*, 'adding node', geometry%nnode, 'via add3 in coordinate'
           print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
           ! Assign physical properties to new node as average of its endmembers
           
           geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.
           geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.
           geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.           
           network%receiver(geometry%nnode)=j
           network%receiver(i)=geometry%nnode
           geometry%fix(geometry%nnode)=0
           geometry%surface(geometry%nnode)=0.0
           geometry%discharge(geometry%nnode)=0.0
           geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.
           geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.
           geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.
           geometry%sediment_flux(geometry%nnode)=0.0
        endif
     endif
  enddo
  !  if (nadd3.ne.0) print*,params%time,nadd3,'nodes added to keep neighbor',geometry%nnode

  ! This section is to add nodes on channels. It checks whether nodes on the drainage network
  ! are too far apart or not. I adds nodes when the distance between two nodes is
  ! greater than a critical distance "lmax" (chosen in initialize_parameters.f90)

  nadd1=0
  do i=1,nnold ! loop over the nodes
     !       if (geometry%fix(i).ne.1) then
     if (network%receiver(i).ne.0) then ! only for connections that are part of the drainage network
        j=network%receiver(i)
        l=sqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length of segment betwenn i and j
        if (l.gt.params%lmax) then ! add a node on the connection
           nadd1=nadd1+1
           geometry%nnode=geometry%nnode+1
           if (geometry%nnode.gt.geometry%nnode_max) then
              print*,'problem while adding nodes'
              STOP 'too many nodes added. Increase nnode_max'
           endif
           ! Add a random perturbation in location of new node equal to 50% of the max length between nodes
           call random_number(xx)
           geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2. &
                +(2.*xx-1.)*0.05*params%lmax 
           call random_number(xx)
           geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2. &
                +(2.*xx-1.)*0.05*params%lmax
           geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.
           print*, 'adding node', geometry%nnode, 'via add1 in coordinate'
           print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
           !note every physical properties will have to be updated at each new
           !node....
           ! Assign physical properties to new node as average of its endmembers
           ! added 15.6.2010 sdw  Not sure if these are all the node-local properties
           geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.
           geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.
           geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.
           geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.
           network%receiver(geometry%nnode)=j
           network%receiver(i)=geometry%nnode
           geometry%fix(geometry%nnode)=0
           geometry%surface(geometry%nnode)=0.0
           geometry%discharge(geometry%nnode)=0.0
           geometry%precipitation(geometry%nnode)=geometry%precipitation(i)
           geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.
           !           geometry%nb(geometry%nnode)=2
           !           geometry%nn(1,geometry%nnode)=i
           !           geometry%nn(2,geometry%nnode)=j
           geometry%sediment_flux(geometry%nnode)=0.0
           !           network%ndon(geometry%nnode)=1
           !           network%donors(1,geometry%nnode)=j
        endif
     endif
     !       endif
  enddo




  ! This section is to add boundary nodes when the boundaries are also advected. 
  ! It adds nodes when the distance between two boundary nodes is
  ! greater than a critical distance "lmax" (chosen in initialize_parameters.f90)

  ! for some reason the nb array for boundary nodes is not trustable. Need to generate a new list.
  ! Each node will point to its neighbor in anticlockwise direction.  
  nadd4=0
  bneighbor = 0
  if (params%move_points) then
     do i=1,nnold
        if (geometry%fix(i).eq.1) then
           curr_dist = max(xmax,ymax)
           if (geometry%y(i).eq.ymin) then
              do j = 1,nnold
                 if (j.ne.i.and.geometry%fix(j).eq.1.and.geometry%y(j).eq.ymin.and.geometry%x(j).gt.geometry%x(i)) then
                    l = sqrt((geometry%x(i)-geometry%x(j))**2)
                    if (l.lt.curr_dist) then
                       bneighbor(i) = j
                       curr_dist = l
                    endif
                 endif
              enddo
              if (bneighbor(i).ne.0.and. curr_dist.gt.params%lmax) then
                 nadd4=nadd4+1
                 geometry%nnode=geometry%nnode+1
                 if (geometry%nnode.gt.geometry%nnode_max) then
                    print*,'problem while adding nodes'
                    STOP 'too many nodes added. Increase nnode_max'
                 endif
                 j = bneighbor(i)
                 geometry%y(geometry%nnode)=geometry%y(i)
                 call random_number(xx)
                 geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2. &
                      +(2.*xx-1.)*0.05*params%lmax
                 geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.
                 print*, 'adding boundary node at coordinates', geometry%x(geometry%nnode),geometry%y(geometry%nnode) 
                 print*, 'neighbor to', i, geometry%x(i),geometry%y(i)
                 print*, 'and to', bneighbor(i), geometry%x(bneighbor(i)),geometry%y(bneighbor(i))
                 print*, 'new node index is', geometry%nnode 
                 print*, 'receiver of i', network%receiver(i), 'receiver of j', network%receiver(j)
                 geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.
                 geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.
                 geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.
                 geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.
                 geometry%fix(geometry%nnode)=1
                 geometry%surface(geometry%nnode)=0.0
                 geometry%discharge(geometry%nnode)=0.0
                 geometry%precipitation(geometry%nnode)=geometry%precipitation(i)
                 geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.
                 geometry%sediment_flux(geometry%nnode)=0.0
                 network%receiver(geometry%nnode)=0
              endif
           endif
           if (geometry%x(i).eq.xmax) then
              do j = 1,nnold
                 if (j.ne.i.and.geometry%fix(j).eq.1.and.geometry%x(j).eq.xmax.and.geometry%y(j).gt.geometry%y(i)) then
                    l = sqrt((geometry%y(i)-geometry%y(j))**2)
                    if (l.lt.curr_dist) then
                       bneighbor(i) = j
                       curr_dist = l
                    endif
                 endif
              enddo
              if (bneighbor(i).ne.0.and.curr_dist.gt.params%lmax) then
                 nadd4=nadd4+1
                 geometry%nnode=geometry%nnode+1
                 if (geometry%nnode.gt.geometry%nnode_max) then
                    print*,'problem while adding nodes'
                    STOP 'too many nodes added. Increase nnode_max'
                 endif
                 j = bneighbor(i)
                 geometry%x(geometry%nnode)=geometry%x(i)
                 call random_number(xx)              
                 geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2. &
                      +(2.*xx-1.)*0.05*params%lmax
                 geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.
                 print*, 'adding boundary node at coordinates', geometry%x(geometry%nnode),geometry%y(geometry%nnode) 
                  print*, 'neighbor to', i, geometry%x(i),geometry%y(i)
                 print*, 'and to', bneighbor(i), geometry%x(bneighbor(i)),geometry%y(bneighbor(i))
                 geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.
                 geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.
                 geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.
                 geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.
                 geometry%fix(geometry%nnode)=1
                 geometry%surface(geometry%nnode)=0.0
                 geometry%discharge(geometry%nnode)=0.0
                 geometry%precipitation(geometry%nnode)=geometry%precipitation(i)
                 geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.
                 geometry%sediment_flux(geometry%nnode)=0.0
                 network%receiver(geometry%nnode)=0
              endif
           endif
           if (geometry%y(i).eq.ymax) then
              do j = 1,nnold
                 if (j.ne.i.and.geometry%fix(j).eq.1.and.geometry%y(j).eq.ymax.and.geometry%x(j).lt.geometry%x(i)) then
                    l = sqrt((geometry%x(i)-geometry%x(j))**2)
                    if (l.lt.curr_dist) then
                       bneighbor(i) = j
                       curr_dist = l
                    endif
                 endif
              enddo
              if (bneighbor(i).ne.0.and.curr_dist.gt.params%lmax) then
                 nadd4=nadd4+1
                 geometry%nnode=geometry%nnode+1
                 if (geometry%nnode.gt.geometry%nnode_max) then
                    print*,'problem while adding nodes'
                    STOP 'too many nodes added. Increase nnode_max'
                 endif
                 j = bneighbor(i)
                 geometry%y(geometry%nnode)=geometry%y(i)
                 call random_number(xx)              
                 geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2. &
                      +(2.*xx-1.)*0.05*params%lmax
                 geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.
                 print*, 'adding boundary node at coordinates', geometry%x(geometry%nnode),geometry%y(geometry%nnode) 
                 print*, 'neighbor to', i, geometry%x(i),geometry%y(i)
                 print*, 'and to', bneighbor(i), geometry%x(bneighbor(i)),geometry%y(bneighbor(i))
                 geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.
                 geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.
                 geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.
                 geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.
                 geometry%fix(geometry%nnode)=1
                 geometry%surface(geometry%nnode)=0.0
                 geometry%discharge(geometry%nnode)=0.0
                 geometry%precipitation(geometry%nnode)=geometry%precipitation(i)
                 geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.
                 geometry%sediment_flux(geometry%nnode)=0.0
                 network%receiver(geometry%nnode)=0
              endif
           endif
           if (geometry%x(i).eq.xmin) then
              do j = 1,nnold
                 if (j.ne.i.and.geometry%fix(j).eq.1.and.geometry%x(j).eq.xmin.and.geometry%y(j).lt.geometry%y(i)) then
                    l = sqrt((geometry%y(i)-geometry%y(j))**2)
                    if (l.lt.curr_dist) then
                       bneighbor(i) = j
                       curr_dist = l
                    endif
                 endif
              enddo
              if (bneighbor(i).ne.0.and.curr_dist.gt.params%lmax) then
                 nadd4=nadd4+1
                 geometry%nnode=geometry%nnode+1
                 if (geometry%nnode.gt.geometry%nnode_max) then
                    print*,'problem while adding nodes'
                    STOP 'too many nodes added. Increase nnode_max'
                 endif
                 j = bneighbor(i)
                 geometry%x(geometry%nnode)=geometry%x(i)
                 call random_number(xx)              
                 geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2. &
                      +(2.*xx-1.)*0.05*params%lmax
                 geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.
                 print*, 'adding boundary node at coordinates', geometry%x(geometry%nnode),geometry%y(geometry%nnode) 
                 print*, 'neighbor to', i, geometry%x(i),geometry%y(i)
                 print*, 'and to', bneighbor(i), geometry%x(bneighbor(i)),geometry%y(bneighbor(i))
                 geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.
                 geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.
                 geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.
                 geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.
                 geometry%fix(geometry%nnode)=1
                 geometry%surface(geometry%nnode)=0.0
                 geometry%discharge(geometry%nnode)=0.0
                 geometry%precipitation(geometry%nnode)=geometry%precipitation(i)
                 geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.
                 geometry%sediment_flux(geometry%nnode)=0.0
                 network%receiver(geometry%nnode)=0
              endif
           endif
        endif
     enddo
  endif




  !  if (nadd1.ne.0) print*,params%time,nadd1,'nodes added in channels',geometry%nnode
  nadd2 = 0
  if(nadd1.eq.0.and.nadd3.eq.0)then ! this condition needed to stop potentially crossing rivers
     ! add node on non-channel divides that are longer than a given distance
     !  add node to lower existing node

     do i=1,nnold ! loop over the original nodes
        if(geometry%nb(i).eq.0)print*,' zero neighbor error in add nodes' ! debug line
        do k=geometry%nb(i),1,-1 ! loop over the neighbors
           j=geometry%nn(k,i)
           !          if(network%receiver(i).ne.j.and.network%receiver(j).ne.i.and.geometry%fix(i).ne.1.)then ! check non-channel,non-boundary connections
           if (geometry%fix(i).ne.1.or.geometry%fix(j).ne.1) then !when both are boundary nodes treated with nadd4
              if(network%receiver(i).ne.j.and.network%receiver(j).ne.i)then ! check non-channel
                 if(geometry%surface_share(k,i).gt.0)then  ! only continue if divide is closer to j than i
                    l=sqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length of segment betwenn i and j
                    if (l.gt.params%ldivmax) then ! add node on the connection
                       xd=(geometry%surface_share(k,i)+1)*l/2.
                       if(xd.gt.params%xc)then  !  continue only if there is a fluvial segment to connection
                          nadd2=nadd2+1
                          geometry%nnode=geometry%nnode+1
                          if (geometry%nnode.gt.geometry%nnode_max) then
                             print*,'problem while adding nodes'
                             STOP 'too many nodes added. Increase nnode_max'
                          endif
                          geometry%x(geometry%nnode)=geometry%x(i)+(xd-params%xc)/l*(geometry%x(j)-geometry%x(i))  
                          geometry%y(geometry%nnode)=geometry%y(i)+(xd-params%xc)/l*(geometry%y(j)-geometry%y(i)) 
                          if (params%hmn.eq.0.d0) then
                             fact2=log(params%xc/((xd)))
                             fact3=(geometry%erosion_rate(i)/geometry%k(i)/ &
                                  (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                             geometry%z(geometry%nnode)=geometry%z(i)-fact2*fact3  ! note that this is only the hm/n=1 solution Needs changing
                          else
                             !fact2=log(params%xc/((xd)))
                             fact2=(1./params%hmn)*(xd**(params%hmn)-(params%xc)**(params%hmn))
                             fact3=(geometry%erosion_rate(i)/geometry%k(i)/ &
                                  (geometry%precipitation(i)*params%ka)**params%m)**(1.d0/params%n)
                             geometry%z(geometry%nnode)=geometry%z(i)+fact2*fact3 
                          endif
                          print*, 'adding node', geometry%nnode, 'via add2 in coordinate'
                          print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)

                          ! Assign physical properties to new node from neighbor node i
                          !  sdw  Not sure if these are all the node-local properties
                          geometry%erosion_rate(geometry%nnode)=geometry%erosion_rate(i)
                          geometry%u(geometry%nnode)=geometry%u(i)
                          geometry%v(geometry%nnode)=geometry%v(i)
                          geometry%w(geometry%nnode)=geometry%w(i)
                          geometry%erosion_rate(geometry%nnode)=geometry%erosion_rate(i)
                          geometry%fix(geometry%nnode)=0
                          !                  geometry%surface(geometry%nnode)=params%amin*2. !set artificially large to avoid removal in check below
                          geometry%discharge(geometry%nnode)=0.0
                          geometry%precipitation(geometry%nnode)=geometry%precipitation(i)
                          geometry%k(geometry%nnode)=geometry%k(i)
                          !                  geometry%nb(geometry%nnode)=1
                          !                  geometry%nn(1,geometry%nnode)=i
                          geometry%sediment_flux(geometry%nnode)=0.0
                          ! update donor count for node i
                          network%ndon(i)=network%ndon(i)+1
                          network%ndon(geometry%nnode)=0
                          network%receiver(geometry%nnode)=i
                       endif
                    endif
                 endif
              endif
           endif
        enddo
     enddo
     !if (nadd2.ne.0) print*,params%time,nadd2,'nodes added in channel heads',geometry%nnode
  endif


   



  ! remove nodes based on minimum surface area in headwater
  nrem=0
  !do i=1,geometry%nnode
  !print*,i,network%receiver(i)
  !enddo


  do i=1, nnold
     if(i.le.geometry%nnode-nadd1-nadd2-nadd3-nadd4)then
        ! print*,i, network%receiver(i)
        if(geometry%fix(i).ne.1.and.geometry%surface(i).lt.params%amin.and.network%ndon(i).eq.0)then
           nrem=nrem+1
           geometry%nnode=geometry%nnode-1
           
           print*,'remove node ',i,'too small drainage area'
           print*, 'number of nodes is', geometry%nnode
           !      network%ndon(network%receiver(i))=network%ndon(network%receiver(i))-1
           do j=1,i-1
              if (network%receiver(j).gt.i) network%receiver(j)=network%receiver(j)-1
           enddo
           if(i.le.geometry%nnode)then
              do j=i,geometry%nnode
                 geometry%x(j)=geometry%x(j+1)
                 geometry%y(j)=geometry%y(j+1)
                 geometry%z(j)=geometry%z(j+1)
                 geometry%u(j)=geometry%u(j+1)
                 geometry%v(j)=geometry%v(j+1)
                 geometry%w(j)=geometry%w(j+1)
                 geometry%fix(j)=geometry%fix(j+1)
                 geometry%surface(j)=geometry%surface(j+1)
                 geometry%discharge(j)=geometry%discharge(j+1)
                 geometry%precipitation(j)=geometry%precipitation(j+1)
                 geometry%k(j)=geometry%k(j+1)
                 geometry%nb(j)=geometry%nb(j+1)
                 geometry%erosion_rate(j)=geometry%erosion_rate(j+1)
                 geometry%sediment_flux(j)=geometry%sediment_flux(j+1)
                 network%ndon(j)=network%ndon(j+1)
                 if (network%receiver(j+1).gt.i) then
                    network%receiver(j)=network%receiver(j+1)-1
                 else
                    network%receiver(j)=network%receiver(j+1)
                 endif

                 !          do k=1, geometry%nb(j)
                 !            geometry%nn(k,j)=geometry%nn(k,j+1)
                 !          enddo
              enddo
           endif
        endif
     endif
  enddo


if (nrem.gt.0.or.nadd1.gt.0.or.nadd2.gt.0.or.nadd3.gt.0.or.nadd4.gt.0)then
     network%nnode=geometry%nnode
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




  ! remove nodes that approach the boundaries due to advection (due to advection)

 nrem1 = 0
  if (params%move_points) then 
     
     min_dist =  params%max_adv*params%deltat
     do i=1, geometry%nnode !also new nodes that were just added
554     if (i.le.geometry%nnode) then
           if(geometry%fix(i).ne.1.or.&
                &geometry%y(i).eq.ymin.and.geometry%x(i).ne.xmin.and.geometry%x(i).ne.xmax.or.&
                &geometry%y(i).eq.ymax.and.geometry%x(i).ne.xmin.and.geometry%x(i).ne.xmax) then
              if (abs(geometry%x(i)-xmin).lt.min_dist.or.&
                   &abs(geometry%x(i)-xmax).lt.min_dist.or.&
                   &abs(geometry%y(i)-ymin).lt.min_dist.and.geometry%y(i).ne.ymin.or.&
                   &abs(geometry%y(i)-ymax).lt.min_dist.and.geometry%y(i).ne.ymax) then
                 print*, 'removing node',i,'with coordinates', geometry%x(i), geometry%y(i)
                 nrem1=nrem1+1
                 geometry%nnode=geometry%nnode-1
                 print*, 'updated number of nodes is',geometry%nnode
                 ! first reorganize network
                 k = network%receiver(i)
                 if (k.ne.0) then !i is not a lake
                    do j = 1,network%ndon(i)
                       h = network%donors(j,i) ! h is a donor of i
                       network%receiver(h) = k ! now h is a donor of k
                       !network%ndon(k) = network%ndon(k)+1
                    enddo
                 else! i is a lake
                    do j = 1,network%ndon(i)
                       h = network%donors(j,i) ! h is a donor of i
                       network%receiver(h) = 0 ! now h is a lake
                    enddo
                 endif
                 do j=1,i-1
                    if (network%receiver(j).gt.i) network%receiver(j)=network%receiver(j)-1
                 enddo
                 if(i.le.geometry%nnode)then
                    do j=i,geometry%nnode
                       geometry%x(j)=geometry%x(j+1)
                       geometry%y(j)=geometry%y(j+1)
                       geometry%z(j)=geometry%z(j+1)
                       geometry%u(j)=geometry%u(j+1)
                       geometry%v(j)=geometry%v(j+1)
                       geometry%w(j)=geometry%w(j+1)
                       geometry%fix(j)=geometry%fix(j+1)
                       geometry%surface(j)=geometry%surface(j+1)
                       geometry%discharge(j)=geometry%discharge(j+1)
                       geometry%precipitation(j)=geometry%precipitation(j+1)
                       geometry%k(j)=geometry%k(j+1)
                       geometry%nb(j)=geometry%nb(j+1)
                       geometry%erosion_rate(j)=geometry%erosion_rate(j+1)
                       geometry%sediment_flux(j)=geometry%sediment_flux(j+1)
                       network%ndon(j)=network%ndon(j+1)
                       if (network%receiver(j+1).gt.i) then
                          network%receiver(j)=network%receiver(j+1)-1
                       else
                          network%receiver(j)=network%receiver(j+1)
                       endif
                    enddo
                 endif
                 network%nnode=geometry%nnode
                 network%donors=0
                 network%ndon=0
                 do j=1,network%nnode
                    k=network%receiver(j)
                    if (k.ne.0) then
                       network%ndon(k)=network%ndon(k)+1
                       network%donors(network%ndon(k),k)=j
                    endif
                 enddo
                 go to 554
              endif
           endif
        endif
     enddo
  endif





  ! remove nodes that approach the boundaries due to advection (due to advection)
!!$  nrem1=0
!!$  if (params%move_points) then 
!!$     do i=1, nnold
!!$        if(i.le.geometry%nnode-nadd1-nadd2-nadd3 .and. geometry%fix(i).ne.1)then
!!$           k = network%receiver(i)
!!$           if (geometry%fix(k).eq.1) then
!!$              l = sqrt((geometry%x(i)-geometry%x(k))**2+(geometry%y(i)-geometry%y(k))**2)
!!$              if (l.lt.params%lmax/5) then 
!!$                 print*, 'removing node',i,'because',l,'is too short'
!!$                 nrem1=nrem1+1
!!$                 geometry%nnode=geometry%nnode-1
!!$                 do j = 1,network%ndon(i)
!!$                    h = network%donors(j,i) ! h is a donor of i
!!$                    network%receiver(h) = k ! now h is a donor of k
!!$                    network%ndon(k) = network%ndon(k)+1
!!$                 enddo
!!$                 do j=1,i-1
!!$                    if (network%receiver(j).gt.i) network%receiver(j)=network%receiver(j)-1
!!$                 enddo
!!$                 if(i.le.geometry%nnode)then
!!$                    do j=i,geometry%nnode
!!$                       geometry%x(j)=geometry%x(j+1)
!!$                       geometry%y(j)=geometry%y(j+1)
!!$                       geometry%z(j)=geometry%z(j+1)
!!$                       geometry%u(j)=geometry%u(j+1)
!!$                       geometry%v(j)=geometry%v(j+1)
!!$                       geometry%w(j)=geometry%w(j+1)
!!$                       geometry%fix(j)=geometry%fix(j+1)
!!$                       geometry%surface(j)=geometry%surface(j+1)
!!$                       geometry%discharge(j)=geometry%discharge(j+1)
!!$                       geometry%precipitation(j)=geometry%precipitation(j+1)
!!$                       geometry%k(j)=geometry%k(j+1)
!!$                       geometry%nb(j)=geometry%nb(j+1)
!!$                       geometry%erosion_rate(j)=geometry%erosion_rate(j+1)
!!$                       geometry%sediment_flux(j)=geometry%sediment_flux(j+1)
!!$                       network%ndon(j)=network%ndon(j+1)
!!$                       if (network%receiver(j+1).gt.i) then
!!$                          network%receiver(j)=network%receiver(j+1)-1
!!$                       else
!!$                          network%receiver(j)=network%receiver(j+1)
!!$                       endif
!!$
!!$                       !          do k=1, geometry%nb(j)
!!$                       !            geometry%nn(k,j)=geometry%nn(k,j+1)
!!$                       !          enddo
!!$                    enddo
!!$                 endif
!!$              endif
!!$           endif
!!$        endif
!!$     enddo
!!$  endif






  !do jj=1,geometry%nnode
  !print*,jj,network%receiver(jj)
  !enddo
  !if(nrem.ne.0)print*,params%time, nrem,' nodes removed', geometry%nnode
  !do i=1,geometry%nnode
  !print*,i,network%receiver(i)
  !enddo
  ! update the donors
  !print*, nadd1, nadd2, nadd3, nadd4, nrem, nrem1
!!$  if (nrem.gt.0.or.nadd1.gt.0.or.nadd2.gt.0.or.nadd3.gt.0.or.nrem1.gt.0.or.nadd4.gt.0)then
!!$     network%nnode=geometry%nnode
!!$     network%donors=0
!!$     network%ndon=0
!!$     do i=1,network%nnode
!!$        k=network%receiver(i)
!!$        if (k.ne.0) then
!!$           network%ndon(k)=network%ndon(k)+1
!!$           network%donors(network%ndon(k),k)=i
!!$        endif
!!$     enddo
!!$  endif

  if(network%nnode.ne.geometry%nnode)print*, 'error in addremovenodes network nnode wrong'

  icount=0
  do i=1,geometry%nnode
     if(geometry%fix(i).eq.1)then
        icount=icount+1
        go to 555
     endif
     if(network%receiver(i).eq.0)then
        icount=icount+1
        go to 555
     endif
     j=network%receiver(i)
     if(j.gt.0.and.j.le.geometry%nnode)then
        icount=icount+1
        go to 555
     endif
     print*,' node not in network: ',i, network%receiver(i)
555  continue
  enddo
  !!print*,'number of nodes in network: ', icount
  !
  !
  !  do i=1,geometry%nnode
  !    if (network%receiver(i).eq.0.or.geometry%fix(i).eq.1) then
  !     nnodecheck=stack%nnode
  !     call add_to_stack (i,stack,network)
  !     if (nnodecheck.eq.stack%nnode) then
  !      print*,i,'this node is not in the stack'
  !      print*,geometry%x(i),geometry%y(i),geometry%z(i),network%receiver(i),network%ndon(i)
  !      pause
  !     endif
  !    endif
  !  enddo
  !if (stack%nnode.ne.geometry%nnode) then
  !  print*,stack%nnode,geometry%nnode,'stack nnode and geometry nnode are different'
  !endif
  !   do i=1, geometry%nnode
  !     do j=1, stack%nnode
  !     if(i.eq.stack%order(j))go to 111
  !     enddo
  !     print*,'add_remove node: missing node in stack', i
  !111 continue
  !   enddo
  !

  !      deallocate(stack%order)

  if (params%move_points) deallocate(bneighbor) 

  call time_out ('add_remove_nodes')

  return 
end subroutine add_remove_nodes
