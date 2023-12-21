subroutine add_remove_nodes(geometry,network,params,delaunay)

  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (parm) params
  type (del) delaunay

  double precision l,xx,xd,fact3,fact2
  integer i,j,k,h 
  integer nrem_small_area,nrem_close_to_boundary,nadd_on_channel,nadd_between_channel,nadd_river_no_connection,nadd_boundary
  integer nnold, icount,flag
  double precision xmin,xmax,ymin,ymax,min_dist
  integer num_samples, num_full_bin,samle_in_last_bin 
  double precision tau,ave_erosion 

  call time_in ('add_remove_nodes')


  !when advection nodes there are special add and remove operations for which the following is needed:
  if (params%move_points) then 
     xmin = minval(geometry%x)
     xmax = maxval(geometry%x)
     ymin = minval(geometry%y)
     ymax = maxval(geometry%y)
  endif


  nnold=geometry%nnode !original number of nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! First stage: take care of all the cases where a new node should added !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! This section is to add nodes on channels where the channel connect two nodes that are not neighbors.
  nadd_river_no_connection=0

  
  ! not part of cascade code
!!$  do i=1,nnold ! loop over the nodes
!!$     if(geometry%fix(i).ne.1.and.network%receiver(i).ne.0)then
!!$        flag=0
!!$        do k=geometry%nb(i),1,-1 ! loop over the neighbors
!!$           if(network%receiver(i).eq.geometry%nn(k,i))flag=1  ! receiver is a neighbor
!!$        enddo
!!$        if(flag.eq.0)then !reciever is not a neighbor
!!$           j=network%receiver(i)
!!$           nadd_river_no_connection=nadd_river_no_connection+1
!!$           geometry%nnode=geometry%nnode+1
!!$           if (geometry%nnode.gt.geometry%nnode_max) then
!!$              print*,'problem while adding nodes'
!!$              STOP 'too many nodes added. Increase nnode_max'
!!$           endif
!!$           geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0
!!$           geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0
!!$           geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
!!$          ! print*, 'adding node', geometry%nnode, 'to prevent channel between non-neighboring node in:'
!!$          ! print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
!!$           ! Assign physical properties to new node as average of its endmembers
!!$           geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
!!$           geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
!!$           geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
!!$           geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
!!$           network%receiver(geometry%nnode)=j
!!$           network%receiver(i)=geometry%nnode
!!$           geometry%fix(geometry%nnode)=0
!!$           geometry%surface(geometry%nnode)=0.0d0
!!$           geometry%discharge(geometry%nnode)=0.0d0
!!$           geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
!!$           geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
!!$           geometry%sediment_flux(geometry%nnode)=0.0d0
!!$        endif
!!$     endif
!!$  enddo

  !print*,  'nadd_river_no_connection', nadd_river_no_connection

  ! This section is to add nodes on channels. It checks whether nodes on the drainage network
  ! are too far apart or not. I adds nodes when the distance between two nodes is
  ! greater than a critical distance "lmax" (chosen in initialize_parameters.f90) - SDW

  nadd_on_channel=0
  do i=1,nnold ! loop over the nodes
     if (network%receiver(i).ne.0) then ! only for connections that are part of the drainage network
        j=network%receiver(i)
        l=dsqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length of segment betwenn i and j
        if (l.gt.params%lmax) then ! add a node on the connection
           nadd_on_channel=nadd_on_channel+1
           geometry%nnode=geometry%nnode+1
           if (geometry%nnode.gt.geometry%nnode_max) then
              print*,'problem while adding nodes'
              STOP 'too many nodes added. Increase nnode_max'
           endif
           ! Add a random perturbation in location of new node equal to 50% of the max length between nodes
           call random_number(xx)
           geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0 &
                +(2.d0*xx-1.d0)*0.3d0*params%lmax 
           call random_number(xx)
           geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0 &
                +(2.d0*xx-1.d0)*0.3d0*params%lmax
           geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
           !print*, 'adding node', geometry%nnode, 'in a too long channel in:'
           !print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
           !note every physical properties will have to be updated at each new
           !node....
           ! Assign physical properties to new node as average of its endmembers
           ! added 15.6.2010 sdw  Not sure if these are all the node-local properties
           geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
           geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
           geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
           geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
           network%receiver(geometry%nnode)=j
           network%receiver(i)=geometry%nnode
           geometry%fix(geometry%nnode)=0
           geometry%surface(geometry%nnode)=0.0d0
           geometry%discharge(geometry%nnode)=0.0d0
           geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
           geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0           
           geometry%sediment_flux(geometry%nnode)=0.0d0 
        endif
     endif
  enddo

!print*, 'nadd_on_channel',nadd_on_channel


  ! add node on non-channel divides that are longer than a given distance
  ! add node to lower existing node
  ! note that this loop does not take into acount boundary nodes.

  nadd_between_channel = 0
  ! the if statement below is needed to stop potentially crossing rivers -SDW
  ! (LG - maybe because the neighbor list was not updated)
  if(nadd_river_no_connection.eq.0.and.nadd_on_channel.eq.0)then
     do i=1,nnold ! loop over the original nodes
        do k=geometry%nb(i),1,-1 ! loop over the neighbors
           j=geometry%nn(k,i)
           if (geometry%fix(i).ne.1.or.geometry%fix(j).ne.1) then !when both are boundary nodes treated with nadd_boundary
              if(network%receiver(i).ne.j.and.network%receiver(j).ne.i) then ! check non-channel
                 if(geometry%z(j).gt.geometry%z(i))then  ! only continue if j is higher than i
                    l=dsqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length of segment betwenn i and j
                    if (l.gt.params%ldivmax) then ! add node on the connection
                       nadd_between_channel=nadd_between_channel+1                          
                       geometry%nnode=geometry%nnode+1
                       !print*, 'adding node', geometry%nnode,'between',i,j
                       if (geometry%nnode.gt.geometry%nnode_max) then
                          print*,'problem while adding nodes'
                          STOP 'too many nodes added. Increase nnode_max'
                       endif
                       call random_number(xx)
                       geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0 &
                            +(2.d0*xx-1.d0)*0.3d0*params%lmax 
                       call random_number(xx)
                       geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0 &
                            +(2.d0*xx-1.d0)*0.3d0*params%lmax
                       geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
                       geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
                       geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                       geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                       geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
                       network%ndon(i)=network%ndon(i)+1
                       network%ndon(geometry%nnode)=0
                       network%receiver(geometry%nnode)=i
                       geometry%fix(geometry%nnode)=0
                       geometry%surface(geometry%nnode)=0.0d0
                       geometry%discharge(geometry%nnode)=0.0d0
                       geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                       geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0           
                       geometry%sediment_flux(geometry%nnode)=0.0d0 
                    endif
                 endif
              endif
           endif
        enddo
     enddo
  endif


  !print*, 'nadd_between_channel',nadd_between_channel

  ! This section is to add boundary nodes when the boundaries are advected. 
  ! It adds nodes when the distance between two boundary nodes is
  ! greater than a critical distance "lmax" (chosen in initialize_parameters.f90)

  ! for some reason the nb array for boundary nodes is not trustable. Need to generate a new list.
  ! Each node will point to its neighbor in anticlockwise direction.

  !divided into four section according tobottom, right, top and left boundaries

  ! try again with the neighbors list 1/4/11
  nadd_boundary=0
  if (params%move_points) then
     do i=1,nnold
        if (geometry%fix(i).eq.1) then
           if (geometry%y(i).eq.ymin) then
              do k = 1,geometry%nb(i)
                 j = geometry%nn(k,i)
                 if (geometry%fix(j).eq.1.and.geometry%y(j).eq.ymin.and.geometry%x(j).gt.geometry%x(i)) then !this is the closesent neighbor to the right
                    l = dsqrt((geometry%x(i)-geometry%x(j))**2)
                    if (l.gt.params%lmax) then
                       nadd_boundary=nadd_boundary+1
                       geometry%nnode=geometry%nnode+1
                       if (geometry%nnode.gt.geometry%nnode_max) then
                          print*,'problem while adding nodes'
                          STOP 'too many nodes added. Increase nnode_max'
                       endif
                       geometry%y(geometry%nnode)=geometry%y(i)
                       call random_number(xx)
                       geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0 &
                            +(2.d0*xx-1.d0)*0.05d0*params%lmax
                       geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0 !should be zero
                       !print*, 'adding node', geometry%nnode, 'because distance between boundary nodes is too large in:'
                       !print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
                       geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.
                       geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                       geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                       geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
                       geometry%fix(geometry%nnode)=1
                       geometry%surface(geometry%nnode)=0.0d0
                       geometry%discharge(geometry%nnode)=0.0d0
                       geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                       geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
                       geometry%sediment_flux(geometry%nnode)=0.0d0
                       network%receiver(geometry%nnode)=0                      
                    endif
                    exit !because the neighbor to the right was found.
                 endif
              enddo
           endif
           if (geometry%x(i).eq.xmax) then
              do k = 1,geometry%nb(i)
                 j = geometry%nn(k,i)
                 if (geometry%fix(j).eq.1.and.geometry%x(j).eq.xmax.and.geometry%y(j).gt.geometry%y(i)) then !this is the closesent neighbor to the top
                    l = dsqrt((geometry%y(i)-geometry%y(j))**2)
                    if (l.gt.params%lmax) then
                       nadd_boundary=nadd_boundary+1
                       geometry%nnode=geometry%nnode+1
                       if (geometry%nnode.gt.geometry%nnode_max) then
                          print*,'problem while adding nodes'
                          STOP 'too many nodes added. Increase nnode_max'
                       endif
                       geometry%x(geometry%nnode)=geometry%x(i)
                       call random_number(xx)              
                       geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0 &
                            +(2.d0*xx-1.d0)*0.05d0*params%lmax
                       geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
                       !print*, 'adding node', geometry%nnode, 'because distance between boundary nodes is too large'
                       !print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
                       geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
                       geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                       geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                       geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
                       geometry%fix(geometry%nnode)=1
                       geometry%surface(geometry%nnode)=0.0d0
                       geometry%discharge(geometry%nnode)=0.0d0
                       geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                       geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
                       geometry%sediment_flux(geometry%nnode)=0.0d0
                       network%receiver(geometry%nnode)=0
                    endif
                    exit !because the neighbor to the right was found.
                 endif
              enddo
           endif
           if (geometry%y(i).eq.ymax) then
              do k = 1,geometry%nb(i)
                 j = geometry%nn(k,i)
                 if (geometry%fix(j).eq.1.and.geometry%y(j).eq.ymax.and.geometry%x(j).lt.geometry%x(i)) then !this is the closesent neighbor to the left
                    l = dsqrt((geometry%x(i)-geometry%x(j))**2)
                    if (l.gt.params%lmax) then
                       nadd_boundary=nadd_boundary+1
                       geometry%nnode=geometry%nnode+1
                       if (geometry%nnode.gt.geometry%nnode_max) then
                          print*,'problem while adding nodes'
                          STOP 'too many nodes added. Increase nnode_max'
                       endif
                       geometry%y(geometry%nnode)=geometry%y(i)
                       call random_number(xx)
                       geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0 &
                            +(2.d0*xx-1.d0)*0.05d0*params%lmax
                       geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0 !should be zero
                       !print*, 'adding node', geometry%nnode, 'because distance between boundary nodes is too large'
                       !print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
                       geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
                       geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                       geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                       geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
                       geometry%fix(geometry%nnode)=1
                       geometry%surface(geometry%nnode)=0.0d0
                       geometry%discharge(geometry%nnode)=0.0d0
                       geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                       geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
                       geometry%sediment_flux(geometry%nnode)=0.0d0
                       network%receiver(geometry%nnode)=0                   
                    endif
                    exit !because the neighbor to the right was found.
                 endif
              enddo
           endif
           if (geometry%x(i).eq.xmin) then
              do k = 1,geometry%nb(i)
                 j = geometry%nn(k,i)
                 if (geometry%fix(j).eq.1.and.geometry%x(j).eq.xmin.and.geometry%y(j).lt.geometry%y(i)) then !this is the closesent neighbor to the bottom
                    l = dsqrt((geometry%y(i)-geometry%y(j))**2)
                    if (l.gt.params%lmax) then
                       nadd_boundary=nadd_boundary+1
                       geometry%nnode=geometry%nnode+1
                       if (geometry%nnode.gt.geometry%nnode_max) then
                          print*,'problem while adding nodes'
                          STOP 'too many nodes added. Increase nnode_max'
                       endif
                       geometry%x(geometry%nnode)=geometry%x(i)
                       call random_number(xx)              
                       geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0 &
                            +(2.d0*xx-1.d0)*0.05d0*params%lmax
                       geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
                       !print*, 'adding node', geometry%nnode, 'because distance between boundary nodes is too large'
                       !print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
                       geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
                       geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                       geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                       geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
                       geometry%fix(geometry%nnode)=1
                       geometry%surface(geometry%nnode)=0.0d0
                       geometry%discharge(geometry%nnode)=0.0d0
                       geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                       geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
                       geometry%sediment_flux(geometry%nnode)=0.0d0
                       network%receiver(geometry%nnode)=0
                    endif
                    exit !because the neighbor to the right was found.
                 endif
              enddo
           endif
        endif
     enddo
  endif


  !print*, 'nadd_boundary',nadd_boundary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Second stage: take care of all the cases where a new node should removed !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! this section remove nodes based on minimum surface area in headwater
  nrem_small_area=0
!!$  do i=1, nnold
!!$     !because we remove nodes, need to check if we are still within the bounds of ald nodes. 
!!$     if(i.le.geometry%nnode-nadd_on_channel-nadd_between_channel-nadd_river_no_connection-nadd_boundary)then 
!!$        if(geometry%fix(i).ne.1.and.network%ndon(i).eq.0.and.geometry%surface(i).lt.params%amin)then
!!$           nrem_small_area=nrem_small_area+1
!!$           geometry%nnode=geometry%nnode-1  
!!$           !print*, 'removing node', i, 'because its a leaf with too small drainage area'
!!$           !print*, geometry%x(i), geometry%y(i), geometry%z(i)
!!$           do j=1,i-1
!!$              if (network%receiver(j).gt.i) network%receiver(j)=network%receiver(j)-1
!!$           enddo
!!$           if(i.le.geometry%nnode)then
!!$              do j=i,geometry%nnode
!!$                 geometry%x(j)=geometry%x(j+1)
!!$                 geometry%y(j)=geometry%y(j+1)
!!$                 geometry%z(j)=geometry%z(j+1)
!!$                 geometry%u(j)=geometry%u(j+1)
!!$                 geometry%v(j)=geometry%v(j+1)
!!$                 geometry%w(j)=geometry%w(j+1)
!!$                 geometry%fix(j)=geometry%fix(j+1)
!!$                 geometry%surface(j)=geometry%surface(j+1)
!!$                 geometry%discharge(j)=geometry%discharge(j+1)
!!$                 geometry%precipitation(j)=geometry%precipitation(j+1)
!!$                 geometry%k(j)=geometry%k(j+1)
!!$                 geometry%nb(j)=geometry%nb(j+1)
!!$                 geometry%erosion_rate(j)=geometry%erosion_rate(j+1)
!!$                 geometry%sediment_flux(j)=geometry%sediment_flux(j+1)
!!$                 network%ndon(j)=network%ndon(j+1)
!!$                 if (network%receiver(j+1).gt.i) then
!!$                    network%receiver(j)=network%receiver(j+1)-1
!!$                 else
!!$                    network%receiver(j)=network%receiver(j+1)
!!$                 endif
!!$              enddo
!!$           endif
!!$        endif
!!$     endif
!!$  enddo


!!! the next remove section operates also on nodes that were just now added. 
!!!For that reason the donor receiver array needs to be updated
  if (nrem_small_area.gt.0.or.nadd_on_channel.gt.0.or.nadd_between_channel.gt.0&
       &.or.nadd_river_no_connection.gt.0.or.nadd_boundary.gt.0)then
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


  ! This section removes nodes that are too close to the boundary due to advection
  ! Both internal modes and boundary nodes

  nrem_close_to_boundary = 0
  if (params%move_points) then 
     min_dist =  params%max_adv*params%deltat !if a distance of a node to the boundary < min_dist-->remove it
     do i=1, geometry%nnode !also new nodes that were just added
554     if (i.le.geometry%nnode) then ! i can grow above geometry%nnode because geometry%nnode can become smaller duringt this loop
           if(geometry%fix(i).ne.1.or.& !internal nodes
                &geometry%y(i).eq.ymin.and.geometry%x(i).ne.xmin.and.geometry%x(i).ne.xmax.or.& !boundary nodes but not coreners
                &geometry%y(i).eq.ymax.and.geometry%x(i).ne.xmin.and.geometry%x(i).ne.xmax.or.& !boundary nodes but not coreners
                &geometry%x(i).eq.xmin.and.geometry%y(i).ne.ymin.and.geometry%y(i).ne.ymax.or.& !boundary nodes but not coreners
                &geometry%x(i).eq.xmax.and.geometry%y(i).ne.ymin.and.geometry%y(i).ne.ymax) then!boundary nodes but not coreners
              if (abs(geometry%x(i)-xmin).lt.min_dist.and.geometry%x(i).ne.xmin.or.& !close to left
                   &abs(geometry%x(i)-xmax).lt.min_dist.and.geometry%x(i).ne.xmax.or.& !close to right
                   &abs(geometry%y(i)-ymin).lt.min_dist.and.geometry%y(i).ne.ymin.or.& !close to bottom
                   &abs(geometry%y(i)-ymax).lt.min_dist.and.geometry%y(i).ne.ymax) then !close to top
                 !print*, 'removing node', i, 'because its too close to a boundary'
                 !print*, geometry%x(i), geometry%y(i), geometry%z(i) 
                 !print*, 'updated number of nodes is',geometry%nnode
                 nrem_close_to_boundary = nrem_close_to_boundary +1
                 geometry%nnode=geometry%nnode-1                 
                 ! first reorganize network
                 k = network%receiver(i)
                 if (k.ne.0) then !i is not a lake
                    do j = 1,network%ndon(i)
                       h = network%donors(j,i) ! h is a donor of i
                       network%receiver(h) = k ! now h is a donor of k
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
                 if(i.le.geometry%nnode)then !not the last node
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
                 !need to continuously update donors and receivers
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
                 go to 554 !without increasing i becuase the node array is smaller.
              endif
           endif
        endif
     enddo
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! Short debuging section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


  
  call time_out ('add_remove_nodes')

  return 
end subroutine add_remove_nodes
