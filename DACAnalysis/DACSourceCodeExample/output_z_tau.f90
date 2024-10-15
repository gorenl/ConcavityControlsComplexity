subroutine output_z_tau (params,geometry,network,stack)


  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (stck) stack
  type (parm) params
  integer i,j,k,icount,ii
  double precision l,ks_K, slope
  character cs*8
  double precision,dimension(:),allocatable::contributing_area, tau, faketau
  

  call time_in ('output_z_tau')


  icount=params%istep

  write(cs,'(I8)') icount
  if (icount.lt.10)      cs(1:7)='0000000'
  if (icount.lt.100)     cs(1:6)='000000'
  if (icount.lt.1000)    cs(1:5)='00000'
  if (icount.lt.10000)   cs(1:4)='0000'
  if (icount.lt.100000)  cs(1:3)='000'
  if (icount.lt.1000000)  cs(1:2)='00'
  if (icount.lt.10000000)  cs(1:1)='0'


 




  allocate(contributing_area(geometry%nnode), tau(geometry%nnode),faketau(geometry%nnode))

  contributing_area=geometry%surface
  do i=stack%nnode,1,-1
     j=stack%order(i)
     k=network%receiver(j)
     if (k.ne.0) then
        contributing_area(k) = contributing_area(k) + contributing_area(j)
     endif
  enddo
  ks_K=params%k_scalar1*params%rainfall_height**params%m
  tau=0.
  faketau=0.
  do i=1,stack%nnode
     j=stack%order(i)
     if (network%receiver(j).eq.0) then
        tau(j)=0.
        faketau(j)=0.
     else                 
        k=network%receiver(j)
        l = dsqrt((geometry%x(k)-geometry%x(j))**2.d0 + (geometry%y(k)-geometry%y(j))**2.d0)
        slope = (geometry%z(j)-geometry%z(k))/l
        tau(j)=tau(k)+(1./(ks_K*contributing_area(j)**params%m*(slope**(params%n - 1.))))*l
        faketau(j)=faketau(k)+(1./(ks_K*contributing_area(j)**params%m))*l
     endif
  enddo
  open (75,file='ASCII/tau_z'//cs,status='unknown')
  do i=1,geometry%nnode
     write(75,'(f8.3,f17.5,f17.5,f17.5,i5)') geometry%z(i), tau(i), faketau(i), contributing_area(i), geometry%catchment(i)
  enddo
  close (75)
  open (78,file='ASCII/no_channel_connection'//cs,status='unknown')
  do i=1,geometry%nnode ! loop over the nodes
     do k=geometry%nb(i),1,-1 ! loop over all neighbours
        j=geometry%nn(k,i)
        if(network%receiver(j).ne.i.and.network%receiver(i).ne.j) then ! only for connections that are not part of the drainage network
           write(78,'(i5,i5)') i, j
        endif
     enddo
  enddo
  close(78)
        

  

  
  deallocate(contributing_area, tau, faketau)
  call time_out ('output_z_tau')

  return
end subroutine output_z_tau
