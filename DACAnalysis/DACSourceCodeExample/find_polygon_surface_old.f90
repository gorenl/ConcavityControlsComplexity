subroutine find_polygon_surface (geometry,network,params,delaunay)


  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (parm) params
  type (del) delaunay
  integer i,j,k,h,tri, i1, i2, i3
  double precision mid_x,mid_y,area, tri_area, total_area, total_area_2
 


  call time_in ('find_polygon_surface')
 

 
  total_area = 0.

! do some preprocessing so that the calculation be faster
  
  do i=1,geometry%nnode
     ! a node is in its polygon
     area = 0.
     do k = 1,geometry%nb(i)
        j = geometry%nn(k,i)
        if (i.eq.network%receiver(j)) then !i and j are part of the channel network
           ! find the centroids of the two triangles that share i and j 
           do h = 1,delaunay%ntriangles
              if (delaunay%icon(1,h).eq.i.or.delaunay%icon(2,h).eq.i.or.delaunay%icon(3,h).eq.i) then
                 if (delaunay%icon(1,h).eq.j.or.delaunay%icon(2,h).eq.j.or.delaunay%icon(3,h).eq.j) then
                    ! h is one of the triangles
                     area = area + tri_area(geometry%x(i),geometry%y(i),geometry%x(j),geometry%y(j),&
                          &delaunay%centers(1,h),delaunay%centers(2,h))
                     ! to find the other triangle go over the neighbors of h
                     tri = delaunay%neighbours(1,h)
                     if (delaunay%icon(1,tri).eq.i.or.delaunay%icon(2,tri).eq.i.or.delaunay%icon(3,tri).eq.i) then
                        if (delaunay%icon(1,tri).eq.j.or.delaunay%icon(2,tri).eq.j.or.delaunay%icon(3,tri).eq.j) then
                           area = area + tri_area(geometry%x(i),geometry%y(i),geometry%x(j),geometry%y(j),&
                                &delaunay%centers(1,tri),delaunay%centers(2,tri))
                           exit
                        endif
                     endif
                     tri = delaunay%neighbours(2,h)
                     if (delaunay%icon(1,tri).eq.i.or.delaunay%icon(2,tri).eq.i.or.delaunay%icon(3,tri).eq.i) then
                        if (delaunay%icon(1,tri).eq.j.or.delaunay%icon(2,tri).eq.j.or.delaunay%icon(3,tri).eq.j) then
                           area = area + tri_area(geometry%x(i),geometry%y(i),geometry%x(j),geometry%y(j),&
                                &delaunay%centers(1,tri),delaunay%centers(2,tri))
                           exit
                        endif
                     endif
                     tri = delaunay%neighbours(3,h)
                     if (delaunay%icon(1,tri).eq.i.or.delaunay%icon(2,tri).eq.i.or.delaunay%icon(3,tri).eq.i) then
                        if (delaunay%icon(1,tri).eq.j.or.delaunay%icon(2,tri).eq.j.or.delaunay%icon(3,tri).eq.j) then
                           area = area + tri_area(geometry%x(i),geometry%y(i),geometry%x(j),geometry%y(j),&
                                &delaunay%centers(1,tri),delaunay%centers(2,tri))
                           exit
                        endif
                     endif
                     print*, 'I should not be here - only one neighboring triangle'
                  endif
               endif
            enddo           
        else
           !if i and j are not part of the channel network, must be divide in between or they are both at the base level
           if (j.ne.network%receiver(i)) then 
              if (geometry%fix(i).eq.1.and.geometry%fix(j).eq.1) then !both on the base level no divide in between
                 !find the common triangle
                 mid_x = (geometry%x(i) + geometry%x(j))/2.
                 mid_y = (geometry%y(i) + geometry%y(j))/2.
                 do h = 1,delaunay%ntriangles
                    if (delaunay%icon(1,h).eq.i.or.delaunay%icon(2,h).eq.i.or.delaunay%icon(3,h).eq.i) then
                       if (delaunay%icon(1,h).eq.j.or.delaunay%icon(2,h).eq.j.or.delaunay%icon(3,h).eq.j) then
                          area = area +  tri_area(geometry%x(i),geometry%y(i),mid_x,mid_y,&
                               &delaunay%centers(1,h),delaunay%centers(2,h))
                          tri = delaunay%neighbours(1,h)
                          if (delaunay%icon(1,tri).eq.i.or.delaunay%icon(2,tri).eq.i.or.delaunay%icon(3,tri).eq.i) then
                             if (delaunay%icon(1,tri).eq.j.or.delaunay%icon(2,tri).eq.j.or.delaunay%icon(3,tri).eq.j) then
                                area = area + tri_area(geometry%x(i),geometry%y(i),mid_x,mid_y,&
                                     &delaunay%centers(1,tri),delaunay%centers(2,tri))
                                exit
                             endif
                          endif
                          tri = delaunay%neighbours(2,h)
                          if (delaunay%icon(1,tri).eq.i.or.delaunay%icon(2,tri).eq.i.or.delaunay%icon(3,tri).eq.i) then
                             if (delaunay%icon(1,tri).eq.j.or.delaunay%icon(2,tri).eq.j.or.delaunay%icon(3,tri).eq.j) then
                                area = area + tri_area(geometry%x(i),geometry%y(i),mid_x,mid_y,&
                                     &delaunay%centers(1,tri),delaunay%centers(2,tri))
                                exit
                             endif
                          endif
                          tri = delaunay%neighbours(3,h)
                          if (delaunay%icon(1,tri).eq.i.or.delaunay%icon(2,tri).eq.i.or.delaunay%icon(3,tri).eq.i) then
                             if (delaunay%icon(1,tri).eq.j.or.delaunay%icon(2,tri).eq.j.or.delaunay%icon(3,tri).eq.j) then
                                area = area + tri_area(geometry%x(i),geometry%y(i),mid_x,mid_y,&
                                     &delaunay%centers(1,tri),delaunay%centers(2,tri))
                                exit
                             endif
                          endif
                          exit ! it is possible that the is only one common triangle
                       endif
                    endif
                 enddo
              else !must be a divide between them, find it
                 do h = 1,geometry%ndivide
                    if (geometry%nndivnode(h,1).eq.i.or.geometry%nndivnode(h,2).eq.i) then
                       if (geometry%nndivnode(h,1).eq.j.or.geometry%nndivnode(h,2).eq.j) then
                          !h is a divide between i and j                     
                          !find the coordinates of the centroids of the two adjecent triangles 
                          tri = geometry%nndivtri(h,1)
                          area = area + tri_area(geometry%x(i),geometry%y(i),geometry%xdiv(h),geometry%ydiv(h),&
                               &delaunay%centers(1,tri),delaunay%centers(2,tri))                     
                          tri = geometry%nndivtri(h,2)
                          area = area + tri_area(geometry%x(i),geometry%y(i),geometry%xdiv(h),geometry%ydiv(h),&
                               &delaunay%centers(1,tri),delaunay%centers(2,tri))                      
                          exit !because the divide was found
                       endif
                    endif
                 enddo
              endif
           endif
        endif
     enddo 
     geometry%surface(i) = area
     total_area = total_area + area
     !print*, 'area is', area
  enddo
  call time_out ('find_polygon_surface')
!  print*, 'total area is', total_area
!!$  total_area_2 = 0.
!!$  do i=1,delaunay%ntriangles
!!$     i1 = delaunay%icon(1,i)
!!$     i2 = delaunay%icon(2,i)
!!$     i3 = delaunay%icon(3,i)
!!$     total_area_2 = total_area_2 + tri_area(geometry%x(i1),geometry%y(i1),geometry%x(i2),geometry%y(i2),geometry%x(i3),geometry%y(i3))
!!$  enddo
!!$  print*, 'area according to triangulation is', total_area_2
 

  return

end subroutine find_polygon_surface


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function tri_area (x1,y1,x2,y2,x3,y3)
  implicit none
  double precision x1,y1,x2,y2,x3,y3,l1,l2,l3,s
  
  l1 = sqrt((x1-x2)**2 + (y1-y2)**2)
  l2 = sqrt((x2-x3)**2 + (y2-y3)**2)
  l3 = sqrt((x3-x1)**2 + (y3-y1)**2)
 ! s = (l1+l2+l3)/2.
  !tri_area =sqrt(s*(s-l1)*(s-l2)*(s-l3))  
  tri_area = 0.5*abs(x1*(y3-y2)+x2*(y1-y3) + x3*(y2-y1))
  
  return 
end function tri_area
