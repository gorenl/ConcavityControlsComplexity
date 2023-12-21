MODULE definitions

   type parm
      integer n,istep,freq,num_restart
      double precision deltat,time,tfinal,h,ka,xc,tanthetac,hmn,m
      double precision rainfall_available,rainfall_minimum,rainfall_height
      double precision lmax,diffusivity,amin,min_erosion_rate,uplift_scalar1
      double precision uplift_scalar2,k_scalar1,ldivmax, max_adv, min_tan_head_slope
      logical plot_triangles,plot_receiver,plot_donors,plot_no_receiver
      logical plot_no_donors,write_discharge,write_stack_order,plot_catchment
      logical plot_height,plot_precipitation
      logical move_points,capture,divide,small_divide,diffusion
      logical add_nodes
      logical read_restart, ascii, transient_divide
      integer num_bins, sample_per_bin
   end type parm

   type geom
      double precision,dimension(:),pointer::x,y,z,u,v,w 
      double precision,dimension(:),pointer::xdiv,ydiv,zdiv
      double precision,dimension(:),pointer::surface,discharge,precipitation
      double precision,dimension(:),pointer::k,erosion_rate,sediment_flux
      integer,dimension(:),pointer::strahler
      double precision,dimension(:,:),pointer::surface_share
      integer,dimension(:),pointer::fix,nb,catchment
      integer,dimension(:,:),pointer::nn
      integer,dimension(:,:),pointer::nndivtri !for each divide two triangles  
      integer,dimension(:,:),pointer::nndivnode !for each divide two nodes
      double precision, dimension(:,:),pointer::erosion_rate_history
      integer nnode,nnode_max,nnmax,ndivide,ndivide_max,ncapture
   end type geom

   type del
      integer ntriangles
      integer,dimension(:,:),pointer::icon,neighbours,numdivides
      double precision,dimension(:,:),pointer::centers
   end type del

   type netw
      integer nnode,ndonmax,nlake
      integer,dimension(:),pointer::receiver,ndon,lakes,lakes_catch
      integer,dimension(:,:),pointer::donors
   end type netw

   type stck
      integer nnode
      integer,dimension(:),allocatable::order
   end type stck

   type timr
      sequence
      integer ntimer
      character*256 name(1024)
      real time_spent(1024)
      character*256 namein
      real timein
      real time_start,time_stop
   end type timr

end MODULE definitions
