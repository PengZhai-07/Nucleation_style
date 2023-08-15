using FortranFiles
using OffsetArrays
using Parameters
using Printf

module spec_grid

#  use stdio, ONLY : IO_abort
#  use constants
#  use fem_grid
#  use bnd_grid

#  implicit none
  # private           # function is not available for execution outside module

#-----------------------------------------------------------------------
#-- Spectral element grid type
#
#------ Topology
#     nelem        = total number of elements
#     ngll         = number of Gauss-Lobato-Legendre nodes per unit segment
#     npoin        = total number of GLL nodes in the mesh
#     ibool        = local to global numbering for the GLL mesh
#                    (igll,jgll,element) -> bulk node index
#     coord        = coordinates of the GLL nodes
#
#------ Working data
#     shape        = control nodes to GLL nodes shape functions
#     dshape       = derivatives of shape functions
#     hprime       = lagrange derivatives on the unit square
#     hTprime      = transpose of hprime
#     wgll         = 1D GLL integration weights
#     wgll2        = 2D GLL integration weights
#     xgll         = GLL points on the unit segment
#
#------ Other
#     tag          = element -> domain tags
#     fmax         = Highest frequency to be resolved by the grid
#     W            = The seismogenic width. Infinity means 2D problem 
#                    and finite W means 2.5D problem (use for elastic material)

#
#
# NOTE: After initialization, the FEM mesh database is only used for plotting
#       purposes, so if no plots are needed you can save them on disk and 
#       clean it from memory.
#
# NOTE: COORD is used in the solver only for computation of the incident
#       wavefield, the boundary coordinates could be stored instead in the 
#       bc structure (boundary_condition).


  #type sem_grid_type
  #    integer :: ngll=0,nelem=0,npoin=0

  mutable struct sem_grid_type
    ngll::Int32 = 0
    nelem::Int32 = 0
    npoin::Int32 = 0
    fem_grid_type(fem)      # one kind of struct
    fmax :: Float64=0e0
    W :: Float64= typemax(Int128)   
    coord = Array{Float64,2}       
    hprime = Array{Float64,2} 
    hTprime = Array{Float64,2} 
    wgll= Vector{Float64,1} 
    xgll= Vector{Float64,1}
    wgll2= Array{Float64,2} 
    shape = Array{Float64,3}
    dshape = Array{Float64,4}
    ibool = Array{Int64, 3}
    tag = Vector{Int64, 1}
    bnd_grid_type(bounds)  
  end
end

# #-----------------------------------------------------------------------

#   begin # interface SE_Jacobian
#     module procedure SE_Jacobian_e, SE_Jacobian_eij
#   end
#       function NAMELESSPROGRAM()

#   begin # interface SE_InverseJacobian
#     module procedure SE_InverseJacobian_e, SE_InverseJacobian_eij
#   end

#   begin # interface SE_VolumeWeights
#     module procedure SE_VolumeWeights_e, SE_VolumeWeights_eij
#   end

#   begin # interface SE_node_belongs_to
#     module procedure SE_node_belongs_to_1, SE_node_belongs_to_2
#   end

#   # (i,j) GLL for element corner nodes
#   #icorner(1) = 1    ; jcorner(1) = 1
#   #icorner(2) = ngll ; jcorner(2) = 1
#   #icorner(3) = ngll ; jcorner(3) = ngll
#   #icorner(4) = 1    ; jcorner(4) = ngll


#   public :: sem_grid_type ,
#          SE_init ,
#          SE_init_interpol ,
#          SE_find_nearest_node ,
#          SE_node_belongs_to,
#          SE_find_point ,
#          SE_get_edge_nodes ,
#          SE_inquire ,
#          SE_isFlat,
#          SE_firstElementTagged,
#          SE_VolumeWeights ,
#          SE_InverseJacobian ,
#          SE_Jacobian ,
#          SE_elem_coord,
#          BC_inquire ,
#          BC_tag_exists,
#          BC_get_normal_and_weights,
#            edge_D, edge_R, edge_U, edge_L

# # ======================================================================
# #
#   function SE_isFlat(se)
#   type(sem_grid_type), intent(in) :: se
#   SE_isFlat = se.fem%flat
#   end

# # ======================================================================
# #
#   function SE_init(se)

# #  use stdio, ONLY : IO_new_unit
#   type(sem_grid_type), intent(inout) :: se
# #  integer :: ounit

#   SE_init_gll(se) # GLL quadrature points, weights 
#                             # and lagrange coefficients,
#                             # shape functions and their derivatives
#   SE_init_numbering(se) # generate the global node numbering
#   SE_init_coord(se) # get the coordinates of the GLL nodes,
#                               # and export the grid to a file
#   SE_BcTopoInit(se) # initialize boundaries generic data

#   # export grid parameters
#   ounit = IO_new_unit()
#   OPEN(ounit, file = "grid_sem2d.hdr", status = "replace")
#   println(ounit, "NELEM  NPGEO  NGNOD  NPOIN  NGLL")
#   println(ounit, se.nelem, FE_GetNbNodes(se.fem), FE_GetNodesPerElement(se.fem),
#                  se.npoin, se.ngll)
#   CLOSE(ounit)

#   end

# # ======================================================================
# #
# # get coordinates and weights of the Gauss-Lobatto-Legendre points
# # and derivatives of Lagrange polynomials  H_ij = h'_i(xgll(j))
#   function SE_init_gll(se)

# #  use gll, ONLY : get_GLL_info, print_GLL

#   type(sem_grid_type), intent(inout) :: se
# #  DOUBLEPRECISION :: xi,eta
# #  integer :: ngll,i,j, ngnod

#   ngll = se.ngll

#   allocate(se.xgll(ngll))
#   allocate(se.wgll(ngll))
#   allocate(se.hprime(ngll, ngll))
#   allocate(se.hTprime(ngll, ngll))
#   allocate(se.wgll2(ngll, ngll))

#   get_GLL_info(ngll, se.xgll, se.wgll, se.hprime)
#   print_GLL(ngll, se.xgll, se.wgll, se.hprime)
#   se.hTprime = TRANSPOSE(se.hprime)

#  # wgll2(i,j) = wgll(i) * wgll(j)
#   for j = 1:ngll
#     se.wgll2(:, j) = se.wgll * se.wgll(j)
#   end

# #-----------------------------------------------------------------------
# #-- set up the shape functions and their local derivatives
#   ngnod = FE_GetNodesPerElement(se.fem)

#   allocate(se.shape(ngnod, ngll, ngll))
#   allocate(se.dshape(ngnod, NDIME, ngll, ngll))

#   for j = 1:ngll
#     eta = se.xgll(j)
#     for i = 1:ngll
#       xi = se.xgll(i)
#       se.shape(:, i, j) = FE_getshape(xi, eta, ngnod)
#       se.dshape(:, :, i, j) = FE_getdershape(xi, eta, ngnod)
#     end
#   end

#   end



# # ======================================================================
# #
# #-- generate the global numbering
# #
#   function SE_init_numbering(se)

# #  use memory_info
# #  use echo, ONLY : echo_init,iout,fmt1,fmtok
# #  use stdio, ONLY : IO_new_unit

#   type(sem_grid_type), intent(inout) :: se

# #  integer, POINTER :: ibool[:,:,:]
# #  integer, dimension(se.ngll,4) :: iedg,jedg,iedgR,jedgR
# #  integer, dimension(4) :: ivtx,jvtx
# #  integer :: i,j,k,e,n,ii,jj,ee,nn,ngll,npoin,iol,ounit
# #  integer, POINTER :: ees[:],nns[:]

# #-----------------------------------------------------------------------

#   if echo_init
#     println(iout, )
#     @printf(iout, "%s\n", " S p e c t r a l   e l e m e n t s   g r i d")
#     @printf(iout, "%s\n", " ===========================================")
#     println(iout, )
#   end

#   ngll  = se.ngll
#   se.nelem = FE_GetNbElements(se.fem)
#   se.tag =Ref().fem%tag

# # global numbering table
#   allocate(se.ibool[ngll,ngll,se.nelem])
#   storearray("ibool", SIZE(se.ibool), iinteg)
#   ibool =Ref().ibool
#   ibool = 0

# #---- start numbering
#   if (echo_init) println(iout, "no", "Numbering GLL points") end
#   npoin = 0

#  # GLL index (i,j) of edge nodes, counterclockwise
#   for n = 1:4
#     SE_inquire(se, edge = n, itab = iedg[:,n], jtab = jedg[:,n])
#   end
#  # reverse order for matching edges
#   iedgR = iedg[ngll:1:-1,:]
#   jedgR = jedg[ngll:1:-1,:]

#  # GLL index (i,j) for vertex nodes
#   ivtx[1] = 1;    jvtx[1] = 1
#   ivtx[2] = ngll; jvtx[2] = 1
#   ivtx[3] = ngll; jvtx[3] = ngll
#   ivtx[4] = 1;    jvtx[4] = ngll

#   for e = 1:se.nelem


#    #-- interior nodes are unique
#     for j = 2:ngll-1
#     for i = 2:ngll-1
#       npoin = npoin + 1
#       ibool[i,j,e] = npoin
#     end
#     end

#    #-- edge nodes
#     for n = 1:4
#       if (ibool[iedg[2,n],jedg[2,n],e]>0) continue end # skip if already processed
#       FE_GetEdgeConn(se.fem, e, n, ee, nn)
#       for k = 2:ngll-1
#         npoin = npoin+1
#         ibool[iedg[k,n],jedg[k,n],e] = npoin
#         if (ee>0) ibool[iedgR[k,nn],jedgR[k,nn],ee] = npoin end
#       end
#     end


#    #-- vertex nodes
#     for n = 1:4
#       i = ivtx[n]
#       j = jvtx[n]
#       if (ibool[i,j,e] > 0) continue end
#       npoin = npoin +1
#      # scan Vertex_Conn for shared elements
#       FE_GetVertexConn(se.fem, e, n, ees, nns)
#       for k = 1:SIZE(ees)
#         ii = ivtx[nns[k]]
#         jj = jvtx[nns[k]]
#         ibool[ii,jj,ees[k]] = npoin
#       end
#     end

#   end

#   FE_UnsetConnectivity(se.fem)

#   se.npoin = npoin

#   if echo_init
#     println(iout, fmtok, )
#     @printf(iout, "     %s%0i\n", "Total number of elements. . . . . . . . = ", se.nelem)
#     @printf(iout, "     %s%0i\n", "Total number of GLL points. . . . . . . = ", npoin)
#     println(iout, )
#     println(iout, "no", "Saving element/node table in binary file ibool_sem2d.dat")
#   end

#   inquire( IOLENGTH = iol ) se.ibool[:,:,1]
#   ounit = IO_new_unit()
#   OPEN(ounit, file = "ibool_sem2d.dat", status = "replace", access = "direct", recl = iol)
#   for e = 1:se.nelem
#     println(ounit, se.ibool[:,:,e])
#   end
#   CLOSE(ounit)
#   if (echo_init) println(iout, fmtok, ) end

#   return nothing

#     @label L100
#     format(5X, A, I0)

#   end


# # ======================================================================
# #
# #  set the global nodal coordinates
# #
#   function SE_init_coord(se)

# #  use echo, ONLY : echo_check,echo_init,iout,fmt1,fmtok
# #  use memory_info
# #  use stdio, ONLY : IO_new_unit

#   type(sem_grid_type), intent(inout) :: se

# #  integer :: i,j,e,ounit,iol
# #  DOUBLEPRECISION, POINTER :: coorg[:,:]

#   allocate(se.coord(NDIME, se.npoin))
#   storearray("coord", SIZE(se.coord), idouble)

# #---- Coordinates of the global points 
#   if (echo_init) println(iout, "no", "Defining nodes coordinates") end
#   for e = 1:se.nelem
#     coorg =Ref()(se.fem, e)
#     for j = 1:se.ngll
#     for i = 1:se.ngll
#       se.coord(:, se.ibool(i, j, e)) = matmul( coorg, se.shape(:, i, j) )
#     end
#     end
#     deallocate(coorg)
#   end
#   if (echo_init) println(iout, fmtok, ) end

#   if echo_check

# #----  Save the grid in a text file
#     println(iout, )
#     println(iout, "no", "Saving the grid coordinates (coord) in a text file")
#     ounit = IO_new_unit()
#     OPEN(ounit, file = "coord_sem2d.tab", status = "unknown")
#     println(ounit, se.npoin)
#     for i = 1:se.npoin
#       println(ounit, i, se.coord(:, i))
#     end
#     CLOSE(ounit)
#     println(iout, fmtok, )
#   end

# #----  Always save the grid in a binary file
#   if (echo_init) println(iout, "no", "Saving the grid coordinates (coord) in a binary file") end
#   inquire( IOLENGTH = iol ) Float32(se.coord(:, 1))
#   ounit = IO_new_unit()
#   OPEN(ounit, file = "coord_sem2d.dat", status = "replace", access = "direct", recl = iol)
#   for i = 1:se.npoin
#     println(ounit, Float32(se.coord(:, i)))
#   end
#   CLOSE(ounit)
#   if (echo_init) println(iout, fmtok, ) end

#   end



# # ======================================================================
# #
# # lagrange interpolation functions at (xi,eta)
# # interp(ngll*ngll) is then used as
# #   interpolated_value = matmul(interp,nodal_values)
# #   interpolated_vector(1:ndof) = matmul(interp,nodal_vectors(:,1:ndof))
# #
#   function SE_init_interpol(xi, eta, interp, grid)

# #  use gll, ONLY : hgll

# #  DOUBLEPRECISION, intent(in) :: xi,eta
#   type(sem_grid_type), intent(in) :: grid
# #  DOUBLEPRECISION, intent(out) :: interp[grid.ngll*grid.ngll]

# #  DOUBLEPRECISION :: fi,fj
# #  integer :: i,j,k

#   k = 0
#   for j = 1:grid.ngll
#     fj = hgll(j-1, eta, grid.xgll, grid.ngll)
#     for i = 1:grid.ngll
#       k = k+1
#       fi = hgll(i-1, xi, grid.xgll, grid.ngll)
#       interp[k] = fi*fj
#     end
#   end

#   end


# # ====================================================================
# #
#   function SE_find_nearest_node(coord_in, grid, iglob, coord, distance)

#   type(sem_grid_type), intent(in)  :: grid
# #  DOUBLEPRECISION   , intent(in)  :: coord_in[:]
# #  DOUBLEPRECISION, optional, intent(out) :: coord[:],distance
# #  integer, intent(out) :: iglob

# #  DOUBLEPRECISION :: d2min,d2
# #  integer :: ip

#   iglob = 0
#   d2min = huge(d2min)
#   for ip = 1:grid.npoin
#     d2 = (coord_in[1]-grid.coord[1,ip])^2 + (coord_in[2]-grid.coord[2,ip])^2
#     if d2 <= d2min
#       d2min  = d2
#       iglob = ip
#     end
#   end

#   if (present(coord)) coord = grid.coord[:,iglob] end
#   if (present(distance)) distance = sqrt(d2min) end

#   end

# # ====================================================================
# # Find (e,i,j), element and local indices associated to a global node iglob
# # Version 1: first element found
# # Version 2: all elements
# #
#   function SE_node_belongs_to_1(iglob, e, i, j, grid)

# #  integer, optional, intent(in) :: iglob
# #  integer, intent(out) :: i,j,e
#   type(sem_grid_type), intent(in)  :: grid

#   if (iglob>grid.npoin) IO_abort("SE_node_belongs_to: iglob out of range") end

#   for e = 1:grid.nelem
#   for j = 1:grid.ngll
#   for i = 1:grid.ngll
#     if (grid.ibool(i, j, e)==iglob) return end
#   end
#   end
#   end

#   end

# #---------------------------------------------------------------------

#   function SE_node_belongs_to_2(iglob, etab, itab, jtab, grid)

# #  integer, optional, intent(in) :: iglob
# #  integer, POINTER, dimension(:) :: itab,jtab,etab
#   type(sem_grid_type), intent(in)  :: grid

# #  integer :: i,j,e,nel
#   # WARNING: dimension should be larger than the maximum expected
#   # number of neighbour elements 
# #  integer, dimension(10) :: etmp,jtmp,itmp

#   if (iglob>grid.npoin) IO_abort("SE_node_belongs_to: iglob out of range") end

#   nel = 0
#   for e = 1:grid.nelem
#   for j = 1:grid.ngll
#   for i = 1:grid.ngll
#     if grid.ibool(i, j, e)==iglob
#       nel = nel +1
#       etmp[nel] = e
#       itmp[nel] = i
#       jtmp[nel] = j
#     end
#   end
#   end
#   end

#   if (@isdefined(etab) && !==(etab, nothing) #= ASSOCIATED(etab) =#) deallocate(etab) end
#   if (@isdefined(itab) && !==(itab, nothing) #= ASSOCIATED(itab) =#) deallocate(itab) end
#   if (@isdefined(jtab) && !==(jtab, nothing) #= ASSOCIATED(jtab) =#) deallocate(jtab) end
#   allocate(etab[nel], itab[nel], jtab[nel])
#   etab = etmp[1:nel]
#   itab = itmp[1:nel]
#   jtab = jtmp[1:nel]

#   end

# # ====================================================================
# #
#   function SE_find_point(coord, grid, e, xi, eta, new_coord)

# #  DOUBLEPRECISION, intent(in)  :: coord[NDIME]
#   type(sem_grid_type), intent(in) :: grid
# #  integer, intent(out) :: e
# #  DOUBLEPRECISION, intent(out) :: xi,eta,new_coord[NDIME]

# #  integer :: iglob,k,istat
# #  integer, POINTER, dimension(:) :: itab,jtab,etab

#   itab = nothing; jtab = nothing; etab = nothing #= NULLIFY(...) =#
#   SE_find_nearest_node(coord, grid, iglob)
#   SE_node_belongs_to(iglob, etab, itab, jtab, grid)
#   for k = 1:SIZE(etab)
#     xi = grid.xgll(itab[k])
#     eta = grid.xgll(jtab[k])
#     e = etab[k]
#     FE_find_point(coord, grid.fem, e, xi, eta, new_coord, istat)
#     if (istat==0) break end
#   end

#   deallocate(itab, jtab, etab)
#   if (istat>0) IO_abort("SE_find_point: point not found") end

#   end

# # ====================================================================
# # get the list of bulk indices of the EDGE nodes of an ELEMENT, 
# # assume counterclockwise numbering for the output list NODES

#   function SE_get_edge_nodes(grid, element, edge, nodes)

#   type(sem_grid_type), intent(in) :: grid
# #  integer, intent(in) :: element,edge
# #  integer, intent(out) :: nodes[grid.ngll]

#   #SELECTCASE (edge)
#     if edge == edge_D); nodes = grid.ibool(: || edge == 1 || edge == element
#     elseif edge == edge_R); nodes = grid.ibool(grid.ngll || edge == : || edge == element
#     elseif edge == edge_U); nodes = grid.ibool(grid.ngll:1:-1 || edge == grid.ngll || edge == element
#     elseif edge == edge_L); nodes = grid.ibool(1 || edge == grid.ngll:1:-1 || edge == element
#   end

#   end


# # ======================================================================
# #
# # Jacobian matrix  =  | dx/dxi   dx/deta |
# #                     | dz/dxi   dz/deta |
# #
#   function SE_Jacobian_e(sem, e) result(jac)

#   type(sem_grid_type), intent(in) :: sem
# #  integer, intent(in) :: e
# #  DOUBLEPRECISION :: jac[2,2,sem.ngll,sem.ngll]

# #  DOUBLEPRECISION, POINTER :: coorg[:,:]
# #  integer :: i,j

#   coorg =Ref()(sem.fem, e)
#   for j = 1:sem.ngll
#   for i = 1:sem.ngll
#     jac[:,:,i,j] = matmul(coorg, sem.dshape(:, :, i, j) )
#   end
#   end
#   deallocate(coorg)

#   end

# #-----------------------------------------------------------------------

#   function SE_Jacobian_eij(sem, e, i, j) result(jac)

#   type(sem_grid_type), intent(in) :: sem
# #  integer, intent(in) :: e,i,j
# #  DOUBLEPRECISION :: jac[2,2]
# #  DOUBLEPRECISION, POINTER :: coorg[:,:]

#   coorg =Ref()(sem.fem, e)
#   jac = matmul(coorg, sem.dshape(:, :, i, j) )
#   deallocate(coorg)

#   end


# # ======================================================================
# #
# #-- Compute weights(i,j) = wgll(i)*wgll(j)*dvol(i,j) for one element
# #
#   function SE_VolumeWeights_e(sem, e) result(dvol)

#   type(sem_grid_type), intent(in) :: sem
# #  integer, intent(in) :: e
# #  DOUBLEPRECISION :: dvol[sem.ngll,sem.ngll]

# #  DOUBLEPRECISION :: jac[2,2,sem.ngll,sem.ngll]

#   jac = SE_Jacobian_e(sem, e)
#   dvol = jac[1,1,:,:]*jac[2,2,:,:] - jac[1,2,:,:]*jac[2,1,:,:]
#   dvol = dvol*sem.wgll2

#   end

# #-----------------------------------------------------------------------

#   function SE_VolumeWeights_eij(sem, e, i, j) result(dvol)

#   type(sem_grid_type), intent(in) :: sem
# #  integer, intent(in) :: e,i,j
# #  DOUBLEPRECISION :: dvol

# #  DOUBLEPRECISION :: jac[2,2]

#   jac = SE_Jacobian_eij(sem, e, i, j)
#   dvol = ( jac[1,1]*jac[2,2] - jac[1,2]*jac[2,1] )*sem.wgll2(i, j)

#   end


# # ======================================================================
# #
# # Inverse Jacobian matrix  =  | dxi/dx    dxi/dz |
# #                             | deta/dx  deta/dz |
# #
#   function SE_InverseJacobian_e(sem, e) result(jaci)

# #  use utils, ONLY : invert2

#   type(sem_grid_type), intent(in) :: sem
# #  integer, intent(in) :: e
# #  DOUBLEPRECISION :: jaci[2,2,sem.ngll,sem.ngll]

# #  DOUBLEPRECISION :: jac[2,2,sem.ngll,sem.ngll]
# #  integer :: i,j

#   jac = SE_Jacobian_e(sem, e)
#   for j = 1:sem.ngll
#   for i = 1:sem.ngll
#     jaci[:,:,i,j] = invert2(jac[:,:,i,j])
#   end
#   end

#   end

# #----------------------------------------------------------------------

#   function SE_InverseJacobian_eij(sem, e, i, j) result(jaci)

# #  use utils, ONLY : invert2

#   type(sem_grid_type), intent(in) :: sem
# #  integer, intent(in) :: e,i,j
# #  DOUBLEPRECISION :: jaci[2,2]

# #  DOUBLEPRECISION :: jac[2,2]

#   jac = SE_Jacobian_eij(sem, e, i, j)
#   jaci = invert2(jac)

#   end


# # ======================================================================
# #
# function SE_BcTopoInit(grid)
#   type(sem_grid_type), intent(inout) :: grid
# #  integer :: i
#   grid.bounds =Ref().fem%bnds
#   for i = 1:SIZE(grid.bounds)
#     BC_set_bulk_node(grid.bounds(i), grid)
#   end
# end

# # ======================================================================
# #
# function BC_set_bulk_node(bc, grid)

# #  use generic_list, ONLY : Link_Ptr_Type,Link_Type,List_Type 
# #     ,LI_Init_List,LI_Add_To_Head,LI_Get_Head 
# #     ,LI_Get_Next,LI_Associated,LI_Remove_Head

# #  use utils, ONLY: drank

#   type(bnd_grid_type), intent(inout) :: bc
#   type(sem_grid_type), intent(in) :: grid

#   #type BC_Node_Data_Type
# #    integer :: bulk_node,bc_node
#   mutable struct bulk_node
#   end

#   #type BC_Node_Type
#     TYPE(Link_Type) :: Link
#   mutable struct Link
#     type(BC_Node_Data_Type), POINTER :: Data
#   end

#   #type BC_Node_Ptr_Type
#     type(BC_Node_Type), POINTER :: P
#   mutable struct P
#   end

# #  integer :: bulk_node,bulk_node_vect[grid.ngll],i,n,kloc,bc_npoin,bc_inode
# #  logical :: node_at_corner,new_node
#   type(List_Type)      :: BC_Nodes_List, BC_Corners_List
#   TYPE(Link_Ptr_Type)  :: Link, Sublink
#   TYPE(BC_Node_Ptr_Type)  :: BC_node, BC_corner

# # For reordering:
# #  DOUBLEPRECISION, allocatable :: BcCoord[:,:]
# #  integer, allocatable :: Isort[:],Iback[:]
# #  DOUBLEPRECISION :: Lx,Lz

#   bc.ngnod = grid.ngll
#   allocate(bc.ibool(bc.ngnod, bc.nelem))

# #-----------------------------------------------------------------------
# #  Build the database as a linked list

#   bc_npoin = 0
#   LI_Init_List(BC_Nodes_List)
#   LI_Init_List(BC_Corners_List)

#   for n = 1:bc.nelem

# #   get edge nodes, counterclockwise by default
#     SE_get_edge_nodes(grid, bc.elem(n), bc.edge(n), bulk_node_vect)

#     for kloc = 1:bc.ngnod

#       bulk_node = bulk_node_vect[kloc]
#       new_node = true

# #  When the node is an element corner (vertex) 
# #  check if it is already in the BC list.
# #  Note that the search is done in a BC_corners sublist.

#       node_at_corner = kloc==1 || kloc==bc.ngnod
#       if node_at_corner
#         SubLink = LI_Get_Head(BC_Corners_List)
#         while (LI_Associated(SubLink))
#           BC_corner = transfer(SubLink, BC_corner)
#           if BC_corner.P%Data.bulk_node == bulk_node
#             new_node = false
#             bc_inode = BC_corner.P%Data.bc_node
#             break
#           end
#           Sublink = LI_Get_Next(Sublink)
#         end
#       end

# #  If it is a new node, add it to the list of BC nodes ...
#       if new_node

#         bc_npoin = bc_npoin +1
#         bc_inode = bc_npoin
#         allocate(BC_node.P)
#         allocate(BC_node.P%data)
#         BC_node.P%data.bulk_node = bulk_node
#         BC_node.P%data.bc_node = bc_inode
#         Link = transfer(BC_node, Link)
#         LI_Add_To_Head(Link, BC_Nodes_List)

# #  ... and possibly to the sublist of BC corners.
#         if node_at_corner
#           allocate(BC_corner.P)
#           BC_corner.P%data =Ref().P%data
#           SubLink = transfer(BC_corner, SubLink)
#           LI_Add_To_Head(SubLink, BC_Corners_List)
#         end

#       end

#      # Set the [ (gll,bc_element) -> bc_node ] table
#       bc.ibool(kloc, n) = bc_inode

#     end
#   end

#  # clean up the list of corners
#   while true
#     SubLink = LI_Remove_Head(BC_Corners_List)
#     if (!LI_Associated(SubLink)) break end
#     BC_corner = transfer(SubLink, BC_corner)
#     deallocate(BC_corner.P)
#   end

# #-----------------------------------------------------------------------
# #  Translate BC database from linked list to array storage

#   bc.npoin = bc_npoin # = LI_Get_Len(BC_Nodes_List)
#   allocate(bc.node(bc.npoin))
#   for i = 1:bc.npoin
#     Link = LI_Remove_Head(BC_Nodes_List)
#     BC_node = transfer(Link, BC_node)
#     bc.node(BC_node.P%data.bc_node) = BC_node.P%data.bulk_node
#     deallocate(BC_node.P%data)
#     deallocate(BC_node.P)
#   end

# #-----------------------------------------------------------------------
# #  Sort BC nodes by increasing coordinate
# #  Use the coordinate with the largest range
#   allocate(BcCoord[2,bc.npoin])
#   allocate(Isort[bc.npoin])
#   BcCoord = grid.coord(:, bc.node)
#   Lx = maxval(BcCoord[1,:]) - minval(BcCoord[1,:])
#   Lz = maxval(BcCoord[2,:]) - minval(BcCoord[2,:])
#   if Lx>Lz
#     drank(BcCoord[1,:], Isort)
#   else
#     drank(BcCoord[2,:], Isort)
#   end
# #    write(51,'(I5,1X,I3)') ( bc%node(n),Isort(n),n=1,bc%npoin) 
#   bc.node = bc.node(Isort)
# #    write(51,'(I5)') bc%node
#   allocate(Iback[bc.npoin])
#   Iback[Isort] = [(ifori=1:bc.npoin)]
#   for n = 1:bc.nelem
#     bc.ibool(:, n) = Iback[bc.ibool(:end,n)]
#   end

#   deallocate(BcCoord, Isort, Iback)

# end


# # ======================================================================
# #
# # example:
# #
# #  integer, pointer :: ptr(:)
# #  ptr => SE_firstElementTagged(grid)
# #  ...
# #  deallocate(ptr)
# #
#   function SE_firstElementTagged(grid) result(elist)

# #  use utils, ONLY : UNIQUE

#   type(sem_grid_type), intent(in) :: grid
# #  integer, POINTER :: elist[:]

# #  integer, POINTER :: tags[:]
# #  integer :: ntags,k,e

#   tags =Ref()(grid.tag)
#   ntags = SIZE(tags)

#   allocate(elist[ntags])
#   for k = 1:ntags
#     for e = 1:grid.nelem
#       if (grid.tag(e)==tags[k]) break end
#     end
#     elist[k] = e
#   end

#   deallocate(tags)

#   end

# # ======================================================================
# #
#   function SE_inquire(grid, element, edge,
#                        itab, jtab, dim_t,
#                        size_min, size_max)

#   type(sem_grid_type), intent(in) :: grid
# #  integer, optional, intent(in) :: element,edge
# #  DOUBLEPRECISION, optional, intent(out) :: size_min,size_max
# #  integer, dimension(grid.ngll), optional, intent(out) :: itab,jtab
# #  integer, optional, intent(out) :: dim_t

# #  integer :: k
# #  DOUBLEPRECISION :: dmin,dmax

#   #-- give the GLL local numbering ITAB,JTAB for a given EDGE
#   #   assuming counterclockwise orientation
#   if present(edge) && present(itab) && present(jtab)
#     #SELECTCASE(edge)
#       case(edge_D); itab = [(kfork=1:grid.ngll)]   ; jtab = 1
#       elseif edge == edge_R); itab = grid.ngll                  ; jtab = [ (k for k = 1:grid.ngll) /
#       case(edge_U); itab = [(k,k=grid.ngll,1,-1)]; jtab = grid.ngll
#       elseif edge == edge_L); itab = 1                          ; jtab = [ (k || edge == k = grid.ngll || edge == 1 || edge == -1) /
#     end
#   end

#   #-- give local dimension (xi,eta) of tangent direction to an EDGE
#   if present(dim_t)
#     #SELECTCASE(edge)
#       case(edge_D, edge_U); dim_t = 1
#       case(edge_R, edge_L); dim_t = 2
#     end
#   end

#   if present(size_min) || present(size_max)
#     FE_GetElementSizes(dmax, dmin, element, grid.fem)
#     if (present(size_min)) size_min = dmin end
#     if (present(size_max)) size_max = dmax end
#   end

#   end

# # ======================================================================
# #
#   function SE_elem_coord(grid, e) result(ecoord)

#   type(sem_grid_type), intent(in) :: grid
# #  integer, intent(in) :: e

# #  DOUBLEPRECISION :: ecoord[2,grid.ngll,grid.ngll]
# #  integer :: i,j

#   for j = 1:grid.ngll
#   for i = 1:grid.ngll
#     ecoord[:,i,j] = grid.coord(:, grid.ibool(i, j, e))
#   end
#   end

#   end

# # ======================================================================
# #
#   function BC_inquire(bounds, tag, bc_topo_ptr)

#   type(bnd_grid_type), POINTER :: bounds[:], bc_topo_ptr
# #  integer, intent(in) :: tag

# #  integer :: i

#  # if no boundary exists with the requested tag a null pointer is returned
#   bc_topo_ptr = nothing #= NULLIFY(bc_topo_ptr) =#
#   for i = 1:SIZE(bounds)
#     if bounds[i]%tag == tag
#       bc_topo_ptr =Ref()[i]
#       return nothing
#     end
#   end

#   end

# # ======================================================================
# #
#   function BC_tag_exists(bounds, tag) result(exists)

#   type(bnd_grid_type), intent(in) :: bounds[:]
# #  integer, intent(in) :: tag

# #  integer :: i

#   exists = false
#   for i = 1:SIZE(bounds)
#     if bounds[i]%tag == tag
#       exists = true
#       break
#     end
#   end

#   end

# # ======================================================================
# # Normal to a boundary, pointing out of the element
# # Normals are assembled --> "average" normal between elements
#   function BC_get_normal_and_weights(bc, grid, NORM, W, periodic)

#   type(bnd_grid_type), intent(inout) :: bc
#   type(sem_grid_type), intent(in) :: grid
# #  DOUBLEPRECISION, intent(out) :: NORM[bc.npoin,2]
# #  DOUBLEPRECISION, intent(out) :: W[bc.npoin]
# #  logical, intent(in) :: periodic

# #  DOUBLEPRECISION :: DGlobDLoc[2,2],DxzDt[2],Tang[2],Jac1D,SignTang
# #  integer :: iGLLtab[bc.ngnod],jGLLtab[bc.ngnod]
# #  integer :: BcElmt,kGLL,BcNode,LocDimTanToEdge

#   NORM = 0e0
#   W = 0e0
#   for BcElmt = 1:bc.nelem
#     SE_inquire(grid, edge = bc.edge(BcElmt),
#                    itab = iGLLtab, jtab = jGLLtab, dim_t = LocDimTanToEdge)
#    # LocDimTanToEdge = local dimension (1=xi, 2=eta) tangent to current edge                    
#    # SignTang enforces the counterclockwise convention for tangent vector
#    # (iGLLtab,jGLLtab) are the GLL indices of the nodes of boundary element BcElmt
#    # in the same order (counterclockwise) as they appear in bc%ibool(:,BcElmt)
#     if bc.edge(BcElmt)==edge_U || bc.edge(BcElmt)==edge_L
#       SignTang = -1e0
#     else
#       SignTang = 1e0
#     end
#     for kGLL = 1:bc.ngnod
#       DGlobDLoc = SE_Jacobian(grid, bc.elem(BcElmt), iGLLtab[kGLL], jGLLtab[kGLL])
#       DxzDt = DGlobDLoc[:,LocDimTanToEdge]
#       Jac1D = sqrt( DxzDt[1]^2 + DxzDt[2]^2 )
#       Tang = SignTang*DxzDt/Jac1D
#       BcNode = bc.ibool(kGLL, BcElmt)
#                             # assembled, outwards normal
#       NORM[BcNode,:] = NORM[BcNode,:] + [Tang[2],-Tang[1]]
#       W[BcNode] = W[BcNode] + grid.wgll(kGLL)*Jac1D
#     end
#   end
#   if periodic
#     NORM[1,:] = NORM[1,:] + NORM[bc.npoin,:]
#     NORM[bc.npoin,:] = NORM[1,:]
#     W[1] = W[1] + W[bc.npoin]
#     W[bc.npoin] = W[1]
#   end
#   for BcElmt = 1:bc.nelem # fix interelement normals:
#     BcNode = bc.ibool(1, BcElmt)
#     NORM[BcNode,:] = NORM[BcNode,:] / sqrt( NORM[BcNode,1]^2 + NORM[BcNode,2]^2 )
#     BcNode = bc.ibool(bc.ngnod, BcElmt)
#     NORM[BcNode,:] = NORM[BcNode,:] / sqrt( NORM[BcNode,1]^2 + NORM[BcNode,2]^2 )
#   end

#   end


# end
