module write_netcdf

   use grid_params
 
   integer :: wr_ncid
   integer :: wrDimIDTime
   integer :: wrDimIDnCells
   integer :: wrDimIDnEdges
   integer :: wrDimIDnVertices
   integer :: wrDimIDmaxEdges
   integer :: wrDimIDmaxEdges2
   integer :: wrDimIDTWO
   integer :: wrDimIDvertexDegree
   integer :: wrDimIDnVertLevels
   integer :: wrDimIDnTracers
   integer :: wrVarIDlatCell
   integer :: wrVarIDlonCell
   integer :: wrVarIDmeshDensity
   integer :: wrVarIDxCell
   integer :: wrVarIDyCell
   integer :: wrVarIDzCell
   integer :: wrVarIDindexToCellID
   integer :: wrVarIDlatEdge
   integer :: wrVarIDlonEdge
   integer :: wrVarIDxEdge
   integer :: wrVarIDyEdge
   integer :: wrVarIDzEdge
   integer :: wrVarIDindexToEdgeID
   integer :: wrVarIDlatVertex
   integer :: wrVarIDlonVertex
   integer :: wrVarIDxVertex
   integer :: wrVarIDyVertex
   integer :: wrVarIDzVertex
   integer :: wrVarIDindexToVertexID
   integer :: wrVarIDcellsOnEdge
   integer :: wrVarIDnEdgesOnCell
   integer :: wrVarIDnEdgesOnEdge
   integer :: wrVarIDedgesOnCell
   integer :: wrVarIDedgesOnEdge
   integer :: wrVarIDweightsOnEdge
   integer :: wrVarIDdvEdge
   integer :: wrVarIDdv1Edge
   integer :: wrVarIDdv2Edge
   integer :: wrVarIDdcEdge
   integer :: wrVarIDangleEdge
   integer :: wrVarIDareaCell
   integer :: wrVarIDareaTriangle
   integer :: wrVarIDcellsOnCell
   integer :: wrVarIDverticesOnCell
   integer :: wrVarIDverticesOnEdge
   integer :: wrVarIDedgesOnVertex
   integer :: wrVarIDcellsOnVertex
   integer :: wrVarIDkiteAreasOnVertex
   integer :: wrVarIDfEdge
   integer :: wrVarIDfVertex
   integer :: wrVarIDh_s
   integer :: wrVarIDu
   integer :: wrVarIDv
   integer :: wrVarIDh
   integer :: wrVarIDvh
   integer :: wrVarIDcirculation
   integer :: wrVarIDvorticity
   integer :: wrVarIDke
   integer :: wrVarIDtracers
 
   integer :: wrLocalnCells
   integer :: wrLocalnEdges
   integer :: wrLocalnVertices
   integer :: wrLocalmaxEdges
   integer :: wrLocalnVertLevels
   integer :: wrLocalnTracers
 
   contains
 
   subroutine write_netcdf_init( &
                               nCells, &
                               nEdges, &
                               nVertices, &
                               maxEdges, &
                               nVertLevels, &
                               nTracers &
                               )
 
      implicit none
 
      include 'netcdf.inc'
 
      integer, intent(in) :: nCells
      integer, intent(in) :: nEdges
      integer, intent(in) :: nVertices
      integer, intent(in) :: maxEdges
      integer, intent(in) :: nVertLevels
      integer, intent(in) :: nTracers
 
      integer :: nferr
      integer, dimension(10) :: dimlist
      real (kind=8) :: sphere_radius
      character (len=16) :: on_a_sphere
 
 
      wrLocalnCells = nCells
      wrLocalnEdges = nEdges
      wrLocalnVertices = nVertices
      wrLocalmaxEdges = maxEdges
      wrLocalnVertLevels = nVertLevels
      wrLocalnTracers = nTracers
 
      on_a_sphere = 'YES             '
      sphere_radius = 1.0

      nferr = nf_create('grid.nc', IOR(NF_CLOBBER,NF_64BIT_OFFSET), wr_ncid)
 
      !
      ! Write Namlist information
      !

      nferr = nf_put_att_text(wr_ncid, NF_GLOBAL, 'on_a_sphere', 16, on_a_sphere)
      nferr = nf_put_att_double(wr_ncid, NF_GLOBAL, 'sphere_radius', NF_DOUBLE, 1, sphere_radius)

      !
      ! Define dimensions
      !
      nferr = nf_def_dim(wr_ncid, 'nCells', nCells, wrDimIDnCells)
      nferr = nf_def_dim(wr_ncid, 'nEdges', nEdges, wrDimIDnEdges)
      nferr = nf_def_dim(wr_ncid, 'nVertices', nVertices, wrDimIDnVertices)
      nferr = nf_def_dim(wr_ncid, 'maxEdges', maxEdges, wrDimIDmaxEdges)
      nferr = nf_def_dim(wr_ncid, 'vertexDegree', 3, wrDimIDvertexDegree)
 
      !
      ! Define variables
      !
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'latCell', NF_DOUBLE,  1, dimlist, wrVarIDlatCell)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'lonCell', NF_DOUBLE,  1, dimlist, wrVarIDlonCell)
      dimlist( 1) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'latVertex', NF_DOUBLE,  1, dimlist, wrVarIDlatVertex)
      dimlist( 1) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'lonVertex', NF_DOUBLE,  1, dimlist, wrVarIDlonVertex)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'nEdgesOnCell', NF_INT,  1, dimlist, wrVarIDnEdgesOnCell)
      dimlist( 1) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'areaCell', NF_DOUBLE,  1, dimlist, wrVarIDareaCell)
      dimlist( 1) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'areaTriangle', NF_DOUBLE,  1, dimlist, wrVarIDareaTriangle)
      dimlist( 1) = wrDimIDmaxEdges
      dimlist( 2) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'cellsOnCell', NF_INT,  2, dimlist, wrVarIDcellsOnCell)
      dimlist( 1) = wrDimIDmaxEdges
      dimlist( 2) = wrDimIDnCells
      nferr = nf_def_var(wr_ncid, 'verticesOnCell', NF_INT,  2, dimlist, wrVarIDverticesOnCell)
      dimlist( 1) = wrDimIDvertexDegree
      dimlist( 2) = wrDimIDnVertices
      nferr = nf_def_var(wr_ncid, 'cellsOnVertex', NF_INT,  2, dimlist, wrVarIDcellsOnVertex)
 
      nferr = nf_enddef(wr_ncid)
 
   end subroutine write_netcdf_init
 
 
   subroutine write_netcdf_fields( &
                                  time, &
                                  latCell, &
                                  lonCell, &
                                  meshDensity, &
                                  xCell, &
                                  yCell, &
                                  zCell, &
                                  indexToCellID, &
                                  latEdge, &
                                  lonEdge, &
                                  xEdge, &
                                  yEdge, &
                                  zEdge, &
                                  indexToEdgeID, &
                                  latVertex, &
                                  lonVertex, &
                                  xVertex, &
                                  yVertex, &
                                  zVertex, &
                                  indexToVertexID, &
                                  cellsOnEdge, &
                                  nEdgesOnCell, &
                                  nEdgesOnEdge, &
                                  edgesOnCell, &
                                  edgesOnEdge, &
                                  weightsOnEdge, &
                                  dvEdge, &
                                  dv1Edge, &
                                  dv2Edge, &
                                  dcEdge, &
                                  angleEdge, &
                                  areaCell, &
                                  areaTriangle, &
                                  cellsOnCell, &
                                  verticesOnCell, &
                                  verticesOnEdge, &
                                  edgesOnVertex, &
                                  cellsOnVertex, &
                                  kiteAreasOnVertex, &
                                  fEdge, &
                                  fVertex, &
                                  h_s, &
                                  u, &
                                  v, &
                                  h, &
                                  vh, &
                                  circulation, &
                                  vorticity, &
                                  ke, &
                                  tracers &
                                 )
 
      implicit none
 
      include 'netcdf.inc'
 
      integer, intent(in) :: time
      real (kind=RKIND), dimension(:), intent(in) :: latCell
      real (kind=RKIND), dimension(:), intent(in) :: lonCell
      real (kind=RKIND), dimension(:), intent(in) :: meshDensity
      real (kind=RKIND), dimension(:), intent(in) :: xCell
      real (kind=RKIND), dimension(:), intent(in) :: yCell
      real (kind=RKIND), dimension(:), intent(in) :: zCell
      integer, dimension(:), intent(in) :: indexToCellID
      real (kind=RKIND), dimension(:), intent(in) :: latEdge
      real (kind=RKIND), dimension(:), intent(in) :: lonEdge
      real (kind=RKIND), dimension(:), intent(in) :: xEdge
      real (kind=RKIND), dimension(:), intent(in) :: yEdge
      real (kind=RKIND), dimension(:), intent(in) :: zEdge
      integer, dimension(:), intent(in) :: indexToEdgeID
      real (kind=RKIND), dimension(:), intent(in) :: latVertex
      real (kind=RKIND), dimension(:), intent(in) :: lonVertex
      real (kind=RKIND), dimension(:), intent(in) :: xVertex
      real (kind=RKIND), dimension(:), intent(in) :: yVertex
      real (kind=RKIND), dimension(:), intent(in) :: zVertex
      integer, dimension(:), intent(in) :: indexToVertexID
      integer, dimension(:,:), intent(in) :: cellsOnEdge
      integer, dimension(:), intent(in) :: nEdgesOnCell
      integer, dimension(:), intent(in) :: nEdgesOnEdge
      integer, dimension(:,:), intent(in) :: edgesOnCell
      integer, dimension(:,:), intent(in) :: edgesOnEdge
      real (kind=RKIND), dimension(:,:), intent(in) :: weightsOnEdge
      real (kind=RKIND), dimension(:), intent(in) :: dvEdge
      real (kind=RKIND), dimension(:), intent(in) :: dv1Edge
      real (kind=RKIND), dimension(:), intent(in) :: dv2Edge
      real (kind=RKIND), dimension(:), intent(in) :: dcEdge
      real (kind=RKIND), dimension(:), intent(in) :: angleEdge
      real (kind=RKIND), dimension(:), intent(in) :: areaCell
      real (kind=RKIND), dimension(:), intent(in) :: areaTriangle
      integer, dimension(:,:), intent(in) :: cellsOnCell
      integer, dimension(:,:), intent(in) :: verticesOnCell
      integer, dimension(:,:), intent(in) :: verticesOnEdge
      integer, dimension(:,:), intent(in) :: edgesOnVertex
      integer, dimension(:,:), intent(in) :: cellsOnVertex
      real (kind=RKIND), dimension(:,:), intent(in) :: kiteAreasOnVertex
      real (kind=RKIND), dimension(:), intent(in) :: fEdge
      real (kind=RKIND), dimension(:), intent(in) :: fVertex
      real (kind=RKIND), dimension(:), intent(in) :: h_s
      real (kind=RKIND), dimension(:,:,:), intent(in) :: u
      real (kind=RKIND), dimension(:,:,:), intent(in) :: v
      real (kind=RKIND), dimension(:,:,:), intent(in) :: h
      real (kind=RKIND), dimension(:,:,:), intent(in) :: vh
      real (kind=RKIND), dimension(:,:,:), intent(in) :: circulation
      real (kind=RKIND), dimension(:,:,:), intent(in) :: vorticity
      real (kind=RKIND), dimension(:,:,:), intent(in) :: ke
      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: tracers
 
      integer :: nferr
      integer, dimension(1) :: start1, count1
      integer, dimension(2) :: start2, count2
      integer, dimension(3) :: start3, count3
      integer, dimension(4) :: start4, count4
 
      start1(1) = 1
 
      start2(1) = 1
      start2(2) = 1
 
      start3(1) = 1
      start3(2) = 1
      start3(3) = 1
 
      start4(1) = 1
      start4(2) = 1
      start4(3) = 1
      start4(4) = 1
 
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDlatCell, start1, count1, latCell)
 
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDlonCell, start1, count1, lonCell)
 
      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDlatVertex, start1, count1, latVertex)
 
      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDlonVertex, start1, count1, lonVertex)
 
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_int(wr_ncid, wrVarIDnEdgesOnCell, start1, count1, nEdgesOnCell)
 
      start1(1) = 1
      count1( 1) = wrLocalnCells
      nferr = nf_put_vara_double(wr_ncid, wrVarIDareaCell, start1, count1, areaCell)
 
      start1(1) = 1
      count1( 1) = wrLocalnVertices
      nferr = nf_put_vara_double(wr_ncid, wrVarIDareaTriangle, start1, count1, areaTriangle)
 
      start2(2) = 1
      count2( 1) = wrLocalmaxEdges
      count2( 2) = wrLocalnCells
      nferr = nf_put_vara_int(wr_ncid, wrVarIDcellsOnCell, start2, count2, cellsOnCell)
 
      start2(2) = 1
      count2( 1) = wrLocalmaxEdges
      count2( 2) = wrLocalnCells
      nferr = nf_put_vara_int(wr_ncid, wrVarIDverticesOnCell, start2, count2, verticesOnCell)
 
      start2(2) = 1
      count2( 1) = 3
      count2( 2) = wrLocalnVertices
      nferr = nf_put_vara_int(wr_ncid, wrVarIDcellsOnVertex, start2, count2, cellsOnVertex)
 
   end subroutine write_netcdf_fields
 
 
   subroutine write_netcdf_finalize()
 
      implicit none
 
      include 'netcdf.inc'
 
      integer :: nferr
 
      nferr = nf_close(wr_ncid)
 
   end subroutine write_netcdf_finalize
 
end module write_netcdf
 
