!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Program: permute_mesh
!
! This program may be used to permute the cells, edges, and vertices in an MPAS
!    grid file. 
!
! Example compilation command:
!    mpif90 -ffree-form -ffree-line-length-none -o permute_mesh permute_mesh.F -I${PNETCDF}/include -I${NETCDF}/include -L${PNETCDF}/lib -L${NETCDF}/lib -lnetcdff -lnetcdf -lpnetcdf
!
! Initial version: 07 November 2014, Michael Duda, NCAR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define NETCDF4_INPUT
!define PNETCDF_INPUT

#define NETCDF4_OUTPUT
!define PNETCDF_OUTPUT

!define TRANSFER_MESH
#define TRANSFER_STATIC

#define R4_INPUT
!define R8_INPUT

!define R4_OUTPUT
#define R8_OUTPUT

program permute_mesh

    use netcdf

    implicit none

#include "mpif.h"
#include "pnetcdf.inc"

#ifdef PNETCDF_INPUT
    integer, parameter :: noerr_input = NF_NOERR
#endif
#ifdef PNETCDF_OUTPUT
    integer, parameter :: noerr_output = NF_NOERR
#endif
#ifdef NETCDF4_INPUT
    integer, parameter :: noerr_input = NF90_NOERR
#endif
#ifdef NETCDF4_OUTPUT
    integer, parameter :: noerr_output = NF90_NOERR
#endif

    integer, parameter :: max_dictionary_size = 16500000   ! should be at least (nEdges/32)
    integer, parameter :: rw_chunk_size = 8000000

    integer :: int_size = 32
    integer, dimension(max_dictionary_size) :: dictionary_array


    integer :: i, ierr
    integer :: ncid_in, ncid_out, dim_out
    character (len=256) :: grid_filename, partition_filename

    integer, dimension(:), pointer :: cell_order, edge_order, vertex_order
    integer, dimension(:), pointer :: cell_lookup, edge_lookup, vertex_lookup
    integer :: nparts

    character (len=512) :: outname


#ifdef TRANSFER_MESH
    character (len=32), dimension(7) :: dims = [ character(len=32) :: &
                                                'nCells', &
                                                'nEdges', &
                                                'nVertices', &
                                                'maxEdges', &
                                                'maxEdges2', &
                                                'TWO', &
                                                'vertexDegree' &
                                               ]

    character (len=32), dimension(36) :: fields = [ character(len=32) :: &
                                                   'latCell', &
                                                   'lonCell', &
                                                   'meshDensity', &
                                                   'xCell', &
                                                   'yCell', &
                                                   'zCell', &
                                                   'indexToCellID', &
                                                   'latEdge', &
                                                   'lonEdge', &
                                                   'xEdge', &
                                                   'yEdge', &
                                                   'zEdge', &
                                                   'indexToEdgeID', &
                                                   'latVertex', &
                                                   'lonVertex', &
                                                   'xVertex', &
                                                   'yVertex', &
                                                   'zVertex', &
                                                   'indexToVertexID', &
                                                   'cellsOnEdge', &
                                                   'nEdgesOnCell', &
                                                   'nEdgesOnEdge', &
                                                   'edgesOnCell', &
                                                   'edgesOnEdge', &
                                                   'weightsOnEdge', &
                                                   'dvEdge', &
                                                   'dcEdge', &
                                                   'angleEdge', &
                                                   'areaCell', &
                                                   'areaTriangle', &
                                                   'cellsOnCell', &
                                                   'verticesOnCell', &
                                                   'verticesOnEdge', &
                                                   'edgesOnVertex', &
                                                   'cellsOnVertex', &
                                                   'kiteAreasOnVertex' &
                                                  ]

    logical, dimension(36) :: reindex = [ &
                                         .false., & !  'latCell'
                                         .false., & !  'lonCell'
                                         .false., & !  'meshDensity'
                                         .false., & !  'xCell'
                                         .false., & !  'yCell'
                                         .false., & !  'zCell'
                                         .true.,  & !  'indexToCellID'
                                         .false., & !  'latEdge'
                                         .false., & !  'lonEdge'
                                         .false., & !  'xEdge'
                                         .false., & !  'yEdge'
                                         .false., & !  'zEdge'
                                         .true.,  & !  'indexToEdgeID'
                                         .false., & !  'latVertex'
                                         .false., & !  'lonVertex'
                                         .false., & !  'xVertex'
                                         .false., & !  'yVertex'
                                         .false., & !  'zVertex'
                                         .true.,  & !  'indexToVertexID'
                                         .true.,  & !  'cellsOnEdge'
                                         .false., & !  'nEdgesOnCell'
                                         .false., & !  'nEdgesOnEdge'
                                         .true.,  & !  'edgesOnCell'
                                         .true.,  & !  'edgesOnEdge'
                                         .false., & !  'weightsOnEdge'
                                         .false., & !  'dvEdge'
                                         .false., & !  'dcEdge'
                                         .false., & !  'angleEdge'
                                         .false., & !  'areaCell'
                                         .false., & !  'areaTriangle'
                                         .true.,  & !  'cellsOnCell'
                                         .true.,  & !  'verticesOnCell'
                                         .true.,  & !  'verticesOnEdge'
                                         .true.,  & !  'edgesOnVertex'
                                         .true.,  & !  'cellsOnVertex'
                                         .false.  & !  'kiteAreasOnVertex'
                                                  ]
#endif


#ifdef TRANSFER_STATIC
    character (len=32), dimension(12) :: dims = [ character(len=32) :: &
                                                'nCells', &
                                                'nEdges', &
                                                'nVertices', &
                                                'maxEdges', &
                                                'maxEdges2', &
                                                'TWO', &
                                                'vertexDegree', &
                                                'FIFTEEN', &
                                                'nMonths', &
                                                'R3', &
                                                'TWENTYONE',  &
                                                'StrLen'  &
                                               ]

    character (len=32), dimension(69) :: fields = [ character(len=32) :: &
                                                   'latCell', &
                                                   'lonCell', &
                                                   'meshDensity', &
                                                   'xCell', &
                                                   'yCell', &
                                                   'zCell', &
                                                   'indexToCellID', &
                                                   'latEdge', &
                                                   'lonEdge', &
                                                   'xEdge', &
                                                   'yEdge', &
                                                   'zEdge', &
                                                   'indexToEdgeID', &
                                                   'latVertex', &
                                                   'lonVertex', &
                                                   'xVertex', &
                                                   'yVertex', &
                                                   'zVertex', &
                                                   'indexToVertexID', &
                                                   'cellsOnEdge', &
                                                   'nEdgesOnCell', &
                                                   'nEdgesOnEdge', &
                                                   'edgesOnCell', &
                                                   'edgesOnEdge', &
                                                   'weightsOnEdge', &
                                                   'dvEdge', &
                                                   'dcEdge', &
                                                   'angleEdge', &
                                                   'areaCell', &
                                                   'areaTriangle', &
                                                   'cellsOnCell', &
                                                   'verticesOnCell', &
                                                   'verticesOnEdge', &
                                                   'edgesOnVertex', &
                                                   'cellsOnVertex', &
                                                   'kiteAreasOnVertex', &
                                                   'edgeNormalVectors', &
                                                   'localVerticalUnitVectors', &
                                                   'cellTangentPlane', &
                                                   'fEdge', &
                                                   'fVertex', &
                                                   'ter', &
                                                   'landmask', &
                                                   'ivgtyp', &
                                                   'mminlu', &
                                                   'isltyp', &
                                                   'soilcat_bot', &
                                                   'snoalb', &
                                                   'soiltemp', &
                                                   'greenfrac', &
                                                   'shdmin', &
                                                   'shdmax', &
                                                   'albedo12m', &
                                                   'varsso', &
                                                   'var2d', &
                                                   'con', &
                                                   'oa1', &
                                                   'oa2', &
                                                   'oa3', &
                                                   'oa4', &
                                                   'ol1', &
                                                   'ol2', &
                                                   'ol3', &
                                                   'ol4', &
                                                   'deriv_two', &
                                                   'advCells', &
                                                   'defc_a', &
                                                   'defc_b', &
                                                   'coeffs_reconstruct' &
                                                  ]

    logical, dimension(69) :: reindex = [ &
                                         .false., & !  'latCell'
                                         .false., & !  'lonCell'
                                         .false., & !  'meshDensity'
                                         .false., & !  'xCell'
                                         .false., & !  'yCell'
                                         .false., & !  'zCell'
                                         .true.,  & !  'indexToCellID'
                                         .false., & !  'latEdge'
                                         .false., & !  'lonEdge'
                                         .false., & !  'xEdge'
                                         .false., & !  'yEdge'
                                         .false., & !  'zEdge'
                                         .true.,  & !  'indexToEdgeID'
                                         .false., & !  'latVertex'
                                         .false., & !  'lonVertex'
                                         .false., & !  'xVertex'
                                         .false., & !  'yVertex'
                                         .false., & !  'zVertex'
                                         .true.,  & !  'indexToVertexID'
                                         .true.,  & !  'cellsOnEdge'
                                         .false., & !  'nEdgesOnCell'
                                         .false., & !  'nEdgesOnEdge'
                                         .true.,  & !  'edgesOnCell'
                                         .true.,  & !  'edgesOnEdge'
                                         .false., & !  'weightsOnEdge'
                                         .false., & !  'dvEdge'
                                         .false., & !  'dcEdge'
                                         .false., & !  'angleEdge'
                                         .false., & !  'areaCell'
                                         .false., & !  'areaTriangle'
                                         .true.,  & !  'cellsOnCell'
                                         .true.,  & !  'verticesOnCell'
                                         .true.,  & !  'verticesOnEdge'
                                         .true.,  & !  'edgesOnVertex'
                                         .true.,  & !  'cellsOnVertex'
                                         .false., & !  'kiteAreasOnVertex'
                                         .false., & !  'edgeNormalVectors', &
                                         .false., & !  'localVerticalUnitVectors', &
                                         .false., & !  'cellTangentPlane', &
                                         .false., & !  'fEdge', &
                                         .false., & !  'fVertex', &
                                         .false., & !  'ter', &
                                         .false., & !  'landmask', &
                                         .false., & !  'ivgtyp', &
                                         .false., & !  'mminlu', &
                                         .false., & !  'isltyp', &
                                         .false., & !  'soilcat_bot', &
                                         .false., & !  'snoalb', &
                                         .false., & !  'soiltemp', &
                                         .false., & !  'greenfrac', &
                                         .false., & !  'shdmin', &
                                         .false., & !  'shdmax', &
                                         .false., & !  'albedo12m', &
                                         .false., & !  'varsso', &
                                         .false., & !  'var2d', &
                                         .false., & !  'con', &
                                         .false., & !  'oa1', &
                                         .false., & !  'oa2', &
                                         .false., & !  'oa3', &
                                         .false., & !  'oa4', &
                                         .false., & !  'ol1', &
                                         .false., & !  'ol2', &
                                         .false., & !  'ol3', &
                                         .false., & !  'ol4', &
                                         .false., & !  'deriv_two', &
                                         .false., & !  'advCells', &
                                         .false., & !  'defc_a', &
                                         .false., & !  'defc_b', &
                                         .false.  & !  'coeffs_reconstruct', &
                                        ]
#endif

    character (len=32), dimension(2) :: atts = [ character(len=32) :: &
                                                'on_a_sphere', &
                                                'sphere_radius' &
                                               ]


    !
    ! Requires two command-line arguments: the name of the grid file and the name of the Metis graph partition file
    !
    if (command_argument_count() /= 2) then
        write(0,*) ' '
        write(0,*) 'Usage: permute_mesh <grid file> <partition file>'
        write(0,*) ' '
        stop
    end if

    call get_command_argument(1, grid_filename)
    call get_command_argument(2, partition_filename)

    call MPI_Init(ierr)

#ifdef PNETCDF_INPUT
    ierr = nfmpi_open(MPI_COMM_WORLD, trim(grid_filename), NF_NOWRITE, MPI_INFO_NULL, ncid_in)
#endif
#ifdef NETCDF4_INPUT
    ierr = nf90_open(trim(grid_filename), NF90_NOWRITE, ncid_in)
#endif
    if (ierr /= noerr_input) then
        write(0,*) '********************************************************************************'
        write(0,*) 'Error: Could not open input grid file '//trim(grid_filename)
        write(0,*) '********************************************************************************'
        stop
    end if

    !
    ! Get permutation arrays for cells, edges, and vertices
    !
    call calc_permutation(ncid_in, partition_filename, &
                          cell_order, cell_lookup, &
                          edge_order, edge_lookup, &
                          vertex_order, vertex_lookup, &
                          nparts)

    if (nparts < 10) then
        write(outname, '(a,i1.1,a)') 'permuted_grid.', nparts, '.nc'
    else if (nparts < 100) then
        write(outname, '(a,i2.2,a)') 'permuted_grid.', nparts, '.nc'
    else if (nparts < 1000) then
        write(outname, '(a,i3.3,a)') 'permuted_grid.', nparts, '.nc'
    else if (nparts < 10000) then
        write(outname, '(a,i4.4,a)') 'permuted_grid.', nparts, '.nc'
    else if (nparts < 100000) then
        write(outname, '(a,i5.5,a)') 'permuted_grid.', nparts, '.nc'
    else if (nparts < 1000000) then
        write(outname, '(a,i6.6,a)') 'permuted_grid.', nparts, '.nc'
    else if (nparts < 10000000) then
        write(outname, '(a,i7.7,a)') 'permuted_grid.', nparts, '.nc'
    end if 

#ifdef NETCDF4_OUTPUT
    ierr = nf90_create(trim(outname), ior(NF90_NETCDF4,NF90_MPIIO), ncid_out)
    if (ierr == NF90_EEXIST) then
#endif
#ifdef PNETCDF_OUTPUT
    ierr = nfmpi_create(MPI_COMM_WORLD, trim(outname), NF_64BIT_DATA, MPI_INFO_NULL, ncid_out)
    if (ierr == NF_EEXIST) then
#endif
        write(0,*) '********************************************************************************'
        write(0,*) 'Error: File '''//trim(outname)//''' already exists.'
        write(0,*) '********************************************************************************'
        stop
    else if (ierr /= noerr_output) then
        write(0,*) '********************************************************************************'
        write(0,*) 'Error: Could not create permuted grid file '''//trim(outname)//'''.'
        write(0,*) 'nferr = ', ierr
        write(0,*) '********************************************************************************'
        stop
    end if

    !
    ! Transfer dimensions
    !
#ifdef NETCDF4_OUTPUT
    ierr = nf90_def_dim(ncid_out, 'Time', NF90_UNLIMITED, dim_out)
#endif
#ifdef PNETCDF_OUTPUT
    ierr = nfmpi_def_dim(ncid_out, 'Time', NFMPI_UNLIMITED, dim_out)
#endif
    if (ierr /= noerr_output) then
        write(0,*) '********************************************************************************'
        write(0,*) 'Error: Creating Time dimension in output file.'
        write(0,*) '********************************************************************************'
        stop
    end if

    do i=1,size(dims)
        call xfer_dim(ncid_in, ncid_out, trim(dims(i)))
    end do

    !
    ! Transfer field definitions
    !
    do i=1,size(fields)
        call xfer_var(ncid_in, ncid_out, trim(fields(i)))
    end do

    !
    ! Transfer global attributes
    !
    do i=1,size(atts)
        call xfer_att(ncid_in, ncid_out, trim(atts(i)))
    end do


#ifdef NETCDF4_OUTPUT
    ierr = nf90_enddef(ncid_out)
#endif
#ifdef PNETCDF_OUTPUT
    ierr = nfmpi_enddef(ncid_out)
    ierr = nfmpi_begin_indep_data(ncid_out)
#endif


    !
    ! Transfer fields themselves
    !
    do i=1,size(fields)
        write(0,*) 'Transferring '//trim(fields(i))
        call copy_var(ncid_in, ncid_out, trim(fields(i)), reindex(i), &
                      cell_order, cell_lookup, &
                      edge_order, edge_lookup, &
                      vertex_order, vertex_lookup)
    end do

#ifdef PNETCDF_INPUT
    ierr = nfmpi_close(ncid_in)
#endif
#ifdef NETCDF4_INPUT
    ierr = nf90_close(ncid_in)
#endif
#ifdef NETCDF4_OUTPUT
    ierr = nf90_close(ncid_out)
#endif
#ifdef PNETCDF_OUTPUT
    ierr = nfmpi_close(ncid_out)
#endif

    write(0,*) ' '
    write(0,*) 'Successfully created permuted grid file ''permuted_grid.nc'''
    write(0,*) ' '

    deallocate(cell_order)
    deallocate(edge_order)
    deallocate(vertex_order)
    deallocate(cell_lookup)
    deallocate(edge_lookup)
    deallocate(vertex_lookup)

    call MPI_Finalize(ierr)

    stop


contains


    subroutine calc_permutation(ncid_in, partition_filename, &
                                cell_order, cell_lookup_out, &
                                edge_order, edge_lookup_out, &
                                vertex_order, vertex_lookup_out, &
                                nparts)

        implicit none

        integer, intent(in) :: ncid_in
        character (len=*), intent(in) :: partition_filename
        integer, dimension(:), pointer :: cell_order
        integer, dimension(:), pointer :: cell_lookup_out
        integer, dimension(:), pointer :: edge_order
        integer, dimension(:), pointer :: edge_lookup_out
        integer, dimension(:), pointer :: vertex_order
        integer, dimension(:), pointer :: vertex_lookup_out
        integer, intent(out) :: nparts

        integer :: i, j, k
        integer :: nCellsID, nEdgesID, nVerticesID, maxEdgesID, ierr
        integer :: nCells, nEdges, nVertices, maxEdges
        integer :: nEdgesOnCellID, cellsOnEdgeID, cellsOnVertexID, edgesOnCellID, verticesOnCellID
        integer :: ntasks
        integer, dimension(:), pointer :: cell_owner, edge_owner, vertex_owner
        integer, dimension(:), pointer :: task_ncells, task_nedges, task_nvertices
        integer, dimension(:), pointer :: task_cell_start, task_edge_start, task_vertex_start
        integer, dimension(:), pointer :: nEdgesOnCell
        integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, edgesOnCell, verticesOnCell
        character(len=32) :: graph_filename

        integer, dimension(:,:), pointer :: cell_lookup
        integer, dimension(:,:), pointer :: edge_lookup
        integer, dimension(:,:), pointer :: vertex_lookup

#ifdef PNETCDF_INPUT
        integer (kind=MPI_OFFSET_KIND), dimension(1) :: start1, count1
        integer (kind=MPI_OFFSET_KIND), dimension(2) :: start2, count2
#endif
#ifdef NETCDF4_INPUT
        integer, dimension(1) :: start1, count1
        integer, dimension(2) :: start2, count2
#endif

#ifdef PNETCDF_INPUT
        ierr = nfmpi_inq_dimid(ncid_in, 'nCells', nCellsID)
        ierr = nfmpi_inq_dimlen(ncid_in, nCellsID, nCells)
        write(0,*) 'nCells=',nCells

        ierr = nfmpi_inq_dimid(ncid_in, 'nEdges', nEdgesID)
        ierr = nfmpi_inq_dimlen(ncid_in, nEdgesID, nEdges)
        write(0,*) 'nEdges=',nEdges

        ierr = nfmpi_inq_dimid(ncid_in, 'nVertices', nVerticesID)
        ierr = nfmpi_inq_dimlen(ncid_in, nVerticesID, nVertices)
        write(0,*) 'nVertices=',nVertices

        ierr = nfmpi_inq_dimid(ncid_in, 'maxEdges', maxEdgesID)
        ierr = nfmpi_inq_dimlen(ncid_in, maxEdgesID, maxEdges)
        write(0,*) 'maxEdges=',maxEdges
#endif
#ifdef NETCDF4_INPUT
        ierr = nf90_inq_dimid(ncid_in, 'nCells', nCellsID)
        ierr = nf90_inquire_dimension(ncid_in, nCellsID, len=nCells)
        write(0,*) 'nCells=',nCells

        ierr = nf90_inq_dimid(ncid_in, 'nEdges', nEdgesID)
        ierr = nf90_inquire_dimension(ncid_in, nEdgesID, len=nEdges)
        write(0,*) 'nEdges=',nEdges

        ierr = nf90_inq_dimid(ncid_in, 'nVertices', nVerticesID)
        ierr = nf90_inquire_dimension(ncid_in, nVerticesID, len=nVertices)
        write(0,*) 'nVertices=',nVertices

        ierr = nf90_inq_dimid(ncid_in, 'maxEdges', maxEdgesID)
        ierr = nf90_inquire_dimension(ncid_in, maxEdgesID, len=maxEdges)
        write(0,*) 'maxEdges=',maxEdges
#endif

        ntasks = 0 
        allocate(cell_owner(nCells))
        allocate(edge_owner(nEdges))
        allocate(vertex_owner(nVertices))

        !
        ! Read partition file to determine ownership of cells
        !
        open(42, file=trim(partition_filename), form='formatted', iostat=ierr)
        if (ierr /= 0) then
            write(0,*) 'Error opening partition file '''//trim(partition_filename)//''''
            stop
        end if
        do i=1,nCells
            read(42,*,iostat=ierr) cell_owner(i)
            if (ierr /= 0) then
                write(0,*) 'Error reading partition file '''//trim(partition_filename)//''''
                stop
            end if
            if (cell_owner(i) > ntasks) then
                ntasks = cell_owner(i)
            end if 
        end do
        close(42)
        ntasks = ntasks + 1
        nparts = ntasks

        write(0,*) 'Found ', ntasks, ' tasks in partition file'
        
#ifdef PNETCDF_INPUT
        ierr = nfmpi_begin_indep_data(ncid_in)
#endif

        !
        ! Determine ownership of edges
        !
        allocate(nEdgesOnCell(nCells))
        start1(1) = 1
        count1(1) = nCells
#ifdef PNETCDF_INPUT
        ierr = nfmpi_inq_varid(ncid_in, 'nEdgesOnCell', nEdgesOnCellID)
        ierr = nfmpi_get_var_int(ncid_in, nEdgesOnCellID, nEdgesOnCell)
#endif
#ifdef NETCDF4_INPUT
        ierr = nf90_inq_varid(ncid_in, 'nEdgesOnCell', nEdgesOnCellID)
        ierr = nf90_get_var(ncid_in, nEdgesOnCellID, nEdgesOnCell)
#endif
        if (ierr /= noerr_input) then
            write(0,*) 'Error reading nEdgesOnCell ', ierr
        end if

        allocate(cellsOnEdge(2,nEdges))
#ifdef PNETCDF_INPUT
        ierr = nfmpi_inq_varid(ncid_in, 'cellsOnEdge', cellsOnEdgeID)
#endif
#ifdef NETCDF4_INPUT
        ierr = nf90_inq_varid(ncid_in, 'cellsOnEdge', cellsOnEdgeID)
#endif

        start2(1) = 1
        count2(1) = int(2,MPI_OFFSET_KIND)
        start2(2) = 1
        count2(2) = min(rw_chunk_size, nEdges)

        do while (start2(2) <= int(nEdges,MPI_OFFSET_KIND)) 
#ifdef PNETCDF_INPUT
            ierr = nfmpi_get_vara_int(ncid_in, cellsOnEdgeID, start2, count2, &
                                      cellsOnEdge(:,start2(2):(start2(2)+count2(2)-1)))
#endif
#ifdef NETCDF4_INPUT
            ierr = nf90_get_var(ncid_in, cellsOnEdgeID, &
                                cellsOnEdge(:,start2(2):(start2(2)+count2(2)-1)), start=start2, count=count2)
#endif
            if (ierr /= noerr_input) then
                write(0,*) 'Error reading cellsOnEdge ', ierr
            end if

            start2(2) = start2(2) + count2(2)
            count2(2) = min(rw_chunk_size, nEdges - start2(2) + 1)
        end do


        allocate(edgesOnCell(maxEdges,nCells))
#ifdef PNETCDF_INPUT
        ierr = nfmpi_inq_varid(ncid_in, 'edgesOnCell', edgesOnCellID)
#endif
#ifdef NETCDF4_INPUT
        ierr = nf90_inq_varid(ncid_in, 'edgesOnCell', edgesOnCellID)
#endif

        start2(1) = 1
        count2(1) = int(maxEdges,MPI_OFFSET_KIND)
        start2(2) = 1
        count2(2) = min(rw_chunk_size, nCells)

        do while (start2(2) <= int(nCells,MPI_OFFSET_KIND)) 
#ifdef PNETCDF_INPUT
            ierr = nfmpi_get_vara_int(ncid_in, edgesOnCellID, start2, count2, &
                                      edgesOnCell(:,start2(2):(start2(2)+count2(2)-1)))
#endif
#ifdef NETCDF4_INPUT
            ierr = nf90_get_var(ncid_in, edgesOnCellID, &
                                edgesOnCell(:,start2(2):(start2(2)+count2(2)-1)), start=start2, count=count2)
#endif
            if (ierr /= noerr_input) then
                write(0,*) 'Error reading edgesOnCell ', ierr
            end if

            start2(2) = start2(2) + count2(2)
            count2(2) = min(rw_chunk_size, nCells - start2(2) + 1)
        end do

        do i=1,nEdges
            edge_owner(i) = cell_owner(cellsOnEdge(1,i))
        end do


        !
        ! Determine ownership of vertices
        !
        allocate(cellsOnVertex(3,nVertices))
#ifdef PNETCDF_INPUT
        ierr = nfmpi_inq_varid(ncid_in, 'cellsOnVertex', cellsOnVertexID)
#endif
#ifdef NETCDF4_INPUT
        ierr = nf90_inq_varid(ncid_in, 'cellsOnVertex', cellsOnVertexID)
#endif

        start2(1) = 1
        count2(1) = int(3,MPI_OFFSET_KIND)
        start2(2) = 1
        count2(2) = min(rw_chunk_size, nVertices)

        do while (start2(2) <= int(nVertices,MPI_OFFSET_KIND)) 
#ifdef PNETCDF_INPUT
            ierr = nfmpi_get_vara_int(ncid_in, cellsOnVertexID, start2, count2, &
                                      cellsOnVertex(:,start2(2):(start2(2)+count2(2)-1)))
#endif
#ifdef NETCDF4_INPUT
            ierr = nf90_get_var(ncid_in, cellsOnVertexID, &
                                cellsOnVertex(:,start2(2):(start2(2)+count2(2)-1)), start=start2, count=count2)
#endif
            if (ierr /= noerr_input) then
                write(0,*) 'Error reading cellsOnVertex ', ierr
            end if

            start2(2) = start2(2) + count2(2)
            count2(2) = min(rw_chunk_size, nVertices - start2(2) + 1)
        end do


        allocate(verticesOnCell(maxEdges,nCells))
#ifdef PNETCDF_INPUT
        ierr = nfmpi_inq_varid(ncid_in, 'verticesOnCell', verticesOnCellID)
#endif
#ifdef NETCDF4_INPUT
        ierr = nf90_inq_varid(ncid_in, 'verticesOnCell', verticesOnCellID)
#endif

        start2(1) = 1
        count2(1) = int(maxEdges,MPI_OFFSET_KIND)
        start2(2) = 1
        count2(2) = min(rw_chunk_size, nCells)

        do while (start2(2) <= int(nCells,MPI_OFFSET_KIND)) 
#ifdef PNETCDF_INPUT
            ierr = nfmpi_get_vara_int(ncid_in, verticesOnCellID, start2, count2, &
                                      verticesOnCell(:,start2(2):(start2(2)+count2(2)-1)))
#endif
#ifdef NETCDF4_INPUT
            ierr = nf90_get_var(ncid_in, verticesOnCellID, &
                                verticesOnCell(:,start2(2):(start2(2)+count2(2)-1)), start=start2, count=count2)
#endif
            if (ierr /= noerr_input) then
                write(0,*) 'Error reading verticesOnCell ', ierr
            end if

            start2(2) = start2(2) + count2(2)
            count2(2) = min(rw_chunk_size, nCells - start2(2) + 1)
        end do

        do i=1,nVertices
            vertex_owner(i) = cell_owner(cellsOnVertex(1,i))
        end do

        !
        ! Determine how many cells/edges/vertices are owned by each task
        !
        allocate(task_ncells(ntasks))
        allocate(task_nedges(ntasks))
        allocate(task_nvertices(ntasks))

        task_ncells(:) = 0
        do i=1,nCells
            task_ncells(cell_owner(i)+1) = task_ncells(cell_owner(i)+1) + 1
        end do

        task_nedges(:) = 0
        do i=1,nEdges
            task_nedges(edge_owner(i)+1) = task_nedges(edge_owner(i)+1) + 1
        end do

        task_nvertices(:) = 0
        do i=1,nVertices
            task_nvertices(vertex_owner(i)+1) = task_nvertices(vertex_owner(i)+1) + 1
        end do

        !
        ! Write out a new graph decomposition file
        !
        if (ntasks < 10) then
            write(graph_filename, '(a,i1.1)') 'permuted_graph.info.part.', ntasks
        else if (ntasks < 100) then
            write(graph_filename, '(a,i2.2)') 'permuted_graph.info.part.', ntasks
        else if (ntasks < 1000) then
            write(graph_filename, '(a,i3.3)') 'permuted_graph.info.part.', ntasks
        else if (ntasks < 10000) then
            write(graph_filename, '(a,i4.4)') 'permuted_graph.info.part.', ntasks
        else if (ntasks < 100000) then
            write(graph_filename, '(a,i5.5)') 'permuted_graph.info.part.', ntasks
        else if (ntasks < 1000000) then
            write(graph_filename, '(a,i6.6)') 'permuted_graph.info.part.', ntasks
        else if (ntasks < 10000000) then
            write(graph_filename, '(a,i7.7)') 'permuted_graph.info.part.', ntasks
        end if 

        open(43, file=trim(graph_filename), form='formatted')
        do i=1,ntasks
            do j=1,task_ncells(i)
                write(43,*) i-1
            end do
        end do
        close(43)


        !
        ! Determine starting cell/edge/vertex index in permuted mesh for each task 
        !
        allocate(task_cell_start(ntasks))
        allocate(task_edge_start(ntasks))
        allocate(task_vertex_start(ntasks))

        task_cell_start(1) = 1
        task_edge_start(1) = 1
        task_vertex_start(1) = 1
        do i=2,ntasks
            task_cell_start(i) = task_cell_start(i-1) + task_ncells(i-1)
            task_edge_start(i) = task_edge_start(i-1) + task_nedges(i-1)
            task_vertex_start(i) = task_vertex_start(i-1) + task_nvertices(i-1)
        end do
        deallocate(task_ncells)
        deallocate(task_nedges)
        deallocate(task_nvertices)

        !
        ! Determine permuted cell/edge/vertex order
        ! (Where a new cell/edge/vertex originates from)
        !
        allocate(cell_order(nCells)) 
        allocate(edge_order(nEdges)) 
        allocate(vertex_order(nVertices)) 
        
        do i=1,nCells
            cell_order(task_cell_start(cell_owner(i)+1)) = i
            task_cell_start(cell_owner(i)+1) = task_cell_start(cell_owner(i)+1) + 1
        end do

        call dictionary_reset()
        do i=1,nCells
            do j=1,nEdgesOnCell(i)
                k = edgesOnCell(j,i)
                if (cell_owner(cellsOnEdge(1,k)) == cell_owner(i)) then
                    if (.not. dictionary_search(k)) then
                        call dictionary_insert(k)
                        edge_order(task_edge_start(cell_owner(i)+1)) = k
                        task_edge_start(cell_owner(i)+1) = task_edge_start(cell_owner(i)+1) + 1
                    end if
                end if
            end do
        end do

        call dictionary_reset()
        do i=1,nCells
            do j=1,nEdgesOnCell(i)
                k = verticesOnCell(j,i)
                if (cell_owner(cellsOnVertex(1,k)) == cell_owner(i)) then
                    if (.not. dictionary_search(k)) then
                        call dictionary_insert(k)
                        vertex_order(task_vertex_start(cell_owner(i)+1)) = k
                        task_vertex_start(cell_owner(i)+1) = task_vertex_start(cell_owner(i)+1) + 1
                    end if
                end if
            end do
        end do

        deallocate(cellsOnVertex)
        deallocate(cellsOnEdge)

        deallocate(edgesOnCell)
        deallocate(verticesOnCell)
        deallocate(nEdgesOnCell)

        deallocate(cell_owner)
        deallocate(edge_owner)
        deallocate(vertex_owner)

        deallocate(task_cell_start)
        deallocate(task_edge_start)
        deallocate(task_vertex_start)


        !
        ! Compute reverse-lookup table to convert original indices to permuted indices
        ! (Where an old cell/edge/vertex maps to)
        !
        allocate(cell_lookup(2,nCells)) 
        allocate(cell_lookup_out(nCells)) 
        do i=1,nCells
            cell_lookup(1,i) = cell_order(i)
            cell_lookup(2,i) = i
        end do
        call mpas_quicksort_2dint(nCells, cell_lookup)
        cell_lookup_out(:) = cell_lookup(2,:)
        deallocate(cell_lookup)

        allocate(edge_lookup(2,nEdges)) 
        allocate(edge_lookup_out(nEdges)) 
        do i=1,nEdges
            edge_lookup(1,i) = edge_order(i)
            edge_lookup(2,i) = i
        end do
        call mpas_quicksort_2dint(nEdges, edge_lookup)
        edge_lookup_out(:) = edge_lookup(2,:)
        deallocate(edge_lookup)

        allocate(vertex_lookup(2,nVertices)) 
        allocate(vertex_lookup_out(nVertices)) 
        do i=1,nVertices
            vertex_lookup(1,i) = vertex_order(i)
            vertex_lookup(2,i) = i
        end do
        call mpas_quicksort_2dint(nVertices, vertex_lookup)
        vertex_lookup_out(:) = vertex_lookup(2,:)
        deallocate(vertex_lookup)

    end subroutine calc_permutation


    subroutine mpas_quicksort_2dint(nArray, array)

        implicit none
  
        integer, intent(in) :: nArray
        integer, dimension(2,nArray), intent(inout) :: array
  
        integer :: i, top, l, r, pivot, s
        integer :: pivot_value
        integer, dimension(2) :: temp
        integer, dimension(1000) :: lstack, rstack
  
        if (nArray < 1) return
  
        top = 1
        lstack(top) = 1
        rstack(top) = nArray
  
        do while (top > 0)
  
            l = lstack(top)
            r = rstack(top)
            top = top - 1
   
            pivot = (l+r)/2
   
            pivot_value = array(1,pivot)
            temp(:) = array(:,pivot)
            array(:,pivot) = array(:,r)
            array(:,r) = temp(:)
   
            s = l
            do i=l,r-1
                if (array(1,i) <= pivot_value) then
                    temp(:) = array(:,s)
                    array(:,s) = array(:,i)
                    array(:,i) = temp(:)
                    s = s + 1
                end if
            end do
   
            temp(:) = array(:,s)
            array(:,s) = array(:,r)
            array(:,r) = temp(:)
   
            if (s-1 > l) then
                top = top + 1
if (top > 1000) write(0,*) 'Error: Quicksort exhausted its stack.'
                lstack(top) = l
                rstack(top) = s-1
            end if
   
            if (r > s+1) then
                top = top + 1
if (top > 1000) write(0,*) 'Error: Quicksort exhausted its stack.'
                lstack(top) = s+1
                rstack(top) = r
            end if
        end do

    end subroutine mpas_quicksort_2dint


    subroutine xfer_dim(ncid_in, ncid_out, dimname)

        implicit none

        integer, intent(in) :: ncid_in, ncid_out
        character(len=*), intent(in) :: dimname

        integer :: dim_in, dim_out, ierr
        integer :: dimlen

#ifdef PNETCDF_INPUT
        ierr = nfmpi_inq_dimid(ncid_in, dimname, dim_in)
        ierr = nfmpi_inq_dimlen(ncid_in, dim_in, dimlen)
#endif
#ifdef NETCDF4_INPUT
        ierr = nf90_inq_dimid(ncid_in, dimname, dim_in)
        ierr = nf90_inquire_dimension(ncid_in, dim_in, len=dimlen)
#endif
#ifdef NETCDF4_OUTPUT
        ierr = nf90_def_dim(ncid_out, dimname, dimlen, dim_out)
#endif
#ifdef PNETCDF_OUTPUT
        ierr = nfmpi_def_dim(ncid_out, dimname, dimlen, dim_out)
#endif

    end subroutine xfer_dim


    subroutine xfer_att(ncid_in, ncid_out, attname)

        implicit none

        integer, intent(in) :: ncid_in, ncid_out
        character(len=*), intent(in) :: attname

        integer :: att_in, att_out, ierr
        integer :: attlen, atttype
        integer :: iatt
        double precision :: datt
        real :: ratt
        character(len=128) :: catt

#ifdef PNETCDF_INPUT
        ierr = nfmpi_inq_att(ncid_in, NF_GLOBAL, attname, atttype, attlen)
#endif
#ifdef NETCDF4_INPUT
        ierr = nf90_inquire_attribute(ncid_in, NF90_GLOBAL, attname, xtype=atttype, len=attlen)
#endif

        if (atttype == NF_INT) then
#ifdef PNETCDF_INPUT
            ierr = nfmpi_get_att_int(ncid_in, NF_GLOBAL, attname, iatt)
#endif
#ifdef NETCDF4_INPUT
            ierr = nf90_get_att(ncid_in, NF90_GLOBAL, attname, iatt)
#endif
#ifdef PNETCDF_OUTPUT
            ierr = nfmpi_put_att_int(ncid_out, NF_GLOBAL, attname, NF_INT, int(1,MPI_OFFSET_KIND), iatt)
#endif
#ifdef NETCDF4_OUTPUT
            ierr = nf90_put_att(ncid_out, NF90_GLOBAL, attname, iatt)
#endif
#ifdef PNETCDF_INPUT
        else if (atttype == NF_DOUBLE .or. atttype == NF_REAL) then
#endif
#ifdef NETCDF4_INPUT
        else if (atttype == NF90_DOUBLE .or. atttype == NF90_REAL) then
#endif
#ifdef PNETCDF_INPUT
#ifdef R4_INPUT
            ierr = nfmpi_get_att_real(ncid_in, NF_GLOBAL, attname, ratt)
            datt = ratt
#endif
#ifdef R8_INPUT
            ierr = nfmpi_get_att_double(ncid_in, NF_GLOBAL, attname, datt)
#endif
#endif
#ifdef NETCDF4_INPUT
#ifdef R4_INPUT
            ierr = nf90_get_att(ncid_in, NF90_GLOBAL, attname, ratt)
            datt = ratt
#endif
#ifdef R8_INPUT
            ierr = nf90_get_att(ncid_in, NF90_GLOBAL, attname, datt)
#endif
#endif
#ifdef NETCDF4_OUTPUT
#ifdef R4_OUTPUT
            ratt = datt
            ierr = nf90_put_att(ncid_out, NF_GLOBAL, attname, ratt)
#endif
#ifdef R8_OUTPUT
            ierr = nf90_put_att(ncid_out, NF_GLOBAL, attname, datt)
#endif
#endif
#ifdef PNETCDF_OUTPUT
#ifdef R4_OUTPUT
            ratt = datt
            ierr = nfmpi_put_att_real(ncid_out, NF_GLOBAL, attname, NF_REAL, int(1,MPI_OFFSET_KIND), ratt)
#endif
#ifdef R8_OUTPUT
            ierr = nfmpi_put_att_double(ncid_out, NF_GLOBAL, attname, NF_DOUBLE, int(1,MPI_OFFSET_KIND), datt)
#endif
#endif
        else if (atttype == NF_CHAR) then
#ifdef PNETCDF_INPUT
            ierr = nfmpi_get_att_text(ncid_in, NF_GLOBAL, attname, catt)
#endif
#ifdef NETCDF4_INPUT
            ierr = nf90_get_att(ncid_in, NF90_GLOBAL, attname, catt)
#endif
#ifdef NETCDF4_OUTPUT
            ierr = nf90_put_att(ncid_out, NF_GLOBAL, attname, catt)
#endif
#ifdef PNETCDF_OUTPUT
            ierr = nfmpi_put_att_text(ncid_out, NF_GLOBAL, attname, int(attlen,MPI_OFFSET_KIND), catt(1:attlen))
#endif
        end if

    end subroutine xfer_att


    subroutine xfer_var(ncid_in, ncid_out, varname)

        implicit none

        integer, intent(in) :: ncid_in, ncid_out
        character(len=*), intent(in) :: varname

        integer :: var_in, var_out, ierr, i
        integer :: vartype
        integer :: varndims
        integer, dimension(5) :: vardimids
        character (len=30), dimension(5) :: dimnames

#ifdef PNETCDF_INPUT
        ierr = nfmpi_inq_varid(ncid_in, varname, var_in)
        ierr = nfmpi_inq_vartype(ncid_in, var_in, vartype)
        ierr = nfmpi_inq_varndims(ncid_in, var_in, varndims)
        ierr = nfmpi_inq_vardimid(ncid_in, var_in, vardimids)
        do i=1,varndims
            ierr = nfmpi_inq_dimname(ncid_in, vardimids(i), dimnames(i))
        end do
#endif
#ifdef NETCDF4_INPUT
        ierr = nf90_inq_varid(ncid_in, varname, var_in)
        ierr = nf90_inquire_variable(ncid_in, var_in, xtype=vartype)
        ierr = nf90_inquire_variable(ncid_in, var_in, ndims=varndims)
        ierr = nf90_inquire_variable(ncid_in, var_in, dimids=vardimids)
        do i=1,varndims
            ierr = nf90_inquire_dimension(ncid_in, vardimids(i), name=dimnames(i))
        end do
#endif

#ifdef NETCDF4_OUTPUT
        do i=1,varndims
            ierr = nf90_inq_dimid(ncid_out, trim(dimnames(i)), vardimids(i))
        end do

#ifdef R4_OUTPUT
        if (vartype == NF90_DOUBLE) then
            vartype = NF90_REAL
        end if
#endif
#ifdef R8_OUTPUT
        if (vartype == NF90_REAL) then
            vartype = NF90_DOUBLE
        end if
#endif

        ierr = nf90_def_var(ncid_out, varname, vartype, vardimids(1:varndims), var_out)
#endif
#ifdef PNETCDF_OUTPUT
        do i=1,varndims
            ierr = nfmpi_inq_dimid(ncid_out, trim(dimnames(i)), vardimids(i))
        end do

#ifdef R4_OUTPUT
        if (vartype == NF_DOUBLE) then
            vartype = NF_REAL
        end if
#endif
#ifdef R8_OUTPUT
        if (vartype == NF_REAL) then
            vartype = NF_DOUBLE
        end if
#endif

        ierr = nfmpi_def_var(ncid_out, varname, vartype, varndims, vardimids(1:varndims), var_out)
#endif

    end subroutine xfer_var


    subroutine copy_var(ncid_in, ncid_out, varname, reindex, &
                        cell_order, cell_lookup, &
                        edge_order, edge_lookup, &
                        vertex_order, vertex_lookup)

        implicit none

        integer, intent(in) :: ncid_in, ncid_out
        character(len=*), intent(in) :: varname
        logical, intent(in) :: reindex
        integer, dimension(:), target, intent(in) :: cell_order
        integer, dimension(:), target, intent(in) :: cell_lookup
        integer, dimension(:), target, intent(in) :: edge_order
        integer, dimension(:), target, intent(in) :: edge_lookup
        integer, dimension(:), target, intent(in) :: vertex_order
        integer, dimension(:), target, intent(in) :: vertex_lookup

        integer :: var_in, var_out, ierr, i, j, k
        integer :: vartype
        integer :: varndims
        integer, dimension(5) :: vardimids
        integer, dimension(5) :: vardimlens
        integer, dimension(:), allocatable :: int1_in, int1_out
        integer, dimension(:,:), allocatable :: int2_in, int2_out
#ifdef R4_INPUT
        real, dimension(:), allocatable :: dbl1_in
        real, dimension(:,:), allocatable :: dbl2_in
        real, dimension(:,:,:), allocatable :: dbl3_in
#endif
#ifdef R8_INPUT
        double precision, dimension(:), allocatable :: dbl1_in
        double precision, dimension(:,:), allocatable :: dbl2_in
        double precision, dimension(:,:,:), allocatable :: dbl3_in
#endif
#ifdef R4_OUTPUT
        real, dimension(:), allocatable :: dbl1_out
        real, dimension(:,:), allocatable :: dbl2_out
        real, dimension(:,:,:), allocatable :: dbl3_out
#endif
#ifdef R8_OUTPUT
        double precision, dimension(:), allocatable :: dbl1_out
        double precision, dimension(:,:), allocatable :: dbl2_out
        double precision, dimension(:,:,:), allocatable :: dbl3_out
#endif
        character (len=:), allocatable :: char1_in 
        integer, dimension(:), pointer :: element_order
        integer, dimension(:), pointer :: element_lookup
        character (len=32) :: dimname 

        integer, dimension(1) :: istart1, icount1
        integer, dimension(2) :: istart2, icount2
        integer, dimension(3) :: istart3, icount3
#ifdef PNETCDF_INPUT
        integer(kind=MPI_OFFSET_KIND), dimension(1) :: start1, count1
        integer(kind=MPI_OFFSET_KIND), dimension(2) :: start2, count2
        integer(kind=MPI_OFFSET_KIND), dimension(3) :: start3, count3
#endif
#ifdef NETCDF4_INPUT
        integer, dimension(1) :: start1, count1
        integer, dimension(2) :: start2, count2
        integer, dimension(3) :: start3, count3
#endif

#ifdef PNETCDF_INPUT
        ierr = nfmpi_inq_varid(ncid_in, varname, var_in)
#endif
#ifdef NETCDF4_INPUT
        ierr = nf90_inq_varid(ncid_in, varname, var_in)
#endif
#ifdef NETCDF4_OUTPUT
        ierr = nf90_inq_varid(ncid_out, varname, var_out)
#endif
#ifdef PNETCDF_OUTPUT
        ierr = nfmpi_inq_varid(ncid_out, varname, var_out)
#endif

#ifdef PNETCDF_INPUT
        ierr = nfmpi_inq_vartype(ncid_in, var_in, vartype)
        ierr = nfmpi_inq_varndims(ncid_in, var_in, varndims)
        ierr = nfmpi_inq_vardimid(ncid_in, var_in, vardimids)
        do i=1,varndims
            ierr = nfmpi_inq_dimlen(ncid_in, vardimids(i), vardimlens(i))
        end do
        ierr = nfmpi_inq_dimname(ncid_in, vardimids(varndims), dimname)
#endif
#ifdef NETCDF4_INPUT
        ierr = nf90_inquire_variable(ncid_in, var_in, xtype=vartype)
        ierr = nf90_inquire_variable(ncid_in, var_in, ndims=varndims)
        ierr = nf90_inquire_variable(ncid_in, var_in, dimids=vardimids)
        do i=1,varndims
            ierr = nf90_inquire_dimension(ncid_in, vardimids(i), len=vardimlens(i))
        end do
        ierr = nf90_inquire_dimension(ncid_in, vardimids(varndims), name=dimname)
#endif

        if (trim(dimname) == 'nCells') then
            element_order => cell_order
        else if (trim(dimname) == 'nEdges') then
            element_order => edge_order
        else if (trim(dimname) == 'nVertices') then
            element_order => vertex_order
#ifdef NETCDF4_INPUT
        else if (vartype == NF90_INT .or. vartype == NF90_DOUBLE .or. vartype == NF90_REAL) then
#endif
#ifdef PNETCDF_INPUT
        else if (vartype == NF_INT .or. vartype == NF_DOUBLE .or. vartype == NF_REAL) then
#endif
            write(0,*) 'Error: No idea how to permute field with dimension '//trim(dimname)
            stop
        end if

        if (reindex) then
            if (varname(1:5) == 'cells') then
                element_lookup => cell_lookup
            else if (varname(1:5) == 'edges') then
                element_lookup => edge_lookup
            else if (varname(1:8) == 'vertices') then
                element_lookup => vertex_lookup
            else if (trim(varname) == 'indexToCellID') then
                element_lookup => cell_lookup
            else if (trim(varname) == 'indexToEdgeID') then
                element_lookup => edge_lookup
            else if (trim(varname) == 'indexToVertexID') then
                element_lookup => vertex_lookup
            else
                write(0,*) 'Error: No idea how to reindex field '//trim(varname)
                stop
            end if
        end if

        if (vartype == NF_INT) then

            if (varndims == 1) then

                allocate(int1_in(vardimlens(1)))
                allocate(int1_out(rw_chunk_size))

                start1(1) = 1
                count1(1) = min(rw_chunk_size, vardimlens(1))

                istart1(1) = start1(1)
                icount1(1) = count1(1)

                do while (istart1(1) <= vardimlens(1)) 
#ifdef PNETCDF_INPUT
                    ierr = nfmpi_get_vara_int(ncid_in, var_in, start1, count1, int1_in(start1(1):(start1(1)+count1(1)-1)))
#endif
#ifdef NETCDF4_INPUT
                    ierr = nf90_get_var(ncid_in, var_in, int1_in(start1(1):(start1(1)+count1(1)-1)), start=start1, count=count1)
#endif
                    if (ierr /= noerr_input) then
                        write(0,*) 'Error reading '//trim(varname)
                    end if

                    start1(1) = start1(1) + count1(1)
                    count1(1) = min(rw_chunk_size, vardimlens(1) - start1(1) + 1)

                    istart1(1) = start1(1)
                    icount1(1) = count1(1)
                end do

                start1(1) = 1
                count1(1) = min(rw_chunk_size, vardimlens(1))

                istart1(1) = start1(1)
                icount1(1) = count1(1)

                do while (istart1(1) <= vardimlens(1)) 

                    do i=start1(1),start1(1)+count1(1)-1
                        int1_out(i-start1(1)+1) = int1_in(element_order(i))
                        if (reindex) then
                            if (int1_out(i-start1(1)+1) /= 0) then
                                int1_out(i-start1(1)+1) = element_lookup(int1_out(i-start1(1)+1))
                            end if
                        end if
                    end do

#ifdef NETCDF4_OUTPUT
                    ierr = nf90_put_var(ncid_out, var_out, int1_out(1:count1(1)), start=istart1, count=icount1)
#endif
#ifdef PNETCDF_OUTPUT
                    ierr = nfmpi_put_vara_int(ncid_out, var_out, start1, count1, int1_out(1:count1(1)))
#endif
                    if (ierr /= noerr_output) then
                        write(0,*) 'Error writing '//trim(varname)
                    end if

                    start1(1) = start1(1) + count1(1)
                    count1(1) = min(rw_chunk_size, vardimlens(1) - start1(1) + 1)

                    istart1(1) = start1(1)
                    icount1(1) = count1(1)
                end do

                deallocate(int1_in)
                deallocate(int1_out)

            else if (varndims == 2) then

                allocate(int2_in(vardimlens(1),vardimlens(2)))
                allocate(int2_out(vardimlens(1),rw_chunk_size))

                start2(1) = 1
                count2(1) = vardimlens(1)
                start2(2) = 1
                count2(2) = min(rw_chunk_size, vardimlens(2))

                istart2(1) = start2(1)
                icount2(1) = count2(1)
                istart2(2) = start2(2)
                icount2(2) = count2(2)

                do while (istart2(2) <= vardimlens(2)) 
#ifdef PNETCDF_INPUT
                    ierr = nfmpi_get_vara_int(ncid_in, var_in, start2, count2, int2_in(:,start2(2):(start2(2)+count2(2)-1)))
#endif
#ifdef NETCDF4_INPUT
                    ierr = nf90_get_var(ncid_in, var_in, int2_in(:,start2(2):(start2(2)+count2(2)-1)), start=start2, count=count2)
#endif
                    if (ierr /= noerr_input) then
                        write(0,*) 'Error reading '//trim(varname)
                    end if

                    start2(2) = start2(2) + count2(2)
                    count2(2) = min(rw_chunk_size, vardimlens(2) - start2(2) + 1)

                    istart2(2) = start2(2)
                    icount2(2) = count2(2)
                end do

                start2(1) = 1
                count2(1) = vardimlens(1)
                start2(2) = 1
                count2(2) = min(rw_chunk_size, vardimlens(2))

                istart2(1) = start2(1)
                icount2(1) = count2(1)
                istart2(2) = start2(2)
                icount2(2) = count2(2)

                do while (istart2(2) <= vardimlens(2)) 

                    do i=start2(2),start2(2)+count2(2)-1
                    do j=1,vardimlens(1)
                        int2_out(j,i-start2(2)+1) = int2_in(j,element_order(i))
                        if (reindex) then
                            if (int2_out(j,i-start2(2)+1) /= 0) then
                                int2_out(j,i-start2(2)+1) = element_lookup(int2_out(j,i-start2(2)+1))
                            end if
                        end if
                    end do
                    end do

#ifdef NETCDF4_OUTPUT
                    ierr = nf90_put_var(ncid_out, var_out, int2_out(:,1:count2(2)), start=istart2, count=icount2)
#endif
#ifdef PNETCDF_OUTPUT
                    ierr = nfmpi_put_vara_int(ncid_out, var_out, start2, count2, int2_out(:,1:count2(2)))
#endif
                    if (ierr /= noerr_output) then
                        write(0,*) 'Error writing '//trim(varname)
                    end if

                    start2(2) = start2(2) + count2(2)
                    count2(2) = min(rw_chunk_size, vardimlens(2) - start2(2) + 1)

                    istart2(2) = start2(2)
                    icount2(2) = count2(2)
                end do

                deallocate(int2_in)
                deallocate(int2_out)
            end if

#ifdef NETCDF4_OUTPUT
        else if (vartype == NF90_DOUBLE .or. vartype == NF90_REAL) then
#endif
#ifdef PNETCDF_OUTPUT
        else if (vartype == NF_DOUBLE .or. vartype == NF_REAL) then
#endif

            if (varndims == 1) then

                allocate(dbl1_in(vardimlens(1)))
                allocate(dbl1_out(rw_chunk_size))

                start1(1) = 1
                count1(1) = min(rw_chunk_size, vardimlens(1))

                istart1(1) = start1(1)
                icount1(1) = count1(1)

                do while (istart1(1) <= vardimlens(1)) 
#ifdef PNETCDF_INPUT
#ifdef R4_INPUT
                    ierr = nfmpi_get_vara_real(ncid_in, var_in, start1, count1, dbl1_in(start1(1):(start1(1)+count1(1)-1)))
#endif
#ifdef R8_INPUT
                    ierr = nfmpi_get_vara_double(ncid_in, var_in, start1, count1, dbl1_in(start1(1):(start1(1)+count1(1)-1)))
#endif
#endif
#ifdef NETCDF4_INPUT
                    ierr = nf90_get_var(ncid_in, var_in, dbl1_in(start1(1):(start1(1)+count1(1)-1)), start=start1, count=count1)
#endif
                    if (ierr /= noerr_input) then
                        write(0,*) 'Error reading '//trim(varname)
                    end if

                    start1(1) = start1(1) + count1(1)
                    count1(1) = min(rw_chunk_size, vardimlens(1) - start1(1) + 1)

                    istart1(1) = start1(1)
                    icount1(1) = count1(1)
                end do

                start1(1) = 1
                count1(1) = min(rw_chunk_size, vardimlens(1))

                istart1(1) = start1(1)
                icount1(1) = count1(1)

                do while (istart1(1) <= vardimlens(1)) 

                    do i=start1(1),start1(1)+count1(1)-1
                        dbl1_out(i-start1(1)+1) = dbl1_in(element_order(i))
                    end do

#ifdef NETCDF4_OUTPUT
                    ierr = nf90_put_var(ncid_out, var_out, dbl1_out(1:count1(1)), start=istart1, count=icount1)
#endif
#ifdef PNETCDF_OUTPUT
#ifdef R4_OUTPUT
                    ierr = nfmpi_put_vara_real(ncid_out, var_out, start1, count1, dbl1_out(1:count1(1)))
#endif
#ifdef R8_OUTPUT
                    ierr = nfmpi_put_vara_double(ncid_out, var_out, start1, count1, dbl1_out(1:count1(1)))
#endif
#endif
                    if (ierr /= noerr_output) then
                        write(0,*) 'Error writing '//trim(varname)
                    end if

                    start1(1) = start1(1) + count1(1)
                    count1(1) = min(rw_chunk_size, vardimlens(1) - start1(1) + 1)

                    istart1(1) = start1(1)
                    icount1(1) = count1(1)
                end do

                deallocate(dbl1_in)
                deallocate(dbl1_out)

            else if (varndims == 2) then

                allocate(dbl2_in(vardimlens(1),vardimlens(2)))
                allocate(dbl2_out(vardimlens(1),rw_chunk_size))

                start2(1) = 1
                count2(1) = vardimlens(1)
                start2(2) = 1
                count2(2) = min(rw_chunk_size, vardimlens(2))

                istart2(1) = start2(1)
                icount2(1) = count2(1)
                istart2(2) = start2(2)
                icount2(2) = count2(2)

                do while (istart2(2) <= vardimlens(2)) 
#ifdef PNETCDF_INPUT
#ifdef R4_INPUT
                    ierr = nfmpi_get_vara_real(ncid_in, var_in, start2, count2, dbl2_in(:,start2(2):(start2(2)+count2(2)-1)))
#endif
#ifdef R8_INPUT
                    ierr = nfmpi_get_vara_double(ncid_in, var_in, start2, count2, dbl2_in(:,start2(2):(start2(2)+count2(2)-1)))
#endif
#endif
#ifdef NETCDF4_INPUT
                    ierr = nf90_get_var(ncid_in, var_in, dbl2_in(:,start2(2):(start2(2)+count2(2)-1)), start=start2, count=count2)
#endif
                    if (ierr /= noerr_input) then
                        write(0,*) 'Error reading '//trim(varname)
                    end if

                    start2(2) = start2(2) + count2(2)
                    count2(2) = min(rw_chunk_size, vardimlens(2) - start2(2) + 1)

                    istart2(2) = start2(2)
                    icount2(2) = count2(2)
                end do

                start2(1) = 1
                count2(1) = vardimlens(1)
                start2(2) = 1
                count2(2) = min(rw_chunk_size, vardimlens(2))

                istart2(1) = start2(1)
                icount2(1) = count2(1)
                istart2(2) = start2(2)
                icount2(2) = count2(2)

                do while (istart2(2) <= vardimlens(2)) 

                    do i=start2(2),start2(2)+count2(2)-1
                    do j=1,vardimlens(1)
                        dbl2_out(j,i-start2(2)+1) = dbl2_in(j,element_order(i))
                    end do
                    end do

#ifdef NETCDF4_OUTPUT
                    ierr = nf90_put_var(ncid_out, var_out, dbl2_out(:,1:count2(2)), start=istart2, count=icount2)
#endif
#ifdef PNETCDF_OUTPUT
#ifdef R4_OUTPUT
                    ierr = nfmpi_put_vara_real(ncid_out, var_out, start2, count2, dbl2_out(:,1:count2(2)))
#endif
#ifdef R8_OUTPUT
                    ierr = nfmpi_put_vara_double(ncid_out, var_out, start2, count2, dbl2_out(:,1:count2(2)))
#endif
#endif
                    if (ierr /= noerr_output) then
                        write(0,*) 'Error writing '//trim(varname)
                    end if

                    start2(2) = start2(2) + count2(2)
                    count2(2) = min(rw_chunk_size, vardimlens(2) - start2(2) + 1)

                    istart2(2) = start2(2)
                    icount2(2) = count2(2)
                end do

                deallocate(dbl2_in)
                deallocate(dbl2_out)

            else if (varndims == 3) then

                allocate(dbl3_in(vardimlens(1),vardimlens(2),vardimlens(3)))
                allocate(dbl3_out(vardimlens(1),vardimlens(2),rw_chunk_size))

                start3(1) = 1
                count3(1) = vardimlens(1)
                start3(2) = 1
                count3(2) = vardimlens(2)
                start3(3) = 1
                count3(3) = min(rw_chunk_size, vardimlens(3))

                istart3(1) = start3(1)
                icount3(1) = count3(1)
                istart3(2) = start3(2)
                icount3(2) = count3(2)
                istart3(3) = start3(3)
                icount3(3) = count3(3)

                do while (istart3(3) <= vardimlens(3)) 
#ifdef PNETCDF_INPUT
#ifdef R4_INPUT
                    ierr = nfmpi_get_vara_real(ncid_in, var_in, start3, count3, dbl3_in(:,:,start3(3):(start3(3)+count3(3)-1)))
#endif
#ifdef R8_INPUT
                    ierr = nfmpi_get_vara_double(ncid_in, var_in, start3, count3, dbl3_in(:,:,start3(3):(start3(3)+count3(3)-1)))
#endif
#endif
#ifdef NETCDF4_INPUT
                    ierr = nf90_get_var(ncid_in, var_in, dbl3_in(:,:,start3(3):(start3(3)+count3(3)-1)), start=start3, count=count3)
#endif
                    if (ierr /= noerr_input) then
                        write(0,*) 'Error reading '//trim(varname)
                        write(0,*) ierr
                        return
                    end if

                    start3(3) = start3(3) + count3(3)
                    count3(3) = min(rw_chunk_size, vardimlens(3) - start3(3) + 1)

                    istart3(3) = start3(3)
                    icount3(3) = count3(3)
                end do

                start3(1) = 1
                count3(1) = vardimlens(1)
                start3(2) = 1
                count3(2) = vardimlens(2)
                start3(3) = 1
                count3(3) = min(rw_chunk_size, vardimlens(3))

                istart3(1) = start3(1)
                icount3(1) = count3(1)
                istart3(2) = start3(2)
                icount3(2) = count3(2)
                istart3(3) = start3(3)
                icount3(3) = count3(3)

                do while (istart3(3) <= vardimlens(3)) 

                    do i=start3(3),start3(3)+count3(3)-1
                    do j=1,vardimlens(2)
                    do k=1,vardimlens(1)
                        dbl3_out(k,j,i-start3(3)+1) = dbl3_in(k,j,element_order(i))
                    end do
                    end do
                    end do

#ifdef NETCDF4_OUTPUT
                    ierr = nf90_put_var(ncid_out, var_out, dbl3_out(:,:,1:count3(3)), start=istart3, count=icount3)
#endif
#ifdef PNETCDF_OUTPUT
#ifdef R4_OUTPUT
                    ierr = nfmpi_put_vara_real(ncid_out, var_out, start3, count3, dbl3_out(:,:,1:count3(3)))
#endif
#ifdef R8_OUTPUT
                    ierr = nfmpi_put_vara_double(ncid_out, var_out, start3, count3, dbl3_out(:,:,1:count3(3)))
#endif
#endif
                    if (ierr /= noerr_output) then
                        write(0,*) 'Error writing '//trim(varname)
                        write(0,*) ierr
                        return
                    end if

                    start3(3) = start3(3) + count3(3)
                    count3(3) = min(rw_chunk_size, vardimlens(3) - start3(3) + 1)

                    istart3(3) = start3(3)
                    icount3(3) = count3(3)
                end do

                deallocate(dbl3_in)
                deallocate(dbl3_out)
            else
                write(0,*) '*** cannot transfer real field '//trim(varname)//'***'
            end if

        else if (vartype == NF_CHAR) then

            if (varndims == 1) then

                allocate(character (len=vardimlens(1)) :: char1_in)
#ifdef PNETCDF_INPUT
                ierr = nfmpi_get_var_text(ncid_in, var_in, char1_in)
#endif
#ifdef NETCDF4_INPUT
                ierr = nf90_get_var(ncid_in, var_in, char1_in)
#endif
                if (ierr /= noerr_input) then
                    write(0,*) 'Error reading '//trim(varname)
                    write(0,*) ierr
                    return
                end if
#ifdef NETCDF4_OUTPUT
                ierr = nf90_put_var(ncid_out, var_out, char1_in)
#endif
#ifdef PNETCDF_OUTPUT
                ierr = nfmpi_put_var_text(ncid_out, var_out, char1_in)
#endif
                if (ierr /= noerr_output) then
                    write(0,*) 'Error writing '//trim(varname)
                    write(0,*) ierr
                    return
                end if
                deallocate(char1_in)
            end if
        end if

    end subroutine copy_var


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Insert an integer into the dictionary
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dictionary_insert(i)
    
        implicit none
    
        integer, intent(in) :: i
    
        integer :: n_integer
        integer :: n_bit
    
        n_integer = ((i-1) / int_size) + 1
        n_bit = mod((i-1), int_size)
    
        if (n_integer > max_dictionary_size) then
            write(0,*) 'Error: dictionary insert out of bounds'
            return
        end if
    
        dictionary_array(n_integer) = ibset(dictionary_array(n_integer), n_bit)
    
    end subroutine dictionary_insert
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Search for an integer in the dictionary
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical function dictionary_search(i)
    
        implicit none
    
        integer, intent(in) :: i
    
        integer :: n_integer
        integer :: n_bit
    
        n_integer = ((i-1) / int_size) + 1
        n_bit = mod((i-1), int_size)
    
        if (n_integer > max_dictionary_size) then
            write(0,*) 'Error: dictionary search out of bounds'
            dictionary_search = .false.
            return
        end if
    
        dictionary_search = btest(dictionary_array(n_integer), n_bit)
    
    end function dictionary_search
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Reset the dictionary to an empty state
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dictionary_reset()
    
        implicit none
    
        dictionary_array(:) = 0
    
    end subroutine dictionary_reset


end program permute_mesh
