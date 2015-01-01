module voronoi_utils

   use grid_constants
   use stripack

   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! SUBROUTINE COMPUTE_DT
   !
   ! Compute the Delaunay triangulation of a set of lat/lon locations.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   subroutine compute_dt(x, y, z, n, ltri, nrow, ntmx, nt)

      implicit none

      integer, intent(in) :: n, nrow, ntmx
      integer, intent(inout) :: nt
      integer, dimension(nrow, ntmx), intent(in) :: ltri
      real, dimension(n), intent(in) :: x, y, z

      integer :: ierr, lnew
      integer, dimension(n) :: near, next
      integer, dimension(n) :: lend
      integer, dimension(6*n+12) :: list, lptr
      real, dimension(n) :: dist
      integer :: i, prev_j, j

      write(0,*) 'started TRMESH'
      call trmesh(n, x, y, z, list, lptr, lend, lnew, near, next, dist, ierr)
      if (ierr /= 0) then
         write(0,*) 'Error: TRMESH returned error code ',ierr
      end if
      write(0,*) 'finished TRMESH ', lnew

      open(88,file='triangles.dat',form='formatted',status='replace')
      do i=1,n
         prev_j = lend(i)
         j = lptr(prev_j)
         do while (j /= lend(i))
            if (i < list(j) .and. i < list(prev_j)) write(88,*) i, list(j), list(prev_j)
            prev_j = j
            j = lptr(j)
         end do
         if (i < list(j) .and. i < list(prev_j)) write(88,*) i, list(j), list(prev_j)
      end do
      close(88)
      write(0,*) 'finished LISTING TRIANGLES '

   end subroutine compute_dt


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! SUBROUTINE COMPUTE_VC
   !
   ! Compute the Voronoi corners of a set of lat/lon locations.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   subroutine compute_vc(rlat, rlon, n, nrow, ntmx, list, lptr, lend, listc, vclat, vclon, nvc)

      implicit none

      integer, intent(in) :: n, nrow, ntmx, nvc
      integer, dimension(nvc), intent(inout) :: list, lptr, listc
      real, dimension(nvc), intent(inout) :: vclat, vclon
      integer, dimension(n), intent(inout) :: lend
      real, dimension(n), intent(in) :: rlat, rlon

      integer :: ierr, lnew, nb
      integer, dimension(n) :: near, next
      integer, dimension(nrow, ntmx) :: ltri
      real, dimension(n) :: x, y, z, dist
      real, dimension(nvc) :: xc, yc, zc, rc

      if (nvc < 6*n-12) then
         write(0,*) 'Error: Argument nvc to COMPUTE_VC must be at least 6*n+12'
         return
      end if

      call trans(n, rlat, rlon, x, y, z)

      call trmesh(n, x, y, z, list, lptr, lend, lnew, near, next, dist, ierr)
      if (ierr /= 0) then
         write(0,*) 'Error: TRMESH returned error code ',ierr
      end if

      call crlist(n, ntmx, x, y, z, list, lend, lptr, lnew, ltri, listc, nb, xc, yc, zc, rc, ierr)
      if (ierr /= 0) then
         write(0,*) 'Error: CRLIST returned error code ',ierr
      end if

      call trans_inv(xc, yc, zc, vclat, vclon, nvc)

   end subroutine compute_vc


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! SUBROUTINE TRANS_INV
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   subroutine trans_inv(x, y, z, lat, lon, n)

      implicit none

      integer, intent(in) :: n
      real, dimension(n), intent(in) :: x, y, z
      real, dimension(n), intent(out) :: lat, lon

      integer :: i

      do i=1,n
         lat(i) = (pii/2.0 - acos(z(i)))
         lon(i) = atan2(y(i),x(i))
      end do

   end subroutine trans_inv

end module voronoi_utils
