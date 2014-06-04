module mgd

   use stripack, only : trmesh, crlist

   private

   real :: pii, width, trans_center, min_val, norm

   public :: mgd_init, density_for_point_mgd, center_of_mass_mgd, triangle_area_mgd, divide_triangle_mgd, compute_vc_mgd, convert_xl_mgd, convert_lx_mgd


   contains

   
   subroutine mgd_init()

      implicit none

      pii = 2.0*asin(1.0)
      width = 0.10
      width = 1.0 / width
      trans_center = pii/7.0
      min_val = (1.0/5.0)**4
      norm = (1.0-min_val)

   end subroutine mgd_init


   real function density_for_point_mgd(x, y, z)

      implicit none

      real, intent(in) :: x, y, z

      real :: r

      r = acos(z)

!      density_for_point_mgd = ((tanh((trans_center-r)*(width))+1.0)*0.5)*norm + min_val
      density_for_point_mgd = 1.0

   end function density_for_point_mgd


   subroutine center_of_mass_mgd(p1, p2, p3, xc, yc, zc)

      implicit none

      real, dimension(3), intent(in) :: p1, p2, p3
      real, intent(out) :: xc, yc, zc

      real :: x1, x2, x3
      real :: y1, y2, y3
      real :: z1, z2, z3

      x1 = p1(1)
      x2 = p2(1)
      x3 = p3(1)

      y1 = p1(2)
      y2 = p2(2)
      y3 = p3(2)

      z1 = p1(3)
      z2 = p2(3)
      z3 = p3(3)

      xc = (x1+x2+x3)/3.
      yc = (y1+y2+y3)/3.
      zc = (z1+z2+z3)/3.

   end subroutine center_of_mass_mgd


   real function triangle_area_mgd(p1, p2, p3)

      implicit none

      real, dimension(3), intent(in) :: p1, p2, p3

      real :: a, b, c, s, tanqe

      a = sphere_distance_mgd(p1,p2)
      b = sphere_distance_mgd(p2,p3)
      c = sphere_distance_mgd(p3,p1)
      s = 0.5*(a+b+c)

      tanqe = sqrt(tan(0.5*s)*tan(0.5*(s-a))*tan(0.5*(s-b))*tan(0.5*(s-c)))
      triangle_area_mgd = 4.0*atan(tanqe)

   end function triangle_area_mgd


   real function sphere_distance_mgd(p1, p2)

      implicit none

      real, dimension(3), intent(in) :: p1, p2

      sphere_distance_mgd = acos(p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3))

   end function sphere_distance_mgd


   subroutine divide_triangle_mgd(p1, p2, p3, n, p)
   
      implicit none
   
      real, dimension(3), intent(in) :: p1, p2, p3
      integer, intent(in) :: n
      real, dimension(3,3,n), intent(out) :: p
   
      integer :: i, j, k
      integer :: glevel      ! Level of decomposition
      real, allocatable, dimension(:,:) :: p1p2, p1p3
      real, allocatable, dimension(:,:,:) :: line
      real :: x, y, z
   
      glevel = nint(log(real(n)) / log(4.0))  ! Each subdivision gives four times the number of 
                                              !   triangles, so log4(n) gives the level decomposition
   
      glevel = (2 ** glevel) + 1
      allocate(line(3,glevel,glevel))
      allocate(p1p2(3,glevel))
      allocate(p1p3(3,glevel))

      call great_circle_points_mgd(p1, p2, p1p2, glevel)
      call great_circle_points_mgd(p1, p3, p1p3, glevel)
      
      line(:,1,1) = p1(:)
      line(:,1,2) = p1p2(:,2)
      line(:,2,2) = p1p3(:,2)
   
      do i = 3,glevel
         call great_circle_points_mgd(p1p2(:,i), p1p3(:,i), line(:,:,i), i)
      end do
      
      k = 1
      do i = 1,glevel-1
         do j = 1,i
            p(:,1,k) = line(:,j,i)
            p(:,2,k) = line(:,j,i+1)
            p(:,3,k) = line(:,j+1,i+1)
            k = k + 1
         end do
      end do
   
      do i = glevel,3,-1
         do j = 2,i-1
            p(:,1,k) = line(:,j,i)
            p(:,2,k) = line(:,j,i-1)
            p(:,3,k) = line(:,j-1,i-1)
            k = k + 1
         end do
      end do
   
      deallocate(line)
      deallocate(p1p2)
      deallocate(p1p3)
   
   end subroutine divide_triangle_mgd


   subroutine great_circle_points_mgd(p1, p2, pl, n)
   
      implicit none
   
      integer, intent(in) :: n
      real, dimension(3), intent(in) :: p1, p2
      real, dimension(3,n), intent(inout) :: pl
   
      real :: dx, dl
      real :: x, y, z
      integer :: i
      real :: dtheta, dinc, dt
   
      if (n == 2) then
         pl(:,1) = p1(:)
         pl(:,2) = p2(:)
      end if
   
      dtheta = sphere_distance_mgd(p1, p2)
      dinc = dtheta / real(n-1)
   
   !  set the end points
   
      pl(:,1) = p1(:)
      pl(:,n) = p2(:)
   
   
   ! compute the interior points.  see notes for derivation
   
      do i=2,n-1
         dt = real(i-1)*dinc
   
         if (dt <= 0.5*dtheta) then
            dx = 1.-tan(0.5*dtheta-dt)/tan(0.5*dtheta)
            pl(1,i) = p1(1)+0.5*dx*(p2(1)-p1(1))
            pl(2,i) = p1(2)+0.5*dx*(p2(2)-p1(2))
            pl(3,i) = p1(3)+0.5*dx*(p2(3)-p1(3))
         else
            dt = dtheta-dt
            dx = 1.-tan(0.5*dtheta-dt)/tan(0.5*dtheta)
            pl(1,i) = p2(1)+0.5*dx*(p1(1)-p2(1))
            pl(2,i) = p2(2)+0.5*dx*(p1(2)-p2(2))
            pl(3,i) = p2(3)+0.5*dx*(p1(3)-p2(3))
         end if

      end do
   
   end subroutine great_circle_points_mgd


   subroutine compute_vc_mgd(x, y, z, n, nrow, ntmx, list, lptr, lend, listc, xc, yc, zc, nvc)

      implicit none

      integer, intent(in) :: n, nrow, ntmx, nvc
      integer, dimension(nvc), intent(inout) :: list, lptr, listc
      real, dimension(nvc), intent(inout) :: xc, yc, zc
      integer, dimension(n), intent(inout) :: lend
      real, dimension(n), intent(in) :: x, y, z

      integer :: ierr, lnew, nb
      integer, dimension(n) :: near, next
      integer, dimension(nrow, ntmx) :: ltri
      real, dimension(n) :: dist
      real, dimension(nvc) :: rc

      call trmesh(n, x, y, z, list, lptr, lend, lnew, near, next, dist, ierr)
      if (ierr /= 0) then
         write(0,*) 'Error: TRMESH returned error code ',ierr
      end if

      call crlist(n, ntmx, x, y, z, list, lend, lptr, lnew, ltri, listc, nb, xc, yc, zc, rc, ierr)
      if (ierr /= 0) then
         write(0,*) 'Error: CRLIST returned error code ',ierr
      end if

   end subroutine compute_vc_mgd


   subroutine convert_lx_mgd(lat, lon, x, y, z)

      implicit none

      real, intent(in) :: lat, lon
      real, intent(out) :: x, y, z

      z = sin(lat)
      x = cos(lon) * cos(lat)
      y = sin(lon) * cos(lat)

   end subroutine convert_lx_mgd


   subroutine convert_xl_mgd(x, y, z, lat, lon)
   
      implicit none
   
      real, intent(in) :: x, y, z
      real, intent(out) :: lat, lon
   
      real :: dl, clat, pii, rtod
      real :: eps
      parameter (eps=1.e-10)
   
      pii = 2.*asin(1.0)
      rtod=180./pii
      dl = sqrt(x*x + y*y + z*z)
   
      lat = asin(z/dl)
   
   !  check for being close to either pole
   
      if (abs(x) > eps) then
   
         if (abs(y) > eps) then
   
            lon = atan(abs(y/x))
   
            if ((x <= 0.) .and. (y >= 0.)) then
               lon = pii-lon
            else if ((x <= 0.) .and. (y < 0.)) then
               lon = lon+pii
            else if ((x >= 0.) .and. (y <= 0.)) then
               lon = 2*pii-lon
            end if
   
         else ! we're either on longitude 0 or 180
   
            if (x > 0) then
               lon = 0.
            else
               lon = pii
            end if
   
         end if
   
      else if (abs(y) > eps) then
   
         if (y > 0) then
            lon = pii/2.
         else
            lon = 3.*pii/2.
         end if
   
      else  ! we are at a pole
   
         lon = 0.
   
      end if
   
   end subroutine convert_xl_mgd

end module mgd
