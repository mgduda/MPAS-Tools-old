module scvt

   use data_types
   use sphere_utilities
   use voronoi_utils
   use grid_constants
   use grid_params
   use mgd


   contains


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  SUBROUTINE SCVT_SOLVE
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine scvt_solve(n, lend, x, y, z, nvc, list, lptr, fn)
   
      implicit none

      integer, intent(in) :: n, nvc, fn
      integer, dimension(n), intent(inout) :: lend
      integer, dimension(nvc), intent(inout) :: list, lptr 
      real, dimension(n), intent(inout) :: x, y, z

      integer :: maxitr

      integer :: i, j, k, iter
      integer :: ntmax, nrow, nptri
      integer, allocatable, dimension(:) :: listc
      real :: area, density, tot_mass
      real :: xm, ym, zm, new_ctr_x, new_ctr_y, new_ctr_z 
      real, allocatable, dimension(:) :: x_2, y_2, z_2, xc, yc, zc
      real, dimension(3) :: p1_cart, p2_cart, p3_cart, pc_cart
      real, dimension(3,3,256) :: ptri_cart
      real :: m
      real :: avg_movement, maxmovement, movement
      logical converged
integer :: itimer, sec, nsec

      maxitr = n_scvt_iterations

      maxmovement = 100000


      ntmax = 6*n
      nrow = 6
      nptri = 256
      allocate(listc(nvc))
      allocate(xc(nvc))
      allocate(yc(nvc))
      allocate(zc(nvc))
      allocate(x_2(n))
      allocate(y_2(n))
      allocate(z_2(n))

      iter = 1
      converged = .false.

      call mgd_init()

      do while (iter <= maxitr .and. .not.converged)

!         do i=1,n
!            m = x(i) * x(i) + y(i) * y(i) + z(i) * z(i)
!            m = sqrt(m)
!            x(i) = x(i) / m
!            y(i) = y(i) / m
!            z(i) = z(i) / m
!         end do

!         write(0,*) 'scvt iteration ',iter

         !
         ! Compute Voronoi corners
         !
itimer = 1
call start_timer(itimer)
         call compute_vc_mgd(x, y, z, n, nrow, ntmax, list, lptr, lend, listc, xc, yc, zc, nvc)
call stop_timer(itimer, sec, nsec)
if (mod(iter,report_interval) == 0) then
   write(0,'(a,i3,a1,i9.9)') 'compute_vc: ', sec, '.', nsec
end if

   
         !
         ! Loop over vertices
         !   Within the loop, p0 always refers to the current vertex being processed
         !
call start_timer(itimer)
!$OMP PARALLEL DO PRIVATE(I, J, K, NEW_CTR_X, NEW_CTR_Y, NEW_CTR_Z, TOT_MASS, P1_CART, P2_CART, P3_CART, AREA, DENSITY, XM, YM, ZM, PTRI_CART) SHARED(X, Y, Z, X_2, Y_2, Z_2, LPTR, LEND, LISTC, XC, YC, ZC, NPTRI)
         do i=1,n

            new_ctr_x = 0.0
            new_ctr_y = 0.0
            new_ctr_z = 0.0
            tot_mass = 0.0

   
            !
            ! Compute center of mass of Voronoi cell
            !
            p1_cart(1) = x(i)
            p1_cart(2) = y(i)
            p1_cart(3) = z(i)
   
            k = lend(i)
            p2_cart(1) = xc(listc(k))
            p2_cart(2) = yc(listc(k))
            p2_cart(3) = zc(listc(k))
            k = lptr(lend(i))
            p3_cart(1) = xc(listc(k))
            p3_cart(2) = yc(listc(k))
            p3_cart(3) = zc(listc(k))

            pc_cart(:) = abs(p3_cart(:) - p2_cart(:))
   
            if (sum(pc_cart) > 0.00001) then
               call divide_triangle_mgd(p1_cart, p2_cart, p3_cart, nptri, ptri_cart)
               do j=1,nptri
                  area = triangle_area_mgd(ptri_cart(:,1,j), ptri_cart(:,2,j), ptri_cart(:,3,j))
                  call center_of_mass_mgd(ptri_cart(:,1,j), ptri_cart(:,2,j), ptri_cart(:,3,j), xm, ym, zm)
                  density = density_for_point_mgd(xm, ym, zm)
                  tot_mass = tot_mass + area * density

                  new_ctr_x = new_ctr_x + xm * area * density
                  new_ctr_y = new_ctr_y + ym * area * density
                  new_ctr_z = new_ctr_z + zm * area * density
               end do
            end if
   
            do while (k /= lend(i))
               k = lptr(k)
               p2_cart(:) = p3_cart(:)

               p3_cart(1) = xc(listc(k))
               p3_cart(2) = yc(listc(k))
               p3_cart(3) = zc(listc(k))

               pc_cart(:) = abs(p3_cart(:) - p2_cart(:))
   
               if (sum(pc_cart) > 0.00001) then
                  call divide_triangle_mgd(p1_cart, p2_cart, p3_cart, nptri, ptri_cart)
                  do j=1,nptri
                     area = triangle_area_mgd(ptri_cart(:,1,j), ptri_cart(:,2,j), ptri_cart(:,3,j))
                     call center_of_mass_mgd(ptri_cart(:,1,j), ptri_cart(:,2,j), ptri_cart(:,3,j), xm, ym, zm)
                     density = density_for_point_mgd(xm, ym, zm)
                     tot_mass = tot_mass + area * density

                     new_ctr_x = new_ctr_x + xm * area * density
                     new_ctr_y = new_ctr_y + ym * area * density
                     new_ctr_z = new_ctr_z + zm * area * density
                  end do
               end if
            end do

            if (tot_mass > 0.0) then
               x_2(i) = new_ctr_x / tot_mass
               y_2(i) = new_ctr_y / tot_mass
               z_2(i) = new_ctr_z / tot_mass
            else
               x_2(i) = x(i)
               y_2(i) = y(i)
               z_2(i) = z(i)
            end if

            m = x_2(i) * x_2(i) + y_2(i) * y_2(i) + z_2(i) * z_2(i)
            m = sqrt(m)
            x_2(i) = x_2(i) / m
            y_2(i) = y_2(i) / m
            z_2(i) = z_2(i) / m

         end do
!$OMP END PARALLEL DO
itimer = 1
call stop_timer(itimer, sec, nsec)
if (mod(iter,report_interval) == 0) then
   write(0,'(a,i3,a1,i9.9)') 'centroids: ', sec, '.', nsec
end if


         !Compute movement
         if (mod(iter,report_interval) == 0) then
           maxmovement = 0.0
           avg_movement = 0.0
           do i=1,n
            
              movement = sqrt((x(i) - x_2(i))**2 + (y(i) - y_2(i))**2 + (z(i) - z_2(i))**2)
  
              if(movement > maxmovement) maxmovement = movement
              avg_movement = avg_movement + movement
  
           end do
           avg_movement = avg_movement / n
           if (avg_movement < eps .and. l2_conv) converged=.true.
           if (avg_movement < eps .and. inf_conv) converged=.true.
           write(0,*) n, iter, maxmovement, avg_movement
         end if

         x(:) = x_2(:)
         y(:) = y_2(:)
         z(:) = z_2(:)
         iter = iter + 1

         if (mod(iter,restart_interval) == 0) then
            write(0,*) 'Writing restart file'
            open(44,file='restart.dat',form='formatted',status='unknown')
            write(44,*) n
            do i=1,n
               write(44,'(f18.15,1x,f18.15,1x,f18.15)') x(i), y(i), z(i)
            end do
            close(44)
         end if
    
      end do


      open(89,file='density.dat',form='formatted',status='replace')
      do i=1,n
         density = density_for_point_mgd(x(i), y(i), z(i))
         write(89,'(f18.16)') density
      end do
      close(89)


      deallocate(listc)
      deallocate(xc)
      deallocate(yc)
      deallocate(zc)
      deallocate(x_2)
      deallocate(y_2)
      deallocate(z_2)

      if (maxitr > 0) write(0,*) 'Finished SCVT solve'

   end subroutine scvt_solve


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  SUBROUTINE RANDOM_POINT
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine random_point(p)

      type (geo_point), intent(inout) :: p
      real :: x, y, z, m
      real :: pi

      pi = 4.0*atan(1.0)

      x = 0.0
      y = 0.0
      z = 0.0
      m = 2.0

      do while (m > 1.0 .or. (x == 0.0 .and. y == 0.0 .and. z == 0.0))
         call random_number(x)
         call random_number(y)
         call random_number(z)
         x = x * 2.0 - 1.0
         y = y * 2.0 - 1.0
         z = z * 2.0 - 1.0
         m = x**2 + y**2 + z**2
      end do

      m = 1.0 / sqrt(m)
      x = x * m
      y = y * m
      z = z * m

      call convert_xl(x, y, z, p)

   end subroutine random_point


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  FUNCTION DENSITY_FOR_POINT
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real function density_for_point(p)

      implicit none

      type (geo_point), intent(in) :: p

      character (len=256) :: fname
      real :: rx, ry, rz, prx, pry, prz
      type (geo_point) :: p_local
      real :: hgt
      real :: r, norm, t_cent
      real :: r1
      real :: pi
      real :: width, trans_center, min_val

      pi = 4.0*atan(1.0)

      !density_for_point = 1.0 + (1.19*cos(p%lat-3.141592654/4.0))**16.0

      ! Uniform Density Function
!      density_for_point = 1.0

      !Target Density Function based on hyperbolic tangent
      !		p_local%lat = latitude (radians) center of high-resolution region
      !         p_local%lon = longitude (radians) center of high-resolution region
      !         width = width of transition zone 
      !         trans_center = width (radians) of high resolution zone
      !         minval = minimum density value. to have grid spacing vary by a factor of 8
      !                  set minval = (1.0 / 8.0)**4

      call convert_lx(rx, ry, rz, 1.0, p)
!      call convert_lx(prx, pry, prz, 1.0, p_local)
      prx = 0.0
      pry = 0.0
      prz = 1.0
      r = acos(rx*prx + ry*pry + rz*prz)

      width = 0.40
      trans_center = pi/7.0
      min_val = (1.0/5.0)**4
      norm = 1.0/(1.0-min_val)
      density_for_point = ((tanh((trans_center-r)*(1.0/width))+1.0)/2)/norm + min_val

   end function density_for_point

end module scvt
