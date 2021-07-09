subroutine evaluation_points
!calculation of evaluation points
!T.Bisbas
use definitions
use healpix_types
use healpix_module
use maincode_module
#ifdef OPENMP
use omp_lib
#endif

double precision::adaptivemin
logical::killray(0:nrays-1)
integer::j
real(kind=dp) :: psi_z, psi_x !angles of rotations for matching UV-source's direction and healpix ray with coordinates (0 1 0)
real(kind=dp) :: vel_vec(1:3),temp_vec(1:3)
real(kind=dp) :: test_angle,test_angle2, test_module, test_module2

if (dark_ptot.gt.0) then

write(6,*) 'Creating evaluation points for the Dark Molecular element.'
!BUILDING HEALPIX VECTORS FOR THE ONE DARK MOLECULAR ELEMENT ======================================
!SERIAL PROCESS -----------------
!defines the origin to transfer all the domain in the original co-ordinate system
killray=.false.
origin(1:3) = pdrpoint(1:3,0)

IF (fieldchoice.EQ."PNT") THEN
call calc_rotation_angles(UV_source_coord-origin,psi_z,psi_x)
ENDIF

allocate(ra(0:grand_ptot-1)) !needs one extra place for sorting in heapsort
allocate(rb(1:grand_ptot-1)) !-1 to avoid overlapping origin & pdrpoint
allocate(ep(1:3,0:nrays-1))
!calculating distances from the origin(1:3)
kk=0
pdr(IDlist_dark(1))%epray = 0
do i=1,grand_ptot
   if (i.eq.IDlist_dark(1)) cycle
   kk=kk+1
   !locates the grid point in the new computational domain
   rvec(1)=pdr(i)%x-origin(1)
   rvec(2)=pdr(i)%y-origin(2)
   rvec(3)=pdr(i)%z-origin(3)

IF (fieldchoice.EQ."PNT") THEN
   call rotate_z(rvec(1),rvec(2),rvec(3),psi_z)
   call rotate_x(rvec(1),rvec(2),rvec(3),psi_x)
ENDIF

   !next two lines return the ipix ray that the rvec(1:3) point belongs to.
   call vec2ang(rvec,theta,phi)
   call ang2pix_nest_id(nside,theta,phi,ipix)
   ra(kk)=sqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)
   rb(kk)=i !stores the identifier of each grid point
enddo

ktot=kk !ktot should be grand_ptot-1
if (ktot.ne.(grand_ptot-1)) then 
   write(6,*) 'ktot = ',ktot,' grand_ptot-1 = ',grand_ptot-1
   stop 'ktot is not equal to grand_ptot-1 !!'
endif
!calling heapsort and sorting with increasing the distance from the origin(1:3)
call heapsort(ktot,rb,ra)
!maximum distance from origin(1:3). This is the radius at 
!which the HEALPix vectors should expand.
radius=ra(ktot)
!gives values for the first evaluation point which is the origin (1:3)
ep=0.
do j=0,nrays-1
   pdr(IDlist_dark(1))%epoint(1:3,j,0) = origin(1:3)
end do
!loops over all the domain and finds evaluation points. [straight N loop]
do k=1,ktot
   !locates the grid point in the new computational domain
   rvec(1)=pdr(rb(k))%x-origin(1)
   rvec(2)=pdr(rb(k))%y-origin(2)
   rvec(3)=pdr(rb(k))%z-origin(3)

   vel_vec(1) = pdr(rb(k))%vx
   vel_vec(2) = pdr(rb(k))%vy
   vel_vec(3) = pdr(rb(k))%vz

IF (fieldchoice.EQ."PNT") THEN
   temp_vec = rvec + vel_vec

   call rotate_z(temp_vec(1),temp_vec(2),temp_vec(3),psi_z)
   call rotate_x(temp_vec(1),temp_vec(2),temp_vec(3),psi_x)

   call rotate_z(rvec(1),rvec(2),rvec(3),psi_z)
   call rotate_x(rvec(1),rvec(2),rvec(3),psi_x)

   vel_vec = temp_vec - rvec   
ENDIF

   !next two lines return the ipix ray that the rvec(1:3) point belongs to.
   call vec2ang(rvec,theta,phi)
   call ang2pix_nest_id(nside,theta,phi,ipix)
   if (killray(ipix)) cycle
   healpixvector(1:3) = 1.1_DP*radius*vectors(1:3,ipix) !expand unit healpix vectors
   !calculates the angle along the line of sight of EVALUATION POINT -- HEALPIXVECTOR
   angle_los=acos(dot_product(rvec(1:3)-ep(1:3,ipix),healpixvector(1:3)-ep(1:3,ipix))/&
       &(sqrt((rvec(1)-ep(1,ipix))**2+(rvec(2)-ep(2,ipix))**2+(rvec(3)-ep(3,ipix))**2) * &
       &sqrt((healpixvector(1)-ep(1,ipix))**2+(healpixvector(2)-ep(2,ipix))**2+&
       &(healpixvector(3)-ep(3,ipix))**2)))
   !if the angle is less than the critical theta (user defined), then we have a new
   !evaluation point. This point is the projection of the grid point in the above line of sight.
   !if ((angle_los.le.theta_crit).and.(angle_los.ge.0D0)) then 
   if (angle_los.le.theta_crit) then
      !All next if-statements are conditions to avoid division by zero (i.e. x-plane, y-plane, z-plane)
      if (healpixvector(3).ne.0.0_dp) then
         ep(3,ipix) = (healpixvector(1)*healpixvector(3)*rvec(1) + healpixvector(2)*healpixvector(3)*&
                  &rvec(2) + (healpixvector(3)**2)*rvec(3))/(healpixvector(1)**2+healpixvector(2)**2+&
                  &healpixvector(3)**2)
         ep(2,ipix) = ep(3,ipix)*healpixvector(2)/healpixvector(3)
         ep(1,ipix) = ep(3,ipix)*healpixvector(1)/healpixvector(3)
      else
         if (healpixvector(1).eq.0.0_dp) then
            ep(1,ipix) = 0.0_dp
            ep(2,ipix) = rvec(2)
            ep(3,ipix) = 0.0_dp
         else if (healpixvector(2).eq.0.0_dp) then
            ep(1,ipix) = rvec(1)
            ep(2,ipix) = 0.0_dp
            ep(3,ipix) = 0.0_dp
         else
            ep(3,ipix) = 0.0_dp
            ep(1,ipix) = ((healpixvector(1)**2)*rvec(1) + healpixvector(1)*healpixvector(2)*rvec(2))/&
                  &(healpixvector(1)**2 + healpixvector(2)**2)
            ep(2,ipix) = ep(1,ipix) * healpixvector(2)/healpixvector(1)
         endif
      endif !healpixvector(3)
    !updates memory and stores the evaluation point in the original computational domain (so evaluation point + origin)
    pdr(IDlist_dark(1))%epray(ipix) = pdr(IDlist_dark(1))%epray(ipix)+1
    id = pdr(IDlist_dark(1))%epray(ipix)
    if (pdr(IDlist_dark(1))%epray(ipix).gt.maxpoints) STOP 'Increase maxpoints!'

    call project(vel_vec,healpixvector,pdr(IDlist_dark(1))%velocity(ipix,id) )

    pdr(IDlist_dark(1))%epoint(1:3,ipix,id)=ep(1:3,ipix)+origin(1:3)
    pdr(IDlist_dark(1))%projected(ipix,id)=rb(k)
    if (pdr(rb(k))%etype.eq.2) killray(ipix)=.true. !if the projected is ionized, stop propagating the ray (it has hit the HII region)
    endif !angle_los
enddo !k=1,ktot
deallocate(ra)
deallocate(rb)
deallocate(ep)
!==========================================================================================
endif
#ifdef OPENMP
!$OMP PARALLEL
!$OMP MASTER
CPUs = OMP_GET_NUM_THREADS()
write(6,*) "Proceeding for the PDR (PARALLEL)..."
write(6,*) "Number of CPUs: ", OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL
#else
write(6,*) 'Proceeding for the PDR (SERIAL)...'
#endif

!BUILDING HEALPIX VECTORS FOR ALL PDR ELEMENTS. 
!PARALLEL PROCESS-----------------------------
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(origin, p, ra, rb, ep) &
!$OMP PRIVATE(i, rvec, vel_vec, temp_vec, theta, phi, ipix, ktot, radius) &
!$OMP PRIVATE(j, healpixvector, angle_los, id, killray,psi_x,psi_z ) REDUCTION (+ : kk)
#endif 

do p=1,pdr_ptot
    !defines the origin to transfer all the domain in the original co-ordinate system
    killray=.false.
    origin(1:3) = pdrpoint(1:3,p)
    allocate(ra(0:grand_ptot-1)) !needs one extra place for sorting in heapsort
    allocate(rb(1:grand_ptot-1)) !-1 to avoid overlapping origin & pdrpoint
    allocate(ep(1:3,0:nrays-1))

    !calculating distances from the origin(1:3)
    kk=0
    pdr(IDlist_pdr(p))%epray = 0

IF (fieldchoice.EQ."PNT") THEN
    call calc_rotation_angles(UV_source_coord-origin,psi_z,psi_x)
ENDIF

    do i=1,grand_ptot
      if (i.eq.IDlist_pdr(p)) cycle
      kk=kk+1
      !locates the grid point in the new computational domain   
      rvec(1)=pdr(i)%x-origin(1)
      rvec(2)=pdr(i)%y-origin(2)
      rvec(3)=pdr(i)%z-origin(3) 

IF (fieldchoice.EQ."PNT") THEN
      call rotate_z(rvec(1),rvec(2),rvec(3),psi_z)
      call rotate_x(rvec(1),rvec(2),rvec(3),psi_x)
ENDIF
       
      !next two lines return the ipix ray that the rvec(1:3) point belongs to.
      call vec2ang(rvec,theta,phi)
      call ang2pix_nest_id(nside,theta,phi,ipix)
      ra(kk)=sqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)
      rb(kk)=i !stores the identifier of each grid point
    enddo
    ktot=kk !ktot should be grand_ptot-1
    if (ktot.ne.(grand_ptot-1)) then 
      write(6,*) 'ktot = ',ktot,' grand_ptot-1 = ',grand_ptot-1
      stop 'ktot is not equal to grand_ptot-1 !!'
    endif

    !calling heapsort and sorting with increasing the distance from the origin(1:3)
    call heapsort(ktot,rb,ra)
    !maximum distance from origin(1:3). This is the radius at 
    !which the HEALPix vectors should expand.
    radius=ra(ktot) 

    !gives values for the first evaluation point which is the origin (1:3)
    ep=0.
    do j=0,nrays-1
      pdr(IDlist_pdr(p))%epoint(1:3,j,0) = origin(1:3)
    end do
    !loops over all the domain and finds evaluation points. [straight N loop]
    do k=1,ktot
      !locates the grid point in the new computational domain
      rvec(1)=pdr(rb(k))%x-origin(1)
      rvec(2)=pdr(rb(k))%y-origin(2)
      rvec(3)=pdr(rb(k))%z-origin(3)

      vel_vec(1) = pdr(rb(k))%vx
      vel_vec(2) = pdr(rb(k))%vy
      vel_vec(3) = pdr(rb(k))%vz

!   test_module = sqrt(pdr(rb(k))%vx**2 + pdr(rb(k))%vy**2 + pdr(rb(k))%vz**2)
!   test_angle = acos(dot_product(rvec,vel_vec)/test_module/sqrt(rvec(1)**2 +rvec(2)**2 + rvec(3)**2))

IF (fieldchoice.EQ."PNT") THEN
      temp_vec = rvec + vel_vec
      
      call rotate_z(temp_vec(1),temp_vec(2),temp_vec(3),psi_z)
      call rotate_x(temp_vec(1),temp_vec(2),temp_vec(3),psi_x)

      call rotate_z(rvec(1),rvec(2),rvec(3),psi_z)
      call rotate_x(rvec(1),rvec(2),rvec(3),psi_x)

      vel_vec = temp_vec - rvec

ENDIF

      !next two lines return the ipix ray that the rvec(1:3) point belongs to.
      call vec2ang(rvec,theta,phi)
      call ang2pix_nest_id(nside,theta,phi,ipix)
      if (killray(ipix)) cycle
      healpixvector(1:3) = 1.1_DP*radius*vectors(1:3,ipix) !expand unit healpix vectors
      !calculates the angle along the line of sight of EVALUATION POINT -- HEALPIXVECTOR
      angle_los=acos(dot_product(rvec(1:3)-ep(1:3,ipix),healpixvector(1:3)-ep(1:3,ipix))/&
          &(sqrt((rvec(1)-ep(1,ipix))**2+(rvec(2)-ep(2,ipix))**2+(rvec(3)-ep(3,ipix))**2) * &
          &sqrt((healpixvector(1)-ep(1,ipix))**2+(healpixvector(2)-ep(2,ipix))**2+&
          &(healpixvector(3)-ep(3,ipix))**2)))
      !if the angle is less than the critical theta (user defined), then we have a new
      !evaluation point. This point is the projection of the grid point in the above line of sight.
      !if ((angle_los.le.theta_crit).and.(angle_los.ge.0D0)) then 
      if (angle_los.le.theta_crit) then
         !All next if-statements are conditions to avoid division by zero (i.e. x-plane, y-plane, z-plane)
         if (healpixvector(3).ne.0.0_dp) then
            ep(3,ipix) = (healpixvector(1)*healpixvector(3)*rvec(1) + healpixvector(2)*healpixvector(3)*&
                         &rvec(2) + (healpixvector(3)**2)*rvec(3))/(healpixvector(1)**2+healpixvector(2)**2+&
                         &healpixvector(3)**2)
            ep(2,ipix) = ep(3,ipix)*healpixvector(2)/healpixvector(3)
            ep(1,ipix) = ep(3,ipix)*healpixvector(1)/healpixvector(3)
         else
            if (healpixvector(1).eq.0.0_dp) then
               ep(1,ipix) = 0.0_dp
               ep(2,ipix) = rvec(2)
               ep(3,ipix) = 0.0_dp
            else if (healpixvector(2).eq.0.0_dp) then
               ep(1,ipix) = rvec(1)
               ep(2,ipix) = 0.0_dp
               ep(3,ipix) = 0.0_dp
            else
               ep(3,ipix) = 0.0_dp
               ep(1,ipix) = ((healpixvector(1)**2)*rvec(1) + healpixvector(1)*healpixvector(2)*rvec(2))/&
                            &(healpixvector(1)**2 + healpixvector(2)**2)
               ep(2,ipix) = ep(1,ipix) * healpixvector(2)/healpixvector(1)
            endif
         endif !healpixvector(3)
!       !updates memory and stores the evaluation point in the original computational domain (so evaluation point + origin)


       pdr(IDlist_pdr(p))%epray(ipix) = pdr(IDlist_pdr(p))%epray(ipix)+1
       id = pdr(IDlist_pdr(p))%epray(ipix)
       if (pdr(IDlist_pdr(p))%epray(ipix).gt.maxpoints) STOP 'Increase maxpoints!'

       call project(vel_vec,healpixvector,pdr(IDlist_pdr(p))%velocity(ipix,id-1) )       

       pdr(IDlist_pdr(p))%epoint(1:3,ipix,id)=ep(1:3,ipix)+origin(1:3)
       pdr(IDlist_pdr(p))%projected(ipix,id)=rb(k)
       if (pdr(rb(k))%etype.eq.2) killray(ipix)=.true.
       endif !angle_los
     enddo !k=1,ktot
    deallocate(ra)
    deallocate(rb)
    deallocate(ep)
enddo !pp/p=1,pdr_ptot
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif


suma=0
do pp=1,pdr_ptot
  p=IDlist_pdr(pp)
  suma = suma + sum(pdr(p)%epray(:))
enddo
if (dark_ptot.gt.0) then
!Include the Dark Molecular element
suma = suma + sum(pdr(IDlist_dark(1))%epray(:))
endif
write(6,*) 'No. evaluation points:',suma
write(6,*) 'Done!';write(6,*) ''

write(6,*) 'Checking for negative steps...'
adaptivemin=100.0D0
do pp=1,pdr_ptot
   p=IDlist_pdr(pp)
  do j=0,nrays-1
    if (pdr(p)%epray(j).gt.0) then
    do i=1,pdr(p)%epray(j)
       adaptive_step = sqrt((pdr(p)%epoint(1,j,0)-pdr(p)%epoint(1,j,i))**2+&
                           &(pdr(p)%epoint(2,j,0)-pdr(p)%epoint(2,j,i))**2+&
                           &(pdr(p)%epoint(3,j,0)-pdr(p)%epoint(3,j,i))**2)-&
                      &sqrt((pdr(p)%epoint(1,j,0)-pdr(p)%epoint(1,j,i-1))**2+&
                           &(pdr(p)%epoint(2,j,0)-pdr(p)%epoint(2,j,i-1))**2+&
                           &(pdr(p)%epoint(3,j,0)-pdr(p)%epoint(3,j,i-1))**2)
       if (adaptive_step.lt.0) stop 'found negative adaptive step!'
       if (adaptive_step.lt.adaptivemin) adaptivemin = adaptive_step
    enddo
    endif
   enddo
enddo
if (dark_ptot.gt.0) then
!Checking for the Dark Molecular element
p=IDlist_dark(1)
do j=0,nrays-1
  if (pdr(p)%epray(j).gt.0) then
     do i=1,pdr(p)%epray(j)
        adaptive_step = sqrt((pdr(p)%epoint(1,j,0)-pdr(p)%epoint(1,j,i))**2+&
                         &(pdr(p)%epoint(2,j,0)-pdr(p)%epoint(2,j,i))**2+&
                         &(pdr(p)%epoint(3,j,0)-pdr(p)%epoint(3,j,i))**2)-&
                    &sqrt((pdr(p)%epoint(1,j,0)-pdr(p)%epoint(1,j,i-1))**2+&
                         &(pdr(p)%epoint(2,j,0)-pdr(p)%epoint(2,j,i-1))**2+&
                         &(pdr(p)%epoint(3,j,0)-pdr(p)%epoint(3,j,i-1))**2)
        if (adaptive_step.lt.0) stop 'found negative adaptive step!'
        if (adaptive_step.lt.adaptivemin) adaptivemin = adaptive_step
     enddo
  endif
enddo
endif

write(6,*) 'No negative steps found'
write(6,*) 'Minimum adaptive step = ',adaptivemin


write(6,*) 'Assigning raytypes'
do pp=1,pdr_ptot
   p=IDlist_pdr(pp)
   allocate(pdr(p)%raytype(0:nrays-1))
   do j=0,nrays-1
     if (pdr(p)%epray(j).gt.0) then
       pdr(p)%raytype(j) = -pdr(pdr(p)%projected(j,pdr(p)%epray(j)))%etype
    else
       pdr(p)%raytype(j) = -pdr(p)%etype
     endif
  enddo
enddo

!Assigning raytype for the Dark Molecular element
if (dark_ptot.gt.0) then
  p=IDlist_dark(1)
  allocate(pdr(p)%raytype(0:nrays-1))
  do j=0,nrays-1
    if (pdr(p)%epray(j).gt.0) then
       pdr(p)%raytype(j) = -pdr(pdr(p)%projected(j,pdr(p)%epray(j)))%etype
    else
       pdr(p)%raytype(j) = -pdr(p)%etype
    endif
  enddo
endif

return
end subroutine evaluation_points

subroutine calc_rotation_angles(UV_vector,psi_z,psi_x)
use definitions
use maincode_module
real(kind=dp), intent(in) :: UV_vector(1:3)
real(kind=dp), intent(out) :: psi_z, psi_x

if ((UV_vector(2).eq.0).and.(UV_vector(1).gt.0)) then
    psi_z = pi/2.0D0
elseif ((UV_vector(2).eq.0).and.(UV_vector(1).lt.0)) then
    psi_z = -pi/2.0D0
elseif ((UV_vector(2).eq.0).and.(UV_vector(1).eq.0)) then
    psi_z = 0.0D0
elseif (UV_vector(2).gt.0) then
     psi_z = atan(UV_vector(1)/UV_vector(2))
else
     psi_z = -pi+atan(UV_vector(1)/UV_vector(2))
endif

if ((UV_vector(2).eq.0).and.(UV_vector(1).eq.0).and.(UV_vector(3).gt.0)) then
   psi_x = -pi/2.0D0
elseif ((UV_vector(2).eq.0).and.(UV_vector(1).eq.0).and.(UV_vector(3).lt.0)) then
   psi_x = pi/2.0D0
elseif ((UV_vector(2).eq.0).and.(UV_vector(1).eq.0).and.(UV_vector(3).eq.0)) then
   psi_x = 0.0D0
else
   psi_x = atan(-UV_vector(3)/(UV_vector(1)*sin(psi_z) + UV_vector(2)*cos(psi_z) )  )
endif
end subroutine calc_rotation_angles

subroutine rotate_z(x,y,z,psi_z)
use definitions
use maincode_module
real(kind=dp), intent(inout) :: x, y, z, psi_z
real(kind=dp) :: x_prev, y_prev, z_prev
   x_prev = x
   y_prev = y
   z_prev = z
   x = cos(psi_z)*x_prev - sin(psi_z)*y_prev
   y = sin(psi_z)*x_prev + cos(psi_z)*y_prev
   z = z_prev
end subroutine rotate_z

subroutine rotate_x(x,y,z,psi_x)
use definitions
use maincode_module
real(kind=dp), intent(inout) :: x, y, z, psi_x
real(kind=dp) :: x_prev, y_prev, z_prev
   x_prev = x
   y_prev = y
   z_prev = z
   x = x_prev
   y = cos(psi_x)*y_prev - sin(psi_x)*z_prev
   z = sin(psi_x)*y_prev + cos(psi_x)*z_prev
end subroutine rotate_x

subroutine project(vec,hvector,proj_length)
use definitions
real(kind=dp), intent(in) :: vec(1:3), hvector(1:3)
real(kind=dp), intent(out) :: proj_length
real(kind=dp) :: cos_angle, projection(1:3)
      if (hvector(3).ne.0.0_dp) then
         projection(3) = (hvector(1)*hvector(3)*vec(1) + hvector(2)*hvector(3)*&
                  &vec(2) + (hvector(3)**2)*vec(3))/(hvector(1)**2+hvector(2)**2+&
                  &hvector(3)**2)
         projection(2) = projection(3)*hvector(2)/hvector(3)
         projection(1) = projection(3)*hvector(1)/hvector(3)
      else
         if (hvector(1).eq.0.0_dp) then
            projection(1) = 0.0_dp
            projection(2) = vec(2)
            projection(3) = 0.0_dp
         else if (hvector(2).eq.0.0_dp) then
            projection(1) = vec(1)
            projection(2) = 0.0_dp
            projection(3) = 0.0_dp
         else
            projection(3) = 0.0_dp
            projection(1) = ((hvector(1)**2)*vec(1) + hvector(1)*hvector(2)*vec(2))/&
                  &(hvector(1)**2 + hvector(2)**2)
            projection(2) = projection(1) * hvector(2)/hvector(1)
         endif
      endif
     proj_length = sqrt(projection(1)**2 + projection(2)**2 + projection(3)**2)
     if (proj_length > 1.0D-8) then
         cos_angle = dot_product(projection,hvector)/proj_length/&
                 &sqrt(hvector(1)**2 + hvector(2)**2 + hvector(3)**2)
         if (abs(cos_angle + 1.0D0) < 1.0D-4) then !if cos = -1.0
             proj_length = -proj_length
         endif
     endif

end subroutine




