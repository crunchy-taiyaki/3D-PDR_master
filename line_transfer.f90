module line_transfer
    use constants
    implicit none

type line
 character(len=20) :: spec
 real(8) :: freq0
 real(8) :: g_j
 real(8) :: g_i
 real(8) :: A !Einstein A coefficient
 integer :: proton_number
contains
 procedure :: solve_rt
end type

 contains

 subroutine solve_rt(this,directory,x,av,v_gas,abun,Tgas,Tdust,rho,pop_i,pop_j,metallicity,gas_to_dust,&
                        &nfreq,min_gas_velocity, max_gas_velocity)
 class (line), intent(in) :: this
 character(len=80), intent(in) :: directory
 real(8), intent(in) :: x(1:),av(1:),v_gas(1:),abun(1:),Tgas(1:),Tdust(1:),rho(1:),pop_i(1:),pop_j(1:)
 real(8), intent(in) :: metallicity, gas_to_dust
 integer, intent(in) :: nfreq
 real(8), intent(in) :: min_gas_velocity, max_gas_velocity
 integer :: ptot !total point number
 real(8), parameter :: rho_grain = 2.0D0
 real(8), parameter :: v_turb = 1.0d5 !cm/s
 real(8) :: tmp
 real(8), allocatable :: ngrain(:), emissivity(:)
 real(8), allocatable :: BB(:),BB_dust(:)
 real(8), allocatable :: tpop(:), S(:)
 integer :: p, ifreq !counters
 real(8) :: mh
 real(8) :: sigma, sigma_ptot
 real(8), allocatable :: tau(:,:)
 real(8), allocatable :: freq(:), phi(:), velocities(:), Tex(:)
 real(8) :: frac
 real(8), allocatable :: dtau(:),tau_incr(:)
 real(8), allocatable :: current_intensity(:),intensity(:,:),antenna_temp(:,:),N(:)
 character(len=100) :: fileout, fileout_radial
 ptot = size(x)
 allocate(ngrain(1:ptot), emissivity(1:ptot))
 allocate(BB(1:ptot), BB_dust(1:ptot))
 allocate(tpop(1:ptot), S(1:ptot))
 allocate(tau(1:ptot,0:nfreq-1))
 allocate(freq(0:nfreq-1), phi(0:nfreq-1), velocities(0:nfreq-1))
 allocate(Tex(1:ptot))
 allocate(tau_incr(0:nfreq-1),dtau(0:nfreq-1))
 allocate(current_intensity(0:nfreq-1),intensity(1:ptot,0:nfreq-1),antenna_temp(1:ptot,0:nfreq-1), N(1:ptot))

!optical depth calculation
mh = this%proton_number*mhp
sigma_ptot=(this%freq0/c)*sqrt(kb*Tgas(ptot-1)/mh+v_turb**2/2.)
  do ifreq=0,nfreq-1
    velocities(ifreq) = min_gas_velocity + ifreq*(max_gas_velocity-min_gas_velocity)/(nfreq-1)
  enddo
  freq = this%freq0/(1.+velocities*1e5/c)
  tau_incr(:) = 0.
  do p=1,ptot-1
  sigma=(this%freq0/c)*sqrt(kb*Tgas(p)/mh+v_turb**2/2.)
  phi=exp(-((1+v_gas(p)*1.0d5/c)*freq(:)-this%freq0)**2/sigma**2/2.)/sigma/sqrt(2.*pi)
  !phi = 1./sigma/sqrt(2.*pi)
  frac=0.5*((pop_j(p)+pop_j(p+1))*this%g_i/this%g_j-(pop_i(p)+pop_i(p+1)))
  tau_incr(:) = tau_incr(:)+phi(:)*(this%A*c**2/8./pi/this%freq0**2)*frac*abs(x(p+1)-x(p))*pc
  !where (tau_incr(:).le.-5) ! Prevent exploding beta values caused by strong masing (tau < -5)
  !  tau_incr = -5.
  !end where
  tau(p,:) = tau_incr(:)
  enddo
  

!velocities(:) = c*(freq(:)/this%freq0 - 1.0D0)*1d-5
Tex(:)=(hp*this%freq0/kb)/log(this%g_i*pop_j(:)/pop_i(:)/this%g_j)

!source function calculation      
  tmp=2.0D0*hp*(this%freq0**3)/(c**2)
  BB(:) = tmp*(1.0D0/(exp(hp*this%freq0/kb/2.7D0)-1.0D0))!Planck function !2.7D0 is the CMBR temperature
  ngrain=2.0D-12*rho(:)*metallicity*100./gas_to_dust
  emissivity=(rho_grain*ngrain(:))*(0.01*(1.3*this%freq0/3.0D11))
  BB_dust(:) = tmp*(1.0D0/(exp(hp*this%freq0/kb/Tdust(:))-1.D0)*emissivity(:))
  BB = BB + BB_dust
  where (pop_i(:).eq.0)
    S(:)=0.0D0
  elsewhere
    tpop=(pop_j(:)*this%g_i)/(pop_i(:)*this%g_j)-1.0D0
    where (abs(tpop).lt.1.0D-50)
      S(:)=hp*this%freq0*pop_i(:)*this%A/4./pi
    elsewhere
      S(:)=tmp/tpop
    end where
  end where

!radiation transfer solving
 current_intensity(:) = 0.
do p=1,ptot-2
   dtau = tau(p+1,:)-tau(p,:)
    where (dtau(:).gt.1d10)
      current_intensity(:)=S(p)
    elsewhere (abs(dtau(:)).gt.1d-6)
      current_intensity(:)=current_intensity(:)*exp(-dtau)+&
             &S(p)*((1-exp(-dtau))/dtau-exp(-dtau))+&
             &S(p+1)*(1.-(1.-exp(-dtau))/dtau)
    elsewhere (abs(dtau(:)).le.1d-6)
      current_intensity(:)=current_intensity(:)*(1-dtau)+(S(p)+S(p+1))*dtau/2.
    end where
      intensity(p,:)=current_intensity(:)
    if (current_intensity(nfreq/2).gt.intensity(p-1,nfreq/2)) then
      intensity(p,:)=current_intensity(:)
    else
      intensity(p,:)=intensity(p-1,:)
    endif

enddo
antenna_temp = intensity*c**2/2./kb/this%freq0**2
N(1) = 0.5*(rho(1)*abun(1)+rho(2)*abun(2))*abs(x(2)-x(1))*pc
do p=2,ptot-2
  N(p)=N(p-1)+0.5*(rho(p)*abun(p)+rho(p+1)*abun(p+1))*abs(x(p+1)-x(p))*pc
enddo

fileout = trim(adjustl(directory))//'/'//trim(adjustl(this%spec))//'.dat'
open(unit=100,file=fileout,status='replace')
  write(100,*) velocities(:)
  write(100,*) tau(ptot-2,:)
  write(100,*) antenna_temp(ptot-2,:)
 close(100)

fileout_radial = trim(adjustl(directory))//'/'//trim(adjustl(this%spec))//'.radial.dat'
open(unit=101,file=fileout_radial,status='replace')
do p=1,ptot-2
  write(101,*) p, av(p), abun(p), N(p), Tgas(p), Tdust(p)
enddo
 close(101)
write(6,'(A10,2X,5ES11.3)') this%spec,N(ptot-2),Tex(ptot-2),tau(ptot-2,nfreq/2),antenna_temp(ptot-2,nfreq/2)

 deallocate(ngrain, emissivity)
 deallocate(BB, BB_dust)
 deallocate(tpop, S)
 deallocate(tau)
 deallocate(freq, phi, velocities)
 deallocate(Tex)
 deallocate(tau_incr,dtau)
 deallocate(current_intensity,intensity,antenna_temp, N)
end subroutine solve_rt	

end module line_transfer

