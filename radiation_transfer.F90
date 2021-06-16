subroutine radiation_transfer(pdr_ptot,nlev,nfreq,&
                     &frequencies,metallicity,density,dust_temperature,&
                     &s_pop_array,x_array,weights,A_COEFFS,&
                     &tau_profile_array,bright_temperature,velocities,coolant,Tguess,v_turb)

use definitions
use maincode_module, only : p,pdr,vectors
use healpix_types
use healpix_module
use global_module, only: g2d
implicit none

integer(kind=i4b), intent(in) :: pdr_ptot
integer(kind=i4b), intent(in) :: nlev
integer(kind=i4b), intent(in) :: nfreq
real(kind=dp), intent(in) :: frequencies(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: metallicity
real(kind=dp), intent(in) :: density(1:pdr_ptot),dust_temperature(1:pdr_ptot)
real(kind=dp), intent(in) :: s_pop_array(1:nlev,1:pdr_ptot)
real(kind=dp), intent(in) :: x_array(1:pdr_ptot)
real(kind=dp), intent(in) :: weights(1:nlev)
real(kind=dp), intent(in) :: A_COEFFS(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: tau_profile_array(1:nlev,1:nlev,0:nfreq-1,1:pdr_ptot)
integer(kind=i4b), intent(in) :: coolant
real(kind=dp), intent(in) :: Tguess, v_turb
integer(kind=i4b) :: ilevel, jlevel
integer(kind=i4b) :: pp
integer(kind=i4b) :: ifreq
real(kind=dp) :: thermal_velocity, TMP
real(kind=dp) :: v_gas = 0.0D0 !in sm*s^-1
real(kind=dp) :: frequency(0:nfreq-1)
real(kind=dp) :: doppler_profile(0:nfreq-1)
real(kind=dp) :: BB_ij(1:pdr_ptot), TPOP(1:pdr_ptot), S_ij(1:pdr_ptot)
real(kind=dp) :: rho_grain
real(kind=dp) :: ngrain(1:pdr_ptot), emissivity(1:pdr_ptot), BB_ij_dust(1:pdr_ptot)
real(kind=dp) :: frac1, frac2, rhs2
real(kind=dp) :: tau(0:nfreq-1,1:pdr_ptot)
real(kind=dp) :: dtau(0:nfreq-1),beta(0:nfreq-1)
real(kind=dp) :: intensity_profile(1:nlev,1:nlev,0:nfreq-1,1:pdr_ptot)
real(kind=dp) :: velocity_limit=0.1D0
real(kind=dp), intent(out) :: bright_temperature(1:nlev,1:nlev,0:nfreq-1)
real(kind=dp), intent(out) :: velocities(1:nlev,1:nlev,0:nfreq-1)

thermal_velocity=sqrt(8*KB*Tguess/PI/MP+v_turb**2)
    do ilevel=1,nlev
       do jlevel=1,nlev !i>j
         if (jlevel.ge.ilevel) exit

         !frequancies and velocities calcuation
         do ifreq=0,nfreq-1
           frequency(ifreq)=frequencies(ilevel,jlevel)-3*thermal_velocity+ifreq*2*3*thermal_velocity/(nfreq-1)
         enddo
         velocities(ilevel,jlevel,:) = C*(frequency(:)/frequencies(ilevel,jlevel) - 1.0D0)*1d-5

         !source function calculation
         TMP=2.0D0*hp*(frequencies(ilevel,jlevel)**3)/(c**2)
         BB_ij = TMP*(1.0D0/(EXP(hp*frequencies(ilevel,jlevel)/kb/2.7D0)-1.0D0))!Planck function !2.7D0 is the CMBR temperature
         ngrain(1:pdr_ptot)=2.0D-12*density(:)*metallicity*100./g2d
         rho_grain=2.0D0
         emissivity(1:pdr_ptot)=(rho_grain*ngrain(:))*(0.01*(1.3*frequencies(ilevel,jlevel)/3.0D11))
         BB_ij_dust(1:pdr_ptot) = TMP*(1.0D0/(EXP(hp*frequencies(ilevel,jlevel)/kb/dust_temperature(:))-1.D0)*emissivity(:))
         BB_ij(1:pdr_ptot) = BB_ij + BB_ij_dust
         where (s_pop_array(ilevel,:).eq.0)
           S_ij(1:pdr_ptot)=0.0D0
         elsewhere
           TPOP(1:pdr_ptot)=(s_pop_array(jlevel,:)*weights(ilevel))/(s_pop_array(ilevel,:)*weights(jlevel))-1.0D0
           where(abs(TPOP).lt.1.0D-50)
             S_ij(1:pdr_ptot)=hp*frequencies(ilevel,jlevel)*s_pop_array(ilevel,:)*A_COEFFS(ilevel,jlevel)/4./pi
           elsewhere
	     S_ij(1:pdr_ptot)=TMP/TPOP
           end where
         end where

         !optical depth calculation
         tau(:,1)=0.0D0
	 frac1=(A_COEFFS(ilevel,jlevel)*(C**3))/(8.0*pi*(frequencies(ilevel,jlevel)**3))
	 doppler_profile(0:nfreq-1)=exp(-((1+v_gas/C)*frequency(:)-frequencies(ilevel,jlevel))**2/thermal_velocity**2)/thermal_velocity
	 do pp=1,pdr_ptot-1
             frac2 = 0.5*((s_pop_array(jlevel,pp)+s_pop_array(jlevel,pp+1))*weights(ilevel)/weights(jlevel)-&
                     &(s_pop_array(ilevel,pp)+s_pop_array(ilevel,pp+1)))
             tau(:,pp+1)=tau(:,pp)+frac1*doppler_profile(:)*frac2*abs(x_array(pp+1)-x_array(pp))*PC
         enddo

         !radiation transfer solving
         intensity_profile(ilevel,jlevel,:,1) = 0.0D0
	 do pp=1,pdr_ptot-2    
          dtau(0:nfreq-1) = tau(:,pp)
          where((dtau(:).ge.1d-6).and.(dtau(:).le.1d10))
           beta(:) = (1-EXP(-dtau))/dtau
           intensity_profile(ilevel,jlevel,:,pp+1) = intensity_profile(ilevel,jlevel,:,pp)*EXP(-dtau)+&
                                  &S_ij(pp+1)*(1-beta(:))+S_ij(pp)*(beta(:)-EXP(-dtau))
	   elsewhere(dtau(:).gt.1d10)
             intensity_profile(ilevel,jlevel,:,pp+1) = BB_ij(pp+1)
           elsewhere(dtau(:).lt.1d-6)
             intensity_profile(ilevel,jlevel,:,pp+1) = (1-dtau)*intensity_profile(ilevel,jlevel,:,pp)+&
                                                      &dtau*(BB_ij(pp)+BB_ij(pp+1))/2
           end where
         enddo
       
       !convert intensity to brightness temperature                          
       bright_temperature(ilevel,jlevel,:) = intensity_profile(ilevel,jlevel,:,pdr_ptot-1)*c**2/(2*kb*frequencies(ilevel,jlevel)**2)
       enddo !jlevel=1,nlev
     enddo !ilevel=1,nlev

  return

end subroutine radiation_transfer
