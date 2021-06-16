!T.Bisbas, T.Bell

subroutine escape_probability(transition, dust_temperature, nrays, nlev,nfreq, &
                   &A_COEFFS, B_COEFFS, C_COEFFS, &
                   &frequencies,s_evalpop, maxpoints,&
                   & T_evalpoint, vel_evalpoint, Tguess,velocity_flag,mode,v_turb, v_gas, &
                   &s_jjr, s_pop, s_evalpoint, weights,cooling_rate,line,tau,&
                   &coolant,density,metallicity,bbeta)


use definitions
use maincode_module, only : p,pdr,vectors
use healpix_types
use healpix_module
use global_module, only: g2d

implicit none

integer(kind=i4b), intent(in) :: nrays
integer(kind=i4b), intent(in) :: nlev
integer(kind=i4b), intent(in) :: nfreq
integer(kind=i4b), intent(in) :: maxpoints
integer(kind=i4b), intent(in) :: s_jjr(0:nrays-1)
integer(kind=i4b), intent(in) :: coolant
real(kind=dp), intent(in) :: A_COEFFS(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: B_COEFFS(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: C_COEFFS(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: frequencies(1:nlev, 1:nlev)
real(kind=dp), intent(in) :: s_evalpop(0:nrays-1,0:maxpoints,1:nlev)
real(kind=dp), intent(in) :: s_evalpoint(1:3,0:nrays-1,0:maxpoints)
real(kind=dp), intent(in) :: T_evalpoint(0:nrays-1,0:maxpoints)
character(len=1),intent(in) :: velocity_flag
character(len=4),intent(in) :: mode
real(kind=dp), intent(in) :: vel_evalpoint(0:nrays-1,0:maxpoints)
real(kind=dp) :: Tguess, v_turb, v_gas !, intent(in)
real(kind=dp), intent(in) :: weights(1:nlev)
real(kind=dp), intent(in) :: s_pop(1:nlev)
real(kind=dp), intent(in) :: dust_temperature,density,metallicity

real(kind=dp), intent(out) :: line(1:nlev,1:nlev)
real(kind=dp), intent(out) :: cooling_rate
real(kind=dp), intent(out) :: tau(1:nlev,1:nlev,0:nrays-1)
real(kind=dp),intent(out) :: bbeta(1:nlev,1:nlev,0:nrays-1)
real(kind=dp), intent(inout) :: transition(1:nlev,1:nlev)

integer(kind=i4b) :: i, j
integer(kind=i4b) :: ifreq
integer(kind=i4b) :: ilevel, jlevel
real(kind=dp) :: beta_ij, beta_ij_sum
real(kind=dp) :: frac1, frac2, frac3, rhs2 
real(kind=dp) :: sigma,sigma_p, thermal_vel
real(kind=dp) :: phi
real(kind=dp) :: tpop, tmp2
real(kind=dp) :: S_ij, BB_ij, J_ij
real(kind=dp) :: tau_increment
real(kind=dp) :: freq(0:nfreq-1),phi_p(0:nfreq-1)
real(kind=dp) :: delta_vel
real(kind=dp) :: BB_CMBR, BB_ij_gas
real(kind=dp) :: tau_ij(0:nrays-1)
real(kind=dp) :: tau_nu(0:nfreq-1),beta_nu(0:nfreq-1),I_nu(0:nfreq-1),J_ij_ray(0:nrays-1)
real(kind=dp) :: beta_ij_ray(0:nrays-1)
real(kind=dp) :: field(1:nlev,1:nlev)
real(kind=dp) :: emissivity, bb_ij_dust, ngrain, rho_grain
character(len=1):: dummystop


line=0.0D0
cooling_rate = 0.0D0
field=0.0D0
frac2=1.0D0/sqrt(8.0*KB*Tguess/PI/MP + v_turb**2)


    do ilevel=1,nlev
       do jlevel=1,nlev !i>j
         if (jlevel.ge.ilevel) exit

	 if (mode.eq.'full') then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!START OF MODE:FULL	 
	 !init frequency array
	 sigma_p=sqrt(8.0*KB*Tguess/PI/MP + v_turb**2)*frequencies(ilevel,jlevel)/C
	 do ifreq=0,nfreq-1
	   freq(ifreq)=frequencies(ilevel,jlevel)*C/(C+v_gas*1d5)-3*sigma_p+ifreq*3*sigma_p*2/(nfreq-1)
	 enddo !ifreq=0,nfreq-1	

         tau_nu = 0.
         tau_ij = 0.
         beta_ij_ray = 0.
         beta_ij = 0.
         J_ij = 0.

         do j=0,nrays-1

	 do ifreq=0,nfreq-1
	    !CMBR emission
            BB_CMBR = 2.0D0*HP*(freq(ifreq)**3)*(1.0D0/(EXP(HP*freq(ifreq)/KB/2.7D0)-1.0D0))&
                            &/(C**2) !Planck function !2.7D0 is the CMBR temperature
            I_nu(ifreq) = BB_CMBR
	    do i=1,s_jjr(j)

#ifdef PSEUDO_1D
             if (j.ne.6) then
               tau_nu(j) = 1.0D50
             else
#endif
#ifdef PSEUDO_2D
             if (abs(vectors(3,j).gt.1d-10) then
	       tau_nu(j) = 1.0D50 !Not in Equator
             else
#endif
            !calculations of tau_ij
            frac1=(A_COEFFS(ilevel,jlevel)*(C**3))/(8.0*pi*(frequencies(ilevel,jlevel)**3))
            frac3=((s_evalpop(j,i-1,jlevel)*weights(ilevel)-s_evalpop(j,i-1,ilevel)*weights(jlevel))+&
            &(s_evalpop(j,i,jlevel)*weights(ilevel)-s_evalpop(j,i,ilevel)*weights(jlevel)))/2./weights(jlevel)
            rhs2=sqrt((s_evalpoint(1,j,i-1)-s_evalpoint(1,j,i))**2+&
              &(s_evalpoint(2,j,i-1)-s_evalpoint(2,j,i))**2+&
              &(s_evalpoint(3,j,i-1)-s_evalpoint(3,j,i))**2) !adaptive step
            sigma = sqrt(8.0*KB*T_evalpoint(j,i)/PI/MP + v_turb**2)*frequencies(ilevel,jlevel)/C
            phi = exp(-((1+vel_evalpoint(j,i)*1d5/C)*freq(ifreq)-frequencies(ilevel,jlevel))**2/&
                        &(2.0*sigma**2))/&
                            &(sigma*sqrt(2.*pi))
            tau_increment=frac1*phi*frac3*rhs2*PC
            tau_nu(ifreq)=tau_nu(ifreq)+tau_increment !optical depth
#ifdef PSEUDO_1D
            endif !j.ne.6
#endif
#ifdef PSEUDO_2D
            endif !(abs(vectors(3,j).gt.1d-10)
#endif
	    !Dust emission
            NGRAIN=2.0D-12*density*metallicity*100./g2d
            rho_grain=2.0D0
            EMISSIVITY=(RHO_GRAIN*NGRAIN)*(0.01*(1.3*freq(ifreq)/3.0D11))
            BB_ij_dust = (1.0D0/(EXP(HP*freq(ifreq)/KB/DUST_TEMPERATURE)-1.D0)*EMISSIVITY)*&
                               &2.0D0*HP*(freq(ifreq)**3)/(C**2)

            !Blackbody emission
            BB_ij_gas = (1.0D0/(EXP(HP*frequencies(ilevel,jlevel)/KB/Tguess)-1.0D0))*&
                            &2.0D0*HP*(frequencies(ilevel,jlevel)**3)/(C**2)

            !calculation of source function (taken from UCL_PDR)
            TMP2=2.0D0*HP*(FREQUENCIES(ilevel,jlevel)**3)/(C**2)
            if (s_evalpop(j,i,ilevel).eq.0) then
              S_ij=0.0D0
            else	    
              TPOP=(s_evalpop(j,i,jlevel)*WEIGHTS(ilevel))/(s_evalpop(j,i,ilevel)*WEIGHTS(jlevel))-1.0D0
              if(abs(TPOP).lt.1.0D-50) then
                S_ij=HP*FREQUENCIES(ilevel,jlevel)*s_evalpop(j,i,ilevel)*A_COEFFS(ilevel,jlevel)/4./pi
              else
                S_ij=TMP2/TPOP
              endif
            endif

            !radiation transfer in integral form
            if (tau_increment.gt.1d10) then
              I_nu(ifreq) = S_ij
            elseif (tau_increment.gt.1d-6) then
	      I_nu(ifreq) = exp(-tau_increment)*I_nu(ifreq) + (1 + exp(-tau_increment))*(S_ij + BB_ij_dust)
            else
              I_nu(ifreq) = I_nu(ifreq)
            endif

         enddo!i=1,s_jjr(j) 

            !Doppler profile for element p
	    phi_p(ifreq) = exp(-((1+v_gas*1d5/C)*freq(ifreq)-frequencies(ilevel,jlevel))**2/&
                        &(2.0*sigma_p**2))/&
                            &(sigma_p*sqrt(2.*pi))

            !monochromatic escape probability
            if (tau_nu(ifreq).lt.-5.0D0) then
              beta_nu(ifreq)=(1.0D0-EXP(5.0D0))/(-5.0D0)*phi_p(ifreq)
            else if (abs(tau_nu(ifreq)).lt.1.0D-8) then
              beta_nu(ifreq)=1.0D0*phi_p(ifreq)
            else
              beta_nu(ifreq)=(1.0D0-EXP(-tau_nu(ifreq)))/tau_nu(ifreq)*phi_p(ifreq)
             endif
             
             enddo !ifreq=0,nfreq-1

             ! mean integrated intensity calculation
             do ifreq=0,nfreq-2
               J_ij_ray(j) = J_ij_ray(j)+&
                 &(I_nu(ifreq+1)*phi_p(ifreq+1)+I_nu(ifreq)*phi_p(ifreq))*&
                 &abs(freq(ifreq+1)-freq(ifreq))/2.
             enddo !ifreq=0,nfreq-2

           ! tau and beta integration
             if (j.ne.6) then
               tau_ij(j) = 1.0D50
               beta_ij_ray(j) = 0.
             else

              do ifreq=0,nfreq-2
                tau_ij(j) = tau_ij(j) +&
                 &(tau_nu(ifreq+1)+tau_nu(ifreq))*&
                 &abs(freq(ifreq+1)-freq(ifreq))/2.
                beta_ij_ray(j) = beta_ij_ray(j)+&
                 &(beta_nu(ifreq+1)+beta_nu(ifreq))*&
                 &abs(freq(ifreq+1)-freq(ifreq))/2.       
              enddo !ifreq=0,nfreq-2
             endif !j.ne.6

         tau(ilevel,jlevel,j) = tau_ij(j)
         bbeta(ilevel,jlevel,j) = beta_ij_ray(j)
         enddo !j=0,nrays-1

#ifdef PSEUDO_1D
         beta_ij = sum(beta_ij_ray)
         field(ilevel,jlevel) = sum(J_ij_ray)
#elif PSEUDO_2D
         beta_ij = sum(beta_ij_ray) / 4.
         field(ilevel,jlevel) = sum(J_ij_ray) / 4.
#else
         beta_ij = sum(beta_ij_ray) / real(nrays,kind=DP) 
         field(ilevel,jlevel) = sum(J_ij_ray) / real(nrays,kind=DP) 
#endif

         bbeta(ilevel,jlevel,6) = beta_ij
         !line(ilevel,jlevel) = A_COEFFS(ilevel,jlevel)*HP*frequencies(ilevel,jlevel) * &
         !                    & s_pop(ilevel)*beta_ij*(S_ij-BB_ij)/S_ij
         line(ilevel,jlevel) = field(ilevel,jlevel)/BB_ij_gas
         cooling_rate = cooling_rate + line(ilevel,jlevel)
         goto 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 else !mode.eq.'lvg'
         tau_ij=0.0D0
         beta_ij=0.0D0; beta_ij_ray=0.0D0
	 beta_ij_sum=0.0D0

         frac1=(A_COEFFS(ilevel,jlevel)*(C**3))/(8.0*pi*(frequencies(ilevel,jlevel)**3))
         TMP2=2.0D0*HP*(FREQUENCIES(ilevel,jlevel)**3)/(C**2)
         BB_ij = TMP2*(1.0D0/(EXP(HP*frequencies(ilevel,jlevel)/KB/2.7D0)-1.0D0)) !Planck function !2.7D0 is the CMBR temperature
         NGRAIN=2.0D-12*density*metallicity*100./g2d
         rho_grain=2.0D0
         EMISSIVITY=(RHO_GRAIN*NGRAIN)*(0.01*(1.3*FREQUENCIES(ilevel,jlevel)/3.0D11))
         BB_ij_dust = TMP2*(1.0D0/(EXP(HP*frequencies(ilevel,jlevel)/KB/DUST_TEMPERATURE)-1.D0)*EMISSIVITY)
         BB_ij = BB_ij + BB_ij_dust
         if (s_pop(ilevel).eq.0) then
            S_ij=0.0D0
            beta_ij=1.0D0
            goto 2
         endif
         TPOP=(s_pop(jlevel)*WEIGHTS(ilevel))/(s_pop(ilevel)*WEIGHTS(jlevel))-1.0D0
         if(abs(TPOP).lt.1.0D-50) then
              S_ij=HP*FREQUENCIES(ilevel,jlevel)*s_pop(ilevel)*A_COEFFS(ilevel,jlevel)/4./pi
              beta_ij=1.0D0
              goto 1
         else
         !calculation of source function (taken from UCL_PDR)
              S_ij=TMP2/TPOP
         endif
         do j=0,nrays-1
#ifdef PSEUDO_1D
         if (j.ne.6) then
           tau_ij(j) = 1.0D50
         else
#endif
#ifdef PSEUDO_2D
         if (abs(vectors(3,j).gt.1d-10) then
	     tau_ij(j) = 1.0D50 !Not in Equator
#endif

         if (velocity_flag.eq.'n') then
           do i=1,s_jjr(j)
             !calculations of tau_ij
           frac3=((s_evalpop(j,i-1,jlevel)*weights(ilevel)-s_evalpop(j,i-1,ilevel)*weights(jlevel))+&
           &(s_evalpop(j,i,jlevel)*weights(ilevel)-s_evalpop(j,i,ilevel)*weights(jlevel)))/2./weights(jlevel)
           rhs2=sqrt((s_evalpoint(1,j,i-1)-s_evalpoint(1,j,i))**2+&
              &(s_evalpoint(2,j,i-1)-s_evalpoint(2,j,i))**2+&
              &(s_evalpoint(3,j,i-1)-s_evalpoint(3,j,i))**2) !adaptive step
           tau_increment=frac1*frac2*frac3*rhs2*PC
           tau_ij(j)=tau_ij(j)+tau_increment !optical depth
           enddo !i=1,jr(j)

         else
           i=1
           delta_vel = abs(vel_evalpoint(j,i)-v_gas)*1d5
           thermal_vel=sqrt(8.0*KB*Tguess/PI/MP + v_turb**2)

           do i=1,s_jjr(j)
             !calculations of tau_ij
           frac3=((s_evalpop(j,i-1,jlevel)*weights(ilevel)-s_evalpop(j,i-1,ilevel)*weights(jlevel))+&
           &(s_evalpop(j,i,jlevel)*weights(ilevel)-s_evalpop(j,i,ilevel)*weights(jlevel)))/2./weights(jlevel)
           rhs2=sqrt((s_evalpoint(1,j,i-1)-s_evalpoint(1,j,i))**2+&
              &(s_evalpoint(2,j,i-1)-s_evalpoint(2,j,i))**2+&
              &(s_evalpoint(3,j,i-1)-s_evalpoint(3,j,i))**2) !adaptive step
           if (delta_vel.lt.thermal_vel) then
	     tau_increment = frac1*frac2*frac3*rhs2*PC
           else
	     tau_increment = frac1*frac3*rhs2*PC/delta_vel
           endif
           tau_ij(j)=tau_ij(j)+tau_increment !optical depth
           enddo !i=1,jr(j)
         endif	

#ifdef PSEUDO_1D
         endif
#endif
#ifdef PSEUDO_2D
         endif
#endif
          ! Prevent exploding beta values caused by strong masing (tau < -5)
          ! Assume tau = -5 and calculate the escape probability accordingly
          if (tau_ij(j).lt.-5.0D0) then
             beta_ij_ray(j)=(1.0D0-EXP(5.0D0))/(-5.0D0)
          ! Treat weak masing using the standard escape probability formalism
          ! else if (tau_ij(j).lt.0.0D0) then
          !     beta_ij_ray(j)=(1.0D0-EXP(-tau_ij(j)))/tau_ij(j)
          ! Prevent floating point overflow caused by very low opacity (tau < 1e-8)
          else if (abs(tau_ij(j)).lt.1.0D-8) then
            beta_ij_ray(j)=1.0D0
          ! For all other cases use the standard escape probability formalism
          else
            beta_ij_ray(j)=(1.0D0-EXP(-tau_ij(j)))/tau_ij(j)
          endif
	 !=============
	 tau(ilevel,jlevel,j)=tau_ij(j)
	 bbeta(ilevel,jlevel,j)=beta_ij_ray(j)
	 !=============
         enddo !j=0,nrays-1

         beta_ij_sum=sum(beta_ij_ray)
         !calculation of average beta_ij in the origin grid point
#ifdef PSEUDO_1D
         beta_ij = beta_ij_sum
#elif PSEUDO_2D
         beta_ij = beta_ij_sum / 4.
#else
         beta_ij = beta_ij_sum / real(nrays,kind=DP) 
#endif


1 continue
         bbeta(ilevel,jlevel,6) = beta_ij
         line(ilevel,jlevel) = A_COEFFS(ilevel,jlevel)*HP*frequencies(ilevel,jlevel) * &
                             & s_pop(ilevel)*beta_ij*(S_ij-BB_ij)/S_ij
         cooling_rate = cooling_rate + line(ilevel,jlevel)
2 continue
         !<J_ij>
         field(ilevel,jlevel) = (1.0D0-beta_ij)*S_ij + beta_ij*BB_ij
         field(jlevel,ilevel) = field(ilevel,jlevel)
         !J_ij(p)

endif !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF MODE: LVG
       enddo !jlevel=1,nlev
     enddo !ilevel=1,nlev
 
3 continue
    !R_IJ CALCULATIONS
    !Update the transition matrix: Rij = Aij + Bij.<J> + Cij				    	 
    DO ilevel=1,NLEV
      DO jlevel=1,NLEV
        TRANSITION(ilevel,jlevel)=A_COEFFS(ilevel,jlevel)&
        & +B_COEFFS(ilevel,jlevel)*FIELD(ilevel,jlevel)&
        & +C_COEFFS(ilevel,jlevel)
        IF(ABS(TRANSITION(ilevel,jlevel)).LT.1.0D-50) TRANSITION(ilevel,jlevel)=0.0D0

      ENDDO !jlevel=1,nlev
    ENDDO !ilevel=1,nlev

  return

end subroutine escape_probability
