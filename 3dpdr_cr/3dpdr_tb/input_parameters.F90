
SUBROUTINE readparams
 
!T.Bisbas, T.Bell
  use definitions
  use healpix_types
  use maincode_module, only : input, level, Tguess, v_turb_inp, iterstep,&
                              & theta_crit, ITERTOT, output, rho_min, rho_max,&
                              & fieldchoice, Gext, AV_fac, UV_fac, nspec, nreac, &
                              & numCRsources, crNene, crfieldchoice, crInputFiles, & !B. Gaches CR input
                              & crSpectrum, crEnergy, crRadialScaling, & !B. Gaches CR input
															& crEnergy_A, crSpectrum_A, crDirection, crindex, & !B. Gaches CR input
                              & CII_NLEV, CII_NTEMP, CI_NLEV, CI_NTEMP, &
                              & OI_NLEV, OI_NTEMP, C12O_NLEV, C12O_NTEMP, maxpoints, &
                              & Tlow0, Thigh0, Tmin, Tmax, Fcrit, Tdiff, dust_temperature,writeiterations,&
                              & chemiterations,  end_time, C12Oinput, CIIinput, CIinput, OIinput
  use global_module, only : g2d, metallicity, omega, grain_radius

  integer::i, j
  real(kind=dp):: dummy
	real(kind=dp):: dum1, dum2
  
  open(unit=12,file='params.dat',status='old')

  read(12,*); read(12,*); read(12,*)
  read(12,*) input
  read(12,*) output
  read(12,*) level
  read(12,*) theta_crit
  read(12,*) chemiterations
  read(12,*) ITERTOT
  read(12,*) rho_min
  read(12,*) rho_max
  read(12,*) v_turb_inp
  read(12,*) dust_temperature
  read(12,*) end_time
  read(12,*) g2d
  read(12,*) metallicity
  read(12,*) omega
  read(12,*) grain_radius
  read(12,*); read(12,*); read(12,*)
  read(12,*) C12Oinput
  read(12,*) CIIinput
  read(12,*) CIinput
  read(12,*) OIinput
  read(12,*); read(12,*); read(12,*)
  read(12,*) Tguess
  read(12,*) Tlow0
  read(12,*) Thigh0
  read(12,*) Tmin
  read(12,*) Tmax
  read(12,*) Fcrit
  read(12,*) Tdiff
  read(12,*); read(12,*); read(12,*)
  read(12,*) fieldchoice
  read(12,*) Gext
  !B. Gaches begin Cosmic ray input
  !Format:
  !===========
  !CR Input
  !===========
  !1                !Number of CR Sources
  !40               !Number of energy bins
	!2								!CR transport index
  !ISO              !ISO-tropic or Source (SRC)
  !input_ISO.txt    !Input file
  !SRC              !ISO-tropic or SOurce (SRC)
  !input_SRC.txt    !Input file
  !1d-2             !Radial scaling for flux (in pc)
  read(12,*); read(12,*); read(12,*)
	read(12, *) numCRsources
	read(12, *) crNene
	read(12, *) crindex
	allocate(crfieldchoice(numCRsources))
	allocate(crInputFiles(numCRsources))
	allocate(crRadialScaling(numCRsources))
	allocate(crDirection(numCRsources))
	do i=1, numCRsources
		read(12, *) crfieldchoice(i)
		read(12, *) crInputFiles(i)
		if (crfieldchoice(i).EQ."ISO") crRadialScaling(i) = 0.0
		if (crfieldchoice(i).EQ.'SRC') then
				read(12,*) dum1
				crRadialScaling(i) = dum1
		endif
		read(12,*) crDirection(i)
	end do
	allocate(crSpectrum(numCRsources, crNene))
	allocate(crEnergy(numCRsources, crNene))
	allocate(crEnergy_A(crNene))
	allocate(crSpectrum_A(crNene))
	
	do i=1, numCRsources
		open(unit=96, file=crInputFiles(i), status='old')
		do j=1, crNene
			read(96, *) dum1, dum2
			crEnergy(i, j) = dum1
			crSpectrum(i, j) = dum2
		end do
	end do
	!========== B. Gaches - END Cosmic Ray read in and initialization ===========
	
  maxpoints = 1000
  AV_fac = 6.289E-22*metallicity
  UV_fac = 3.02

#ifdef REDUCED
close(1)
open(unit=1,file='species_reduced.d',status='old')
nspec=0
do 
 read(1,*,end=100) dummy
 nspec=nspec+1
enddo
100 continue
close(1)
open(unit=1,file='rates_reduced.d',status='old')
nreac=0
do 
 read(1,*,end=101) dummy
 nreac=nreac+1
enddo
101 continue
close(1)
#endif

#ifdef FULL
close(1)
open(unit=1,file='species_full.d',status='old')
nspec=0
do 
 read(1,*,end=100) dummy
 nspec=nspec+1
enddo
100 continue
close(1)
open(unit=1,file='rates_full.d',status='old')
nreac=0
do 
 read(1,*,end=101) dummy
 nreac=nreac+1
enddo
101 continue
close(1)
#endif

#ifdef MYNETWORK
close(1)
open(unit=1,file='species_mynetwork.d',status='old')
nspec=0
do 
 read(1,*,end=100) dummy
 nspec=nspec+1
enddo
100 continue
close(1)
open(unit=1,file='rates_mynetwork.d',status='old')
nreac=0
do 
 read(1,*,end=101) dummy
 nreac=nreac+1
enddo
101 continue
close(1)
#endif

CII_NLEV = 5
CII_NTEMP = 18
CI_NLEV = 5
CI_NTEMP = 29
OI_NLEV = 5
OI_NTEMP = 27
C12O_NLEV = 41
C12O_NTEMP = 25

close(12)

return
END SUBROUTINE readparams
