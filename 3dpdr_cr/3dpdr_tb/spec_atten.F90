module spec_module
		use definitions
    implicit none
    integer(kind=4),parameter :: ERES = 40
    real(kind = DP), parameter :: EMIN = 100.0
    real(kind = DP), parameter :: EMAX = 2.0e9
    real(kind = DP), dimension(ERES) :: rangeFunc
    real(kind = DP), dimension(ERES) :: eArr
    real(kind = DP), allocatable, dimension(:) :: eloss
		real(kind = DP), allocatable, dimension(:) :: lloss
    !character, allocatable, dimension(:) :: lossFile


    contains
        function logspace(xmin, xmax, n) !Create a log-spaced vector
            implicit none
            integer, intent(in) :: n
            integer :: i
            real(kind = DP), intent(in) :: xmin, xmax
            real(kind = DP), dimension(n) :: logspace
            do i=1,n
                logspace(i) = (xmax - xmin)*real(i-1) / real(n-1) + xmin
								logspace(i) = 10.0**logspace(i)
            end do
        end function

        subroutine CRsetup() !Setup for the attenuation
            implicit none
            integer :: i, nloss, j, k
            real(kind = DP) :: logEmin, logEmax
            real(kind = DP), allocatable, dimension(:) :: eprime
						real(kind = DP), allocatable, dimension(:) :: lprime
            logEmin = log10(EMIN)
            logEmax = log10(EMAX)

            eArr = logspace(logEmin, logEmax, ERES)

						!WRITE(*,*) ""
						!WRITE(*,*) "TESTING SETUP"

            OPEN(10, action='read', file="LE_loss_p_2.txt")
            READ(10, *) nloss
            allocate(eloss(nloss))
            allocate(lloss(nloss))
            do i=1,nloss
                READ(10, *) eloss(i), lloss(i)
            end do
						CLOSE(10)
            lloss = lloss * 1e-16
            do i=1,ERES
                j = MINLOC(eloss, DIM = 1, mask = (eloss > eArr(i)))
                allocate(eprime(j))
                allocate(lprime(j))
                eprime = eloss(1:j)
                lprime = 1.0/lloss(1:j)
                !integrate here
								!WRITE(*,*) eprime(1), eprime(j)
								!WRITE(*,*) "EPRIME = "
								!do k=1,j
									!WRITE(*,*) eprime(k), lprime(k)
								!end do
                call avint(SIZE(eprime), eprime, lprime, eprime(1), eprime(j), rangeFunc(i))
								deallocate(eprime)
								deallocate(lprime)
            end do
            OPEN(12, action='write', file="range_test.txt")
						WRITE(12,*) "!Energy, density*Range\n"
            do i=1,ERES
                WRITE(12, *) eArr(i), rangeFunc(i)
            end do
						CLOSE(12)
        end subroutine CRsetup

        subroutine spec_atten(espec, j0, ja, eret, NCOL)
            real(kind = DP),intent(in), dimension(:) :: j0 !Input spectrum
            real(kind = DP),intent(in), dimension(:) :: espec
            real(kind = DP),intent(in) :: NCOL !Column density
            real(kind = DP),intent(out), dimension(:) :: eret
            real(kind = DP),intent(out), dimension(:) :: ja !Attenuated spectrum
            real(kind = DP) :: RI, DR, E0, L1, L2
            real(kind = DP), dimension(2) :: LTMP, ETMP
            integer :: i, n, j
            n = SIZE(j0)

						!allocate(eret(ERES))
						!allocate(ja(ERES))

						do i=1, ERES
            	eret(i) = eArr(i)
						end do
            do i=1,ERES
                call interp_linear ( 1, SIZE(rangeFunc), log10(eArr), log10(rangeFunc), 1, log10(eArr(i)), RI )
                RI = 10.0**RI
                DR = NCOL + RI
                call interp_linear ( 1, SIZE(rangeFunc), log10(rangeFunc), log10(eArr), 1, log10(DR), E0 )
                E0 = 10.0**E0
								ETMP(1) = log10(eArr(i))
								ETMP(2) = log10(E0)
                call interp_linear ( 1, SIZE(eloss), log10(eloss), log10(lloss), 2, ETMP, LTMP )
								L1 = 10.0**LTMP(1)
                L2 = 10.0**LTMP(2)
                j = MINLOC(espec, DIM=1, MASK=(espec > E0)) - 1
                ja(i) = j0(j)*(L2/L1)
								!WRITE(*,*) eArr(i), RI, DR, E0, L1, L2, j, j0(j), ja(i)
            end do

        end subroutine spec_atten
				subroutine ionRate(e, j, zeta)
					implicit none
					real(kind = DP), dimension(:), intent(in) :: e,j
					real(kind = DP), intent(out) :: zeta !the CRIR

					real(kind = DP), parameter :: a0 = 5.29177211e-9 !cm^2 - bohm radius
					real(kind = DP), parameter :: IP = 13.598 !eV - H ionization potential
					real(kind = DP), parameter :: memp = 5.44617e-4 !me/mp
					real(kind = DP), parameter :: A = 0.71
					real(kind = DP), parameter :: B = 1.63
					real(kind = DP), parameter :: C = 0.51
					real(kind = DP), parameter :: D = 1.24
					real(kind = DP), parameter :: mpc2 = 0.938 !GeV - mass of proton
					real(kind = DP), parameter :: me = 9.109383e-28 !g - mass of electron
					real(kind = DP), parameter :: mh = 1.67372e-24 !g - mass of hydrogen/proton
					real(kind = DP), parameter :: ev_to_ergs = 1.602e-12 !ev to ergs
					real(kind = DP), parameter :: pi = 3.14159265358979323
					real(kind = DP), parameter :: cl = 2.998e10 !cm/s - speed of light
					integer :: n, i
					real(kind = DP), dimension(:), allocatable :: xarr, sigmal, sigmah, sigmap, vrel2, integrand

					n = SIZE(e)
					allocate(xarr(n))
					allocate(sigmal(n))
					allocate(sigmah(n))
					allocate(sigmap(n))
					allocate(vrel2(n))
					allocate(integrand(n))

					!Cross section WITH relativistic correction
					vrel2 = (1.0 - (mpc2*1E9/(earr + mpc2*1E9))**2)*cl**2
					xarr = (mh/2.0)*(me/(mh*IP*ev_to_ergs))*vrel2
					sigmal = 4.0*pi*a0**2*C*xarr**D
					sigmah = 4.0*pi*a0**2*(A*log(1.0+xarr) + B)*(1.0/xarr)
					sigmap = 1.0/((1.0/sigmal) + (1.0/sigmah))

					!OPEN(13, action='write', file='outCross.txt')
					!do i=1, n
					!	WRITE(13, *), e(i), sigmap(i)
					!end do
					!CLOSE(13)

					integrand = j*sigmap
					!OPEN(14, action='write', file='outIntegrand.txt')
					!do i=1, n
					!	WRITE(14, *), e(i), integrand(i)
					!end do
					!CLOSE(14)
					!OPEN(15, action='write', file='outTmp.txt')
					!do i =1, n
					!	WRITE(15, *), e(i), j(i)
					!end do
					!CLOSE(15)
					call avint(n, e, integrand, e(1), e(n), zeta)
					zeta = 2.0*pi*zeta

				end subroutine ionRate
end module spec_module
