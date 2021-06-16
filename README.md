# 3D-PDR

escape_probability subroutine

escape_probability(transition, dust_temperature, nrays, nlev,nfreq, A_COEFFS, B_COEFFS, C_COEFFS,
	frequencies,s_evalpop, maxpoints, Tguess, v_turb, s_jjr, s_pop, s_evalpoint,
	weights,cooling_rate,line,line_profile,tau,coolant,density,metallicity,bbeta)

Level population calculation

	Parameters:
	...

	nlev : integer(kind=i4b), intent(in)
	Number of levels

	nfreq : integer(kind=i4b), intent(in)
	Number of ticks on frequency axis defines the resolution of line profiles

	line_profile : real(kind=dp), intent(out), shape(nlev,nlev,nfreq)
	The line profile intensity values for frequencies range from -3*sigma to +3*sigma


