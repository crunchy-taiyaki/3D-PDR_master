## Cosmic Ray Attenuation Module
### Introduction
In one-dimension, the code is now able to calculate the cosmic ray ionizate rate in-situ from user defined input cosmic ray spectrum. The attenuation follows the prescription of Padovani+2009 using the 'continuous-slowing-down approximation (CSDA)' also known as the 'continuous energy loss regime'. This approximation will break down towards very high column densities (>100 g/cm^2) (see Padovani+2018 for details). The full details of the implementation are in Gaches+2019a (submitted). In one-dimension, the user defines cosmic ray spectra on either of the available surfaces. The code calculates the cosmic ray attenuation using the in-situ calculated molecular hydrogen from a user-given loss function.

### Outputs
The chemistry code includes two main cosmic-ray related output files:
- zeta.txt : This text file contains three columns: (x, CRIR/1.3E-17, NCOL) where CRIR is the cosmic ray ionization rate and NCOL is the column density of H2
- OUT.CR.fin : This file is the main output file for the cosmic ray spectrum. The first line is (NENE, Energy_array), where NENE is the number of energy bins and the Energy_array is NENE entries long. The following lines are (x, CR_spectrum) where the CR_spectrum is the spectrum at point x. If NENE = 40 and the 1D domain has 1000 points, the file has 1001 rows and 41 columns.

### User-defined Inputs
The user can define several different inputs. The first of which is their own loss function. The public distribution contains two pre-defined loss function input files: LE_loss_p.txt and LE_loss_p_2.txt. They are the same loss function but at different energy resoutions. The use can define their own where with the format:
bins
E, LE/1E-16
These functions are used to interpolate over, so the number of bins can be (and should be) higher than NENE. LE_loss_p2.txt has 256 bins. If you define your own loss function, change the input file name in spec_atten.F90 in the CRsetup() function.

The next inputs are all defined in the params.dat file. An example is given below for two sources:
```
===========================|
Cosmic Ray Sources         |
===========================|
2	 		   !Number of CR Sources
40			   !Number of CR energy bins
1			   !Index for CR transport - 1 = diffusive, 2 = rectilinear
SRC			   !ISO-tropic or Source (SRC)
Taurus_spec.dat		   !Input CR Spectrum as (E, F(E))
0.1			   !Radial scaling in (pc)
-1			   !CR Surface - (1 for +x direction, -1 for -x direction)
ISO			   !ISO-tropic or Source (SRC)
IS_L_spectrum.dat	   !Input CR Spectrum as (E, F(E))
1			   !CR Surface - (1 for +x direction, -1 for -x direction)
```
In this case, Taurus_spec.dat is a spectrum for cosmic rays from a Taurus-like protocluster which is embedded in the center (-1) of the cloud. The other source is IS_L_spectrum.dat which is an external interstellar spectrum at the external surface (1). The index for CR transport is for sources only, where the physical transport scales as (rscale/r)^a. The public distribution contains three example spectra:
- IS_L_spectrum.dat and IS_H_spectrum.dat are the Low and High interstellar spectrum from Ivlev+2015
- Taurus_spec.dat is an example protocluster Spectrum

User-defined spectrum must have the format: (e, je) where je is in units of particles/eV/s/cm^2/sr. The file IS_Spectra.pdf shows a plot of the two different provided interstellar spectra.
