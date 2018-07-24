# Absorption coefficient curve fitter (ACCF)

ACCF is a simple R script that fits experimental optical-absorption coefficient curves to a theoretical model for indirect band-gap semiconductor materials, returns an estimation of the band gap energy and generates a plot of the experimental data and the fitted curves. It can fit data from multiple experiments simultaneously.

## Usage

After downloading the script, you can run it as any other R script. 

The data file must be a csv file containing the incident radiation energy (in eV) in odd columns and the measured absorption coefficient (in 1/cm) in even columns. Columns 1 and 2 must correspond to data of the fist optical-absorption experiment, columns 3 and 4 correspond  to the second experiment and so on.

### From the R console:

	source(file="ACCF.R")

The experimental data should be located in the same directory as the ACCF.R script in a file called "data.csv".

### From the linux terminal:

	Rscript ACCF.R /PATH/TO/DATAFILE.csv

If no arguments are given, the experimental data is assumed to be in the same directory as ACCF.R.

## Theoretical model

For low energy radiation (lower than the band gap energy) the absorption coefficient is characterized by the Urbach's rule.

	alpha1*(exp((e-eg)/en))  for  e < ec  	(1)

The absorption at high energies for indirect band-gap semiconductors is described by equation (2). This is an approximation that doesn't take into account the energy of the indirect transition photon, which is usually a goop approximation at ambient temperature.

	alpha2*(e-eg)^2  for  e > ec  		(2)

Since both equations are describing two regions of the same curve, the must fulfill the condition for continuity at one point (*ec*). In other words, both equations and their derivatives must have the same value at *ec*. After a little bit of algebra we obtain the following conditions:

	ec = 2*en + eg				(3)
	en = sqrt(alpha1/(4*alpha2))*exp(1)	(4)


## The algorithm

1. A first approximation of *eg* is obtained by a line fitting of the high energy region data. The intersection of this line and the x-axis is our first estimation of *eg*.
2. A small value of en is assumed (0.1 eV) since *ec* and *eg* have a similar value.
3. *alpha2* and *eg* are estimated with a non-linear fit to equation (2).
4. *alpha1* is estimated with a non-linear fit to equation (1).
5. A new estimation of en and *ec* are obtained with equations (3) and (4).
6. We return to the third step and repeat the process until the difference between estimated values of consecutive  iterations is sufficiently small.


