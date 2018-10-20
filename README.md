# particleDetection

>> Gabriel Moreira, 2018

This research project was proposed by ISAE-Supaero in partnership with Nuclétudes.
The goal was to use CMOS Image sensors (CIS) to identify ionizing particles (alpha, X-rays, cosmic-rays).

>> subDeviations.m
This function calculates the standard deviations of the submatrices of a given matrix. Since the fixed pattern noise has a small standard deviation, this criterion can be used to distinguish between particle events and the former.

>> main.m
Opens a .JPG image with a specified integration time (see images provided) and retrives the RAW data. The image processing algorithm is comprised of a Wiener filter followed by a convolution with an event kernel whose matrix is also provided. After these steps the frame obtained can be easily used for event detection using just a threshold or the subDeviations values. The results plotted include:

1.) initial frame
2.) frame resulting from the wiener filter
3.) the convoluted frameaverage 
4.) Histogram of the average energy per event in eV
5.) The representations of all the events identified

>> NOTE: All the constants used in the calculations were determined and calibrated for a SONY IMX219 image sensor.

>> For more information please refer to the article provided: CIS_for_Particle_Detection_Bielawski_Li_Moreira.pdf
