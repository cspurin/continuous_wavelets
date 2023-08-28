#!/usr/bin/env python3.8

########################################################################

import numpy as np
import sys
import os
import inspect
import matplotlib.pyplot as plt
import pycwt
import shutil

from processing_wav import pad, autoscales, recon, fourier_from_scales, icwt_fixed

##### EXAMPLE VARIABLES TO PASS IN #####
#	# input transect
#	infile = './amphibia_richness_americas.dat' or './carnivora_richness_americas.dat'
#	dt = 10000.
#	mirror = True
#	cut1 = 100000.
#	cut2 = 1000000.
#	write_output = True
#	om0 = 6
#	dj=0.1
#	wf = 'morlet'

def run_full_wavelet_analysis(infile, dt=10000., mirror=True, cut1=None, cut2=None, write_output=True, wf='morlet', dj=0.1, om0=6, normmean=True, mirrormethod=1):

	outdir = '.'

	p = int(om0)
	if wf == 'morlet':
		fullwavelet = pycwt.wavelet.Morlet(p)
	elif wf == 'dog':
		fullwavelet = pycwt.wavelet.DOG(p)
	elif wf == 'paul':
		fullwavelet = pycwt.wavelet.Paul(p)
	else:
		sys.exit('Could not recognise desired mother wavelet, exiting...')

	## read in data ##
	xtemp = np.loadtxt(infile)
	if len(np.shape(xtemp)) > 1:
		x = [ i[-1] for i in xtemp ]
	else:
		x = xtemp
	x_orig = x

	## mirror/flip and save mean of result ##
	if mirror == True :
		xmax = np.amax(x)

		if mirrormethod == 1:
			xstart = x[0]
			xend = x[-1]
			x_rev = x[::-1]
			x_rev_flip = np.multiply(x_rev,-1.) + xend*2.
			x_flip = np.multiply(x,-1.) + xstart*2
			x = np.concatenate((x_orig,x_rev_flip,x_flip,x_rev,x_orig,x_rev_flip,x_flip,x_rev),axis=0)

		if mirrormethod == 2:
			x_rev = x[::-1]
			x_flip = (x * -1.)
			x_rev_flip = (x_rev * -1.) 
			x = np.concatenate((x,x_rev,x_flip,x_rev_flip,x,x_rev,x_flip,x_rev_flip),axis=0)

		xmirrormean = np.mean(x)

		if normmean == True:
			xmname = outdir + '/xmirror.mean'
			np.savetxt(xmname, np.array([xmirrormean]))
			x = x - xmirrormean
	else:
		if normmean == True:
			x = x - np.mean(x)

	## padding, save signal to be analysed ##
	[x_pad, x_pad_orig] = pad(x,method='zeros')
	x_pad = x # no padding!
	xlen = len(x_pad)
	if write_output == True :
		signame = outdir + '/signal.h'
		outfile = open(signame, 'w')
		np.savetxt(signame, x_pad)

	fft_signal = np.fft.rfft(x_pad)
	fkinv = np.fft.irfft(fft_signal)
	fft_freq = np.fft.rfftfreq(xlen, d=dt)

	fft_freq = fft_freq[1:int((xlen/2)+1)] 
	fft_signal = (2. / xlen) * abs(fft_signal[1:int(xlen/2)+1])
	fft_power  = abs(fft_signal[0:int(xlen/2)+1]) ** 2.

	fftave = np.convolve(fft_power, np.ones(5)/5)

	### saving np fft results separate to those done internally by pycwt (by pyfftw.interfaces.scipy_fftpack - see helpers.py from pycwt)
	fft_freq_np = fft_freq
	fft_power_np = fft_power
	fft_ave_np = fftave
	fkinv_np = fkinv

	if write_output == True:
		fftpowname = outdir + '/fft.pow'
		fftinvname = outdir + '/fft.inv'
		outfile = open(fftpowname, 'w')
		for i in range (0,len(fft_power)):
			print(fft_freq_np[i], fft_power_np[i] , fft_ave_np[i], file=open(fftpowname, 'a'))
		outfile.close()
		
		outfile = open(fftinvname, 'w')
		for i in range (0,len(fkinv_np.real)):
			print(fkinv_np.real[i],file=open(fftinvname, 'a'))
		outfile.close()
 
	## calculate wavelet scales ##
	N=int(x_pad.shape[0])
	N_orig = int(len(x_orig))

	s0 = 2.*dt

	j_full = (1/dj) * np.log(N*dt /s0) / np.log(2.)
	j_orig = (1/dj) * np.log(N_orig*dt /s0) / np.log(2.)

	## perform continuous wavelet transform ##
	cwtX = pycwt.cwt(x_pad, dt, dj=dj, s0=s0, J=j_full, wavelet=fullwavelet)
	cwtX_orig = pycwt.cwt(x_orig, dt, dj=dj, s0=s0, J=j_orig, wavelet=fullwavelet)
	X = cwtX[0]
	sj = cwtX[1]
	sj_orig = cwtX_orig[1]
	coi = cwtX[3]
	fft_power_pycwt = cwtX[4] ** 2.
	fft_freq_pycwt = cwtX[5]
	fftave_pycwt = np.convolve(fft_power_pycwt, np.ones(5)/5)
	power = (np.abs(X))**2.
	sumpow = np.sum(power, axis=1)

	scales = sj
	scales_orig = sj_orig

	## convert scales to fourier periods ##
	period = fourier_from_scales(scales,wf=wf,p=p)
	scale_len = len(scales)
	if write_output == True :
		spname = outdir + '/scales.periods'
		outfile = open(spname, 'w')
		for i in range (0, scale_len):
			print(i, scales[i], period[i], 1./period[i], file=open(spname, 'a'))
		outfile.close()

	period_orig = fourier_from_scales(scales_orig,wf=wf,p=p)
	scale_len_orig = len(scales_orig)

	if write_output == True :
		sporigname = outdir + '/scales_orig.periods'
		outfile = open(sporigname, 'w')
		for i in range (0, scale_len_orig):
			print(i, scales_orig[i], period_orig[i], 1./period_orig[i], file=open(sporigname, 'a'))
		outfile.close()

	sumpowname = outdir + '/sumpower.fp'
	if write_output == True :
		outfile = open(sumpowname, 'w')
		for i in range (0, len(sj)):
			print(sumpow[i] / xlen , period[i], 1./period[i], scales[i], sumpow[i] / (scales[i] * xlen),file=open(sumpowname, 'a'))
		outfile.close()

	if write_output == True:
		fftpowname = outdir + '/fft_pycwt.pow'
		outfile = open(fftpowname, 'w')
		for i in range (0,len(fft_power_pycwt)):
			print(fft_freq_pycwt[i], fft_power_pycwt[i] , fftave_pycwt[i], file=open(fftpowname, 'a'))
		outfile.close()


	## Gabor limit check ##
	x_array = dt*np.arange(1, len(x_pad)+1)
	freqs = 1/period
	s1 = np.std(x_array)
	s2 = np.std(freqs)
	gtest = s1*s2
	if gtest < 1/(4*np.pi):
		print('\nWarning, signal spacing/frequency choice may not conform to Gabor limit...\n')

	# perform inverse transform (to check correct scales have been used)...
	## NOTE there was an error in previous versions of pycwt - check wavelet.py line 170 in icwt - sj should be square rooted on bottom of iW = ...
	## AND brackets should be added... should read:
	# iW = (dj * np.sqrt(dt) / (wavelet.cdelta * wavelet.psi(0)) *
	#          (np.real(W) / np.sqrt(sj)).sum(axis=0))
	# and then will work fine
	# INCLUDES calculation of recon factor (from empirical cdelta) unlike mlpy

#	x_icwt = pycwt.icwt(X, sj, dt, dj=dj, wavelet=wf).real
	x_icwt = icwt_fixed(X, sj, dt, dj=dj, wavelet=fullwavelet).real

	icwtname = outdir + '/x_icwt.x'
	if write_output == True :
		np.savetxt(icwtname, x_icwt)

	# calculate mean squared error of icwt... not used currently
	diff = np.sqrt((x_pad - x_icwt)**2.)
	diffmean = np.mean(diff)

	if cut1 is not None and cut2 is None:
		print('Please specify two cut-off wavelengths to calculate ICWTs...')
	elif cut1 is not None and cut2 is not None:
		## use only scales above cut1 or cut2 for reconstruction of signal...
		X_numcols = len(X[0])
		scales_cut1 = []
		scales_cut2 = []
		for i in scales:
			if i > cut1:
				scales_cut1.append(i)
			if i > cut2:
				scales_cut2.append(i)

		ncut1 = len(scales_cut1)
		ncut2 = len(scales_cut2)
		X_cut1 = X[len(scales)-ncut1:len(scales),0:X_numcols]
		X_cut2 = X[len(scales)-ncut2:len(scales),0:X_numcols]

		x_icwt_part1 = icwt_fixed(X_cut1, np.array(scales_cut1), dt, dj=dj, wavelet=fullwavelet).real
		x_icwt_part2 = icwt_fixed(X_cut2, np.array(scales_cut2), dt, dj=dj, wavelet=fullwavelet).real

		if write_output == True :
			op1name = outdir + '/x_icwt_part1.x'
			np.savetxt(op1name, x_icwt_part1.real)	
			op2name = outdir + '/x_icwt_part2.x'
			np.savetxt(op2name, x_icwt_part2.real)


		# calculate mean squared error of both filtered icwts... not used currently
		diff_cut1 = np.sqrt((x_pad - x_icwt_part1.real)**2.)
		diffmean_cut1 = np.mean(diff_cut1)
		diff_cut2 = np.sqrt((x_pad - x_icwt_part2.real)**2.)
		diffmean_cut2 = np.mean(diff_cut2)
	else:
		print('Not calculating filtered ICWTs...')

	## output power normalised by scale [e.g. Liu et al., 2007] as list of matrix elements and their values
	numrows = len(power)
	numcols = len(power[0])
    
    
    # attempting matplotlib contour plot
    
       
	if write_output == True :
		print('Writing out to file.txt (can be large!)...')
		print('Format = n, t, power, period, 1/period, rectified power')
		filename = outdir + '/file.txt'
		outfile = open(filename, 'w')
		for i in range (0,numrows):
			for j in range (0,numcols):
				print(i, j*dt, power[i,j], period[i], 1./period[i], power[i,j]/scales[i], file=open(filename, 'a'))
		outfile.close()
		print("Wrote to file.txt succesfully.")

	return cwtX, scales, x_pad, xlen, period, power



def run_double_wavelet_analysis(infile1, infile2, dt=10000., mirror=True, cut1=None, cut2=None, write_output=True, wf='morlet', dj=0.1, om0=6, normmean=True, mirrormethod=1):

	print('\nWarning, currently run_double_wavelet_analysis will overwrite any single wavelet analysis results in the same directory!\n')
	print('Also, not saving cut-off inverse transforms...')
	try:
		dat1 = np.loadtxt(infile1)
		dat2 = np.loadtxt(infile2)
		len1 = len(dat1)
		len2 = len(dat2)
		if len1 != len2:
			raise ValueError('The two input signals need to be the same length!')
		if len1 == 0. or len2 == 0.:
			raise ValueError('The signals cannot have zero length.')
	except (ValueError):
		exit('Exited.')

	cwtX, scales, x_pad, xlen, period = run_full_wavelet_analysis(infile1, dt=dt, mirror=mirror, cut1=cut1, cut2=cut2, write_output=write_output, wf=wf, dj=dj, om0=om0, normmean=normmean, mirrormethod=mirrormethod)

	shutil.copyfile('./xmirror.mean', './xmirror1.mean')
	shutil.copyfile('./signal.h', './signal1.h')
	shutil.copyfile('./fft.pow', './fft1.pow')
	shutil.copyfile('./fft_pycwt.pow', './fft1_pycwt.pow')
	shutil.copyfile('./fft.inv', './fft1.inv')
	shutil.copyfile('./scales.periods', './scales1.periods')
	shutil.copyfile('./scales_orig.periods', './scales_orig1.periods')
	shutil.copyfile('./sumpower.fp', './sumpower1.fp')
	shutil.copyfile('./x_icwt.x', './x_icwt1.x')
	shutil.copyfile('./file.txt', './file1.txt')

	cwtY, scales, y_pad, ylen, period = run_full_wavelet_analysis(infile2, dt=dt, mirror=mirror, cut1=cut1, cut2=cut2, write_output=write_output, wf=wf, dj=dj, om0=om0, normmean=normmean, mirrormethod=mirrormethod)

	shutil.copyfile('./xmirror.mean', './xmirror2.mean')
	shutil.copyfile('./signal.h', './signal2.h')
	shutil.copyfile('./fft.pow', './fft2.pow')
	shutil.copyfile('./fft_pycwt.pow', './fft2_pycwt.pow')
	shutil.copyfile('./fft.inv', './fft2.inv')
	shutil.copyfile('./scales.periods', './scales2.periods')
	shutil.copyfile('./scales_orig.periods', './scales_orig2.periods')
	shutil.copyfile('./sumpower.fp', './sumpower2.fp')
	shutil.copyfile('./x_icwt.x', './x_icwt2.x')
	shutil.copyfile('./file.txt', './file2.txt')

	X = cwtX[0]
	Y = cwtY[0]
	sjx = cwtX[1]
	sjy = cwtY[1]
	coix = cwtX[3]
	coiy = cwtY[3]
	powerx = (np.abs(X))**2.
	powery = (np.abs(Y))**2.
	sumpowx = np.sum(powerx, axis=1)
	sumpowy = np.sum(powery, axis=1)
	rectpowx = []
	distavgrectpowx = []
	rectpowy = []
	distavgrectpowy = []
	for i in range(len(scales)):
		rectpowx.append(powerx[i]/scales[i])
		distavgrectpowx.append(sumpowx[i] / (scales[i] * xlen ))
		rectpowy.append(powery[i]/scales[i])
		distavgrectpowy.append(sumpowy[i] / (scales[i] * ylen ))

	## cross wavelets and coherence
	Wxy = X * np.conjugate(Y)
	xypower = np.abs(Wxy)
	Wyy = Y * np.conjugate(Y)
	yypower = np.abs(Wyy)
	phasexy = np.angle(Wxy,deg=True)
	
	if wf == 'morlet':
		fullwavelet = pycwt.wavelet.Morlet(om0)
	elif wf == 'dog':
		fullwavelet = pycwt.wavelet.DOG(om0)
	elif wf == 'paul':
		fullwavelet = pycwt.wavelet.Paul(om0)
	else:
		sys.exit('Could not recognise desired mother wavelet, exiting...')

	print("\nCalculating coherence with pycwt...\n")
	R2ns, R_phase, coi, freq, sig = pycwt.wct(x_pad,y_pad,dt,dj=dj,s0=scales[0],J=len(scales)-1,sig=False, wavelet=fullwavelet)



	# save XWT, coherence and phase #
	if write_output == True:
		xwsumpow = np.sum(xypower, axis=1)
		oname7 = "./sumpower_xw.fp"
		outfile = open(oname7, 'w')
		for i in range (0, len(scales)):
			print(xwsumpow[i] / xlen , period[i], 1./period[i], scales[i], xwsumpow[i] / (scales[i] * xlen),file=open(oname7, 'a'))
		outfile.close()

		sumcoh = np.sum(R2ns, axis=1)
		oname7 = "./sumcoh.fp"
		outfile = open(oname7, 'w')
		for i in range (0, len(scales)):
			print(sumcoh[i] / xlen , period[i], 1./period[i], scales[i], sumcoh[i] / (scales[i] * xlen ),file=open(oname7, 'a'))
		outfile.close()

		oname = "./file_xw.txt"
		oname2 = "./phase.txt"
		oname3 = "./file_coh.txt"
		outfile = open(oname, 'w')
		outfile2 = open(oname2, 'w')
		outfile3 = open(oname3, 'w')
		numrows = len(powery)
		numcols = len(powery[0])
		for i in range (0,numrows):
			for j in range (0,numcols):
				print(i, j*dt, xypower[i,j], period[i], 1./period[i], file=open(oname, 'a'))
				print(i, j*dt, phasexy[i,j], period[i], 1./period[i], file=open(oname2, 'a'))
				print(i, j*dt, R2ns[i,j], period[i], 1./period[i], file=open(oname3, 'a'))
		print("Wrote to file_xw.txt, phase.txt and file_coh.txt succesfully...")
		outfile.close()
		outfile2.close()
		outfile3.close()



