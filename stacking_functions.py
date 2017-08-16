import sys
import subprocess
import numpy as np
import healpy as hp

def check_dir(dirname):
	if subprocess.check_output("if test -d {}; then echo 'exist'; fi".format(dirname), shell=True) == "exist\n":
		return True
	else:
		return False

def check_file(fname):
	if subprocess.check_output("if [ -f {} ]; then echo 'exist'; fi".format(fname), shell=True) == "exist\n":
		return True
	else:
		return False

def read_it(rtype, fname, nside):
	if rtype == hpix:
		if fname[-3:] == 'npy' or fname[-3:] == 'npz':
			data = np.load(fname)
		else:
			try:
				data = np.loadtxt(fname)
			except (ValueError):
				print "Error loading file. Ensure the file consists of a single column of floats with comments commented out using '#'s."
				sys.exit()
		pix = hp.ring2nest(nside, data)

	elif rtype == xyz:
		if fname[-3:] == 'npy' or fname[-3:] == 'npz':
			data = np.load(fname)
		else:
			try:
				data = np.loadtxt(fname)
			except (ValueError):
				print "Error loading file. Ensure the file consists of three, space delimited columns of floats with comments commented out using '#'s."	
				sys.exit()
		if data.shape[1] == 3:
			x, y, z = data[:,0], data[:,1], data[:,2]
			pix = hp.vec2pix(nside, x, y, z, nest=True)
		else:
			print "Error loading file: incorrect number of columns for xyz file."
			sys.exit()
	
	elif rtype == lonlat:
		if fname[-3:] == 'npy' or fname[-3:] == 'npz':
			data = np.load(fname)
		else:
			try:
				data = np.loadtxt(fname)
			except (ValueError):
				print "Error loading file. Ensure the file consists of two, space delimited columns of floats with comments commented out using '#'s."	
				sys.exit()
		if data.shape[1] == 2:
			theta, phi = data[:,0], data[:,1]
			if max(theta) <= np.pi and max(phi) <= 2*np.pi:
				pix = hp.ang2pix(nside, theta, phi, nest=True)
			else:
				pix = hp.ang2pix(nside, theta, phi, nest=True, lonlat=True)
		else:
			print "Error loading file: incorrect number of columns for thetaphi file."
			sys.exit()
	return theta, phi, pix



def txt2bin_file_writer(fname, outdir, header, theta, phi, pix, peak_type, nu, masked=False, mask=None):
	if peak_type == 'hotPeaks':
		genre = 1
	elif peak_type == 'coldPeaks':
		genre = 2
	elif peak_type == 'hotPeaksOrient':
		genre = 4
	elif peak_type == 'coldPeaksOrient':
		genre = 5
	elif peak_type == 'randomHotSpots':
		genre = 13
	elif peak_type == 'randomHotSpotsOrient':
		genre = 14
	elif peak_type == 'randomColdSpots':
		genre = 15
	elif peak_type == 'randomColdSpotsOrient':
		genre = 16
	
	print "\nReading in map...\n"

	m = hp.read_map(fname, verbose=False, field=None)
	
	count1 = 1
	while fname[-count1] != '.':
		count1 += 1
	count2 = count1
	while fname[-count2] != '/':
		count2 += 1
	bname = outdir+fname[-count2:-count1]

	if masked:
		mask_int = True
		mask = hp.read_map(mask, verbose=False)
	else:
		mask_int = False
	mask_pol = False
	nmaps = header['TFIELDS']
	nside = header['NSIDE']
	order = header['ORDERING']
	nside2 = 0

	print "Obtaining peak values...\n"

	if nmaps == 3:
		index_I = 1
		index_Q = 2
		index_U = 3
		if order == 'RING':
			rpix = hp.nest2ring(nside, pix)
			if masked:
				rpix = rpix[mask[rpix]!=0]
			I, Q, U = m[0][rpix], m[1][rpix], m[2][rpix]
		else:
			if masked:
				pix = pix[mask[pix]!=0]
			I, Q, U = m[0][pix], m[1][pix], m[2][pix]
		threshold_option = 6
		P = np.sqrt(Q**2 + U**2)
		sigma_P = np.std(P)
	elif nmaps == 1:
		index_I = 1
		index_Q = 0
		index_U = 0
		if order == 'RING':
			rpix = hp.nest2ring(nside, pix)
			I = m[rpix]
		else:
			I = m[pix]
		threshold_option = 4
		sigma_P = 0
		
	index_L = 0
	index_peak = 1
	caption = "T maxima, $Q_{nabla^2T}U_nabla^2T}$ oriented"
	sigma_I = np.std(I)
	sigma_L = 0
	I_lower = 0
	I_upper = 1.e20
	L_lower = -1.e20
	L_upper = 1.e20
	P2_lower = 0.
	P2_upper = np.inf
	I_lower_nu = nu
	I_upper_nu = 1.e20
	L_lower_nu = -1.e20
	L_upper_nu = 1.e20
	P_lower_nu = 0.
	P_upper_nu = 1.e20
	P2byI2_lower = 0.
	P2byI2_upper = 1.e15
	P2byL2_lower = 0.
	P2byL2_upper = 1.e15
	abs_threshold = False
	norm_to_corr = False
	norm_power = 0
	addpi = True  
	angzero = False
	nested = True
	n = len(I)

	print "Writing {} peaks to {}.dat\n".format(n, bname)

	f = open('{}.txt'.format(bname), 'w')
	if nmaps == 3:
		for j in range(n-1):
			f.write('{} {} {} {} {} {}\n'.format(pix[j], theta[j], phi[j], I[j], Q[j], U[j]))
		f.write('{} {} {} {} {} {}'.format(pix[-1], theta[-1], phi[-1], I[-1], Q[-1], U[-1]))
	f.close()

	h = open('{}_header.txt'.format(bname), 'w')
	h.write('{} {}\n'.format(mask_int, mask_pol))
	h.write('{} {} {} {} {} {} {} {} {}\n'.format(genre, nmaps, nside, nside2, index_peak, index_I, index_Q, index_U, index_L))
	h.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(I_lower, I_upper, L_lower, L_upper, P2_lower, P2_upper, I_lower_nu, I_upper_nu, L_lower_nu, L_upper_nu, P_lower_nu, P_upper_nu, P2byI2_lower, P2byI2_upper, P2byL2_lower, P2byL2_upper, norm_power))
	h.write('{}\n'.format(caption))
	h.write('{} {} {} {} {} {}\n'.format(threshold_option, abs_threshold, norm_to_corr, addpi, nested, angzero))
	h.write('{}\n'.format(n))
	h.write('{} {} {}'.format(sigma_I, sigma_L, sigma_P))
	h.close()

	subprocess.Popen("export OASTACK_FNAME={}; ./custom_peaks".format(bname), stdout=subprocess.PIPE, shell=True)



	
	 
	







