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
	return pix



def txt2bin_file_writer(fname, header, pix, peak_type, nu, masked=False):
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
	m = hp.read_map(fname, verbose=False, field=None)
	if masked:
		mask_int = 'T'
	else:
		mask_int = 'F'
	mask_pol = 'F'
	genre = genre
	nmaps = header['TFIELDS']
	nside = header['NSIDE']
	nside2 = 0
	if nmaps == 3:
		index_I = 1
		index_Q = 2
		index_U = 3
		I, Q, U = m[0][pix], m[1][pix], m[2][pix]
		threshold_option = 6
	elif nmaps == 1:
		index_I = 1
		index_Q = 0
		index_U = 0
		I = m[pix]
		threshold_option = 4
	index_L = 0
	index_peak = 1
	caption = 'T'
	sigma_I = np.std(I)
	sigma_L = 0
	sigma_P = 0
	I_lower = np.min(I)
	I_upper = 1.e20
	L_lower = -1.e20
	L_upper = 1.e20
	P2_lower = 0.
	P2_upper = 1.e15
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
	abs_threshold = 'F'
	norm_to_corr = 'F'
	norm_power = 0
	addpi = 'T'  
	angzero = 'F'
	nested = 'T'
	







