import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.io import fits

# import stacking functions
from stacking_functions import *

# find location of COOP and writeto locations
try:
	settings = np.loadtxt('settings.txt', dtype='str', delimiter='\t')
except (IOError):
	print "Cannot find OAStack/settings.txt. Either setup.py has not been run, or OAStack.py has been moved from where setup.py was run"
	sys.exit()
COOP_dir	= settings[0,1]
map_dir		= settings[1,1]
stack_dir	= settings[2,1]
ftcomp		= settings[3,1]

# main menu
print "\nOAStack, an oriented asymmetric stacking manager for COOP. Select an option below:\n"
print "\t[1] Unoriented (random rotation) stacking"
print "\t[2] Oriented stacking"
print "\t[3] Oriented, symmetry breaking stacking"
print "\t[4] Help"
print "\t[0] exit\n"
case1 = raw_input()
while not (case1 == "0" or case1 == "1" or case1 == "2" or case1 == "3" or case1 == "4"):
	print "\nPlease select a valid option:\n"
	case1 = raw_input()

# exit program
if case1 == "0":
	sys.exit()

# unoriented stacking
elif case1 == "1":
	stack_type = 'unoriented'
	orient = 'RANDOM'
	print "\nSelect type of map to be stacked:\n"
	print "\t[1] thermal Sunyaev-Zeldovich (tSZ)"
	print "\t[2] CO emission"
	print "\t[3] weak lensing"
	print "\t[4] dust"
	print "\t[0] exit\n"
	case2 = raw_input()
	while not (case2 == "0" or case2 == "1" or case2 == "2" or case2 == "3" or case3 == "4"):
		print "\nPlease select a valid option:\n"
		case2 = raw_input()

	# exit
	if case2 == "0":
		sys.exit()
	
	# tsZ maps
	elif case2 == "1":
		map_type = 'tSZ'
		print "\nMaps are set to be found in {}".format(map_dir)
		print "Maps are assumed to be in Healpix format. List of potential maps to stack:\n"
		print subprocess.check_output("cd {}; ls".format(map_dir), stderr=subprocess.STDOUT, shell=True)
		print "Select a map to stack with, or change map locations by running setup.py:\n"
		map_name = raw_input()
		while not check_file(map_dir+'/'+map_name):
			print "\nERROR: No such file. Try Again.\n"
			map_name = raw_input()

		print "\nChecking if valid temperature/IQU map..."

		raw_map = fits.open('{}'.format(map_dir+'/'+map_name))
		header	= raw_map[1].header
		
		try:
			pixtype = header['PIXTYPE']
		except (KeyError):
			print "Invalid map: no PIXTYPE specified. Exiting."
			sys.exit()
		if pixtype != 'HEALPIX':
			print "Invalid map: must be Healpix format. Exiting."
			sys.exit()
		try:
			nside = header['NSIDE']
		except (KeyError):
			print "Invalid map: no NSIDE specified. Exiting."
			sys.exit()
		try:
			numcols = header['TFIELDS']
		except (KeyError):
			print "Invalid map: Unknown number of columns. Exiting."

		if numcols == 3:
			print "3 columns detected."
			I, Q, U = header['TTYPE1'].upper(), header['TTYPE2'].upper(), header['TTYPE3'].upper()
			if (I=='I' or I=='ISTOKES' or I=='INTENSITY' or I=='T' or I=='TEMPERATURE') and (Q=='Q' or Q=='QT' or Q=='QSTOKES') and (U=='U' or U=='UT' or U=='USTOKES'):
				print "IQU map detected. Proceeding...\n"
			else:
				print "WARNING: unknown IQU columns names:\n"
				print "\t"+I
				print "\t"+Q
				print "\t"+U
				print "\nProceed? [y/n]\n"
				iqu_proceed = raw_input()
				while not (iqu_proceed == "y" or iqu_proceed == "n"):
					print "\nPlease enter either 'y' or 'n'\n"
					iqu_proceed = raw_input()
				if iqu_proceed == 'n':
					print "\nValid I column names: 'I' 'ISTOKES' 'INTENSITY' 'T' 'TEMPERATURE'"
					print "Valid Q column names: 'Q' 'QT' 'QSTOKES'"
					print "Valid U column names: 'U' 'UT' 'USTOKES'"
					time.sleep(3)
					print "\nExiting."
					sys.exit()
		elif numcols == 2:
			pass


		elif numcols == 1:
			I = header['TTYPE1'].upper()
			if I=='I' or I=='ISTOKES' or I=='INTENSITY' or I=='T' or I=='TEMPERATURE':
				print "Temperature map detected. Proceeding...\n"
			else:
				print "WARNING: unknown columns temperature column name:\n"
				print "\t"+I
				print "\nProceed? [y/n]\n"
				i_proceed = raw_input()
				while not (i_proceed == "y" or i_proceed == "n"):
					print "\nPlease enter either 'y' or 'n'\n"
					i_proceed = raw_input()
				if i_proceed == 'n':
					print "\nValid temperature column names: 'I' 'ISTOKES' 'INTENSITY' 'T' 'TEMPERATURE'"
					time.sleep(3)
					print "\nExiting."
					sys.exit()
		else:
			print "Invalid map: require either 3 columns for an IQU map or 1 column for temperature map. Exiting."
			sys.exit()
			

		print "\nMasking the map? [y/n]\n"
		masked = raw_input()
		while not (masked == "y" or masked == "n"):
			print "\nPlease enter either 'y' or 'n'\n"
			masked = raw_input()
		if masked == "y":
			print "\nList of potential masks:\n"
			print subprocess.check_output("cd {}; ls".format(map_dir), stderr=subprocess.STDOUT, shell=True)
			print "Select a mask:\n"
			mask_name = raw_input()
			while not check_file(map_dir+'/'+mask_name):
				print "\nERROR: No such file. Try Again."
				mask_name = raw_input()
		else:
			mask_name = None

		print "\nStacking on predetermined or COOP-determined peaks?\n"
		print "\t[1] predetermined"
		print "\t[2] COOP-determined"
		print "\t[0] exit\n"
		case3 = raw_input()
		while not (case3 == "0" or case3 == "1" or case3 == "2"):
			print "\nPlease select a valid option:\n"
			case3 = raw_input()
		
		# exit
		if case3 == "0":
			sys.exit()
	
		# TODO: custom peaks
		elif case3 == "1":
			peak_def = 'customPeaks'
			print "\nCustom peaks can be supplied in the form of an array of Healpix (ring) pixels or peak positions. Please select an option:\n"
			print "\t[1] Healpix pixels"
			print "\t[2] (x,y,z) positions"
			print "\t[3] (theta, phi) positions"
			print "\t[0] exit\n"
			case4 = raw_input()
			while not (case4 == "0" or case4 == "1" or case4 == "2" or case4 == "3"):
				print "\nPlease select a valid option:\n"
				case4 = raw_input()
	
			# exit
			if case4 == "0":
				sys.exit()

			# healpix pixel peaks
			elif case4 == "1":
				read_type = hpix
			# xyz peak positions
			elif case4 == "2":
				read_type = xyz

			# thetaphi peak positions
			elif case4 == "3":
				read_type = thetaphi

			print "\nPeak files must be readable with numpy.loadtxt() or numpy.load(). Enter peaks full directory filename:"
			peak_dir = raw_input()
			while not check_file(peak_dir):
				print "\nERROR: No such file. Try Again."
				peak_dir = raw_input()
			pix = read_it(read_type, peak_dir, nside, peak_type)
			#if max_name:
			#	txt2bin_file_writer(map_dir+'/'+map_name, header, pix, 1, masked=True)
			#else:
			#	txt2bin_file_writer(map_dir+'/'+map_name, header, pix, 1)
				
			


		# COOP peak interface
		elif case3 == "2":
			peak_def = 'coopPeaks'
			print "\nSelect a COOP-defined peak:\n"
			print "\t[1] temperature maxima (hot peaks)"
			print "\t[2] temperature minima (cold peaks)"
			print "\t[3] random hot field points (hot spots)"
			print "\t[4] random cold field points (cold spots)"
			print "\t[5] random hot/cold field points (hot/cold spots)"
			print "\t[0] exit\n"
			case4 = raw_input()
			while not (case4 == "0" or case4 == "1" or case4 == "2" or case4 == "3" or case4 == "4" or case4 == "5"):
				print "\nPlease select a valid option:\n"
				case4 = raw_input()
			if case4 == "0":
				sys.exit()
			elif case4 == "1":
				hot, cold, peak = 'T', 'F', 'T'
				peak_type = 'hotPeaks'
			elif case4 == "2":
				hot, cold, peak = 'F', 'T', 'T'
				peak_type = 'coldPeaks'
			elif case4 == "3":
				hot, cold, peak = 'T', 'F', 'RANDOM'
				peak_type = 'randomHotSpots'
			elif case4 == "4":
				hot, cold, peak = 'F', 'T', 'RANDOM'
				peak_type = 'randomColdSpots'
			elif case4 == "5":
				hot, cold, peak = 'T', 'T', 'RANDOM'
				peak_type = 'randomHotColdSpots'

			print "\nThreshold value (nu) for peaks? [0 < nu <= 10 or 0 for no threshold]\n"
			nu = raw_input()
			while not (nu >= 0 or nu <= 10):
				print "\nPlease specify a value 0 < nu <= 10 or 0 for no threshold:\n"
				nu = raw_input()

			print "\nMap resolution in pixels: [10-200]\n"
			res = raw_input()
			while not (res >= 10 or nu <= 200):
				print "\nPlease specify a value 10 <= res <= 200:\n"
				res = raw_input()


			print "\nStacking options:"
			print "\tmap type:\t {}".format(map_type)
			print "\torient:\t False"
			print "\tmap dir:\t {}".format(map_dir+'/'+map_name)
			if mask_name:
				print "\tmask dir:\t {}".format(mask_name)
			else:
				print "\tmask:\t None"
			print "\tpeak type:\t {}".format(peak_type)
			if nu > 0:
				print "\tthreshold:\t {}".format(nu)
			else:
				print "\tthreshold:\t None"
			print "\tres (pix):\t {}".format(res)

			print "\nConfirm? [y/n]\n"
			confirm = raw_input()
			while not (confirm == "y" or confirm == "n"):
				print "\nPlease enter either 'y' or 'n'\n"
				confirm = raw_input()
			if confirm == "n":
				print "\nPlease start again. Exiting."
				sys.exit()

			if mask_name:
				fname = "{}_stack_{}_masked_{}_nu{}_res{}".format(map_type, map_name[:-5], peak_type, int(nu), res)
			else:
				fname = "{}_stack_{}_{}_nu{}_res{}".format(map_type, map_name[:-5], peak_type, int(nu), res)

			print "\nProposed output filename: {}".format(fname)
			print "Type 'y' to confirm or otherwise enter a custom filename:\n"
			new_fname = raw_input()
			if new_fname != 'y':
				fname = new_fname.strip()
			
			print "\nFinding peaks...\n"
			if mask_name:
				print subprocess.check_output("{}/mapio/GetPeaks -map {} -out {} -hot {} -cold {} -peak {} -mask {} -nu {}".format(COOP_dir, map_dir+'/'+map_name, stack_dir+'/out/'+fname+'.dat', hot, cold, peak, map_dir+'/'+mask_name, float(nu)), stderr=subprocess.STDOUT, shell=True)
			else:
				print subprocess.check_output("{}/mapio/GetPeaks -map {} -out {} -hot {} -cold {} -peak {} -nu {}".format(COOP_dir, map_dir+'/'+map_name, stack_dir+'/out/'+fname+'.dat', hot, cold, peak, float(nu)), stderr=subprocess.STDOUT, shell=True)
			print "\nStacking...\n"
			if mask_name:
				print subprocess.check_output("{}/mapio/Stack -map {} -out {} -peaks {} -field T -mask {} -res {} -colortable 'jet'".format(COOP_dir, map_dir+'/'+map_name, stack_dir+'/out/'+fname, stack_dir+'/out/'+fname+'.dat', map_dir+'/'+mask_name, res), stderr=subprocess.STDOUT, shell=True)
			else:
				print subprocess.check_output("{}/mapio/Stack -map {} -out {} -peaks {} -field T -res {} -colortable 'jet'".format(COOP_dir, map_dir+'/'+map_name, stack_dir+'/out/'+fname, stack_dir+'/out/'+fname+'.dat'), stderr=subprocess.STDOUT, shell=True)
			print "\nGenerating PDF...\n"
			print subprocess.check_output("{}/utils/pypl.py {}.txt {}.pdf".format(COOP_dir, stack_dir+'/out/'+fname, stack_dir+'/'+fname, res), stderr=subprocess.STDOUT, shell=True)

			print "\nPDF available at {}.pdf".format(stack_dir+'/'+fname)

			f = open(stack_dir+'/'+fname+'.log', 'w')
			f.write("map type:\t {}\n".format(map_type))
			f.write("orient: \t False")
			f.write("map dir:\t {}\n".format(map_dir+'/'+map_name))
			if mask_name:
				f.write("mask dir:\t {}\n".format(mask_name))
			else:
				f.write("mask:\t None\n")
			f.write("peak type:\t {}\n".format(peak_type))
			f.write("threshold:\t {}\n".format(nu))
			f.write("resolution:\t {}\n".format(res))
			f.close()
			print "Logfile available at {}.log".format(stack_dir+'/'+fname)

		
	# TODO: CO maps
	elif case2 == "2":
		map_type = 'CO'
		print "\nUnder construction..."

	# TODO: weak lensing maps
	elif case2 == "3":
		map_type = 'weakLensing'
		print "\nUnder construction..."

	# TODO: dust maps
	elif case2 == "4":
		map_type = 'dust'
		print "\nUnder construction..."
		
# TODO: oriented stacking
elif case1 == "2":
	stack_type = 'oriented'
	orient = r'Q_{\nabla^2T},U_{\nabla^2T}'
	print "\nSelect type of map to be stacked:\n"
	print "\t[1] thermal Sunyaev-Zeldovich (tSZ)"
	print "\t[2] CO emission"
	print "\t[3] weak lensing"
	print "\t[4] dust"
	print "\t[0] exit\n"
	case2 = raw_input()
	while not (case2 == "0" or case2 == "1" or case2 == "2" or case2 == "3" or case3 == "4"):
		print "\nPlease select a valid option:\n"
		case2 = raw_input()

	# exit
	if case2 == "0":
		sys.exit()
	
	# tsZ maps
	elif case2 == "1":
		map_type = 'tSZ'
		print "\nMaps are set to be found in {}".format(map_dir)
		print "Maps are assumed to be in Healpix format. List of potential maps to stack:\n"
		print subprocess.check_output("cd {}; ls".format(map_dir), stderr=subprocess.STDOUT, shell=True)
		print "Select a map to stack with, or change map locations by running setup.py:\n"
		map_name = raw_input()
		while not check_file(map_dir+'/'+map_name):
			print "\nERROR: No such file. Try Again.\n"
			map_name = raw_input()

		print "\nChecking if valid temperature/IQU map..."

		raw_map = fits.open('{}'.format(map_dir+'/'+map_name))
		header	= raw_map[1].header
		
		try:
			pixtype = header['PIXTYPE']
		except (KeyError):
			print "Invalid map: no PIXTYPE specified. Exiting."
			sys.exit()
		if pixtype != 'HEALPIX':
			print "Invalid map: must be Healpix format. Exiting."
			sys.exit()
		try:
			nside = header['NSIDE']
		except (KeyError):
			print "Invalid map: no NSIDE specified. Exiting."
			sys.exit()
		try:
			numcols = header['TFIELDS']
		except (KeyError):
			print "Invalid map: Unknown number of columns. Exiting."

		if numcols == 3:
			print "3 columns detected."
			I, Q, U = header['TTYPE1'].upper(), header['TTYPE2'].upper(), header['TTYPE3'].upper()
			if (I=='I' or I=='ISTOKES' or I=='INTENSITY' or I=='T' or I=='TEMPERATURE') and (Q=='Q' or Q=='QT' or Q=='QSTOKES') and (U=='U' or U=='UT' or U=='USTOKES'):
				print "IQU map detected. Proceeding...\n"
			else:
				print "WARNING: unknown IQU columns names:\n"
				print "\t"+I
				print "\t"+Q
				print "\t"+U
				print "\nProceed? [y/n]\n"
				iqu_proceed = raw_input()
				while not (iqu_proceed == "y" or iqu_proceed == "n"):
					print "\nPlease enter either 'y' or 'n'\n"
					iqu_proceed = raw_input()
				if iqu_proceed == 'n':
					print "\nValid I column names: 'I' 'ISTOKES' 'INTENSITY' 'T' 'TEMPERATURE'"
					print "Valid Q column names: 'Q' 'QT' 'QSTOKES'"
					print "Valid U column names: 'U' 'UT' 'USTOKES'"
					time.sleep(3)
					print "\nExiting."
					sys.exit()

		else:
			print "Invalid map: require 3 column IQU map for oriented stacking. Exiting."
			sys.exit()
			

		print "\nMasking the map? [y/n]\n"
		masked = raw_input()
		while not (masked == "y" or masked == "n"):
			print "\nPlease enter either 'y' or 'n'\n"
			masked = raw_input()
		if masked == "y":
			print "\nList of potential masks:\n"
			print subprocess.check_output("cd {}; ls".format(map_dir), stderr=subprocess.STDOUT, shell=True)
			print "Select a mask:\n"
			mask_name = raw_input()
			while not check_file(map_dir+'/'+mask_name):
				print "\nERROR: No such file. Try Again."
				mask_name = raw_input()
		else:
			mask_name = None

		print "\nStacking on predetermined or COOP-determined peaks?\n"
		print "\t[1] predetermined"
		print "\t[2] COOP-determined"
		print "\t[0] exit\n"
		case3 = raw_input()
		while not (case3 == "0" or case3 == "1" or case3 == "2"):
			print "\nPlease select a valid option:\n"
			case3 = raw_input()
		
		# exit
		if case3 == "0":
			sys.exit()
	
		# TODO: custom peaks
		elif case3 == "1":
			peak_def = 'customPeaks'
			print "\nCustom peaks can be supplied in the form of an array of Healpix (ring) pixels or peak positions. Please select an option:\n"
			print "\t[1] Healpix pixels"
			print "\t[2] (x,y,z) positions"
			print "\t[3] (theta, phi) positions"
			print "\t[0] exit\n"
			case4 = raw_input()
			while not (case4 == "0" or case4 == "1" or case4 == "2" or case4 == "3"):
				print "\nPlease select a valid option:\n"
				case4 = raw_input()
	
			# exit
			if case4 == "0":
				sys.exit()

			# healpix pixel peaks
			elif case4 == "1":
				read_type = hpix
			# xyz peak positions
			elif case4 == "2":
				read_type = xyz

			# thetaphi peak positions
			elif case4 == "3":
				read_type = thetaphi

			print "\nPeak files must be readable with numpy.loadtxt() or numpy.load(). Enter peaks full directory filename:"
			peak_dir = raw_input()
			while not check_file(peak_dir):
				print "\nERROR: No such file. Try Again."
				peak_dir = raw_input()
			pix = read_it(read_type, peak_dir, nside, peak_type)
			#if max_name:
			#	txt2bin_file_writer(map_dir+'/'+map_name, header, pix, 1, masked=True)
			#else:
			#	txt2bin_file_writer(map_dir+'/'+map_name, header, pix, 1)
				
			


		# COOP peak interface
		elif case3 == "2":
			peak_def = 'coopPeaks'
			print "\nSelect a COOP-defined peak:\n"
			print "\t[1] temperature maxima (hot peaks)"
			print "\t[2] temperature minima (cold peaks)"
			print "\t[3] random hot field points (hot spots)"
			print "\t[4] random cold field points (cold spots)"
			print "\t[5] random hot/cold field points (hot/cold spots)"
			print "\t[0] exit\n"
			case4 = raw_input()
			while not (case4 == "0" or case4 == "1" or case4 == "2" or case4 == "3" or case4 == "4" or case4 == "5"):
				print "\nPlease select a valid option:\n"
				case4 = raw_input()
			if case4 == "0":
				sys.exit()
			elif case4 == "1":
				hot, cold, peak = 'T', 'F', 'T'
				peak_type = 'hotPeaksOrient'
			elif case4 == "2":
				hot, cold, peak = 'F', 'T', 'T'
				peak_type = 'coldPeaksOrient'
			elif case4 == "3":
				hot, cold, peak = 'T', 'F', 'RANDOM'
				peak_type = 'randomHotSpotsOrient'
			elif case4 == "4":
				hot, cold, peak = 'F', 'T', 'RANDOM'
				peak_type = 'randomColdSpotsOrient'
			elif case4 == "5":
				hot, cold, peak = 'T', 'T', 'RANDOM'
				peak_type = 'randomHotColdSpotsOrient'

			print "\nThreshold value (nu) for peaks? [0 < nu <= 10 or 0 for no threshold]\n"
			nu = raw_input()
			while not (nu >= 0 or nu <= 10):
				print "\nPlease specify a value 0 < nu <= 10 or 0 for no threshold:\n"
				nu = raw_input()

			print "\nMap resolution in pixels: [10-200]\n"
			res = raw_input()
			while not (res >= 10 or nu <= 200):
				print "\nPlease specify a value 10 <= res <= 200:\n"
				res = raw_input()

			print "\nStacking options:"
			print "\tmap type:\t {}".format(map_type)
			print "\torient:\t True"
			print "\tmap dir:\t {}".format(map_dir+'/'+map_name)
			if mask_name:
				print "\tmask dir:\t {}".format(mask_name)
			else:
				print "\tmask:\t None"
			print "\tpeak type:\t {}".format(peak_type)
			if nu > 0:
				print "\tthreshold:\t {}".format(nu)
			else:
				print "\tthreshold:\t None"
			print "\tres (pix):\t {}".format(res)

			print "\nConfirm? [y/n]\n"
			confirm = raw_input()
			while not (confirm == "y" or confirm == "n"):
				print "\nPlease enter either 'y' or 'n'\n"
				confirm = raw_input()
			if confirm == "n":
				print "\nPlease start again. Exiting."
				sys.exit()

			if mask_name:
				fname = "{}_stack_orient_{}_masked_{}_nu{}".format(map_type, map_name[:-5], peak_type, nu)
			else:
				fname = "{}_stack_orient_{}_{}_nu{}".format(map_type, map_name[:-5], peak_type, nu)

			print "\nProposed output filename: {}".format(fname)
			print "Type 'y' to confirm or otherwise enter a custom filename:\n"
			new_fname = raw_input()
			if new_fname != 'y':
				fname = new_fname.strip()
			
			print "\nFinding peaks...\n"
			if mask_name:
				print subprocess.check_output("{}/mapio/GetPeaks -map {} -out {} -hot {} -cold {} -peak {} -mask {} -orient {} -nu {}".format(COOP_dir, map_dir+'/'+map_name, stack_dir+'/out/'+fname+'.dat', hot, cold, peak, map_dir+'/'+mask_name, orient, float(nu)), stderr=subprocess.STDOUT, shell=True)
			else:
				print subprocess.check_output("{}/mapio/GetPeaks -map {} -out {} -hot {} -cold {} -peak {} -orient {} -nu {}".format(COOP_dir, map_dir+'/'+map_name, stack_dir+'/out/'+fname+'.dat', hot, cold, peak, orient, float(nu)), stderr=subprocess.STDOUT, shell=True)
			print "\nStacking...\n"
			if mask_name:
				print subprocess.check_output("{}/mapio/Stack -map {} -out {} -peaks {} -field T -mask {} -res {} -colortable 'jet'".format(COOP_dir, map_dir+'/'+map_name, stack_dir+'/out/'+fname, stack_dir+'/out/'+fname+'.dat', map_dir+'/'+mask_name, res), stderr=subprocess.STDOUT, shell=True)
			else:
				print subprocess.check_output("{}/mapio/Stack -map {} -out {} -peaks {} -field T -res {} -colortable 'jet'".format(COOP_dir, map_dir+'/'+map_name, stack_dir+'/out/'+fname, stack_dir+'/out/'+fname+'.dat'), stderr=subprocess.STDOUT, shell=True)
			print "\nGenerating PDF...\n"
			print subprocess.check_output("{}/utils/pypl.py {}.txt {}.pdf".format(COOP_dir, stack_dir+'/out/'+fname, stack_dir+'/'+fname, res), stderr=subprocess.STDOUT, shell=True)

			print "\nPDF available at {}.pdf".format(stack_dir+'/'+fname)

			f = open(stack_dir+'/'+fname+'.log', 'w')
			f.write("map type:\t {}\n".format(map_type))
			f.write("orient: \t True")
			f.write("map dir:\t {}\n".format(map_dir+'/'+map_name))
			if mask_name:
				f.write("mask dir:\t {}\n".format(mask_name))
			else:
				f.write("mask:\t None\n")
			f.write("peak type:\t {}\n".format(peak_type))
			f.write("threshold:\t {}\n".format(nu))
			f.write("resolution:\t {}\n".format(res))
			f.close()
			print "Logfile available at {}.log".format(stack_dir+'/'+fname)

		
	# TODO: CO maps
	elif case2 == "2":
		map_type = 'CO'
		print "\nUnder construction..."

	# TODO: weak lensing maps
	elif case2 == "3":
		map_type = 'weakLensing'
		print "\nUnder construction..."

	# TODO: dust maps
	elif case2 == "4":
		map_type = 'dust'
		print "\nUnder construction..."

# TODO: oriented symmetry breaking stacking
elif case1 == "3":
	print "\nUnder construction..."
	"""
	stack_type = 'orientedSymmetryBreaking'
	orient = 'RANDOM'
	print "\nSelect type of map to be stacked:\n"
	print "\t[1] thermal Sunyaev-Zeldovich (tSZ)"
	print "\t[2] CO emission"
	print "\t[3] weak lensing"
	print "\t[4] dust"
	print "\t[0] exit\n"
	case2 = raw_input()
	while not (case2 == "0" or case2 == "1" or case2 == "2" or case2 == "3" or case3 == "4"):
		print "\nPlease select a valid option:\n"
		case2 = raw_input()
	"""
# TODO: help documentation
elif case1 == "4":
	print "\nUnder construction..."
