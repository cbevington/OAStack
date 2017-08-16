import time
import subprocess
import numpy as np
from stacking_functions import check_dir

wd = subprocess.check_output('pwd', shell=True)[:-1]

print "\nBeginning OAStack setup. Make sure you are running this in the OAStack directory."
time.sleep(3)
print "\nBy default raw maps are assumed to be located in OAStack/maps and outputted stacks in OAStack/stacks. Do you want to use the default settings? [y/n]\n"
default = raw_input()

while not (default == "y" or default == "n"):
	print "\nPlease enter either 'y' or 'n'"
	default = raw_input('\n')

COOP_dir = raw_input('\nEnter location of COOP:\n')
while not check_dir(COOP_dir):
	print "\nERROR: No such directory. Try Again."
	COOP_dir = raw_input('\n')

if default == "y":
	if not check_dir(COOP_dir+'/maps'):
		_ = subprocess.Popen("mkdir {}/maps".format(wd), shell=True),
	if not check_dir(COOP_dir+'/stacks'):
		_ = subprocess.Popen("mkdir {}/stacks".format(wd), shell=True),
		_ = subprocess.Popen("mkdir {}/stacks/out".format(wd), shell=True),
	map_dir, stack_dir = wd+"/maps", wd+"/stacks"
else:
	map_dir	= raw_input('\nEnter location where raw maps will be located:\n')
	while not check_dir(map_dir):
		print "\nERROR: No such directory. Try Again."
		map_dir = raw_input('\n')

	stack_dir = raw_input('\nEnter location for outputted stacks:\n')
	while not check_dir(stack_dir):
		print "\nERROR: No such directory. Try Again."
		stack_dir = raw_input('\n')

print "\nRaw maps to be stored at {}".format(map_dir)
print "Stacks to be stored at {}".format(stack_dir)
time.sleep(1)

print "\nChecking for required libraries..."

try:
	import healpy
	print "Healpy installed."
except (ImportError):
	print "Healpy not installed. Beginning installation...\n"
	subprocess.check_output("pip install healpy", stderr=subprocess.STDOUT, shell=True)

try:
	import astropy
	print "Astropy installed."
except (ImportError):
	print "Astropy not installed. Beginning installation...\n"
	subprocess.check_output("pip install astropy", stderr=subprocess.STDOUT, shell=True)

ftcomp = raw_input('\nEnter name of Fortran compiler used to build COOP: [gfortran/ifort]\n')
while not (ftcomp == 'gfortran' or ftcomp == 'ifort'):
	print "\nOnly gfortran or ifort supported. Try Again."
	ftcomp = raw_input('\n')

settings = np.array([['COOP', COOP_dir], ['maps', map_dir], ['stacks', stack_dir], ['fortran', ftcomp]])
np.savetxt('{}/settings.txt'.format(wd), settings, delimiter='\t', fmt='%s')

print "\nSetup Complete! Settings have been stored in {}/settings.txt and can be edited at any time directly or by running this script again.".format(wd)


