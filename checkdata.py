#!/usr/bin/env python
#
# PdV 2017
# Fran Beltran 2017. Moved the display part to checkdata_vdifNew_display.py

import sys
import os,math
import string
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pylab
import numpy
import time

# Copy the extracted data from the recorder to the FS and use mk5access tools to plot whatever is important. To be used with checkdata_vdifNew_display.py

## Tuning here ##_________________________

debug = False
sleepsecs = 20				# Windows is shown for some seconds
recdirectory = '/home/oper/data'	# Directory where the data are stored at the Flexbuff, Mark5C or Mark6
recfile = 'systest.vdif'		# File name
fsdirectory = '/usr2/log'		# Directory where autocorrelation spectrum will be copied
specfile = 'spec.out'			# File name with autocorrelation spectrum
pylab.rcParams["figure.figsize"] = [12,8]	# Sise of the window

##########################################
def getRecorderIP():
	'''Get IP address from mk5ad.ctl
	'''

	fileRec = "/usr2/control/mk5ad.ctl"	
	fInr = open(fileRec, 'r')
	content = fInr.readlines()
	fInr.close()
	for line in content:
		if line[0] == '*':
			pass
		else:
			ipaddress = line.split()[0]
	return ipaddress

##########################################
def getVDIFformat(logfile):
	'''Get VDIF format, number of channels, BW per channel
	'''

	try: 
		fIn = open(logfile)
		fileContent = fIn.readlines()
		fIn.close()
	except:
		print "File not found"
		sys.exit()
  
	#start=0
	#if len(fileContent)>4000:
	#	start=len(fileContent)-4000

	for line in fileContent:
		# Find mode 
		# jive5ab/!mode? 0 : VDIF_8000-128-16-2 : VDIF : 32 : 4000000.000 : 8000 ;
		if "/jive5ab/!mode?" in line:
			newline=line.split(":")
			vdifformat = newline[3].strip()

	# VDIF_8000-128-16-2 => 128 Mbs, 16 channels, 2 bits => 128 / (2 bands * 2 bits * channels)
	vdiflist = vdifformat.split('-')
	aux = int(vdiflist[0].split('_')[1])
	mbps = int(vdiflist[1])
	bbcs = int(vdiflist[2])
	bits = int(vdiflist[3])
	sampling = 2			# Twice the BW
	BW = mbps/(bbcs*bits*sampling)

	if debug:
		print vdifformat, mbps,bbcs, bits, BW

	return vdifformat, mbps, bbcs, BW
	
##########################################
def main(args):

	if len(args) < 2:
		print "\n Usage: %s lognm vdiffile\n" % args[0]
		sys.exit()

	PrcDir = "/usr2/proc/"
	LogDir = "/usr2/log/"

	# Recorder: flexbuff or Mark5C or Mark6. Get it from mk5ad.ctl
	recIP = getRecorderIP()

	#prcfile = PrcDir+args[1]+'.prc'
	logfile = LogDir+args[1]+'.log'

	fileVDIF=args[2]

	vdifformat, mbps, bbcs, BW = getVDIFformat(logfile)

	# get the file from the recorder ...
	#linecmd = "scp -p oper@%s:%s/%s %s/"  % (recIP, recdirectory, recfile, fsdirectory)
	linecmd = "scp -p oper@%s:%s %s/%s"  % (recIP, fileVDIF, fsdirectory, recfile)
	os.system(linecmd)

	linecmd = "m5bstate %s/%s %s 200 | tee -a %s.m5bstat" % (fsdirectory, recfile, vdifformat, args[1])
	os.system(linecmd)

	# Get the spectrum with DifX
	linecmd = "m5spec %s/%s %s 4000 1024 '%s/%s' -dbbc > m5spec.tmp" % (fsdirectory, recfile, vdifformat, fsdirectory, specfile)
	#linecmd = "m5spec %s/%s %s 32000 1024 '%s/%s' -dbbc > m5spec.tmp" % (fsdirectory, recfile, vdifformat, fsdirectory, specfile)
	#print linecmd
	os.system(linecmd)

if __name__ == "__main__":
        main(sys.argv)
