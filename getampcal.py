#!/usr/bin/env python

#----------------------------------------------------------------------------------
# This script gets the amplitude corrections applied to a priori calibrated and fringe-fitted 
# data by self-calibration for all experiments observed at a given year and month and 
# for the specified station.
#
# According tot the description from JIVE:
# Amplitude corrections applied to a priori calibrated and fringe-fitted data 
# by self-calibration:
# These are the amplitude corrections derived from amplitude-and-phase self-calibration 
# of the sources (done after an initial phase-only self-calibration). Corrections should 
# be small if the a priori amplitude calibration was good. Please note that the error 
# given here is the square root of the error in the original Tsys value.
#
# An example of the pipeline results is:
# http://archive.jive.nl/exp/N15M2_151021/pipe/n15m2.html
#
# Explanations of the content here:
# http://www.evlbi.org/pipeline/pipe_desc.html
#
# Usage example: Get data from YS station in October 2016
# get_ampcal.py 2016 10 YS
#
#---------------------------------------------------------------------------------
# Pablo de Vicente  (p.devicente@oan.es) 2009

import sys
from ftplib import FTP
import urllib2

debug = True

def main():

	# Go to http://archive.jive.nl/exp/ and get the list of directories
	# and store it in projectList

	try:
		url = 'http://archive.jive.nl/exp/'
	        r = urllib2.urlopen(url)
		fileContent = r.read()
		fileList = fileContent.split('\n')
		projectList = []
		for line in fileList:
			if 'DIR' in line:
				if '_' in line:
					projectList.append(line.split('>')[2].split('/')[0])

		#if debug:
		#	print projectList
	except Exception, e:
		print "Caught an exception. Connection to JIVE may have failed"
		print str(e)
		sys.exit(-1)

	# Ask the user for the input parameters
        if len(sys.argv) > 3:

                year = sys.argv[1]
                month = sys.argv[2]
		station = sys.argv[3]
		station = station.upper()
		if len(year) != 4 or len(month) != 2:
			print "Year is a 4 digit number and month a two digit number"
			sys.exit(-1)
	else:
		print "Usage is: %s year month stationCode" % (sys.argv[0])
		sys.exit(-1)

	dateObs = "_%s%s" %(year[-2:], month)
	dateObsLong = "%s%s" %(year, month)
	outputFileName = 'ampcal_' + dateObsLong + '_' + station.lower() + '.dat' 

	ampcalList = []

	# Go through all experiments for the specified month and year and make a list
	# with the name of the project plus the name of the .ampcal file in that project
	# There may be more than one ampcal per project. Indeed there is one ampcal per
	# source

	print "Getting all experiments names for the specified date from http://archive.jive.nl/exp ..."
	print " "
	for project in projectList:
		if dateObs in project:
			url = "http://archive.jive.nl/exp/%s/pipe/" % (project)
			dirpipe = urllib2.urlopen(url)
			dirpipeContent = dirpipe.read()
			dirpipeList = dirpipeContent.split('\'n')
			for line in dirpipeList:
				line2List = line.split('\n')
				for line2 in line2List: 
					if '.ampcal' in line2:
						ampcalfile = line2.split('>')[2].split('/')[0]
						ampcalfile = ampcalfile.split('<')[0]
						ampcalList.append([project, ampcalfile])

	# For each .ampcal file in the list get the amplitude calibration information 
	# for the specified station. Only the mean calibration is obtained. Individual
	# channels are dropped.
	# A list is constructed with station, project name, frequency, observed source, and results

	resultsList = []

	print "Retrieving and parsing ampcal files from http://archive.jive.nl/exp/*/pipe ..."
	print " "
	if debug:
		print ampcalList
	for ampcalfile in ampcalList:
		fileWeb = "http://archive.jive.nl/exp/%s/pipe/%s" % (ampcalfile[0],ampcalfile[1])
		r = urllib2.urlopen(fileWeb)
		fileContent = r.read()
		fileList = fileContent.split('\n')
		foundStation = False
		foundFrequency = False
		observedSource = ampcalfile[1].split('_')[-3]
		for lines in fileList:
			try:
				allAmpCalList = lines.split()
				if station in lines:
					foundStation = True
				if foundFrequency:
					foundFrequency = False
					freqGHz = float(lines.split()[2])
				if 'FQID' in lines:
					foundFrequency = True
				if foundStation and 'All:' in lines:
					auxList = [freqGHz, station, ampcalfile[0], observedSource]
					auxList.extend(allAmpCalList)
					resultsList.append(auxList)
					foundStation = False
					if debug:
						print "(%s) \t %s \t %s \t %s \t %s" %(station, ampcalfile[0], freqGHz, observedSource, lines)
			except Exception, ex:
				print "Error in the loop: %s" % (str(ex))

	# Open the file
	outputFile = open(outputFileName, 'w')

	# Sort the results by frequency and print them
	print "Sorting the results by frequency ..."
	print " "
	resultsList.sort()

	firstElement = resultsList[0]
	freqList = [[]]
	freqPrevious = firstElement[0]
	i = 0
	for res in resultsList:
		if res[0] == freqPrevious:
			freqPrevious = res[0]
			freqList[i].append([res[6], res[0], res[1], res[2], res[3], res[4], res[5], res[7], res[8], res[9]])
		else:
			freqList.append([[res[6], res[0], res[1], res[2], res[3], res[4], res[5], res[7], res[8], res[9]]])
			i = i + 1
			freqPrevious = res[0]

	print "Found %d different frequency setups" %(i+1)

	for felist in freqList:
	
		avg_res6 = 0
		avg_res0 = 0
		counter = 0

		felist.sort()
		print "----------------"
		for res in felist:
			strLine = "%3.5f \t (%s) \t %s \t %12s %s %s %s[0;35m %s%s[0m %s %s %s" %(res[1], res[2], res[3], res[4], res[5], res[6], chr(27), res[0],chr(27), res[7], res[8], res[9])
			strLineFile = "%3.5f \t (%s) \t %s \t %12s %s %s %s %s %s %s" %(res[1], res[2], res[3], res[4], res[5], res[6], res[0], res[7], res[8], res[9])
			print strLine
			strLineFile = strLineFile + '\n'
			outputFile.write(strLineFile)

			# get average per frequency
			# Amplituce correction:
			avg_res6 = avg_res6 + float(res[6])
			# Amplitude median error
			avg_res0 = avg_res0 + float(res[0])
			counter = counter + 1

		avg_res6 = avg_res6 / counter
		avg_res0 = avg_res0 / counter
		strLine = "\t \t \t \t \t \t \t   %3.4f  %3.5f" %(avg_res6, avg_res0)
		print strLine
		strLineFile = strLine + '\n'
		outputFile.write(strLineFile)


	# Old way without ordering by median error
	'''
	for res in resultsList:
		#strLine = "%3.5f \t (%s) \t %s \t %12s %s" %(res[0], res[1], res[2], res[3], res[4])
		strLine = "%3.5f \t (%s) \t %s \t %12s %s %s %s[0;35m %s%s[0m %s %s %s" %(res[0], res[1], res[2], res[3], res[4], res[5], chr(27), res[6],chr(27), res[7], res[8], res[9])
		strLineFile = "%3.5f \t (%s) \t %s \t %12s %s %s %s %s %s %s" %(res[0], res[1], res[2], res[3], res[4], res[5], res[6], res[7], res[8], res[9])
		print strLine
		strLineFile = strLineFile + '\n'
		outputFile.write(strLineFile)
	'''

	outputFile.close()

	print " "
	print "Results are stored at %s for plotting !!" % (outputFileName)
	print 'You can use gnuplot for plotting. Do inside gnuplot: plot "%s" using 6' % (outputFileName)
	print " "

			
if __name__ == "__main__":
	main()

