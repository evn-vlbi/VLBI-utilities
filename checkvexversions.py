#!/usr/bin/env python

import sys
from ftplib import FTP

schedsDirectory = "/usr2/oper/scheds/"
stationCode = "Ys"


def main():

	global serverFile, localFile, localVersion, experimentArray, stationPresent, stationCode

	if len(sys.argv) == 3:
                monthyear = sys.argv[1]
		stationCode = sys.argv[2]
	else:
		print "Usage: %s monthyear stationCode   (for example %s oct15 Ys)" % (sys.argv[0], sys.argv[0])
		sys.exit()

	ftp = None
	#try block?:
	ftp = FTP('vlbeer.ira.inaf.it')
	ftp.login('evn','morevlbeer!')
	directory = 'vlbi_arch/%s' % monthyear
	ftp.cwd(directory)
	listOut = ftp.nlst("*.vex")

	experimentArray = []
	experimentLatestArray = []

	print " Checking directory %s ..." % monthyear

	for serverFile in listOut:
		lastestServerVersion = '?'
		stationPresent = False
		localVersion = extractLocalFileVersion(schedsDirectory, serverFile)
		getFile = "RETR %s" % serverFile
		ftp.retrlines(getFile, extractvexFileVersion)
		if stationPresent:
			experimentArray.append([serverFile, serverVersion, localVersion]) 
	

	#printInfoVersions(experimentArray)

	print " Checking %s/.latest ..." % monthyear
	listOutLatest = ftp.nlst(".latest/*.vex")
	for serverFile in listOutLatest:
		stationPresent = False
		localList = serverFile.split('/')
		localFile = localList[-1] 
		localVersion = extractLocalFileVersion(schedsDirectory, localFile)
		getFile = "RETR %s" % serverFile
		ftp.retrlines(getFile, extractvexFileVersion)
		if stationPresent:
			experimentLatestArray.append([localFile, serverVersion, localVersion]) 
	ftp.close()

	finalArray = mixArrays(experimentArray, experimentLatestArray)

	printInfoVersions(finalArray, monthyear)

def getPIVersion(line):
	'''get the PI version from the line. We get the thrid element from the back. We assume this line already has a string like
	PI revision version
	@param line
	@return PI version as a string
	'''
	listaLine = line.split(' ')
	return listaLine[-3]

def addLatestServerVersion(experimentArray, serverFile, serverVersion):
	'''Being built ....
	'''
	for element in experimentArray:
		if element[0] == serverFile:
			pass
			 
	
def printInfoVersions(finalArray, monthyear):
	'''Prints information with the name of the experiment, the version in the server and the local version 
	@param experimentArray  Array with information on VEX file name, server version and localversion. Like this:
		['f15m1', '2.0000', '2.000']
	'''
	color = "\033[91m"
	resetColor = "\033[0m"

	strLine = " Experiment \t %s / .latest / local \t Update VEX?" % (monthyear)
	print strLine
	strLine = " --------------------------------------------------------- "
	print strLine

	for element in finalArray:
		paddedName = element[0].rjust(10,' ')
		if element[3] == -1:
			element[3] = '?     '

		if element[1] == element[2] and element[1] == element[3]:
			updateFile = "No"
			strLine = " %s \t %s / %s / %s \t %s" % (paddedName, element[1], element[2], element[3], updateFile)
		else:
			updateFile = "Yes"
			strLine = " %s \t %s / %s / %s \t %s%s%s%s%s" % (paddedName, element[1], element[2], element[3], chr(27), color, updateFile, chr(27), resetColor)

		print strLine
		

def mixArrays(experimentArray, experimentLatestArray):
	'''
	'''

	finalArray = []

	i = 0
	# Loop trough the difrectoy and add information from .latest. If the experiment
	# is not in .latest, add a ?

	for experiment in experimentArray:
		finalArray.append([experiment[0], experiment[1], '?     ', experiment[2]])
		for latestexperiment in experimentLatestArray:
			if latestexperiment[0] == experiment[0]:
				finalArray[i][2] = latestexperiment[1]
		i = i + 1

	auxArray = []
	for felement in finalArray:
		auxArray.append(felement)

	latestOnly = False
	for latestexperiment in experimentLatestArray:
		'''
		for experiment in experimentArray:
			if latestexperiment[0] == experiment[0]:
				latestOnly = True
		'''
		for experiment in auxArray:
			if latestexperiment[0] == experiment[0]:
				latestOnly = True
		if not latestOnly:
			finalArray.append([latestexperiment[0], '?', latestexperiment[1], latestexperiment[2]])
		latestOnly = False

	return finalArray


def extractvexFileVersion(line):
	'''This is a callback, it is called with every line read from the FTP
	It requires global variables to get and return values, since the inpout parameter is
	always line.
	Goal: get the PI version from the file in the FTP server and set a variable to True if it finds
	the station code in the VEX file

	@param line Read from the file in the fTP

	@globalparam, stationCode
	@globalreturn stationPresent: True or False (default vaule)
	@globalreturn serverVersion: PI version from the VEX file in the FTP server
	'''
	global serverVersion, stationCode, stationPresent


	if "PI revision number:" in line:
		serverVersion = getPIVersion(line)
	stLine = "def %s;" % stationCode
	if stLine in line:
		stationPresent = True


def extractLocalFileVersion(schedsDirectory, fileName):
	'''Get the PI version from the vex file in our local host. If the file does not exist
	return -1
	@param schedsDirectory. Direcotry where we look for the VEX file. Could be read from /usr2/control/sked.ctl
	@fileName Name of the VEX file
	'''

	fileIn = "%s%s" % (schedsDirectory, fileName)
	try:
		fIn = open(fileIn)
		content = fIn.readlines()
		fIn.close()
	except Exception, ex:
		return -1

	for line in content:
		if "PI revision number:" in line:
			return getPIVersion(line)
	return -1



if __name__ == "__main__":
        main()
