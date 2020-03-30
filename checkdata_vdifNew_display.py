#!/usr/bin/env python

#	FJB
import sys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy
import math
if sys.version_info[0] < 3:
   import Tkinter as Tk
else:
   import tkinter as Tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import os
import signal
from time import sleep, ctime


def handler(signum, frame):
   print '\nSignal handler called with signal', signum
   print "CTRL-C pressed..."
   aux.shutdown()

class specDisplay:

   def __init__(self):

	self.root = Tk.Tk()
        self.root.wm_title("SPEC DISPLAY")
        self.root.geometry('800x600')
        self.root.protocol("WM_DELETE_WINDOW", self.shutdown)

	self.fig = plt.figure(1)
	plt.tick_params(axis='x', labelbottom='off')
        plt.tick_params(axis='y', labelleft='off')

	self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
	self.canvas.show()
	self.canvas.get_tk_widget().pack(fill=Tk.BOTH, expand=1)

	toolbar = NavigationToolbar2TkAgg(self.canvas, self.root)
	toolbar.update()

	# Plot the results with matplotlib. Show them for some seconds and close the plot.
        self.fullspecfile = "/usr2/log/spec.out"
	self.fileTime = ""
	self.root.after(100,self.printCanvas)

   def startMainloop(self):
	self.root.mainloop()

   def printCanvas(self):

	#Get file's last modified date
	fileTime = os.path.getmtime(self.fullspecfile)

	if fileTime != self.fileTime:

		#Read file
	        data = numpy.loadtxt(self.fullspecfile)
	
	        # First column contains frequency in MHz.
		if len(data) == 0:
			self.root.after(100,self.printCanvas)
			return
	        ncols = len(data[0])-1
	        ncols = ncols / 2
	        #if ncols != bbcs:
	                #print "Columns (%d) != Channels (%d)" % (ncols, bbcs)
	
	        # First column is the X with already the correct frequency
	        x = data[:,0]           
	
	        if ncols > 4:
	                plines =int(math.ceil(ncols/4.0))
	                pcols = 4
	        elif ncols == 4:
	                plines = 2
	                pcols = 2
	        else:
	                plines = 1
	                pcols = ncols
	                
		#Clear current figure
		self.fig.clear()

		plt.suptitle("Output Spectrum: %s" % ctime(fileTime))

	        # Now plot
	        i = 0
	        j = 0
	        while i < ncols:
	                j = j + 1
	                for k in range(1, pcols+1):
				if i >= ncols:
					break
	                        #print "i: ", i, "j, k", j, k
				axarr = self.fig.add_subplot(plines, pcols, i+1)
	        		# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
				if j < plines:
					plt.setp(axarr.get_xticklabels(), visible=False)
				#if k == pcols:
					#plt.setp(axarr.get_yticklabels(), visible=False)
				axarr.plot(x, data[:,1+i])
				axarr.grid()
	                        strline = "Channel %d" % (1+i)
				axarr.set_title(strline)
	                        if j == plines:
					axarr.set_xlabel("BW MHz")
	                        i = i + 1

		#Draw new figure	
		self.canvas.draw()

	#Save file's last modified time
	self.fileTime = fileTime

	#Check file to draw a new figure every 100 microseconds
	self.root.after(100,self.printCanvas)

   def shutdown(self):
	print "Disconnecting ..."
        self.root.destroy()
        sys.exit(0)

if __name__ == "__main__":
	signal.signal(signal.SIGTERM, handler)
	signal.signal(signal.SIGINT, handler)
	aux = specDisplay()
	aux.startMainloop()
