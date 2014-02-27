#!/usr/bin/python3

import glob
import numpy as np

import glob
import numpy as np

files = glob.glob('run_*.data')

for filename in files:
	data = np.genfromtxt(filename)
	print(filename)
	
	try :
		hocr    = data[:,7]
		fixetal = data[:,17]
		optimal = data[:,9]
		heur    = data[:,10]
		
		lab_hocr = data[:,3]
		lab_fixetal = data[:,15]
		lab_optimal = data[:,5]
		lab_heur = data[:,6]
		
		time_hocr    = data[:,11]*1000
		time_optimal = data[:,13]*1000
		time_heur    = data[:,14]*1000
		
		def printinfo(data) :
			print('%.2f %.2f %.2f' % (np.min(data), np.median(data), np.max(data)), end="")
		
		print('    HOCR      : labeled ', end="")
		printinfo(lab_hocr)
		print('   rel.bound  ', end="") 
		printinfo((optimal - hocr)/np.abs(optimal))
		print('   time  ', end="")
 
		printinfo(time_hocr)
		print('')
		
		print('    Fix et al.: labeled ', end="")
		printinfo(lab_fixetal)
		print('   rel.bound  ', end="")
		printinfo((optimal - fixetal)/np.abs(optimal))
		print('   time  ', end="") 
		printinfo(time_hocr)
		print('')
		
		print('    Optimal   : labeled ', end="")
		printinfo(lab_optimal)
		print('   rel.bound  ', end="") 
		printinfo((optimal - optimal)/np.abs(optimal))
		print('   time  ', end="") 
		printinfo(time_optimal)
		print('')
		
		print('    Heur      : labeled ', end="")
		printinfo(lab_heur)
		print('   rel.bound  ', end="") 
		printinfo((optimal - heur)/np.abs(optimal))
		print('   time  ', end="") 
		printinfo(time_heur)
		print('')
		
	except :
		pass
	
	print(' ')
	