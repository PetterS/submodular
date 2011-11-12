import glob
import numpy as np

import glob
import numpy as np

files = glob.glob('run_*.data')

for filename in files :
	data = np.genfromtxt(filename)
	print filename
	
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
			print  '%.2f %.2f %.2f' % (np.min(data), np.median(data), np.max(data)),
		
		print '    HOCR      : labeled ',
		printinfo(lab_hocr)
		print '   rel.bound  ', 
		printinfo((optimal - hocr)/np.abs(optimal))
		print '   time  ', 
		printinfo(time_hocr)
		print ''
		
		print '    Fix et al.: labeled ',
		printinfo(lab_fixetal)
		print '   rel.bound  ', 
		printinfo((optimal - fixetal)/np.abs(optimal))
		print '   time  ', 
		printinfo(time_hocr)
		print ''
		
		print '    Optimal   : labeled ',
		printinfo(lab_optimal)
		print '   rel.bound  ', 
		printinfo((optimal - optimal)/np.abs(optimal))
		print '   time  ', 
		printinfo(time_optimal)
		print ''
		
		print '    Heur      : labeled ',
		printinfo(lab_heur)
		print '   rel.bound  ', 
		printinfo((optimal - heur)/np.abs(optimal))
		print '   time  ', 
		printinfo(time_heur)
		print ''
		
	except :
		pass
	
	print ' '
	