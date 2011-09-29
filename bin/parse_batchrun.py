import glob
import numpy as np

import glob
import numpy as np

files = glob.glob('run_*.data')

for filename in files :
	data = np.genfromtxt(filename)
	print filename
	
	try :
		hocr = data[:,7]
		lp   = data[:,9]
		heur = data[:,10]
		
		lab_hocr = data[:,3]
		lab_lp   = data[:,5]
		lab_heur = data[:,6]
		
		time_hocr = data[:,11]*1000
		time_lp   = data[:,13]*1000
		time_heur = data[:,14]*1000
		
		def printinfo(data) :
			print  '%.2f %.2f %.2f' % (np.min(data), np.median(data), np.max(data)),
		
		print '    HOCR      : labeled ',
		printinfo(lab_hocr)
		print '   rel.bound  ', 
		printinfo((lp - hocr)/np.abs(lp))
		print '   time  ', 
		printinfo(time_hocr)
		print ''
		
		print '    Optimal   : labeled ',
		printinfo(lab_lp)
		print '   rel.bound  ', 
		printinfo((lp - lp)/np.abs(lp))
		print '   time  ', 
		printinfo(time_lp)
		print ''
		
		print '    Heur      : labeled ',
		printinfo(lab_heur)
		print '   rel.bound  ', 
		printinfo((lp - heur)/np.abs(lp))
		print '   time  ', 
		printinfo(time_heur)
		print ''
		
	except :
		pass
	
	print ' '
	