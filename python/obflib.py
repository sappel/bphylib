import numpy
   
def smooth(x,window_len=10,window='hanning'):
	
	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."
 
	if x.size < window_len:
		raise ValueError, "Input vector needs to be bigger than window size."
  
	if window_len<3:
		return x
 
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
  
	s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
     #print(len(s))
	if window == 'flat': #moving average
		w=ones(window_len,'d')
	else:
		w=eval('numpy.'+window+'(window_len)')
 
	y=numpy.convolve(w/w.sum(),s,mode='same')
	return y[window_len-1:-window_len+1]
