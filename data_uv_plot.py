import time,aipy,astropy.time as ap,sys,numpy as np,matplotlib,matplotlib.pyplot as plt

def plot_all(f_min,f_max,data):
	plt.ion()	
	fig=plt.figure()
	freq=np.linspace(f_min,f_max,512)

	ax3=fig.add_subplot(1,1,1)	
	line1,=ax3.plot(freq,10*np.log10(np.abs(data['00'][-1])),'-b',label='0x-0x')
	line2,=ax3.plot(freq,10*np.log10(np.abs(data['11'][-1])),'-r',label='1x-1x')
	line3,=ax3.plot(freq,10*np.log10(np.abs(data['01'][-1])),'-g',label='0x-1x')
	ax3.legend(loc='upper right')
	plt.xlabel('Frequency(MHz)')
	plt.ylabel('Power level(dBm)')
		
	#print np.shape(data['00'])
	for i in range(len(data['00'])):	
		
		line1.set_ydata(10*np.log10(np.abs(data['00'][i])))		
		line2.set_ydata(10*np.log10(np.abs(data['11'][i])))
		line3.set_ydata(10*np.log10(np.abs(data['01'][i])))
		plt.title('Correlation sample # %d' %i)		
		plt.pause(0.1)
		fig.canvas.draw()
		
	plt.ioff()
def get_seconds(t1,t2):
	t_1=ap.Time(t1,format='jd')
	t_2=ap.Time(t2,format='jd')
	return (t_2.unix-t_1.unix)
		

    
		
if __name__ == '__main__':
	uv=aipy.miriad.UV(sys.argv[1])
	times=[]
	data={}
	data['00']=[]
	data['01']=[]
	data['11']=[]
	freq_MHz=[]
	#for (uvw,t,(i,j)),d,f in uv.all(raw=True):
	for (uvw,t,(i,j)),d in uv.all():
		if i==0 and j==0:
			data['00']+=[d]
			times+=[t]
			
		if i==0 and j==1:	
			data['01']+=[d]
		if i==1 and j==1:
			data['11']+=[d]
		
		freq_MHz=uvw
	#print len(times)
	#print len(data['00'][0])
	

	f_min=freq_MHz[0]
	f_max=freq_MHz[-1]
	t_min=times[0]
	t_max=times[-1]
	

	
	
	sec=get_seconds(t_min,t_max)
	
	plot_all(f_min,f_max,data)
	plt.figure()
	#plt.ion()
	plt.imshow(10*np.log10(np.abs(data['01'])),cmap='coolwarm',aspect='auto',extent=[f_min,f_max,sec,0],interpolation='none')
	plt.colorbar()
			
	plt.xlabel('Frequency(MHz)')
	plt.ylabel('Time(seconds)')
	plt.title('Waterfall:Cross-Correlation')
	plt.show()
	
	
	#plt.pause(30)
	#plt.ioff()
