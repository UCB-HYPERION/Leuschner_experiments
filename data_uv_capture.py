import casperfpga,casperfpga.snapadc,aipy,astropy.time, numpy as np, time,struct,sys,logging,matplotlib.pyplot as plt

katcp_port=7147



def get_data(fpga):
     #get the data...    
 d={}
 d['xx']=np.fromstring(fpga.read('xx', 8*512), dtype='>i8')
 d['yy']=np.fromstring(fpga.read('yy', 8*512), dtype='>i8')
 d['xy']=np.fromstring(fpga.read('xy_r',8*512), dtype = '>i8')+1j*np.fromstring(fpga.read('xy_i',8*512), dtype = '>i8')
 
 return d

def return_jd(t):
 t_u = astropy.time.Time(t, format= 'unix')
 return t_u.jd

def write_file(c,n,times,dat,swp):
 fmt="%m-%d-%Y_%H-%M-%S"
 t0=times[0]
 t0_str = time.strftime(fmt,time.localtime(t0))
 print 'writing obs_dat %s.uv....' % t0_str
 sys.stdout.flush()
 uv = aipy.miriad.UV('obs_dat_%s.uv' % t0_str,'new')
 uv['history'] = 'this observation started at %s ' % t0_str

 pols = [-5]
 uv._wrhd('obstype','auto-cross')
 #uv._wrhd('history','MDLVIS: created file.\nMDLVIS: ' + ' '.join(sys.argv) + '\n')
 uv.add_var('telescop','a'); uv['telescop'] = 'AIPY'
 uv.add_var('operator','a'); uv['operator'] = 'AIPY'
 uv.add_var('version' ,'a'); uv['version'] = '0.0.1'
 uv.add_var('epoch'   ,'r'); uv['epoch'] = 2000.
 uv.add_var('source'  ,'a'); uv['source'] = 'zenith'
 uv.add_var('latitud' ,'d'); uv['latitud'] = 0.0
 uv.add_var('dec'     ,'d'); uv['dec'] = 0.0
 uv.add_var('obsdec'  ,'d'); uv['obsdec'] = 0.0
 uv.add_var('longitu' ,'d'); uv['longitu'] = 0.0
 uv.add_var('npol'    ,'i'); uv['npol'] = len(pols)
 uv.add_var('nspect'  ,'i'); uv['nspect'] = 1
 uv.add_var('nants'   ,'i'); uv['nants'] = 2
 uv.add_var('antpos'  ,'d')
 antpos = np.array(([0,0,0],[0,0,1]), dtype=np.double)
 uv['antpos'] = antpos.transpose().flatten()
 uv.add_var('sfreq'   ,'d'); uv['sfreq'] = 0.08
 uv.add_var('freq'    ,'d'); uv['freq'] = 0.16
 uv.add_var('restfreq','d'); uv['restfreq'] = 0.
 uv.add_var('sdf'     ,'d'); uv['sdf'] = .00015625
 uv.add_var('nchan'   ,'i'); uv['nchan'] = 512
 uv.add_var('nschan'  ,'i'); uv['nschan'] = 0
 uv.add_var('inttime' ,'r'); uv['inttime'] = 0.
 # These variables just set to dummy values
 uv.add_var('vsource' ,'r'); uv['vsource'] = 0.
 uv.add_var('ischan'  ,'i'); uv['ischan'] = 1
 uv.add_var('tscale'  ,'r'); uv['tscale'] = 0.
 uv.add_var('veldop'  ,'r'); uv['veldop'] = 0.
 # These variables will get updated every spectrum
 #uv.add_var('coord'   ,'d')
 uv.add_var('time'    ,'d')
 uv.add_var('lst'     ,'d')
 uv.add_var('ra'      ,'d')
 uv.add_var('obsra'   ,'d')
 #uv.add_var('baseline','r')
 uv.add_var('pol'     ,'i')

 uv.add_var('switch','i');

 c2=float(c)/1024
 uvw=np.array([(c*0.5*(n-1)),c2,(c*0.5*n)],dtype=np.double)

 for i in range(len(times)):
  uv['pol']=pols[0]
  uv['lst']=0.
  uv['ra']=0.
  uv['obsra']=0.

  for key in dat[i].keys():
    if key == 'xx':
     preamble = (uvw,return_jd(times[i]),(0,0))
     if (n%2)==0:
      dat[i]['xx']=dat[i]['xx'][::-1]
     data=np.ma.array(dat[i]['xx'],mask=np.array(5*[1]+(len(dat[i]['xx'])-10)*[0]+5*[1]))
     uv.write(preamble,data)
    elif key == 'xy':
     preamble = (uvw,return_jd(times[i]),(0,1))
     if (n%2)==0:
      dat[i]['xy']=dat[i]['xy'][::-1]
     #data=np.ma.array(dat[i]['xy'],mask=np.zeros(len(dat[i]['xy'])))
     data=np.ma.array(dat[i]['xy'],mask=np.array(5*[1]+(len(dat[i]['xy'])-10)*[0]+5*[1]))
     uv.write(preamble,data)
    else:
     preamble = (uvw,return_jd(times[i]),(1,1))
     if (n%2)==0:
      dat[i]['yy']=dat[i]['yy'][::-1]
     #data=np.ma.array(dat[i]['yy'],mask=np.zeros(len(dat[i]['yy'])))
     data=np.ma.array(dat[i]['yy'],mask=np.array(5*[1]+(len(dat[i]['yy'])-10)*[0]+5*[1]))
     uv.write(preamble,data)
 del(uv)

 time.sleep(0.5)
 t1=times[-1]
 t_d=(t1-t0)
 print 'done in %.2f seconds' % t_d
 


 

#START OF MAIN:
if __name__ == '__main__':
 
 from optparse import OptionParser
 p = OptionParser()
 p.set_usage('snap_init.py <ROACH_HOSTNAME_or_IP> [options]')
 p.set_description(__doc__)
 p.add_option('-l', '--acc_len', dest='acc_len', type='int',default=(2**28)//512,
    help='Set the number of vectors to accumulate between dumps. default is 2*(2^27)/512, or just under 4 seconds.')
 p.add_option('-s', '--skip', dest='skip', action='store_true',
    help='Skip reprogramming the FPGA and configuring EQ.')
 p.add_option('-b', '--fpg', dest='fpgfile',type='str', default='dual_input_test_2.fpg',
    help='Specify the fpg file to load. Default: dual_input_test_2.fpg')
 p.add_option('-f', '--fftshift', dest='fftshift',type='int', default=0xffff,
    help='FFT shift schedule as an integer. Default:0xffff')
 p.add_option('-t', '--filetime', dest='filetime', type=int, default=300,
    help='Time in seconds of each data file. Default:300')
 p.add_option('-n', '--nnumber', dest='nnumber', type=int, default=1,
    help='enter 1 for first nyquist zone,2 for second etc...')
 p.add_option('-p', '--port_numbers', dest='inputs', type='str', default='1,9',
    help='enter the input ports on the SNAP. Default: 1,9 ')
 
 opts, args = p.parse_args(sys.argv[1:])

 if args==[]:
  print 'Please specify a SNAP board. Run with the -h flag to see all options.\nExiting.'        
  exit()
 else:
  snap = args[0]
 
 if opts.fpgfile != '':
  bitstream = opts.fpgfile

 inp = opts.inputs.split(",")
 for i,v in enumerate(inp):
  inp[i]=int(inp[i])-1



 print('Connecting to server %s on port %i... '%(snap,katcp_port))
 fpga = casperfpga.CasperFpga(sys.argv[1],transport=casperfpga.transport_katcp.KatcpTransport)
 time.sleep(0.2)
 
 if fpga.is_connected():
  print 'ok\n'
 else:
  print 'ERROR connecting to server %s on port %i.\n'%(snap,katcp_port)
  exit_fail()

 print '------------------------'
 print 'Programming FPGA with %s...' %bitstream,
 sys.stdout.flush()
 
 if not opts.skip:
  fpga.upload_to_ram_and_program(bitstream)
  print 'done'
 else:
  print 'Skipped.'


 fpga.write_int('adc_0_select',inp[0])
 fpga.write_int('adc_1_select',inp[1])

 adc = casperfpga.snapadc.SNAPADC(fpga, ref = None) # reference at clock input provided

 if not opts.skip:
  print 'Attempting to initialize ADC chips...'
  sys.stdout.flush()
        # try initializing a few times for good measure in case it fails...
  done = False
  for i in range(3):
   if adc.init(samplingRate= 250, numChannel=4, resolution=8) == 0:
    done = True
    break
  print 'done (took %d attempts)' % (i+1)
  
  if not done:
   print 'Failed to calibrate after %d attempts' % (i+1)
   exit_clean()

 c1 = fpga.estimate_fpga_clock()
 print 'FPGA board clock is %d MHz' %c1


 adc.selectADC() # send commands to all the ADC chips
 adc.adc.selectInput([1,2,3,4]) #No  Interleaving



  #Configure registers
 print 'Setting FFT shift to %x' % opts.fftshift
 fpga.write_int('fft_shift', opts.fftshift & 0xffff)
 print 'Checking for FFT overflows...'
 oflow = False
 
 for i in range(5):
  oflow = oflow or bool(fpga.read_int('fft_of'))
  time.sleep(1)
 
 if oflow:
  print 'Overflows detected -- consider increasing FFT shift'
 else:
  print 'No overflows detected'

 print 'Setting accumulation length to %.2f spectra' % opts.acc_len #* 2 * 512/ (c1*10e6)
 fpga.write_int('acc_len', opts.acc_len)
  

 print 'Triggering sync'
 fpga.write_int('cnt_rst', 0)
 fpga.write_int('sw_sync', 0)
 fpga.write_int('sw_sync', 1)
    
 trig_time = time.time()
 fpga.write_int('sw_sync', 0)
 fpga.write_int('cnt_rst', 1)
 fpga.write_int('cnt_rst', 0)

 this_acc = 0
 this_acc_time = trig_time
 file_start_time = time.time()

 data  = []
 times = []
 n1=opts.nnumber
 while(True):
  try:
   latest_acc = fpga.read_int('acc_cnt')
   latest_acc_time = time.time()
   if latest_acc == this_acc:
      time.sleep(0.05)
   elif latest_acc == this_acc + 1:
      print 'Got %d accumulation after %.2f seconds' % (latest_acc, (latest_acc_time - this_acc_time))
      data += [get_data(fpga)]
      times += [latest_acc_time]
      this_acc = latest_acc
      this_acc_time = latest_acc_time
      if time.time() > (file_start_time + opts.filetime):
        write_file(c1,n1,times,data)
        file_start_time = time.time()
        data  = []
        times = []
   else:
        print 'Last accumulation was number %d' % this_acc,
        print 'Next accumulation is number %d' % latest_acc,
        print 'Bad!'
        this_acc = latest_acc
        this_acc_time = latest_acc_time
  except KeyboardInterrupt:
   print 'Exiting'
   write_file(c1,n1,times,data)
   exit()
