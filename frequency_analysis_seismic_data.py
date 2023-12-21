## ObsPy packages

import obspy.signal
from obspy.core import read, UTCDateTime, AttribDict, Stream
from obspy.signal.trigger import plot_trigger
from obspy.signal.trigger import carl_sta_trig
from obspy.signal.freqattributes import spectrum, welch
from obspy.signal.invsim import corn_freq_2_paz
from obspy.signal.tf_misfit import plot_tfr, plot_tf_gofs, plot_tf_misfits
from obspy.signal.cross_correlation import _pad_zeros
from obspy.imaging.spectrogram import spectrogram
from obspy.io.ascii.core import _is_slist, _read_slist, _read_tspair

## Numpy

import numpy as np
from numpy import shape

## Scipy

import scipy.signal
from scipy import pi
from scipy.fftpack import fft, rfft
from scipy.signal import get_window

## Matplotlib

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import gridspec

## Rest

import os
import glob
import linecache
import math
from geopy import distance
from pylab import *
from sys import exit
from PIL import Image

#---------------------------------------------------------------------------------

## Specific dataset files

cluster1_data = open('locations_2017_DD_AMT_clust1.txt', 'r')
cluster2_data = open('locations_2017_DD_AMT_clust2.txt', 'r')
cluster3_data = open('locations_2017_DD_AMT_clust3.txt', 'r')
cluster4_data = open('locations_2017_DD_AMT_clust4.txt', 'r')

data_merged = open('locations_2017_merged_list.txt', 'r')
stat_catalog = open('reykjanet_stanice_new.txt', 'r')

## Creating empty lists for categorising data

event_list =[]
mag_list = []
#lines_list = []
event_lat = []
event_lon = []
depth_list = []
cluster_list = []

for line in data_merged:
#for line in cluster1_data:
    line_elements = line.split()
    event_n = line_elements[11]
    lat = line_elements[7]
    lon = line_elements[8]
    depth = line_elements[9]
    event_mag = line_elements[10]
    cluster_nr = line_elements[18]
    
    event_list.append(event_n)
    depth_list.append(depth)
    mag_list.append(event_mag)
    event_lat.append(lat)
    event_lon.append(lon)
    cluster_list.append(cluster_nr)
    
#print(event_list[73])
#print(depth_list[73])
#print(mag_list[73])
#print(event_lat[73])
#print(event_lon[73])
#print(cluster_list[73])
#exit()
#event_list.pop(0)
#mag_list.pop(0)
#event_lat.pop(0)
#event_lon.pop(0)
#depth_list.pop(0)

## Seismic stations naming

##Stations
ASH = 'ASH'
FAF = 'FAF'
HRG = 'HRG'
KLV = 'KLV'
LSF = 'LSF'
SEA = 'SEA'
ELB = 'ELB'
LAT = 'LAT'
ISS = 'ISS'
GEI = 'GEI'
LHL = 'LHL'
LAG = 'LAG'
HDV = 'HDV'
STH = 'STH'
MOH = 'MOH'
SEL = 'SEL'

#stat_list = ('ASH', 'FAF', 'HRG', 'KLV', 'LSF', 'SEA', 'ELB', 'LAT', 'ISS', 'GEI', 'LHL', 'LAG', 'HDV',
#             'STH','MOH', 'SEL')

## Lists for stations parameters

stat_list = []
stat_lat = []
stat_lon = []

## Categorizing data for distance calculation (between hypocenter and geophone)

for line1 in stat_catalog:
    stat_name = line1[7:10]
    st_lat = line1[21:28]
    st_lon = line1[29:37]
    stat_list.append(stat_name)
    stat_lat.append(st_lat)
    stat_lon.append(st_lon)

## The range can be specified based on current analysis

for a in range (0, len(event_list)):
#for a in range (0, 1):

 if cluster_list[a] == '1':  
    print('Event number: ', event_list[a])
    mag_str = str(mag_list[a])
    event_str = str(event_list[a])
    event_lat_str = str(event_lat[a])
    event_lon_str = str(event_lon[a])
    depth_fl = float(depth_list[a])
    print('Depth =' + str(depth_fl) + ' km')
    
    for b in range (0, len(stat_list)):
    #for b in range (9, 10):
        try:
            stat_str = str(stat_list[b])
            stat_lat_str = str(stat_lat[b])
            stat_lon_str = str(stat_lon[b])
            print('Station: ', stat_list[b])

            name_combi = 'figures/cluster_1/P_wave_record_1s/' + event_str + '/' + event_str + '_' + stat_str + '_combi_plot.png'
            #name_combi = 'figures/P_wave_record_1s_1.5-15_Hz/' + event_str + '/' + event_str + '_' + stat_str + '_combi_plot.png'

            data_e = 'Data3_asc/A' + event_str + '/' + event_str + '_' + stat_str + '_CHE.ascii'
            data_n = 'Data3_asc/A' + event_str + '/' + event_str + '_' + stat_str + '_CHN.ascii'
            data_z = 'Data3_asc/A' + event_str + '/' + event_str + '_' + stat_str + '_CHZ.ascii'

            ## Get time of current event

            with open(data_e) as f:
                first_line = f.readline()
            time_orig = UTCDateTime(first_line[51:74])
            stime = time_orig + 1.5
            etime_short = stime + 1.0
            etime_long = stime + 1.0
            

            ## Calculate distance between drill and hypocenter

            event_coord = (event_lat_str, event_lon_str)
            stat_coord = (stat_lat_str, stat_lon_str)
            
            epi_dist = distance.distance(event_coord, stat_coord).km
            hyp_dist = math.sqrt(epi_dist**2 + depth_fl**2)

            ep = round(epi_dist, 3)
            hp = round(hyp_dist, 3)

            print('Epi.dist =' + str(ep) + ' km')
            print('Hyp.dist =' + str(hp) + ' km')
                
            
            st_z=read(data_z)
            st_e=read(data_e)
            st_n=read(data_n)

            f_min = 1.5
            f_max = 15.0
            str_f_min = str(f_min)
            str_f_max = str(f_max)
           

            st = st_z    ### only z for spectrograms

            ## Specifying lenght of the record for analysis

            if hp > 10:
                st=st.trim(starttime=stime, endtime=etime_long)
            if hp <= 10:
                st=st.trim(starttime=stime, endtime=etime_short)

            tr = st[0]
            for tr in st:
                tr.detrend()
                tr.taper(max_percentage=0.1)
                tr.filter('bandpass', freqmin=f_min, freqmax=f_max, corners=2, zerophase=True)

            ## Saving generating images

            plot_tfr(tr.data, dt=0.004, fmin=1.5, fmax=50, w0=6, nf=120, fft_zero_pad_fac=6, show=False)
            #plt.show()
            #exit()
            plt.savefig('tfr_plot.png')
            plt.close("all")

            #st[0].spectrogram(samp_rate=250, per_lap = 0.95, wlen=1.0, log=True, mult=100.0, title='Spectrogram', show=False)
            #plt.xlim(0, 7)
            #plt.ylim(2, 50)
            #plt.savefig('tfr_plot.png')
            #plt.close("all")
            #plt.show()
            #exit()

            ##### adding components unless it's STH station or LAT

                ### rewriting Z component to avoid double filtering
            
            st_z2=read(data_z)
            
            st = st_z2
            
            if stat_str in ('ASH', 'FAF', 'HRG', 'KLV', 'LSF', 'SEA', 'LAT', 'ELB', 'ISS', 'GEI', 'LHL', 'LAG', 'HDV', 'MOH', 'SEL'):
                st += st_n

            if stat_str in ('ASH', 'FAF', 'HRG', 'KLV', 'LSF', 'SEA', 'ELB', 'ISS', 'GEI', 'LHL', 'LAG', 'HDV', 'MOH', 'SEL'):
                st += st_e

            ## Specifying lenght of the record for analysis

            if hp > 10:
                st=st.trim(starttime=stime, endtime=etime_long)
            if hp <= 10:
                st=st.trim(starttime=stime, endtime=etime_short)
                

            tr = st[0]
            for tr in st:
                tr.detrend()
                tr.taper(max_percentage=0.1)
                tr.filter('bandpass', freqmin=f_min, freqmax=f_max, corners=2, zerophase=True)

            #exit()

            ## Plotting relative data
            
            st.plot(type='relative', size=(1050, 900), outfile='st_plot.png')
            #st.plot(type='relative')
            #st.plot(type='relative', size=(1050, 300), outfile='st_plot.png')
            #plt.locator_params(axis='x', nbins=2)
            #plt.show()
            plt.close('all')
            #exit()


            ## Combined figure

            suptit = str('Event number: ' + event_str + ', Station: ' + stat_str +
                         ', Ml=' + mag_str + ', fmin = 1.5 Hz, fmax = 15.0 Hz')

            img1 = mpimg.imread('st_plot.png')
            img2 = mpimg.imread('tfr_plot.png')
            #print(img)

            
            fig = plt.figure(figsize = (15,17))
            #gs = gridspec.GridSpec(2, 1, width_ratios=[1, 1.2])
            fig.suptitle(suptit, fontsize=20)
            
            
            #a = fig.add_subplot(2,1,1)
            a = fig.add_axes([0.25, 0.55, 0.47, 0.4])
            imgplot = plt.imshow(img1)
            #plt.subplots_adjust(left=0.0125, bottom=0.14, right=0.9, top=0.9, wspace=0.01, hspace=0.01)
            plt.axis('off')
            
            #a = fig.add_subplot(2,1,2)
            a = fig.add_axes([0.07, 0.001, 0.7, 0.6])
            #plt.tight_layout()
            imgplot = plt.imshow(img2)
            plt.axis('off')

            plt.savefig(name_combi)
            #plt.savefig('extra.png')
            plt.close('all')
            #plt.show()
            #exit()

            ## Spectrograms
            
            #st[0].spectrogram(samp_rate=250, per_lap = 0.95, wlen=1.0, log=True, mult=100.0, title='Spectrogram', show=False)
            #plt.xlim(0, 7)
            #plt.ylim(2, 50)
            #plt.show()
            #plt.savefig('name3.fig')
        except:
            print('Station ' + stat_list[b] + ' out of order')
            #exit()
 else:
     print('out of cluster')

exit()

#---------------------------------------------------------------------

## Counting the lines

count = 0
content = file1.read()
colist = content.split("\n")

for i in colist:
    if i:
        count += 1

 ###count is 323 lines

list_time = []
list_lat = []
list_lon = []

## Lists of events for further analysis

## List 2 ---with well seen phases

list2 = [5, 21, 22, 23, 24, 25, 26, 59, 64, 65, 67, 106, 109, 112, 123, 124, 126,
         145, 158, 164, 211, 212, 214, 237, 247, 276]

len_list2 = len(list2)

## START LOOP 1 for specific range

for i in range(1, 2):
        d=linecache.getline(data_catalog, i)

        ## START LOOP 2 for list events

#for i in range(0, len(list2)):
#        a = list2[i]
#        d=linecache.getline(data_catalog, a)
        
## Rest of the loop  
#       
        year = int(d[0:4])
        month = int(d[5:7])
        day = int(d[8:10])
        day_str = str(d[8:10])
        hour = int(d[22:24])
        minute = int(d[25:27])
        second = int(d[28:30])

        event_no = int(d[11:20])
        event_no_str = str(d[11:20])
        name1=('figures/visual_top/' + event_no_str + '_top_zoom1_5Hz')
        name2=('figures/visual_top/' + event_no_str + '_borehole_zoom1_5Hz')
        name3=('figures/fft_visual_top_Z/' + event_no_str + '_Z_top_freq_time_spec_5Hz')
        name4=('figures/fft_visual_top_Z/' + event_no_str + '_Z_borehole_freq_time_spec')
        name5=('figures/fft_visual_top_Z/' + event_no_str + '_Z_top_fft')
        name6=('figures/fft_visual_top_Z/' + event_no_str + '_Z_borehole_fft')
        name7=('figures/visual_top/' + event_no_str + '_Z_fft_spec_borehole')

        event_h = float(d[68:76])
        magnitude = float(d[90:97])

        lat = d[38:47]
        lon = d[52:62]
        list_lat.append(lat)
        list_lon.append(lon)

        ## Calculate distance between drill and hypocenter

        coords_1 = (F2_bot_lat, F2_bot_lon)
        coords_2 = (lat, lon)
        coords_3 = (F2_surf_lat, F2_surf_lon)
        
        epi_dist = distance.distance(coords_1, coords_2).km
        hyp_dist = math.sqrt(epi_dist**2 + event_h**2)

        ## Print info about event

        print('Event number, depth, epicentral distance, hypocentral distance in [km] and magnitude:')
        print(i)
        print(event_no_str)
        #print(event_h)
        #print(round(epi_dist, 3))
        #print(round(hyp_dist, 3))
        #print(magnitude)

        ## Defining the start time

        origin_t = UTCDateTime(year, month, day, hour, minute, second)

        if hyp_dist <= 20:
            stime = origin_t
            list_time.append(stime)
            etime = origin_t + 10

        if hyp_dist > 20:
            #continue
            dist_koeff = (hyp_dist/4.2)
            #stime = origin_t - 2 + dist_koeff
            stime = origin_t
            list_time.append(stime)
            #etime = stime + 10
            etime = stime + 8 + dist_koeff

        #if magnitude < 0:
        #    continue

        ## Read data

        data1_z = 'mseed_b6p/c0b6p200' + str(month) + day_str + '??????.pri0'
        data1_x = 'mseed_b6p/c0b6p200' + str(month) + day_str + '??????.pri1'
        data1_y = 'mseed_b6p/c0b6p200' + str(month) + day_str + '??????.pri2'

        data2_z = 'mseed_b6w/c0b6w200' + str(month) + day_str + '??????.pri0'
        data2_x = 'mseed_b6w/c0b6w200' + str(month) + day_str + '??????.pri1'
        data2_y = 'mseed_b6w/c0b6w200' + str(month) + day_str + '??????.pri2'

        ## Fixing the labels
        
        b=0
        
        st_z=read(data1_z)
        for tr in st_z:
            if tr.stats.channel == 'p0':
                st_z[b].stats.channel = 'Z'
        st_x=read(data1_x)
        for tr in st_x:
            if tr.stats.channel == 'p1':
                st_x[b].stats.channel = 'N'
        st_y=read(data1_y)
        for tr in st_y:
            if tr.stats.channel == 'p2':
                st_y[b].stats.channel = 'E'
        
        st = st_z
        #st += st_x
        #st += st_y

        c=0
        
        st2_z=read(data2_z)
        for tr in st2_z:
            if tr.stats.channel == 'p0':
                st2_z[c].stats.channel = 'Z'
        st2_x=read(data2_x)
        for tr in st2_x:
            if tr.stats.channel == 'p1':
                st2_x[c].stats.channel = 'N'
        st2_y=read(data2_y)
        for tr in st2_y:
            if tr.stats.channel == 'p2':
                st2_y[c].stats.channel = 'E'

        st2 = st2_z
        #st2 += st2_x
        #st2 += st2_y

        ####Trim
        
        st=st.trim(starttime=stime, endtime=etime)
        st2=st2.trim(starttime=stime, endtime=etime)
        
        ## Downsample/Upsample data
        
        st_copy=st.copy()
        st_int=st_copy.interpolate(sampling_rate=100.0)
        st_res=st_copy.resample(sampling_rate=100.0)
        #st=st_res
        #st=st
        #print(st)
        #print(st_int)
        #exit()


        ## Finding local maximum aplitude for equal scaling up to 1
        
        list_max = []
        #st=st.trim(starttime=stime, endtime=etime)
        tr = st[0]
        for tr in st:
            tr.detrend()
            tr.taper(max_percentage=0.1)
            if event_no == 20202127:
                tr.filter('bandpass', freqmin=5.0, freqmax=49.0, corners=2, zerophase=True)
            else:
                tr.filter('bandpass', freqmin=2.0, freqmax=49.0, corners=2, zerophase=True)
            max1 = max(tr.data)
            min1 = sqrt((min(tr.data)**2))
            list_max.append(max1)
            list_max.append(min1)
            

        #st2=st2.trim(starttime=stime, endtime=etime)
        tr2 = st2[0]
        for tr2 in st2:
            tr2.detrend()
            tr2.taper(max_percentage=0.1)
            tr2.filter('bandpass', freqmin=2.0, freqmax=49.0, corners=2, zerophase=True)
            max2 = max(tr2.data)
            min2 = sqrt((min(tr2.data)**2))
            list_max.append(max2)
            list_max.append(min2)

        max_equal = max(list_max)
       
        ## EXTRAS -------------------------------------------

        ### Plotting with use of the scaling
        
        #for tr in st:
        #    tr.data = tr.data / max_equal

        #for tr2 in st2:
        #    tr2.data = tr2.data / max_equal

        #---------------------------------------------------
            ### Z component for FFT

        #st=st.trim(starttime=stime, endtime=etime)
        #tr = st[0]
        #tr_z.plot()

        #st2=st2.trim(starttime=stime, endtime=etime)
        #tr2 = st2[0]
        #tr2_z.plot()

        #a = np.array(tr)
        #print(a.shape)
        #exit()


        #---------------------------------------------------

        ## Downsampling

        #tr2_copy = tr2.copy()
        #tr2_copy.decimate(factor=4, strict_length=False)
        #print(tr2_copy.stats)
        #exit()


        #---------------------------------------------------
            ## FFT

            ## Use of spectrum function (inside ObsPy package)


        #sp1 = spectrum(tr.data, win=1, nfft = 1001)
        #sp2 = spectrum(tr2.data, win=1, nfft = 4001)

            ### cut in hlaf to get rid of mirror effect
        #sp1_half = sp1[:len(sp1)//2]
        #sp2_half = sp2[:len(sp2)//2]

        #max_spec1 = max(sp1_half)
        #max_spec2 = max(sp2_half)
        #list_max = [max_spec1, max_spec2]
        #max_spec = max(list_max)
        #print(max_spec)
        

        #x_freq1 = np.linspace(0.0, 50.0, 500)
        #x_freq2 = np.linspace(0.0, 200.0, 500)



        #plt.plot(x_freq2, sp2_half, label="borehole")
        #plt.plot(x_freq1, sp1_half, label="surface")
        #plt.title('Z Frequency domain compare' + event_no_str)
        #plt.xlabel('Frequency [Hz]')
        #plt.xlim(0, 80)
        #plt.ylim(-10, max_spec+20)
        #plt.ylabel('Amplitude')
        #plt.legend(loc='upper right')
        #plt.savefig(name7)
        #plt.show()

        #---------------------------------------
            
        #st.plot(size=(1500, 900), outfile=name1)
        #st2.plot(size=(1500, 900), outfile=name2)
            
        #st.plot(type='relative', size=(1500, 900), outfile=name1)
        #st2.plot(type='relative', size=(1500, 900), outfile=name2)

        #st.plot(type='relative', size=(1200, 650))
        #st2.plot(type='relative', size=(1200, 650))

        #st_z.plot(type='relative', size=(1200, 650))
        #st2_z.plot(type='relative', size=(1200, 650))

        #plot_tfr(tr.data, dt=0.01, fmin=2.0, fmax=75.0, w0=8, nf=64, fft_zero_pad_fac=4, show=False)
        #plt.show()
        #plt.savefig(name3)

        #plot_tfr(tr2_copy.data, dt=0.01, fmin=2.0, fmax=75.0, w0=8, nf=64, fft_zero_pad_fac=4, show=False)
        #plt.show()
        #plt.savefig(name4)


        ## Additional code snippets --------------------------------------------------

        
        #st_z=st_z.trim(starttime=origin_minus, endtime=etime)
        #tr_z= st_z[0]
        #for tr_z in st_z:
        #    tr_z.detrend()
        #    tr_z.taper(max_percentage=0.1)
        #    tr_z.filter('bandpass', freqmin=2.0, freqmax=49.0, corners=2, zerophase=True)
            # scale amplitudes
        #    tr_z.data = tr_z.data / max(tr_z.data)
            #tr_z.plot()

        #st2_z=st2_z.trim(starttime=origin_minus, endtime=etime)
        #tr2_z= st2_z[0]
        #for tr2_z in st2_z:
        #    tr2_z.detrend()
        #    tr2_z.taper(max_percentage=0.1)
        #    tr2_z.filter('bandpass', freqmin=2.0, freqmax=49.0, corners=2, zerophase=True)
            # scale amplitudes
        #    tr2_z.data = tr2_z.data / max(tr2_z.data)
            #tr2_z.plot()

        ## Ploting freqency spectrum

        #plot_tfr(tr.data, dt=0.0025, fmin=2.0, fmax=75.0, w0=8, nf=64, fft_zero_pad_fac=4, show=False)
        #plt.show()
        #plt.savefig(name3)

        #plot_tfr(tr2.data, dt=0.0025, fmin=2.0, fmax=49.0, w0=6, nf=64, fft_zero_pad_fac=4, show=False)
        #plt.show()
        #plt.savefig(name4)


        ## Testing coordinates (has to return 279.35 km)
        
        #coords_1 = (52.2296756, 21.0122287)
        #coords_2 = (52.406374, 16.9251681)
        #print('Event number: ', event_no_str)
        #print(stime)
        #print('Distance from F2: ', distance.distance(coords_1, coords_2).km, 'Km')

        ## Calculating spectrum

        #spc1 = spectrum(tr.data, win = 10, nfft = 10)
        #spc2 = spectrum(tr2.data, win = 10, nfft = 10)
        #print(spc1)
        #print(spc2)

        #---------------------------------------
            ## FFT and plotting frequency spectrum

        #tr2.decimate(factor=4)
        
        sr1 = tr.stats.sampling_rate
        N1 = len(tr.data)
        dt1 = tr.stats.delta

        sr2 = tr2.stats.sampling_rate
        N2 = len(tr2.data)
        dt2 = tr2.stats.delta

        print(sr1, N1, dt1, sr2, N2, dt2)


        #a = np.array(tr)
        #print(len(a.shape))
        #exit()

        Fsig = np.fft.rfft(tr.data, n=N1)
        Fsig2 = np.fft.rfft(tr2.data, n=N2)

        #### zero padding and inverse fft

        N = 1500
        Fsig_pad1 = np.pad(Fsig, (0, N), 'constant')
        Fsig_pad2 = _pad_zeros(Fsig, 0, 1500)

        st_pad = np.fft.irfft(Fsig, n=4001)
        t = np.arange(4001)
        #tr_pad = st_pad[0]
        N3 = len(st_pad)
        print(N3)
        #sr3 = tr_pad.stats.sampling_rate
        #print(sr3)
        plt.plot(t, st_pad)
        plt.show()
        exit()

        #a2 = np.array(tr2)
        #print(len(a2.shape))
        #exit()
        
        freq1 = np.linspace(0.0, int(sr1/2), int((N1/2)+1))
        freq2 = np.linspace(0.0, int(sr2/2), int((N2/2)+1))


        #plt.plot(freq2, y2, label="borehole")
        plt.plot(freq1, dt1 * np.abs(Fsig), label="surface")
        plt.plot(freq2, dt1 * np.abs(Fsig_pad1), label = "surface, pad")
        #plt.plot(freq2, dt2 * np.abs(Fsig2), label="borehole")
        plt.title('Frequency domain borehole' + event_no_str)
        plt.xlabel('Frequency [Hz]')
        plt.xlim(0, 200)
        #plt.ylim(0, 0.23)
        plt.ylabel('Amplitude')
        plt.legend(loc='upper right')
        #plt.savefig(name7)
        plt.show()
        plt.clf()

#plt.show()