import scipy.signal
import scipy.interpolate
import numpy as np
import matplotlib.pyplot as plt
import csv
import datetime as dt

data = []
t = []
fs_vlf = 1.42857e7 # Hz
fs = 131072.  

with open('/home/rileyannereid/Downloads/scope_50usdiv_5mVdiv.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for ri, row in enumerate(csv_reader):
        if ri > 15:
            t.append((ri-16)*1/fs_vlf)
            data.append(float(row[1]))
data = np.array(data)

# resample
n_samples = 8192
data_dt_vlf = dt.timedelta(microseconds=1e6/fs_vlf) # convert time delta to datetime obj.
time_vec_vlf = [data_dt_vlf*i for i in range(int(n_samples))] # create time vec

# create a timevec for the data at desired sampling freq.
data_dt = dt.timedelta(microseconds=1e6/fs) # convert time delta to datetime obj. - NOT WORKING roundoff error
time_vec = [data_dt*i for i in range(int(fs * n_samples / fs_vlf))] # create time ve

# interpolate w a linear func 
t_vlf = np.linspace(0, len(time_vec_vlf), num=len(time_vec_vlf), endpoint=True)
t_fs = np.linspace(0, len(time_vec_vlf), num=len(time_vec), endpoint=True)

f_x = scipy.interpolate.interp1d(t_vlf, data)
data_rs = f_x(t_fs)

print(len(data_rs))
#overlap = 0.5
nfft = 1024
window = 'hanning'

ff, tt, FE = scipy.signal.spectrogram(data, fs=fs_vlf, window=window,
                                      nperseg=nfft, mode='psd', scaling='density')

E_mag = np.log10(FE)
pcm = plt.pcolormesh(tt, ff/1e3, E_mag, cmap='jet')
cbar = plt.colorbar(pcm)
cbar.set_label('log(uV^2/m^2/Hz)')
plt.title('signal test 1-17')
plt.xlabel('sec')
plt.ylabel('Frequency [kHz]')
plt.show()
plt.close()