import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
"""
STEP 2: Smooth, Despike, and Segment Data

This script regularly samples data time-series at near-acquisition-equivalent sampling rates, regularizing
sampling times across the experiment and filling minor data gaps that arise in the course of processing.

An iterative smoothing routine using a boxcar averaging kernel is then used to 

Data are saved in CSV format with date-time strings that parse quickly using the Pandas package
in successive processing steps

ALSO GENERATES SUPPLEMENTARY FIG. S1

"""
# Relative path to repository root directory
ROOT = os.path.join('..','..','..')
IDIR = os.path.join(ROOT,'processed_data','timeseries','STEP_0')
ODIR = os.path.join(ROOT,'processed_data','timeseries','STEP_1')
OFDIR = os.path.join(ROOT,'results','figures','supplement')

# Load Data
df_NTP = pd.read_csv(os.path.join(IDIR,'NTP_Local_PD_DateTime.csv'),parse_dates=True,index_col=[0])
df_LVDT = pd.read_csv(os.path.join(IDIR,'LVDT_Local_PD_DateTime_FULL.csv'),parse_dates=True,index_col=[0])
print('Data loaded')


def evenly_sample(df):
	# Make a copy & sort index

	df_ = df.copy().sort_index()
	# Get time index sampling raate
	dt_sec = (df_.index[1:] - df_.index[:-1]).total_seconds()
	dt_bar = np.round(np.median(dt_sec),decimals=0)
	# Resample to equally sampled representation
	df_r = df_.resample(pd.Timedelta(dt_bar,unit='sec')).median()
	return df_r


def z_despike(df,threshscale = 3, spikewin=pd.Timedelta(5,unit='min')):
	df_ = df.copy()
	df_O = pd.DataFrame()
	df_O.index = df_.index
	for i_,C_ in enumerate(df_.columns):
		print('Running %s'%(C_))
		# Get rolling median of data & use for despiking
		S_ = df_[C_].copy().rolling(spikewin).median()
		S_.index -= spikewin/2.
		S_.name = '%s_med'%(C_)
		_df_ = pd.concat([df_[C_].copy(),S_],axis=1,ignore_index=False)
		S_Z = (_df_[C_] - _df_[S_.name])
		IND = np.abs(S_Z) <= threshscale*S_Z.std()
		df_O = pd.concat([df_O,_df_[C_][IND]],axis=1,ignore_index=False)
	return df_O

# Trim out flat spot in LVDT
IND = (df_LVDT.index >= pd.Timestamp("2021-11-06T19:30")) & (df_LVDT.index <= pd.Timestamp("2021-11-07T18"))
df_LVDT = df_LVDT[~IND]

# Resample
df_NTPr = evenly_sample(df_NTP)
# Remove mean pressure from normal stress
S_Nc = df_NTPr['N_kPa'] - df_NTPr[['Pw1_kPa','Pw2_kPa']].mean(axis=1)
S_Nc.name = 'N_kPa'
df_NTr = pd.concat([S_Nc,df_NTPr['T_kPa']],axis=1,ignore_index=False)

df_LVDTr = evenly_sample(df_LVDT)


# Despike
df_NTrd = z_despike(df_NTr,threshscale=2)
df_LVDTrd = z_despike(df_LVDTr,threshscale=2)


issave = True
DPI = 300
FMT = 'PNG'

fig, axs = plt.subplots(ncols=1,nrows=4,figsize=(8,10))
axs[0].plot(df_NTPr['N_kPa'],label='raw')
axs[0].plot(df_NTrd['N_kPa'],label='$P_w$ removed + despiked')
axs[0].set_ylabel('Vertical pressure (kPa)')
axs[0].legend()
axs[1].plot(df_NTPr['Pw1_kPa'],label='Transducer 1')
axs[1].plot(df_NTPr['Pw2_kPa'],label='Transducer 2')
axs[1].legend(loc='upper center')
axs[1].set_ylabel('Water pressure (kPa)')
axs[2].plot(df_NTPr['T_kPa'],label='raw')
axs[2].plot(df_NTrd['T_kPa'],label='despiked')
axs[2].set_ylabel('Shear stress (kPa)')
axs[2].legend()
axs[3].plot(df_LVDTr['LVDT_mm'],label='shift corrected')
axs[3].plot(df_LVDTrd['LVDT_mm'],label='despiked')
axs[3].set_ylabel('LVDT Position (mm)')
axs[3].set_xlabel('Date / Time (UTC - 5 h)')
axs[3].legend()
for i_,l_ in enumerate(['a','b','c','d']):
	if i_ < 3:
		axs[i_].xaxis.set_visible(False)
	ylim = axs[i_].get_ylim()
	axs[i_].text(df_NTPr.index[0],(np.max(ylim) - np.min(ylim))*0.9 + np.min(ylim),l_,\
				 fontweight='extra bold',fontstyle='italic',fontsize=16,va='center',ha='right')
plt.subplots_adjust(wspace=0,hspace=0)
plt.show()
if issave:
	OFILE = os.path.join(OFDIR,'CNSB_FigS1_Raw_Timeseries_diss_v2_%ddpi.%s'%(DPI,FMT.lower()))
	plt.savefig(OFILE,dpi=DPI,format=FMT.lower())


#Smooth data
siter = 3
smoothwin_min =10
smoothwin = pd.Timedelta(smoothwin_min,unit='min') 
df_NTrds = df_NTrd.copy().rolling(smoothwin).mean()
df_LVDTrds = df_LVDTrd.copy().rolling(smoothwin).mean()
# Window-center data
df_NTrds.index -= smoothwin/2
df_LVDTrds.index -= smoothwin/2
# Trim off leading nubbin
df_NTrds = df_NTrds[(df_NTrds.index >= df_NTr.index.min()) & \
					 (df_NTrds.index <= df_NTr.index.max())]
df_LVDTrds = df_LVDTrds[(df_LVDTrds.index >= df_LVDTr.index.min()) & \
					    (df_LVDTrds.index <= df_LVDTr.index.max())]


if siter > 1:
	for i_ in range(siter - 1):
		df_NTrds = df_NTrds.rolling(smoothwin).mean()
		df_LVDTrds = df_LVDTrds.rolling(smoothwin).mean()

		# Window-center data
		df_NTrds.index -= smoothwin/2
		df_LVDTrds.index -= smoothwin/2
		# Trim off leading nubbin
		df_NTrds = df_NTrds[(df_NTrds.index >= df_NTr.index.min()) & \
							  (df_NTrds.index <= df_NTr.index.max())]
		df_LVDTrds = df_LVDTrds[(df_LVDTrds.index >= df_LVDTr.index.min()) & \
							    (df_LVDTrds.index <= df_LVDTr.index.max())]


### WRITE FULL TIME-SERIES TO DISK ###
df_NTrds.to_csv(os.path.join(ODIR,'FULL_NT_%smin_x%s_smoothed.csv')%(smoothwin_min,siter))
df_LVDTrds.to_csv(os.path.join(ODIR,'FULL_LVDT_%smin_x%s_smoothed.csv')%(smoothwin_min,siter))


### CONDUCT SEGMENTATION ###
# Cycle start and end times 24 hr and 6 hr
TSV = [pd.Timestamp("2021-10-26T02:00:00"),pd.Timestamp('2021-10-26T13:56'),\
	   pd.Timestamp("2021-10-31T15:40"),pd.Timestamp('2021-11-1T11:09:15'),\
	   pd.Timestamp("2021-11-03T15:52")]
TEV = [pd.Timestamp("2021-10-26T13"),pd.Timestamp("2021-10-31T14:03"),\
	   pd.Timestamp("2021-11-01T10:30"),pd.Timestamp("2021-11-02T17:26:30"),\
	   pd.Timestamp("2021-11-22T12:09")]
DT_PAD = pd.Timedelta(6,unit='hour')

df_NTS1 = df_NTrds[(df_NTrds.index >= TSV[0] - DT_PAD) & (df_NTrds.index <= TEV[0] + DT_PAD)]
df_NTS1.to_csv(os.path.join(ODIR,'SS1_NT_%smin_x%s_smoothed.csv'%(smoothwin_min,siter)))

df_NT24 = df_NTrds[(df_NTrds.index >= TSV[1] - DT_PAD) & (df_NTrds.index <= TEV[1] + DT_PAD)]
df_NT24.to_csv(os.path.join(ODIR,'T24_NT_%smin_x%s_smoothed.csv'%(smoothwin_min,siter)))

df_NTS2 = df_NTrds[(df_NTrds.index >= TSV[2] - DT_PAD) & (df_NTrds.index <= TEV[2] + DT_PAD)]
df_NTS2.to_csv(os.path.join(ODIR,'SS2_NT_%smin_x%s_smoothed.csv'%(smoothwin_min,siter)))

df_NT06 = df_NTrds[(df_NTrds.index >= TSV[3] - DT_PAD) & (df_NTrds.index <= TEV[3] + DT_PAD)]
df_NT06.to_csv(os.path.join(ODIR,'T06_NT_%smin_x%s_smoothed.csv'%(smoothwin_min,siter)))

df_NT96 = df_NTrds[(df_NTrds.index >= TSV[4] - DT_PAD) & (df_NTrds.index <= TEV[4] + DT_PAD)]
df_NT96.to_csv(os.path.join(ODIR,'T96_NT_%smin_x%s_smoothed.csv'%(smoothwin_min,siter)))


df_LVDTS1 = df_LVDTrds[(df_LVDTrds.index >= TSV[0] - DT_PAD) & (df_LVDTrds.index <= TEV[0] + DT_PAD)]
df_LVDTS1.to_csv(os.path.join(ODIR,'SS1_LVDT_%smin_x%s_smoothed.csv'%(smoothwin_min,siter)))

df_LVDT24 = df_LVDTrds[(df_LVDTrds.index >= TSV[1] - DT_PAD) & (df_LVDTrds.index <= TEV[1] + DT_PAD)]
df_LVDT24.to_csv(os.path.join(ODIR,'T24_LVDT_%smin_x%s_smoothed.csv'%(smoothwin_min,siter)))

df_LVDTS2 = df_LVDTrds[(df_LVDTrds.index >= TSV[2] - DT_PAD) & (df_LVDTrds.index <= TEV[2] + DT_PAD)]
df_LVDTS2.to_csv(os.path.join(ODIR,'SS2_LVDT_%smin_x%s_smoothed.csv'%(smoothwin_min,siter)))

df_LVDT06 = df_LVDTrds[(df_LVDTrds.index >= TSV[3] - DT_PAD) & (df_LVDTrds.index <= TEV[3] + DT_PAD)]
df_LVDT06.to_csv(os.path.join(ODIR,'T06_LVDT_%smin_x%s_smoothed.csv'%(smoothwin_min,siter)))

df_LVDT96 = df_LVDTrds[(df_LVDTrds.index >= TSV[4] - DT_PAD) & (df_LVDTrds.index <= TEV[4] + DT_PAD)]
df_LVDT96.to_csv(os.path.join(ODIR,'T96_LVDT_%smin_x%s_smoothed.csv'%(smoothwin_min,siter)))






