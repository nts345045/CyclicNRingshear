import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from matplotlib.gridspec import GridSpec
sys.path.append('..')
import util.datetimeindex as dtiu

issave = True
DPI = 300
FMT = 'PNG'

ROOT = os.path.join('..','..','..')
DDIR = os.path.join(ROOT,'processed_data','STEP_1')
T6_NT = os.path.join(DDIR,'T06_NT_10min_x3_smoothed.csv')
T6_LV = os.path.join(DDIR,'T06_LVDT_10min_x3_smoothed_QSS_reduced.csv')
FDIR = os.path.join(ROOT,'results','figures','main')
ODIR = os.path.join(ROOT,'processed_data','STEP_2')

# Load data into merged data frame
df_6 = pd.concat([pd.read_csv(T6_NT,parse_dates=True,index_col=[0]),\
				   pd.read_csv(T6_LV,parse_dates=True,index_col=[0])],\
				   axis=1,ignore_index=False)


print('data loaded')
# Trim off lead & lags
DT_MASK = pd.Timedelta(3,unit='hour')
IND6 = (df_6.index >= df_6.index.min() + DT_MASK) & (df_6.index <= df_6.index.max() - DT_MASK)

IND6x = (df_6.index >= df_6.index.min() + DT_MASK*2) & (df_6.index <= df_6.index.max() - DT_MASK*2)


# Upsample LVDT Data
df_6r = df_6[IND6].copy().interpolate(method='quadratic',limit=15,limit_direction='both')
df_6x = df_6[IND6x].copy().interpolate(method='quadratic',limit=15,limit_direction='both')
print('data upsampled')



# fig = plt.figure(figsize=(8,4))
nrows,ncols = 3,1
fig,axs = plt.subplots(ncols=ncols,nrows=nrows,figsize=(8,10))
# kick-in shared axes
axs[0].get_shared_x_axes().join(axs[0],*axs[1:])

ylims = {'N_kPa':[180,520],'LVDT_mm red':[-0.5,2.5],'T_kPa':[110,210]}
yticks = {'N_kPa':np.arange(200,550,50),'LVDT_mm red':np.arange(0,3,1),'T_kPa':np.arange(120,210,20)}
labels = {'LVDT_mm red':'Relative ice-bed separation\n(mm)',\
		  'N_kPa':'Effective pressure (kPa)',\
		  'T_kPa':'Shear stress (kPa)'}
TLAB = [6]

T0_6 = df_6r.index.min() + DT_MASK

### ACTUAL PLOTTING ###
for i_,C_ in enumerate(df_6r.columns):
	# Plot data
	axs[i_].plot((df_6r.index - T0_6).total_seconds()/3600,df_6r[C_].values,'k')
	axs[i_].set_ylabel(labels[C_])
	axs[i_].set_xlim([(df_6r.index.min() - T0_6).total_seconds()/3600,\
						(df_6r.index.max() - T0_6).total_seconds()/3600])
	axs[i_].set_ylim(ylims[C_])
	axs[i_].set_yticks(yticks[C_])



# Overlay minima and maxima
for k_, C_ in enumerate(df_6r.columns):
	ext6 = dtiu.pick_extrema_indices(df_6x[C_],T=pd.Timedelta(5,unit='hour'))
	if C_ == 'N_kPa':
		for t_ in ext6['I_min']:
			axs[k_].plot(np.ones(2,)*(t_ - T0_6).total_seconds()/3600,ylims[C_],'b--',alpha=0.25)
			axs[1].plot(np.ones(2,)*(t_ - T0_6).total_seconds()/3600,ylims['T_kPa'],'b--',alpha=0.25)
		for t_ in ext6['I_max']:
			axs[k_].plot(np.ones(2,)*(t_ - T0_6).total_seconds()/3600,ylims[C_],'b-',alpha=0.25)
			axs[1].plot(np.ones(2,)*(t_ - T0_6).total_seconds()/3600,ylims['T_kPa'],'b-',alpha=0.25)


	elif C_ == 'LVDT_mm red':
		for t_ in ext6['I_min']:
			axs[k_].plot(np.ones(2,)*(t_ - T0_6).total_seconds()/3600,ylims[C_],'r--',alpha=0.25)
			axs[1].plot(np.ones(2,)*(t_ - T0_6).total_seconds()/3600,ylims['T_kPa'],'r--',alpha=0.25)
		for t_ in ext6['I_max']:
			axs[k_].plot(np.ones(2,)*(t_ - T0_6).total_seconds()/3600,ylims[C_],'r-',alpha=0.25)
			axs[1].plot(np.ones(2,)*(t_ - T0_6).total_seconds()/3600,ylims['T_kPa'],'r-',alpha=0.25)
#### SHADE IN SELECTED DATA FOR FIGURE 8 #####
		for i_ in range(3):
			axs[i_].fill_between([(ext6['I_max'][1] - T0_6).total_seconds()/3600,\
				 				  (ext6['I_max'][2] - T0_6).total_seconds()/3600],\
								 np.ones(2,)*ylims[df_6r.columns[i_]][0],\
								 np.ones(2,)*ylims[df_6r.columns[i_]][1],\
								 color='yellow',alpha=0.25)
		if k_ == 0:
			df_6r_out = df_6r[(df_6r.index >= ext6['I_max'][1]) & (df_6r.index <= ext6['I_max'][2])]
			df_6r_out.to_csv(os.path.join(ODIR,'T6_hysteresis_subset.csv'),header=True,index=True)

T0v = [T0_6]
for k_ in range(1):
	# axs[0].set_title('%d hour cycling'%(TLAB[k_]))
	axs[2].set_xlabel('Elapsed time since %s\n(hours)'%(T0v[k_]))
plt.subplots_adjust(wspace=0,hspace=0)


for i_,l_ in enumerate(['a','b','c']):
	C_ = df_6r.columns[i_]
	# if i_%2 == 0:
	# 	pct = 0.1
	# else:
	pct = 0.90
	yloc = (ylims[C_][1] - ylims[C_][0])*pct + ylims[C_][0]
	axs[i_].text(-2,yloc,l_,fontweight='extra bold',fontstyle='italic',ha='center',va='center',fontsize=16)

# Save Figure
if issave:
	SAVE_FILE = 'CNSB_Fig6_T6_TimeSeries_v9_%ddpi.%s'%(DPI,FMT.lower())
	plt.savefig(os.path.join(FDIR,SAVE_FILE),dpi=DPI,format=FMT.lower())

#### TODO ####
# Get minimum and maxium cavity geometry timings from videos into row 3
# Outline of Fig. 2 timeframe into both columns









plt.show()