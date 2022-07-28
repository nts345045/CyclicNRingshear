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
include_mu = False

ROOT = os.path.join('..','..','..')
DDIR = os.path.join(ROOT,'processed_data','STEP_1')
T96_NT = os.path.join(DDIR,'T96_NT_10min_x3_smoothed.csv')
T96_LV = os.path.join(DDIR,'T96_LVDT_10min_x3_smoothed_min_reduced.csv')
FDIR = os.path.join(ROOT,'results','figures','main')
ODIR = os.path.join(ROOT,'processed_data','STEP_2')

# Load data into merged data frame
df_96 = pd.concat([pd.read_csv(T96_NT,parse_dates=True,index_col=[0]),\
				   pd.read_csv(T96_LV,parse_dates=True,index_col=[0])],\
				   axis=1,ignore_index=False)

print('data loaded')
# Trim off lead & lags
DT_MASK = pd.Timedelta(3,unit='hour')
IND96 = (df_96.index >= df_96.index.min() + DT_MASK) & (df_96.index <= df_96.index.max() - DT_MASK)

IND96x = (df_96.index >= df_96.index.min() + DT_MASK*2) & (df_96.index <= df_96.index.max() - DT_MASK*2)


# Upsample LVDT Data
df_96r = df_96[IND96].copy().interpolate(method='quadratic',limit=15,limit_direction='both')

df_96x = df_96[IND96x].copy().interpolate(method='quadratic',limit=15,limit_direction='both')
print('data upsampled')



# fig = plt.figure(figsize=(8,4))
nrows,ncols = 3,1
fig,axs = plt.subplots(ncols=ncols,nrows=nrows,figsize=(8,10))
# kick-in shared axes
axs[0].get_shared_x_axes().join(axs[0],*axs[1:])

ylims = {'N_kPa':[180,520],'LVDT_mm red':[-1,19],'T_kPa':[110,210]}
yticks = {'N_kPa':np.arange(200,550,50),'LVDT_mm red':np.arange(0,19,6),'T_kPa':np.arange(120,210,20)}
labels = {'LVDT_mm red':'Relative ice-bed separation\n(mm)',\
		  'N_kPa':'Effective pressure (kPa)',\
		  'T_kPa':'Shear stress (kPa)'}
TLAB = [96]

T0_96 = df_96r.index.min() + DT_MASK

### ACTUAL PLOTTING ###

for i_,C_ in enumerate(df_96r.columns):

	axs[i_].plot((df_96r.index - T0_96).total_seconds()/3600,df_96r[C_].values,'k')
	axs[i_].set_ylabel(labels[C_])
	axs[i_].set_xlim([(df_96r.index.min() - T0_96).total_seconds()/3600,\
						(df_96r.index.max() - T0_96).total_seconds()/3600])
	axs[i_].set_ylim(ylims[C_])
	axs[i_].set_yticks(yticks[C_])

if include_mu:
	axs[1].plot((df_96r.index - T0_96).total_seconds()/3600,100*(df_96r['T_kPa'].values/df_96r['N_kPa'].values)+100,'m')


# Overlay minima and maxima
for k_, C_ in enumerate(df_96r.columns):
	ext96 = dtiu.pick_extrema_indices(df_96x[C_],T=pd.Timedelta(36,unit='hour'))
	N_min_ind = [0,1,2,3,4,5,6,7,8,9,10]
	N_max_ind = N_min_ind
	Z_min_ind = N_min_ind
	Z_max_ind = N_min_ind

	N_min_ind = [1,6,9]
	N_max_ind = [0,2,3,5,10]
	Z_min_ind = [4,5]
	Z_max_ind = [0,4]
	if C_ == 'N_kPa':
		for i_,t_ in enumerate(ext96['I_min']):
			if i_ in N_min_ind:
				axs[k_].plot(np.ones(2,)*(t_ - T0_96).total_seconds()/3600,ylims[C_],'b--',alpha=0.25)
				axs[1].plot(np.ones(2,)*(t_ - T0_96).total_seconds()/3600,ylims['T_kPa'],'b--',alpha=0.25)
		for i_,t_ in enumerate(ext96['I_max']):
			if i_ in N_max_ind:
				axs[k_].plot(np.ones(2,)*(t_ - T0_96).total_seconds()/3600,ylims[C_],'b-',alpha=0.25)
				axs[1].plot(np.ones(2,)*(t_ - T0_96).total_seconds()/3600,ylims['T_kPa'],'b-',alpha=0.25)
##### SHADE IN SELECTED DATA FOR FIGURE 8 #####
		for i_ in range(3):
			axs[i_].fill_between([(ext96['I_max'][N_max_ind[1]] - T0_96).total_seconds()/3600,\
				 				  (ext96['I_max'][N_max_ind[2]] - T0_96).total_seconds()/3600],\
								 np.ones(2,)*ylims[df_96r.columns[i_]][0],\
								 np.ones(2,)*ylims[df_96r.columns[i_]][1],\
								 color='yellow',alpha=0.25)

		# Kick over data for hysteresis calculation
		if k_ == 0:
			df_96r_out = df_96r[(df_96r.index >= ext96['I_max'][2]) & (df_96r.index <= ext96['I_max'][3])]
			df_96r_out.to_csv(os.path.join(ODIR,'T96_hysteresis_subset.csv'),header=True,index=True)
	elif C_ == 'LVDT_mm red':
		for i_,t_ in enumerate(ext96['I_min']):
			if i_ in Z_min_ind:
				axs[k_].plot(np.ones(2,)*(t_ - T0_96).total_seconds()/3600,ylims[C_],'r--',alpha=0.25)
				axs[1].plot(np.ones(2,)*(t_ - T0_96).total_seconds()/3600,ylims['T_kPa'],'r--',alpha=0.25)
		for i_,t_ in enumerate(ext96['I_max']):
			if i_ in Z_max_ind:
				axs[k_].plot(np.ones(2,)*(t_ - T0_96).total_seconds()/3600,ylims[C_],'r-',alpha=0.25)
				axs[1].plot(np.ones(2,)*(t_ - T0_96).total_seconds()/3600,ylims['T_kPa'],'r-',alpha=0.25)

	# elif C_ == 'T_kPa':
	# 	axs[k_].plot((ext96['I_min'] - T0_96).total_seconds()/3600,ext96['V_min'],'ko')
	# 	axs[k_].plot((ext96['I_max'] - T0_96).total_seconds()/3600,ext96['V_max'],'ks')


T0v = [T0_96]
for k_ in range(1):
	# axs[0].set_title('%d hour cycling'%(TLAB[k_]))
	axs[2].set_xlabel('Elapsed time since %s\n(hours)'%(T0v[k_]))
plt.subplots_adjust(wspace=0,hspace=0)




for i_,l_ in enumerate(['a','b','c']):
	C_ = df_96r.columns[i_]
	# if i_%2 == 0:
	# 	pct = 0.1
	# else:
	pct = 0.925
	yloc = (ylims[C_][1] - ylims[C_][0])*pct + ylims[C_][0]
	axs[i_].text(5,yloc,l_,fontweight='extra bold',fontstyle='italic',ha='center',va='center',fontsize=16)


if issave:
	SAVE_FILE = 'CNSB_Fig7_T96_TimeSeries_v9_%ddpi.%s'%(DPI,FMT.lower())
	plt.savefig(os.path.join(FDIR,SAVE_FILE),dpi=DPI,format=FMT.lower())

#### TODO ####
# Get minimum and maxium cavity geometry timings from videos into row 3
# Outline of Fig. 2 timeframe into both columns









plt.show()