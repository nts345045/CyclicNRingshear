import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from matplotlib.gridspec import GridSpec
sys.path.append('../..')
import util.datetimeindex as dtiu

issave = True
DPI = 300
FMT = 'PNG'
include_mu = False

ROOT = os.path.join('..','..','..','..')
DDIR = os.path.join(ROOT,'processed_data','STEP_1')
T24_NT = os.path.join(DDIR,'T24_NT_10min_x3_smoothed.csv')
T24_LV = os.path.join(DDIR,'T24_LVDT_10min_x3_smoothed_QSS_reduced.csv')
T24_S = os.path.join(DDIR,'Postprocessed_Contact_Areas.csv')
FDIR = os.path.join(ROOT,'results','figures','main')
ODIR = os.path.join(ROOT,'processed_data','STEP_2')

# Load data into merged data frame
df_24 = pd.concat([pd.read_csv(T24_NT,parse_dates=True,index_col=[0]),\
				   pd.read_csv(T24_LV,parse_dates=True,index_col=[0])],\
				  axis=1,ignore_index=False)
df_S24 = pd.read_csv(T24_S,parse_dates=True,index_col=[0])

# Trim off lead & lags
DT_MASK = pd.Timedelta(3,unit='hour')
IND24 = (df_24.index >= df_24.index.min() + DT_MASK) & (df_24.index <= df_24.index.max() - DT_MASK)

IND24x = (df_24.index >= df_24.index.min() + DT_MASK*2) & (df_24.index <= df_24.index.max() - DT_MASK*2)


# Upsample LVDT Data
df_24r = df_24[IND24].copy().interpolate(method='quadratic',limit=15,limit_direction='both')
df_24x = df_24[IND24x].copy().interpolate(method='quadratic',limit=15,limit_direction='both')

##########################
#### PLOTTING SECTION ####
##########################

# fig = plt.figure(figsize=(8,4))
nrows,ncols = 4,1
fig,axs = plt.subplots(ncols=ncols,nrows=nrows,figsize=(8,10))
# kick-in shared axes
axs[0].get_shared_x_axes().join(axs[0],*axs[1:])
# axs[0,1].get_shared_x_axes().join(axs[0,1],*axs[1:,1])

ylims = {'N_kPa':[180,520],'LVDT_mm red':[-1,7],'T_kPa':[110,210]}
yticks = {'N_kPa':np.arange(200,550,50),'LVDT_mm red':np.arange(0,7,1),'T_kPa':np.arange(120,210,20)}
labels = {'LVDT_mm red':'Relative ice-bed\nseparation(mm)',\
		  'N_kPa':'Effective pressure (kPa)',\
		  'T_kPa':'Shear stress (kPa)'}
TLAB = [24]

T0_24 = df_24r.index.min() + DT_MASK

### ACTUAL PLOTTING ###
for i_,C_ in enumerate(df_24r.columns):
	# Plot data
	axs[i_].plot((df_24r.index - T0_24).total_seconds()/3600,df_24r[C_].values,'k',zorder=1,linewidth=2)
	axs[i_].set_ylabel(labels[C_])
	axs[i_].set_xlim([(df_24r.index.min() - T0_24).total_seconds()/3600,\
						(df_24r.index.max() - T0_24).total_seconds()/3600])
	axs[i_].set_ylim(ylims[C_])
	axs[i_].set_yticks(yticks[C_])

### PLOT CONTACT AREA STUFF ###
lbda = 0.3*2*np.pi/4 # [m] Bed wavelength for outer edge of experimental chamber

S_Packs = [('Lee',df_S24['Lee X Piecewise']/lbda,'blue'),\
		   ('Stoss',df_S24['Stoss X Piecewise']/lbda,'darkorange'),\
		   ('Total',(df_S24['Lee X Piecewise'] + df_S24['Stoss X Piecewise'])/lbda,'firebrick')]
axb = axs[3]
for i_,P_ in enumerate(S_Packs):
	l_,S_,c_ = P_

	axb.plot((S_.index - T0_24).total_seconds()/3600,(S_ - S_.min())/(S_.max() - S_.min()),\
			 '.-',color=c_,label=l_,zorder=1000,alpha=0.66)
	axb.set_ylabel('Normalized\ncontact fraction')
	axb.legend(ncol=3,loc='upper right')


axb.set_ylim([-0.1,1.25])
# axc.set_ylim([-0.1,1.25])
axb.set_yticks(np.linspace(0,1,6))
# axc.set_yticks(np.linspace(0,1,6))



# Overlay minima and maxima
for k_, C_ in enumerate(df_24r.columns):
	ext24 = dtiu.pick_extrema_indices(df_24x[C_],T=pd.Timedelta(23,unit='hour'))
	if C_ == 'N_kPa':
		for t_ in ext24['I_min']:
			axs[k_].plot(np.ones(2,)*(t_ - T0_24).total_seconds()/3600,ylims[C_],'b--',alpha=0.25)
			axs[1].plot(np.ones(2,)*(t_ - T0_24).total_seconds()/3600,ylims['T_kPa'],'b--',alpha=0.25)
			axs[3].plot(np.ones(2,)*(t_ - T0_24).total_seconds()/3600,[-.1,1.25],'b--',alpha=0.25)

		for t_ in ext24['I_max']:
			axs[k_].plot(np.ones(2,)*(t_ - T0_24).total_seconds()/3600,ylims[C_],'b-',alpha=0.25)
			axs[1].plot(np.ones(2,)*(t_ - T0_24).total_seconds()/3600,ylims['T_kPa'],'b-',alpha=0.25)
			axs[3].plot(np.ones(2,)*(t_ - T0_24).total_seconds()/3600,[-.1,1.25],'b-',alpha=0.25)

	elif C_ == 'LVDT_mm red':
		for t_ in ext24['I_min']:
			axs[k_].plot(np.ones(2,)*(t_ - T0_24).total_seconds()/3600,ylims[C_],'r--',alpha=0.25)
			axs[1].plot(np.ones(2,)*(t_ - T0_24).total_seconds()/3600,ylims['T_kPa'],'r--',alpha=0.25)
			axs[3].plot(np.ones(2,)*(t_ - T0_24).total_seconds()/3600,[-.1,1.25],'r--',alpha=0.25)
		for t_ in ext24['I_max']:
			axs[k_].plot(np.ones(2,)*(t_ - T0_24).total_seconds()/3600,ylims[C_],'r-',alpha=0.25)
			axs[1].plot(np.ones(2,)*(t_ - T0_24).total_seconds()/3600,ylims['T_kPa'],'r-',alpha=0.25)
			axs[3].plot(np.ones(2,)*(t_ - T0_24).total_seconds()/3600,[-.1,1.25],'r-',alpha=0.25)

#### SHADE IN SELECTED DATA FOR FIGURE 8 #####
		for i_ in range(3):
			axs[i_].fill_between([(ext24['I_max'][4] - T0_24).total_seconds()/3600,\
				 				  (ext24['I_max'][5] - T0_24).total_seconds()/3600],\
								 np.ones(2,)*ylims[df_24r.columns[i_]][0],\
								 np.ones(2,)*ylims[df_24r.columns[i_]][1],\
								 color='yellow',alpha=0.25)		
		if k_ == 0:
			df_24r_out = df_24r[(df_24r.index >= ext24['I_max'][4]) & (df_24r.index <= ext24['I_max'][5])]
			df_24r_out.to_csv(os.path.join(ODIR,'T24_hysteresis_subset.csv'),header=True,index=True)



T0v = [T0_24]
for k_ in range(1):
	# axs[0].set_title('%d hour cycling'%(TLAB[k_]))
	axs[3].set_xlabel('Elapsed time since %s\n(hours)'%(T0v[k_]))
# for k_ in range(4):
# 	axs[k_].text(120,)
for i_,l_ in enumerate(['a','b','c']):
	C_ = df_24r.columns[i_]
	if i_%2 == 0:
		pct = 0.1
	else:
		pct = 0.90
	yloc = (ylims[C_][1] - ylims[C_][0])*pct + ylims[C_][0]
	axs[i_].text(0,yloc,l_,fontweight='extra bold',fontstyle='italic',ha='center',va='center',fontsize=16)

axs[3].text(0,.9*1.25,'d',fontweight='extra bold',fontstyle='italic',ha='center',va='center',fontsize=16)
# axs[3].set_ylim([-.1,1.1])
plt.subplots_adjust(wspace=0,hspace=0)


if issave:
	SAVE_FILE = 'CNSB_Fig3_T24_TimeSeries_diss_v2_%ddpi.%s'%(DPI,FMT.lower())
	plt.savefig(os.path.join(FDIR,SAVE_FILE),dpi=DPI,format=FMT.lower())


### PLOT CONTACT AREA STUFF (FIGURE 4) ###
lbda = 0.3*2*np.pi/4 # [m] Bed wavelength for outer edge of experimental chamber
nrows,ncols = 3,1
fig,axs = plt.subplots(ncols=ncols,nrows=nrows,figsize=(8,10))
SPL = ['a','b','c']
S_Packs = [('Lee',df_S24['Lee X Piecewise']/lbda,'blue'),\
		   ('Stoss',df_S24['Stoss X Piecewise']/lbda,'darkorange'),\
		   ('Total',(df_S24['Lee X Piecewise'] + df_S24['Stoss X Piecewise'])/lbda,'firebrick')]

for i_,P_ in enumerate(S_Packs):
	l_,S_,c_ = P_
	axs[i_].plot((S_.index - T0_24).total_seconds()/3600,S_ ,\
			 '.-',color=c_,label=l_,zorder=1000,alpha=1)
	axs[i_].set_ylabel('%s contact fraction'%(l_))
	if i_ == 0:
		pct = 0.1
	else:
		pct = 0.90
	yloc = (S_.max() - S_.min())*pct + S_.min()
	axs[i_].text(0,yloc,SPL[i_],fontweight='extra bold',fontstyle='italic',ha='center',va='center',fontsize=16)


	# axs[i_].legend(ncol=3,loc='upper right')
	axb = axs[i_].twinx()
	axb.plot((S_.index - T0_24).total_seconds()/3600,S_*lbda*1e3,\
			 '.-',color=c_,zorder=1000,alpha=0.66)
	axb.set_ylabel('%s contact length (mm)'%(l_),rotation=270,labelpad=15)
	if i_ < 2:
		axs[i_].xaxis.set_visible(False)
	else:
		axs[2].set_xlabel('Elapsed time since %s\n(hours)'%(T0v[k_]))

	# axs[i_].set_xlim([])


plt.subplots_adjust(wspace=0,hspace=0)

# axs[i_].set_ylim([-0.1,1.25])
# axc.set_ylim([-0.1,1.25])
# axs[i_].set_yticks(np.linspace(0,1,6))
# axc.set_yticks(np.linspace(0,1,6))


if issave:
	SAVE_FILE = 'CNSB_Fig4_T24_ContactAreas_diss_v2_%ddpi.%s'%(DPI,FMT.lower())
	plt.savefig(os.path.join(FDIR,SAVE_FILE),dpi=DPI,format=FMT.lower())









plt.show()