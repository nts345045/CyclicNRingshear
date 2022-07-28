import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
sys.path.append(os.path.join('..','..'))
import core.model.steadystate_model as ssm

ROOT = os.path.join('..','..','..','..')
DDIR = os.path.join(ROOT,'processed_data','STEP_1')
T24_NT = os.path.join(DDIR,'T24_NT_10min_x3_smoothed.csv')
T24_LV = os.path.join(DDIR,'T24_LVDT_10min_x3_smoothed_QSS_reduced.csv')
T6_NT = os.path.join(DDIR,'T06_NT_10min_x3_smoothed.csv')
T6_LV = os.path.join(DDIR,'T06_LVDT_10min_x3_smoothed_QSS_reduced.csv')
ODIR = os.path.join(ROOT,'results','figures','main')

# LOAD EXPERIMENTAL DATA #
df_NT6 = pd.read_csv(T6_NT,parse_dates=True,index_col=[0])
df_Z6 = pd.read_csv(T6_LV,parse_dates=True,index_col=[0])
df_NT24 = pd.read_csv(T24_NT,parse_dates=True,index_col=[0])
df_Z24 = pd.read_csv(T24_LV,parse_dates=True,index_col=[0])

issave = True
DPI = 300
FMT = 'PNG'

# CALCULATE MODEL VALUES #
# Populate \tau(N) and S(N) vectors from steady-state model
model_vectors = ssm.calc_TS_vectors(100e3,600e3)

# INTERPOLATE STEADY-STATE VALUES
T_hat_24 = np.interp(df_NT24['N_kPa'].values*1e3,model_vectors['N_Pa'],model_vectors['T_Pa'])
S_hat_24 = np.interp(df_NT24['N_kPa'].values*1e3,model_vectors['N_Pa'],model_vectors['S'])

T_hat_6 = np.interp(df_NT6['N_kPa'].values*1e3,model_vectors['N_Pa'],model_vectors['T_Pa'])
S_hat_6 = np.interp(df_NT6['N_kPa'].values*1e3,model_vectors['N_Pa'],model_vectors['S'])



# NORMALIZE MODEL AND EXPERIMENTAL VALUES #
NT_hat_24 = (T_hat_24 - T_hat_24.min())/(T_hat_24.max() - T_hat_24.min())
NS_hat_24 = (S_hat_24 - S_hat_24.min())/(S_hat_24.max() - S_hat_24.min())
NT_hat_6 = (T_hat_6 - np.nanmin(T_hat_6))/(np.nanmax(T_hat_6) - np.nanmin(T_hat_6))
NS_hat_6 = (S_hat_6 - np.nanmin(S_hat_6))/(np.nanmax(S_hat_6) - np.nanmin(S_hat_6))

NT_obs_24 = (df_NT24['T_kPa'].values - df_NT24['T_kPa'].min())/(df_NT24['T_kPa'].max() - df_NT24['T_kPa'].min())
NS_obs_24 = 1 -1*(df_Z24["LVDT_mm red"].values - df_Z24["LVDT_mm red"].min())/(df_Z24["LVDT_mm red"].max() - df_Z24["LVDT_mm red"].min())
NT_obs_6 = (df_NT6['T_kPa'].values - df_NT6['T_kPa'].min())/(df_NT6['T_kPa'].max() - df_NT6['T_kPa'].min())
NS_obs_6 = 1 -1*(df_Z6["LVDT_mm red"].values - df_Z6["LVDT_mm red"].min())/(df_Z6["LVDT_mm red"].max() - df_Z6["LVDT_mm red"].min())

# PLOT #
fig,axs = plt.subplots(ncols=2,nrows=2,figsize=(8,10))

OBS = [[NT_obs_24,NT_obs_6],[NS_obs_24,NS_obs_6]]
HAT = [[NT_hat_24,NT_hat_6],[NS_hat_24,NS_hat_6]]
IND = [[(df_NT24.index - df_NT24.index.min() - pd.Timedelta(6,unit='hour')).total_seconds()/3600,\
		(df_NT6.index - df_NT6.index.min() - pd.Timedelta(6,unit='hour')).total_seconds()/3600],\
	   [(df_Z24.index - df_Z24.index.min() - pd.Timedelta(6,unit='hour')).total_seconds()/3600,\
	    (df_Z6.index - df_Z6.index.min() - pd.Timedelta(6,unit='hour')).total_seconds()/3600]]
COL = ['k','royalblue']
LAB = ['Observed','Modeled']
PAR = ['shear stress','contact size']
SPL = [['a','c'],['b','d']]

for i_ in range(2):
	for j_ in range(2):
		ijIND = IND[i_][j_]
		ijOBS = OBS[i_][j_]
		ijHAT = HAT[i_][j_]
		iIND = IND[0][j_]
		axs[i_,j_].plot(ijIND,ijOBS,color=COL[0],label=LAB[0])
		axs[i_,j_].plot(iIND,HAT[i_][j_],color=COL[1],label=LAB[1])

		if i_ == 0:
			axs[i_,j_].xaxis.set_visible(False)

		else:
			axs[i_,j_].set_xlabel('Elapsed time (hours)')
		if j_ == 1:
			axs[i_,j_].yaxis.set_visible(False)
			if i_ == 0:
				axs[i_,j_].legend(ncol=2,bbox_to_anchor=(0.4,1.10))
		else:
			axs[i_,j_].set_ylabel('Normalized %s'%(PAR[i_]))

		axs[i_,j_].set_xlim([np.nanmin(ijIND),np.nanmax(ijIND)])
		axs[i_,j_].set_ylim([-0.075,1.075])

		axs[i_,j_].text((np.nanmax(ijIND) - np.nanmin(ijIND))*0.95 + np.nanmin(ijIND),0.95,SPL[i_][j_],\
						fontweight='extra bold',fontstyle='italic',fontsize=14,\
						ha = 'right')


plt.subplots_adjust(wspace=0,hspace=0)



# SAVE FIGURE #
if issave:
	SAVE_FILE = 'CNSB_Fig9_NormModel_v9_%ddpi.%s'%(DPI,FMT.lower())
	print('Saving figure to disk at: %s'%(os.path.join(ODIR,SAVE_FILE)))
	plt.savefig(os.path.join(ODIR,SAVE_FILE),dpi=DPI,format=FMT.lower())


plt.show()