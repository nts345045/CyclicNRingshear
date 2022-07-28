import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def calc_hysteresis(df,fldX,fldY):
	"""
	Trapezoidal method for estimating area inside a discretized loop
	Calculate hysteresis based on equation 2 in Stevens (2022) Chapter 3
	"""
	# Filter out non-finite data
	IND = np.isfinite(df[fldX].values) & np.isfinite(df[fldY].values)
	Xdata = df[fldX][IND].values
	Ydata = df[fldY][IND].values
	# Get area inside the curve
	A = arb_area(Xdata,Ydata)
	# Get curve centroid
	Xbar = np.nanmean(Xdata)
	Ybar = np.nanmean(Ydata)
	# Get maximum radius
	rx = np.nanmax(np.abs(Xdata - Xbar))
	ry = np.nanmax(np.abs(Ydata - Ybar))
	rcirc = np.sqrt(rx**2 + ry**2)
	# Get area of maximum circle 
	Acirc = np.pi*rcirc**2
	# Normalize measured area by maximum circle area
	H = np.abs(A/Acirc)
	return H


issave = True
DPI = 300
FMT = 'PNG'

ROOT = os.path.join('..','..','..')
DDIR = os.path.join(ROOT,'processed_data','STEP_1')
# Normal and shear stresses
T06_NT = os.path.join(DDIR,'T06_NT_10min_x3_smoothed.csv')
T24_NT = os.path.join(DDIR,'T24_NT_10min_x3_smoothed.csv')
T96_NT = os.path.join(DDIR,'T96_NT_10min_x3_smoothed.csv')
# Processed LVDT Data
T06_LV = os.path.join(DDIR,'T06_LVDT_10min_x3_smoothed_QSS_reduced.csv')
T24_LV = os.path.join(DDIR,'T24_LVDT_10min_x3_smoothed_QSS_reduced.csv')
T96_LV = os.path.join(DDIR,'T96_LVDT_10min_x3_smoothed_min_reduced.csv')
# Results Directory
HDIR = os.path.join(ROOT,'processed_data','STEP_2')
# Selected data for Hysteresis
H6_DAT = os.path.join(HDIR,'T6_hysteresis_subset.csv')
H24_DAT = os.path.join(HDIR,'T24_hysteresis_subset.csv')
H96_DAT = os.path.join(HDIR,'T96_hysteresis_subset.csv')

ODIR = os.path.join(ROOT,'results','figures','main')


# Load data into merged data frames
df_24 = pd.concat([pd.read_csv(T24_NT,parse_dates=True,index_col=[0]),\
				   pd.read_csv(T24_LV,parse_dates=True,index_col=[0])],\
				   axis=1,ignore_index=False)
df_06 = pd.concat([pd.read_csv(T06_NT,parse_dates=True,index_col=[0]),\
				   pd.read_csv(T06_LV,parse_dates=True,index_col=[0])],\
				   axis=1,ignore_index=False)
df_96 = pd.concat([pd.read_csv(T96_NT,parse_dates=True,index_col=[0]),\
				   pd.read_csv(T96_LV,parse_dates=True,index_col=[0])],\
				   axis=1,ignore_index=False)


# Load Hysteresis defining cycle data and downsample to match LVDT sampling
H06 = pd.read_csv(H6_DAT,parse_dates=True,index_col=[0]).resample(pd.Timedelta(15,unit='sec')).median()
H24 = pd.read_csv(H24_DAT,parse_dates=True,index_col=[0]).resample(pd.Timedelta(15,unit='sec')).median()
H96 = pd.read_csv(H96_DAT,parse_dates=True,index_col=[0]).resample(pd.Timedelta(15,unit='sec')).median()

# Set up period indexing
H06p = {0:H06[H06.index == H06.index[0]],\
	    1:H06[H06['N_kPa'] == H06['N_kPa'].max()],\
 	    2:H06[H06['LVDT_mm red'] == H06['LVDT_mm red'].min()],\
	    3:H06[H06['N_kPa'] == H06['N_kPa'].min()],\
		4:H06[H06.index == H06.index[0]]}

H24p = {0:H24[H24.index == H24.index[0]],\
	    1:H24[H24['N_kPa'] == H24['N_kPa'].max()],\
 	    2:H24[H24['LVDT_mm red'] == H24['LVDT_mm red'].min()],\
	    3:H24[H24['N_kPa'] == H24['N_kPa'].min()],\
		4:H24[H24.index == H24.index[0]]}

H96p = {1:H96[H96.index == H96.index[0]],\
	    2:H96[H96['T_kPa'] == H96['T_kPa'].max()],\
 	    3:H96[H96['N_kPa'] == H96['N_kPa'].min()],\
	    4:H96[H96['LVDT_mm red'] == H96['LVDT_mm red'].max()],\
		0:H96[H96.index == H96.index[-1]]}




# Trim off lead & lags
DT_MASK = pd.Timedelta(6,unit='hour')
IND06 = (df_06.index >= df_06.index.min() + DT_MASK + pd.Timedelta(0.1,unit='hour')) &\
	 	(df_06.index <= df_06.index.max() - DT_MASK - pd.Timedelta(1,unit='hour'))
IND24 = (df_24.index >= df_24.index.min() + DT_MASK + pd.Timedelta(0.2,unit='hour')) &\
		(df_24.index <= df_24.index.max() - DT_MASK - pd.Timedelta(0.2,unit='hour'))
IND96 = (df_96.index >= df_96.index.min() + DT_MASK + pd.Timedelta(0.2,unit='hour')) &\
		(df_96.index <= df_96.index.max() - DT_MASK - pd.Timedelta(0.2,unit='hour'))

# Downsample Data to Match LVDT S/R
df_06r = df_06[IND06].copy().resample(pd.Timedelta(15,unit='sec')).median()
#interpolate(method='quadratic',limit=15,limit_direction='both')
df_24r = df_24[IND24].copy().resample(pd.Timedelta(15,unit='sec')).median()
#interpolate(method='quadratic',limit=15,limit_direction='both')
df_96r = df_96[IND96].copy().resample(pd.Timedelta(15,unit='sec')).median()
#interpolate(method='quadratic',limit=15,limit_direction='both')


# Get Cycle Boundaries
T24o = df_24r.index.min()
dt24 = pd.Timedelta(1,unit='day')
T6o = df_06r.index.min()
dt6 = pd.Timedelta(6,unit='hour')
T24_bounds = [T24o]
T6_bounds = [T6o]
for i_ in range(4):
	T24_bounds.append(T24o + (i_+1)*dt24)
	T6_bounds.append(T6o + (i_+1)*dt6)
T24_bounds.append(df_24r.index.max())
T6_bounds.append(df_06r.index.max())

T96_bounds = [df_96r.index.min(),\
			  pd.Timestamp("2021-11-05T14:32"),\
			  pd.Timestamp("2021-11-05T14:32") + pd.Timedelta(4,unit='day'),\
			  pd.Timestamp("2021-11-13T22:30"),\
			  pd.Timestamp("2021-11-13T22:30") + pd.Timedelta(4,unit='day'),\
			  df_96r.index.max()]


### DO HYSTERESIS CALCULATION ON SELECT CYCLES
# Get field names
C_ = df_06r.columns
# Define combinations to assess
XSET = ([C_[0],C_[1]],[C_[0],C_[2]],[C_[2],C_[1]])
H06_hat = [] ; H24_hat = []; H96_hat = []
for i_ in range(3):
	H06_hat.append(calc_hysteresis(H06,XSET[i_][0],XSET[i_][1]))
	H24_hat.append(calc_hysteresis(H24,XSET[i_][0],XSET[i_][1]))
	H96_hat.append(calc_hysteresis(H96,XSET[i_][0],XSET[i_][1]))
H06_hat = np.array(H06_hat)
H24_hat = np.array(H24_hat)
H96_hat = np.array(H96_hat)




#### PLOTTING ####
ncol,nrow = 3,3#,2
fig, axs = plt.subplots(ncols=ncol,nrows=nrow,figsize=(9.5,10))
# fig = plt.figure(figsize)
cmap = 'viridis'
cmaps = ['Blues_r','Purples_r','RdPu_r','Oranges_r','Greys_r']
chs = []
C_ = df_06r.columns
# XSET = ([C_[0],C_[1]],[C_[0],C_[2]],[C_[1],C_[2]])
DSET = [df_06r,df_24r,df_96r]
HD_SET = [H06,H24,H96]
PSET = [H06p,H24p,H96p]
HSET = [H06_hat,H24_hat,H96_hat]
TSET = [T6_bounds,T24_bounds,T96_bounds]
TT = [6,24,96]


labels = {'LVDT_mm red':'Relative ice-bed separation (mm)',\
		  'N_kPa':'Effective pressure (kPa)',\
		  'T_kPa':'Shear stress (kPa)'}
LBL = ['a','b','c','d']
MKR = ['ks','ko','kv','k^']
SLBL = [['a','b','c'],['d','e','f'],['g','h','i']]
for i_ in range(ncol):
	DATA = DSET[i_]
	for j_ in range(nrow):
		S_X = DATA[XSET[j_][0]]
		S_X = S_X[(S_X.index >= TSET[i_][0]) & (S_X.index <= TSET[i_][-1])]
		S_Y = DATA[XSET[j_][1]]
		S_Y = S_Y[(S_Y.index >= TSET[i_][0]) & (S_Y.index <= TSET[i_][-1])]
		# Highlight hysteresis data
		axs[j_,i_].scatter(HD_SET[i_][XSET[j_][0]],HD_SET[i_][XSET[j_][1]],s=16,color='yellow')
		
		# Place labels for phases a--d from Fig 4
		Hp = PSET[i_]
		for p_ in range(4):
			# Plot Interfaces
			axs[j_,i_].plot(Hp[p_][XSET[j_][0]],Hp[p_][XSET[j_][1]],MKR[p_])

		for k_ in range(5):
			# Filter data by time
			XIND = (S_X.index >= TSET[i_][k_]) & (S_X.index <= TSET[i_][k_+1])
			YIND = (S_Y.index >= TSET[i_][k_]) & (S_Y.index <= TSET[i_][k_+1])
			# Extract relevant data
			XD = S_X[XIND]
			YD = S_Y[YIND]
			# Scale color ramp to period length
			CD = (0.66*(XD.index - TSET[i_][k_])/pd.Timedelta(TT[i_],unit='hour'))
			# PLOT!
			cbl = axs[j_,i_].scatter(XD,YD,c=CD,cmap=cmaps[k_],s=1,vmin=0,vmax=1)
			if i_ == 0 & j_ == 0:
				chs.append(cbl)

		# Label all x-axes
		axs[j_,i_].set_xlabel(labels[XSET[j_][0]])
		# Insert hysteresis measurements
		if j_ == 0:
			# Hysteresis measure
			axs[j_,i_].text(0.95*(S_X.max() - S_X.min()) + S_X.min(),\
					0.025*(S_Y.max() - S_Y.min()) + S_Y.min(),\
					'H: %.3f'%(np.abs(HSET[i_][j_])),ha='right')
			# Subplot label
			axs[j_,i_].text(0.05*(S_X.max() - S_X.min()) + S_X.min(),\
					0.95*(S_Y.max() - S_Y.min()) + S_Y.min(),\
					SLBL[j_][i_],ha='left',va='center',fontweight='extra bold',fontstyle='italic',fontsize=16)
		elif j_ == 2 and i_ == 2:
			axs[j_,i_].text(0,120,\
					'H: %.3f'%(np.abs(HSET[i_][j_])),ha='left')
						# Subplot label
			axs[j_,i_].text(0.90*(S_X.max() - S_X.min()) + S_X.min(),\
					0.95*(S_Y.max() - S_Y.min()) + S_Y.min(),\
					SLBL[j_][i_],ha='left',va='center',fontweight='extra bold',fontstyle='italic',fontsize=16)
		else:
			axs[j_,i_].text(0.05*(S_X.max() - S_X.min()) + S_X.min(),\
					0.025*(S_Y.max() - S_Y.min()) + S_Y.min(),\
					'H: %.3f'%(np.abs(HSET[i_][j_])),ha='left')
			# Subplot label
			axs[j_,i_].text(0.90*(S_X.max() - S_X.min()) + S_X.min(),\
					0.95*(S_Y.max() - S_Y.min()) + S_Y.min(),\
					SLBL[j_][i_],ha='left',va='center',fontweight='extra bold',fontstyle='italic',fontsize=16)
		if i_ == 0:
			axs[j_,i_].set_ylabel(labels[XSET[j_][1]])

for k_ in range(5):
	cax = fig.add_axes([.925,.1 + (.8/5)*k_,.025,.8/5])
	chb = plt.colorbar(chs[k_],cax=cax,orientation='vertical',ticks=[0])
	if k_ == 2:
		cax.text(2,0.5,'Cycle number (dt/T)',rotation=270,ha='center',va='center')
	chb.ax.set_yticklabels(str(k_+1))

sax = fig.add_axes([.333,.9,.333,.025])
for i_ in range(4):
	sax.plot(i_,0,MKR[i_])
	sax.text(i_+0.5,0,LBL[i_],ha='center',va='center',fontsize=12,fontstyle='italic')
sax.plot(4,0,MKR[0])
sax.xaxis.set_visible(False)
sax.yaxis.set_visible(False)

if issave:
	SAVE_FILE = 'CNSB_Fig8_Hysteresis_v9_%ddpi.%s'%(DPI,FMT.lower())
	plt.savefig(os.path.join(ODIR,SAVE_FILE),dpi=DPI,format=FMT.lower())

plt.show()



