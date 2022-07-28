import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('..')
import util.datetimeindex as dtiu

def reorg_picks(df,camno=None):
	"""
	Reorganize stoss, crest, and lee pick geometries into a data-frame with 
	a DatetimeIndex'd DataFrame
	"""
	odf = pd.DataFrame()
	for C_ in ['Crest','Stoss','Lee']:
		idf = df.copy()[df['Location']==C_]
		Sx = idf['X']
		Sx.name = '%s X'%(C_)
		Sy = idf['Y']
		Sy.name = '%s Y'%(C_)
		Si = idf['id']
		Si.name = '%s id'%(C_)

		iodf = pd.concat([Si,Sx,Sy],axis=1,ignore_index=False)
		odf = pd.concat([odf,iodf],axis=1,ignore_index=False)
	if camno is not None:
		cam_idx = []
		for i_ in range(len(odf)):
			cam_idx.append(camno)
		odf = pd.concat([odf,pd.Series(cam_idx,odf.index,name='cam')],axis=1,ignore_index=False)
	return odf


	
ROOT = os.path.join('..','..','..')
IDIR = os.path.join(ROOT,'data','images','image_processing','picks')
PICKS = os.path.join(IDIR,'master_picks_Exported_Points.csv')
ODIR = os.path.join(ROOT,'processed_data','timeseries','STEP_1')
# Organize picks into organized sets w/ DateTime indexing
df_main = reorg_picks(pd.read_csv(PICKS,parse_dates=['DateTime'],index_col='DateTime'))

# Put in camera indexing
IND_c4 = (df_main.index >= pd.Timestamp("2021-10-28T12")) & (df_main.index <= pd.Timestamp("2021-10-30T13"))

cam = []
for i_ in IND_c4:
	if i_:
		cam.append('cam4')
	else:
		cam.append('cam2')
df_main = pd.concat([df_main,pd.Series(cam,index=df_main.index,name='cam')],axis=1,ignore_index=False)

# Create stoss & lee lengths
df_stoss = pd.concat([df_main[['Stoss X','Crest X','cam']],\
							  pd.Series(df_main['Crest X'] - df_main['Stoss X'],name='Stoss dx')],\
							 axis=1,ignore_index=False)
df_lee = pd.concat([df_main[['Lee X','Crest X','cam']],\
							  pd.Series(df_main['Lee X'] - df_main['Crest X'],name='Lee dx')],\
							 axis=1,ignore_index=False)

#### ADDITIONAL FINE-SCALE ADJUSTMENTS ####
"""
Observations from extended viewing of time-lapse images and camera transfers
1) Camera 4's vantage point provides a greater about of fitting surfaces to the flattened model
   so preference should be given to the general scaling of stoss contact areas for absolute values
2) The transition from cam2 to cam4 a the subsequent transitions occurs when cavity dilation is at maximum, 
   which generally stabilizes cavity geometries. As such, the values of stoss & lee contact areas should
   match across these transitions
3) Stoss contact areas are observed to be slightly larger for the first peak observed by cam2 (c. 2021-10-28) 
   compared to the second peak (c. 2021-10-31)
4) Stoss contact areas are observed to be slightly smaller for the first peak observed by cam4 (c. 2021-10-29)
   compared to the second peak (c. 2021-10-30)
5) Due to the proximity of detachment points (lee contact areas) compared to reattachement points (stoss contact
   areas) lee contact length estimates are likely less sensitive to re-projection artifacts.

"""
# Fetch stoss length values at each side of camera transfers
for i_ in range(len(IND_c4) - 1):
	if not IND_c4[i_] and IND_c4[i_ + 1]:
		loc1 = [i_, i_ + 1]
	elif IND_c4[i_] and not IND_c4[i_ + 1]:
		loc2 = [i_, i_ + 1]

S_dxs = df_stoss['Stoss dx']
L_dxs = df_lee['Lee dx']

sref0 = S_dxs.values[loc1[0]]
sref1 = S_dxs.values[loc1[1]]
sref2 = S_dxs.values[loc2[0]]
sref3 = S_dxs.values[loc2[1]]
S2_dxs = S_dxs[~IND_c4]
S4_dxs = S_dxs[IND_c4]
u4 = S4_dxs.mean()
u2b = np.mean([sref0,sref3])


lref0 = L_dxs.values[loc1[0]]
lref1 = L_dxs.values[loc1[1]]
lref2 = L_dxs.values[loc2[0]]
lref3 = L_dxs.values[loc2[1]]
L2_dxs = L_dxs[~IND_c4]
L4_dxs = L_dxs[IND_c4]




##### START WITH CAM2 CORRECTIONS BASED ON POINTS 1--3 FROM ABOVE #####
# More likely scenario because we can get more definitive features from cam4's coverage
# i.e., 1 crest, 1 trough, 1 full stoss face, and 2 partial lee faces
# camera 2 only shows a crest and halves of 2 faces

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(S_dxs,'.',label="Crest Corrected")

# Shift cam2 data to line up at first transition
# This transition corresponds with a period of max dilation when
# contact areas should be relatively low and stable
DC = sref1 - sref0 
# DC = u4 - u2b

## Apply DC correction to S2 data to match average values
S_dxsc1 = pd.concat([S4_dxs,S2_dxs + DC],axis=0).sort_index()
S_dxsc1.name = 'Stoss X DC'
ax1.plot(S_dxsc1,'o',label='Offset corrected')

##### NOW CORRECT CAM4 BASED ON POINTS 1--4 FROM ABOVE #####
# Get the slope of cam4 data - this assumes a linear time-dependent adjustment
# While somewhat nonphysical, this is a simpler model for correcting artifacts from
# the nonlinear image projection onto a focal-plane used

# Get slope of average linear fit to cam4 stoss contact lengths
mod = dtiu.fit_dateline(S4_dxs.index,S4_dxs.values)
# Define any scalar modifications here and print out model line being subtracted
# from cam4 stoss contact length observations
mod *= 1
ref_v = sref1
print('Correction model for camera 4 data: %.3e t + %.3e '%(mod,ref_v))
# Reduce data by the defined linear model
S4_dxstc = dtiu.reduce_series_by_slope(S4_dxs,mod*1.2,S4_dxs.index.min(),sref1)

# Recombine into fully corrected series
S_dxsc2 = pd.concat([S4_dxstc + ref_v,S2_dxs + DC],axis=0).sort_index()
S_dxsc2.name = 'Stoss X DC+M'
ax1.plot(S_dxsc2,'v-',label='Offset + C4 Slope Correction',alpha=0.5)


ax1b = ax1.twinx()
ax1b.plot(df_lee['Lee X'],'ro-',alpha=0.5)


ax1.set_xlabel("Date Time (UTC - 5)")
ax1.set_ylabel("Stoss contact length (m)")
ax1b.set_ylabel("Lee contact length (m)",color='red')
ax1.legend()





##### ALTERNATE APPROACH - PIECEWISE STITCHGING #####
# Hypothesis is that mismatches in edge fits arise from some adjustment between cam2 
# scene 1 and cam2 scene 2
# This seeks to preserve 

DC1 = sref0 - sref1
DC2 = sref3 - sref2

S_c2s1 = S_dxs[S_dxs.index <= S_dxs.index[loc1[0]]] - DC1
S_c2s2 = S_dxs[S_dxs.index >= S_dxs.index[loc2[1]]] - DC2

S_DC = pd.concat([S_c2s1,S4_dxs,S_c2s2],axis=0).sort_index()
S_DC.name = 'Stoss X Piecewise'
ax1.plot(S_DC,'k^',label='24-mean edge offsets')

### LEE PROCESSING ###
DC1 = lref0 - lref1
DC2 = lref3 - lref2

L_c2s1 = L_dxs[L_dxs.index <= L_dxs.index[loc1[0]]] - DC1
L_c2s2 = L_dxs[L_dxs.index >= L_dxs.index[loc2[1]]] - DC2

L_DC = pd.concat([L_c2s1,L4_dxs,L_c2s2],axis=0).sort_index()
L_DC.name = 'Lee X Piecewise'


plt.show()

#### SAVE DATA TO DISK ####

df_out = pd.concat([df_lee['Lee X'],L_DC,S_DC,S_dxs,S_dxsc1,S_dxsc2],axis=1,ignore_index=False)
# 
df_out.to_csv(os.path.join(ODIR,'Postprocessed_Contact_Areas.csv'),header=True,index=True)

