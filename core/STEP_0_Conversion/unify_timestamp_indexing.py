import pandas as pd
import sys
import os
from tqdm import tqdm
"""
This module shifts all data to local Wisconsin time for October/November
and saves the CSV data with strings formatted for quicker read/write using
the Pandas library.

"""


#### FILE PATHS & DATA LOADING WRAPPERS ####
ROOT = os.path.join('..','..','..')
# Raw Files
D_LVDT_UTC = os.path.join(ROOT,'data','physical_units','LVDT','UTC_LVDT_data_PD_thresh0.05mm.csv')
D_NTP_UTC = os.path.join(ROOT,'data','physical_units','Pressures_Stresses','UTC_N_T_P.csv')

ODIR = os.path.join(ROOT,'processed_data','timeseries','STEP0')


### LOAD DATA ###
DTL = pd.Timedelta(-5,unit='hour')
df_NTP = pd.read_csv(D_NTP_UTC)
# Reformat indices (speeds up read compared to calling parse_dates=True in the above line)
IDX = []
for i_ in tqdm(range(len(df_NTP))):
	IDX.append(pd.Timestamp(df_NTP['TimeUTC'].values[i_]) + DTL)
df_NTP.index = IDX
df_NTP = df_NTP[['N_kPa','T_kPa','Pw1_kPa','Pw2_kPa']]
print('Stress data loaded')
print('Writing in Pandas Friendly Format (drastically accelerates I/O with pd.read_csv)')
df_NTP.to_csv(os.path.join(ODIR,'NTP_Local_PD_DateTime.csv'),header=True,index=True)

df_LVDT = pd.read_csv(D_LVDT_UTC,parse_dates=True,index_col=[0])
df_LVDT_out = df_LVDT['LVDT_mm_stitched']
df_LVDT_out.index += DTL
df_LVDT_out.to_csv(os.path.join(ODIR,'LVDT_Local_PD_DateTime_FULL.csv'),header=True,index=True)



