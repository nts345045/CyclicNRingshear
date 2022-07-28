import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob

issave = True
DPI = 300
FMT = 'PNG'

ROOT = os.path.join('..','..','..')
DDIR = os.path.join(ROOT,'processed_data','STEP_1')
D_NT = os.path.join(DDIR,'FULL_NT_10min_x3_smoothed.csv')
ODIR = os.path.join(ROOT,'results','figures','main')

df = pd.read_csv(D_NT,parse_dates=True,index_col=[0])


LBV = ['Steady State','T24','T6','T96']
TSV = [pd.Timestamp("2021-10-25T13:00:00"),pd.Timestamp('2021-10-26T13:56'),\
	   pd.Timestamp('2021-11-1T11:09:15'),pd.Timestamp("2021-11-03T15:52")]
TEV = [pd.Timestamp("2021-10-26T13:56"),pd.Timestamp("2021-10-31T14:03"),\
	   pd.Timestamp("2021-11-02T17:26:30"),pd.Timestamp("2021-11-22T12:09")]

xlims = [pd.Timestamp('2021-10-25'),TEV[-1] + pd.Timedelta(3,unit='hour')]
fig = plt.figure(figsize=(12,5.3))
ax1 = fig.add_subplot(111)

ax1.plot(df['N_kPa'])

ylims = ax1.get_ylim()

# Insert experiment labels
for i_,l_ in enumerate(LBV):
	ax1.plot([TSV[i_],TSV[i_]],ylims,'k-.')
	ax1.plot([TEV[i_],TEV[i_]],ylims,'k-.')
	if i_ == 0:
		ax1.text(TSV[i_] + (TEV[i_] - TSV[i_])/2,348,LBV[i_],ha='center',va='center',rotation=90)
	else:
		ax1.text(TSV[i_] + (TEV[i_] - TSV[i_])/2,180,LBV[i_],ha='center')

# Insert cycle numbering
T96C = [pd.Timestamp("2021-11-04T12"),pd.Timestamp("2021-11-06T12"),pd.Timestamp("2021-11-10T12"),\
		pd.Timestamp("2021-11-14T20"),pd.Timestamp("2021-11-19T13")]

for i_ in range(5):
	ax1.text(TSV[1] + pd.Timedelta(6,unit='hour') + i_*pd.Timedelta(24,unit='hour'),510,i_+1,ha='center')
	# ax1.text(TSV[2] + pd.Timedelta())
	ax1.text(T96C[i_],510,i_ + 1, ha='center')
	# ax1.text(TSV[2] + pd.Timedelta(3,unit='hour') + i_*pd.Timedelta(6,unit='hour'),\
	# 		 190 + ((i_+1)%2)*310,i_ + 1, ha='center')
	ax1.text(TSV[2] + pd.Timedelta(3,unit='hour') + i_*pd.Timedelta(6,unit='hour'),510,i_ + 1, ha='center')
# ax1.text(pd.Timestamp("2021-10-21T20"),510,'Cycle Number',ha='left')
# ax1.text(pd.Timestamp("2021-10-21T20"),180,'Experiment',ha='left')

ax1.set_xlim(xlims)
ax1.set_ylim([170,530])
ax1.set_ylabel('Effective pressure (kPa)')
ax1.set_xlabel('Local date/time (UTC - 5)')


if issave:
	OFILE = os.path.join(ODIR,'CNSB_Fig2_Effective_Pressure_Labeled_v9_%ddpi.%s'%(DPI,FMT.lower()))
	plt.savefig(OFILE,dpi=DPI,format=FMT.lower())

plt.show()