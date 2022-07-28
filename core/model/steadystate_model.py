import os
import numpy as np
from tqdm import tqdm
"""
This module calculates discretized solutions for cavity geometry and resultant
drag
"""



def acot(x):
	return np.pi/2 - np.arctan(x)

def defined_bed():
	# 	   Step height, Viscosity, Flow exp, Wlength. ,  Slide speed (m/a)
	out = {'h':0.0253*2,'B':6.3e7,'nn':3,'lbda':.31425,'US':15}
	return out

def bedmodel(lbda=.3*np.pi*2.*0.25,hh=.078,npts=5001,pi_offset=0.5,lbda_offset=0.5,ncycles=1.5):
	xv = np.linspace(0,ncycles*lbda,npts)
	# Calculate wavenumber
	kk = 2.*np.pi/lbda
	# Calculate bed amplitude
	aa = hh/2.
	# Calculate bed elevation vector (h(x) - defines a sinusoidal bed with wavelength lambda)
	# hx = (h/2.)*np.sin((2.*np.pi*xv/lbda) + np.pi/2.) + h/2.
	hx = aa*(np.cos(kk*xv - pi_offset*np.pi) + 1.)
	# Calculate get normlized locatins
	iv = xv/lbda
	df = pd.DataFrame({'x_m':xv - lbda_offset*lbda,'y_m':hx,'x_ind':iv})
	return df



def lvel2avel(r_ref=.2,lvel=15):
	C_ref = 2.*np.pi*r_ref
	vv = C_ref/lvel
	return vv

def avel2lvel(c_ref,vv):
	return c_ref*vv


def bedslice(r_ref):
	ao = 0.04; ai = 0.01375
	ro = 0.3;  ri = 0.1
	kk = 4
	mm = (ao - ai)/(ro - ri)
	a_ref = ai + (r_ref-ri)*mm 
	l_ref = 2.*np.pi*r_ref/kk
	return a_ref,l_ref




def calc_SS(NN,US=15,lbda=.31425,hh=0.0253*2,BB=6.3e7,nn=3,npts=5001):
	# Convert sliding velocity into m/sec
	Uss = US/(365*24*3600)
	# Discretize model domain
	xv = np.linspace(0,lbda,npts)
	# Calculate length-scale of the cavity in the absence of subsequent bumps (lc) (Equation A2 in Zoet & Iverson, 2015)
	lc = np.sqrt(8.*Uss*hh/np.pi*(BB/NN)**nn)
	# Calculate cavity roof elevation vector (g(x) - equation 4 in Kamb, 1987)
	gx = np.real(hh*(0.5 - (1./np.pi)*np.arcsin((2.*xv-lc)/lc) - (2.*(xv - lc) * np.sqrt(xv*(lc - xv))/(np.pi*lc**2))))
	# Calculate wavenumber
	kk = 2.*np.pi/lbda
	# Calculate bed amplitude
	aa = hh/2.
	# Calculate bed elevation vector (h(x) - defines a sinusoidal bed with wavelength lambda)
	# hx = (h/2.)*np.sin((2.*np.pi*xv/lbda) + np.pi/2.) + h/2.
	hx = aa*(np.cos(kk*xv) + 1.)
	# Calculate cavity space
	dd = gx - hx
	# breakpoint()
	# Get indices of cavity edges
	try:
		sa = np.max(np.argwhere(dd > 0))
		# Get first index of positive cavity volume
		sb = np.min(np.argwhere(dd > 0))
	# If model says no cavity, assign cavity start and end as same point
	except:
		sa = 0
		sb = 0
		# SS = 1.
	# Calculate bed contact fraction
	SS = 1. - ((xv[sa] - xv[sb])/lbda)

	return SS

def calc_from_SN(NN,SS,hh=0.0253*2,lbda=.31425):
	kk = 2.*np.pi/lbda
	# Calculate bed amplitude
	aa = hh/2.
	# Calculate wave-number scaled cavity critical length (k*x_c)
	kxc = acot((2.*np.pi*(1. - SS) + np.sin(2.*np.pi*SS))/(1. - np.cos(2.*np.pi*SS)))
	# Calculate fo factor (Coefficient of friction)
	# fo = ((h/2.)/lbda)*np.pi*N
	fo = aa*kk*0.5*NN
	# Calculate geometric factor
	PHI = ((np.pi*SS - 0.5*np.sin(2.*np.pi*SS))*np.sin(np.pi*SS - kxc))/(np.sin(np.pi*SS) - np.pi*SS*np.cos(np.pi*SS))
	# Calculate drag
	TAU = fo*PHI
	return TAU

def calc_drag(NN,hh=0.0253*2,BB=6.3e7,nn=3,lbda=0.31425,US=15,npts=5001):
	SS = calc_SS(NN,BB=BB,nn=nn,lbda=lbda,US=US,npts=npts)
	TT = calc_from_SN(NN,SS,hh=hh,lbda=lbda)
	return SS,TT


def bedslice_fit(NN,nvals=31,**kwargs):
	r_vect = np.linspace(0.1,0.3,nvals)
	U_vect = avel2lvel(2.*np.pi*r_vect,lvel2avel())
	SSo = []; TTo = []; lbdo = []; aao = []
	for i_,r_ in enumerate(r_vect):
		a_i,l_i = bedslice(r_)
		S_i,T_i = calc_drag(NN,hh=a_i*2,lbda=l_i,US=U_vect[i_],**kwargs)
		SSo.append(S_i)
		TTo.append(T_i)
		lbdo.append(l_i)
		aao.append(a_i)
	return {'r':r_vect,'Us':U_vect,'S':SSo,'T':TTo,'lbda':lbdo,'a':aao}





def model_TS_single(N,h=0.0253*2,B=6.3e7,nn=3,lbda=.31425,US=15,npts=5001,output='tau'):
	"""
	Model cavity geometry (S) from specified sinusoidal bed
	geometry (assuming 1-D obstacles) based on theory from Lliboutry (1968)

	:: INPUTS ::
	:type N: float
	:param N: effective stress in [Pa]
	:type h: float
	:param h: obstacle height in [m] (double step amplitude)
	:type B: float
	:param B: effective viscosity of ice in [Pa a^1/nn]
	:type nn: float
	:param nn: flowlaw exponent associated with dimensionality of B, 
			   from Glen (1955)
	:type lbda: float
	:param lbda: bed wavelength in [m]
	:type US: float
	:param US: sliding velocity in [m/a]
	:type npts: int
	:param npts: number of points to discretize x-domain with, default is 1000
	:type output: str
	:param output: output format
					'tau': (default) just output estimate of shear stress (TAU)
					'params': TAU, SS
					'model': TAU, gx, hx
					'full': TAU, SS, gx, hx
	:: RETURN ::
	:rtype TAU: float
	:
	:rtype SS: float
	:return SS: ice-bed contact fraction [length / wavelength]
	"""
	# Convert sliding velocity into m/sec
	Uss = US/(365*24*3600)
	# Discretize model domain
	xv = np.arange(0,lbda,1/int(npts))
	# Calculate length-scale of the cavity in the absence of subsequent bumps (lc) (Equation A2 in Zoet & Iverson, 2015)
	lc = np.sqrt(8.*Uss*h/np.pi*(B/N)**nn)
	# Calculate cavity roof elevation vector (g(x) - equation 4 in Kamb, 1987)
	gx = np.real(h*(0.5 - (1./np.pi)*np.arcsin((2.*xv-lc)/lc) - (2.*(xv - lc) * np.sqrt(xv*(lc - xv))/(np.pi*lc**2))))
	# Calculate wavenumber
	kk = 2.*np.pi/lbda
	# Calculate bed amplitude
	aa = h/2.
	# Calculate bed elevation vector (h(x) - defines a sinusoidal bed with wavelength lambda)
	# hx = (h/2.)*np.sin((2.*np.pi*xv/lbda) + np.pi/2.) + h/2.
	hx = aa*(np.cos(kk*xv) + 1.)
	# Calculate cavity space
	dd = gx - hx
	# breakpoint()
	# Get indices of cavity edges
	try:
		sa = np.max(np.argwhere(dd > 0))
		# Get first index of positive cavity volume
		sb = np.min(np.argwhere(dd > 0))
	# If model says no cavity, assign cavity start and end as same point
	except:
		sa = 0
		sb = 0
		# SS = 1.
	# Calculate bed contact fraction
	SS = 1. - ((xv[sa] - xv[sb])/lbda)
	# Calculate wave-number scaled cavity critical length (k*x_c)
	kxc = acot((2.*np.pi*(1. - SS) + np.sin(2.*np.pi*SS))/(1. - np.cos(2.*np.pi*SS)))
	# Calculate fo factor (Coefficient of friction)
	# fo = ((h/2.)/lbda)*np.pi*N
	fo = aa*kk*0.5*N
	# Calculate geometric factor
	PHI = ((np.pi*SS - 0.5*np.sin(2.*np.pi*SS))*np.sin(np.pi*SS - kxc))/(np.sin(np.pi*SS) - np.pi*SS*np.cos(np.pi*SS))
	# Calculate drag
	TAU = fo*PHI


	if output == 'tau':
		return TAU
	elif output == 'model':
		return TAU, gx, hx
	elif output == 'full':
		return TAU, SS, gx, hx
	elif output == 'params':
		return TAU, SS
	else:
		return TAU


def calc_TS_vectors(N_min,N_max,nnods=5001,**kwargs):
	TV = []; SV = []
	for i_,N_ in tqdm(enumerate(np.linspace(N_min,N_max,nnods))):
		iT,iS = model_TS_single(N_,output='params',npts=nnods,**kwargs)
		TV.append(iT)
		SV.append(iS)
	return {'T_Pa':np.array(TV),'S':np.array(SV),'N_Pa':np.linspace(N_min,N_max,nnods)}

# CODE VALIDATION SECTION

def run_ISU_test2(Uval):
	CLR = .4
	CLC = 2.*np.pi*CLR
	CLW = .183
	CLa = .0153
	BB = 6.3e7
	nn = 3.
	NN = 500e3
	SS = calc_SS(NN=NN,BB=BB,nn=nn,lbda=CLW,US=Uval,npts=10001)
	TT = calc_from_SN(NN,SS,hh=2*CLa,lbda=CLW)
	return SS, TT


def run_ISU_test(Uval):
	CLR = .4
	CLC = 2.*np.pi*CLR
	CLW = .183
	CLa = .0153
	DC_Tau = 9e3 # [Pa] Wall drag
	T_ = model_TS_single(500e3,h=2*CLa,B=6.3e7,nn=3,lbda=CLW,US=Uval,npts=10001)
	return T_ + DC_Tau


def run_ISU_tests():
	Uvals = np.linspace(2,200,100)
	Tvals = []
	for U_ in Uvals:
		Tvals.append(run_ISU_test(U_))
	return (Uvals,Tvals)



### These are computationally expensive. Just use linear interpolation
# def loop_N2Tau(N_vect,output='tau',**kwargs):
# 	TV = np.zeros(len(N_vect),)
# 	for i_,N_ in enumerate(N_vect):
# 		TV[i_] += model_TS_single(N_,output=output,**kwargs)
# 	return TV

# def loop_N2TauS(N_vect,output='params',**kwargs):
# 	# TV = np.zeros(len(N_vect),)
# 	# SV = np.zeros(len(N_vect),)
# 	for i_, N_ in tqdm(enumerate(N_vect)):
# 		iT,iS = model_TS_single(N_,output=output,**kwargs)
# 		TV[i_] += iT
# 		SV[i_] += iS
# 	return TV, SV




# S = np.linspace(0,1,1000)
