"""
python module for calculating integrals with 1 Bessel functions.
This module contains:
-- FFTLog method for integrals with 1 spherical Bessel function;
-- integrals with 1 (cylindrical) Bessel function, i.e. Hankel transform;
-- window function (optional) for smoothing Fourier coefficients

by Xiao Fang
Apr 10, 2019
"""

import numpy as np
from scipy.special import gamma
from numpy.fft import rfft, irfft

class fftlog(object):

	def __init__(self, x, fx, ell, nu, N_extrap_low=0, N_extrap_high=0, c_window_width=0.25):

		self.x_origin = x # x is logarithmically spaced
		# self.lnx = np.log(x)
		self.dlnx = np.log(x[1]/x[0])
		self.fx_origin= fx # f(x) array
		self.ell = ell
		self.nu = nu
		self.N_extrap_low = N_extrap_low
		self.N_extrap_high = N_extrap_high
		self.c_window_width = c_window_width

		# extrapolate x and f(x) linearly in log(x), and log(f(x))
		self.x = log_extrap(x, N_extrap_low, N_extrap_high)
		self.fx = log_extrap(fx, N_extrap_low, N_extrap_high)
		self.N = self.x.size

		if(self.N%2==1): # Make sure the array sizes are even
			self.x= self.x[:-1]
			self.fx=self.fx[:-1]
			self.N -= 1
			self.N_extrap_high -=1

		self.m, self.c_m = self.get_c_m()
		self.eta_m = 2*np.pi/(float(self.N)*self.dlnx) * self.m


	def get_c_m(self):
		"""
		return m and c_m
		c_m: the smoothed FFT coefficients of "biased" input function f(x): f_b = f(x) / x^\nu
		number of x values should be even
		c_window_width: the fraction of c_m elements that are smoothed,
		e.g. c_window_width=0.25 means smoothing the last 1/4 of c_m elements using "c_window".
		"""

		f_b=self.fx * self.x**(-self.nu)
		c_m=rfft(f_b)
		m=np.arange(0,self.N//2+1) 
		c_m = c_m*c_window(m, int(self.c_window_width*self.N//2.) )
		return m, c_m

	def fftlog(self):
		"""
		Calculate F(y) = \int_0^\infty dx / x * f(x) * j_\ell(xy),
		where j_\ell is the spherical Bessel func of order ell.
		array y is set as y[:] = (ell+0.5)/x[::-1]
		"""
		x0 = self.x[0]
		z_ar = self.nu + 1j*self.eta_m
		y = (self.ell+1.) / self.x[::-1]
		h_m = self.c_m * (self.x[0]*y[0])**(-1j*self.eta_m) * g_l(self.ell, z_ar)

		Fy = irfft(np.conj(h_m)) * y**(-self.nu) * np.sqrt(np.pi)/4.
		return y[self.N_extrap_high:self.N-self.N_extrap_low], Fy[self.N_extrap_high:self.N-self.N_extrap_low]

class hankel(object):
	def __init__(self, x, fx, n, nu, N_extrap_low=0, N_extrap_high=0, c_window_width=0.25):
		print('nu is required to be between (0.5-n) and 2.')
		self.myfftlog = fftlog(x, np.sqrt(x)*fx, n-0.5, nu, N_extrap_low, N_extrap_high, c_window_width)
	
	def hankel(self):
		y, Fy = self.myfftlog.fftlog()
		Fy *= np.sqrt(2*y/np.pi)
		return y, Fy


### Utility functions ####################

def log_extrap(x, N_extrap_low, N_extrap_high):

	low_x = high_x = []
	if(N_extrap_low):
		dlnx_low = np.log(x[1]/x[0])
		low_x = x[0] * np.exp(dlnx_low * np.arange(-N_extrap_low, 0) )
	if(N_extrap_high):
		dlnx_high= np.log(x[-1]/x[-2])
		high_x = x[-1] * np.exp(dlnx_high * np.arange(1, N_extrap_high+1) )
	x_extrap = np.hstack((low_x, x, high_x))
	return x_extrap

def c_window(n,n_cut):
	"""
	One-side window function of c_m,
	Adapted from Eq.(C1) in McEwen et al. (2016), arXiv:1603.04826
	"""
	n_right = n[-1] - n_cut
	# n_left = n[0]+ n_cut 
	n_r=n[ n[:]  > n_right ] 
	# n_l=n[ n[:]  <  n_left ] 
	theta_right=(n[-1]-n_r)/float(n[-1]-n_right-1) 
	# theta_left=(n_l - n[0])/float(n_left-n[0]-1) 
	W=np.ones(n.size)
	W[n[:] > n_right]= theta_right - 1/(2*np.pi)*np.sin(2*np.pi*theta_right)
	# W[n[:] < n_left]= theta_left - 1/(2*pi)*sin(2*pi*theta_left)
	return W

def g_m_vals(mu,q):
	'''
	g_m_vals function is adapted from FAST-PT
	'''
	imag_q= np.imag(q)
	g_m=np.zeros(q.size, dtype=complex)
	cut =200
	asym_q=q[np.absolute(imag_q) >cut]
	asym_plus=(mu+1+asym_q)/2.
	asym_minus=(mu+1-asym_q)/2.

	q_good=q[ (np.absolute(imag_q) <=cut) & (q!=mu + 1 + 0.0j)]

	alpha_plus=(mu+1+q_good)/2.
	alpha_minus=(mu+1-q_good)/2.
	
	g_m[(np.absolute(imag_q) <=cut) & (q!= mu + 1 + 0.0j)] =gamma(alpha_plus)/gamma(alpha_minus)

	# asymptotic form 								
	g_m[np.absolute(imag_q)>cut] = np.exp( (asym_plus-0.5)*np.log(asym_plus) - (asym_minus-0.5)*np.log(asym_minus) - asym_q \
		+1./12 *(1./asym_plus - 1./asym_minus) +1./360.*(1./asym_minus**3 - 1./asym_plus**3) +1./1260*(1./asym_plus**5 - 1./asym_minus**5) )

	g_m[np.where(q==mu+1+0.0j)[0]] = 0.+0.0j
	return g_m

def g_l(l,z_array):
	'''
	gl = 2.**z_array * gamma((l+z_array)/2.) / gamma((3.+l-z_array)/2.)
	'''
	gl = 2.**z_array * g_m_vals(l+0.5,z_array-1.5)
	return gl

if __name__ == '__main__':
	import matplotlib.pyplot as plt

	print('This is a test of fftlog module written by Xiao Fang.')
	print('nu is required to be between -ell to 2.')
	k, pk = np.loadtxt('Pk_test', usecols=(0,1), unpack=True)
	N = k.size
	print('number of input data points: '+str(N))
	ell = 1
	nu = 1.
	myfftlog = fftlog(k, pk, ell=ell, nu=nu, N_extrap_low=1500, N_extrap_high=1500, c_window_width=0.25)
	r, Fr = myfftlog.fftlog()

	fig = plt.figure(figsize=(8,4))
	fig.suptitle(r'$F(y) = \int_0^\infty f(x)j_{\ell}(xy) dx/x, \ell=$%.1f'%(ell))

	subfig1 = fig.add_subplot(1,2,1)
	subfig1.set_xscale('log')
	subfig1.set_yscale('log')
	subfig1.set_xlabel('x')
	subfig1.set_ylabel('f(x)')
	subfig1.plot(k, pk)
	plt.tight_layout()

	subfig2 = fig.add_subplot(1,2,2)
	subfig2.set_title(r'$\nu=$%.2f'%(nu))
	subfig2.set_xscale('log')
	subfig2.set_yscale('log')
	subfig2.set_xlabel('y')
	subfig2.set_ylabel('F(y)')
	subfig2.plot(r, Fr)

	r_bf, Fr_bf = np.loadtxt('test_bruteforce.txt', usecols=(0,1), unpack=True)
	subfig2.plot(r_bf, Fr_bf)

	plt.tight_layout()
	plt.show()

	#################

	n = 0
	nu = 1.
	myhankel = hankel(k, pk, n=n, nu=nu, N_extrap_low=1500, N_extrap_high=1500, c_window_width=0.25)
	r, Fr = myhankel.hankel()

	fig = plt.figure(figsize=(8,4))
	fig.suptitle(r'$F(y) = \int_0^\infty f(x)J_{n}(xy) dx/x, n=$%.1f'%(n))

	subfig1 = fig.add_subplot(1,2,1)
	subfig1.set_xscale('log')
	subfig1.set_yscale('log')
	subfig1.set_xlabel('x')
	subfig1.set_ylabel('f(x)')
	subfig1.plot(k, pk)
	plt.tight_layout()

	subfig2 = fig.add_subplot(1,2,2)
	subfig2.set_title(r'$\nu=$%.2f'%(nu))
	subfig2.set_xscale('log')
	subfig2.set_yscale('log')
	subfig2.set_xlabel('y')
	subfig2.set_ylabel('F(y)')
	subfig2.plot(r, Fr)
	plt.tight_layout()
	plt.show()
