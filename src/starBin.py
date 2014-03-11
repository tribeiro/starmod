
'''
	Class for stars in binary. The surface of the star is defined by Roche equipontential surfaces.
'''

from baseAstroObj import *
import healpy as H

################################################################################

class StarBinary(BaseAstroObj):

	############################################################################

	def __init__(self,i,q,ffac=1.0,m=1.0,nside=16,**kwargs):

		BaseAstroObj.__init__(self,	i=i,
									q=q,
									ffac=ffac,
									m=m,
									nside=nside,
									**kwargs)

		self.__MAX_RL_ITER__ = 100
		# Make surface
		self.npix = H.nside2npix(self.nside)
		self.sresize(self.npix)
		self.RL1 = self.rl1()
		self.r = self.RL1
		
		for ipix in range(self.npix):
			self.phi[ipix],self.theta[ipix] = H.pix2ang(self.nside,ipix)
			self.phi[ipix] -= np.pi/2.0

		#self.nx = np.sin(self.phi)
		#self.ny = np.cos(self.phi)*np.cos(self.theta)
		#self.nz = np.cos(self.phi)*np.sin(self.theta)

		#rlob = np.zeros(self.npix)
		#niter = np.zeros(self.npix,dtype=int)
		
		#for ipix in range(self.npix):
		#	rlob[ipix],niter[ipix] = self.rocheLobe(ipix)
		
		rlob,niter = self.rocheLobeAll()
		
		print '# - NITER = %i'%(niter)
		#print len(rlob)
		#print '# - Std. NITER = %i'%(np.std(niter))
		#print '# - Min. NITER = %i'%(np.min(niter))
		#print '# - Max. NITER = %i'%(np.max(niter))
		
		self.px = rlob*np.sin(self.phi)
		self.py = rlob*np.cos(self.phi)*np.cos(self.theta)
		self.pz = rlob*np.cos(self.phi)*np.sin(self.theta)

		nr = np.sqrt(self.px*self.px+self.py*self.py+self.pz*self.pz)
		
		n_r = np.sqrt((1.-self.px)*(1.-self.px)+self.py*self.py+self.pz*self.pz)
		nr = nr*nr*nr
		n_r = n_r*n_r*n_r
		
		self.nx = -self.px/nr+(self.q*(1.-self.px)/n_r)+(self.q+1)*(self.px)-self.q;
		self.ny = -self.py*(1./nr+self.q/n_r-(self.q+1.));
		self.nz = -self.pz*(1./nr+self.q/n_r);

		#np = np.sqrt(nx*nx+ny*ny+nz*nz);

		self.vx =  np.cos(self.phi)
		self.vy = -np.sin(self.phi)*np.cos(self.theta)
		self.vz = -np.sin(self.phi)*np.sin(self.theta)

		self.I+=1

	############################################################################

	def rocheLobe(self,ipix):
		'''
		Find radius of the Roche lobe for specified angles.
		'''

		v0 = self.theta[ipix]
		u0 = self.phi[ipix]
		
		R1 = 0.
		R2 = self.RL1

		a = 2./ ( ( 1.+self.q )*R2 )
		b = 2.*self.q/( ( 1.+self.q ) * np.sqrt( (R2*R2) - (2.*R2)+1.0 ) )
		c = ( R2 - ( self.q /( 1.0+self.q ) ) )**2.0
		
		fimax = (a+b+c)*self.ffac

		for i in range(self.__MAX_RL_ITER__):

			R = (R2+R1)/2.
			r_2 = np.sqrt(R**2.-2.*R*np.sin(u0)+1.)
			x = R*np.sin(u0)
			y = R*np.cos(u0)*np.cos(v0)
			
			a = 2./( (1.+self.q) * R )
			b = 2.*self.q / ( ( 1.+self.q )*r_2 )
			c = (x - self.q / ( 1.+self.q ) )**2.
			d = y**2.
			
			fi = a+b+c+d;

			#print '%i %f %f %f %f %f %f'%(i,R,R1,R2,fi,fimax,fi-fimax)

			if(fi < fimax):
				R2 = R
			elif(fi > fimax):
				R1 = R
				
			if(np.abs(fi-fimax) < 1e-5):
				return R,i

		raise RuntimeError('''Maximum number of interations reached (NITER=%i).
Could not find Roche lobe radius (pix = %i). Pot = %+8.2e'''%(self.__MAX_RL_ITER__,ipix,fi))
		
	############################################################################

	def rocheLobeAll(self):
		'''
		Find radius of the Roche lobe for specified angles.
		'''

		v0 = self.theta
		u0 = self.phi
		
		R1 = np.zeros(self.npix)
		R2 = self.RL1+np.zeros(self.npix)
		rNotFound = np.zeros(self.npix)==0 # Mask for pixels that were not found yet

		a = 2./ ( ( 1.+self.q )*R2 )
		b = 2.*self.q/( ( 1.+self.q ) * np.sqrt( (R2*R2) - (2.*R2)+1.0 ) )
		c = ( R2 - ( self.q /( 1.0+self.q ) ) )**2.0
		
		fimax = (a+b+c)*self.ffac

		for i in range(self.__MAX_RL_ITER__):

			R = (R2+R1)/2.
			r_2 = np.sqrt(R**2.-2.*R*np.sin(u0)+1.)
			x = R*np.sin(u0)
			y = R*np.cos(u0)*np.cos(v0)
			
			a = 2./( (1.+self.q) * R )
			b = 2.*self.q / ( ( 1.+self.q )*r_2 )
			c = (x - self.q / ( 1.+self.q ) )**2.
			d = y**2.
			
			fi = a+b+c+d;

			#print '%i %f %f %f %f %f %f'%(i,R,R1,R2,fi,fimax,fi-fimax)
			m2 = fi <= fimax
			m1 = fi > fimax
			R2[m2] = R[m2]
			R1[m1] = R[m1]
			
			rNotFound = np.abs(fi-fimax) > 1e-5
			
			if not rNotFound.any() :
				return R,i

		raise RuntimeError('''Maximum number of interations reached (NITER=%i).
Could not find Roche lobe radius.'''%(self.__MAX_RL_ITER__))
		

	############################################################################

	def rl1(self):
		'''
		Find distance to internal lagrange point (RL1).
		'''

		x1 = 0.
		x2 = 1.

		for i in range(self.__MAX_RL_ITER__):

			x = ((x2-x1)/2.0)+x1
			a = -2.*x/((1.+self.q)*(x*x)**(3./2.))
			b = -2.*self.q*(x-1.)/((1.+self.q)*((x-1)**2.)**(3./2.))
			c = 2.*(x-(self.q/(1.0+self.q)))
			dfi = a+b+c
			
			if dfi >= 0.:
				x2 = x
			elif dfi < 0.:
				x1 = x
			
			if (np.abs(dfi) < 1.e-5 ):
				return x
				
		raise RuntimeError('Could not find RL1...')

	############################################################################

################################################################################








