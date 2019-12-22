from disklab.diskradial import *
from disklab.natconst import *
import matplotlib.pyplot as plt

Rc = 30.0*au
Sigc = 10**2.5
gam = 1.0
alpha = 0.001
d    = DiskRadialModel(alpha=alpha)

d.make_disk_from_simplified_lbp(Sigc,Rc,gam)

# agrain = np.array([1.6e-4, 8.7e-2, 0.13, 0.9])	# observed wavelength

n 		    = 5
amin		= 1.6e-4
amax		= 0.9
agraini 	= amin * (amax/amin)**np.linspace(0.,1.,n+1)
agrain      = 0.5 * ( agraini[1:] + agraini[:-1] )
xi		    = 2.
mgrain 		= (4.*np.pi/3.)*xi*agrain**3
gamma 		= -3.5
dtgr		= 0.01  				# dust to gas ratio
abun 		= agrain**(gamma+4.)
abun		/= dtgr*abun.sum()

d.make_disk_from_simplified_lbp(Sigc,Rc,gam)
for ia in range(n):
    d.add_dust(agrain=agrain[ia], xigrain=xi, dtg=abun[ia])
sigmadust = d.dust[0].join_multi_array(d.dust)  # Make a 2-D array: Sigma(r,a)
dlna = np.log(agraini[1]) - np.log(agraini[0])  # Only valid for log grid in a

ntime = 10
tstart = 1e3 * year     # Initial time at 1000 year (for log-spaced timesteps)
tend = 5e6 * year       # Final time at 5 Myr
time = tstart * (tend / tstart)**(np.linspace(0., 1., ntime+1))

plt.figure()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('r [au]')
plt.ylabel(r'$\Sigma$')
plt.plot(d.r / au, d.sigma)
print('Disk mass = {} Msun'.format(d.mass/MS))


for na in range(n):
    for itime in range(1, ntime + 1):
        d.sigma = d.get_viscous_evolution_next_timestep(time[itime] - time[itime - 1])
        d.dust[na].sigma = d.dust[na].get_dust_radial_drift_next_timestep(time[itime] - time[itime - 1],fixgas=False)
        d.dust[na].compute_mass()
    plt.plot(d.r / au, d.sigma)

plt.show()

d.compute_mass()
print('Evolved disk mass = {} Msun'.format(d.mass/MS))
