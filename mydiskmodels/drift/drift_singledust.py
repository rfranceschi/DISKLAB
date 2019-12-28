from snippet_header import DiskRadialModel, np, plt, MS, year, au, finalize
tstart = 1e3 * year
tend = 5e6 * year
ntime = 10
nr = 100
time = tstart * (tend / tstart)**(np.linspace(0., 1., ntime + 1))

Rc = 30.0*au
Sigc = 10**2.5
gam = 1.0
alpha = 0.001
d = DiskRadialModel(rout=1000 * au, nr=nr)
d.make_disk_from_simplified_lbp(Sigc, Rc, gam)
d.add_dust(agrain=1.6e-4)
# d.Sc = 1e10    # Switch off mixing by putting Schmidt number to 'infinity'

# plotting

plt.plot(d.r / au, d.dust[0].sigma)
for itime in range(1, ntime + 1):
    d.dust[0].sigma = d.dust[0].get_dust_radial_drift_next_timestep(time[itime] - time[itime - 1], fixgas=True)
    plt.plot(d.r / au, d.dust[0].sigma)
    d.dust[0].compute_mass()
    print('{}'.format(d.dust[0].mass / MS))
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-6, 1e3)
plt.xlabel('r [au]')
plt.ylabel(r'$\Sigma_d$')

finalize(results=d.dust[0].sigma)
