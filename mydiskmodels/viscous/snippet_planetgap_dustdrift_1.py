from snippet_header import DiskRadialModel, np, plt, MS, au, finalize, year, Mea
import copy
tstart    = 1e4 * year
tend      = 2e6 * year
ntime     = 100
nr        = 2000
alpha     = 1e-4
agrain    = 1e-1
mplanet   = 100 * Mea
aplanet   = 5 * au
mdisk0    = 1e-2 * MS
rdisk0    = 10 * au

d = DiskRadialModel(rin=0.1 * au, rout=1000 * au, nr=nr)
d.make_disk_from_lbp_alpha(mdisk0, rdisk0, alpha, tstart)
d.add_dust(agrain=agrain)
d.alphamix = d.alpha      # Make separate alpha for mixing
d.add_planet_gap(aplanet, 'duffell', mpl=mplanet, smooth=2., log=True, innu=True)

# iteration

time = tstart * (tend / tstart)**(np.linspace(0., 1., ntime + 1))
dlist = [copy.deepcopy(d)]
for itime in range(1, ntime + 1):
    dt = time[itime] - time[itime - 1]
    d.compute_viscous_evolution_and_dust_drift_next_timestep(dt)
    dlist.append(copy.deepcopy(d))

# plotting

plt.figure()
for itime in range(0, ntime + 1, 20):
    s = '{0:8.2e} years'.format(time[itime] / year)
    plt.plot(dlist[itime].r / au, dlist[itime].dust[0].sigma, label=s)
plt.xscale('log')
plt.yscale('log')
plt.ylim(bottom=1e-6)
plt.xlabel('r [au]')
plt.ylabel(r'$\Sigma_d$')
plt.legend()

finalize(results=(dlist[-1].dust[0].sigma))
