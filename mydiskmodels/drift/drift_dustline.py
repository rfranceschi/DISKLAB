from snippet_header import DiskRadialModel, np, plt, MS, year, au, finalize
from disklab.natconst import pi

tstart = 1e3 * year
tend = 5e6 * year
ntime = 100
nr = 1000
time = tstart * (tend / tstart)**(np.linspace(0., 1., ntime + 1))  # log spacing for disk age

d = DiskRadialModel(rout=1000 * au, nr=nr)

# Simplified Lynden-Bell & Pringle density distribution
Rc = 30.0 * au
Sigc = 175
gam = 1.0
d.make_disk_from_simplified_lbp(Sigc, Rc, gam)

alpha = 0.1  # viscous alpha
wl = 0.087  # observed wavelength
agrain = wl/(2*pi)  # corresponding grain size
dtg = 0.001  # dust to gas ratio
d.add_dust(agrain=agrain, dtg=dtg)
# d.Sc = 1e10    # switch off mixing by putting Schmidt number to 'infinity'

# Plot initial profile

fig = plt.figure()

ax = plt.axes()

ax_start, = ax.loglog(d.r / au, d.dust[0].sigma, label=r'$10^3$ yr')

# Initial dust line at an arbitrary total mass threshold
th = 0.95
d.dust[0].compute_mass()
m_int = d.dust[0].mass
print(f'Initial dust mass = {m_int/MS:.2e} M_S')

integrand = d.dust[0].sigma*2*pi*d.r
i = 0
while m_int > th*d.dust[0].mass:
    i = i+1
    m_int = np.trapz(integrand[:-i], d.r[:-i])
print(f'Initial dustline at {d.r[-i]/au:.0f} au')
plt.axvline(x=d.r[-i] / au, color='b')

# Drift evolution
for itime in range(1, ntime + 1):
    d.dust[0].sigma = d.dust[0].get_dust_radial_drift_next_timestep(time[itime] - time[itime - 1], fixgas=True)


d.dust[0].compute_mass()
m_int = d.dust[0].mass
print(f'Final dust mass = {m_int/MS:.2e} M_S')

integrand = d.dust[0].sigma*2*pi*d.r
i = 0
while m_int > th*d.dust[0].mass:
    i = i+1
    m_int = np.trapz(integrand[:-i], d.r[:-i])
print(f'Final dustline at {d.r[-i]/au:.0f} au')

ax_end, = ax.loglog(d.r / au, d.dust[0].sigma, label=r'$5$ Myr')
plt.axvline(x=d.r[-i] / au, color='orange')

ax.legend()
plt.title(f'Dustline drift evolution, a_grain = {agrain:.2e} cm')
plt.show()
