from snippet_header import np, plt, DiskRadialModel, au, year, MS, LS, finalize
from disklab.interactive_plot import *
#
# Model function in a form that can be used by the interactive_plot widget
#
def modelfunc(rau,param,fixedpar=None):
    #
    # Fixed parameters
    #
    mstar     = 1*MS  # Default stellar mass: Solar mass
    lstar     = 1*LS  # Default stellar luminosity: Solar luminosity
    if fixedpar is not None:
        if 'mstar' in fixedpar: mstar=fixedpar['mstar']
        if 'lstar' in fixedpar: lstar=fixedpar['lstar']
    #
    # Variable parameters (sliders of widget)
    #
    tend    = param[0]
    Sigc    = param[1]
    Rc      = param[2]
    agrain  = param[3]
    #
    # Model setup
    #
    disk      = DiskRadialModel(mstar=mstar,lstar=lstar)
    disk.make_disk_from_simplified_lbp(Sigc,Rc,1)
    disk.add_dust(agrain=agrain)
    #
    # Run the model
    #
    ntime     = 100
    time = tend**(np.linspace(0., 1., ntime+1))
    for itime in range(1,ntime+1):
        dt = time[itime]-time[itime-1]
        disk.compute_viscous_evolution_and_dust_drift_next_timestep(dt)
    #
    # Return the result (you can return whatever you want to plot; here
    # we plot two results simultaneously: the gas and dust surface
    # densities)
    #
    return np.array([disk.sigma, disk.dust[0].sigma])

#
# Create the plot we wish to make interactive
#
xmin     = 1.0
xmax     = 1e3
rau      = xmin * (xmax/xmin)**np.linspace(0.,1.,100)
par      = [5e6*year,300,30.0*au,1e-2]
fixedpar = {'mstar':MS,'lstar':LS}
result   = modelfunc(rau,par,fixedpar=fixedpar)
sigmagas = result[0]
sigmadust= result[1]
ymin     = 1e-5
ymax     = 1e+3
fig      = plt.figure()
ax       = plt.axes(xlim=(xmin,xmax),ylim=(ymin,ymax))
axgas,   = ax.loglog(rau,sigmagas,linewidth=2,label='Gas')
axdust,  = ax.loglog(rau,sigmadust,'--',linewidth=2,label='Dust')
axmodel  = [axgas,axdust]
plt.xlabel(r'$r [\mathrm{au}]$')
plt.ylabel(r'$\Sigma [\mathrm{g}/\mathrm{cm}^2]$')
plt.legend()
#
# Now make the plot interactive with sliders. We have to
# specify the slider names, the possible values and (optionally)
# the units. Then call interactive_plot() to bring it to life,
# with the function modelfunction() (see above) being the life-giver.
# Type interactive_plot? to find out more about interactive_plot():
# there are numerous examples in the document string.
#
parnames = ['tend =','Sigc =','Rc =','agrain = ']
params   = [1e6*1e1**np.linspace(0.,1.,100)*year,
            1e1*1e2**np.linspace(0.,1.,100),
            1e1*1e1**np.linspace(0.,1.,100)*au,
            1.6e-4*(0.9/1.6e-4)**np.linspace(0.,1.,100)]
parunits = [year,1,au,1.]
ipar     = interactive_plot(rau,modelfunc,params,parnames=parnames,
                            parunits=parunits,fixedpar=fixedpar,
                            parstart=par,plotbutton=False,
                            fig=fig,ax=ax,axmodel=axmodel,returnipar=True)

finalize([])
