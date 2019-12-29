from disklab import natconst as nc

# Modify the built-in drift function to go back in time
def get_dust_radial_drift_previous_timestep(self, dt, alphamodel=True, fixgas=False, extracond=None,
                                        bcin='zerogradient', bcout='zerogradient'):
    """
    Advance the dust component one time step into the future.
    Radial drift and turbulent mixing included, as well as the
    gas drag as the gas is moving inward.

    ARGUMENTS:
    dt          = Time step in seconds
    alphamodel  = If True, then recompute self.nu from alpha-recipe (default)
    fixgas      = If True, then do *not* include the inward gas motion in dust drift.
    extracond   = (for special purposes only) List of extra internal conditions
    bcin,bcout  = String denoting the kind of boundary condition for the dust
                  'zerogradient'   : Default: simply set gradient of Sigma_dust to zero
                  'closed'         : Do not allow dust to cross this border

    Note: If self.diskradialmodel.alphamix is present, then this alpha will be used (instead of the
    usual self.alpha) for the turbulent mixing.

    ** BEWARE: **
    Always make sure to have updated the midplane density and temperature,
    and then call the compute_stokes_from_agrain() before calling this subroutine,
    if you have evolved the gas beforehand.
    """
    #
    # If requested, compute nu and dmix from alpha
    #
    if alphamodel:
        self.diskradialmodel.compute_nu()
        if hasattr(self.diskradialmodel, 'alphamix'):
            self.dmix = self.diskradialmodel.alphamix * self.diskradialmodel.cs * self.diskradialmodel.cs / self.omk / self.diskradialmodel.Sc
        else:
            self.dmix = self.diskradialmodel.nu / self.diskradialmodel.Sc
        self.dmix[:] *= 1.0 / (1.0 + self.St ** 2)
    #
    # Cast into diffusion equation form
    #
    x = self.r
    y = 2 * nc.pi * self.r * self.sigma  # Dust
    g = 2 * nc.pi * self.r * self.diskradialmodel.sigma  # Gas
    d = self.dmix
    di = 0.5 * (d[1:] + d[:-1])
    self.compute_dustvr_at_interfaces(alphamodel=alphamodel, fixgas=fixgas)
    vi = self.vr
    #
    # Source term only if self.sigdot is given
    #
    if hasattr(self, 'sigdot'):
        s = 2 * nc.pi * self.r * self.sigdot
    else:
        s = np.zeros(len(x))  # For now no source term
    #
    # Set boundary conditions
    #
    if bcin == 'zerogradient':
        bcl = (1, 0, 0, 1)  # Simply set dy/dx=0 at inner edge
    elif bcin == 'closed':
        bcl = (1, 0, 0, -1)  # Closed boundary at inner edge
    else:
        raise (ValueError('Inner boundary condition for dust invalid'))
    if bcout == 'zerogradient':
        bcr = (1, 0, 0, 1)  # Simply set dy/dx=0 at outer edge
    elif bcout == 'closed':
        bcr = (1, 0, 0, -1)  # Closed boundary at outer edge
    else:
        raise (ValueError('Outer boundary condition for dust invalid'))
    #
    # Get the new value of y after one time step dt
    #
    y = solvediffonedee(x, y, vi, di, g, s, bcl, bcr, dt=dt,
                        int=True, upwind=True, extracond=extracond)
    #
    # Obtain new sigdust
    #
    sigma = y / (2 * nc.pi * self.r)
    #
    # Return
    #
    return sigma
