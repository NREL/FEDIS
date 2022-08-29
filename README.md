   '''
    
    The “Fresnel Equations” for Diffuse radiation on Inclined photovoltaic Surfaces (FEDIS)
    is an analytical solution of diffuse transmission based on the rigorous integration of
    an alternate form of the Fresnel equations. The approach leads to a simple yet accurate
    relative transmittance model (often referred to as incident angle modifier or angle of 
    incidence model) that reconciles the solar energy sensed by pyranometers and PV panels.

    Parameters
    ----------
    aoi : Angle of incidence in degrees.

    rfnt: the refractive index of the pyranometer cover
          For a fused silica dome over a CMP22, the rfnt is 1.4585

    rfn: the refractive index of the PV cover
         The suggested value is 1.5

    surface_tilt : numeric
        Surface tilt angles in decimal degrees.
        The tilt angle is defined as degrees from horizontal
        (e.g. surface facing up = 0, surface facing horizon = 90).

    Returns
    -------
    cd : the incidence angle modifier (IAM) for direct radiation
    cuk: the IAM for diffuse radiation from the sky
    cug: the IAM for diffuse radiation from the ground reflection

    Usage
    ----------
    The solar energy received by a PV can be given by
    F = cd*Fd + cuk*Fuk + cug*Fug
    Fd, Fuk, and Fug are the POA irradiances, observed by the pyranometer, 
    that are associated with direct radiation, diffuse radiation from the sky, and 
    diffuse radiation from ground reflection, respectively.
    An example can be found in Eq.(1-2) from the reference

    rfnt = 1.4585
    rfn = np.arange(nPV)
    aoi = np.arange(n)
    surface_tilt = np.arange(n)
    cd, cuk, cug = FEDIS(aoi, surface_tilt, rfn, rfnt )
        Will return 2 dimensional arrays (nPV, n)
        n represents the number of total scenarios
        nPV represents the number of PV cover types

    Example
    -------
    belta = np.arange(100)*0.9
    rfn = np.array([1.3, 2.0]) 
    theta0p = belta*np.pi/180.0
    aoi = theta0p*180.0/np.pi
    surface_tilt = belta

    cd, cuk, cug = FEDIS(aoi, surface_tilt, rfn, rfnt )
    print( cd[0][:] )
    print( cuk[0][:] )
    print( cug[0][:] )
    This generate the results of Fig.3 in the reference

    Reference
    ----------
    Xie, Y., M. Sengupta, A. Habte, A. Andreas, The "Fresnel Equations for Diffuse 
    radiation on Inclined photovoltaic Surfaces (FEDIS), Rew. Sus. En. Rev., 161, 112362

    '''

