


def FEDIS(aoi, surface_tilt, rfn, rfnt=1.4585):

    '''
    The “Fresnel Equations” for Diffuse radiation on Inclined photovoltaic Surfaces (FEDIS)
    is an analytical solution of diffuse transmission based on the rigorous integration of
    an alternate form of the Fresnel equations. The approach leads to a simple yet accurate
    relative transmittance model that reconciles the solar energy sensed by pyranometers and PV panels.

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

    surface_tilt[surface_tilt<0.01] = 0.01
    aoi[aoi==0.0] = 0.001

    r0 = ( (rfnt-1.0)/(rfnt+1.0) )**2.0
    aoi = aoi*np.pi/180.0
    theta0tp = [] 
    for rfn1 in rfn:
        theta0tp1 =  np.arcsin( np.sin(aoi)/rfn1 ) 
        theta0tp.append( theta0tp1 )
    theta0tp = np.asarray(theta0tp, dtype=np.float)

    rd = [ ]
    for i in range( rfn.shape[0] ):
        rd1 = ( np.sin(aoi-theta0tp[i,:])/np.sin(aoi+theta0tp[i,:]) )**2.0 + \
             ( np.tan(aoi-theta0tp[i,:])/np.tan(aoi+theta0tp[i,:]) )**2.0
        rd1 = rd1*0.5
        rd.append( rd1 )
    rd = np.asarray(rd, dtype=np.float)

    cd = [ ]
    for j in range( surface_tilt.shape[0] ):
        cd1 = (1.0-rd[:,j])/(1.0-r0)
        cd.append( cd1 )
    cd = np.asarray(cd, dtype=np.float)
    cd = cd.T

    cuk = [ ]
    cug = [ ]
    w = 2.77526e-09  +3.74953*rfn -5.18727*rfn**2.0 +3.41186*rfn**3.0 -1.08794*rfn**4.0 +0.136060*rfn**5.0
    w = w*(rfn*(1.0+rfnt)**2.0)/( rfnt*(1.0+rfn)**2.0 )

    for i in range( rfn.shape[0] ):
        w1 = w[i]
        cuk1 = 30.0*np.pi/7.0 - (160.0/21.0)*(surface_tilt*np.pi/180.0) - (10.0*np.pi/3.0)*np.cos(surface_tilt*np.pi/180.0) + \
            (160.0/21.0)*np.cos(surface_tilt*np.pi/180.0)*np.sin(surface_tilt*np.pi/180.0) \
            - (5.0*np.pi/3.0)*np.cos(surface_tilt*np.pi/180.0)*( np.sin(surface_tilt*np.pi/180.0) )**2.0 \
            + (20.0/7.0)*np.cos(surface_tilt*np.pi/180.0)*( np.sin(surface_tilt*np.pi/180.0) )**3.0 \
            - (5.0*np.pi/16.0)*np.cos(surface_tilt*np.pi/180.0)*( np.sin(surface_tilt*np.pi/180.0) )**4.0 \
            + (16.0/105.0)*np.cos(surface_tilt*np.pi/180.0)*( np.sin(surface_tilt*np.pi/180.0) )**5.0
        cuk1 = cuk1*w1*2.0/(np.pi*(1.0+np.cos(surface_tilt*np.pi/180.0)))
        cuk.append( cuk1 )
        cug1 = 40.0*w1/(21.0*(1.0-np.cos(surface_tilt*np.pi/180.0))) - \
               (1.0+np.cos(surface_tilt*np.pi/180.0))*cuk1/(1.0-np.cos(surface_tilt*np.pi/180.0))
        cug.append( cug1 )
    cuk = np.asarray( cuk, dtype=np.float)
    cug = np.asarray( cug, dtype=np.float)

    return cd, cuk, cug




