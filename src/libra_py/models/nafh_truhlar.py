def model_nafh_adiab_nD(ndof,q,nsurf,model_params):
    """Returns the value of the (adiabatic) NaFH potential computed from the PES of Truhlar et. al. (CITE)

    Args:
        ndof (integer): The number of degrees of freedom.

        q (MATRIX): List of basis positions extracted from a qpas list.

        model_params (dictionary): Dictionary containing the potential parameters.

        nsurf (integer): Parameter specifying the surface under consideration.

    Returns:
        energy (complex): Potential energy V at the specified point on the specified surface.

        denergy (complex): 1st derivative of the potential energy V at the specified point on the specified surface.

        d2energy (complex): 2nd derivative of the potential energy V at the specified point on the specified surface.
    """

    r=np.zeros([1,3],dtype='d')
    e=np.zeros([1],dtype='d')
    de=np.zeros([3,1],dtype='d')

    if (ndof == 1):
        rFH=q.get(0);rNaF=3.779
    elif (ndof == 2):
        rFH=q.get(0);rNaF=q.get(1)

    r[0,0]=rNaF+rFH
    r[0,1]=rFH
    r[0,2]=rNaF

    if nsurf == 2:
        nsurf_conv=1
    elif nsurf == 3:
        nsurf_conv=2
    else:
        nsurf_conv=3

    nafh2dpy.prepot()
    nafh2dpy.pot(r,e,de,nsurf_conv,1)

    energy=e[0]
    if (ndof == 1):
        denergy=[de[1,0]]
    elif (ndof == 2):
        denergy=[de[1,0],de[2,0]]

    d2energy=0.0

    return(energy,denergy,d2energy)
