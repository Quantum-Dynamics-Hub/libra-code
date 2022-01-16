def exact_gauss_cpl(ndof,qpasi,qpasj,model_params,*args,**kwargs):
    """Returns the value *v*=<gi|Vcpl|gj> for a Gaussian coupling potential computed in exact fashion between basis functions on different surfaces defined by the parameters of *qpasi* and *qpasj*.

    Args:
        ndof (integer): The number of degrees of freedom.

        qpasi (list): List containing the 1-by-4 basis parameter matrices for the basis function at point i.

        qpasj (list): List containing the 1-by-4 basis parameter matrices for the basis function at point j.

        model_params (dictionary): Dictionary containing the potential parameters.

    Returns:
        v (complex): The exact value of the Gaussian coupling potential between basis functions i and j.

    """

    v=complex(1.0,0.0)
    q1,p1,a1=qpasi[0], qpasi[1], qpasi[2]
    q2,p2,a2=qpasj[0], qpasj[1], qpasj[2]

    for i in range(ndof):
    	et1=-(q1.get(i)-model_params['d3'][i])**2*a1.get(i)**2/2.0-(q2.get(i)-model_params['d3'][i])**2*a2.get(i)**2/2.0+(p1.get(i)-p2.get(i))**2/2.0
    	et2=(q1.get(i)-model_params['d3'][i])*((model_params['d3'][i]-q2.get(i))*a2.get(i)+1j*(p1.get(i)-p2.get(i)))*a1.get(i)
    	et3=1j*(q2.get(i)-model_params['d3'][i])*(p1.get(i)-p2.get(i))*a2.get(i)
    	et4=(a1.get(i)+2*model_params['d2'][i]+a2.get(i))*(a1.get(i)+a2.get(i))

	v*=np.exp(2.0*model_params['d2'][i]*(et1+et2+et3)/et4)*model_params['d1'][i]* \
	np.sqrt(2.0*a1.get(i)+2.0*a2.get(i))/np.sqrt(2.0*a1.get(i)+4.0*model_params['d2'][i]+2.0*a2.get(i))

    return(v)