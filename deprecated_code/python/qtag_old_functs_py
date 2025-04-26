
def cls_force(univ,mss,mom_calc,props,model_params,qpas1,c1_new,qpas2,c2_new,norm2,beta):
    """Returns the values for the new basis parameter matrices on surfaces 1 (*qpas1n*) and 2 (*qpas2n*), 
       as well as their corresponding projection vectors *b1* and *b2*, where the motion of both sets of 
       functions are calculated via classical forces computed at their centers. Also necessary are the functions 
       for calculating momentum (*mom_calc*) and basis updates (*props*).

    Args:
        univ (dictionary): Dictionary containing various system parameters.

        mss (dictionary): Dictionary containing multi-surface scheme parameters.

        mom_calc (function object): The function object needed to calculate the momentum, as defined by qtag_config.

        props (list): A list of function objects for propagating basis parameters {q,p,a,s}.

        model_params (dictionary): Dictionary containing the potential parameters.

        qpas1 (list): List of {q,p,a,s} MATRIX objects for surface 1.

        c1_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 1.

        qpas2 (list): List of {q,p,a,s} MATRIX objects for surface 2.

        c2_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 2.

        norm2 (float): The population on surface 2.

        beta (float): Parameter determining tolerance in the momentum linear fitting algorithm.

    Returns:
        qpas1n (list): List of updated {q,p,a,s} MATRIX objects for surface 1.

        qpas2n (list): List of updated {q,p,a,s} MATRIX objects for surface 2.

        b1 (CMATRIX): The updated projection vector of the old basis (defined by qpas1) onto the new basis (defined by qpas1n).

        b2 (CMATRIX): The updated projection vector of the old basis (defined by qpas2) onto the new basis (defined by qpas2n).
    """

    ndof,ntraj=univ['ndof'],univ['ntraj']
    decpl=mss['decpl']
    qprop=props[0];pprop=props[1];aprop=props[2];sprop=props[3];cls_prop=props[4]
    qvals,pvals,avals,svals=MATRIX(qpas1[0]),MATRIX(qpas1[1]),MATRIX(qpas1[2]),MATRIX(qpas1[3])
    qvalsn,pvalsn,avalsn,svalsn=MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj)

    mom,r,gmom,gr=mom_calc(univ,beta,qpas1,c1_new)

    for i in range(ndof):
        qn=qprop(univ,i,qvals.row(i),mom.row(i));
        pn=pprop(mom.row(i));
        an=aprop(univ,i,avals.row(i),gmom.row(i));
        sn=sprop(univ,i,svals.row(i))

        for j in range(ntraj):
            qvalsn.set(i,j,qn.get(j))
            pvalsn.set(i,j,pn.get(j))
            avalsn.set(i,j,an.get(j))
            svalsn.set(i,j,sn.get(j))

    qpas1n=[MATRIX(qvalsn),MATRIX(pvalsn),MATRIX(avalsn),MATRIX(svalsn)]
    st=qtag_calc.overlap(ntraj,qpas1n,qpas1)  ## AVA: should this be time_overlap maybe?
    b1=st*c1_new

    if norm2 < decpl:
        qvals,pvals=MATRIX(qpas2[0]),MATRIX(qpas2[1])
        qvals_cls,pvals_cls=cls_prop(univ,qvals,pvals,model_params)
        qpas2n=[MATRIX(qvals_cls),MATRIX(pvals_cls),MATRIX(avalsn),MATRIX(svalsn)]
    else:
        qvals,pvals,avals,svals=MATRIX(qpas2[0]),MATRIX(qpas2[1]),MATRIX(qpas2[2]),MATRIX(qpas2[3])
        mom,r,gmom,gr=mom_calc(univ,beta,qpas2,c2_new)

        for i in range(ndof):
            qn2=qprop(univ,i,qvals.row(i),mom.row(i));
            pn2=pprop(mom.row(i));
            an2=aprop(univ,i,avals.row(i),gmom.row(i));
            sn2=sprop(univ,i,svals.row(i))

            for j in range(ntraj):
                qvalsn.set(i,j,qn2.get(i,j))
                pvalsn.set(i,j,pn2.get(i,j))
                avalsn.set(i,j,an2.get(i,j))
                svalsn.set(i,j,sn2.get(i,j))

        qpas2n=[MATRIX(qvalsn),MATRIX(pvalsn),MATRIX(avalsn),MATRIX(svalsn)]

    st=qtag_calc.overlap(ntraj,qpas2n,qpas2)  ### AVA: again - maybe time-overlap?
    b2=st*c2_new

    return(qpas1n,qpas2n,b1,b2)

def fixed():
    return ()

def mean_field(univ,mss,mom_calc,props,model_params,qpas1,c1_new,qpas2,c2_new,norm2,beta):
    """Returns the values for the new basis parameter matrices on surfaces 1 (*qpas1n*) and 2 (*qpas2n*), 
       as well as their corresponding projection vectors *b1* and *b2*, where the basis parameters are calculated 
       according to an "average" surface, defined by the average momentum calculed in the *mom_avg* function. 
       Also necessary are the functions for calculating momentum (*mom_calc*, although specifically mom_avg here) 
       and basis updates (*props*).

    Args:
        univ (dictionary): Dictionary containing various system parameters.

        mss (dictionary): Dictionary containing multi-surface scheme parameters.

        mom_calc (function object): The function object needed to calculate the momentum, as defined by qtag_config.

        props (list): A list of function objects for propagating basis parameters {q,p,a,s}.

        model_params (dictionary): Dictionary containing potential parameters.

        qpas1 (list): List of {q,p,a,s} MATRIX objects for surface 1.

        c1_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 1.

        qpas2 (list): List of {q,p,a,s} MATRIX objects for surface 2.

        c2_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 2.

        norm2 (float): The population on surface 2

        beta (float): Parameter determining tolerance in the momentum linear fitting algorithm.

    Returns:
        qpas1n (list): List of updated {q,p,a,s} MATRIX objects for surface 1.

        qpas2n (list): List of updated {q,p,a,s} MATRIX objects for surface 2.

        b1 (CMATRIX): The updated projection vector of the old basis (defined by qpas1) onto the new basis (defined by qpas1n).

        b2 (CMATRIX): The updated projection vector of the old basis (defined by qpas2) onto the new basis (defined by qpas2n).
    """

    ndof,ntraj=univ['ndof'],univ['ntraj']
    decpl=mss['decpl']
    qvalsn,pvalsn,avalsn,svalsn=MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj)
    qprop=props[0];pprop=props[1];aprop=props[2];sprop=props[3]
    qvals,pvals,avals,svals=MATRIX(qpas1[0]),MATRIX(qpas1[1]),MATRIX(qpas1[2]),MATRIX(qpas1[3])

    mom,r,gmom,gr=mom_calc(univ, beta, qpas1, c1_new, qpas2, c2_new)
    for i in range(ndof):
        qn=qprop(univ,i,qvals.row(i),mom.row(i),mss);
        pn=pprop(mom.row(i));
        an=aprop(univ,i,avals.row(i),gmom.row(i),mss);
        sn=sprop(univ,i,svals.row(i))

        for j in range(ntraj):
            qvalsn.set(i,j,qn.get(i,j))
            pvalsn.set(i,j,pn.get(i,j))
            avalsn.set(i,j,an.get(i,j))
            svalsn.set(i,j,sn.get(i,j))

    qpas1n=[MATRIX(qvalsn),MATRIX(pvalsn),MATRIX(avalsn),MATRIX(svalsn)]
    qpas2n=[MATRIX(qvalsn),MATRIX(pvalsn),MATRIX(avalsn),MATRIX(svalsn)]

    st=qtag_calc.overlap(ntraj,qpas1n,qpas1)
    b1=st*c1_new

    st=qtag_calc.overlap(ntraj,qpas2n,qpas2)
    b2=st*c2_new

    return(qpas1n,qpas2n,b1,b2)


def two_surf(univ,mss,mom_calc,props,model_params,qpas1,c1_new,qpas2,c2_new,norm2,beta):
    """Returns the values for the new basis parameter matrices on surfaces 1 (*qpas1n*) and 2 (*qpas2n*), 
       as well as their corresponding projection vectors *b1* and *b2*, where the basis parameters are 
       calculated in pairs, although the surface assignment alternates with each function (i.e. even basis 
       are ground surface, odd basis are excited surface). Also necessary are the functions for calculating 
       momentum (*mom_calc*, although specifically mom_avg here) and basis updates (*props*).

    Args:
        univ (dictionary): Dictionary containing various system parameters.

        mss (dictionary): Dictionary containing multi-surface scheme parameters.

        mom_calc (function object): The function object needed to calculate the momentum, as defined by qtag_config.

        props (list): A list of function objects for propagating basis parameters {q,p,a,s}.

        model_params (dictionary): Dictionary containing potential parameters.

        qpas1 (list): List of {q,p,a,s} MATRIX objects for surface 1.

        c1_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 1.

        qpas2 (list): List of {q,p,a,s} MATRIX objects for surface 2.

        c2_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 2.

        norm2 (float): The population on surface 2

        beta (float): Parameter determining tolerance in the momentum linear fitting algorithm.

    Returns:
        qpas1n (list): List of updated {q,p,a,s} MATRIX objects for surface 1.

        qpas2n (list): List of updated {q,p,a,s} MATRIX objects for surface 2.

        b1 (CMATRIX): The updated projection vector of the old basis (defined by qpas1) onto the new basis (defined by qpas1n).

        b2 (CMATRIX): The updated projection vector of the old basis (defined by qpas2) onto the new basis (defined by qpas2n).

    """

    ndof,ntraj=univ['ndof'],univ['ntraj']
    qprop=props[0];pprop=props[1];aprop=props[2];sprop=props[3]

    qvals1,pvals1,avals1,svals1=MATRIX(qpas1[0]),MATRIX(qpas1[1]),MATRIX(qpas1[2]),MATRIX(qpas1[3])
    qvals2,pvals2,avals2,svals2=MATRIX(qpas2[0]),MATRIX(qpas2[1]),MATRIX(qpas2[2]),MATRIX(qpas2[3])
    qvals1n,pvals1n,avals1n,svals1n=MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj)
    qvals2n,pvals2n,avals2n,svals2n=MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj)
    b1,b2=CMATRIX(ntraj,1),CMATRIX(ntraj,1)
    dummy=CMATRIX(ntraj,ntraj)

    mom1,r1,gmom1,gr1=mom_calc(univ,beta,qpas1,c1_new)
    mom2,r2,gmom2,gr2=mom_calc(univ,beta,qpas2,c2_new)

    for i in range(ndof):
        qn1=qprop(univ,i,qvals1.row(i),mom1.row(i),mss);pn1=pprop(mom1.row(i));an1=aprop(univ,i,avals1.row(i),gmom1.row(i),mss);sn1=sprop(univ,i,svals1.row(i))
        qn2=qprop(univ,i,qvals2.row(i),mom2.row(i),mss);pn2=pprop(mom2.row(i));an2=aprop(univ,i,avals2.row(i),gmom2.row(i),mss);sn2=sprop(univ,i,svals2.row(i))
        for j in range(ntraj):
            if (j%2==0):
                qvals1n.set(i,j,qn1.get(j)); qvals2n.set(i,j,qn1.get(j))
                pvals1n.set(i,j,pn1.get(j)); pvals2n.set(i,j,pn1.get(j))
                avals1n.set(i,j,an1.get(j)); avals2n.set(i,j,an1.get(j))
                svals1n.set(i,j,sn1.get(j)); svals2n.set(i,j,sn1.get(j))
            else:
                qvals1n.set(i,j,qn2.get(j)); qvals2n.set(i,j,qn2.get(j))
                pvals1n.set(i,j,pn2.get(j)); pvals2n.set(i,j,pn2.get(j))
                avals1n.set(i,j,an2.get(j)); avals2n.set(i,j,an2.get(j))
                svals1n.set(i,j,sn2.get(j)); svals2n.set(i,j,sn2.get(j))

    qpas1n=[MATRIX(qvals1n),MATRIX(pvals1n),MATRIX(avals1n),MATRIX(svals1n)]
    qpas2n=[MATRIX(qvals2n),MATRIX(pvals2n),MATRIX(avals2n),MATRIX(svals2n)]

    st=qtag_calc.overlap(ntraj,qpas1n,qpas1)
    b1=st*c1_new

    st=qtag_calc.overlap(ntraj,qpas2n,qpas2)
    b2=st*c2_new
	
    return(qpas1n,qpas2n,b1,b2)
