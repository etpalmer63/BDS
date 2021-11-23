
using namespace amres;

/**
 * Compute sedgex without transverse corrections.
 *
 */

Real function sedgexWOTransCorr (int i, int j, int k,
                Array3D<Real, 1, > uadv,
                Real isign,
                int ioff,
                Real dt,
                Array1D<Real, 1,3> del,
                Array3D<Real, 1, > sedgex,
                Real hx,
                Array4D<Real,    > slope)
{

    if (uadv(i,j,k) > 0.0)
    {
        isign = 1.0;
        ioff = -1;
    } else {
        isign = -1.0;
        ioff = 0;
    }

    // centroid of rectangular volume
    del(1) = isign*0.5*hx - 0.5*uadv(i,j,k)*dt;
    del(2) = 0.0;
    del(3) = 0.0;
    
    sedgex(i,j,k) = eval(s(i+ioff,j,k),slope(i+ioff,j,k,:),del);

    // source term
    return sedgex(i,j,k) - dt2*sedgex(i,j,k)*ux(i+ioff,j,k);
    //sedgex(i,j,k) = sedgex(i,j,k) - dt2*sedgex(i,j,k)*ux(i+ioff,j,k)

}




/**
 * Compue \Gamma^{y+} without corner corrections
 *
 */

Real function gammaYpWOCrnCorr ( int i, int j, int k,
                Array3D<Real, 1, > uadv,
                Array3D<Real, 1, > vadv,
                Real isign, Real jsign,
                int ioff, int joff,
                Real dt,
                Real dt3,
                Array1D<Real, 1,3> p1,
                Array1D<Real, 1,3> p2,
                Array1D<Real, 1,3> p3,
                Array1D<Real, 1,3> del,
                Array3D<Real, 1, > sedgex,
                Real hx, Real hy,
                Array4D<Real,    > slope)

{
    Real val1;
    Real val2;
    Real val3;

    Real gamma;

    if (vadv(i+ioff,j+1,k) < 0.0)
    {
       jsign = 1.0; 
       joff = 0; 
    } else {
       jsign = -1.0; 
       joff = 1; 
    }
    
    u = 0.0; 
    if (uadv(i,j,k)*uadv(i,j+joff,k) .gt. 0) then
       u = uadv(i,j+joff,k); 
    endif
    
    p1(1) = isign*0.5*hx; 
    p1(2) = jsign*0.5*hy; 
    p1(3) = 0.0; 
    
    p2(1) = isign*0.5*hx - uadv(i,j,k)*dt; 
    p2(2) = jsign*0.5*hy; 
    p2(3) = 0.0; 
    
    p3(1) = isign*0.5*hx - u*dt; 
    p3(2) = jsign*0.5*hy - vadv(i+ioff,j+1,k)*dt; 
    p3(3) = 0.0; 
    
    //do ll=1,3
    for (int ll=1; ll<=3; ++ll)
    {
       del(ll) = (p2(ll)+p3(ll))/2.0; 
    }
    //call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)
    val1 = eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del);
    
    //do ll=1,3
    for (int ll=1; ll<=3; ++ll)
    {
       del(ll) = (p1(ll)+p3(ll))/2.0; 
    } 
    //call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)
    val2 = eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del);
    
    for (int ll=1; ll<=3; ++ll)
    {
       del(ll) = (p1(ll)+p2(ll))/2.0; 
    } 

    //call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)
    val3 = eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del);
    
    // average these centroid values to get the average value
    gamma = (val1+val2+val3)/3.0; 
    
    // source term
    gamma = gamma - dt3*(gamma*ux(i+ioff,j+joff,k) + gamma*vy(i+ioff,j+joff,k)); 

    return gamma;
}


/**
 * Correct \Gamma^{y+} with \Gamma^{y+,z+}
 */

Real function corrGammaYpWGamaYpZp (int i, int j, int k,
                int joff, int koff,
                Real jsign, Real ksign,
                Array3D<Real, 1, > uadv,
                Array3D<Real, 1, > vadv,
                Real dt,
                Real dt3,
                Array1D<Real, 1,3> p1,
                Array1D<Real, 1,3> p2,
                Array1D<Real, 1,3> p3,
                Array1D<Real, 1,3> del,
                Array3D<Real, 1, > sedgex,
                Real hx, Real hy,
                Array4D<Real,    > slope)

{

    if (wadv(i+ioff,j+joff,k+1) .gt. 0) then
       ksign = 1.d0
       koff = 0
    else
       ksign = -1.d0
       koff = 1
    endif
    
    uu = 0.d0
    if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
       uu = uadv(i,j+joff,k+koff)
    endif
    
    vv = 0.d0
    if (vadv(i+ioff,j+1,k)*vadv(i+ioff,j+1,k+koff) .gt. 0) then
       vv = vadv(i+ioff,j+1,k+koff)
    endif
    
    p1(1) = isign*0.5d0*hx
    p1(2) = jsign*0.5d0*hy
    p1(3) = ksign*0.5d0*hz
    
    p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
    p2(2) = jsign*0.5d0*hy
    p2(3) = ksign*0.5d0*hz
    
    p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
    p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k)*dt
    p3(3) = ksign*0.5d0*hz
    
    p4(1) = isign*0.5d0*hx - uu*dt
    p4(2) = jsign*0.5d0*hy - vv*dt
    p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k+1)*dt
    
    do ll=1,3
       del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
    end do
    call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)
    
    do ll=1,3
       del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
    end do
    call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)
    
    do ll=1,3
       del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
    end do
    call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)
    
    do ll=1,3
       del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
    end do
    call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)
    
    do ll=1,3
       del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
    end do
    call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)
    
    gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)
    
    ! source term
    gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                             +gamma2*vy(i+ioff,j+joff,k+koff) &
                             +gamma2*wz(i+ioff,j+joff,k+koff))
    
    gamma2 = gamma2 * wadv(i+ioff,j+joff,k+1)
    
    gamma = gamma - dt*gamma2/(3.d0*hz)


}
