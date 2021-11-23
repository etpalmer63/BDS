#include <AMReX_Array.H>





using namespace amrex;

// figure out ml_layout data

void function bds (mla, MultiFab s, MultiFab s_new, MultiFab umac, Real dx, Real dt, is_converv)
{

    //Array4


}


















//Does this only work for 3d? 

Real function eval (Real s,
                    const Array1D<const Real, 1, 7> slope,     //why 7 
                    const Array1D<const Real, 1, 3> del)       //3 - space_dim
{

    val = s + del(1)*slope(1)
            + del(2)*slope(2)
            + del(3)*slope(3)
            + del(1)*del(2)*slope(4)
            + del(1)*del(3)*slope(5)
            + del(2)*del(3)*slope(6)
            + del(1)*del(2)*del(3)*slope(7);

    return val;

}

