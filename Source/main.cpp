#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <math.h>

#include "bds.H"

constexpr amrex::Real PI = 3.141592653589793238;


using namespace amrex;
void main_main();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();
    return 0;
}

void main_main ()
{

    // **********************************
    // SIMULATION PARAMETERS

    // number of cells on each side of the domain
    int n_cell;

    // size of each box (or grid)
    int max_grid_size;

    // total steps in simulation
    int nsteps;

    // how often to write a plotfile
    int plot_int;

    // time step
    Real dt;

    // input of constant velocities
    Real u_val;
    Real v_val;
    Real w_val;

    int is_conserv = 0; //default to non-conservative update
    //Array1D<const bool, 1, AMREX_SPACEDIM> is_conserv = {true, true, true};

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        // pp.get means we require the inputs file to have it
        // pp.query means we optionally need the inputs file to have it - but we must supply a default here
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of
        //   a square (or cubic) domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size",max_grid_size);

        // Default nsteps to 10, allow us to set it to something else in the inputs file
        nsteps = 10;

        pp.query("nsteps",nsteps);

        // Default plot_int to -1, allow us to set it to something else in the inputs file
        //  If plot_int < 0 then no plot files will be written
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // time step
        pp.get("dt",dt);

        pp.get("is_conserv", is_conserv);

        // Inputs, S old, S new, U mac, dx, dt, is_conservative

        pp.get("u_val", u_val);
        pp.get("v_val", v_val);
        pp.get("w_val", w_val);

    }

    // **********************************
    // SIMULATION SETUP

    // make BoxArray and Geometry
    // ba will contain a list of boxes that cover the domain
    // geom contains information such as the physical domain size,
    //               number of points in the domain, and periodicity
    BoxArray ba;
    Geometry geom;

    // AMREX_D_DECL means "do the first X of these, where X is the dimensionality of the simulation"
    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // This defines the physical box, [0,1] in each direction.
    RealBox real_box({AMREX_D_DECL( 0.0, 0.0, 0.0)},
                     {AMREX_D_DECL( 1.0, 1.0, 1.0)});

    // periodic in all direction
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

    // This defines a Geometry object
    geom.define(domain, real_box, CoordSys::cartesian, is_periodic);

    // extract dx from the geometry object
    GpuArray<Real,AMREX_SPACEDIM> const dx = geom.CellSizeArray();

    // Nghost = number of ghost cells for each array
    int Nghost = 3;

    // Ncomp = number of components for each array
    int Ncomp = 1;

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate two phi multifabs; one will store the old state, the other the new.


    //initialize data on the multifab
    MultiFab s_old_mf(ba, dm, Ncomp, Nghost);
    MultiFab s_new_mf(ba, dm, Ncomp, Nghost);

    std::array<MultiFab, AMREX_SPACEDIM> umac_mf{AMREX_D_DECL(MultiFab(convert(ba,IntVect::TheDimensionVector(0)), dm, 3, Nghost),
                                                              MultiFab(convert(ba,IntVect::TheDimensionVector(1)), dm, 3, Nghost),
                                                              MultiFab(convert(ba,IntVect::TheDimensionVector(2)), dm, 3, Nghost))};
    
    // time = starting time in the simulation
    Real time = 0.0;

    // **********************************
    // INITIALIZE DATA

            //Cylindrical shear layer -- save for later
            //umac(i,j,k) = std::tanh(amrex::Math::abs(0.15 - std::sqrt( std::pow(y-0.5,2) ))/0.5);
            //vmac(i,j,k) = 0.25;
            //
            //Cylindrical shear layer
            //umac(i,j,k) = std::tanh(amrex::Math::abs(0.15 - std::sqrt( std::pow(y-0.5,2) + std::pow(z-0.5,2) ))/0.333);
            //vmac(i,j,k) = 0.25;
            //wmac(i,j,k) = 0.05*std::exp(-15*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));



    //umac is constant even in variable velocity scheme
    //umac_mf[0].setVal(1.0);
    //umac_mf[0].setVal(u_val);
    //umac_mf[1].setVal(v_val);
    //umac_mf[2].setVal(w_val);

    for (MFIter mfi(umac_mf[0]); mfi.isValid(); ++mfi){

        const Box& bx = mfi.grownnodaltilebox(0,Nghost);
        //const Box& bx = mfi.validbox(); // HACK -- why didn't this work
        Array4<Real> const& umac = umac_mf[0].array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Real y = (j+0.5) * dx[1];
            Real z = (k+0.5) * dx[2];
            umac(i,j,k) = tanh((0.15 - std::sqrt( std::pow(y-0.5,2) + std::pow(z-0.5,2)))/0.333);

        });
    }

#if (AMREX_SPACEDIM > 1)
    umac_mf[1].setVal(0.25);

    //for (MFIter mfi(umac_mf[1]); mfi.isValid(); ++mfi){

    //    const Box& bx = mfi.grownnodaltilebox(1,Nghost);
    //    //const Box& bx = mfi.validbox(); // HACK -- why didn't this work
    //    Array4<Real> const& vmac = umac_mf[1].array(mfi);

    //    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    //    {
    //        Real x = (i+0.5) * dx[0];
    //        vmac(i,j,k) = 0.5 + 0.5  * std::sin(2*PI*x);
    //    });
    //}

#endif
#if (AMREX_SPACEDIM > 2)

    for (MFIter mfi(umac_mf[2]); mfi.isValid(); ++mfi){

        const Box& bx = mfi.grownnodaltilebox(2,Nghost);
        //const Box& bx = mfi.validbox(); // HACK -- why didn't this work
        Array4<Real> const& wmac = umac_mf[2].array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Real x = (i+0.5) * dx[0];
            Real y = (j+0.5) * dx[1];
            //wmac(i,j,k) = 0.25 + 0.25 * std::cos(2*PI*x);
            wmac(i,j,k) = 0.05 * std::exp(-15*(std::pow(x-0.5,2) + std::pow(y-0.5,2)));
        });
    }

#endif

    // Constant velocity in all directions
    //AMREX_D_DECL( umac_mf[0].setVal(1.0),
    //              umac_mf[1].setVal(0.5),
    //              umac_mf[2].setVal(0.25) );


    AMREX_D_DECL( VisMF::Write(umac_mf[0], "umac.mf"),
                  VisMF::Write(umac_mf[1], "vmac.mf"),
                  VisMF::Write(umac_mf[2], "wmac.mf") );


    Real umac_max = amrex::max(AMREX_D_DECL(umac_mf[0].norm0(0,0),
                                            umac_mf[1].norm0(0,0),
                                            umac_mf[2].norm0(0,0)));

    Print() << "umac_max " << umac_max << std::endl;



    // Set time step size dt = dx / (largest velocities) * CFL
    // For constant velocities this is:
    constexpr Real CFLnum = 0.9;
    dt = dx[0] / umac_max * CFLnum;  //assumes dx constant across all dimensions

    Print() << "dt: " << dt << std::endl;


    // loop over boxes
    for (MFIter mfi(s_old_mf); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        Array4< Real> const& S_old = s_old_mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {

#if (AMREX_SPACEDIM == 2)

            Real x = (i+0.5) * dx[0];
            Real y = (j+0.5) * dx[1];

            Real r = std::sqrt(std::pow(x-0.5,2) + std::pow(y-0.5,2));

            //Real r = std::sqrt(std::pow(x-0.5,2) + std::pow(y-1.0,2));

            //if ( r <= 0.1 ) {
            //  S_old(i,j,k) = 1.0;
            //} else {
            //  S_old(i,j,k) = 0.0;
            //}

            //r = std::sqrt(std::pow(x-0.5,2) + std::pow(y-0.0,2));

            //if ( r <= 0.1 ) {
            //  S_old(i,j,k) += 1.0;
            //} else {
            //  S_old(i,j,k) += 0.0;
            //}



#elif (AMREX_SPACEDIM == 3)

            Real x = (i+0.5) * dx[0];
            Real y = (j+0.5) * dx[1];
            Real z = (k+0.5) * dx[2];

            
            Real r = std::sqrt(std::pow(x-0.375,2) + std::pow(y-0.5,2) + std::pow(z-0.5,2));
            //Real r = std::sqrt(std::pow(x-0.5,2) + std::pow(y-0.5,2) + std::pow(z-0.5,2));

            // Exact solution hack -- shperical step
            //Real r = std::sqrt(std::pow(x-0.5,2) + std::pow(y-1.0,2) + std::pow(z-0.75,2));
            //
            //if ( r <= 0.1 ) {
            //  S_old(i,j,k) = 1.0;
            //} else {
            //  S_old(i,j,k) = 0.0;
            //}
            //
            //r = std::sqrt(std::pow(x-0.5,2) + std::pow(y-0.0,2) + std::pow(z-0.75,2));
            //
            //if ( r <= 0.1 ) {
            //  S_old(i,j,k) += 1.0;
            //} else {
            //  S_old(i,j,k) += 0.0;
            //}

            // Exact solution hack -- Gaussian Bump
            //Real r = std::sqrt(std::pow(x-0.5,2) + std::pow(y-1.0,2) + std::pow(z-0.75,2));
            //S_old(i,j,k) = std::exp(-300.0*std::pow(r,2));
            //r = std::sqrt(std::pow(x-0.5,2) + std::pow(y-0.0,2) + std::pow(z-0.75,2));
            //S_old(i,j,k) += std::exp(-300.0*std::pow(r,2));

#endif
            //Gaussian bump
            //S_old(i,j,k) = std::exp(-300.0*std::pow(r,2));
            
            // Spherical step function
Real sum = 0.;
for (int kk=0; kk<10; ++kk) {
for (int jj=0; jj<10; ++jj) {
for (int ii=0; ii<10; ++ii) {
    Real zz = (k + (kk+0.5)/10.) * dx[2];
    Real yy = (j + (jj+0.5)/10.) * dx[1];
    Real xx = (i + (ii+0.5)/10.) * dx[0];

    Real dist = std::sqrt((xx-0.375)*(xx-0.375) + (yy-0.5)*(yy-0.5) + (zz-0.5)*(zz-0.5));

    if (dist < 0.1) {
        sum = sum + 1.0;
    } else if (dist == 0.1) {
        sum = sum + 0.5;
    }

}
}
}

sum /= 1000.;


//		if ( r <= 0.1 ) {
//              S_old(i,j,k) = 1.0;
//            } else {
//              S_old(i,j,k) = 0.0;
//            }
S_old(i,j,k) = sum;
        });
    }

    Real S_max = s_old_mf.norm0(0,0);
    Print() << "S max " << S_max << std::endl;

    Print() << "Write step 0" << std::endl;

    // Write a plotfile of the initial data if plot_int > 0
    if (plot_int > 0)
    {
        int step = 0;
        const std::string& pltfile = amrex::Concatenate("plt",step,5);
        WriteSingleLevelPlotfile(pltfile, s_old_mf, {"S"}, geom, time, 0);

        Print() << std::fixed << std::setprecision(8)
                << "Time: " << std::setw(16) << time 
		<< " Sum of S: " << std::setw(16) << s_old_mf.sum()
          	<< " Max: " << std::setw(16) << s_old_mf.max(0) 
          	<< " Min: " << std::setw(16) << s_old_mf.min(0) << std::endl;
    }



    int comp = 0; //HACK figure out what to do with this later

    constexpr Real EndTime = 1.0;

    Print() << "Begin time step loop" << std::endl;

    for (int step = 1; step <= nsteps; ++step)
    {
        // Land at exactly 1 sec.
        if ( time + dt > EndTime) { dt = EndTime - time; }



        // fill periodic ghost cells
        s_old_mf.FillBoundary(geom.periodicity());

        bds(s_old_mf, geom, s_new_mf, umac_mf, dt, comp, is_conserv);

        // update time
        time = time + dt;

        // copy new solution into old solution
        MultiFab::Copy(s_old_mf, s_new_mf, 0, 0, 1, 0);

        // calculate norm0 of divergence
        //computeDivergence(umac_div, umac_mf, geom);

        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if (plot_int > 0 && step%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",step,5);
            WriteSingleLevelPlotfile(pltfile, s_new_mf, {"S"}, geom, time, step);
        }

        Real S_sum = s_old_mf.sum();
        Print() << std::fixed << std::setprecision(8)
                << "Time: " << std::setw(16) << time 
		<< " Sum of S: " << std::setw(16) << S_sum
          	<< " Max: " << std::setw(16) << s_old_mf.max(0) 
          	<< " Min: " << std::setw(16) << s_old_mf.min(0) << std::endl;

        if ( time >= EndTime ) { return; }
    }
}
