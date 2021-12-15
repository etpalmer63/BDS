#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>

#include "bds.H"

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

    //std::array<MultiFab, AMREX_SPACEDIM> umac_mf(ba, dm, 3, Nghost);
    Array<MultiFab, AMREX_SPACEDIM> umac_mf{AMREX_D_DECL(MultiFab(convert(ba,IntVect::TheDimensionVector(0)), dm, 3, Nghost),
                                                              MultiFab(convert(ba,IntVect::TheDimensionVector(1)), dm, 3, Nghost),
                                                              MultiFab(convert(ba,IntVect::TheDimensionVector(2)), dm, 3, Nghost))};
    //MultiFab uVel(ba, dm, Ncomp, Nghost);
    //MultiFab vVel(ba, dm, Ncomp, Nghost);
    //MultiFab wVel(ba, dm, Ncomp, Nghost);

    // Track Divergence
    //MultiFab umac_div(ba, dm, 3, Nghost);


    // time = starting time in the simulation
    Real time = 0.0;

    // **********************************
    // INITIALIZE DATA

            // Save this for later
            // 2D
            //umac(i,j,k) = std::tanh(amrex::Math::abs(0.15 - std::sqrt( std::pow(y-0.5,2) ))/0.5);
            //vmac(i,j,k) = 0.25;
            // 3D
            //umac(i,j,k) = std::tanh(amrex::Math::abs(0.15 - std::sqrt( std::pow(y-0.5,2) + std::pow(z-0.5,2) ))/0.333);
            //vmac(i,j,k) = 0.25;
            //wmac(i,j,k) = 0.05*std::exp(-15*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));


    AMREX_D_DECL( umac_mf[0].setVal(1.0),
                  umac_mf[1].setVal(0.5),
                  umac_mf[2].setVal(0.25) );

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

//        Array4< Real> const& umac = umac_mf[0].array(mfi);  //HACK -- later may need to adjusted.
//        Array4< Real> const& vmac = umac_mf[1].array(mfi);
//#if (AMREX_SPACEDIM == 3)
//        Array4< Real> const& wmac = umac_mf[2].array(mfi);
//#endif
        Array4< Real> const& S_old = s_old_mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {

#if (AMREX_SPACEDIM == 2)

            Real x = (i+0.5) * dx[0];
            Real y = (j+0.5) * dx[1];

            Real r = std::sqrt(std::pow(x-0.375,2) + std::pow(y-0.5,2));

#elif (AMREX_SPACEDIM == 3)

            Real x = (i+0.5) * dx[0];
            Real y = (j+0.5) * dx[1];
            Real z = (k+0.5) * dx[2];

            Real r = std::sqrt(std::pow(x-0.375,2) + std::pow(y-0.5,2) + std::pow(z-0.5,2));

#endif
            if ( r <= 0.1 ) {
              S_old(i,j,k) = 1.0;
            } else {
              S_old(i,j,k) = 0.0;
            }

        });
    }



    // Write a plotfile of the initial data if plot_int > 0
    if (plot_int > 0)
    {
        int step = 0;
        const std::string& pltfile = amrex::Concatenate("plt",step,5);
        WriteSingleLevelPlotfile(pltfile, s_old_mf, {"S"}, geom, time, 0);
        Print() << "Time: " << time << " Sum of S: " << s_old_mf.sum() << std::endl;
    }

    int comp = 0; //HACK figure out what to do with this later

    constexpr Real EndTime = 1.0;

    for (int step = 1; step <= nsteps; ++step)
    {
        // Land at exactly 1 sec.
        if ( time + dt > EndTime) { dt = EndTime - time; }



        // fill periodic ghost cells
        s_old_mf.FillBoundary(geom.periodicity());

        bds(s_old_mf, geom, s_new_mf, umac_mf, dt, comp, is_conserv);
        // new_phi = old_phi + dt * Laplacian(old_phi)
        // loop over boxes
    /*
        for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            const Array4<Real>& phiOld = phi_old.array(mfi);
            const Array4<Real>& phiNew = phi_new.array(mfi);


            //Call BDS Fortran routine






            // advance the data by dt
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                phiNew(i,j,k) = phiOld(i,j,k) + dt *
                    ( (phiOld(i+1,j,k) - 2.*phiOld(i,j,k) + phiOld(i-1,j,k)) / (dx[0]*dx[0])
                     +(phiOld(i,j+1,k) - 2.*phiOld(i,j,k) + phiOld(i,j-1,k)) / (dx[1]*dx[1])
#if (AMREX_SPACEDIM == 3)
                     +(phiOld(i,j,k+1) - 2.*phiOld(i,j,k) + phiOld(i,j,k-1)) / (dx[2]*dx[2])
#endif
                        );
            });
        }
*/
        // update time
        time = time + dt;

        // copy new solution into old solution
        MultiFab::Copy(s_old_mf, s_new_mf, 0, 0, 1, 0);

        // calculate norm0 of divergence
        //computeDivergence(umac_div, umac_mf, geom);

        // Tell the I/O Processor to write out which step we're doing
        //amrex::Print() << "Advanced step " << step << "\n";

        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if (plot_int > 0 && step%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",step,5);
            WriteSingleLevelPlotfile(pltfile, s_new_mf, {"S"}, geom, time, step);
        }

        printf("Time: %.6f Sum of S: %.6f\n", time, s_old_mf.sum());
        //Print() << "Time: " << time << " Sum of S: " << s_old_mf.sum() << std::endl;

        if ( time >= EndTime ) { return; }
    }
}
