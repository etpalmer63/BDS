#ifndef BDS_FUNCS_
#define BDS_FUNCS_

#include<AMReX_REAL.H>
#include<AMReX_MultiFab.H>
#include<AMReX_Array.H>

//using namespace amrex; //HACK --this has to go

void bds (amrex::BoxArray& ba,
          const amrex::Geometry& geom,
          const amrex::DistributionMapping& dmap,
          int nlevs,
          const amrex::MultiFab& s_mf,
          amrex::MultiFab& sn_mf,
          const amrex::MultiFab& umac_mf,
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx, amrex::Real dt,
          amrex::Array1D<const bool, 1, AMREX_SPACEDIM> const& is_conserv);

void bdsslope_2d (amrex::Array1D<const int > const& lo,
                  amrex::Array1D<const int> const& hi,
                  amrex::Array2D<const amrex::Real> const& s,      const int ng_s,
                  amrex::Array3D<      amrex::Real> const& slope,  const int ng_c,
                  amrex::Array1D<const amrex::Real> const& dx);

void bdsslope_3d (amrex::GpuArray<amrex::Real, AMREX_SPACEDIM > const& lo,
                  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM > const& hi,
                  const amrex::MultiFab& s_mf,     const int ng_s,
                  const amrex::MultiFab& slope_mf, const int ng_c,
                  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx, 
                  amrex::BoxArray& ba, 
                  amrex::Geometry& geom, 
                  amrex::DistributionMapping& dmap);


void bdsconc_2d (amrex::Array1D<const int > const& lo,
                 amrex::Array1D<const int > const& hi,
                 amrex::Array2D<const amrex::Real> const& s,
                 amrex::Array2D<      amrex::Real> const& sn,      const int ng_s,
                 amrex::Array3D<const amrex::Real> const& slope,   const int ng_c,
                 amrex::Array2D<const amrex::Real> const& uadv,
                 amrex::Array2D<const amrex::Real> const& vadv,    const int ng_u,
                 amrex::Array1D<const int > const& dx,      const amrex::Real dt,
                 const bool is_conserv);

void bdsconc_3d (amrex::Array1D<const int > const& lo, 
                 amrex::Array1D<const int > const& hi, 
                 amrex::Array3D<const amrex::Real> const& s,
                 amrex::Array3D<      amrex::Real> const& sn,      const int ng_s,
                 amrex::Array4< const amrex::Real> const& slope,   const int ng_c,
                 amrex::Array3D<const amrex::Real> const& uadv,
                 amrex::Array3D<const amrex::Real> const& vadv,
                 amrex::Array3D<const amrex::Real> const& wadv,    const int ng_u,
                 amrex::Array1D<const int > const& dx,      const amrex::Real dt,
                 const bool is_conserv);

Real eval ( const amrex::Real s,
            amrex::Array1D<const Real> const& slope,
            amrex::Array1D<const Real> const& del );

#endif
