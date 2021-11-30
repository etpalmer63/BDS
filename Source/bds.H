#ifndef BDS_FUNCS_
#define BDS_FUNCS_

#include<AMReX_REAL.H>
#include<AMReX_MultiFab.H>
#include<AMReX_Array.H>



void bds (const amrex::MultiFab& s_mf,
          const amrex::Geometry& geom,      
          const amrex::MultiFab& sn_mf,
          std::array<amrex::MultiFab, AMREX_SPACEDIM>& umac_mf,
          amrex::Real dt,
          int comp,
          int is_conserv);

void bdsslope (const amrex::MultiFab& s_mf,   
               const amrex::Geometry& geom,
               const amrex::MultiFab& slope_mf,
               int comp);

/*
void bdsconc_2d (amrex::Array1D<const int > const& lo,
                 amrex::Array1D<const int > const& hi,
                 amrex::Array2D<const amrex::Real> const& s,
                 amrex::Array2D<      amrex::Real> const& sn,      const int ng_s,
                 amrex::Array3D<const amrex::Real> const& slope,   const int ng_c,
                 amrex::Array2D<const amrex::Real> const& uadv,
                 amrex::Array2D<const amrex::Real> const& vadv,    const int ng_u,
                 amrex::Array1D<const int > const& dx,      const amrex::Real dt,
                 const bool is_conserv);
*/

void bdsconc_3d (const amrex::MultiFab& s_mf,   
                 const amrex::Geometry& geom,
                 const amrex::MultiFab& sn_mf,
                 amrex::MultiFab& slope_mf,
                 amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& umac_mf,
                 const amrex::Real dt,
                 int comp,
                 int is_conserv);

amrex::Real eval ( const amrex::Real s,
                   amrex::Array1D<amrex::Real,1,7>& slope,
                   amrex::Array1D<amrex::Real,1,3>& del );

#endif
