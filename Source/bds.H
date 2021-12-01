#ifndef BDS_FUNCS_
#define BDS_FUNCS_

#include<AMReX_REAL.H>
#include<AMReX_MultiFab.H>
#include<AMReX_Array.H>



void bds (const amrex::MultiFab& s_mf,
          const amrex::Geometry& geom,      
          amrex::MultiFab& sn_mf,
          std::array<amrex::MultiFab, AMREX_SPACEDIM>& umac_mf,
          amrex::Real dt,
          int comp,
          int is_conserv);

void bdsslope (const amrex::MultiFab& s_mf,   
               const amrex::Geometry& geom,
               amrex::MultiFab& slope_mf,
               int comp);

#if (AMREX_SPACEDIM == 2)
void bdsconc_2d (const amrex::MultiFab& s_mf,   
                 const amrex::Geometry& geom,
                 amrex::MultiFab& sn_mf,
                 const amrex::MultiFab& slope_mf,
                 const std::array<amrex::MultiFab, AMREX_SPACEDIM>& umac_mf,
                 const amrex::Real dt,
                 int comp,
                 int is_conserv);

#elif (AMREX_SPACEDIM == 3)
void bdsconc_3d (const amrex::MultiFab& s_mf,   
                 const amrex::Geometry& geom,
                 amrex::MultiFab& sn_mf,
                 const amrex::MultiFab& slope_mf,
                 const std::array<amrex::MultiFab, AMREX_SPACEDIM>& umac_mf,
                 const amrex::Real dt,
                 int comp,
                 int is_conserv);
#endif

amrex::Real eval ( const amrex::Real s,
                   amrex::Array1D<amrex::Real,1,7>& slope,
                   amrex::Array1D<amrex::Real,1,3>& del );

#endif
