#include <PeleLMDeriveFunc.H>
#include <PeleLM_Index.H>
#include <PelePhysics.H>
#include <mechanism.H>
#include <PeleLM.H>

using namespace amrex;

//
// Extract temp
//
void pelelm_dertemp (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                     const FArrayBox& statefab, const FArrayBox& /*pressfab*/,
                     const Geometry& /*geomdata*/,
                     Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    auto const in_dat = statefab.array();
    auto       der = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) = in_dat(i,j,k,TEMP);
    });
}

//
// Extract species mass fractions Y_n
//
void pelelm_dermassfrac (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                         const FArrayBox& statefab, const FArrayBox& /*pressfab*/,
                         const Geometry& /*geomdata*/,
                         Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(statefab.nComp() >= NUM_SPECIES+1);
    AMREX_ASSERT(ncomp == NUM_SPECIES);
    auto const in_dat = statefab.array();
    auto       der = derfab.array(dcomp);
    amrex::ParallelFor(bx, NUM_SPECIES,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        amrex::Real rhoinv = 1.0 / in_dat(i,j,k,DENSITY);
        der(i,j,k,n) = in_dat(i,j,k,FIRSTSPEC+n) * rhoinv;
    });
}

//
// Compute cell averaged pressure from nodes
//
void pelelm_deravgpress (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                         const FArrayBox& /*statefab*/, const FArrayBox& pressfab,
                         const Geometry& /*geomdata*/,
                         Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    auto const in_dat = pressfab.array();
    auto       der = derfab.array(dcomp);
    Real factor = 1.0 / ( AMREX_D_TERM(2.0,*2.0,*2.0) );
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) =  factor * (  in_dat(i+1,j,k)     + in_dat(i,j,k)
#if (AMREX_SPACEDIM >= 2 )
                                + in_dat(i+1,j+1,k)   + in_dat(i,j+1,k)
#if (AMREX_SPACEDIM == 3 )
                                + in_dat(i+1,j,k+1)   + in_dat(i,j,k+1)
                                + in_dat(i+1,j+1,k+1) + in_dat(i,j+1,k+1)
#endif
#endif
                                );
    });
}

//
// Compute vorticity magnitude
//
void pelelm_dermgvort (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                       const FArrayBox& statefab, const FArrayBox& /*pressfab*/,
                       const Geometry& geomdata,
                       Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{

    AMREX_D_TERM(const amrex::Real idx = geomdata.InvCellSize(0);,
                 const amrex::Real idy = geomdata.InvCellSize(1);,
                 const amrex::Real idz = geomdata.InvCellSize(2););

    auto const& dat_arr = statefab.const_array();
    auto const&vort_arr = derfab.array(dcomp);

    // TODO : EB
    // TODO : BCs

    {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if ( AMREX_SPACEDIM == 2 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            vort_arr(i,j,k) = vx-uy;

#elif ( AMREX_SPACEDIM == 3 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;

            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            amrex::Real wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;

            amrex::Real uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
            amrex::Real vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;

            vort_arr(i,j,k) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
#endif
        });
    }

}

//
// Compute the kinetic energy
//
void pelelm_derkineticenergy (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                              const FArrayBox& statefab, const FArrayBox& /*pressfab*/,
                              const Geometry& /*geomdata*/,
                              Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    auto const rho = statefab.array(DENSITY);
    auto const vel = statefab.array(VELX);
    auto       der = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) = 0.5 * rho(i,j,k)
                         * ( AMREX_D_TERM(  vel(i,j,k,0)*vel(i,j,k,0),
                                          + vel(i,j,k,1)*vel(i,j,k,1),
                                          + vel(i,j,k,2)*vel(i,j,k,2)) );
    });
}

//
// Compute mixture fraction
//
void pelelm_dermixfrac (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                        const FArrayBox& statefab, const FArrayBox& /*pressfab*/,
                        const Geometry& /*geomdata*/,
                        Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(ncomp == 1);

    if (a_pelelm->Zfu < 0.0) amrex::Abort("Mixture fraction not initialized");

    auto const density   = statefab.array(DENSITY);
    auto const rhoY      = statefab.array(FIRSTSPEC);
    auto       mixt_frac = derfab.array(dcomp);

    amrex::Real Zox_lcl   = a_pelelm->Zox;
    amrex::Real Zfu_lcl   = a_pelelm->Zfu;
    amrex::Real denom_inv = 1.0 / (Zfu_lcl - Zox_lcl);
    amrex::GpuArray<amrex::Real,NUM_SPECIES> fact_Bilger;
    for (int n=0; n<NUM_SPECIES; ++n) {
        fact_Bilger[n] = a_pelelm->spec_Bilger_fact[n];
    }

    amrex::ParallelFor(bx,
    [density, rhoY, mixt_frac, fact_Bilger, Zox_lcl, denom_inv] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real rho_inv = 1.0_rt / density(i,j,k);
        mixt_frac(i,j,k) = 0.0_rt;
        for (int n = 0; n<NUM_SPECIES; ++n) {
            mixt_frac(i,j,k) += ( rhoY(i,j,k,n) * fact_Bilger[n] ) * rho_inv;
        }
        mixt_frac(i,j,k) = ( mixt_frac(i,j,k) - Zox_lcl ) * denom_inv;
    });
}

//
// Compute progress variable
//
void pelelm_derprogvar (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                        const FArrayBox& statefab, const FArrayBox& /*pressfab*/,
                        const Geometry& /*geomdata*/,
                        Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(ncomp == 1); 

    if (a_pelelm->m_C0 < 0.0) amrex::Abort("Progress variable not initialized");

    auto const density  = statefab.array(DENSITY);
    auto const rhoY     = statefab.array(FIRSTSPEC);
    auto const temp     = statefab.array(TEMP);
    auto       prog_var = derfab.array(dcomp);

    amrex::Real C0_lcl   = a_pelelm->m_C0;
    amrex::Real C1_lcl   = a_pelelm->m_C1;
    amrex::Real denom_inv = 1.0 / (C1_lcl - C0_lcl);
    amrex::GpuArray<amrex::Real,NUM_SPECIES+1> Cweights;
    for (int n=0; n<NUM_SPECIES+1; ++n) {
        Cweights[n] = a_pelelm->m_Cweights[n];
    }   

    amrex::ParallelFor(bx, [=,revert=a_pelelm->m_Crevert]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {   
        amrex::Real rho_inv = 1.0_rt / density(i,j,k);
        prog_var(i,j,k) = 0.0_rt;
        for (int n = 0; n<NUM_SPECIES; ++n) {
            prog_var(i,j,k) += ( rhoY(i,j,k,n) * Cweights[n] ) * rho_inv;
        }   
        prog_var(i,j,k) += temp(i,j,k) * Cweights[NUM_SPECIES];
        if (revert) {
           prog_var(i,j,k) = 1.0 - ( prog_var(i,j,k) - C0_lcl ) * denom_inv;
        } else {
           prog_var(i,j,k) = ( prog_var(i,j,k) - C0_lcl ) * denom_inv;
        }   
    });
}

