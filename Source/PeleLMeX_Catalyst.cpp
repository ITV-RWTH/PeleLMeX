#ifndef PeleLMeX_Catalyst_H
#define PeleLMeX_Catalyst_H

#include <PeleLMeX.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>
#include <array>
#include <string>
#include <AMReX_Reduce.H> 
#ifdef PELE_USE_CATALYST
#include <catalyst.hpp>
#include <conduit_cpp_to_c.hpp>
#endif

#ifdef PELE_USE_CATALYST
using namespace amrex;
void PeleLM::CatalystInit() {
    amrex::Print() << "Running Catalyst Initialize script...\n";
    ParmParse const pp_catalyst("catalyst");
    std::string scriptPaths;
    std::string implementation {"paraview"};
    std::string searchPaths;
    std::string proxyPaths;
    pp_catalyst.query("script_paths", scriptPaths);
    pp_catalyst.query("implementation", implementation);
    pp_catalyst.query("implementation_search_paths", searchPaths);
    pp_catalyst.query("proxy_paths", proxyPaths);
    next_time = -1;
    steering_dt = (prob_parm->wall_height/(PeleLM::prob_parm->V_mean*1.5))/10;
    phiSteering.push_back(prob_parm->phi);
    phiSteering.push_back(10000);
    phiSteering.push_back(0.006);
    
    conduit::Node node;

    size_t scriptNumber = 0;
    size_t pos = 0;
    std::string subpath;
    while (scriptPaths.find(':') != std::string::npos || scriptPaths.find(';') != std::string::npos) {
        pos = std::min(scriptPaths.find(':'), scriptPaths.find(';'));
        subpath = scriptPaths.substr(0, pos);

        node["catalyst/scripts/script" + std::to_string(scriptNumber)].set_string(subpath);

        scriptNumber++;
        scriptPaths.erase(0, pos + 1);
    }
    // Prevent empty end paths
    if (scriptPaths.length() != 0) {
    //    node["catalyst/scripts/script" + std::to_string(scriptNumber)].set_string(scriptPaths);
     node["catalyst/scripts/script"].set_string(scriptPaths);
    }

    if (do_inSitu_Steering) {
        node["catalyst/proxies/proxy"].set_string(proxyPaths);
        amrex::Print() << "Proxy Path is : " << proxyPaths;
    }
    node["catalyst_load/implementation"].set_string(implementation);
    node["catalyst_load/search_paths/" + implementation].set_string(searchPaths);

    catalyst_status err = catalyst_initialize(conduit::c_node(&node));
    if (err != catalyst_status_ok)
    {
        std::string message = " Error: Failed to initialize Catalyst!\n";
        std::cerr << message << err << std::endl;
        amrex::Print() << message;
        amrex::Abort(message);
    }
}

void PeleLM::AddDummyZAxis (conduit::Node &mesh_data) {
    // We assume that the input node has multiple children named domain_XXXXXX
    // (as produced by amrex::MultiLevelToBlueprint).  Each child is a Blueprint
    // mesh domain. We'll iterate over them:
    const std::vector<std::string> domain_names = mesh_data.child_names();
    for (const auto &dom_name : domain_names)
    {
        conduit::Node &dom = mesh_data[dom_name];

        // If this domain is already 3D, skip. We check if dims/k exists:
        conduit::Node &coords = dom["coordsets/coords"];

        // At this point, we assume it's 2D and we want to add the dummy third dimension
        // 1) Add dims/k = 2
        coords["dims/k"] = 1;

        // 2) Add spacing/dz
        coords["spacing/dz"] = 0;  

        // 3) Add origin/z
        coords["origin/z"] = 0;

        // 4) Update the topology
        if (dom.has_path("topologies/topo/elements/origin"))
        {
            conduit::Node &orig = dom["topologies/topo/elements/origin"];
            // just set k0 = 0
            orig["k0"] = (conduit::int32)0;
        }

        // 5) If nestsets exist, update them to reflect 3D
        if (dom.has_path("nestsets/nest"))
        {
            conduit::Node &nest = dom["nestsets/nest"];
            if (nest.has_path("windows"))
            {
                conduit::Node &windows = nest["windows"];
                const std::vector<std::string> w_names = windows.child_names();
                for (const auto &w : w_names)
                {
                    conduit::Node &wnd = windows[w];
                    // add 'k' = 0, dims/k = 1, ratio/k = 1
                    wnd["origin/k"] = 0;
                    wnd["dims/k"]   = 1;
                    wnd["ratio/k"]  = 1;
                }
            }
        }
        // Fields do not need changes:  Nx*Ny stays Nx*Ny for Nx*Ny*1 3D cells.
    }
}

void PeleLM::EmptyFieldData(const std::vector<std::string>& field_names, conduit::Node& node){
    const std::string domain_name = amrex::Concatenate("domain_", 0, 6);
    conduit::Node &dom = node[domain_name];

    dom["state/domain_id"] = 0;
    dom["state/cycle"] = 0;
    dom["state/time"] = 0;
    dom["state/level"] = 0;

    dom["coordsets/coords/type"] = "uniform";
    dom["coordsets/coords/dims/i"] = 0;
    dom["coordsets/coords/dims/j"] = 0;
    if(BL_SPACEDIM < 3){
        dom["coordsets/coords/dims/k"] = 0;
    }

    dom["coordsets/coords/spacing/dx"] = 0;
    dom["coordsets/coords/spacing/dy"] = 0;
    if(BL_SPACEDIM < 3){
        dom["coordsets/coords/spacing/dz"] = 0;
    }


    dom["coordsets/coords/origin/x"] = 0;
    dom["coordsets/coords/origin/y"] = 0;
    if(BL_SPACEDIM < 3){
        dom["coordsets/coords/origin/z"] = 0;
    }

    dom["topologies/topo/coordset"] = "coords";
    dom["topologies/topo/type"]     = "uniform";
    dom["topologies/topo/elements/origin/i0"] = 0;
    dom["topologies/topo/elements/origin/j0"] = 0;
     if(BL_SPACEDIM < 3){
        dom["topologies/topo/elements/origin/k0"] = 0;
    }

    auto &fields = dom["fields"]; 
    for(auto &fname : field_names){
        auto &fnode = fields[fname];

        fnode["association"] = "element";
        fnode["topology"] = "topo";
        fnode["values"].set(conduit::DataType::c_double(0));
    }
}

void PeleLM::CatalystExecute () {
    amrex::Print() << "Running Catalyst pipeline scripts... \n";
    BL_PROFILE("PeleLM::CatalystExecute()");
    //----------------------------------------------------------------
    // Blueprint : Mesh data
    std::string copy_conduit_node = "copy_conduit_node";
    BL_PROFILE_VAR("copy_conduit_node", copy_conduit_node);
    
    conduit::Node node;
    auto & state = node["catalyst/state"];
    state["timestep"].set(m_nstep);
    state["time"].set(m_cur_time);

    auto& meshChannel = node["catalyst/channels/mesh"];
    meshChannel["type"].set_string("amrmesh");
    auto& meshData = meshChannel["data"];

    //----------------------------------------------------------------
    // Blueprint : MultiFabs
    // Blueprint : find number of components for multifabs
      //----------------------------------------------------------------
  // Number of components
  int ncomp = 0;

  // State
  if (m_incompressible != 0) {
    // Velocity + pressure gradients
    ncomp = 2 * AMREX_SPACEDIM;
  } else {
    // State + pressure gradients
    if (m_plot_grad_p != 0) {
      ncomp = NVAR + AMREX_SPACEDIM;
    } else {
      ncomp = NVAR;
    }
    // Make the plot lighter by dropping species by default
    if (m_plotStateSpec == 0) {
      ncomp -= NUM_SPECIES;
    }
    if (m_has_divu != 0) {
      ncomp += 1;
    }
  }

  // Reactions
  if ((m_do_react != 0) && (m_skipInstantRR == 0) && (m_plot_react != 0)) {
    // Cons Rate
    ncomp += nCompIR();
    // FunctCall
    ncomp += 1;
    // Extras:
    if (m_plotHeatRelease != 0) {
      ncomp += 1;
    }
  }

#ifdef AMREX_USE_EB
  // Include volume fraction in plotfile
  ncomp += 1;
#endif

#ifdef PELE_USE_RADIATION
  if (do_rad_solve) {
    ncomp += 3;
  }
#endif

  // Derive
  int deriveEntryCount = 0;
  for (int ivar = 0; ivar < m_derivePlotVarCount; ivar++) {
    const PeleLMDeriveRec* rec = derive_lst.get(m_derivePlotVars[ivar]);
    deriveEntryCount += rec->numDerive();
  }
  ncomp += deriveEntryCount;
#ifdef PELE_USE_SPRAY
  if (do_spray_particles) {
    ncomp += SprayParticleContainer::NumDeriveVars();
    if (SprayParticleContainer::plot_spray_src) {
      ncomp += AMREX_SPACEDIM + 2 + SPRAY_FUEL_NUM;
    }
  }
#endif

#ifdef PELE_USE_EFIELD
  if (m_do_extraEFdiags) {
    ncomp += NUM_IONS * AMREX_SPACEDIM;
  }
#endif

  if (m_do_les && m_plot_les) {
    ncomp += 1;
  }

    Vector<MultiFab> mf_plt(finest_level + 1);
    Vector<int> level_steps(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        mf_plt[lev].define(grids[lev], dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
        level_steps[lev] = m_nstep;
    }
    // Blueprint : fill the multifabs
for (int lev = 0; lev <= finest_level; ++lev) {
    int cnt = 0;
    if (m_incompressible != 0) {
        MultiFab::Copy(
        mf_plt[lev], m_leveldata_new[lev]->state, 0, cnt, AMREX_SPACEDIM, 0);
        cnt += AMREX_SPACEDIM;
    } else {
        // Velocity and density
        MultiFab::Copy(
        mf_plt[lev], m_leveldata_new[lev]->state, 0, cnt, AMREX_SPACEDIM + 1,
        0);
        cnt += AMREX_SPACEDIM + 1;
        // Species only if requested
        if (m_plotStateSpec != 0) {
            MultiFab::Copy(
                mf_plt[lev], m_leveldata_new[lev]->state, FIRSTSPEC, cnt, NUM_SPECIES,
                0);
            cnt += NUM_SPECIES;
        }
        MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->state, RHOH, cnt, 3, 0);
        cnt += 3;
#ifdef PELE_USE_EFIELD
        MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->state, NE, cnt, 2, 0);
        cnt += 2;
#endif
#ifdef PELE_USE_SOOT
        MultiFab::Copy(
        mf_plt[lev], m_leveldata_new[lev]->state, FIRSTSOOT, cnt, NUMSOOTVAR,
        0);
        cnt += NUMSOOTVAR;
#endif
#ifdef PELE_USE_RADIATION
        if (do_rad_solve) {
            MultiFab::Copy(mf_plt[lev], rad_model->G()[lev], 0, cnt, 1, 0);
            cnt += 1;
            MultiFab::Copy(mf_plt[lev], rad_model->kappa()[lev], 0, cnt, 1, 0);
            cnt += 1;
            MultiFab::Copy(mf_plt[lev], rad_model->emis()[lev], 0, cnt, 1, 0);
            cnt += 1;
        }
#endif
        if (m_has_divu != 0) {
            MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->divu, 0, cnt, 1, 0);
            cnt += 1;
        }
    }
    if (m_plot_grad_p != 0) {
        MultiFab::Copy(
        mf_plt[lev], m_leveldata_new[lev]->gp, 0, cnt, AMREX_SPACEDIM, 0);
        cnt += AMREX_SPACEDIM;
    }

    if ((m_do_react != 0) && (m_skipInstantRR == 0) && (m_plot_react != 0)) {
        MultiFab::Copy(
        mf_plt[lev], m_leveldatareact[lev]->I_R, 0, cnt, nCompIR(), 0);
        cnt += nCompIR();

        MultiFab::Copy(mf_plt[lev], m_leveldatareact[lev]->functC, 0, cnt, 1, 0);
        cnt += 1;

        if (m_plotHeatRelease != 0) {
            std::unique_ptr<MultiFab> mf;
            mf = std::make_unique<MultiFab>(grids[lev], dmap[lev], 1, 0);
            getHeatRelease(lev, mf.get());
            MultiFab::Copy(mf_plt[lev], *mf, 0, cnt, 1, 0);
            cnt += 1;
        }
    }

#ifdef AMREX_USE_EB
    MultiFab::Copy(mf_plt[lev], EBFactory(lev).getVolFrac(), 0, cnt, 1, 0);
    cnt += 1;
#endif

    for (int ivar = 0; ivar < m_derivePlotVarCount; ivar++) {
        std::unique_ptr<MultiFab> mf;
        mf = derive(m_derivePlotVars[ivar], m_cur_time, lev, 0);
        MultiFab::Copy(mf_plt[lev], *mf, 0, cnt, mf->nComp(), 0);
        cnt += mf->nComp();
    }
#ifdef PELE_USE_SPRAY
    if (SprayParticleContainer::NumDeriveVars() > 0) {
        const int num_spray_derive = SprayParticleContainer::NumDeriveVars();
        mf_plt[lev].setVal(0., cnt, num_spray_derive);
        SprayPC->computeDerivedVars(mf_plt[lev], lev, cnt);
        if (lev < finest_level) {
        MultiFab tmp_plt(
            grids[lev], dmap[lev], num_spray_derive, 0, MFInfo(), Factory(lev));
        tmp_plt.setVal(0.);
        VirtPC->computeDerivedVars(tmp_plt, lev, 0);
        MultiFab::Add(mf_plt[lev], tmp_plt, 0, cnt, num_spray_derive, 0);
        }
        cnt += num_spray_derive;
    }
    if (do_spray_particles && SprayParticleContainer::plot_spray_src) {
        SprayComps scomps = SprayParticleContainer::getSprayComps();
        MultiFab::Copy(
        mf_plt[lev], *m_spraysource[lev], scomps.rhoSrcIndx, cnt++, 1, 0);
        MultiFab::Copy(
        mf_plt[lev], *m_spraysource[lev], scomps.engSrcIndx, cnt++, 1, 0);
        MultiFab::Copy(
        mf_plt[lev], *m_spraysource[lev], scomps.momSrcIndx, cnt,
        AMREX_SPACEDIM, 0);
        cnt += AMREX_SPACEDIM;
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
            MultiFab::Copy(
            mf_plt[lev], *m_spraysource[lev], scomps.specSrcIndx + spf, cnt++, 1,
            0);
        }
    }
#endif
#ifdef PELE_USE_EFIELD
    if (m_do_extraEFdiags) {
        MultiFab::Copy(
        mf_plt[lev], *m_ionsFluxes[lev], 0, cnt, m_ionsFluxes[lev]->nComp(), 0);
        cnt += m_ionsFluxes[lev]->nComp();
    }
#endif
#if NUM_ODE > 0
    MultiFab::Copy(
        mf_plt[lev], m_leveldata_new[lev]->state, FIRSTODE, cnt, NUM_ODE, 0);
    cnt += NUM_ODE;
#endif

    if (m_do_les && m_plot_les) {
        constexpr amrex::Real fact = 0.5 / AMREX_SPACEDIM;
        auto const& plot_arr = mf_plt[lev].arrays();
        AMREX_D_TERM(auto const& mut_arr_x =
                        m_leveldata_old[lev]->visc_turb_fc[0].const_arrays();
                    , auto const& mut_arr_y =
                        m_leveldata_old[lev]->visc_turb_fc[1].const_arrays();
                    , auto const& mut_arr_z =
                        m_leveldata_old[lev]->visc_turb_fc[2].const_arrays();)
        // interpolate turbulent viscosity from faces to centers
        amrex::ParallelFor(
        mf_plt[lev],
        [plot_arr, AMREX_D_DECL(mut_arr_x, mut_arr_y, mut_arr_z),
            cnt] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
            plot_arr[box_no](i, j, k, cnt) =
            fact *
            (AMREX_D_TERM(
                mut_arr_x[box_no](i, j, k) + mut_arr_x[box_no](i + 1, j, k),
                +mut_arr_y[box_no](i, j, k) + mut_arr_y[box_no](i, j + 1, k),
                +mut_arr_z[box_no](i, j, k) + mut_arr_z[box_no](i, j, k + 1)));
        });
        Gpu::streamSynchronize();
    }

#ifdef AMREX_USE_EB
    if (m_plot_zeroEBcovered != 0) {
        EB_set_covered(mf_plt[lev], 0.0);
    }
#endif
    }

    //----------------------------------------------------------------
    // Blueprint : Components names
    Vector<std::string> names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
    names, &(eos_parms.host_parm()));

    Vector<std::string> plt_VarsName;
    AMREX_D_TERM(plt_VarsName.push_back("x_velocity");
                , plt_VarsName.push_back("y_velocity");
                , plt_VarsName.push_back("z_velocity"));
    if (m_incompressible == 0) {
        plt_VarsName.push_back("density");
        if (m_plotStateSpec != 0) {
            for (int n = 0; n < NUM_SPECIES; n++) {
            plt_VarsName.push_back("rho.Y(" + names[n] + ")");
            }
        }
    plt_VarsName.push_back("rhoh");
    plt_VarsName.push_back("temp");
    plt_VarsName.push_back("RhoRT");
#ifdef PELE_USE_EFIELD
    plt_VarsName.push_back("nE");
    plt_VarsName.push_back("phiV");
#endif
#ifdef PELE_USE_SOOT
    for (int mom = 0; mom < NUMSOOTVAR; mom++) {
        std::string sootname = soot_model->sootVariableName(mom);
        plt_VarsName.push_back(sootname);
    }
#endif
#ifdef PELE_USE_RADIATION
    if (do_rad_solve) {
        plt_VarsName.push_back("rad.G");
        plt_VarsName.push_back("rad.kappa");
        plt_VarsName.push_back("rad.emis");
    }
#endif
    if (m_has_divu != 0) {
        plt_VarsName.push_back("divu");
    }
    }

    if (m_plot_grad_p != 0) {
        AMREX_D_TERM(plt_VarsName.push_back("gradpx");
                    , plt_VarsName.push_back("gradpy");
                    , plt_VarsName.push_back("gradpz"));
    }

    if ((m_do_react != 0) && (m_skipInstantRR == 0) && (m_plot_react != 0)) {
        for (int n = 0; n < NUM_SPECIES; n++) {
            plt_VarsName.push_back("I_R(" + names[n] + ")");
    }
#ifdef PELE_USE_EFIELD
    plt_VarsName.push_back("I_R(nE)");
#endif
    plt_VarsName.push_back("FunctCall");
    // Extras:
    if (m_plotHeatRelease != 0) {
        plt_VarsName.push_back("HeatRelease");
    }
    }

#ifdef AMREX_USE_EB
    plt_VarsName.push_back("volFrac");
#endif

    for (int ivar = 0; ivar < m_derivePlotVarCount; ivar++) {
        const PeleLMDeriveRec* rec = derive_lst.get(m_derivePlotVars[ivar]);
        for (int dvar = 0; dvar < rec->numDerive(); dvar++) {
            plt_VarsName.push_back(rec->variableName(dvar));
        }
    }
#ifdef PELE_USE_SPRAY
    if (SprayParticleContainer::NumDeriveVars() > 0) {
    // We need virtual particles for the lower levels
        setupVirtualParticles(0);
        for (const auto& spray_derive_name :
                SprayParticleContainer::DeriveVarNames()) {
            plt_VarsName.push_back(spray_derive_name);
        }
    }
    if (do_spray_particles && SprayParticleContainer::plot_spray_src) {
        plt_VarsName.push_back("spray_mass_src");
        plt_VarsName.push_back("spray_energy_src");
        AMREX_D_TERM(plt_VarsName.push_back("spray_momentumX_src");
                        , plt_VarsName.push_back("spray_momentumY_src");
                        , plt_VarsName.push_back("spray_momentumZ_src"));
        for (const auto& spray_fuel_name :
                SprayParticleContainer::m_sprayDepNames) {
            plt_VarsName.push_back("spray_" + spray_fuel_name + "_src");
        }
    }
#endif

#ifdef PELE_USE_EFIELD
    if (m_do_extraEFdiags) {
        for (int ivar = 0; ivar < NUM_IONS; ++ivar) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            std::string dir = (idim == 0) ? "X" : ((idim == 1) ? "Y" : "Z");
            plt_VarsName.push_back(
                "DriftFlux_" + names[NUM_SPECIES - NUM_IONS + ivar] + "_" + dir);
            }
    }
    }
#endif

    if (m_do_les && m_plot_les) {
        plt_VarsName.push_back("viscturb");
    }

#if NUM_ODE > 0
    for (int n = 0; n < NUM_ODE; n++) {
        plt_VarsName.push_back(m_ode_names[n]);
    }
#endif

    //----------------------------------------------------------------
    // Blueprint : level 
    //Vector<int> level_steps(1,m_nstep);

    amrex::MultiLevelToBlueprint(
        finest_level + 1, amrex::GetVecOfConstPtrs(mf_plt), plt_VarsName, 
        Geom(), m_cur_time, level_steps, refRatio(), meshData);
    
    
    auto& root = node["catalyst/channels/mesh/data"];

    if(!root.dtype().is_object())
    {
        EmptyFieldData(plt_VarsName, root);
    }

    BL_PROFILE_VAR_STOP(copy_conduit_node);
    if(AMREX_SPACEDIM<3){
         AddDummyZAxis(root);
    }
    node.print();
    


    // Catalyst Execute
    std::string catalyst_execute_fun = "catalyst_execute_fun";
    BL_PROFILE_VAR("catalyst_execute_fun", catalyst_execute_fun);
    catalyst_status err = catalyst_execute(conduit::c_node(&node));
    BL_PROFILE_VAR_STOP(catalyst_execute_fun);

    if (err != catalyst_status_ok) {
        std::string message = " Error: Failed to execute Catalyst!\n";
        std::cerr << message << err << std::endl;
        amrex::Print() << message;
    }
    else{
        amrex::Print() << "Succesfully execute the catalyst \n";
    }
}

void PeleLM::compute_local_phi(std::vector<double> phi){
    double dx = PeleLM::prob_parm->slot_width / 2.0 / 64;
    for (int i=0; i < 64; i++){
        PeleLM::prob_parm->loc_phi[i] = phi[0] / (1+exp(phi[1]*(i*dx-phi[2]))) ;
    }
}


void PeleLM::CatalystSteering() {

    BL_PROFILE("PeleLM::CatalystSteering()");
    amrex::Print() << "running Catalyst Steering... \n"; 

    bool foundCell = false;
    const amrex::Real x_probe = (prob_parm->slot_width / 2.0) * 0.9;
    const amrex::Real y_probe = prob_parm->wall_height;
    amrex::Real localTemp = -9999.0;

    for (int lev = finest_level-1; lev >=finest_level-1 ; --lev)
        {

        const auto geomdata = geom[lev].data();
        const amrex::Real* dx_lev = geomdata.CellSize();
        const amrex::Real* plo_lev = geomdata.ProbLo();


        int iWanted_lev = static_cast<int>((x_probe - plo_lev[0]) / dx_lev[0]);
        int jWanted_lev = static_cast<int>((y_probe - plo_lev[1]) / dx_lev[1]);

        amrex::IntVect iv_lev(AMREX_D_DECL(iWanted_lev, jWanted_lev, 0));


        MultiFab& state_mf = m_leveldata_new[lev]->state;

        for (amrex::MFIter mfi(state_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            if (bx.contains(iv_lev))
            {
                auto const& arr = state_mf.const_array(mfi);
                localTemp = arr(iWanted_lev, jWanted_lev, 0, TEMP);

                foundCell = true;
                break;
            }
        }
        }

    Real tmax = -std::numeric_limits<Real>::max();
    for (int lev = 0; lev <= finest_level; ++lev) {
        const MultiFab& mf = m_leveldata_new[lev]->state;
        tmax = std::max(tmax, mf.max(TEMP, 0, true));  // local max
    }
    amrex::ParallelDescriptor::ReduceRealMax(tmax);

    

    // if(prev_steering_time == 0) {
    //     PrevTemp = localTemp;
    //     prev_steering_time = m_cur_time;
    //     return;
    // }
    // else {
    //     steering_dt = m_cur_time - prev_steering_time;
    // }
 
    // 2) Blueprint 
    conduit::Node node;
    auto & state = node["catalyst/state"];
    state["timestep"].set(m_nstep);
    state["time"].set(m_cur_time);

    // fill  node

    // fill steering node
    auto &steerable  = node["catalyst/channels/steerable"];
    steerable ["type"].set_string("mesh");
    auto &steerable_data = steerable["data"];
    steerable_data["coordsets/coords/type"].set_string("explicit");
    steerable_data["coordsets/coords/values/x"].set_float64_vector({ phiSteering[0] });
    steerable_data["coordsets/coords/values/y"].set_float64_vector({ phiSteering[1] });
    steerable_data["coordsets/coords/values/z"].set_float64_vector({ phiSteering[2] });

    steerable_data["topologies/mesh/type"].set("unstructured");
    steerable_data["topologies/mesh/coordset"].set("coords");
    steerable_data["topologies/mesh/elements/shape"].set("point");
    steerable_data["topologies/mesh/elements/connectivity"].set_int32_vector({ 0 });
    
    // steerable_data["fields/T_center/association"].set("vertex");
    // steerable_data["fields/T_center/topology"].set("mesh");
    // steerable_data["fields/T_center/volume_dependent"].set("false");
    // steerable_data["fields/T_center/values"].set_float64_vector(
    //     { PeleLM::prob_parm->T_center });

    // steerable_data["fields/V_mean/association"].set("vertex");
    // steerable_data["fields/V_mean/topology"].set("mesh");
    // steerable_data["fields/V_mean/values"].set_float64_vector(
    //     { PeleLM::prob_parm->V_mean });

    
    // steerable_data["fields/error/association"].set("vertex");
    // steerable_data["fields/error/topology"].set("mesh");
    // steerable_data["fields/error/volume_dependent"].set("false");
    // steerable_data["fields/error/values"].set_float64_vector({ error });

    steerable_data["fields/next_time/association"].set("vertex");
    steerable_data["fields/next_time/topology"].set("mesh");
    steerable_data["fields/next_time/volume_dependent"].set("false");
    steerable_data["fields/next_time/values"].set_float64_vector({ next_time });

    steerable_data["fields/temperature/association"].set("vertex");
    steerable_data["fields/temperature/topology"].set("mesh");
    steerable_data["fields/temperature/volume_dependent"].set("false");
    steerable_data["fields/temperature/values"].set_float64_vector({ localTemp });

    // steerable_data["fields/pv_temp_prev/association"].set("vertex");
    // steerable_data["fields/pv_temp_prev/topology"].set("mesh");
    // steerable_data["fields/pv_temp_prev/volume_dependent"].set("false");
    // steerable_data["fields/pv_temp_prev/values"].set_float64_vector({ PrevTemp });

    steerable_data["fields/steering_dt/association"].set("vertex");
    steerable_data["fields/steering_dt/topology"].set("mesh");
    steerable_data["fields/steering_dt/volume_dependent"].set("false");
    steerable_data["fields/steering_dt/values"].set_float64_vector({ steering_dt });

    steerable_data["fields/max_temp/association"].set("vertex");
    steerable_data["fields/max_temp/topology"].set("mesh");
    steerable_data["fields/max_temp/volume_dependent"].set("false");
    steerable_data["fields/max_temp/values"].set_float64_vector({ tmax });

    // steerable_data["fields/phi_global/association"].set("vertex");
    // steerable_data["fields/phi_global/topology"].set("mesh");
    // steerable_data["fields/phi_global/volume_dependent"].set("false");
    // steerable_data["fields/phi_global/values"].set_float64_vector({ PeleLM::prob_parm->phi });

    // steerable_data["fields/P00/association"].set("vertex");
    // steerable_data["fields/P00/topology"].set("mesh");
    // steerable_data["fields/P00/volume_dependent"].set("false");
    // steerable_data["fields/P00/values"].set_float64_vector({ p00 });

    // steerable_data["fields/P01/association"].set("vertex");
    // steerable_data["fields/P01/topology"].set("mesh");
    // steerable_data["fields/P01/volume_dependent"].set("false");
    // steerable_data["fields/P01/values"].set_float64_vector({ p01 });

    // steerable_data["fields/P10/association"].set("vertex");
    // steerable_data["fields/P10/topology"].set("mesh");
    // steerable_data["fields/P10/volume_dependent"].set("false");
    // steerable_data["fields/P10/values"].set_float64_vector({ p10 });

    // steerable_data["fields/P11/association"].set("vertex");
    // steerable_data["fields/P11/topology"].set("mesh");
    // steerable_data["fields/P11/volume_dependent"].set("false");
    // steerable_data["fields/P11/values"].set_float64_vector({ p11 });

    // steerable_data["fields/A/association"].set("vertex");
    // steerable_data["fields/A/topology"].set("mesh");
    // steerable_data["fields/A/volume_dependent"].set("false");
    // steerable_data["fields/A/values"].set_float64_vector({ A });

    // steerable_data["fields/B/association"].set("vertex");
    // steerable_data["fields/B/topology"].set("mesh");
    // steerable_data["fields/B/volume_dependent"].set("false");
    // steerable_data["fields/B/values"].set_float64_vector({ B });


    // for visualization && find maximum temperature
    /*
    auto& meshChannel = node["catalyst/channels/mesh"];
    meshChannel["type"].set_string("amrmesh");
    auto& meshData = meshChannel["data"];

    
    Vector<MultiFab*> mf_plt(finest_level + 1);
    Vector<int> level_steps(finest_level + 1);
    Vector<std::string> plt_VarsName;

    //plt_VarsName.push_back("temp");
    plt_VarsName.push_back("HR");

    // Blueprint : fill the multifabs
    for (int lev = 0; lev <= finest_level; ++lev) {

        level_steps[lev] = m_nstep;

        mf_plt[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), Factory(lev));
    
        //MultiFab::Copy(*mf_plt[lev], m_leveldata_new[lev]->state, TEMP, 0, 1, 0);

        if (m_plotHeatRelease != 0) {
            std::unique_ptr<MultiFab> mf;
            mf = std::make_unique<MultiFab>(grids[lev], dmap[lev], 1, 0);
            getHeatRelease(lev, mf.get());
            MultiFab::Copy(*mf_plt[lev], *mf, 0, 0, 1, 0);
        }
        //amrex::EB_set_covered(*mf_plt[lev], 0.0);
    }

    amrex::MultiLevelToBlueprint(
    finest_level + 1, amrex::GetVecOfConstPtrs(mf_plt), plt_VarsName, 
    Geom(), m_cur_time, level_steps, refRatio(), meshData);
    
    auto& root = node["catalyst/channels/mesh/data"];

    if (!root.dtype().is_object())
    {
        conduit::Node empty_node;
       //EmptyFieldData(plt_VarsName, root);
       catalyst_execute(conduit::c_node(&empty_node));
    }

    
    if(AMREX_SPACEDIM<3){
        AddDummyZAxis(root);
    }

    if(ParallelDescriptor::IOProcessor()){
        node.print();

    }

    //amrex::Abort("just stop here");
    */
    catalyst_status err = catalyst_execute(conduit::c_node(&node));
    if (err != catalyst_status_ok)
    {
        std::string message = " Error: Failed to execute Catalyst!\n";
        std::cerr << message << err << std::endl;
        amrex::Print() << message;
        amrex::Abort(message);
    } 
    else{
        amrex::Print() << "Succesfully execute the catalyst \n";
    }
    
    //for (int lev = 0; lev <= finest_level; ++lev) {
    //    delete mf_plt[lev];
    //}

    conduit::Node result;
    catalyst_status err_result = catalyst_results(conduit::c_node(&result));
    if (err_result != catalyst_status_ok)
    {
        std::cerr << "Failed to execute Catalyst-results: " << err_result << std::endl;
    }
    else
    {
        int localFlag   = (foundCell ? ParallelDescriptor::MyProc() : -1);
        ParallelDescriptor::ReduceIntMax(localFlag);
        int rootRank    = localFlag; 
        
        double localSteering_dt = steering_dt;
        // double localError           = error;
        // double localIntegral        = integral;
        //double localVmean           = PeleLM::prob_parm->V_mean;    
        double localnext_time = next_time;
        std::array<double,3> localPhi;
        for(int i = 0 ; i <3; i ++) {
            localPhi[i] = phiSteering[i];
        }
        //double localTcenter         = PeleLM::prob_parm->T_center;
        // int    localSteeringInt     = inSitu_Steering_int;
        // double localmaxTemp         = max_temp;
        // double localPrevTemp      = PrevTemp;
        // double localP00           = p00;
        // double localP01           = p01;
        // double localP10           = p10;
        // double localP11           = p11;
        // double localA             = A;
        // double localB             = B;

        if (rootRank != -1) {

            if (ParallelDescriptor::MyProc() == rootRank)
            {   
                auto &fields  = result["catalyst/steerable/fields"];
                //localError    = fields["error/values"].to_double();
                //localIntegral = fields["integral/values"].to_double();
                localPhi[0]   = fields["phi_local/values/0"].to_double();
                localPhi[1]   = fields["phi_local/values/1"].to_double();
                localPhi[2]   = fields["phi_local/values/2"].to_double();
                //localSteering_dt = fields["steering_dt/values"].to_double();
                localnext_time = fields["next_time/values"].to_double();
                if (localnext_time == 0){
                    localnext_time = next_time;
                }
                if (localPhi[0]==0 && localPhi[1]==0 && localPhi[2]==0){
                    localPhi[0] = phiSteering[0];
                    localPhi[1] = phiSteering[1];
                    localPhi[2] = phiSteering[2];
                }
                // localSteeringInt = m_nstep + static_cast<int>(
                //                     (y_probe/PeleLM::prob_parm->V_mean)/m_dt*5
                //                 );
                // localSteeringInt = 1;
                // localPrevTemp = fields["pv_temp_prev/values"].to_double();
                // localmaxTemp  = fields["max_temp/values"].to_double();
                // localP00      = fields["P00/values"].to_double();
                // localP01      = fields["P01/values"].to_double();
                // localP10      = fields["P10/values"].to_double();
                // localP11      = fields["P11/values"].to_double();
                // localA        = fields["A/values"].to_double();
                // localB        = fields["B/values"].to_double();

            }

            // ParallelDescriptor::Bcast(&localError,      1, rootRank);
            // ParallelDescriptor::Bcast(&localIntegral,   1, rootRank);
            //ParallelDescriptor::Bcast(&localVmean,      1, rootRank);
            ParallelDescriptor::Bcast(localPhi.data(),  3, rootRank);
            //ParallelDescriptor::Bcast(&localSteering_dt,1, rootRank);
            ParallelDescriptor::Bcast(&localnext_time,  1, rootRank);
            //ParallelDescriptor::Bcast(&localTcenter,    1, rootRank);
            //ParallelDescriptor::Bcast(&localSteeringInt,1, rootRank);
            // ParallelDescriptor::Bcast(&localmaxTemp,1, rootRank);
            // ParallelDescriptor::Bcast(&localPrevTemp,1, rootRank);
            // ParallelDescriptor::Bcast(&localP00,1, rootRank);
            // ParallelDescriptor::Bcast(&localP01,1, rootRank);
            // ParallelDescriptor::Bcast(&localP10,1, rootRank);
            // ParallelDescriptor::Bcast(&localP11,1, rootRank);
            // ParallelDescriptor::Bcast(&localA,1, rootRank);
            // ParallelDescriptor::Bcast(&localB,1, rootRank);


            // error                          = localError;
            // integral                       = localIntegral;
            // PrevTemp                       = localPrevTemp;
            //PeleLM::prob_parm->V_mean      = localVmean;
            //PeleLM::prob_parm->T_center    = localTcenter;
            for(int i = 0 ; i <3; i ++) {
                phiSteering[i] = localPhi[i];
            }
            //steering_dt = localSteering_dt;
            next_time = localnext_time;
            // inSitu_Steering_int            = localSteeringInt;
            // p00                            = localP00;
            // p01                            = localP01;
            // p10                            = localP10;
            // p11                            = localP11;
            // A                               = localA;
            // B                               = localB;


            compute_local_phi(phiSteering);
            /*
            inSitu_Steering_int = m_nstep + static_cast<int>(
                                 (y_probe/PeleLM::prob_parm->V_mean)*5);

            amrex::Print() << "[next steering time " << m_cur_time +  (y_probe/PeleLM::prob_parm->V_mean)*5 << "approximately at steps  "<< static_cast<int>(
                (y_probe/PeleLM::prob_parm->V_mean)*5/m_dt) << "\n";
            */
            //inSitu_Steering_int = 1;
            
            // amrex::Print() << 
            //            "[steering] next steering time :"
            //            << next_time
            //            << "\n";
            
            
            Gpu::copy(Gpu::hostToDevice, prob_parm, prob_parm + 1, prob_parm_d);
        }
        
    }
    // prev_steering_time = m_cur_time;
}
    


void PeleLM::CatalystFinalize() {
    conduit::Node node;
    catalyst_status err = catalyst_finalize(conduit::c_node(&node));
    if (err != catalyst_status_ok)
    {
        std::string message = " Error: Failed to finalize Catalyst!\n";
        std::cerr << message << err << std::endl;
        amrex::Print() << message;
        amrex::Abort(message);
    } else {
        amrex::Print() << "Successfully finalized Catalyst\n";
    }
}
#endif
#endif