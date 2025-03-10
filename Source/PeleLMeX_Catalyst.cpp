#ifndef PeleLMeX_Catalyst_H
#define PeleLMeX_Catalyst_H

#include <PeleLMeX.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>

#include <string>

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
    error = 0;
    integral = 0;
    

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

        // Look for a uniform coordset at dom["coordsets/coords"]
        if (!dom.has_path("coordsets/coords/type")) continue;
        if (std::string(dom["coordsets/coords/type"].as_string()) != "uniform")
            continue;

        // If this domain is already 3D, skip. We check if dims/k exists:
        conduit::Node &coords = dom["coordsets/coords"];
        if (coords.has_path("dims/k"))
        {
            // This domain is already 3D, no need to modify
            continue;
        }

        // At this point, we assume it's 2D and we want to add the dummy third dimension
        // 1) Add dims/k = 2
        coords["dims/k"] = (conduit::int32)2;

        // 2) Add spacing/dz
        coords["spacing/dz"] = 1.0;  // or your chosen dummy value

        // 3) Add origin/z
        coords["origin/z"] = 0.0;

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


void PeleLM::CatalystExecute () {
    amrex::Print() << "Running Catalyst pipeline scripts... \n";
    BL_PROFILE("PeleLM::CatalystExecute()");
    //----------------------------------------------------------------
    // Blueprint : Mesh data
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
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf_plt[lev].define(grids[lev], dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
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
    Vector<int> level_steps(1,m_nstep);

    amrex::MultiLevelToBlueprint(
        finest_level + 1, amrex::GetVecOfConstPtrs(mf_plt), plt_VarsName, 
        Geom(), m_cur_time, level_steps, refRatio(), meshData);
    
    if (AMREX_SPACEDIM == 2) {
        //AddDummyZAxis(meshData);
    }
    // Catalyst Execute
    catalyst_status err = catalyst_execute(conduit::c_node(&node));
    if (err != catalyst_status_ok) {
        std::string message = " Error: Failed to execute Catalyst!\n";
        std::cerr << message << err << std::endl;
        amrex::Print() << message;
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
        if (foundCell) {
            std::cout << "we found the cell" << std::endl;
            break;
        } 
    }
 
    // 2) Blueprint 
    conduit::Node node;
    auto & state = node["catalyst/state"];
    state["timestep"].set(m_nstep);
    state["time"].set(m_cur_time);

    // fill steering node
    auto &steerable  = node["catalyst/channels/steerable"];
    steerable ["type"].set_string("mesh");
    auto &steerable_data = steerable["data"];
    steerable_data["coordsets/coords/type"].set_string("explicit");
    steerable_data["coordsets/coords/values/x"].set_float64_vector({ 1 });
    steerable_data["coordsets/coords/values/y"].set_float64_vector({ 2 });
    steerable_data["coordsets/coords/values/z"].set_float64_vector({ 3 });

    steerable_data["topologies/mesh/type"].set("unstructured");
    steerable_data["topologies/mesh/coordset"].set("coords");
    steerable_data["topologies/mesh/elements/shape"].set("point");
    steerable_data["topologies/mesh/elements/connectivity"].set_int32_vector({ 0 });
    
    steerable_data["fields/T_center/association"].set("vertex");
    steerable_data["fields/T_center/topology"].set("mesh");
    steerable_data["fields/T_center/volume_dependent"].set("false");
    steerable_data["fields/T_center/values"].set_float64_vector(
        { PeleLM::prob_parm->T_center });

    steerable_data["fields/V_mean/association"].set("vertex");
    steerable_data["fields/V_mean/topology"].set("mesh");
    steerable_data["fields/V_mean/volume_dependent"].set("false");
    steerable_data["fields/V_mean/values"].set_float64_vector(
        { PeleLM::prob_parm->V_mean });

    steerable_data["fields/phi/association"].set("vertex");
    steerable_data["fields/phi/topology"].set("mesh");
    steerable_data["fields/phi/volume_dependent"].set("false");
    steerable_data["fields/phi/values"].set_float64_vector({ PeleLM::prob_parm->phi });

    steerable_data["fields/error/association"].set("vertex");
    steerable_data["fields/error/topology"].set("mesh");
    steerable_data["fields/error/volume_dependent"].set("false");
    steerable_data["fields/error/values"].set_float64_vector({ error });

    steerable_data["fields/integral/association"].set("vertex");
    steerable_data["fields/integral/topology"].set("mesh");
    steerable_data["fields/integral/volume_dependent"].set("false");
    steerable_data["fields/integral/values"].set_float64_vector({ integral });

    steerable_data["fields/temperature/association"].set("vertex");
    steerable_data["fields/temperature/topology"].set("mesh");
    steerable_data["fields/temperature/volume_dependent"].set("false");
    steerable_data["fields/temperature/values"].set_float64_vector({ localTemp });

    catalyst_status err = catalyst_execute(conduit::c_node(&node));
    if (err != catalyst_status_ok)
    {
        std::string message = " Error: Failed to execute Catalyst!\n";
        std::cerr << message << err << std::endl;
        amrex::Print() << message;
        amrex::Abort(message);
    } 
}

    


void PeleLM::CatalystResult() {
    conduit::Node node;
    //amrex::Print() << "catalyst result...\n";

    catalyst_status err = catalyst_results(conduit::c_node(&node));

    if (err != catalyst_status_ok)
    {
        std::cerr << "Failed to execute Catalyst-results: " << err << std::endl;
    }
    else
    {
        auto &timestep_node = node["catalyst/state/timestep"];
        int timestep = timestep_node.to_int();
        if(timestep != 0){
            auto &results = node["catalyst/steerable/fields"];
            error = results["error/values"].to_double();
            integral = results["integral/values"].to_double();
            PeleLM::prob_parm->V_mean = results["V_mean/values"].to_double();
            //inSitu_Steering_int = m_nstep + static_cast<int>((PeleLM::prob_parm->wall_height/PeleLM::prob_parm->V_mean)/m_dt);
            std::cout << "Next Time Step to steering is :" << inSitu_Steering_int ;
        }
    }
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
        std::cout << "Successfully finalized Catalyst\n";
    }
}
#endif
#endif