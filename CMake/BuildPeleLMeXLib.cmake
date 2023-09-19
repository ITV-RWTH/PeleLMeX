function(build_pelelmex_lib pelelmex_lib_name)
  if (NOT (TARGET ${pelelmex_lib_name}))
    add_library(${pelelmex_lib_name} OBJECT)

    set(PELE_PHYSICS_SRC_DIR ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics)
    set(PELE_PHYSICS_TRANSPORT_DIR "${PELE_PHYSICS_SRC_DIR}/Transport")
    set(PELE_PHYSICS_EOS_DIR "${PELE_PHYSICS_SRC_DIR}/Eos")
    set(PELE_PHYSICS_MECHANISM_DIR "${PELE_PHYSICS_SRC_DIR}/Support/Mechanism/Models/${PELELMEX_CHEMISTRY_MODEL}")
    set(AMREX_SUNDIALS_DIR ${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS)

    if(CLANG_TIDY_EXE)
      set_target_properties(${pelelmex_lib_name} PROPERTIES CXX_CLANG_TIDY
                            "${CLANG_TIDY_EXE};--config-file=${CMAKE_SOURCE_DIR}/.clang-tidy")
    endif()

    include(SetPeleLMeXCompileFlags)
    
    target_sources(${pelelmex_lib_name}
      PRIVATE
        ${PELE_PHYSICS_SRC_DIR}/Utility/TurbInflow/turbinflow.cpp
        ${PELE_PHYSICS_SRC_DIR}/Utility/TurbInflow/turbinflow.H
    )
    target_include_directories(${pelelmex_lib_name} PUBLIC ${PELE_PHYSICS_SRC_DIR}/Utility/TurbInflow)

    target_sources(${pelelmex_lib_name}
      PRIVATE
        ${PELE_PHYSICS_SRC_DIR}/Utility/Diagnostics/DiagBase.H
        ${PELE_PHYSICS_SRC_DIR}/Utility/Diagnostics/DiagBase.cpp
        ${PELE_PHYSICS_SRC_DIR}/Utility/Diagnostics/DiagConditional.H
        ${PELE_PHYSICS_SRC_DIR}/Utility/Diagnostics/DiagConditional.cpp
        ${PELE_PHYSICS_SRC_DIR}/Utility/Diagnostics/DiagFilter.H
        ${PELE_PHYSICS_SRC_DIR}/Utility/Diagnostics/DiagFilter.cpp
        ${PELE_PHYSICS_SRC_DIR}/Utility/Diagnostics/DiagFramePlane.H
        ${PELE_PHYSICS_SRC_DIR}/Utility/Diagnostics/DiagFramePlane.cpp
        ${PELE_PHYSICS_SRC_DIR}/Utility/Diagnostics/DiagPDF.H
        ${PELE_PHYSICS_SRC_DIR}/Utility/Diagnostics/DiagPDF.cpp
    )
    target_include_directories(${pelelmex_lib_name} PUBLIC ${PELE_PHYSICS_SRC_DIR}/Utility/Diagnostics)
    
    target_sources(${pelelmex_lib_name}
      PRIVATE
        ${PELE_PHYSICS_SRC_DIR}/Utility/PltFileManager/PltFileManager.cpp
        ${PELE_PHYSICS_SRC_DIR}/Utility/PltFileManager/PltFileManager.H
        ${PELE_PHYSICS_SRC_DIR}/Utility/PltFileManager/PltFileManagerBCFill.H
    )
    target_include_directories(${pelelmex_lib_name} PUBLIC ${PELE_PHYSICS_SRC_DIR}/Utility/PltFileManager)

    target_sources(${pelelmex_lib_name}
      PRIVATE
        ${PELE_PHYSICS_SRC_DIR}/Utility/PMF/PMF.H
        ${PELE_PHYSICS_SRC_DIR}/Utility/PMF/PMFData.cpp
        ${PELE_PHYSICS_SRC_DIR}/Utility/PMF/PMFData.H
    )
    target_include_directories(${pelelmex_lib_name} PUBLIC ${PELE_PHYSICS_SRC_DIR}/Utility/PMF)
    
    target_sources(${pelelmex_lib_name}
      PRIVATE 
        ${AMREX_SUNDIALS_DIR}/AMReX_Sundials.H
        ${AMREX_SUNDIALS_DIR}/AMReX_Sundials_Core.cpp
        ${AMREX_SUNDIALS_DIR}/AMReX_Sundials_Core.H
        ${AMREX_SUNDIALS_DIR}/AMReX_NVector_MultiFab.cpp
        ${AMREX_SUNDIALS_DIR}/AMReX_NVector_MultiFab.H
        ${AMREX_SUNDIALS_DIR}/AMReX_SUNMemory.cpp
        ${AMREX_SUNDIALS_DIR}/AMReX_SUNMemory.H
    )
    target_include_directories(${pelelmex_lib_name} SYSTEM PUBLIC ${AMREX_SUNDIALS_DIR})

    target_include_directories(${pelelmex_lib_name} PUBLIC "${PELE_PHYSICS_SRC_DIR}/Source")

    target_sources(${pelelmex_lib_name}
      PRIVATE
        ${PELE_PHYSICS_TRANSPORT_DIR}/Transport.H
        ${PELE_PHYSICS_TRANSPORT_DIR}/Transport.cpp
        ${PELE_PHYSICS_TRANSPORT_DIR}/TransportParams.H
        ${PELE_PHYSICS_TRANSPORT_DIR}/TransportTypes.H
        ${PELE_PHYSICS_TRANSPORT_DIR}/Constant.H
        ${PELE_PHYSICS_TRANSPORT_DIR}/Simple.H
        ${PELE_PHYSICS_TRANSPORT_DIR}/Sutherland.H
    )
    target_include_directories(${pelelmex_lib_name} PUBLIC ${PELE_PHYSICS_TRANSPORT_DIR})
    if("${PELELMEX_TRANSPORT_MODEL}" STREQUAL "Constant")
      target_compile_definitions(${pelelmex_lib_name} PUBLIC USE_CONSTANT_TRANSPORT)
    endif()
    if("${PELELMEX_TRANSPORT_MODEL}" STREQUAL "Simple")
      target_compile_definitions(${pelelmex_lib_name} PUBLIC USE_SIMPLE_TRANSPORT)
    endif()
    if("${PELELMEX_TRANSPORT_MODEL}" STREQUAL "Sutherland")
      target_compile_definitions(${pelelmex_lib_name} PUBLIC USE_SUTHERLAND_TRANSPORT)
    endif()

    target_sources(${pelelmex_lib_name}
      PRIVATE
        ${PELE_PHYSICS_EOS_DIR}/EOS.cpp
        ${PELE_PHYSICS_EOS_DIR}/EOS.H
        ${PELE_PHYSICS_EOS_DIR}/GammaLaw.H
        ${PELE_PHYSICS_EOS_DIR}/Fuego.H
        ${PELE_PHYSICS_EOS_DIR}/SRK.H
    )
    target_include_directories(${pelelmex_lib_name} PUBLIC ${PELE_PHYSICS_EOS_DIR})
    if("${PELELMEX_EOS_MODEL}" STREQUAL "GammaLaw")
      target_compile_definitions(${pelelmex_lib_name} PUBLIC USE_GAMMALAW_EOS)
    endif()
    if("${PELELMEX_EOS_MODEL}" STREQUAL "Fuego")
      target_compile_definitions(${pelelmex_lib_name} PUBLIC USE_FUEGO_EOS)
    endif()
    if("${PELELMEX_EOS_MODEL}" STREQUAL "Soave-Redlich-Kwong")
      target_compile_definitions(${pelelmex_lib_name} PUBLIC USE_SRK_EOS)
    endif()

    target_sources(${pelelmex_lib_name}
      PRIVATE
        ${PELE_PHYSICS_MECHANISM_DIR}/mechanism.cpp
        ${PELE_PHYSICS_MECHANISM_DIR}/mechanism.H
    )
    target_include_directories(${pelelmex_lib_name} SYSTEM PUBLIC ${PELE_PHYSICS_MECHANISM_DIR})

    target_sources(${pelelmex_lib_name}
      PRIVATE
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorArkode.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorArkode.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorBase.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorBase.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvode.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvode.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeCustomLinSolver.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeCustomLinSolver.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeJacobian.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeJacobian.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodePreconditioner.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodePreconditioner.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeUtils.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeUtils.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorNull.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorNull.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorRK64.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorRK64.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorTypes.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorUtils.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorUtils.cpp
    )
    target_include_directories(${pelelmex_lib_name} PUBLIC ${PELE_PHYSICS_SRC_DIR}/Reactions)

    include(AMReXBuildInfo)
    generate_buildinfo(${pelelmex_lib_name} ${CMAKE_SOURCE_DIR})
    target_include_directories(${pelelmex_lib_name} SYSTEM PUBLIC ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)
    
    target_link_libraries(${pelelmex_lib_name} PUBLIC sundials_arkode sundials_cvode)
    
    if(PELELMEX_ENABLE_CUDA)
      target_link_libraries(${pelelmex_lib_name} PUBLIC sundials_nveccuda sundials_sunlinsolcusolversp sundials_sunmatrixcusparse)
    elseif(PELELMEX_ENABLE_HIP)
      target_link_libraries(${pelelmex_lib_name} PUBLIC sundials_nvechip)
    elseif(PELELMEX_ENABLE_SYCL)
      target_link_libraries(${pelelmex_lib_name} PUBLIC sundials_nvecsycl)
    endif()
    
    if(PELELMEX_ENABLE_MPI)
      target_link_libraries(${pelelmex_lib_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
    endif()
    
    #Link to amrex libraries
    target_link_libraries(${pelelmex_lib_name} PUBLIC AMReX-Hydro::amrex_hydro_api AMReX::amrex)

  endif()
endfunction()