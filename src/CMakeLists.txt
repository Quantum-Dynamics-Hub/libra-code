#
#  Compile PCH first
#
#file(SET PCH_HEADERS pch.h)
#file(SET PCH_SRC pch.cpp pch.h) 
set(PCH_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/pch.h)
set(PCH_SRC ${CMAKE_CURRENT_SOURCE_DIR}/pch.cpp ${PCH_HEADERS})

ADD_LIBRARY(pch SHARED ${PCH_SRC})
ADD_LIBRARY(pch_stat STATIC ${PCH_SRC})

target_precompile_headers(pch PUBLIC ${PCH_SRC}) # With PUBLIC they will be used by targets using this target
target_precompile_headers(pch_stat PUBLIC ${PCH_SRC}) # With PUBLIC they will be used by targets using this target


#add_compile_definitions(USING_PCH)

TARGET_LINK_LIBRARIES(pch      ${ext_libs} )
TARGET_LINK_LIBRARIES(pch_stat ${ext_libs} )


#
#  Building sub-directories
#
MESSAGE("Going into subdirectory util...")
ADD_SUBDIRECTORY("util")

MESSAGE("Going into subdirectory io...")
ADD_SUBDIRECTORY("io")

MESSAGE("Going into subdirectory timer...")
ADD_SUBDIRECTORY("timer")

MESSAGE("Going into subdirectory common_types...")
ADD_SUBDIRECTORY("common_types")




MESSAGE("Going into subdirectory linalg...")
ADD_SUBDIRECTORY("math_linalg")

MESSAGE("Going into subdirectory meigen...")
ADD_SUBDIRECTORY("math_meigen")

MESSAGE("Going into subdirectory specialfunctions...")
ADD_SUBDIRECTORY("math_specialfunctions")

MESSAGE("Going into subdirectory symmetry...")
ADD_SUBDIRECTORY("math_symmetry")

MESSAGE("Going into subdirectory operators...")
ADD_SUBDIRECTORY("math_operators")

MESSAGE("Going into subdirectory random...")
ADD_SUBDIRECTORY("math_random")

MESSAGE("Going into subdirectory graph...")
ADD_SUBDIRECTORY("math_graph")

MESSAGE("Going into subdirectory data...")
ADD_SUBDIRECTORY("math_data")

MESSAGE("Going into subdirectory ann...")
ADD_SUBDIRECTORY("math_ann")



MESSAGE("Going into subdirectory molint...")
ADD_SUBDIRECTORY("molint")

MESSAGE("Going into subdirectory qobjects...")
ADD_SUBDIRECTORY("qobjects")

MESSAGE("Going into subdirectory basis...")
ADD_SUBDIRECTORY("basis")

MESSAGE("Going into subdirectory libint2_wrappers...")
ADD_SUBDIRECTORY("libint2_wrappers")



MESSAGE("Going into subdirectory calculators...")
ADD_SUBDIRECTORY("calculators")


MESSAGE("Going into subdirectory models...")
ADD_SUBDIRECTORY("models")


MESSAGE("Going into subdirectory context...")
ADD_SUBDIRECTORY("context")



MESSAGE("Going into subdirectory rigidbody...")
ADD_SUBDIRECTORY("dyn_rigidbody")



MESSAGE("Going into subdirectory chemobjects...")
ADD_SUBDIRECTORY("chemobjects")


MESSAGE("Going into subdirectory cell...")
ADD_SUBDIRECTORY("cell")


MESSAGE("Going into subdirectory opt...")
ADD_SUBDIRECTORY("opt")


MESSAGE("Going into subdirectory pot...")
ADD_SUBDIRECTORY("pot")


MESSAGE("Going into subdirectory forcefield...")
ADD_SUBDIRECTORY("forcefield")



MESSAGE("Going into subdirectory control_parameters...")
ADD_SUBDIRECTORY("control_parameters")

MESSAGE("Going into subdirectory model_parameters...")
ADD_SUBDIRECTORY("model_parameters")


MESSAGE("Going into subdirectory basis_setups...")
ADD_SUBDIRECTORY("basis_setups")



MESSAGE("Going into subdirectory atomistic...")
ADD_SUBDIRECTORY("atomistic")



MESSAGE("Going into subdirectory nhamiltonian...")
ADD_SUBDIRECTORY("nhamiltonian")


MESSAGE("Going into subdirectory montecarlo...")
ADD_SUBDIRECTORY("montecarlo")


MESSAGE("Going into subdirectory fgr...")
ADD_SUBDIRECTORY("fgr")


MESSAGE("Going into subdirectory ivr...")
ADD_SUBDIRECTORY("ivr")


MESSAGE("Going into subdirectory dyn...")
ADD_SUBDIRECTORY("dyn")


MESSAGE("Going into subdirectory converters...")
ADD_SUBDIRECTORY("converters")


MESSAGE("Going into subdirectory scripts...")
ADD_SUBDIRECTORY("scripts")

MESSAGE("Going into subdirectory qchem_tools...")
ADD_SUBDIRECTORY("qchem_tools")


MESSAGE("Going into subdirectory solvers...")
ADD_SUBDIRECTORY("solvers")


MESSAGE("Going into subdirectory integrators...")
ADD_SUBDIRECTORY("integrators")



#
#  Source files and headers in this directory
#
#file(GLOB LIBRA_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
#file(GLOB LIBRA_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${LIBRA_HEADERS}) 

#file(GLOB LIBRA_HEADERS liblibra_core.h realtype.h Units.h)
#file(GLOB LIBRA_SRC liblibra_core.cpp ${LIBRA_HEADERS}) 
set(LIBRA_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/liblibra_core.h ${CMAKE_CURRENT_SOURCE_DIR}/realtype.h ${CMAKE_CURRENT_SOURCE_DIR}/Units.h)
set(LIBRA_SRC ${CMAKE_CURRENT_SOURCE_DIR}/liblibra_core.cpp ${LIBRA_HEADERS})


ADD_LIBRARY(libra_core_stat STATIC ${LIBRA_SRC})
ADD_LIBRARY(libra_core SHARED ${LIBRA_SRC})

#ADD_LIBRARY(libra_core SHARED ${LIBRA_SRC})


#
#  PCH (Precompiled header) support
# 
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support")
   #target_precompile_headers(libra_core PUBLIC pch.h) # With PUBLIC they will be used by targets using this target
   #target_precompile_headers(libra_core_stat PUBLIC pch.h) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(libra_core REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(libra_core_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
   set(CMAKE_PCH_INSTANTIATE_TEMPLATES ON)
   add_compile_definitions(USING_PCH)
else()
   message(STATUS "No pre-compiled header support")
endif()





TARGET_LINK_LIBRARIES(libra_core_stat
                      integrators_stat
                      solvers_stat
                      qchem_tools_stat
                      atomistic_stat
                      scripts_stat
                      converters_stat
                      montecarlo_stat
                      fgr_stat
                      ivr_stat
                      dyn_stat
                      common_types_stat
                      basis_setups_stat
                      model_parameters_stat
                      models_stat
                      control_parameters_stat
                      nhamiltonian_stat
                      util_stat
                      forcefield_stat
                      opt_stat
                      pot_stat
                      cell_stat
                      chemobjects_stat
                      rigidbody_stat
                      basis_stat
                      qobjects_stat
                      libint2_wrappers_stat
                      molint_stat
                      linalg_stat meigen_stat specialfunctions_stat symmetry_stat operators_stat random_stat graph_stat data_stat ann_stat timer_stat 
                      context_stat io_stat calculators_stat ${ext_libs}  )


TARGET_LINK_LIBRARIES(libra_core libra_core_stat)



#
#  Copy python files
#
add_custom_target(copy-libra-py ALL
    COMMAND cmake -E remove_directory ${CMAKE_BINARY_DIR}/src/libra_py
    COMMAND cmake -E copy_directory ${CMAKE_SOURCE_DIR}/src/libra_py ${CMAKE_BINARY_DIR}/src/libra_py
    COMMAND cmake -E copy ${CMAKE_SOURCE_DIR}/src/util/__init__.py ${CMAKE_BINARY_DIR}/src/util/__init__.py
    COMMAND cmake -E copy ${CMAKE_SOURCE_DIR}/src/__init__.py ${CMAKE_BINARY_DIR}/src/__init__.py
)

add_custom_target(copy-libra-params ALL
    COMMAND cmake -E copy_directory ${CMAKE_SOURCE_DIR}/src/libra_params ${CMAKE_BINARY_DIR}/src/libra_params
)




