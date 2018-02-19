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




MESSAGE("Going into subdirectory calculators...")
ADD_SUBDIRECTORY("calculators")


MESSAGE("Going into subdirectory context...")
ADD_SUBDIRECTORY("context")



MESSAGE("Going into subdirectory rigidbody...")
ADD_SUBDIRECTORY("dyn_rigidbody")



MESSAGE("Going into subdirectory chemobjects...")
ADD_SUBDIRECTORY("chemobjects")


MESSAGE("Going into subdirectory cell...")
ADD_SUBDIRECTORY("cell")

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



MESSAGE("Going into subdirectory hamiltonian...")
ADD_SUBDIRECTORY("hamiltonian")


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



#
#  Source files and headers in this directory
#
file(GLOB LIBRA_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB LIBRA_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${LIBRA_HEADERS}) 

ADD_LIBRARY(libra_core SHARED ${LIBRA_SRC})
ADD_LIBRARY(libra_core_stat STATIC ${LIBRA_SRC})



TARGET_LINK_LIBRARIES(libra_core
                      solvers_stat
                      qchem_tools_stat
                      scripts_stat
                      converters_stat
                      ivr_stat
                      dyn_stat
                      common_types_stat                    
                      basis_setups_stat
                      model_parameters_stat
                      control_parameters_stat
                      hamiltonian_stat
                      util_stat
                      forcefield_stat
                      pot_stat
                      cell_stat
                      chemobjects_stat
                      rigidbody_stat
                      basis_stat
                      qobjects_stat
                      molint_stat
                      linalg_stat meigen_stat specialfunctions_stat symmetry_stat operators_stat random_stat graph_stat data_stat ann_stat timer_stat 
                      context_stat io_stat calculators_stat ${ext_libs}  )


TARGET_LINK_LIBRARIES(libra_core_stat
                      solvers_stat
                      qchem_tools_stat
                      scripts_stat
                      converters_stat
                      ivr_stat
                      dyn_stat
                      common_types_stat
                      basis_setups_stat
                      model_parameters_stat
                      control_parameters_stat
                      hamiltonian_stat
                      util_stat
                      forcefield_stat
                      pot_stat
                      cell_stat
                      chemobjects_stat
                      rigidbody_stat
                      basis_stat
                      qobjects_stat
                      molint_stat
                      linalg_stat meigen_stat specialfunctions_stat symmetry_stat operators_stat random_stat graph_stat data_stat ann_stat timer_stat 
                      context_stat io_stat calculators_stat ${ext_libs}  )

#
#  Copy python files
#
add_custom_target(copy-libra-py ALL
    COMMAND cmake -E copy_directory ${CMAKE_SOURCE_DIR}/src/libra_py ${CMAKE_BINARY_DIR}/src/libra_py
)

add_custom_target(copy-libra-params ALL
    COMMAND cmake -E copy_directory ${CMAKE_SOURCE_DIR}/src/libra_params ${CMAKE_BINARY_DIR}/src/libra_params
)




