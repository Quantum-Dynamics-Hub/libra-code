#
# Set the libraries
# 
# For Cygwin
SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}/libpython2.6.dll )   
# For Linux
#SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ) 



#
#  Building sub-directories
#
MESSAGE("Going into subdirectory calculators...")
ADD_SUBDIRECTORY("calculators")


MESSAGE("Going into subdirectory io...")
ADD_SUBDIRECTORY("io")


MESSAGE("Going into subdirectory context...")
ADD_SUBDIRECTORY("context")


MESSAGE("Going into subdirectory mmath...")
ADD_SUBDIRECTORY("mmath")


MESSAGE("Going into subdirectory chemobjects...")
ADD_SUBDIRECTORY("chemobjects")


MESSAGE("Going into subdirectory qchem...")
ADD_SUBDIRECTORY("qchem")


MESSAGE("Going into subdirectory cell...")
ADD_SUBDIRECTORY("cell")

MESSAGE("Going into subdirectory pot...")
ADD_SUBDIRECTORY("pot")


MESSAGE("Going into subdirectory hamiltonian...")
ADD_SUBDIRECTORY("hamiltonian")


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

MESSAGE("Going into subdirectory util...")
ADD_SUBDIRECTORY("util")


#
#  Source files and headers in this directory
#
file(GLOB LIBRA_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB LIBRA_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${LIBRA_HEADERS}) 

ADD_LIBRARY(libra_core SHARED ${LIBRA_SRC})
ADD_LIBRARY(libra_core_stat STATIC ${LIBRA_SRC})

TARGET_LINK_LIBRARIES(libra_core      util_stat solvers_stat qchem_tools_stat scripts_stat converters_stat dyn_stat hamiltonian_stat pot_stat cell_stat qchem_stat chemobjects_stat mmath_stat context_stat io_stat calculators_stat ${ext_libs}  )
TARGET_LINK_LIBRARIES(libra_core_stat util_stat solvers_stat qchem_tools_stat scripts_stat converters_stat dyn_stat hamiltonian_stat pot_stat cell_stat qchem_stat chemobjects_stat mmath_stat context_stat io_stat calculators_stat ${ext_libs}  )


#
#  Copy python files
#
add_custom_target(copy-libra-py ALL
    COMMAND cmake -E copy_directory ${CMAKE_SOURCE_DIR}/src/libra_py ${CMAKE_BINARY_DIR}/src/libra_py
)

#    COMMAND cmake -E copy ${CMAKE_SOURCE_DIR}/src/__init__.py ${CMAKE_BINARY_DIR}/src/__init__.py

