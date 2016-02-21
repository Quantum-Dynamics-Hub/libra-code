#
# Set the libraries
# 
SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}/libpython2.6.dll )
#SET( EIGEN_INCLUDE /home/Alexey_2/Soft/eigen-eigen-2249f9c22fe8/ )


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










