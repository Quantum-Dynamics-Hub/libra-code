#
# Set the libraries
# 
SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}/libpython2.6.dll )


#
#  Building sub-directories
#
MESSAGE("Going into subdirectory mmath...")
ADD_SUBDIRECTORY("mmath")

MESSAGE("Going into subdirectory mmath_eigen...")
ADD_SUBDIRECTORY("mmath_eigen")

MESSAGE("Going into subdirectory ann...")
ADD_SUBDIRECTORY("ann")

MESSAGE("Going into subdirectory rigidbody...")
ADD_SUBDIRECTORY("rigidbody")

MESSAGE("Going into subdirectory cell...")
ADD_SUBDIRECTORY("cell")

MESSAGE("Going into subdirectory pot...")
ADD_SUBDIRECTORY("pot")

MESSAGE("Going into subdirectory operators...")
ADD_SUBDIRECTORY("operators")

MESSAGE("Going into subdirectory hamiltonian...")
ADD_SUBDIRECTORY("hamiltonian")





