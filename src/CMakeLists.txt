#
# Set the libraries
# 
SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}/libpython2.6.dll )


#
#  Building sub-directories
#
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


MESSAGE("Going into subdirectory hamiltonian...")
ADD_SUBDIRECTORY("hamiltonian")

MESSAGE("Going into subdirectory dyn...")
ADD_SUBDIRECTORY("dyn")




#MESSAGE("Going into subdirectory cell...")
#ADD_SUBDIRECTORY("cell")

#MESSAGE("Going into subdirectory pot...")
#ADD_SUBDIRECTORY("pot")








