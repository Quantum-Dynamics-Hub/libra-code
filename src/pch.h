#ifndef PCH_H
#define PCH_H



#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#ifdef __OPENMP
#include <omp.h>
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime> 

#include <complex>
#include <array>
#include <vector>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <chrono>
#include <thread>
#include <algorithm>
#include <map>
#include <memory> // for std::auto_ptr<>
#include <set>
#include <exception>


// Eigen matrix algebra library
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <Eigen/SVD>


// Libint Gaussian integrals library
//#include <libint2.hpp>

// Boost
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
#include <boost/regex.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>


#define BOOST_PYTHON_MAX_ARITY 30


#endif // PCH_H
