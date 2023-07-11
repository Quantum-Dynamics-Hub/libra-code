/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef CONVERTERS_H
#define CONVERTERS_H

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#endif

#include "../chemobjects/libchemobjects.h"
#include "../dyn/libdyn.h"


/// liblibra namespace
namespace liblibra{

using namespace libchemobjects;
using namespace libdyn;


namespace libconverters{


void system_to_nuclear(System& syst, Nuclear& nucl);
void nuclear_to_system(Nuclear& nucl, System& syst);

void system_to_vector_q(System& syst, vector<double>& q);
void system_to_vector_p(System& syst, vector<double>& p);



typedef std::vector<std::string> StringList;
typedef std::vector<vector<std::string> > StringMap;


typedef std::map<std::string, double> StringDoubleMap;
typedef std::map<std::string, int> StringIntMap;
typedef std::map<std::string, vector<double> > StringVDoubleMap;
//typedef std::map<std::string, double> StringIntMap;
//boost::python::dict map_to_dict(std::map<std::string, double> map_);



template<class T>
vector<T> Py2Cpp(boost::python::list x){
/** 
  Converts the input Python list to vector<T> object
*/

  int sz = boost::python::len(x);
  vector<T> res(sz);

  for(int i=0;i<sz;i++){ 
    res[i] = boost::python::extract<T>(x[i]);
  }

  return res;
}

template<class T>
boost::python::list Cpp2Py(vector<T>& x){
/**
  Converts the vector<T> object into Python list
*/

  int sz = x.size();
  boost::python::list res;

  for(int i=0;i<sz;i++){ res.append(x[i]); }

  return res;
}





/// from : https://wiki.python.org/moin/boost.python/StlContainers
void IndexError();
void KeyError();

template<class T>
struct std_item
{
    typedef typename T::value_type V;
    static V& get(T const& x, int i)
    {
        if( i<0 ) i+=x.size();
        if( i>=0 && i<x.size() ) return x[i];
        IndexError();
    }
    static void set(T const& x, int i, V const& v)
    {
        if( i<0 ) i+=x.size();
        if( i>=0 && i<x.size() ) x[i]=v;
        else IndexError();
    }
    static void del(T const& x, int i)
    {
        if( i<0 ) i+=x.size();
        if( i>=0 && i<x.size() ) x.erase(i);
        else IndexError();
    }
    static void add(T const& x, V const& v)
    {
        x.push_back(v);
    }

    static int index(T const& x, V const& v)
    {
        int i=0;
        for(typename T::const_iterator it=x.begin; it!=x.end(); ++it,++i)
          if( *it == v ) return i;
        return -1;
    }
    static bool in(T const& x, V const& v)
    {
        return find_eq(x.begin, x.end, v) != x.end();
    }


};




template<class T>
struct map_item
{
    typedef typename T::key_type K;
    typedef typename T::mapped_type V;
    static V& get(T const& x, K const& i)
    {
        if( x.find(i) != x.end() ) return x[i];
        KeyError();
    }
    static void set(T const& x, K const& i, V const& v)
    {
        x[i]=v; // use map autocreation feature
    }
    static void del(T const& x, K const& i)
    {
        if( x.find(i) != x.end() ) x.erase(i);
        else KeyError();
    }

    static int index(T const& x, K const& k)
    {
        int i=0;
        for(typename T::const_iterator it=x.begin(); it!=x.end(); ++it,++i){
          if( it->first == k ) return i;        
        }
        return -1;
    }
 
    static bool in(T const& x, K const& i)
    {
        return x.find(i) != x.end();
    }


};

template<class T>
static boost::python::list keys(T const& x){
  boost::python::list t;
  for(typename T::const_iterator it=x.begin(); it!=x.end(); ++it){  t.append(it->first); }

  return t;
}

template<class T>
static boost::python::list values(T const& x){
  boost::python::list t;
 
  for(typename T::const_iterator it=x.begin(); it!=x.end(); ++it){  t.append(it->second); }

  return t;
}

template<class T>
static boost::python::list items(T const& x){
  boost::python::list t;
  for(typename T::const_iterator it=x.begin(); it!=x.end(); ++it){  t.append(make_tuple(it->first,it->second)); }

  return t;
}





}// namespace libconverters
}// liblibra

#endif // CONVERTERS_H

