#ifndef PY_COPY_H
#define PY_COPY_H

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <stdint.h>
#define PYTHON_ERROR(TYPE, REASON) \
{ \
    PyErr_SetString(TYPE, REASON); \
    throw bp::error_already_set(); \
}


namespace bp = boost::python;
//using namespace boost;

template<class T>
inline PyObject * managingPyObject(T *p)
{
    return typename bp::manage_new_object::apply<T *>::type()(p);
}

template<class Copyable>
bp::object
generic__copy__(bp::object copyable)
{
    Copyable *newCopyable(new Copyable(bp::extract<const Copyable
&>(copyable)));
    bp::object
result(bp::detail::new_reference(managingPyObject(newCopyable)));

    bp::extract<bp::dict>(result.attr("__dict__"))().update(
        copyable.attr("__dict__"));

    return result;
}

template<class Copyable>
bp::object
generic__deepcopy__(bp::object copyable, bp::dict memo)
{
    bp::object copyMod = bp::import("copy");
    bp::object deepcopy = copyMod.attr("deepcopy");

    Copyable *newCopyable(new Copyable(bp::extract<const Copyable
&>(copyable)));
    bp::object
result(bp::detail::new_reference(managingPyObject(newCopyable)));

    // HACK: copyableId shall be the same as the result of id(copyable)
    //in Python -
    // please tell me that there is a better way! (and which ;-p)
//    int copyableId = (int)(copyable.ptr());
    uintptr_t copyableId = (uintptr_t)(copyable.ptr()); // This is needed to compile on 64-bit machines
    // see more: http:\//stackoverflow.com/questions/153065/converting-a-pointer-into-an-integer
    memo[copyableId] = result;

    bp::extract<bp::dict>(result.attr("__dict__"))().update(
        deepcopy(bp::extract<bp::dict>(copyable.attr("__dict__"))(),
memo));

    return result;
}


#endif // PY_COPY_H
