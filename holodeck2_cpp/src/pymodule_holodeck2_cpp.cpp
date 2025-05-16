/**
 *
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>
// #include "test.h"


typedef struct {
    PyObject_HEAD
    // Test *test;
} CustomObject;


static PyTypeObject CustomType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "holodeck2_cpp.Test",
    .tp_doc = PyDoc_STR("Test object."),
    .tp_basicsize = sizeof(CustomObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = PyType_GenericNew,
};


// =============================================================================
// ====    Define Functions    =====
// =============================================================================


static PyObject* custom_function(PyObject *self, PyObject *args) {
    int argOne;
    int argTwo;
    double argThree = 0.0;

    if (!PyArg_ParseTuple(args, "ii|f", &argOne, &argTwo, &argThree))
        return NULL;

    return PyFloat_FromDouble(argOne + argTwo + argThree);
}


// =============================================================================
// ====    Define Module    ====
// =============================================================================


static PyMethodDef holodeck2_cpp_methods[] = {
    {
        "custom_function",             // Name (str)
        custom_function,               // Function (pointer)
        METH_VARARGS,                  // Argument flags (typically:  METH_VARARGS or METH_VARARGS | METH_KEYWORDS)
        "A custom function in C."      // Description (str)
    },
    {NULL, NULL, 0, NULL}              // Sentinel to designate end of methods.
};


static PyModuleDef holodeck2_cpp_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "holodeck2cpp",
    .m_doc = "Example module that creates an extension type.",
    .m_size = -1,
    .m_methods = holodeck2_cpp_methods,
};


PyMODINIT_FUNC PyInit_holodeck2_cpp(void) {

    fprintf(stderr, "Running `PyInit_holodeck2_cpp`\n");

    if (PyType_Ready(&CustomType) < 0) return NULL;

    PyObject *m = PyModule_Create(&holodeck2_cpp_module);
    if (m == NULL) return NULL;

    if (PyModule_AddObjectRef(m, "Custom", (PyObject *) &CustomType) < 0) {
        Py_DECREF(m);
        return NULL;
    }

    return m;
}


