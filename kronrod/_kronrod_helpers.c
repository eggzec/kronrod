/**
 * @file _kronrod_helpers.c
 * @brief Helper wrappers for f2py Python bindings.
 *
 * The kronrod() function returns arrays via output pointers; this thin
 * wrapper packs the results into pre-allocated arrays suitable for f2py.
 */

#include "kronrod.h"
