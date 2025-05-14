/**
 *
 */

#pragma once

#include <iostream>
#include <cassert>
#include <stdexcept>


using namespace std;


void check(bool cond, const string& msg) {
    if (!cond) throw std::runtime_error("\n" + msg);
}


inline bool is_almost_equal(double a, double b, double atol = 1e-8f, double rtol = 1e-6f) {
    double rval;
    if ((fabs(a) > 0) && (fabs(b) > 0)) {
        rval = fmin(fabs(a), fabs(b));
    } else {
        rval = fmax(fabs(a), fabs(b));
    }
    // printf("is_almost_equal(%e, %e, %e, %e) :: rval=%e, diff=%e, rdiff=%e\n", a, b, atol, rtol, rval, fabs(a-b), rtol * rval);
    // bool rv = (fabs(a - b) <= (atol + rtol * rval));
    // printf("rv = %i\n", (int)rv);
    // return rv;
    return (fabs(a - b) <= (atol + rtol * rval));
}