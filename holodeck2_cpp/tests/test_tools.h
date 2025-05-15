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

