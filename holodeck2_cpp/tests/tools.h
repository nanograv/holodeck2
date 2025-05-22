/**
 *
 */

#pragma once

#include <iostream>
#include <cassert>
#include <stdexcept>


using namespace std;


template <typename... Args>
void check_throw(bool condition, string_view fmt_str, Args&&... args) {
    if (!condition) {
        auto message = vformat(fmt_str, make_format_args(std::forward<Args>(args)...));
        throw runtime_error(message);
    }
}