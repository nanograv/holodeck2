/**
 *
 */

#pragma once

// #include <array>
#include <format>

#include "../src/cosmology.h"
#include "../src/sam.h"
#include "tools.h"


void test__mmbulge_mstar_from_mbh() {
    printf(" - test_sam::test__mmbulge_mstar_from_mbh()\n");

    Cosmology cosmo;
    SAM sam(&cosmo);

    // ---- Test the units

    double mbh = 1E8;
    double lo = 1E6;
    double hi = 1E15;
    double test;

    const char* name = "`mmbulge_mstar_from_mbh` units";

    printf("Test %s...\n", name);
    test = sam.mmbulge_mstar_from_mbh(mbh);
    check_throw((test > lo), "Test 1 failed: mstar {:.2e} is too low  ( < {:.2e})!", test, lo);
    check_throw((test < hi), "Test 1 failed: mstar {:.2e} is too high ( > {:.2e})!", test, hi);
    printf("✅ %s passed.\n", name);

    //! FIX: once parameters can be changed, test known values

    std::cout << "test__mmbulge_mstar_from_mbh passed all tests.\n";
}


// void test__mmbulge_mstar_from_mbh() {
//     printf(" - test_sam::test__mmbulge_mstar_from_mbh()\n");

//     Cosmology cosmo;
//     SAM sam(&cosmo);

//     // ---- Test the units
//     double mbh = 1E8;
//     double lo = 1E6;
//     double hi = 1E15;
//     double test;

//     const char* name = "`mmbulge_mstar_from_mbh` units";

//     printf("Test %s...\n", name);
//     test = sam.mmbulge_mstar_from_mbh(mbh);
//     check_throw((test > lo), "Test 1 failed: mstar {:.2e} is too low  ( < {:.2e})!", test, lo);
//     check_throw((test < hi), "Test 1 failed: mstar {:.2e} is too high ( > {:.2e})!", test, hi);
//     printf("✅ %s passed.\n", name);

//     std::cout << "test__mmbulge_mstar_from_mbh passed all tests.\n";
// }


void test_sam() {
    printf(" ====    test_sam::test_sam()\n");

    test__mmbulge_mstar_from_mbh();

    std::cout << "All SAM tests passed.\n";

}



// std::string err_msg = std::format(
//     "Interpolation test failed on {}!: test={:.8e} vs. truth={:.8e} at redz={:.8e}",
//     msg, test, truth, redz
// );
// throw std::runtime_error(err_msg);
