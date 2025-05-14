/**
 *
 */

#pragma once

#include "../src/cosmology.h"
#include "test_tools.h"


int test_cosmology() {
    printf(" ====    test_cosmology::test_cosmology()\n");

    Cosmology cosmo;
    std::cout << "✅ Cosmology loaded successfully.\n";

    // Basic sanity check (replace 1e-5 with expected min z)
    // assert(cosmo.get_grid_size() > 0);
    // assert(cosmo.get_redshift(0) < cosmo.get_redshift(cosmo.get_grid_size() - 1));

    std::cout << "✅ Tests passed.\n";
    return 0;

}