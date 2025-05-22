/**
 *
 */

#include <stdio.h>

#include "constants.h"
#include "sam.h"
#include "utils.h"


int main(int argc, char *argv[]) {
    LOG_WARNING(get_logger(), " ====    holodeck2::main(): beginning    ====\n");

    // ---- Setup / Initialize PTA, GravWaves, and SAM

    Cosmology cosmo;
    LOG_DEBUG(get_logger(), "Initialized cosmology instance.\n");

    // PTA* pta = new PTA(20.0, 30);
    PTA pta(20.0, 30);
    LOG_DEBUG(get_logger(),
        "Initialized PTA with {} frequencies between [{:.2e}, {:.2e}] Hz\n",
        pta.num_freq_cents, pta.fobs_cents[0], pta.fobs_edges[pta.num_freq_cents-1]
    );

    // GravWaves* gw = new GravWaves(&pta, 10, 5);
    GravWaves gw(&pta, 10, 5);
    LOG_DEBUG(get_logger(),
        "Initialized GravWaves with {} freqs, {} reals, {} louds\n",
        gw.num_freq_cents, gw.num_reals, gw.num_louds
    );

    // SAM* sam = new SAM(&cosmo);
    SAM sam(&cosmo);
    LOG_DEBUG(get_logger(),
        "Initialized SAM with {} mass bins, {} redshift bins\n",
        sam.num_mass_edges, sam.num_redz_edges
    );

    // ---- Run Calculations

    LOG_INFO(get_logger(), "Calculating gravitational waves...\n");
    sam.grav_waves(pta, gw);

    LOG_WARNING(get_logger(), " ====    holodeck2::main(): completed.    ====\n");

    // ---- Cleanup

}
