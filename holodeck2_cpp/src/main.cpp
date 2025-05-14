/**
 *
 */

//  #include <iostream>
//  #include <deque>
#include <stdio.h>
//  #include <string>

// #include "test.h"

#include "sam.h"
#include "constants.h"


int main(int argc, char *argv[]) {
    printf("main()\n");

    // ---- Setup / Initialize PTA, GravWaves, and SAM

    PTA* pta = new PTA(20.0, 30);
    printf(
        "Initialized PTA with %d frequencies between [%.2e, %.2e] Hz\n",
        pta->num_freqs, pta->fobs_cents[0], pta->fobs_edges[pta->num_freqs-1]
    );

    GravWaves* gw = new GravWaves(pta, 10, 5);
    printf(
        "Initialized GravWaves with %d freqs, %d reals, %d louds\n",
        gw->num_freqs, gw->num_reals, gw->num_louds
    );

    SAM* sam = new SAM();
    printf(
        "Initialized SAM with %d mass bins, %d redshift bins\n",
        sam->num_mass, sam->num_redz
    );

    // ---- Run Calculations

    printf("Calculating gravitational waves...\n");
    grav_waves(pta, gw, sam);

    printf("Done.\n\n");

    // ---- Cleanup

    delete(pta);
    delete(gw);
    delete(sam);
}
