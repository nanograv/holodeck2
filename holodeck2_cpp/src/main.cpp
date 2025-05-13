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

    PTA* pta = new PTA(20.0, 30);
    printf(
        "Initialized PTA with %d frequencies between [%.2e, %.2e] Hz\n",
        pta->numFreqs, pta->fobsCents[0], pta->fobsEdges[pta->numFreqs-1]
    );

    GravWaves* gw = new GravWaves(pta, 10, 5);
    printf(
        "Initialized GravWaves with %d freqs, %d reals, %d louds\n",
        gw->numFreqs, gw->numReals, gw->numLouds
    );

    SAM* sam = new SAM();
    printf(
        "Initialized SAM with %d mass bins, %d redshift bins\n",
        sam->numMass, sam->numRedz
    );

    delete(pta);
    delete(gw);
    delete(sam);
}
