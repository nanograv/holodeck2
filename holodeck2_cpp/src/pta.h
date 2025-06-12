/**
 *
 *
 */

#include "config.h"
#include "utils.h"
#include "constants.h"


class PTA {

public:
    double obs_dur_yr;     // [yr]
    int num_freq_cents;
    double* fobs_cents;    // [Hz]
    double* fobs_edges;    // [Hz]

    PTA(double dur=20.0, double nfreqs=30) : obs_dur_yr(dur), num_freq_cents(nfreqs) {
        fobs_cents = (double*)malloc(num_freq_cents * sizeof(double));
        fobs_edges = (double*)malloc((num_freq_cents+1) * sizeof(double));
        double df = 1.0 / (obs_dur_yr * YR);
        fobs_edges[0] = df * 0.5;
        for (int ii = 0; ii < num_freq_cents; ii++) {
            fobs_cents[ii] = df * (ii + 1);
            fobs_edges[ii+1] = df * (ii + 1.5);
        }
        LOG_DEBUG(get_logger(),
            "Initialized PTA with {} frequencies between [{:.2e}, {:.2e}] Hz = [{:.2e}, {:.2e}] yr^-1\n",
            num_freq_cents, fobs_cents[0], fobs_edges[num_freq_cents-1],
            fobs_cents[0]*YR, fobs_edges[num_freq_cents-1]*YR
        );
    };

    ~PTA() {
        free(fobs_cents);
        free(fobs_edges);
    };

};
