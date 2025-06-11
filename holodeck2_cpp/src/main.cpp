/**
 *
 */

#include <stdio.h>
#include <chrono>

#include "constants.h"
#include "sam.h"
#include "utils.h"


int main(int argc, char *argv[]) {

    LOG_WARNING(get_logger(), " ====    holodeck2::main(): beginning    ====\n");

    // ---- Setup / Initialize PTA, GravWaves, and SAM

    Cosmology cosmo;

    // PTA* pta = new PTA(20.0, 30);
    PTA pta(15.0, 30);

    // GravWaves* gw = new GravWaves(&pta, 10, 5);
    GravWaves gw(&pta, 100, 5);

    // SAM* sam = new SAM(&cosmo);
    SAM sam(&cosmo);


    // ---- Run Calculations

    LOG_INFO(get_logger(), "Calculating gravitational waves...\n");



    int num_runs = 10;
    double dur;
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    double* durations = (double*)malloc(num_runs * sizeof(double));
    double ave_dur = 0.0, min_dur = 1e10, max_dur = 0.0;
    for (int i = 0; i < num_runs; ++i) {
        auto start = std::chrono::steady_clock::now();
        sam.grav_waves(pta, gw);
        auto end = std::chrono::steady_clock::now();
        // auto duration_ms = duration_cast<milliseconds>(end - start).count();
        dur = std::chrono::duration<double>(end - start).count();
        printf("%d: %.4e\n", i, dur);
        durations[i] = dur;
        ave_dur += dur / num_runs;
        min_dur = (dur < min_dur) ? dur : min_dur;
        max_dur = (dur > max_dur) ? dur : max_dur;
    }

    printf("Duration ave: %.2e, min: %.2e, max: %.2e\n", ave_dur, min_dur, max_dur);

    LOG_INFO(get_logger(), "Gravitational waves calculated.\n");


    // ---- Cleanup

    LOG_WARNING(get_logger(), " ====    holodeck2::main(): completed.    ====\n");

    free(durations);

}
