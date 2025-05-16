/**
 *
 */

#pragma once

#include <array>
#include <format>

#include "../src/cosmology.h"
#include "test_tools.h"


void test_interp_values(double redz, double truth, double test, const char* msg) {
    if (!utils::is_almost_equal(test, truth, 0.0, 1E-4)) {
        std::string err_msg = std::format(
            "Interpolation test failed on {}!: test={:.8e} vs. truth={:.8e} at redz={:.8e}",
            msg, test, truth, redz
        );
        throw std::runtime_error(err_msg);
    }
}


void test_interpolation(Cosmology &cosmo) {
    printf(" - test_cosmology::test_interpolation()\n");

    double redz[]  = {
        1.80422089e-05, 2.34300168e-05, 2.45124079e-05, 2.97102841e-05, 8.72639006e-05,
        4.43186134e-03, 9.82815660e-03, 1.40178590e+00, 6.29892108e+00, 8.64871553e+02
    };
    double scafa[] = {
        9.99981958e-01, 9.99976571e-01, 9.99975488e-01, 9.99970291e-01, 9.99912744e-01,
        9.95587693e-01, 9.90267496e-01, 4.16356846e-01, 1.37006551e-01, 1.15490571e-03
    };
    double dcom[]  = {
        7.80166894e-02, 1.01314103e-01, 1.05994463e-01, 1.28470531e-01, 3.77333999e-01,
        1.91456138e+01, 4.24078764e+01, 4.25567300e+03, 8.56800310e+03, 1.39825572e+04
    };
    double vcom[]  = {
        1.98907501e-03, 4.35610486e-03, 4.98813438e-03, 8.88176278e-03, 2.25043474e-01,
        2.93965587e+04, 3.19468606e+05, 3.22844475e+11, 2.63467744e+12, 1.14511320e+13
    };
    double tlook[] = {
        2.54454113e-01, 3.30438539e-01, 3.45703464e-01, 4.19008607e-01, 1.23064521e+00,
        6.23067191e+01, 1.37641198e+02, 9.17562386e+03, 1.28651872e+04, 1.37520453e+04
    };
    double efunc[] = {
        1.00000779e+00, 1.00001012e+00, 1.00001059e+00, 1.00001284e+00, 1.00003770e+00,
        1.00192122e+00, 1.00427848e+00, 2.16845714e+00, 1.06159892e+01, 1.36733979e+04
    };
    double dtdz[]  = {
        4.45059550e+17, 4.45056116e+17, 4.45055426e+17, 4.45052114e+17, 4.45015436e+17,
        4.42257586e+17, 4.38861734e+17, 8.54563246e+16, 5.74394415e+15, 3.75923453e+10
    };

    int num_test = std::size(redz);

    double test;

    for (int ii = 0; ii < num_test; ii++) {
        // `dcom` - comoving distance
        cosmo.dcom_from_redz(redz[ii], &test);
        test_interp_values(redz[ii], dcom[ii], test, "dcom");
    }
    std::cout << "✅ dcom interpolation successful.\n";

    for (int ii = 0; ii < num_test; ii++) {
        // `vcom` - comoving volume
        cosmo.vcom_from_redz(redz[ii], &test);
        test_interp_values(redz[ii], vcom[ii], test, "vcom");
    }
    std::cout << "✅ vcom interpolation successful.\n";

    for (int ii = 0; ii < num_test; ii++) {
        // `tlook` - lookback time
        cosmo.tlook_from_redz(redz[ii], &test);
        test_interp_values(redz[ii], tlook[ii], test, "tlook");
    }
    std::cout << "✅ tlook interpolation successful.\n";

    for (int ii = 0; ii < num_test; ii++) {
        // `tlook` - lookback time
        cosmo.efunc_from_redz(redz[ii], &test);
        test_interp_values(redz[ii], efunc[ii], test, "efunc");
    }
    std::cout << "✅ efunc interpolation successful.\n";

    for (int ii = 0; ii < num_test; ii++) {
        // `tlook` - lookback time
        cosmo.dtdz_from_redz(redz[ii], &test);
        test_interp_values(redz[ii], dtdz[ii], test, "dtdz");
    }
    std::cout << "✅ dtdz interpolation successful.\n";

}


int test_cosmology() {
    printf(" ====    test_cosmology::test_cosmology()\n");

    Cosmology cosmo;
    std::cout << "✅ Cosmology loaded successfully.\n";

    test_interpolation(cosmo);

    std::cout << "✅ Tests passed.\n";
    return 0;

}


