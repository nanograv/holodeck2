/**
 *
 */

#pragma once

// #include <array>
#include <format>

#include "../src/utils.h"
#include "tools.h"


void test_argsort() {
    printf(" - test_utils::test_argsort()\n");

    double a[] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    int truth1[] = { 0, 1, 2, 3, 4 };

    // ---- Test 1: Already-Sorted array

    int* sorted_indices = utils::argsort(a, 5);
    for (int i = 0; i < 5; ++i) {
        // std::string err_msg = std::format(
        //     "Interpolation test failed on {}!: test={:.8e} vs. truth={:.8e} at redz={:.8e}",
        //     msg, test, truth, redz
        // );
        check_throw(sorted_indices[i] == truth1[i], "argsort() failed!");
    }
    std::cout << "✅ Already-sorted array passed.\n";

    // ---- Test 2: reverse-sorted array

    double b[] = { 5.0, 4.0, 3.0, 2.0, 1.0 };
    int truth2[] = { 4, 3, 2, 1, 0 };

    sorted_indices = utils::argsort(b, 5);
    for (int i = 0; i < 5; ++i) {
        check_throw(sorted_indices[i] == truth2[i], "argsort() failed!");
    }
    std::cout << "✅ Reverse-sorted array passed.\n";


    // ---- Test 3: random order array
    double c[] = { 3.0, 1.0, 4.0, 2.0, 5.0 };
    int truth3[] = { 1, 3, 0, 2, 4 };
    sorted_indices = utils::argsort(c, 5);
    for (int i = 0; i < 5; ++i) {
        check_throw(sorted_indices[i] == truth3[i], "argsort() failed!");
    }
    std::cout << "✅ Random array passed.\n";


    delete[] sorted_indices;
}


void test_index_2d_to_1d() {
    printf(" - test_utils::test_index_2d_to_1d()\n");

    int dim1 = 5, dim2 = 3;
    int index, a, b;
    for (int i = 0; i < dim1; ++i) {
        for (int j = 0; j < dim2; ++j) {
            utils::index_2d_to_1d(i, j, dim1, dim2, &index);
            check_throw(index == (i * dim2 + j), "index_2d_to_1d() failed!");
            utils::index_1d_to_2d(index, dim1, dim2, &a, &b);
            check_throw(i == a, "index_1d_to_2d() failed!  i != a");
            check_throw(j == b, "index_1d_to_2d() failed!  j != b");
        }
    }

    std::cout << "✅ index_2d_to_1d passed.\n";
}


void test_index_3d_to_1d() {
    printf(" - test_utils::test_index_3d_to_1d()\n");

    int dim1 = 5, dim2 = 3, dim3 = 4;
    int index, a, b, c;
    for (int i = 0; i < dim1; ++i) {
        for (int j = 0; j < dim2; ++j) {
            for (int k = 0; k < dim3; ++k) {
                utils::index_3d_to_1d(i, j, k, dim1, dim2, dim3, &index);
                check_throw(index == (i * (dim2*dim3) + j*dim3 + j), "index_3d_to_1d() failed!");
                utils::index_1d_to_3d(index, dim1, dim2, dim3, &a, &b, &c);
                check_throw(i == a, "index_1d_to_3d() failed!  i != a");
                check_throw(j == b, "index_1d_to_3d() failed!  j != b");
                check_throw(k == c, "index_1d_to_3d() failed!  k != c");
            }
        }
    }

    std::cout << "✅ index_3d_to_1d passed.\n";
}


void test_is_almost_equal() {
    printf(" - test_utils::test_is_almost_equal()\n");

    // Test the assert_almost_equal function
    double a = 1.0f;
    double diff = 1e-8f;
    double b = a + diff;

    printf("Test atol\n");
    check_throw(utils::is_almost_equal(a, b, 1.01*diff, 0.0), "Test 1 failed: a and b should be almost equal");
    check_throw(!utils::is_almost_equal(a, b, 0.99*diff, 0.0), "Test 2 failed: a and b should not be almost equal");
    std::cout << "✅ assert_almost_equal passed 'atol' test.\n";

    printf("Test rtol\n");
    b = a + diff * a;
    check_throw(utils::is_almost_equal(a, b, 0.0, 1.01*diff), "Test 3 failed: a and b should be almost equal");
    check_throw(!utils::is_almost_equal(a, b, 0.0, 0.99*diff), "Test 4 failed: a and b should not be almost equal");

    std::cout << "✅ assert_almost_equal passed 'rtol' test.\n";
}


void test_quantiles() {
    printf(" - test_utils::test_quantiles()\n");

    double vals0[]   = {0.00, 1.00};
    double percs0[]  = {0.00, 0.25, 0.50, 0.75, 1.00};
    double quants0[] = {0.00, 0.25, 0.50, 0.75, 1.00};
    double* test0 = utils::quantiles(vals0, 2, percs0, 5, nullptr, true);
    for (int i = 0; i < 5; ++i) {
        // printf("test1[%d]=%.8e vs. quants1[%d]=%.8e\n", i, test1[i], i, quants1[i]);
        check_throw(
            utils::is_almost_equal(test0[i], quants0[i], 0.0, 1E-4),
            "`quantiles()` faield test 0!!  perc[{}]={:.6e} ==> {:.6e}, should be {:.6e}!",
            i, percs0[i], test0[i], quants0[i]
        );
    }
    printf("✅ `quantiles()` test 0 passed.\n");

    int size1        = 11;
    double vals1[]   = {0.00, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 10.0};
    double percs1[]  = {0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00};
    double quants1[] = {0.00, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 10.0};
    double* test1 = utils::quantiles(vals1, size1, percs1, size1, nullptr, true);
    for (int i = 0; i < size1; ++i) {
        // printf("test1[%d]=%.8e vs. quants1[%d]=%.8e\n", i, test1[i], i, quants1[i]);
        check_throw(
            utils::is_almost_equal(test1[i], quants1[i], 0.0, 1E-4),
            "quantiles() failed test1!!  perc[{}]={:.6e} ==> {:.6e}, should be {:.6e}!",
            i, percs1[i], test1[i], quants1[i]
        );
    }
    printf("✅ `quantiles()` test 1 passed.\n");

    // Test the quantiles function
    double vals2[]    = {-1.7419e-01, -1.2005e-01, -9.5085e-02, +1.0767e-02, +2.2249e-01};
    double percs2[]   = {+4.8312e-01, +1.3999e-01, +1.1176e-01, +7.4080e-01, +2.9461e-01};
    double quants2[]  = {-9.6771e-02, -1.4387e-01, -1.4999e-01, +6.8729e-03, -1.1560e-01};
    double* test2 = utils::quantiles(vals2, 5, percs2, 5, nullptr, true);
    for (int i = 0; i < 5; ++i) {
        check_throw(
            utils::is_almost_equal(test2[i], quants2[i], 0.0, 1E-3),
            "quantiles() failed test2!!  perc[{}]={:.6e} ==> {:.6e}, should be {:.6e}!",
            i, percs2[i], test2[i], quants2[i]
        );
    }
    printf("✅ `quantiles()` test 2 passed.\n");

    std::cout << "quantiles passed all tests.\n";

}


void test_utils() {
    printf(" ====    test_utils::test_utils()\n");

    test_argsort();

    test_index_2d_to_1d();

    test_is_almost_equal();

    test_quantiles();

    std::cout << "All utils tests passed.\n";

}



// std::string err_msg = std::format(
//     "Interpolation test failed on {}!: test={:.8e} vs. truth={:.8e} at redz={:.8e}",
//     msg, test, truth, redz
// );
// throw std::runtime_error(err_msg);
