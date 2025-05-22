/**
 *
 */

// #include "../src/utils.h"
#include "test_cosmology.h"
#include "test_sam.h"
#include "test_utils.h"
#include "tools.h"

bool EXIT_ON_FAIL = true;

#define TEST_COSMOLOGY
#define TEST_SAM
#define TEST_UTILS


int main() {
    printf(" ====    tester::main()\n");

    #ifdef TEST_COSMOLOGY
    printf("Running test_cosmology...\n");
    try {
        int result = test_cosmology();
        printf("✅ `test_cosmology()` passed.\n");
    } catch (const std::exception& e) {
        std::cerr << "❌ Exception in `test_cosmology!` : " << e.what() << std::endl;
        if (EXIT_ON_FAIL) return 1;
    }
    #endif

    #ifdef TEST_SAM
    printf("Running test_sam...\n");
    try {
        test_sam();
        printf("✅ `test_sam()` passed.\n");
    } catch (const std::exception& e) {
        std::cerr << "❌ Exception in `test_sam!` : " << e.what() << std::endl;
        if (EXIT_ON_FAIL) return 1;
    }
    #endif

    #ifdef TEST_UTILS
    printf("Running test_utils...\n");
    try {
        test_utils();
        printf("✅ `test_utils()` passed.\n");
    } catch (const std::exception& e) {
        std::cerr << "❌ Exception in `test_utils!` : " << e.what() << std::endl;
        if (EXIT_ON_FAIL) return 1;
    }
    #endif

    std::cout << "All tests passed.\n";
    return 0;

}