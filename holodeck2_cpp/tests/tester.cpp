/**
 *
 */

// #include "../src/utils.h"
#include "test_cosmology.h"
#include "test_tools.h"
#include "test_utils.h"

bool EXIT_ON_FAIL = true;


int main() {
    printf(" ====    tester::main()\n");

    printf("Running test_utils...\n");
    try {
        test_utils();
    } catch (const std::exception& e) {
        std::cerr << "❌ Exception in `test_utils!` : " << e.what() << std::endl;
        if (EXIT_ON_FAIL) return 1;
    }

    printf("Running test_cosmology...\n");
    try {
        int result = test_cosmology();
    } catch (const std::exception& e) {
        std::cerr << "❌ Exception in `test_cosmology!` : " << e.what() << std::endl;
        if (EXIT_ON_FAIL) return 1;
    }

    std::cout << "✅ All tests passed.\n";
    return 0;

}