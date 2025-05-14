/**
 *
 */

#include "test_tools.h"
#include "test_cosmology.h"

int main() {
    printf(" ====    tester::main()\n");

    // Test the assert_almost_equal function
    double a = 1.0f;
    double diff = 1e-8f;
    double b = a + diff;

    printf("Test atol\n");
    check(is_almost_equal(a, b, 1.01*diff, 0.0), "Test 1 failed: a and b should be almost equal");
    check(!is_almost_equal(a, b, 0.99*diff, 0.0), "Test 2 failed: a and b should not be almost equal");

    printf("Test rtol\n");
    b = a + diff * a;
    check(is_almost_equal(a, b, 0.0, 1.01*diff), "Test 3 failed: a and b should be almost equal");
    check(!is_almost_equal(a, b, 0.0, 0.99*diff), "Test 4 failed: a and b should not be almost equal");
    std::cout << "✅ assert_almost_equal passed.\n";

    // Add more tests as needed
    printf("Running test_cosmology...\n");
    try {
        int result = test_cosmology();
    } catch (const std::exception& e) {
        std::cerr << "Exception in `test_cosmology!` : " << e.what() << std::endl;
        return 1;
    }


    // std::cout << "✅ All tests passed.\n";
    return 0;

    // try {
    //     test_tools_main();  // put all test logic here
    // } catch (const std::exception& e) {
    //     std::cerr << "❌ Uncaught exception: " << e.what() << "\n";
    //     return 1;

}