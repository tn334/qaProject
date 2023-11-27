#include <deepstate/DeepState.hpp>
#include "qaProject.cpp"

using namespace deepstate;

TEST(QAProject, SequentialVsParallel) {
    int N = DeepState_IntInRange(parameterStart, parameterEnd);
    int M = DeepState_IntInRange(parameterStart, parameterEnd);
    double tolerance = .0005;
    // Run the sequential implementation
    double seqResult = sequentialRun(N, M);

    // Run the parallel implementation
    double paraResult = parrallelOptimizedRun(N, M);

    // Check if the results are equal within a tolerance
    ASSERT_EQ(seqResult, paraResult) << "Sequential and parallel results should be equal.";

    // Add more assertions if needed
    //ASSERT_LT(seqTime, paraTime) << "Sequential time should be less than parallel time.";

    //ASSERT_NEAR(expectedValue, actualValue, tolerance) << "Values are not within tolerance.";

    // You can also check specific conditions based on your implementation
    // For example, you might want to check that the sequential implementation is faster
    // in certain cases or that the parallel implementation produces correct results.
}
