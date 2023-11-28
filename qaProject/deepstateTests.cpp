#include <deepstate/DeepState.hpp>
//#include "qaProject.cpp"
#include "qaProject.h"
#include <fcntl.h>
#include <unistd.h>
#include <string>

using namespace deepstate;

TEST(QAProject, SequentialVsParallel) {
    int N = DeepState_IntInRange(parameterStart, parameterEnd);
    int M = DeepState_IntInRange(parameterStart, parameterEnd);
    //double tolerance = DeepState_Double(timeDiffThreshold);
    double seqTime = 0, paraTime = 0;
    const char* filename = "output.txt";

    ASSERT_GE(N, 50); // Assert that N is greater than or equal to 0
    ASSERT_GE(M, 50); // Assert that N is greater than or equal to 0

    //ASSERT_LT(loseRate, 1.0); // Assert that loseRate is less than 1.0

    //ASSERT_EQ(totalRuns, (parameterEnd / incrementor) * (parameterEnd / incrementor));

    inLoopRun(N, M, seqTime, paraTime, )

    // Run the sequential implementation
    double seqResult = sequentialRun(N, M);

    // Run the parallel implementation
    double paraResult = parrallelOptimizedRun(N, M);

    ASSERT_EQ(seqResult, paraResult) << "Sequential and parallel results should be equal.";



    // Check if the results are equal within a tolerance


    // Add more assertions if needed
    ASSERT_LT(seqTime, paraTime) << "Sequential time should be less than parallel time.";
    seqTimer(N, M, seqTime);
    paraTimer(N, M, paraTime);

    //ASSERT_FLOAT_NEAR(seqTime, paraTime, timeDiffThreshold) << "Values are not within tolerance.";

    // You can also check specific conditions based on your implementation
    // For example, you might want to check that the sequential implementation is faster
    // in certain cases or that the parallel implementation produces correct results.
}
