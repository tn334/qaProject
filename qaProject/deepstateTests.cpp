#include <deepstate/DeepState.hpp>
//#include "qaProject.cpp"
#include "qaProject.h"
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <omp.h>

using namespace deepstate;

TEST(QAProject, SequentialVsParallel) {
    int N = DeepState_IntInRange(PARAMTER_START, PARAMETER_END);
    int M = DeepState_IntInRange(PARAMTER_START, PARAMETER_END);
    int threads = DeepState_IntInRange(MINTHREADS, MAXTHREADS);

    omp_set_num_threads(threads);
    //double tolerance = DeepState_Double(timeDiffThreshold);
    double seqTime = 0, paraTime = 0;
    const char* filename = "output.txt";

    LOG(TRACE) << "N=" << N;
    LOG(TRACE) << "M=" << M;

    ASSERT_GE(N, 50); // Assert that N is greater than or equal to 0
    ASSERT_LE(N, 300);
    ASSERT_GE(M, 50); // Assert that N is greater than or equal to 0
    ASSERT_LE(M, 300);

    // Test Threads is right number
    ASSERT_GE(threads, MINTHREADS);
    ASSERT_LE(threads, MAXTHREADS)
    //ASSERT_LT(loseRate, 1.0); // Assert that loseRate is less than 1.0

    //ASSERT_EQ(totalRuns, (parameterEnd / incrementor) * (parameterEnd / incrementor));

    inLoopRun(N, M, seqTime, paraTime);

    // Run the sequential implementation
    double seqResult = sequentialRun(N, M);

    // Run the parallel implementation
    double paraResult = parrallelOptimizedRun(N, M);



    ASSERT_EQ(seqResult, paraResult) << "Sequential and parallel results should be equal.";


    ASSERT_EQ(firstSeqLoop, firstParallelLoop);
    // Check if the results are equal within a tolerance


    // Add more assertions if needed
    seqTimer(N, M, seqTime);
    paraTimer(N, M, paraTime);
    LOG(TRACE) << "Parallelized Time=" << paraTime;
    LOG(TRACE) << "Sequential Time=" << seqTime;
    //ASSERT_LT(seqTime, paraTime) << "Sequential time should be less than parallel time.";
    
    
    ASSERT_FLOAT_NEAR(seqTime, paraTime, timeDiffThreshold) << "Values are not within tolerance.";

    // You can also check specific conditions based on your implementation
    // For example, you might want to check that the sequential implementation is faster
    // in certain cases or that the parallel implementation produces correct results.
}
