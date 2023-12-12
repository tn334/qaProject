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

    LOG(TRACE) << "N = " << N;
    LOG(TRACE) << "M = " << M;
    LOG(TRACE) << "Number of Threads = " << threads;

    ASSERT_GE(N, 50); // Assert that N is greater than or equal to 50
    ASSERT_LE(N, 300);
    ASSERT_GE(M, 50); // Assert that N is greater than or equal to 50
    ASSERT_LE(M, 300);

    // Test Threads is right number
    ASSERT_GE(threads, MINTHREADS);
    ASSERT_LE(threads, MAXTHREADS);
    //ASSERT_LT(loseRate, 1.0); // Assert that loseRate is less than 1.0

    //ASSERT_EQ(totalRuns, (parameterEnd / incrementor) * (parameterEnd / incrementor));

    inLoopRun(N, M, seqTime, paraTime);

    // Run the sequential implementation
    double seqResult = sequentialRun(N, M);
    LOG(TRACE) << "Sequential Result = " << seqResult;
    // Run the parallel implementation
    double paraResult = parrallelOptimizedRun(N, M);
    LOG(TRACE) << "Parallel Result = " << paraResult;



    ASSERT_EQ(seqResult, paraResult); // "Sequential and parallel results should be equal.";


    //ASSERT_EQ(firstSeqLoop, firstParallelLoop);
    // Check if the results are equal within a tolerance


    // Add more assertions if needed
    seqTimer(N, M, seqTime);
    paraTimer(N, M, paraTime);
    
    LOG(TRACE) << "Sequential Time = " << seqTime;
    LOG(TRACE) << "Parallelized Time = " << paraTime;
    // May make crashes
    ASSERT_GT(seqTime, paraTime); // "Sequential time should be less than parallel time.";
    ASSERT_FLOAT_NEAR(seqTime, paraTime, .005); //Asserting that we are within .005
}
