// qaProject.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include "qaProject.h"
#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <string>
#include <vector>
#include <tuple>
#include <fstream> //ofstream, instream, .eof(), .get()

using namespace std;

int main()
{
    int threads = 2, maxThreads = 5;
    //open output file
    ofstream outFile;
    //try for threads 2-5???

    //cout << "Opening output.txt file" << endl;
    outFile.open("../output.txt");
    assert(outFile.is_open() && "Failed to open output.txt file");

    //find a way to test threads numbers
    for (threads = 2; threads < maxThreads; threads++)
    {
        // you may add some code here to measure the execution time
        double seqTime = 0, paraTime = 0;
        double seqWins = 0, totalRuns = 0;
        double loseRate = 0;
        vector<tuple<int, int, double, double>> sequentialFasterCases;

        omp_set_num_threads(threads);
        //sequential
        for (int N = PARAMTER_START; N <= PARAMETER_END; N += INCREMENTOR)
        {
            for (int M = PARAMTER_START; M <= PARAMETER_END; M += INCREMENTOR)
            {
                inLoopRun(N, M, seqTime, paraTime);
            }
        }
        loseRate = seqWins / totalRuns;

        outFile << "\n\nPercentage of Sequential Wins for " << threads << " is " << loseRate << endl;

        // handle sequential wins
        for (const auto& caseData : sequentialFasterCases)
        {
            int N_case, M_case;
            double Seq_case, Para_case;
            tie(N_case, M_case, Seq_case, Para_case) = caseData;
            outFile << "N+M = " << N_case + M_case << endl;
            outFile << "Sequential is faster for N=" << N_case << ", M=" << M_case << ", by " << Para_case - Seq_case << endl;
        }

    }
    //close and check output file closed
    outFile.close();
    assert(!outFile.is_open() && "Failed to close output.txt file");

    return 0;
}

void inLoopRun(int& N, int& M, double& seqTime, double& paraTime)
{
    seqTimer(N, M, seqTime);

    paraTimer(N, M, paraTime);

    /*if (seqTime < paraTime && (paraTime - seqTime) >= timeDiffThreshold)
    {
        sequentialFasterCases.emplace_back(N, M, seqTime, paraTime);
        seqWins++;
#ifdef DEBUG_ASSERT
        Assert(seqTime < paraTime, "Sequential is faster than parallelized");
#endif
    }

    totalRuns++;*/
}

void paraTimer(int N, int M, double& paraTime)
{
    //para testing
    double paraStart, paraEnd;
    paraStart = omp_get_wtime();
    double resultParallel = parrallelOptimizedRun(N, M);
    paraEnd = omp_get_wtime();
    paraTime = paraEnd - paraStart;
}

void seqTimer(int N, int M, double& seqTime)
{
    //sequential testing
    double seqStart, seqEnd;

    //get start time
    seqStart = omp_get_wtime();

    //run sequential function
    double result = sequentialRun(N, M);
    //get stop time
    seqEnd = omp_get_wtime();
    //get total time
    seqTime = seqEnd - seqStart;
}

double sequentialRun(int N, int M)
{
    long A = 0;
    double B = 0, C = 0, D = 0;

    for (long i = 0; i < N; i++)
    {
        for (long j = 0; j < M; j++)
        {
            A += i * j;
        }
    }

    for (long i = 1; i < (long)sqrt(A); i++)
    {
        B += 1 / i;
    }

    for (long i = 0; i < M * N; i++)
    {
        for (long j = 0; j < M; j++)
        {
            D += pow(0.1, i * j);
        }
    }

    for (long i = 0; i < (long)B * (N + 1); i++)
    {
        for (long j = 1; j < (long)sqrt(D); j++)
        {
            C += i / j;
        }
    }

    return A + B - C / D;
}

double parrallelOptimizedRun(int N, int M)
{

    long A = 0;
    double B = 0, C = 0, D = 0;

    // handle nested loop with data dependency A
#pragma omp parallel for if(N + M < 350) private(j) reduction(+:A) schedule(dynamic)
    for (long i = 0; i < N; i++)
    {
        for (long j = 0; j < M; j++)
        {
            A += i * j;
        }
    }

    // data dependency
    //if(N + M < 350)
#pragma omp parallel for reduction(+:B) schedule(dynamic)
    for (long i = 1; i < (long)sqrt(A); i++)
    {
        B += 1 / i;
    }

    //nested loops data dependency D and B
#pragma omp parallel for if(N + M < 350) private(j) reduction(+:D) schedule(dynamic)
    for (long i = 0; i < M * N; i++)
    {
        for (long j = 0; j < M; j++)
        {
            D += pow(0.1, i * j);
        }
    }

    //nested loops
    //if(N + M < 350)
#pragma omp parallel for reduction(+:C) schedule(dynamic)
    for (long i = 0; i < (long)B * (N + 1); i++)
    {
        for (long j = 1; j < (long)sqrt(D); j++)
        {
            C += i / j;
        }
    }

    //return result
    return A + B - C / D;
}



