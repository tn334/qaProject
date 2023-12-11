// qaProject.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <iostream>
#include <fstream>
#include "qaProject.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <string>
#include <vector>
#include <tuple>
#include <fstream>

using namespace std;

int myMain()
{
    //open output file
    //ofstream outFile("../output.txt");
    //try for threads 2-5???

    //cout << "Opening output.txt file" << endl;
    //outFile.open("../output.txt");
    //assert(outFile.is_open() && "Failed to open output.txt file");

    //int threads;

    for (int threads = MINTHREADS; threads < MAXTHREADS; threads++)
    {
        
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

                if (seqTime < paraTime && (paraTime - seqTime) >= TIME_THRESHOLD)
                {
                    seqWins++;
                    sequentialFasterCases.emplace_back(N, M, seqTime, paraTime);
                }
                totalRuns++;
            }
        }

        loseRate = seqWins / totalRuns;
        
        //outFile << "\n\nPercentage of Sequential Wins for " << threads << " is " << loseRate << endl;
        
        // handle sequential wins
        for (const auto& caseData : sequentialFasterCases)
        {
            int N_case, M_case;
            double Seq_case, Para_case;
            tie(N_case, M_case, Seq_case, Para_case) = caseData;
            //outFile << "N+M = " << N_case + M_case << endl;
            //outFile << "Sequential is faster for N=" << N_case << ", M=" << M_case << ", by " << Para_case - Seq_case << endl;
        }
    }

    //close and check output file closed
    //outFile.close();
    //assert(!outFile.is_open() && "Failed to close output.txt file");

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

    firstSeqLoop(N, M, A);

    secondSeqLoop(A, B);

    thirdSeqLoop(M, N, D);

    fourthSeqLoop(B, N, D, C);

    return A + B - C / D;
}

void firstSeqLoop( int N, int M, long& A)
{
    for (long i = 0; i < N; i++)
    {
        for (long j = 0; j < M; j++)
        {
            A += i * j;
        }
    }
}

void secondSeqLoop(long A, double& B)
{
    for (long i = 1; i < (long)sqrt(A); i++)
    {
        B += 1 / i;
    } 
}

void thirdSeqLoop(int M, int N, double& D)
{
    for (long i = 0; i < M * N; i++)
    {
        for (long j = 0; j < M; j++)
        {
            D += pow(0.1, i * j);
        }
    }   
}

void fourthSeqLoop(double B, int N, double D, double& C)
{
    for (long i = 0; i < (long)B * (N + 1); i++)
    {
        for (long j = 1; j < (long)sqrt(D); j++)
        {
            C += i / j;
        }
    }    
}

double parrallelOptimizedRun(int N, int M)
{
    long A = 0;
    double B = 0, C = 0, D = 0;

    // handle nested loop with data dependency A
    firstParallelLoop(N, M, A);

    // data dependency
    secondParallelLoop(A, B);

    //nested loops data dependency D and B
    thirdParallelLoop(M, N, D);

    //nested loops
    fourthParallelLoop(B, N, D, C);

    //return result
    return A + B - C / D;
}

void firstParallelLoop(int N, int M, long& A)
{
    long j;
#pragma omp parallel for if(N + M < 350) private(j) reduction(+:A) schedule(dynamic)
    for (long i = 0; i < N; i++)
    {
        for ( j = 0; j < M; j++)
        {
            A += i * j;
        }
    }
}

void secondParallelLoop(long A, double& B)
{
    //if(N + M < 350)
#pragma omp parallel for reduction(+:B) schedule(dynamic)
    for (long i = 1; i < (long)sqrt(A); i++)
    {
        B += 1 / i;
    }
        
}

void thirdParallelLoop(int M, int N, double& D)
{
    long j;
#pragma omp parallel for if(N + M < 350) private(j) reduction(+:D) schedule(dynamic)
    for (long i = 0; i < M * N; i++)
    {
        for ( j = 0; j < M; j++)
        {
            D += pow(0.1, i * j);
        }
    }
        
}

void fourthParallelLoop(double B, int N, double D, double& C)
{
    //if(N + M < 350)
#pragma omp parallel for reduction(+:C) schedule(dynamic)
    for (long i = 0; i < (long)B * (N + 1); i++)
    {
        for (long j = 1; j < (long)sqrt(D); j++)
        {
            C += i / j;
        }
    }
        
}
