// qaProject.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <string>
#include <vector>
#include <tuple>
#include <fstream> //ofstream, instream, .eof(), .get()

using namespace std;

//declare consts
const int parameterStart = 50;
const int parameterEnd = 300;
const int incrementor = 50;
const double timeDiffThreshold = .0005;

// functions
double sequentialRun(int N, int M);
void fourthSeqLoop(long& i, double B, int N, long& j, double D, double& C);
void thirdSeqLoop(long& i, int M, int N, long& j, double& D);
void secondSeqLoop(long& i, long A, double& B);
void firstSeqLoop(long& i, int N, long& j, int M, long& A);
double parrallelOptimizedRun(int N, int M);

void firstParallelLoop(long& i, int N, long& j, int M, long& A);
void secondParallelLoop(long& i, long A, double& B);
void thirdParallelLoop(long& i, int M, int N, long& j, double& D);
void fourthParallelLoop(long& i, double B, int N, long& j, double D, double& C);

int main()
{
    int threads = 2, maxThreads = 5;
    //open output file
    ofstream outFile;
    //try for threads 2-5???
    for (threads = 2; threads <= maxThreads; threads++)
    {
        // you may add some code here to measure the execution time
        double seqTime, paraTime;
        double seqWins = 0, totalRuns = 0;
        double loseRate = 0;
        vector<tuple<int, int, double, double>> sequentialFasterCases;

        omp_set_num_threads(threads);
        //sequential
        for (int N = parameterStart; N <= parameterEnd; N += incrementor)
        {
            for (int M = parameterStart; M <= parameterEnd; M += incrementor)
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

                //para testing
                double paraStart, paraEnd;
                paraStart = omp_get_wtime();
                double resultParallel = parrallelOptimizedRun(N, M);
                paraEnd = omp_get_wtime();
                paraTime = paraEnd - paraStart;
                
                if (seqTime < paraTime && (paraTime - seqTime) >= timeDiffThreshold)
                {
                    sequentialFasterCases.emplace_back(N, M, seqTime, paraTime);
                    seqWins++;
                }

                totalRuns++;
            }
        }
        loseRate = seqWins / totalRuns;
        // open output file
        outFile.open("/output.txt");
        outFile << "\n\nPercentage of Sequential Wins for" << threads << "is" << loseRate << endl;
        
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
    outFile.close();
    return 0;
}

double sequentialRun(int N, int M)
{
    long i, j;
    long A = 0;
    double B = 0, C = 0, D = 0;

    firstSeqLoop(i, N, j, M, A);

    secondSeqLoop(i, A, B);

    thirdSeqLoop(i, M, N, j, D);

    fourthSeqLoop(i, B, N, j, D, C);

    return A + B - C / D;
}

void firstSeqLoop(long& i, int N, long& j, int M, long& A)
{
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            A += i * j;
        }
    }
}

void secondSeqLoop(long& i, long A, double& B)
{
    for (i = 1; i < (long)sqrt(A); i++)
    {
        B += 1 / i;
    }
        
}

void thirdSeqLoop(long& i, int M, int N, long& j, double& D)
{
    for (i = 0; i < M * N; i++)
    {
        for (j = 0; j < M; j++)
        {
            D += pow(0.1, i * j);
        }
    }
        
}

void fourthSeqLoop(long& i, double B, int N, long& j, double D, double& C)
{
    for (i = 0; i < (long)B * (N + 1); i++)
    {
        for (j = 1; j < (long)sqrt(D); j++)
        {
            C += i / j;
        }
    }
        
}

double parrallelOptimizedRun(int N, int M)
{
    long i, j;
    long A = 0;
    double B = 0, C = 0, D = 0;

//#pragma omp parallel for if(N + M < 350) private(j) reduction(+:A) schedule(dynamic)
//    for (i = 0; i < N; i++)
//    {
//        for (j = 0; j < M; j++)
//        {
//            A += i * j;
//        }
//    }
//
//#pragma omp parallel for if(N + M < 350) reduction(+:B) schedule(dynamic)
//    for (i = 1; i < (long)sqrt(A); i++)
//        B += 1 / i;
//
//#pragma omp parallel for if(N + M < 350) private(j) reduction(+:D) schedule(dynamic)
//    for (i = 0; i < M * N; i++)
//        for (j = 0; j < M; j++)
//        {
//            D += pow(0.1, i * j);
//        }
//
//#pragma omp parallel for if(N + M < 350) reduction(+:C) schedule(dynamic)
//    for (i = 0; i < (long)B * (N + 1); i++)
//        for (j = 1; j < (long)sqrt(D); j++)
//        {
//            C += i / j;
//        }
    // handle nested loop with data dependency A
        //fast paralellization
    //if(N * M > 100)
    firstParallelLoop(i, N, j, M, A);

    // data dependency
    //if(N > 100)
    secondParallelLoop(i, A, B);

    //nested loops data dependency D and B
    //if(N * M > 100)
    thirdParallelLoop(i, M, N, j, D);

    //nested loops
    //if((long)B * (N + 1) > 100)
    fourthParallelLoop(i, B, N, j, D, C);

    //return result
    return A + B - C / D;
}

void firstParallelLoop(long& i, int N, long& j, int M, long& A)
{
#pragma omp parallel for if(N + M < 350) private(j) reduction(+:A) schedule(dynamic)
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            A += i * j;
        }
    }
}

void secondParallelLoop(long& i, long A, double& B)
{
    //if(N + M < 350)
#pragma omp parallel for reduction(+:B) schedule(dynamic)
    for (i = 1; i < (long)sqrt(A); i++)
    {
        B += 1 / i;
    }
        
}

void thirdParallelLoop(long& i, int M, int N, long& j, double& D)
{
#pragma omp parallel for if(N + M < 350) private(j) reduction(+:D) schedule(dynamic)
    for (i = 0; i < M * N; i++)
    {
        for (j = 0; j < M; j++)
        {
            D += pow(0.1, i * j);
        }
    }
        
}

void fourthParallelLoop(long& i, double B, int N, long& j, double D, double& C)
{
    //if(N + M < 350)
#pragma omp parallel for reduction(+:C) schedule(dynamic)
    for (i = 0; i < (long)B * (N + 1); i++)
    {
        for (j = 1; j < (long)sqrt(D); j++)
        {
            C += i / j;
        }
    }
        
}
