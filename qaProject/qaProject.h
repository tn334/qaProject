// qaProject.h

#ifndef QAPROJECT_H
#define QAPROJECT_H

const int MAXTHREADS = 5;
const int MINTHREADS = 2;
const int PARAMTER_START = 50;
const int PARAMETER_END = 300;
const int INCREMENTOR = 50;
const double TIME_THRESHOLD = .005;
//const int 

int myMain();

void inLoopRun(int& N, int& M, double& seqTime, double& paraTime);
void firstSeqLoop(int N, int M, long& A);
void thirdSeqLoop(int M, int N, double& D);
void fourthSeqLoop(double B, int N, double D, double& C);
void secondSeqLoop(long A, double& B);
double sequentialRun(int N, int M);

void firstParallelLoop(int N, int M, long& A);
void thirdParallelLoop(int M, int N, double& D);
void fourthParallelLoop(double B, int N, double D, double& C);
void secondParallelLoop(long A, double& B);
double parrallelOptimizedRun(int N, int M, int threads);

void seqTimer(int N, int M, double& seqTime);
void paraTimer(int N, int M, double& paraTime);

#endif // QAPROJECT_H