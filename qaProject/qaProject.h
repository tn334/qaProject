// qaProject.h

#ifndef QAPROJECT_H
#define QAPROJECT_H

const int parameterStart = 50;
const int parameterEnd = 300;
const int incrementor = 50;
const double timeDiffThreshold = .0005;

void firstSeqLoop(int N, int M, long& A);
void thirdSeqLoop(int M, int N, double& D);
void fourthSeqLoop(double B, int N, double D, double& C);
void secondSeqLoop(long A, double& B);
double sequentialRun(int N, int M);

void firstParallelLoop(int N, int M, long& A);
void thirdParallelLoop(int M, int N, double& D);
void fourthParallelLoop(double B, int N, double D, double& C);
void secondParallelLoop(long A, double& B);
double parrallelOptimizedRun(int N, int M);

void seqTimer(int N, int M, double& seqTime);
void paraTimer(int N, int M, double& paraTime);

#endif // QAPROJECT_H