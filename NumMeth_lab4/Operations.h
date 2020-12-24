#pragma once
#define ALLMATR -1
#define ISVECT -2

void CopyMatr(double** in, double** from, int stolbNum);
double Find2VectNorma(double** matrix, int stolb, int isVect);
double Find1MatrNorma(double** matrix);
double Find1VectNorma(double* vect);
void ShareVectByNorm(double** matrix, int stolb);
void MinusVect(double** minuend, double** subtrahend, int stolbMin, int stolbSub, int isVect);
void PrintMatrix(double** matrix, int stolbNum);
void TranspMatr(double** matr);
void CompMatr(double** result, double** matr1, double** matr2, int matr1Stolb_c, int matr2Str_c, int isVect);
double getRand();
double getRandRange(double min, double max);
void CreateDiagMatrix(double** diagMatrix, double* eigs);
void Q_MatrixGeneration(double** Matrix_Q, double* Vector_W);
void CreateMatr(double** A, double* eigs);
double* CreateEigVect(double minEig, double step);
double FindMaxElUnderMainDiag(double** A, int* _string, int* _column);
int sign(double val);