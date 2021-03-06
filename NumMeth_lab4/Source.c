#include<stdlib.h>
#include<stdio.h>
#include <math.h>
#include <time.h>
#include "Operations.h"
#pragma warning(disable:4996)

int size = 10;

void FindCS(double* c, double* s, double p, double q) {
	double r;
	double d = sqrt(pow(p, 2) + pow(q, 2));

	if (q != 0) {
		r = fabs(q) / (2 * d);
		*c = sqrt(0.5 + r);
		*s = sqrt(0.5 - r) * sign(p * q);
	}
	else {
		*c = sqrt(2) / 2;
		*s = *c;
	}
}

void JacobiIterB(double** B, const double** B_prev, double c, double s, int strM, int colM) {
	B[strM][strM] = pow(c, 2) * B_prev[strM][strM] + pow(s, 2) * B_prev[colM][colM] + 2 * c * s * B_prev[strM][colM];
	B[colM][colM] = pow(s, 2) * B_prev[strM][strM] + pow(c, 2) * B_prev[colM][colM] - 2 * c * s * B_prev[strM][colM];
	B[strM][colM] = (pow(c, 2) - pow(s, 2)) * B_prev[strM][colM] + c * s * (B_prev[colM][colM] - B_prev[strM][strM]);
	B[colM][strM] = B[strM][colM];

	for (int m = 0; m < size; m++) {
		if (m != strM && m != colM) {
			B[strM][m] = c * B_prev[m][strM] + s * B_prev[m][colM];
			B[m][strM] = B[strM][m];
			B[colM][m] = -s * B_prev[m][strM] + c * B_prev[m][colM];
			B[m][colM] = B[colM][m];
		}
	}
}

void JacobiIterT(double** Tk, double c, double s, int strM, int colM) {
	double** Tcur = malloc(sizeof(double*) * size);
	double** TkCopy = malloc(sizeof(double*) * size);
	for (int i = 0; i < size; i++) {
		Tcur[i] = malloc(sizeof(double) * size);
		TkCopy[i] = malloc(sizeof(double) * size);
	}
	CopyMatr(TkCopy, Tk, ALLMATR);
	CopyEMatr(Tcur);
	Tcur[strM][strM] = c;
	Tcur[colM][colM] = c;
	Tcur[colM][strM] = s;
	Tcur[strM][colM] = -s;
	CompMatr(Tk, TkCopy, Tcur, size, size, 0);
	for (int i = 0; i < size; i++) {
		free(Tcur[i]);
		free(TkCopy[i]);
	}
	free(Tcur);
	free(TkCopy);
}

int JacobiFindEig(const double** A, double** B, double** Tk, double eps) {
	//B - matrix of eigenvalues (name B as in book)
	//Tk - matrix of own vectors (name Tk as in book)
	int iterCount = 0;
	int strM, colM;
	double maxEl, p, q, d, r, c, s;
	double** B_prev = malloc(sizeof(double*) * size);
	for (int i = 0; i < size; i++) {
		B_prev[i] = malloc(sizeof(double) * size);
	}

	CopyEMatr(Tk);
	CopyMatr(B, A, ALLMATR);

	do {
		CopyMatr(B_prev, B, ALLMATR);
		maxEl = FindMax(B, &strM, &colM); // FindMaxElUnderMainDiag(B, &strM, &colM);
		p = 2 * maxEl;
		q = B[strM][strM] - B[colM][colM];
		FindCS(&c, &s, p, q);
		JacobiIterB(B, B_prev, c, s, strM, colM);
		JacobiIterT(Tk, c, s, strM, colM);
		iterCount++;
	} while (fabs(maxEl) >= eps);

	for (int i = 0; i < size; i++) {
		free(B_prev[i]);
	}
	free(B_prev);
	return iterCount;
}

void PrintIterCountForStepEigs(double otd, FILE* file, double** A, double** B, double** Tk) {
	int iterCount = 0;
	double eps = 0;
	double* eigs = CreateEigVect(1000, 1000, otd);
	for (int i = 0; i < size; i++)
		printf("%f\n", eigs[i]);
	CreateMatr(A, eigs);
	for (eps = 1e-5; eps < 10; eps *= 10) {
		iterCount = JacobiFindEig(A, B, Tk, eps);
		printf("%d\t", iterCount);
		fprintf(file, "%d\t", iterCount);
	}
	printf("\n");
	fprintf(file, "\n");
	free(eigs);
}

void main(void) {
	int iterCount = 0;
	double eps = 1e-4;
	double* eigs = NULL;
	double** A = malloc(sizeof(double*) * size);
	double** B = malloc(sizeof(double*) * size);
	double** Tk = malloc(sizeof(double*) * size);
	for (int i = 0; i < size; i++) {
		A[i] = malloc(sizeof(double) * size);
		B[i] = malloc(sizeof(double) * size);
		Tk[i] = malloc(sizeof(double) * size);
	}


	srand(time(0));
	//A[0][0] = 1;
	//A[0][1] = 2;
	//A[1][0] = 3;
	//A[1][1] = 4;
	//eigs = CreateEigVect(1000, 1000, 10);
	//for (int i = 0; i < size; i++) {
	//	printf("%f\t", eigs[i]);
	//}
	//printf("\n");
	//CreateMatr(A, eigs);
	//PrintMatrix(A, ALLMATR);
	//printf("\n");
	//iterCount = JacobiFindEig(A, B, Tk, eps);
	//printf("\nEigs:\n");
	//PrintMatrix(B, ALLMATR);
	//printf("\nOwn vectors:\n");
	//PrintMatrix(Tk, ALLMATR);

	FILE* file = fopen("D:\\MyProgramms\\Num_Meth_Lab4\\IterCount.txt", "w");
	for (double otd = 1e-4; otd < 1e2; otd *= 10) {
		PrintIterCountForStepEigs(otd, file, A, B, Tk);
	}
	
	for (int i = 0; i < size; i++) {
		free(A[i]);
	}
	free(A);
	free(eigs);
}