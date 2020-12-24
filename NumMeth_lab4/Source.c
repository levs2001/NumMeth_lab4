#include<stdlib.h>
#include<stdio.h>
#include <math.h>
#include <time.h>
#include "Operations.h"
#pragma warning(disable:4996)

int size = 3;

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
	B[strM][colM] = 0;
	B[colM][strM] = 0;

	for (int m = 0; m < size; m++) {
		if (m != strM && m != colM) {
			B[strM][m] = c * B_prev[m][strM] + s * B_prev[m][colM];
			B[m][strM] = B[strM][m];
			B[colM][m] = -s * B_prev[m][strM] + c * B_prev[m][colM];
			B[m][colM] = B[colM][m];
		}
	}
}

int JacobiFindEig(const double** A, double** B, double** Tk) {
	int strM, colM;
	double maxEl, p, q, d, r, c, s;
	double** B_prev = malloc(sizeof(double*) * size);
	for (int i = 0; i < size; i++)
		B_prev[i] = malloc(sizeof(double) * size);
	CopyMatr(B_prev, A, ALLMATR);

	CopyMatr(B, B_prev, ALLMATR);
	maxEl = FindMaxElUnderMainDiag(B, &strM, &colM);
	p = 2 * maxEl;
	q = B[strM][strM] - B[colM][colM];
	FindCS(&c, &s, p, q);
	JacobiIterB(B, B_prev, c, s, strM, colM);

	for (int i = 0; i < size; i++) {
		free(B_prev[i]);
	}
	free(B_prev);
}

void main(void) {
	double* eigs = NULL;
	double** A = malloc(sizeof(double*) * size);
	for (int i = 0; i < size; i++) {
		A[i] = malloc(sizeof(double) * size);
	}
	eigs = CreateEigVect(1, 1);
	srand(time(0));
	CreateMatr(A, eigs);
	PrintMatrix(A, ALLMATR);
	for (int i = 0; i < size; i++) {
		free(A[i]);
	}
	free(A);
	free(eigs);
}