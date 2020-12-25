#include<stdlib.h>
#include<stdio.h>
#include <math.h>
#include "Operations.h"

#pragma warning(disable:4996)

extern int size;

void CopyMatr(double** in, double** from, int stolbNum) {
	if (stolbNum == ALLMATR) {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				in[i][j] = from[i][j];
			}
		}
	}
	else {
		for (int i = 0; i < size; i++) {
			in[i][stolbNum] = from[i][stolbNum];
		}
	}
}

double Find1MatrNorma(double** matrix) {
	double maxNorma = 0, norma = 0;
	for (int j = 0; j < size;j++) {
		for (int i = 0; i < size; i++) {
			norma += fabs(matrix[i][j]);
		}
		if (norma > maxNorma) {
			maxNorma = norma;
		}
		norma = 0;
	}
	return maxNorma;
}

double Find1VectNorma(double* vect) {
	double norma = 0;
	for (int i = 0; i < size; i++) {
		norma += fabs(vect[i]);
	}
	return norma;
}

double Find2VectNorma(double** matrix, int stolb, int isVect) {
	double norma = 0;
	if (isVect) {
		for (int i = 0; i < size; i++) {
			norma += pow((*matrix)[i], 2);
		}
	}
	else {
		for (int i = 0; i < size; i++) {
			norma += pow(matrix[i][stolb], 2);
		}
	}
	norma = sqrt(norma);
	return norma;
}

void ShareVectByNorm(double** matrix, int stolb) {
	double norma = Find2VectNorma(matrix, stolb, 0);
	for (int i = 0; i < size; i++) {
		matrix[i][stolb] /= norma;
	}
}

void MinusVect(double** minuend, double** subtrahend, int stolbMin, int stolbSub, int isVect) {
	if (isVect) {
		for (int i = 0; i < size; i++) {
			(*minuend)[i] -= (*subtrahend)[i];
		}
	}
	else {
		for (int i = 0; i < size; i++) {
			if (stolbSub == ISVECT)
				minuend[i][stolbMin] -= (*subtrahend)[i];
			else
				minuend[i][stolbMin] -= subtrahend[i][stolbSub];
		}
	}
}

void PrintMatrix(double** matrix, int stolbNum) {
	if (stolbNum == ALLMATR) {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				printf("%.15lf\t", matrix[i][j]);
			}
			printf("\n");
		}
	}
	else if (stolbNum == ISVECT) {
		for (int i = 0; i < size; i++) {
			printf("%.15lf\n", (*matrix)[i]);
		}
	}
	else {
		for (int i = 0; i < size; i++) {
			printf("%.15lf\n", matrix[i][stolbNum]);
		}
	}
	printf("\n");

}

void TranspMatr(double** matr) {
	double tmp;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < i;j++) {
			tmp = matr[i][j];
			matr[i][j] = matr[j][i];
			matr[j][i] = tmp;
		}
	}
}

void CompMatr(double** result, double** matr1, double** matr2, int matr1Stolb_c, int matr2Str_c, int isVect) {
	if (isVect) {
		for (int i = 0; i < matr2Str_c; i++) {
			for (int j = 0; j < matr1Stolb_c; j++) {
				(*result)[i] = 0;
				for (int k = 0; k < matr2Str_c; k++) {
					(*result)[i] += matr1[i][k] * (*matr2)[k];
				}
			}
		}
	}
	else {
		for (int i = 0; i < matr2Str_c; i++) {
			for (int j = 0; j < matr1Stolb_c; j++) {
				result[i][j] = 0;
				for (int k = 0; k < matr2Str_c; k++) {
					result[i][j] += matr1[i][k] * matr2[k][j];
				}
			}
		}
	}

}

double getRand()
{
	return (double)rand() / (double)RAND_MAX;
}

double getRandRange(double min, double max)
{
	return getRand() * (max - min) + min;
}

void CreateDiagMatrix(double** diagMatrix, double* eigs)
{
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			diagMatrix[i][j] = 0;
		}
	}
	for (int i = 0; i < size; i++)
	{
		diagMatrix[i][i] = eigs[i];
	}

}

void Q_MatrixGeneration(double** Matrix_Q, double* Vector_W)
{
	//
	//Q = E - 2 * V*V^T / ||V||^2
	//

	double VectorNormSquared = 0;
	for (int i = 0; i < size; i++)
		VectorNormSquared += Vector_W[i] * Vector_W[i];

	for (int i = 0; i < size; i++)
		for (int j = 0; j <= i; j++)
		{
			Matrix_Q[i][j] = -2 * Vector_W[i] * Vector_W[j] / VectorNormSquared;
			Matrix_Q[j][i] = Matrix_Q[i][j];
		}

	for (int i = 0; i < size; i++)
		Matrix_Q[i][i] += 1;
}

void CreateMatr(double** A, double* eigs) {
	double** D = malloc(sizeof(double) * size), ** Q = malloc(sizeof(double) * size), ** tempM = malloc(sizeof(double) * size), * vect_W = malloc(sizeof(double) * size);
	for (int i = 0; i < size; i++) {
		D[i] = malloc(sizeof(double) * size);
		Q[i] = malloc(sizeof(double) * size);
		tempM[i] = malloc(sizeof(double) * size);
	}

	//Задаём рандомный вектор W
	for (int i = 0; i < size; i++) {
		vect_W[i] = getRandRange(1, size);
	}
	CreateDiagMatrix(D, eigs);
	Q_MatrixGeneration(Q, vect_W);

	CompMatr(tempM, Q, D, size, size, 0);
	CompMatr(A, tempM, Q, size, size, 0);


	for (int i = 0; i < size; i++) {
		free(D[i]);
		free(Q[i]);
		free(tempM[i]);
	}
	free(D);
	free(Q);
	free(tempM);
}

double* CreateEigVect(double minEig, double step) {
	double* eigs = malloc(sizeof(double) * size);
	for (int i = 0; i < size; i++) {
		eigs[i] = minEig + i * step;
	}
	return eigs;
}

double FindMaxElUnderMainDiag(double** matr, int* strM, int* colM) {
	*strM = 1;
	*colM = 0;
	double max = matr[*strM][*colM];
	
	for (int i = 1; i < size; i++) {
		for (int j = 0; j < i; j++) {
			if (matr[i][j] > max) {
				max = matr[i][j];
				*strM = i;
				*colM = j;
			}
		}
	}
	return max;
}

int sign(double val) {
	return val > 0 ? 1 : val == 0 ? 0 : -1;
}

void CopyEMatr(double** matr) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (i == j)
				matr[i][j] = 1;
			else
				matr[i][j] = 0;
		}
	}
}