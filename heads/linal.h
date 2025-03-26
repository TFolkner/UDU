#include <iostream>
#include <math.h>
#include <fstream>
#include <stdio.h>


// gauss rand
double Gauss_Rnd(double m, double sigma) {
	double r1, r2, z;

	r1 = double(rand()) / double(RAND_MAX);
	r2 = double(rand()) / double(RAND_MAX);

	if (r1 < 1e-10)	r1 = 1e-10;

	z = sqrt(-2 * log(r1)) * cos(2 * 3.14*r2);
	z = sigma * z + m;
	return z;
}

// zeros matrix
void ZerosMatrix(double* matrix, int matrix_len) {
	for (int i = 0; i < matrix_len; i++)
		*(matrix + i) = 0;
	return;
}

// line_matr[M][N] * vector[N]
void line_matr_mull_vec(double* matr, double* vec, double* res, const int M, const int N) {
	int i, j;
	for (i = 0; i < M; i++)
		for (j = 0; j < N; j++)
			*(res+i) += *(matr + N*i + j) * *(vec+j);

	return;
}

// line_matr[M][N] * line_matr[N][L] = res[M][L]
void line_matr_mull_line_matr(double* matr_1, double* matr_2, double* res,
	const int M, const int N, const int L) {

	int i, j, k;
	for (i = 0; i < M; i++)
		for (j = 0; j < L; j++)
			for (k = 0; k < N; k++)
				*(res+L*i + j) += *(matr_1+N*i + k) * *(matr_2+L*k + j);

	return;
}

// line_matr[M][N] + line_matr[M][N] = res[M][N]
void line_matr_plus_line_matr(double* matr_1, double* matr_2, double* res, const int M, const int N) {
	for (int i = 0; i < M*N; i++)
		*(res + i) = *(matr_1+i) + *(matr_1+i);
	return;
}

// line_matr[M][N] - line_matr[M][N] = res[M][N]
void line_matr_minus_line_matr(double* matr_1, double* matr_2, double* res, const int M, const int N) {
	for (int i = 0; i < M*N; i++)
		*(res + i) = *(matr_1 + i) - *(matr_1 + i);
	return;
}

// transpose matrix
void transpose_matrix(double* matr, double* res, const int M, const int N) {
	int line, column;
	for (line = 0; line < M; line++)
		for (column = 0; column < N; column++)
			*(res+N * line + column) = *(matr+N * column + line);

	return;
}

// print matrix[M][N] -> matrix save as line
void prt_matr_line(const int M, const int N, double* array) {
	int i, j;
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++)
			printf("%-10.2f ", *(array+N*i + j));
		printf("\n");
	}
	printf("\n");
}