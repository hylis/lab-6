#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "winbgi2.h"
#include "gauss.h"

int n = 21;
int L = 58;
double h = 1.0/ n ;
double Q(int x) { return pow(10, 4) * sin(x * 3.14); }
void HilberMatrix(int N, double** H);
void displayMatrix(int N, double** H);
void computeVec(int N, double** H, double* b);
void plotVec(int N, double* v);
void computeMatrix(int N, double** K);
void computeVector(int N, double* F);
void displayVector(int n, double* F);
void main() {
	FILE* f;
	f = fopen("dane.txt", "w");
int N = 6;
double** H = (double**)malloc(N * sizeof(double*));
for (int i = 0;i < N;i++) {
	H[i] = (double*)malloc(N * sizeof(double));
}
double* b = (double*)malloc(N * sizeof(double));
double* x = (double*)malloc(N * sizeof(double));

HilberMatrix(N, H);
displayMatrix(N, H);
computeVec(N, H, b);
printf("\n");
plotVec(N, b);
gauss(N, H, x, b);
printf("\n");
plotVec(N, x);
double** K;
double* F, * T;
K = (double**)malloc(n * sizeof(double*));
for (int i = 0; i < n; i++) {
	K[i] = (double*)malloc(n * sizeof(double));
}
F = (double*)malloc(n * sizeof(double));
T = (double*)malloc(n * sizeof(double));
printf("\n");
computeMatrix(n, K);
displayMatrix(n, K);
computeVector(n, F);
plotVec(n, F);
gauss(n, K, F, T);
plotVec(n, T);
for (int i = 0;i < n;i++) {
	fprintf(f, "%lf\n", T[i]);
}
for (int i = 0;i < N;i++) {
	free(H[i]);
}
free(H);
free(b);
free(x);
for (int i = 0;i < N;i++) {
	free(K[i]);
}
free(K);
free(F);
free(T);
fclose(f);
}
void HilberMatrix(int N, double** H) {
	for (int i = 0;i < N;i++) {
		for (int j = 0;j < N;j++) {
			H[i][j] = 1.0 / (1.0 + i + j);
		}
	}
}
void displayMatrix(int N, double** H) {
	for (int i = 0;i < N;i++) {
		for (int j = 0;j < N;j++) {
			printf("%lf\t", H[i][j]);
		}
		printf("\n");
	}
}
void computeVec(int N, double** H, double* b) {
	for (int i = 0;i < N;i++) {
		for (int j = 0;j < N;j++) {
			b[i]=b[i] + H[i][j];
		}
	}
}
void plotVec(int N, double* v) {
	for (int i = 0;i < N ;i++) {
		printf("%lf\n", v[i]);
	}
}
void computeMatrix(int N, double** K) {
	K[0][0] = 1;
	K[N][N] = 1;
	for(int i = 1;i < N;i++) {
		K[i][i - 1] = 1;
		K[i][i] = -2;
		K[i][i + 1] = 1;
	}
}
void computeVector(int N, double* F) {
	F[0] = 273;
	F[N] = 300;
	for (int i = 1;i < N - 1;i++) {
		F[i] = Q(h * i) / L * pow(h, 2);
	}
}