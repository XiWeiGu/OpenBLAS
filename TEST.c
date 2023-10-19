#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"

int main() {
    int m, n, k;

    // 输入矩阵规模
    printf("请输入矩阵的规模 (m n k): ");
    scanf("%d %d %d", &m, &n, &k);

    if (m <= 0 || n <= 0 || k <= 0) {
        printf("矩阵规模必须为正整数。\n");
        return 1;
    }

    // 分配内存并初始化矩阵 A, B, 和 C
    double *A = (double *)malloc(m * k * sizeof(double));
    double *B = (double *)malloc(k * n * sizeof(double));
    double *C = (double *)malloc(m * n * sizeof(double));

    if (A == NULL || B == NULL || C == NULL) {
        printf("内存分配失败。\n");
        return 1;
    }

    // 填充矩阵 A 和 B
    for (int i = 0; i < m * k; i++) {
        A[i] = (double)rand() / RAND_MAX; // 随机数填充 A
    }

    for (int i = 0; i < k * n; i++) {
        B[i] = (double)rand() / RAND_MAX; // 随机数填充 B
    }

    // 调用OpenBLAS的dgemm接口进行矩阵相乘
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, A, k, B, n, 0.0, C, n);

    // 释放内存
    free(A);
    free(B);
    free(C);
}
