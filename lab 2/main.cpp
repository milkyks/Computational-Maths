#include <iostream>
#include <iomanip>
#include "inverse-matrix.hpp"


int main() {
    const int N = 8;
    double norma = 0, max = 0;
    double p_values[] = {1.0, 0.1, 0.01, 0.0001, 0.000001};
    double inverted[N][N];
    double comp[N][N];
    double R[N][N];
    double work[N];
    double cond = 0;
    int ipvt[N];
    double **lu = new double *[N];
    double E[N][N] =
            {{1, 0, 0, 0, 0, 0, 0, 0},
             {0, 1, 0, 0, 0, 0, 0, 0},
             {0, 0, 1, 0, 0, 0, 0, 0},
             {0, 0, 0, 1, 0, 0, 0, 0},
             {0, 0, 0, 0, 1, 0, 0, 0},
             {0, 0, 0, 0, 0, 1, 0, 0},
             {0, 0, 0, 0, 0, 0, 1, 0},
             {0, 0, 0, 0, 0, 0, 0, 1}};

    for (int k = 0; k < 5; k++) {
        std::cout << "p = " << p_values[k] << "\n";

        double matrix[N][N] =
                {{p_values[k] + 12, 3,  2,  0,  -3, -8,  0,  0},
                 {-6,               26, 0,  -3, 8,  7,   -7, 7},
                 {-5,               -2, -3, -4, -6, 1,   0,  7},
                 {-6,               -5, -4, -2, 8,  7,   4,  -8},
                 {-2,               -3, 1,  -2, 10, -8,  6,  0},
                 {7,                -7, 2,  2,  0,  -24, -6, 4},
                 {-3,               -7, 0,  0,  6,  1,   13, 2},
                 {-6,               -3, 4,  0,  2,  -8,  6,  11}};

        for (int i = 0; i < N; ++i) {
            lu[i] = new double[N];

            for (int j = 0; j < N; ++j) {
                lu[i][j] = matrix[i][j];
            }
        }

        decomp(N, lu, &cond, ipvt, work);

        double vector[N];

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                vector[j] = 0;
            }
            vector[i] = 1.0;

            solve(N, lu, vector, ipvt);

            for (int j = 0; j < N; j++) {
                inverted[j][i] = vector[j];
            }
            vector[i] = 0;
        }

        std::cout << "Matrix A^(-1):" << std::endl;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                std::cout << std::left << std::setw(9) << inverted[i][j] << " ";
            }
            std::cout << "\n";
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                comp[i][j] = 0;

                for (int p = 0; p < N; p++) {
                    comp[i][j] += matrix[i][p] * inverted[p][j];
                }
            }
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                R[i][j] = E[i][j] - comp[i][j];
                max += fabs(R[i][j]);
            }

            if (max > norma) {
                norma = max;
            };
        }

        std::cout << "Cond - " << cond << "\n"
                  << "||R|| = ||E - A * A^(-1)|| = " << norma << "\n\n";

        norma = 0;
    }

    return 0;
}

