#include "quanc8.hpp"
#include <cmath>

void quanc8(double (*FUN)(double), double A, double B, double ABSERR,
            double RELERR, double *RESULT, double *ERREST, int *NOFUN, double *FLAG)

{
    double QRIGHT[32], F[17], X[17], FSAVE[9][31], XSAVE[9][31];
    int LEVMIN, LEVMAX, LEVOUT, NOMAX, NOFIN, LEV, NIM, J, I;
    double W0, W1, W2, W3, W4, COR11, AREA, X0, F0, STONE, STEP;
    double QLEFT, QNOW, QDIFF, QPREV, TOLERR, ESTERR;


    LEVMIN = 1;
    LEVMAX = 30;
    LEVOUT = 6;
    NOMAX = 5000;
    NOFIN = NOMAX - (8 * (LEVMAX - LEVOUT + pow(2, (double) (LEVOUT + 1))));

    W0 = 3956.0 / 14175.0;
    W1 = 23552.0 / 14175.0;
    W2 = -3712.0 / 14175.0;
    W3 = 41984.0 / 14175.0;
    W4 = -18160.0 / 14175.0;

    *FLAG = 0.0;
    *RESULT = 0.0;
    COR11 = 0.0;
    *ERREST = 0.0;
    AREA = 0.0;
    *NOFUN = 0;
    if (A == B)
        return;

    LEV = 0;
    NIM = 1;
    X0 = A;
    X[16] = B;
    QPREV = 0.0;
    F0 = FUN(X0);
    STONE = (B - A) / 16.0;
    X[8] = (X0 + X[16]) / 2.0;
    X[4] = (X0 + X[8]) / 2.0;
    X[12] = (X[8] + X[16]) / 2.0;
    X[2] = (X0 + X[4]) / 2.0;
    X[6] = (X[4] + X[8]) / 2.0;
    X[10] = (X[8] + X[12]) / 2.0;
    X[14] = (X[12] + X[16]) / 2.0;
    for (J = 2; J <= 16; J = J + 2) {
        F[J] = FUN(X[J]);
    }
    *NOFUN = 9;

    trenta:
    X[1] = (X0 + X[2]) / 2.0;
    F[1] = FUN(X[1]);
    for (J = 3; J <= 15; J = J + 2) {
        X[J] = (X[J - 1] + X[J + 1]) / 2.0;
        F[J] = FUN(X[J]);
    }
    *NOFUN = *NOFUN + 8;
    STEP = (X[16] - X0) / 16.0;
    QLEFT = (W0 * (F0 + F[8]) + W1 * (F[1] + F[7]) + W2 * (F[2] + F[6]) + W3 * (F[3] + F[5])
             + W4 * F[4]) * STEP;
    QRIGHT[LEV + 1] = (W0 * (F[8] + F[16]) + W1 * (F[9] + F[15]) + W2 * (F[10] + F[14])
                       + W3 * (F[11] + F[13]) + W4 * F[12]) * STEP;
    QNOW = QLEFT + QRIGHT[LEV + 1];
    QDIFF = QNOW - QPREV;
    AREA = AREA + QDIFF;

    ESTERR = fabs(QDIFF) / 1023.0;
    if (ABSERR > (RELERR * fabs(AREA)) * (STEP / STONE))
        TOLERR = ABSERR;
    else
        TOLERR = (RELERR * fabs(AREA)) * (STEP / STONE);
    if (LEV < LEVMIN)
        goto cinquanta;
    if (LEV >= LEVMAX)
        goto sessantadue;
    if (*NOFUN > NOFIN)
        goto sessanta;
    if (ESTERR <= TOLERR)
        goto settanta;

    cinquanta:
    NIM = 2 * NIM;
    LEV = LEV + 1;

    for (I = 1; I <= 8; I++) {
        FSAVE[I][LEV] = F[I + 8];
        XSAVE[I][LEV] = X[I + 8];
    }

    QPREV = QLEFT;
    for (I = 1; I <= 8; I++) {
        J = -I;
        F[2 * J + 18] = F[J + 9];
        X[2 * J + 18] = X[J + 9];
    }
    goto trenta;

    sessanta:
    NOFIN = 2 * NOFIN;
    LEVMAX = LEVOUT;
    *FLAG = *FLAG + ((B - X0) / (B - A));
    goto settanta;

    sessantadue:
    *FLAG = *FLAG + 1.0;

    settanta:
    *RESULT = *RESULT + QNOW;
    *ERREST = *ERREST + ESTERR;
    COR11 = COR11 + QDIFF / 1023.0;

    while (NIM % 2 != 0) {
        NIM = NIM / 2;
        LEV = LEV - 1;
    }
    NIM = NIM + 1;
    if (LEV <= 0)
        goto ottanta;

    QPREV = QRIGHT[LEV];
    X0 = X[16];
    F0 = F[16];
    for (I = 1; I <= 8; I++) {
        F[2 * I] = FSAVE[I][LEV];
        X[2 * I] = XSAVE[I][LEV];
    }
    goto trenta;

    ottanta:
    *RESULT = *RESULT + COR11;
    if (*ERREST == 0.0)
        return;
    while (fabs(*RESULT) + (*ERREST) == fabs(*RESULT))
        *ERREST = 2.0 * (*ERREST);
    return;
}
