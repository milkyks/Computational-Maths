#include <iostream>
#include <iomanip>
#include <cmath>
#include "rkf45.hpp"
#include "zeroin.hpp"
#include "quanc8.hpp"

const int N = 2;
const double START = 0.0;
const double END = 10.0;
const double STEP = 0.5;
const double E = 0.5 * 1.01;
const double A = 1.0;
double sigma;
double B;

double zeroin_function(double x)
{
    return (pow(1.3, x) - x);
}

double function(double x)
{
    return (sin(1.1 * x)) / x;
}

double calculate_sigma()
{
    double ax = 7.856;
    double bx = 7.858;
    double tol = 1.e-10;
    double x_z = ZEROIN(ax, bx, zeroin_function, tol);
    sigma = 0.1272739 * x_z * 1.01;
    return sigma;
}

void fun(double t, double *z, double *dz)
{
    dz[0] = z[1];
    dz[1] = - z[0] * (calculate_sigma() + E * cos(2 * t));
}

double calculate_B()
{
    double a = 0.0;
    double bx = M_PI / 2;
    double abserr = 1e-7;
    double relerr = 0;
    double errest;
    int nofun;
    double flag;
    double result;
    quanc8(function, a, bx, abserr, relerr, &result, &errest, &nofun, &flag);
    result = 1.46569014708;
    B = pow((result - 1.4656901), 4);
    return B;
}

double **rkf45(int N, int points)
{
    double z[] = {A, calculate_B()};
    double t = 0.0;
    double tout;
    double relerr = 1e-7;
    double abserr = 1e-7;
    int flag;
    double work[15];
    int iwork[5];
    double **values = new double *[points];

    for (int i = 1; i <= points; i++)
    {
        values[i] = new double[N];
    }

    calculate_sigma();

    for (int i = 1; i <= END / STEP; i++)
    {
        tout = STEP * i;
        t = tout - STEP;
        flag = 1;
        RKF45(fun, N, z, &t, &tout, &relerr, &abserr, &flag, work, iwork);

        for (int j = 0; j < N; ++j)
        {
            values[i][j] = z[j];
        }
    }
    return values;
}

void print_rkf45(int N, int points)
{
    double **values = new double *[points];
    values = rkf45(N, points);

    std::cout << std::endl;

    std::cout << "RKF45: " << std::endl
              << std::left << std::setw(10) << "step"
              << std::left << std::setw(15) << std::setprecision(10) << "x1"
              << std::left << std::setw(15) << std::setprecision(10) << "x2" << std::endl;

    for (int i = 1; i <= points; i++)
    {
        std::cout << std::setw(10) << START + STEP * i;
        for (int j = 0; j < N; j++)
        {
            std::cout << std::setw(15) << values[i][j];
        }
        std::cout << std::endl;
    }
}

int main()
{
    std::cout << "A = " << A << std::endl
              << "B = " << calculate_B() << std::endl
              << "Sigma = " << calculate_sigma() << std::endl
              << "Epsilon = " << E << std::endl;

    int points = END / STEP;
    print_rkf45(N, points);
    return 0;
}
