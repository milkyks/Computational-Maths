#include <iostream>
#include <iomanip>
#include <cmath>
#include "rkf45.hpp"
#include "functions.hpp"

const int N = 2;
const double START = 0.0;
const double STEP = 0.1;
const double END = 2.0;

void fun(double t, double *x, double *dx)
{
  dx[0] = -4 * x[0] + 23 * x[1] + exp(-t);
  dx[1] = 4 * x[0] - 48 * x[1] + sin(t);
}

double **rkf45(int N, int points)
{
  double x[] = {1, 0};
  double t = 0.0;
  double tout;
  double relerr = 1e-4;
  double abserr = 1e-4;
  int flag;
  double work[15];
  int iwork[5];
  double **values = new double *[points];

  for (int i = 0; i < points; i++)
  {
  values[i] = new double[N];
  }

  for (int i = 0; i < points; i++)
  {
    tout = START + STEP * i;
    flag = 1;
    RKF45(fun, N, x, &t, &tout, &relerr, &abserr, &flag, work, iwork);

    for (int j = 0; j < N; ++j)
    {
      values[i][j] = x[j];
    }
  }

  return values;
}

double **runge(int N, int points)
{
  double **values = new double *[points];
  double x[] = {1, 0};
  double k1[N];
  double k2[N];
  double k3[N];
  double tout;

  for (int i = 0; i < points; i++)
  {
    values[i] = new double[N];
  }

  for (int i = 0; i < N; ++i)
  {
    values[0][i] = x[i];
  }

  for (int i = 0; i < points; i++)
  {
    tout = START + STEP * i;
    k1[0] = STEP * fun1(tout, x[0], x[1]);
    k1[1] = STEP * fun2(tout, x[0], x[1]);

    k2[0] = STEP * fun1(tout + STEP / 2, x[0] + k1[0] / 2, x[1] + k1[1] / 2);
    k2[2] = STEP * fun2(tout + STEP / 2, x[0] + k1[0] / 2, x[1] + k1[1] / 2);

    k3[0] = STEP * fun1(tout + STEP, x[0] - k1[0] + 2 * k2[0], x[1] - k1[1] + 2 * k2[1]);
    k3[1] = STEP * fun2(tout + STEP, x[0] - k1[0] + 2 * k2[0], x[1] - k1[1] + 2 * k2[1]);

    for (int j = 0; j < N; ++j)
    {
    values[i][j] = x[j];
    }

    x[0] += (k1[0] + 4 * k2[0] + k3[0]) / 6;
    x[1] += (k1[1] + 4 * k2[1] + k3[1]) / 6;
  }

  return values;
}

void print_runge(int N, int points)
{
  double **values = new double *[points];
  values = runge(N, points);

  std::cout << "Runge Kutta method 3 degrees of accuracy: " << std::endl
            << std::left << std::setw(10) << "step"
            << std::left << std::setw(15) << "x1"
            << std::left << std::setw(15) << "x2" << std::endl;

  for (int i = 0; i < points; i++)
  {
    std::cout << std::setw(10) << START + STEP * i;
    for (int j = 0; j < N; j++)
    {
      std::cout << std::setw(15) << values[i][j];
    }
    std::cout << std::endl;
  }
}

void print_rkf45(int N, int points)
{
  double **values = new double *[points];
  values = rkf45(N, points);

  std::cout << "RKF45: " << std::endl
            << std::left << std::setw(10) << "step"
            << std::left << std::setw(15) << "x1"
            << std::left << std::setw(15) << "x2" << std::endl;

  for (int i = 0; i < points; i++)
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
  int points = (END - START) / STEP + 1;
  print_rkf45(N, points);
  print_runge(N, points);
  return 0;
}



