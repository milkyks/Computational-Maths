#include <cmath>
#include <iostream>
#include <iomanip>
#include "lagrange.hpp"
#include "quanc8.hpp"
#include "spline.hpp"

double function(double x) 
{
  return sin(x*x);
}

int main()
{
  const int N = 7;
  double nodes[N] = {};
  double node = 1.5;

  for (int i = 0; i < N; i++)
  {
    nodes[i] = node;
    node += 0.2;
  }

  double values[N] = {};
  
  double a = 0.0;
  double b;
  double abserr = 1e-10;
  double relerr = 0;
  double errest;
  int nofun;
  double flag;
  double result;

  for (int i = 0; i < N; i++) // подсчет интеграла с пределами интегрирования от 0 до 1.5, 1.7 и т.д.
  {
    b = nodes[i];
    quanc8(function, a, b, abserr, relerr, &result, &errest, &nofun, &flag);
    values[i] = result;
  }

  std::cout << "Values of integral at first: " << "\n";

  std::cout << std::setw(12) << "x"
            << std::setw(12) << "f(x)\n";

  for (int i = 0; i < N; i++)
  {
    std::cout << std::setw(12) << nodes[i]
              << std::setw(12) << values[i] << "\n";
  }
  
  std::cout << "_______________________________________________________________\n";

  double x_1[N-1] = {};  // массив точек, в которых нужно считать функцию, полином Лагранжа и сплайн-функцию
  node = 1.6;

  for (int i = 0; i < N-1; i++)
  {
    x_1[i] = node;
    node += 0.2;
  }

  std::cout << std::setw(6) << std::left << "x"
            << std::setw(12) << "f(x)"
            << std::setw(12) << "Lagrange"
            << std::setw(12) << "Spline"
            << std::setw(14) << "|f(x)-L(x)|"
            << std::setw(14) << "|f(x)-S(x)|"
            << std::endl;

  double x_spline[N+1] = {}; // для сплайн-функции
  double f_spline[N+1] = {};
  double b_k[N] = {};
  double c_k[N] = {};
  double d_k[N] = {};

  for (int i = 0; i < N; i++)
  {
    x_spline[i+1] = nodes[i];
    f_spline[i+1] = values[i];
  }

  a = 0.0;
  b = 1.6;
  
  for (int i = 0; i < N-1; i++)
  {
    std::cout << std::setw(6) << x_1[i];
    quanc8(function, a, b, abserr, relerr, &result, &errest, &nofun, &flag);
    std::cout << std::setw(12) << result
              << std::setw(12) << calculateLagrangePolynom(N, x_1[i], nodes, values);
    spline(N, x_spline, f_spline, b_k, c_k, d_k);
    std::cout << std::setw(12) << seval(N, x_1[i], x_spline, f_spline, b_k, c_k, d_k)
              << std::setw(14) << fabs(result - calculateLagrangePolynom(N, x_1[i], nodes, values))
              << std::setw(14) << fabs(result - seval(N, x_1[i], x_spline, f_spline, b_k, c_k, d_k)) << "\n";

    b += 0.2;
  }
 
  return 0;
}