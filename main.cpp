#include <cmath>
#include <iostream>
#include <iomanip>
#include "lagrange.hpp"
#include "quanc8.hpp"
#include "spline.hpp"

double function(double x) 
{
  return sin(x * x);
}

int main()
{
  double nodes[7] = {};

  double node = 1.5;
  for (int i = 0; i < 7; i++)
  {
    nodes[i] = node;
    node += 0.2;
  }

  double values[7] = {};
  
  double a = 0.0;
  double b;
  double abserr = 1e-4;
  double relerr = 0;
  double errest;
  int nofun;
  double flag;

  double result;
  for (int i = 0; i < 7; i++) // подсчет интеграла с пределами интегрирования от 0 до 1.5, 1.7 и т.д.
  {
    b = nodes[i];
    quanc8(function, a, b, abserr, relerr, &result, &errest, &nofun, &flag);
    values[i] = result;
  }

  std::cout << "Values of integral at first: " << "\n";
  for (int i = 0; i < 7; i++)
  {
    std::cout << values[i] << "\n";
  }
  
  std::cout << "_______________________________________________________________\n";

  double x_1[7] = {};  // массив точек, в которых нужно считать функцию, полином Лагранжа и сплайн-функцию

  node = 1.6;

  for (int i = 0; i < 7; i++)
  {
    x_1[i] = node;
    node += 0.2;
  }

  std::cout << std::setw(12) << std::left << "f(x)"
            << std::setw(13) << "Lagrange"
            << std::setw(13) << "Spline"
            << std::setw(14) << "|f(x)-L(x)|"
            << std::setw(14) << "|f(x)-S(x)|"
            << std::endl;

  double x_k[8] = {}; // для сплайн-функции
  double f_k[8] = {};
  double b_k[7] = {};
  double c_k[7] = {};
  double d_k[7] = {};

  for (int i = 0; i < 7; i++)
  {
    x_k[i+1] = nodes[i];
    f_k[i+1] = values[i];
  }

  a = 0.0;
  b = 1.6;
  for (int i = 0; i < 7; i++)
  {
    quanc8(function, a, b, abserr, relerr, &result, &errest, &nofun, &flag);
    std::cout << std::setw(12) << result
              << std::setw(13) << calculateLagrangePolynom(x_1[i], nodes, values); 
    spline(7, x_k, f_k, b_k, c_k, d_k);
    std::cout << std::setw(13) << seval(7, x_1[i], x_k, f_k, b_k, c_k, d_k)
              << std::setw(14) << fabs(result - calculateLagrangePolynom(x_1[i], nodes, values))
              << std::setw(14) << fabs(result - seval(7, x_1[i], x_k, f_k, b_k, c_k, d_k)) << "\n";

    b += 0.2;
  }
 
  return 0;
}