#include "lagrange.hpp"

double calculateLagrangePolynom(int N, double x, double nodes[], double values[])
{
  double sum_of_polynoms = 0;
  double polynom = 1;

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      if (j == i)
      {
        continue;
      }
      else
      {
        polynom *= ((x - nodes[j]) / (nodes[i] - nodes[j]));
      }
    }
    sum_of_polynoms += polynom * values[i];
    polynom = 1;
  }

  return sum_of_polynoms;
}
