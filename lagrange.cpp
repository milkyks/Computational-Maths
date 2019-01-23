#include "lagrange.hpp"

double calculateLagrangePolynom(double x, double nodes[], double values[])
{
  double sum_of_polynoms = 0;
  double polynom = 1;

  for (int i = 0; i < 7; i++)
  {
    for (int j = 0; j < 7; j++)
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
