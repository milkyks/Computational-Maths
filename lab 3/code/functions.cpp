#include "functions.hpp"
#include <cmath>

double fun1(double t, double x1, double x2)
{
    return (-4 * x1 + 23 * x2 + exp(-t));
}

double fun2(double t, double x1, double x2)
{
    return (4 * x1 - 48 * x2 + sin(t));
}