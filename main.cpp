#include <algorithm>
#include <cmath>
#include "Problem.h"

int main()
{
    Parameters params("0.txt");
    Problem problem(params);
    problem.Solve();
    return 0;
}