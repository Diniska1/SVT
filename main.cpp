#include "inmost.h"
#include <stdio.h>
#include <cmath>

/*
mkdir build
cd build
cmake ..
make
./main
*/

using namespace INMOST;
using namespace std;

double f(double x, double y)
{
    return sin(5 * x) * sin(5 * y);
}

double d2f(double x, double y)
{
    return 50 * sin(5 * x) * sin(5 * y);
}

void fill_in_A(Sparse::Matrix & A, unsigned step)
{
    for (int i = 0; i < step; ++i) {
        for (int j = 0; j < step; ++j) {
            A[(i * step + j)][(i * step + j)] = 4; 
            if (j + 1 < step) {
                A[(i * step + j)][(i * step + (j + 1))] = -1;
                A[(i * step + j + 1)][(i * step + j)] = -1;
            }
            if (i + 1 < step) {
                A[((i + 1) * step + j)][(i * step + j)] = -1;
                A[(i * step + j)][((i + 1) * step + j)] = -1;
            }
        }
    }
}

void fill_in_b(Sparse::Vector & b, unsigned step, double h)
{
    for(unsigned i = 1; i < step + 1; i++){
        for (int j = 1; j < step + 1; j++) {
            b[(i - 1) * step + (j - 1)] = d2f(i * h, j * h) * h * h;
            if (j == step) b[(i - 1) * step + (j - 1)] += f(i * h, 1);
            if (i == step) b[(i - 1) * step + (j - 1)] += f(1, j * h);
        }
    }
}

void norms(Sparse::Vector & sol, unsigned step)
{
    double max = 0, l2 = 0;
    double h = 1.0 / (step + 1);
    for(unsigned i = 1; i < step + 1; i++){
        for (int j = 1; j < step + 1; j++) {
            double t = abs(f(i * h, j * h) - sol[(i - 1) * step + (j - 1)]);
            l2 += t;
            if (t > max) max = t;
        }
    }
    cout << "l2: " << sqrt(l2) * h << endl;
    cout << "ch: " << max << endl;

}


int main(int argc, char *argv[])
{
    // Get number of nodes
    unsigned N;
    cout << "N: ";
    cin >> N;
    unsigned int step = N - 1;

    double h = 1.0 / (step + 1);

    // Create sparse matrix, RHS vector and solution vector
    Sparse::Matrix A;
    Sparse::Vector b;
    Sparse::Vector sol;
    // Set their size
    unsigned long long size = step * step;
    A.SetInterval(0, size);
    b.SetInterval(0, size);
    sol.SetInterval(0, size);

    
    fill_in_A(A, step);

    fill_in_b(b, step, h);



    Solver S(Solver::INNER_MPTILUC);
    S.SetParameter("absolute_tolerance", "1e-10");
    S.SetParameter("relative_tolerance", "1e-6");

    S.SetMatrix(A);

    // Solve
    bool solved = S.Solve(b, sol);
    cout << "num.iters: " << S.Iterations() << endl;
    cout << "prec.time: " << S.PreconditionerTime() << endl;
    cout << "iter.time: " << S.IterationsTime() << endl;
    if(!solved){
        cout << "Linear solver failure!" << endl;
        cout << "Reason: " << S.ReturnReason() << endl;
    }
    
    norms(sol, step);

	return 0;
}