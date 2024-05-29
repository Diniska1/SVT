#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>

// compile + run:
// g++ NAMEFILE.cpp && ./a.out

double u(double x)
{
    return x + sin(x);
}

double f(double x)
{
    return sin(x);
}


void 
progonka(double * x, double * b, int n)
{
    double * v = new double[n];
    double * u = new double[n];

    v[0] = (-1) / (-2);
    u[0] = -b[0] / (-2);

    for(int i = 1; i < n - 1; ++i) {
        v[i] = (-1) / (-2 - (-1) * v[i - 1]);
        u[i] = ((-1) * u[i - 1] - b[i]) / (-2 - (-1) * v[i - 1]);
    }

    v[n - 1] = 0;
    u[n - 1] = ((-1) * u[n - 2] - b[n - 1]) / (-2 - (-1) * v[n - 2]);

    x[n - 1] = u[n - 1];
    for(int i = n - 1; i >= 1; --i) {
        x[i - 1] = v[i - 1] * x[i] + u[i - 1];
    }
}

double C_norm(double u_accurate(double), double *u_calculated, int N) { 
    double h = 1 / (double) N; 
    double norm = std::abs(u_accurate(h) - u_calculated[0]); 
    for (int i = 2; i < N; i++) { 
        double temp = std::abs(u_accurate(i * h) - u_calculated[i]); 
        if (temp > norm) { 
            norm = temp; 
        } 
    } 
    return norm; 
} 


double L_2_norm(double u_accurate(double), double* u_calculated, int N) { 
    double h = 1 / (double) N; 
    double norm = 0; 
    for (int i = 1; i < N - 1; i++) { 
        norm += std::pow((u_accurate(i * h) - u_calculated[i]), 2); 
    } 
    norm *= h;
    return std::sqrt(norm); 
}

void 
find_x(double * x, int n)
{
    double h = 1.0 / n;
        
    double * b = new double[n - 1];
    b[0] = f(h) * h * h + u(0);
    b[n - 1] = f((n - 1) * h) * h * h + u(1);
    for(int i = 1; i < n - 1; i ++) {
        b[i] = f((i + 1) * h) * h * h;
    }


    progonka(x, b, n);
}

int main()
{
    for(int n = 10; n < 10e6; n *= 10) {
        double * x = new double[n + 1];
        find_x(x, n);

        std::cout << n << "\t" << C_norm(u, x, n) << "    " << L_2_norm(u, x, n) << std::endl;
    }
    
    std::cout << std::endl;

}