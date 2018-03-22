/*
 Ignacio Cordón
 Bartolomé Ortiz
 Julián Pozuelo
*/
#include<cmath>
#include<vector>
#include<numeric>
#include<iostream>
#include<iterator>
#include "gnuplot-iostream.h"
using namespace std;

/*
 Solve problem
 u''(x) = f(x)
*/
double f(double x) {
  return exp(x) * sin(x);
}


/*
True solution of the problem
 u''(x) = f(x)
 u(-1) = 1
 u(-1) = 1
*/
double TrueSolution(double x) {
  double c_1 = 1 + 1.0/4 * (exp(-1) * cos(-1) + exp(1) * cos(1));
  double c_2 = 1.0/4 * (exp(1) * cos(1) - exp(-1) * cos(-1));
  return c_2 * x + c_1 - 0.5 * exp(x) * cos(x);
}

/*
Given a lists of coefficients, it evaluates the the linear
combination of chebyshev functions in x
*/
double EvalChebyshevSeries(vector<double> coeffs, double x) {
  double result = 0;
  double T_nold = 1;
  double T_n = x;
  double new_T_n;
  result += coeffs[0] / 2 * T_nold;
  result += coeffs[1] * T_n;
  
  for (int j = 2; j < coeffs.size(); j++) {
    new_T_n = 2 * x * T_n - T_nold;
    T_nold = T_n;
    T_n = new_T_n;
    
    result += coeffs[j] * new_T_n;
  }

  return result;
}
  
/*
 Computes the coefficients of the linear combination of Chebyshev's
 polynomials to approximate f, using a and b as boundary conditions
 u(-1) = a
 u(+1) = b
 O(n²) implementation
*/
vector<double> ChebyshevCoeffs(int n, double a, double b) {
  vector<double> collocation_point(n);
  vector<double> fs(n);
  vector<vector<double> > chebyshev_values(n, vector<double>(n));
  // what we are trying to compute
  vector<double> coeffs(n + 1, 0);
  double rhs;
  double alternate_coeffs_sum = 0;
  double coeffs_sum = 0;
  int signature = 1;
  
  // Compute collocation points x_j, j = 1, ..., n-1 and f(x_j)
  for (int j = 1; j <= n - 1; j++){
    collocation_point[j] = cos((j - 0.5) * M_PI / (n - 1));
    fs[j] = f(collocation_point[j]);      
  }
  
  // Compute value of T_r(x_j) with T_r the chebyshev polynom
  // using recursion T_{r + 1}(x) = 2x T_{r}(x) - T_{r - 1}(x)
  for (int j = 1; j <= n - 1; j++) {
    chebyshev_values[0][j] = 1;
    chebyshev_values[1][j] = collocation_point[j];
      
    for (int r = 2; r <= n - 2; r++) {
      chebyshev_values[r][j] = 2 * collocation_point[j] * chebyshev_values[r - 1][j] -
                              chebyshev_values[r - 2][j];
    }
  }

  // Solve system backwards to compute coeffs[n]...coeffs[2]
  for (int r = n - 2; r >= 0; r--) {
    rhs = inner_product(chebyshev_values[r].begin() + 1, chebyshev_values[r].begin() + n,
                        fs.begin() + 1, 0.0);
    rhs *= 2.0 / (n - 1);

    for (int k = r + 4; k <= n; k += 2)
      rhs -= coeffs[k] * (k - r) * k * (k + r);
    
    coeffs[r + 2] = rhs / (2.0 * (r + 2) * (2 * r + 2));
  }


  // Compute coeffs[0] and coeffs[1]
  for (int k = 2; k <= n; k++) {
    alternate_coeffs_sum += signature * coeffs[k];
    coeffs_sum += coeffs[k];
    signature *= -1;
  }

  coeffs[0] = (b + a - alternate_coeffs_sum - coeffs_sum);
  coeffs[1] = (b - a + alternate_coeffs_sum - coeffs_sum) / 2.0;

  return coeffs;
}


int main() {
  int n = 10;
  int a = 1;
  int b = 1;
  vector<double> coeffs = ChebyshevCoeffs(n, a, b);
  // Points where to compute approximation
  vector<pair<double, double>> approximation;
  vector<pair<double, double>> utrue;  

  // Equally-spaced points from -1 to 1
  for (double x = -1.0; x <= 1.01; x += 0.05){
    approximation.push_back({x, EvalChebyshevSeries(coeffs, x)});
    utrue.push_back({x, TrueSolution(x)});
  }

  // Print result of discrete approximation in comparison with true solution
  Gnuplot gp;

  gp << "set xrange [-1:1]\nset yrange [0.85:1.05]\n";
  gp << "plot" << gp.file1d(approximation) << "with lines title 'approximation',"
     << gp.file1d(utrue) << "with points title 'true solution'," << endl;
  
  
  return 0;
}
