/*
 Ignacio Cordón
 Bartolomé Ortiz
 Julián Pozuelo
*/
#include<cmath>
#include<vector>
#include<numeric>
using namespace std;

/*
 Solve problem
 u''(x) = f(x)
*/
double f(double x) {
  return exp(x) * sin(x);
}

// u(-1) = a
// u(+1) = b
// O(n²) implementation
vector<double> ChebyshevCoeffs(int n, double a, double b) {
  vector<double> collocation_point(n - 1);
  vector<double> fs(n);
  vector<vector<double> > chebyshev_values(n, n);
  // what we are trying to compute
  vector<double> coeffs(n + 1, 0);
  double rhs;
  double alternate_coeffs_sum = 0;
  int signature = 1;
  double coeffs_sum = 0;
  
  // Compute collocation points x_j, j = 1, ..., n-1 and f(x_j)
  for (int j = 1; j <= n - 1; j++){
    collocation_point[j] = cos((j - 1/2) * M_PI / (n - 1));
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
  for (int r = n - 2; r >= 0; r--){
    rhs = inner_product(chebyshev_values[r].begin() + 1, chebyshev_values[r].begin() + n - 1,
                        fs.begin() + 1, 0);
    rhs *= 2 / (n - 1);

    for (int k = r + 4; k <= n; k += 2){
      rhs -= coeffs[k] * (k - r) * k * (k + r);
    }
    
    coeffs[r + 2] = rhs / (2.0 * (r + 2) * (2 * r + 2));
  }


  // Compute coeffs[0] and coeffs[1]
  for (int k = 2; k <= n; k++) {
    alternate_coeffs_sum += signature * coeffs[k];
    coeffs_sum += coeffs[k];
    signature *= -1;
  }

  coeffs[0] = (b + a - alternate_coeffs_sum - coeffs_sum) / 2.0;
  coeffs[1] = (b - a + alternate_coeffs_sum - coeffs_sum) / 2.0;

  return coeffs;
}

