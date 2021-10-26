
# 
# Seperatable: Rho_S Rho_B
#   https://kaskr.github.io/adcomp/classdensity_1_1SEPARABLE__t.html
#   
#   https://kaskr.github.io/adcomp/_book/Densities.html
#   
#   
#   https://kaskr.github.io/adcomp/_book/Densities.html#separable-construction-of-covariance-precision-matrices
  
  # 
  # using namespace density;
  # int n_s = 10;                   // Number of grid points in space
  # int n_t = 10;                   // Number of grid points in time
  # Type rho_s = 0.2;               // Correlation in space
  # Type rho_t = 0.4;               // Correlation in time
  # 
  # array<Type> x(n_s,n_t);
  # x.setZero();                    // x = 0
  # 
  # res = SEPARABLE(AR1(rho_t),AR1(rho_s))(x);
  # 
  
  