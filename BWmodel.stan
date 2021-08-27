//Stan model for a simple exponential model of dChl/dt based on Borman an Webster 1999
      
      data{
      int < lower = 1 > N; // Sample size
      vector[N] H; // Predictor
      vector[N] conc; // Outcome
}

parameters {
 real alpha; // Intercept
 real beta; // Slope (regression coefficients)
 real < lower = 0 > sigma; // Error SD
}

model {
 conc ~ normal(alpha + H * beta , sigma);
}

generated quantities {
} // The posterior predictive distribution
