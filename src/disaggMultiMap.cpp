// EXACTLY FROM DESKTOP THAT WORKED
// Data: Spatial field mesh and matrices, polygon data, covariate pixel data

#define TMB_LIB_INIT R_init_disaggMultiMap
#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  // ------------------------------------------------------------------------ //
  // Spatial field data
  // ------------------------------------------------------------------------ //
  
  // The A matrices are for projecting the mesh to a point for the pixel and point data respectively.
  DATA_SPARSE_MATRIX(Apixel);
  DATA_STRUCT(spde, spde_t);
  
  // ------------------------------------------------------------------------ //
  // Polygon level data
  // ------------------------------------------------------------------------ //
  
  // Covariate pixel data (combined across time)
  DATA_MATRIX(x);
  
  // Two-column integer array with start and end indices (already computed in R).
  DATA_IARRAY(start_end_index);
  
  // Polygon data: responses and (if applicable) sample sizes.
  DATA_VECTOR(polygon_response_data);
  DATA_VECTOR(response_sample_size);
  
  // Aggregation values (combined across time).
  DATA_VECTOR(aggregation_values);
  
  // ------------------------------------------------------------------------ //
  // Likelihood and link functions
  // ------------------------------------------------------------------------ //
  
  DATA_INTEGER(family);
  DATA_INTEGER(link);
  
  // ------------------------------------------------------------------------ //
  // NEW for multiple maps: Temporal data (one entry per polygon)
  // ------------------------------------------------------------------------ //
  DATA_VECTOR(time);
  
  // ------------------------------------------------------------------------ //
  // Parameters
  // ------------------------------------------------------------------------ //
  
  PARAMETER(intercept);
  PARAMETER_VECTOR(slope);
  
  DATA_SCALAR(priormean_intercept);
  DATA_SCALAR(priorsd_intercept);
  DATA_SCALAR(priormean_slope);
  DATA_SCALAR(priorsd_slope);
  
  // Likelihood for Gaussian likelihood hyperparameter.
  PARAMETER(log_tau_gaussian);
  Type tau_gaussian = exp(log_tau_gaussian);
  Type gaussian_sd = 1 / sqrt(tau_gaussian);
  
  // INLA defines a loggamma prior on log tau.
  Type prior_gamma_shape = 1;
  Type prior_gamma_rate = 5e-05;
  
  PARAMETER_VECTOR(iideffect);
  PARAMETER(iideffect_log_tau);
  Type iideffect_tau = exp(iideffect_log_tau);
  Type iideffect_sd = 1 / sqrt(iideffect_tau);
  Type iideffect_mean = 0.0;
  
  DATA_SCALAR(prior_iideffect_sd_max);
  DATA_SCALAR(prior_iideffect_sd_prob);
  
  // SPDE hyperparameters
  PARAMETER(log_sigma);
  PARAMETER(log_rho);
  Type sigma = exp(log_sigma);
  Type rho = exp(log_rho);
  
  DATA_SCALAR(prior_rho_min);
  DATA_SCALAR(prior_rho_prob);
  DATA_SCALAR(prior_sigma_max);
  DATA_SCALAR(prior_sigma_prob);
  
  DATA_SCALAR(nu);
  Type kappa = sqrt(8.0) / rho;
  
  // Random field parameters
  PARAMETER_VECTOR(nodemean);
  
  // Model component flags
  DATA_INTEGER(field);
  DATA_INTEGER(iid);
  
  // ------------------------------------------------------------------------ //
  // Other dimensions
  // ------------------------------------------------------------------------ //
  
  // Number of polygons (across all time points)
  int n_polygons = polygon_response_data.size();
  // Number of pixels (all time points combined)
  int n_pixels = x.rows();
  
  Type nll = 0.0;
  
  // ------------------------------------------------------------------------ //
  // Priors on fixed effects and hyperparameters
  // ------------------------------------------------------------------------ //
  
  nll -= dnorm(intercept, priormean_intercept, priorsd_intercept, true);
  for (int s = 0; s < slope.size(); s++) {
    nll -= dnorm(slope[s], priormean_slope, priorsd_slope, true);
  }

  if(iid && family != 3) {
    // Likelihood of hyperparameter of polygon iid random effect.
    // From https://projecteuclid.org/euclid.ss/1491465621 (Eqn 3.3)
    Type lambda = -log(prior_iideffect_sd_prob) / prior_iideffect_sd_max;
    Type log_pcdensity_iid = log(lambda / 2) - (3.0/2.0)*iideffect_log_tau - lambda * pow(iideffect_tau, -1.0/2.0);
    nll -= log_pcdensity_iid + iideffect_log_tau;
    for (int p = 0; p < iideffect.size(); p++) {
      nll -= dnorm(iideffect[p], iideffect_mean, iideffect_sd, true);
    }
  }
  
  if(family == 0) {
    nll -= dgamma(tau_gaussian, prior_gamma_shape, prior_gamma_rate, true);
  }
  
  if(field) {
    // Likelihood of hyperparameters for field.
    // From https://www.tandfonline.com/doi/full/10.1080/01621459.2017.1415907 (Theorem 2.6)
    Type lambdatilde1 = -log(prior_rho_prob) * prior_rho_min;
    Type lambdatilde2 = -log(prior_sigma_prob) / prior_sigma_max;
    Type log_pcdensity = log(lambdatilde1) + log(lambdatilde2) - 2 * log_rho - lambdatilde1 * pow(rho, -1) - lambdatilde2 * sigma;
    nll -= log_pcdensity + log_rho + log_sigma;
    
    // build spde matrix
    SparseMatrix<Type> Q = Q_spde(spde, kappa);

    // From Lindgren (2011) https://doi.org/10.1111/j.1467-9868.2011.00777.x, see equation for the marginal variance
    Type scaling_factor = sqrt(exp(lgamma(nu)) / (exp(lgamma(nu + 1)) * 4 * M_PI * pow(kappa, 2 * nu)));

    // Likelihood of the random field.
    nll += SCALE(GMRF(Q), sigma / scaling_factor)(nodemean);
  }

  // if(family==3) { // old PC prior
  //   Type lambda = -log(prior_iideffect_sd_prob)
  //                 / prior_iideffect_sd_max;
  //     // log π(η) = log(λ/2) + 0.5·η − λ·exp(η/2)
  //     Type log_pcdensity_nb =
  //         log(lambda / Type(2.0))
  //       + Type(0.5) * iideffect_log_tau
  //       - lambda * exp(iideffect_log_tau * Type(0.5));
  //     nll -= log_pcdensity_nb;
  // }

  if (family == 3) {  // NB PC prior   AUG 1
    Type tau_nb = exp(iideffect_log_tau);

    // PC‐prior on tau = exp(iideffect_log_tau), tail P(τ>τ_max)=p_tail
    Type lambda_nb = -log(prior_iideffect_sd_prob) / prior_iideffect_sd_max;
    
    nll += log(Type(2.0)) + lambda_nb * sqrt(tau_nb) + Type(0.5) * log(tau_nb) - log(lambda_nb);
  }

  // NEW GAMMA
  // if (family == 3) {  // NB Gamma prior  
  //   Type tau_nb = exp(iideffect_log_tau);

  //   // tau_nb prior Gamma(a=1,b=2):
  //   Type a = Type(1.0);
  //   Type b = Type(2.0);
  //   // log‐density of tau_nb ~ Gamma(a,b):  (dgamma returns density in tau‐space)
  //   nll -= dgamma(tau_nb, a, b, true);
  //   // + Jacobian term:  dτ_nb/dη = τ_nb, so subtract log(τ_nb)
  //   nll -= log(tau_nb);
  // }
  
  Type nll_priors = nll;
  
  // ------------------------------------------------------------------------ //
  // Data likelihood: Pixel-level predictions and polygon aggregation
  // ------------------------------------------------------------------------ //
  
  vector<Type> pixel_linear_pred(n_pixels);
  pixel_linear_pred = intercept + x * slope;
  
  if(field) {
    vector<Type> linear_pred_field(n_pixels);
    linear_pred_field = Apixel * nodemean;
    pixel_linear_pred += linear_pred_field.array();
  }
  
  // Adjust start_end_index: column 2 becomes (end - start + 1)
  start_end_index.col(1) = start_end_index.col(1) - start_end_index.col(0) + 1;
  
  // Variables for polygon-level aggregation
  Type polygon_response;
  Type normalised_polygon_response;
  Type normalisation_total;
  Type pred_polygoncases;
  Type pred_polygonrate;
  Type polygon_sd;
  vector<Type> pixel_pred;
  vector<Type> numerator_pixels;
  vector<Type> normalisation_pixels;
  vector<Type> reportnormalisation(n_polygons);
  vector<Type> reportprediction_cases(n_polygons);
  vector<Type> reportprediction_rate(n_polygons);
  vector<Type> reportnll(n_polygons);
  vector<Type> reportpolygonsd(n_polygons);
  
  // Loop over each polygon (each polygon now carries an associated time index)
  for (int polygon = 0; polygon < n_polygons; polygon++) {
    // Rprintf("Polygon %d at time %f\n", polygon, time[polygon]);
    
    // pixel-level predictions for this polygon
    pixel_pred = pixel_linear_pred.segment(start_end_index(polygon, 0), start_end_index(polygon, 1)).array();
    if(iid && family != 3) {
      pixel_pred += iideffect[polygon];
    }
    // Apply link function
    if(link == 0) {
      pixel_pred = invlogit(pixel_pred);
    } else if(link == 1) {
      pixel_pred = exp(pixel_pred);
    } else if(link == 2) {
      // identity link; do nothing
    } else {
      Rf_error("Link function not implemented.");
    }
    
    // Aggregate pixel predictions using the aggregation weights
    numerator_pixels = pixel_pred * aggregation_values.segment(start_end_index(polygon, 0), start_end_index(polygon, 1)).array();
    normalisation_pixels = aggregation_values.segment(start_end_index(polygon, 0), start_end_index(polygon, 1));
    normalisation_total = sum(normalisation_pixels);
    pred_polygoncases = sum(numerator_pixels);
    pred_polygonrate = pred_polygoncases / normalisation_total;
    
    reportnormalisation[polygon] = normalisation_total;
    reportprediction_cases[polygon] = pred_polygoncases;
    reportprediction_rate[polygon] = pred_polygonrate;
    
    // Compute likelihood contribution based on specified family
    if(family == 0) {
      polygon_sd = gaussian_sd * sqrt((normalisation_pixels * normalisation_pixels).sum()) / normalisation_total;
      reportpolygonsd[polygon] = polygon_sd;
      polygon_response = polygon_response_data(polygon);
      normalised_polygon_response = polygon_response / normalisation_total;
      nll -= dnorm(normalised_polygon_response, pred_polygonrate, polygon_sd, true);
      reportnll[polygon] = -dnorm(normalised_polygon_response, pred_polygonrate, polygon_sd, true);
    } else if(family == 1) {
      nll -= dbinom(polygon_response_data[polygon], response_sample_size[polygon], pred_polygonrate, true);
      reportnll[polygon] = -dbinom(polygon_response_data[polygon], response_sample_size[polygon], pred_polygonrate, true);
    } else if(family == 2) {
      nll -= dpois(polygon_response_data[polygon], pred_polygoncases, true);
      reportnll[polygon] = -dpois(polygon_response_data[polygon], pred_polygoncases, true);
    } else if(family == 3) {
      // OG HERE!
      // Type tau_nb = exp(iideffect_log_tau);

      // // tau_nb prior Gamma(a=1,b=0.5):
      // Type a = Type(1.0);
      // Type b = Type(2.0);
      // // log‐density of tau_nb ~ Gamma(a,b):  (dgamma returns density in τ‐space)
      // nll -= dgamma(tau_nb, a, b, true);
      // // + Jacobian term:  dτ_nb/dη = τ_nb, so subtract log(τ_nb)
      // // nll -= iideffect_log_tau;
      // nll -= log(tau_nb);

      // // NB‐likelihood:
      // Type y_i  = polygon_response_data[polygon];
      // Type mu_i  = pred_polygoncases;


      // Type r = tau_nb;          // e.g. φ = 1/τ_nb or φ = τ_nb depending on your convention
      // Type p = r / (r + mu_i); 

      // nll -= dnbinom(y_i, r, p, true);
      // reportnll[polygon] = -dnbinom(y_i, r, p, true);


      // // NICE PC PRIOR? FROM JULY 24
      // Type lambda = -log(prior_iideffect_sd_prob)
      //             / prior_iideffect_sd_max;
      // // log π(η) = log(λ/2) + 0.5·η − λ·exp(η/2)
      // Type log_pcdensity_nb =
      //     log(lambda / Type(2.0))
      //   + Type(0.5) * iideffect_log_tau
      //   - lambda * exp(iideffect_log_tau * Type(0.5));
      // nll -= log_pcdensity_nb;

      // // ——— Negative‐Binomial likelihood ———
      Type tau_nb = exp(iideffect_log_tau);
      Type r      = Type(1.0) / tau_nb;
      Type p      = Type(1.0) / (Type(1.0) + tau_nb * pred_polygoncases);

      // y_i ~ NB(r, p)
      nll -= dnbinom(polygon_response_data[polygon],
                     r, p,
                     true);
      reportnll[polygon] = 
          -dnbinom(polygon_response_data[polygon],
                   r, p,
                   true);



      
    } else {
      Rf_error("Likelihood not implemented.");
    }
  }
  
  // ------------------------------------------------------------------------ //
  // Reporting
  // ------------------------------------------------------------------------ //
  
  REPORT(reportprediction_cases);
  REPORT(reportprediction_rate);
  REPORT(reportnormalisation);
  REPORT(reportnll);
  REPORT(polygon_response_data);
  REPORT(nll_priors);
  REPORT(nll);
  if(family == 0) {
    REPORT(reportpolygonsd);
  }
  // Report the time vector for debugging and verification
  REPORT(time);
  
  return nll;
}




// NORMAL VERSION
// // Data: Spatial field mesh and matrices, polygon data, covariate pixel data

// #define TMB_LIB_INIT R_init_disaggMultiMap
// #include <TMB.hpp>

// template <class Type>
// Type objective_function<Type>::operator()()
// {
//   using namespace R_inla;
//   using namespace density;
//   using namespace Eigen;
  
//   // ------------------------------------------------------------------------ //
//   // Spatial field data
//   // ------------------------------------------------------------------------ //
  
//   // The A matrices are for projecting the mesh to a point for the pixel and point data respectively.
//   DATA_SPARSE_MATRIX(Apixel);
//   DATA_STRUCT(spde, spde_t);
  
//   // ------------------------------------------------------------------------ //
//   // Polygon level data
//   // ------------------------------------------------------------------------ //
  
//   // Covariate pixel data (combined across time)
//   DATA_MATRIX(x);
  
//   // Two-column integer array with start and end indices (already computed in R).
//   DATA_IARRAY(start_end_index);
  
//   // Polygon data: responses and (if applicable) sample sizes.
//   DATA_VECTOR(polygon_response_data);
//   DATA_VECTOR(response_sample_size);
  
//   // Aggregation values (combined across time).
//   DATA_VECTOR(aggregation_values);
  
//   // ------------------------------------------------------------------------ //
//   // Likelihood and link functions
//   // ------------------------------------------------------------------------ //
  
//   DATA_INTEGER(family);
//   DATA_INTEGER(link);
  
//   // ------------------------------------------------------------------------ //
//   // Temporal data (one entry per polygon)
//   // ------------------------------------------------------------------------ //
//   DATA_VECTOR(time);
  
//   // ------------------------------------------------------------------------ //
//   // Parameters
//   // ------------------------------------------------------------------------ //
  
//   PARAMETER(intercept);
//   PARAMETER_VECTOR(slope);
  
//   DATA_SCALAR(priormean_intercept);
//   DATA_SCALAR(priorsd_intercept);
//   DATA_SCALAR(priormean_slope);
//   DATA_SCALAR(priorsd_slope);
  
//   // Likelihood for Gaussian likelihood hyperparameter.
//   PARAMETER(log_tau_gaussian);
//   Type tau_gaussian = exp(log_tau_gaussian);
//   Type gaussian_sd = 1 / sqrt(tau_gaussian);
  
//   // INLA defines a loggamma prior on log tau.
//   Type prior_gamma_shape = 1;
//   Type prior_gamma_rate = 5e-05;
  
//   PARAMETER_VECTOR(iideffect);
//   PARAMETER(iideffect_log_tau);
//   Type iideffect_tau = exp(iideffect_log_tau);
//   Type iideffect_sd = 1 / sqrt(iideffect_tau);
//   Type iideffect_mean = 0.0;
  
//   DATA_SCALAR(prior_iideffect_sd_max);
//   DATA_SCALAR(prior_iideffect_sd_prob);
  
//   // SPDE hyperparameters
//   PARAMETER(log_sigma);
//   PARAMETER(log_rho);
//   Type sigma = exp(log_sigma);
//   Type rho = exp(log_rho);
  
//   DATA_SCALAR(prior_rho_min);
//   DATA_SCALAR(prior_rho_prob);
//   DATA_SCALAR(prior_sigma_max);
//   DATA_SCALAR(prior_sigma_prob);
  
//   DATA_SCALAR(nu);
//   Type kappa = sqrt(8.0) / rho;
  
//   // Random field parameters
//   PARAMETER_VECTOR(nodemean);
  
//   // Model component flags
//   DATA_INTEGER(field);
//   DATA_INTEGER(iid);
  
//   // ------------------------------------------------------------------------ //
//   // Other dimensions
//   // ------------------------------------------------------------------------ //
  
//   // Number of polygons (across all time points)
//   int n_polygons = polygon_response_data.size();
//   // Number of pixels (all time points combined)
//   int n_pixels = x.rows();
  
//   Type nll = 0.0;
  
//   // ------------------------------------------------------------------------ //
//   // Priors on fixed effects and hyperparameters
//   // ------------------------------------------------------------------------ //
  
//   nll -= dnorm(intercept, priormean_intercept, priorsd_intercept, true);
//   for (int s = 0; s < slope.size(); s++) {
//     nll -= dnorm(slope[s], priormean_slope, priorsd_slope, true);
//   }

//   if(iid && family != 3) {
//     // Likelihood of hyperparameter of polygon iid random effect.
//     // From https://projecteuclid.org/euclid.ss/1491465621 (Eqn 3.3)
//     Type lambda = -log(prior_iideffect_sd_prob) / prior_iideffect_sd_max;
//     Type log_pcdensity_iid = log(lambda / 2) - (3.0/2.0)*iideffect_log_tau - lambda * pow(iideffect_tau, -1.0/2.0);
//     nll -= log_pcdensity_iid + iideffect_log_tau;
//     for (int p = 0; p < iideffect.size(); p++) {
//       nll -= dnorm(iideffect[p], iideffect_mean, iideffect_sd, true);
//     }
//   }
  
//   if(family == 0) {
//     nll -= dgamma(tau_gaussian, prior_gamma_shape, prior_gamma_rate, true);
//   }

//   if (family == 3) {  // NB PC prior: tau ~ Exp(lambda), with eta = log(tau) (log scale)
//     Type eta = iideffect_log_tau;  // eta = log(tau)
//     Type lambda_nb = -log(prior_iideffect_sd_prob) / prior_iideffect_sd_max;

//     // pi(tau) = lambda * exp(-lambda * tau), tau = exp(eta)
//     // pi(eta) = pi(tau) * |dtau/deta| = lambda * exp(-lambda * exp(eta)) * exp(eta)
//     // => -log pi(eta) = lambda * exp(eta) - eta   (drop constant)
//     nll += lambda_nb * exp(eta) - eta;
//   }


//   if (family == 3) {  // Aug 12: NB PC prior: tau ~ Exp(lambda_nb) (non log)
//     Type tau_nb = exp(iideffect_log_tau);
//     Type lambda_nb = -log(prior_iideffect_sd_prob) / prior_iideffect_sd_max;

//     nll += lambda_nb * tau_nb - iideffect_log_tau; // drop the constant -log(lambda_nb)
//   }

//   // if (family == 3) {  // NB PC prior   AUG 1
//   //   Type tau_nb = exp(iideffect_log_tau);

//   //   // PC‐prior on tau = exp(iideffect_log_tau), tail P(τ>τ_max)=p_tail
//   //   Type lambda_nb = -log(prior_iideffect_sd_prob) / prior_iideffect_sd_max;
    
//   //   nll += log(Type(2.0)) + lambda_nb * sqrt(tau_nb) + Type(0.5) * log(tau_nb) - log(lambda_nb);
//   // }

//   // if (family == 3) {  // NB Gamma prior  
//   //   Type tau_nb = exp(iideffect_log_tau);

//   //   // tau_nb prior Gamma(a=1,b=2):
//   //   Type a = Type(1.0);
//   //   Type b = Type(2.0);
//   //   // log‐density of tau_nb ~ Gamma(a,b):  (dgamma returns density in tau‐space)
//   //   nll -= dgamma(tau_nb, a, b, true);
//   //   // + Jacobian term:  dτ_nb/dη = τ_nb, so subtract log(τ_nb)
//   //   nll -= log(tau_nb);
//   // }

//   // if (family == 3) { // OLD VERSION
//   //   Type eta = iideffect_log_tau;                    // eta = log tau
//   //   Type lambda_nb = -log(prior_iideffect_sd_prob) / prior_iideffect_sd_max;
//   //   // log π(eta) = log(lambda/2) + 0.5*eta − lambda*exp(0.5*eta)
//   //   Type log_pcdensity_nb = log(lambda_nb / Type(2.0))
//   //                         + Type(0.5) * eta
//   //                         - lambda_nb * exp(Type(0.5) * eta);
//   //   nll -= log_pcdensity_nb;
//   // }


//   if(field) {
//     // Likelihood of hyperparameters for field.
//     // From https://www.tandfonline.com/doi/full/10.1080/01621459.2017.1415907 (Theorem 2.6)
//     Type lambdatilde1 = -log(prior_rho_prob) * prior_rho_min;
//     Type lambdatilde2 = -log(prior_sigma_prob) / prior_sigma_max;
//     Type log_pcdensity = log(lambdatilde1) + log(lambdatilde2) - 2 * log_rho - lambdatilde1 * pow(rho, -1) - lambdatilde2 * sigma;
//     nll -= log_pcdensity + log_rho + log_sigma;
    
//     // build spde matrix
//     SparseMatrix<Type> Q = Q_spde(spde, kappa);

//     // From Lindgren (2011) https://doi.org/10.1111/j.1467-9868.2011.00777.x, see equation for the marginal variance
//     Type scaling_factor = sqrt(exp(lgamma(nu)) / (exp(lgamma(nu + 1)) * 4 * M_PI * pow(kappa, 2 * nu)));

//     // Likelihood of the random field.
//     nll += SCALE(GMRF(Q), sigma / scaling_factor)(nodemean);
//   }
  
//   Type nll_priors = nll;
  
//   // ------------------------------------------------------------------------ //
//   // Data likelihood: Pixel-level predictions and polygon aggregation
//   // ------------------------------------------------------------------------ //
  
//   vector<Type> pixel_linear_pred(n_pixels);
//   if (slope.size() > 0) { // if there are covariates
//     pixel_linear_pred = intercept + x * slope;
//   } else {
//     // no covariates means all pixels share the intercept
//     for (int i = 0; i < n_pixels; i++) {
//       pixel_linear_pred[i] = intercept;
//     }
//   }
  
//   if(field) {
//     vector<Type> linear_pred_field(n_pixels);
//     linear_pred_field = Apixel * nodemean;
//     pixel_linear_pred += linear_pred_field.array();
//   }
  
//   // Adjust start_end_index: column 2 becomes (end - start + 1)
//   start_end_index.col(1) = start_end_index.col(1) - start_end_index.col(0) + 1;
  
//   // Variables for polygon-level aggregation
//   Type polygon_response;
//   Type normalised_polygon_response;
//   Type normalisation_total;
//   Type pred_polygoncases;
//   Type pred_polygonrate;
//   Type polygon_sd;
//   vector<Type> pixel_pred;
//   vector<Type> numerator_pixels;
//   vector<Type> normalisation_pixels;
//   vector<Type> reportnormalisation(n_polygons);
//   vector<Type> reportprediction_cases(n_polygons);
//   vector<Type> reportprediction_rate(n_polygons);
//   vector<Type> reportnll(n_polygons);
//   vector<Type> reportpolygonsd(n_polygons);


  
//   // Loop over each polygon (each polygon carries an associated time index)
//   for (int polygon = 0; polygon < n_polygons; polygon++) {
//     // pixel-level predictions for this polygon
//     pixel_pred = pixel_linear_pred.segment(start_end_index(polygon, 0), start_end_index(polygon, 1)).array();
//     if(iid && family != 3) {
//       pixel_pred += iideffect[polygon];
//     }
//     // Apply link function
//     if(link == 0) {
//       pixel_pred = invlogit(pixel_pred);
//     } else if(link == 1) {
//       pixel_pred = exp(pixel_pred);
//     } else if(link == 2) {
//       // identity link; do nothing
//     } else {
//       Rf_error("Link function not implemented.");
//     }
    
//     // Aggregate pixel predictions using the aggregation weights
//     numerator_pixels = pixel_pred * aggregation_values.segment(start_end_index(polygon, 0), start_end_index(polygon, 1)).array();
//     normalisation_pixels = aggregation_values.segment(start_end_index(polygon, 0), start_end_index(polygon, 1));
//     normalisation_total = sum(normalisation_pixels);
//     pred_polygoncases = sum(numerator_pixels);
//     pred_polygonrate = pred_polygoncases / normalisation_total;
    
//     reportnormalisation[polygon] = normalisation_total;
//     reportprediction_cases[polygon] = pred_polygoncases;
//     reportprediction_rate[polygon] = pred_polygonrate;
    
//     // Compute likelihood contribution based on specified family
//     if(family == 0) {
//       polygon_sd = gaussian_sd * sqrt((normalisation_pixels * normalisation_pixels).sum()) / normalisation_total;
//       reportpolygonsd[polygon] = polygon_sd;
//       polygon_response = polygon_response_data(polygon);
//       normalised_polygon_response = polygon_response / normalisation_total;
//       nll -= dnorm(normalised_polygon_response, pred_polygonrate, polygon_sd, true);
//       reportnll[polygon] = -dnorm(normalised_polygon_response, pred_polygonrate, polygon_sd, true);
//     } else if(family == 1) {
//       nll -= dbinom(polygon_response_data[polygon], response_sample_size[polygon], pred_polygonrate, true);
//       reportnll[polygon] = -dbinom(polygon_response_data[polygon], response_sample_size[polygon], pred_polygonrate, true);
//     } else if(family == 2) {
//       nll -= dpois(polygon_response_data[polygon], pred_polygoncases, true);
//       reportnll[polygon] = -dpois(polygon_response_data[polygon], pred_polygoncases, true);
//     } else if(family == 3) {
//       // Aug 12 (log version)
//       Type y_i = polygon_response_data[polygon];
//       Type mu_i = pred_polygoncases;

//       Type eta = iideffect_log_tau;       // eta = log(tau)
//       Type tau_nb = exp(eta);             // tau
//       Type alpha_nb = tau_nb * tau_nb;    // alpha = tau^2

//       // r = 1/alpha, p = 1/(1 + alpha * mu)
//       Type r = Type(1.0) / alpha_nb;
//       Type p = Type(1.0) / (Type(1.0) + alpha_nb * mu_i);

//       nll -= dnbinom(y_i, r, p, true);
//       reportnll[polygon] = -dnbinom(y_i, r, p, true);

//       // r = 1/alpha, p = 1/(1 + alpha * mu)
//       // Type r = Type(1.0) / alpha_nb;
//       // Type p = Type(1.0) / (Type(1.0) + alpha_nb * mu_i);

//       // nll -= dnbinom(y_i, r, p, true);
//       // reportnll[polygon] = -dnbinom(y_i, r, p, true);

//       // // AUG 12 version (non log)
//       // Type tau_nb = exp(iideffect_log_tau);   // tau is sd of the gamma mixing
//       // Type y_i  = polygon_response_data[polygon];
//       // Type mu_i = pred_polygoncases;

//       // // alpha = tau^2
//       // Type alpha = tau_nb * tau_nb;
//       // // r = 1/alpha, p = 1/(1 + alpha * mu)
//       // Type r = Type(1.0) / alpha;
//       // Type p = Type(1.0) / (Type(1.0) + alpha * mu_i);

//       // nll -= dnbinom(y_i, r, p, true);


//       // AUG 1
//       // Type tau_nb = exp(iideffect_log_tau);

//       // Type y_i  = polygon_response_data[polygon];
//       // Type mu_i  = pred_polygoncases;


//       // Type r = Type(1.0)/tau_nb;          
//       // Type p = 1/(1+tau_nb*mu_i);


//       // nll -= dnbinom(y_i, r, p, true);
//       // reportnll[polygon] = -dnbinom(y_i, r, p, true);


//       // JUL 25
//       // Type tau_nb = exp(iideffect_log_tau);
//       // // α = r = 1/τ²
//       // Type r = Type(1.0) / (tau_nb * tau_nb);
//       // // p = 1 / (1 + τ² · μ)
//       // Type p = Type(1.0) 
//       //        / (Type(1.0) + tau_nb * tau_nb * pred_polygoncases);

//       // nll -= dnbinom(polygon_response_data[polygon],
//       //                r, p,
//       //                true);
//       // reportnll[polygon] =
//       //   -dnbinom(polygon_response_data[polygon],
//       //            r, p,
//       //            true);



    
      
//     } else {
//       Rf_error("Likelihood not implemented.");
//     }
//   }


//   // ------------------------------------------------------------------------ //
//   // Reporting
//   // ------------------------------------------------------------------------ //
  
//   REPORT(reportprediction_cases);
//   REPORT(reportprediction_rate);
//   REPORT(reportnormalisation);
//   REPORT(reportnll);
//   REPORT(polygon_response_data);
//   REPORT(nll_priors);
//   REPORT(nll);
//   if(family == 0) {
//     REPORT(reportpolygonsd);
//   }
//   // Report the time vector for debugging and verification
//   REPORT(time);
  
//   return nll;
// }
