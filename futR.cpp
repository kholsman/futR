#include <TMB.hpp>

  // ------------------------------------------------------------------------- //
  //                    futR version 1.0.1                                     //
  //       Climate enhanced recruitment for projections                        //
  //           Under climate and environmental change                          //
  //                         ACLIM                                             //
  //                                                                           //
  // AUTHORS:   K Holsman, Aydin, Adams, Punt et al.                           //
  // CITATIONS:                                                                //
  // 1. Holsman, et al. 2020 Climate and trophic controls on groundfish        //
  // recruitment in Alaska.                                                    //
  // ------------------------------------------------------------------------- //
  // 
  //  INDEX (following the fab G. Adams design)
  //  0. Load data
  //  1. Local parameters
  //  2. Estimated parameters
  //  3. Init. calcs
  //  4. Recruitment model
  //  5. Calc obj fun
  //  6. Report section
  //  7. End
  
template<class Type>
Type objective_function<Type>::operator(  ) (  )
{
  // ------------------------------------------------------------------------- //
  // 0. LOAD DATA
  // ------------------------------------------------------------------------- //
  // DATA_ARRAY(), DATA_FACTOR(), DATA_IARRAY(), DATA_IMATRIX(), DATA_INTEGER(), 
  // DATA_IVECTOR(), DATA_MATRIX(), DATA_SCALAR(), DATA_SPARSE_MATRIX(), 
  // DATA_STRING(), DATA_STRUCT(), DATA_UPDATE(), DATA_VECTOR()
  DATA_INTEGER( rectype );         // Based Rec~S_hat functions. Options are 
                                    // 1 ) Linear with biomass ( y-1 )
                                    // 2 ) Beverton holt 
                                    // 3 ) Ricker
                                    // 4 ) exponential
  DATA_VECTOR(  years );            // recruiment years                   ( age class+1 if rec is in age 1 )( 1,nyrs,"years" );
  DATA_VECTOR(  S_obs );            // Spawning biomass (or Number) in previous year  ( 1,nyrs,"S_hat" );  
  DATA_VECTOR(  R_obs );             // Recruitment in a given year        ( 1,nyrs,"Robs" );
  DATA_VECTOR(  sdS );              // stdev of spawning biomass in previous year  ( 1,nyrs,"S_hat" );  
  DATA_VECTOR(  sdR );              // stdev of recruitment in a given year        ( 1,nyrs,"R" );
  //DATA_VECTOR(  Ration_scaled );    // Predator ration in a given year    ( 1,nyrs,"Ration" );
  DATA_MATRIX(  rs_cov );           // environmental covariates           rs_cov.allocate( 1,ncov,1,nyrs,"rs_cov" );
  DATA_VECTOR(  sdrs_cov );          // stdev environmental covariates     sdrs_cov.allocate( 1,ncov,1,nyrs,"rs_cov" );
  DATA_INTEGER( nyrs );
  DATA_IVECTOR( cov_type );         // covariate type if 1 then run ration cov_type.allocate( 1,ncov,"cov_type" );
  DATA_INTEGER( ncov );
  DATA_INTEGER( sigMethod);         // apply bias correction to sigma? 0 = NO; 1 = yes
  DATA_INTEGER( tMethod);           // apply 0-1 transformation to gamma for ricker model; 1 = complementary loglog (default); 2 = logit
  DATA_SCALAR(  tau );
  
  // ------------------------------------------------------------------------- //
  // 1. LOCAL PARAMETERS
  // ------------------------------------------------------------------------- //

  // int nyrs  =   years.size(  );         // integer of years
  // int ncov  =   rs_cov.size(  );        // integer of rs_cov
  // 4.1 Set up phase values
  Type a       = Type(0.0);
  Type b       = Type(0.0);
  Type g       = Type(0.0);
  Type Acovars = Type(0.0);
  Type Bcovars = Type(0.0);
  Type k       = Type(0.0);                // declare number of parms (minus sig)
  
  //vector<Type>    R = R_obs * exp(epsi_r);    // annual random effect on Rec or deterministic estimated Rec error
  vector<Type>    S_hat( nyrs );            // declare local S
  vector<Type>    R_hat( nyrs );            // declare local R
   
  // ------------------------------------------------------------------------- //
  // 2. INIT CALCS
  // ------------------------------------------------------------------------- //
  
  
  // ------------------------------------------------------------------------- //
  // 3. MODEL PARAMETERS
  // ------------------------------------------------------------------------- //
  
	PARAMETER( log_a ); 
	PARAMETER( log_b );   // b doesn't always have to be positive 
	PARAMETER( gamma ); 
	PARAMETER_VECTOR( beta ); 
	PARAMETER_VECTOR( lambda ); 
  //PARAMETER( logit_tau );
  PARAMETER( epsi_s  );  // vector of error around SSB estimates
  //PARAMETER_VECTOR( epsi_S  );  // vector of error around Spawner estimates
  //PARAMETER_VECTOR( epsi_R  );  // vector of error around Rec estimates
	PARAMETER( logsigma ); 
	
	// ------------------------------------------------------------------------- //
	// 4. REC MODEL
	// ------------------------------------------------------------------------- //
    
    
    // 4.1. Precalcs:
    // 4.3 calculate recruitment:
    // 1. Linear (gamma = 0)
    // 2. Beverton Holt (gamma = -1)
    // 3. Ricker (0 < gamma <1 )
    // 4. Exponential (gamma=1, b<0)
    switch ( rectype ) {
      case 1: // Linear with biomass ( y-1 )  BLM
        a = exp( log_a );
        b = log_b;  // isn't used in the model
        g = 0;
        k = 1+ncov;
      break;
      case 2:  // Beverton Holt
        a = exp( log_a );
        b = exp( log_b );
        g = -1;
        k = 2+ncov;
        
      break;
    case 3: // Ricker
      a = exp( log_a );
      b = exp( log_b );
      switch ( tMethod ){
        case 1:
          // gamma = log(-log(1-g));
          g  = 1-exp(-exp(gamma)) ;
        break;
        case 2:
          //gamma = log((g/(1-g)));
          g = exp(gamma)/(1+exp(gamma));
        break;
        default:
          error( "Invalid 'tMethod'" );
      }
      k = 3 + ncov;
      break;
    case 4: // exponential
      a     = exp( log_a );
      g     =   1;
      b     =   -1 * exp( log_b );
      k     =   2+ncov;
      
    break;
    default:
      error( "Invalid 'rectype'" );
    }
    
   // R_hat = 0.0;
   S_hat = S_obs  * exp(epsi_s);  // annual random effect on SSB or deterministic estimated SBB error
    for( int i= 0; i< nyrs; i++ )
    {
      Acovars = Type(0);  // Initialize
      Bcovars = Type(0);  // Initialize
      for ( int c= 0; c< ncov; c++ )
      {
        Acovars+= beta( c )*rs_cov( c,i );
        Bcovars+= lambda( c )*rs_cov( c,i );
      }

     // R_hat(i)  =  a*S_hat(i)*exp(Acovars)*pow( Type(1.0) - (b*g*S_hat(i)*exp(Acovars)) ,Type(1.0)/g  )*exp(Bcovars);
      R_hat(i)  =  a*S_obs(i)*exp(Acovars)*pow( Type(1.0) - (b*g*S_obs(i)*exp(Acovars)) ,Type(1.0)/g  )*exp(Bcovars);
      
    }
    
     
  // ------------------------------------------------------------------------- //
  // 5. CALC OBJ FUN
  // ------------------------------------------------------------------------- //
  
    // 5.1 calculate obj function based on sigMethod
    
    // formulation based on: Porch CE, Lauretta MV (2016) On Making Statistical 
    // Inferences Regarding the Relationship  between Spawners and Recruits and the Irresolute
    // Case of Western Atlantic Bluefin Tuna (Thunnus thynnus). PLoS ONE 11(6): e0156767. 
    // doi:10.1371/journal.pone.0156767
    
    vector<Type> epsi_R     = log(R_obs) - log(R_hat);
    vector<Type> epsi_S     = log(S_obs) - log(S_hat);
   // matrix<Type> epsi_cov   = rs_cov     - rs_cov_hat;
   // nll -= sum(dnorm(epsi_cov,Type(0.0),rs_cov_sd,true));
    
    Type sig   = exp(logsigma);
    Type sig_R = Type(0.0);
    Type sig_S = Type(0.0);
    Type nll   = Type(0.0);
    
    switch ( sigMethod ) {
      case 1: // estimate sigma, random effects on SSB if tau >0
        // observation errors sig_R and sig_S are statistically independent, random normal variables with similar variance sig_R
        // but because of lack of other info sig_R is a function of process error via the tau scalar see eq 4 in 
        // Porch, C. E., and M. V. Lauretta. 2016. 
         sig  = exp(logsigma);  // process error of recruitment
         sig_R = (1+tau)*sig;   // observation error of recruitment
         sig_S = tau*(sig);     // observation error of spawners (can be same as sig_R)
         nll -= sum(dnorm(epsi_R,Type(0.0),sig_R,true));
         nll -= sum(dnorm(epsi_S,Type(0.0),sig_S,true));
         SIMULATE {
           
           for( int i= 0; i< nyrs; i++ ){
            R_hat(i) = exp(rnorm(Type(0.0),sig_R)+log(R_hat(i)));  // Simulate response
            S_hat(i) = exp(rnorm(Type(0.0),sig_S)+log(S_hat(i)));    // Simulate response
           }
                         REPORT(R_hat);          // Report the simulation
                         REPORT(S_hat);          // Report the simulation
                         ADREPORT(R_hat);          // Report the simulation
                         ADREPORT(S_hat);          // Report the simulation
         }
      break;
      case 2:// as in case 1 but using an unbiased sigma (rather than biased estimate from MLE; sensu Ludwig and Walters 1982):
        sig = 0.0;
        for( int i= 0; i< nyrs; i++ )
         sig =  ( ((epsi_R[i]*epsi_R[i])/(1+tau)) + ((epsi_S[i]*epsi_S[i])/(tau)) );
         sig = (1/(nyrs-k))*sig;
         sig_R = (1+tau)*sig;
         sig_S = tau*(sig);
         nll -= sum(dnorm(epsi_R,Type(0.0),sig_R,true));
         nll -= sum(dnorm(epsi_S,Type(0.0),sig_S,true));
         SIMULATE {
           for( int i= 0; i< nyrs; i++ ){
            R_hat(i) = exp(rnorm(Type(0.0),sig_R)+log(R_hat(i)));    // Simulate response
            S_hat(i) = exp(rnorm(Type(0.0),sig_S)+log(S_hat(i)));    // Simulate response
           }
           REPORT(R_hat);          // Report the simulation
           REPORT(S_hat);          // Report the simulation
           ADREPORT(R_hat);          // Report the simulation
           ADREPORT(S_hat);          // Report the simulation
         }
        break;
      case 3:
        // as in 1 but with defined measurement error for rec (indep of random effects on S_hat)
        // tau = 0 defaults to no random effects on Spawners
         sig  = exp(logsigma);
         sig_S = tau*(sig);
         nll -= sum(dnorm(epsi_R,Type(0.0),sig + sdR,true));
         nll -= sum(dnorm(epsi_S,Type(0.0),sig_S,true));
         SIMULATE {
           for( int i= 0; i< nyrs; i++ ){
            R_hat(i) = exp(rnorm(Type(0.0),sig + sdR(i))+log(R_hat(i)));  // Simulate response
            S_hat(i) = exp(rnorm(Type(0.0),sig_S)+log(S_hat(i)));    // Simulate response
           }
           REPORT(R_hat);          // Report the simulation
           REPORT(S_hat);          // Report the simulation
           ADREPORT(R_hat);          // Report the simulation
           ADREPORT(S_hat);          // Report the simulation
         }
      break;
      case 4:// as in 1 but with defined measurement error for rec and SSB; tau = 0 defaults to no random effects on S, and sig for R
        sig  = exp(logsigma);
        sig_R =  ((1+tau)*sig);
        sig_S =  (tau*(sig));
        nll -= sum(dnorm(epsi_R,Type(0.0),sig_R + sdR,true));
        nll -= sum(dnorm(epsi_S,Type(0.0),sig_S + sdS,true));
        SIMULATE {
          for( int i= 0; i< nyrs; i++ ){
            R_hat(i) = exp(rnorm(Type(0.0),sig_R + sdR(i))+log(R_hat(i)));  // Simulate response
            S_hat(i) = exp(rnorm(Type(0.0),sig_S + sdS(i))+log(S_hat(i)));    // Simulate response
          }
          REPORT(R_hat);          // Report the simulation
          REPORT(S_hat);          // Report the simulation
          ADREPORT(R_hat);          // Report the simulation
          ADREPORT(S_hat);          // Report the simulation
        }
      break;
      default:
        error( "Invalid 'sigMethod'" );
    }

    
    
  // ------------------------------------------------------------------------- //
  // 6. Report
  // ------------------------------------------------------------------------- // 

  REPORT( S_hat );
  REPORT( S_obs );
  REPORT( R_hat );
  REPORT( R_obs );
  REPORT( sig_S );
  REPORT( sig_R );
  REPORT( sig );
  REPORT( tau );
  REPORT (gamma);
  REPORT(g);
  REPORT( log_a );
  REPORT( a );
  REPORT( log_b );
  REPORT( b );
  REPORT (logsigma);
  REPORT( epsi_R );
  REPORT( epsi_S );

  REPORT( beta );
  REPORT( lambda );


  ADREPORT( S_hat );
  ADREPORT( S_obs );
  ADREPORT( R_hat );
  ADREPORT( R_obs );
  ADREPORT( sig_S );
  ADREPORT( sig_R );
  ADREPORT( sig );
  ADREPORT( tau );
  ADREPORT (gamma);
  ADREPORT(g);
  ADREPORT( log_a );
  ADREPORT( a );
  ADREPORT( log_b );
  ADREPORT( b );
  ADREPORT (logsigma);
  ADREPORT( epsi_R );
  ADREPORT( epsi_S ); 
  ADREPORT( beta );
  ADREPORT( lambda );
  
  // ------------------------------------------------------------------------- //
  // 7. END
  // ------------------------------------------------------------------------- //
  
  return nll;
  
  
}
