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

  // Based Rec~S_hat functions. rectype options are 
  // 0 -- Linear 
  // 1 -- Linear with biomass ( y-1 )
  // 2 -- Beverton holt 
  // 3 -- Ricker
  // 4 -- DEFUNCT exponential
  
  //0.1 -- switches
  DATA_INTEGER( estMode);          // estimation mode?
  DATA_INTEGER( rectype );         // which type of recruitment model
  
  //0.2 -- RS data attributes
  DATA_VECTOR( years );            // recruiment years ( age class+1 if rec is in age 1 )( 1,nyrs,"years" );
  DATA_VECTOR( S_obs );            // Spawning biomass (or Number) in previous year  ( 1,nyrs,"S_hat" );  
  DATA_VECTOR( R_obs );            // Recruitment in a given year ( 1,nyrs,"Robs" );
  DATA_VECTOR( sdS );              // stdev of spawning biomass in previous year  ( 1,nyrs,"S_hat" );  
  DATA_VECTOR( sdR );              // stdev of recruitment in a given year ( 1,nyrs,"R" );
  
  //0.3 -- covariate attributes  
  DATA_MATRIX( rs_cov );           // environmental covariates        rs_cov.allocate( 1,ncov,1,nyrs,"rs_cov" );
  DATA_MATRIX( sdrs_cov );         // stdev environmental covariates  sdrs_cov.allocate( 1,ncov,1,nyrs,"sdrs_cov" );
  DATA_MATRIX( beta_0 );           // 0,1 matrix of inclusion         beta_0.allocate( 1,ncov,1,nyrs,"beta_0" );
  DATA_MATRIX( lambda_0 );         // binary 0,1 matrix of inclusion         lambda_0.allocate(1,ncov,1,nyrs,"lambda_0" );
  DATA_INTEGER( nyrs );            // number of years
  DATA_IVECTOR( cov_type );        // covariate type if 1 then run ration cov_type.allocate( 1,ncov,"cov_type" );
  DATA_INTEGER( ncov );            // number of environmental covariates
  
  //0.4 -- Observation error options
  DATA_INTEGER( sigMethod);         // which method for estimating sigma
  DATA_SCALAR( tau );               // scalar for observation and process error
  
  // ------------------------------------------------------------------------- //
  // 1. LOCAL PARAMETERS
  // ------------------------------------------------------------------------- //
  Type a       = Type(0.0);
  Type b       = Type(0.0);
  Type Acovars = Type(0.0);
  Type Bcovars = Type(0.0);
  Type k       = Type(0.0);                 // declare number of parms (minus sig)
  vector<Type>    S_hat( nyrs );            // declare local S
  vector<Type>    R_hat( nyrs );            // declare local R
  matrix<Type>    cov( ncov,nyrs );         // declare local cov
   
  // ------------------------------------------------------------------------- //
  // 2. INIT CALCS
  // ------------------------------------------------------------------------- //

  // ------------------------------------------------------------------------- //
  // 3. MODEL PARAMETERS
  // ------------------------------------------------------------------------- //
 
	PARAMETER( log_a ); 
	PARAMETER( log_b );         // b doesn't always have to be positive 
	PARAMETER_VECTOR( beta );   // covariate effects on pre-spawning success
	PARAMETER_VECTOR( lambda ); // covariate effects on post-spawning success
  PARAMETER( epsi_s  );       // error around SSB estimates
	PARAMETER( logsigma );      // process error
	PARAMETER( skipFit);        // used when projecting
	
	// ------------------------------------------------------------------------- //
	// 4. REC MODEL
	// ------------------------------------------------------------------------- //
	  
	 
    S_hat = S_obs * exp(epsi_s);  // annual random effect on SSB or deterministic estimated SBB error

    switch ( rectype ) {
      case 1: // Linear 
        a = exp( log_a ); // must be positive
        b = log_b;
        k = 1+ncov;
        for( int i= 0; i< nyrs; i++ ){
          Acovars = Type(0.0);  // Initialize
          Bcovars = Type(0.0);  // Initialize
         
          for ( int c= 0; c< ncov; c++ ){
            cov(c,i)     = rs_cov( c,i );
            SIMULATE {
              cov( c,i ) = rnorm( rs_cov( c,i ) , sdrs_cov( c,i )) ; // Simulate env
            }
            Acovars+= beta_0( c,i )*beta( c )*cov( c,i );
            Bcovars+= lambda_0( c,i )*lambda( c )*cov( c,i );
          }
          R_hat(i)  =  exp(a+Acovars);
        }
      break;
      case 2: // Linear with biomass ( y-1 )  BLM
        a = exp( log_a ); // must be positive
        b = exp( log_b );   
        k = 2+ncov;
        for( int i= 0; i< nyrs; i++ ){
          Acovars = Type(0.0);  // Initialize
          Bcovars = Type(0.0);  // Initialize
          
          for ( int c= 0; c< ncov; c++ ){
            cov(c,i)     = rs_cov( c,i );
            SIMULATE {
              cov( c,i ) = rnorm( rs_cov( c,i ) , sdrs_cov( c,i )) ; // Simulate env
            }
            Acovars+= beta_0( c,i )*beta( c )*cov( c,i );
            Bcovars+= lambda_0( c,i )*lambda( c )*cov( c,i );
          }
          R_hat(i)  =  exp(a+Acovars+(b+Bcovars)*log(S_hat(i)));
      }
      break;
      case 3:  // Beverton Holt
        a = exp( log_a );// must be positive
        b = exp( log_b );// must be positive
        k = 2+ncov;
        for( int i= 0; i< nyrs; i++ ){
          Acovars = Type(0.0);  // Initialize
          Bcovars = Type(0.0);  // Initialize
          
          for ( int c= 0; c< ncov; c++ ){
            cov(c,i)     = rs_cov( c,i );
            SIMULATE {
              cov( c,i ) = rnorm( rs_cov( c,i ) , sdrs_cov( c,i )) ; // Simulate env
            }
            Acovars+= beta_0( c,i )*beta( c )*cov( c,i );
            Bcovars+= lambda_0( c,i )*lambda( c )*cov( c,i );
          }
          R_hat(i)  =  (a*S_hat(i)*exp(Acovars)) / (1+ b*(S_hat(i)*exp(Bcovars)) );
        }
      break;
    case 4: // Ricker
      a = exp( log_a );
      b = exp( log_b );
      k = 2+ncov;
      for( int i= 0; i< nyrs; i++ ){
        
        Acovars = Type(0.0);  // Initialize
        Bcovars = Type(0.0);  // Initialize
        
        for ( int c= 0; c< ncov; c++ ){
          cov(c,i)     = rs_cov( c,i );
          SIMULATE {
            cov( c,i ) = rnorm( rs_cov( c,i ) , sdrs_cov( c,i )) ; // Simulate env
          }
          Acovars+= beta_0( c,i )*beta( c )*cov( c,i );
          Bcovars+= lambda_0( c,i )*lambda( c )*cov( c,i );
        }
        //modified from mueter 2011 formulation
        R_hat(i)= exp( (a+Acovars)-( (b+Bcovars)*S_hat(i))+log(S_hat(i))  ) ;
        
        //modified from mueter 2011 formulation old way - not right?
        //R_hat(i)= exp((a*exp(Acovars))-(b*exp(Bcovars)*S_hat(i))+log(S_hat(i))) ;
        
      }
      break;
    // case 5: // exponential
    //   a     = exp( log_a );
    //   b     = -.99*exp( log_a-log_b );
    //   k     =  2+ncov;
    //   Acovars = Type(0.0);  // Initialize
    //   Bcovars = Type(0.0);  // Initialize
    //   for( int i= 0; i< nyrs; i++ ){
    //     for ( int c= 0; c< ncov; c++ ){
    // Acovars+= beta_0( c,i )*beta( c )*rs_cov( c,i );
    // Bcovars+= lambda_0( c,i )*lambda( c )*rs_cov( c,i );
    //     }
    //     R_hat(i)  =  exp(a+b*S_hat(i)+exp(Acovars));
    //   }
    //   
    //   
    break;
    default:
      error( "Invalid 'rectype'" );
    }
    
    vector<Type> ols_R     = (log(R_obs) - log(R_hat));
    vector<Type> ols_S     = (log(S_obs) - log(S_hat));
     
  // ------------------------------------------------------------------------- //
  // 5. CALC OBJ FUN
  // ------------------------------------------------------------------------- //
  
    // 5.1 calculate obj function based on sigMethod
    
    // formulation based on: Porch CE, Lauretta MV (2016) On Making Statistical 
    // Inferences Regarding the Relationship  between Spawners and Recruits and the Irresolute
    // Case of Western Atlantic Bluefin Tuna (Thunnus thynnus). PLoS ONE 11(6): e0156767. 
    // doi:10.1371/journal.pone.0156767

    // pre-declare variables
    Type sig    = Type(0.0);
    Type sig_R  = Type(0.0);
    Type sig_S  = Type(0.0);
    Type mn_sdR = Type(0.0);
    Type mn_sdS = Type(0.0);
    
    switch ( sigMethod ) {
    case 1: //  no observation error tau = 0
        sig   =  exp(logsigma);   // -- process error of recruitment
        sig_R =  sig;             // -- observation error of recruitment
        sig_S =  Type(0.00001);   // -- Not used bc epsi_S  = 0
        
      break;
      case 2: // estimate sigma, random effects on SSB if tau >0
        // observation errors sig_R and sig_S are statistically independent, random normal variables with similar variance sig_R
        // but because of lack of other info sig_R is a function of process error via the tau scalar see eq 4 in 
        // Porch, C. E., and M. V. Lauretta. 2016. 
        
         sig   = exp(logsigma);          // -- process error of recruitment
         sig_R = (Type(1.0)+tau)*sig;    // -- observation error of recruitment
         sig_S = tau*(sig);              // -- observation error of spawners (can be same as sig_R)
      break;
      case 3:// as in case 1 but using an unbiased sigma (rather than biased estimate from MLE; sensu Ludwig and Walters 1982):
         sig = Type(0.0);
         for( int i = 0; i< nyrs; i++ )
            sig +=  ( (ols_R[i]*ols_R[i])/(1+tau)) + ((ols_S[i]*ols_S[i])/(tau) );
         sig = (1/(nyrs-k))*sig;
         sig_R = (Type(1.0)+tau)*sig;
         sig_S = tau*(sig);
        break;
      case 4:
        // as in 1 but with defined measurement error for rec (indep of random effects on S_hat)
        // tau = 0 defaults to no random effects on Spawners
         sig  = exp(logsigma);
         for( int i= 0; i< nyrs; i++ )
           mn_sdR += sdR(i)/nyrs;
         sig_R = sig + mn_sdR;
         sig_S = tau*(sig);
      break;
      case 5:// as in 1 but with defined measurement error for rec and SSB; tau = 0 defaults to no random effects on S, and sig for R
        sig  = exp(logsigma);
        for( int i = 0; i< nyrs; i++ ){
          mn_sdR += sdR(i)/nyrs;
          mn_sdS += sdS(i)/nyrs;
        }
        sig_R =  ((Type(1.0)+tau)*sig) + mn_sdR;
        sig_S =  (tau*(sig)) + mn_sdS;
      break;
      default:
        error( "Invalid 'sigMethod'" );
    }
   
    
    
    // ------------------------------------------------------------------------- //
    // 6. END
    // ------------------------------------------------------------------------- //

    //calculate negative log-likelihood
    Type nll_test  = Type(0.0);
    for( int i = 0; i< nyrs; i++ )
      nll_test += ((ols_R[i]/sig_R)*(ols_R[i]/sig_R)) + ((ols_S[i]/sig_S)*(ols_S[i]/sig_S));
    nll_test =  0.5*(  (nyrs*log(2*3.141593*sig_S*sig_S)) + (nyrs*log(2*3.141593*sig_R*sig_R)) + nll_test);


    Type nll = -sum(dnorm(log(R_obs),log(R_hat),sig_R,true));
        nll += -sum(dnorm(log(S_obs),log(S_hat),sig_S,true));
     
    
    if(estMode == 0)  nll  = skipFit;
  
    
    SIMULATE {
      //for( int i= 0; i< nyrs; i++ ){
        R_obs = rnorm( log(R_hat) ,sig_R) ; // Simulate response
        R_obs = exp(R_obs); // Simulate response
        S_obs = rnorm( log(S_hat) ,sig_S) ; // Simulate response
        S_obs = exp(S_obs); // Simulate response
      //}
       REPORT(R_obs);            // Report the simulation
       REPORT(S_obs);            // Report the simulation
    }
    
    
  // ------------------------------------------------------------------------- //
  // 7. Report
  // ------------------------------------------------------------------------- // 
  
  REPORT( nll );
  REPORT( nll_test );
  REPORT( S_hat );
  REPORT( S_obs );
  REPORT( R_hat );
  REPORT( R_obs );
  REPORT( cov );
  REPORT( sig_S );
  REPORT( sig_R );
  REPORT( sig );
  REPORT( tau );
  REPORT( log_a );
  REPORT( a );
  REPORT( log_b );
  REPORT( b );
  REPORT (logsigma);
  REPORT( ols_R );
  REPORT( ols_S );
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
  ADREPORT( log_a );
  ADREPORT( a );
  ADREPORT( log_b );
  ADREPORT( b );
  ADREPORT (logsigma);
  ADREPORT( ols_R );
  ADREPORT( ols_S ); 
  ADREPORT( beta );
  ADREPORT( lambda );
  
  return nll;
  
  
}
