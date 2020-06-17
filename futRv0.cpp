#include <TMB.hpp>

  // ------------------------------------------------------------------------- //
  //                    futR version 1.0.1                                     //
  //       Climate enhanced recruitment for projections                        //
  //           Under climate and environmental change                          //
  //                         ACLIM                                             //
  //                                                                           //
  // AUTHORS:   K Holsman, G Adams, K Aydin, J Ianelli, A Punt                 //
  // CITATIONS:                                                                //
  // 1. Holsman, et al. 2019.                                                  //
  // ------------------------------------------------------------------------- //
  // 
  //  INDEX (following the impeccable G. Adams design)
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
                                    // 2 ) linear 
                                    // 3 ) Beverton holt 
                                    // 4 ) Ricker
  DATA_VECTOR( years );             // recruiment years                   ( age class+1 if rec is in age 1 )( 1,nyrs,"years" );
  DATA_VECTOR( SSB );               // Spawning biomass in previous year  ( 1,nyrs,"S_hat" );  
  DATA_VECTOR( Robs );              // Recruitment in a given year        ( 1,nyrs,"Robs" );
  DATA_VECTOR( sdS );               // Spawning biomass in previous year  ( 1,nyrs,"S_hat" );  
  DATA_VECTOR( sdR );               // Recruitment in a given year        ( 1,nyrs,"R" );
  DATA_VECTOR( Ration_scaled );     // Predator ration in a given year    ( 1,nyrs,"Ration" );
  DATA_MATRIX( rs_cov );            // environmental covariates           rs_cov.allocate( 1,ncov,1,nyrs,"rs_cov" );
  DATA_INTEGER( nyrs );
  DATA_IVECTOR( cov_type );         // covariate type if 1 then run ration cov_type.allocate( 1,ncov,"cov_type" );
  DATA_INTEGER( ncov );
  DATA_INTEGER( sigMethod);     // apply bias correction to sigma? 0 = NO; 1 = yes
  DATA_SCALAR( tau );
  
  // ------------------------------------------------------------------------- //
  // 1. LOCAL PARAMETERS
  // ------------------------------------------------------------------------- //

  // int nyrs  =   years.size(  );         // integer of years
  // int ncov  =   rs_cov.size(  );        // integer of rs_cov
  
  vector<Type> logR_hat( nyrs );            // declare local R
  vector<Type> R_hat( nyrs );               // declare local Rhat
  Type k           = 0.0;                   // declare number of parms (minus sig)

  
  // ------------------------------------------------------------------------- //
  // 2. INIT CALCS
  // ------------------------------------------------------------------------- //
  
  
  // ------------------------------------------------------------------------- //
  // 3. MODEL PARAMETERS
  // ------------------------------------------------------------------------- //
  
	PARAMETER( log_a ); 
	PARAMETER( log_b ); 
  //PARAMETER( logit_tau );
	PARAMETER_VECTOR( rs_parm );  
	PARAMETER( epsi_s  );  // vector of error around SSB estimates
  //PARAMETER_VECTOR(epsi_r);  // vector of error around Rec estimates
	PARAMETER( logsigma ); 
	
	// ------------------------------------------------------------------------- //
	// 4. REC MODEL
	// ------------------------------------------------------------------------- //
    
    // 4.1 Set up phase values
    Type a   = Type(0.0);
    Type b   = Type(0.0);
    Type covars = Type(0.0);
    vector<Type>    S_hat = SSB  * exp(epsi_s);  // annual random effect on SSB or deterministic estimated SBB error
    //vector<Type>    R = Robs * exp(epsi_r);    // annual random effect on Rec or deterministic estimated Rec error
   
    // 4.2 calculate recruitment    
    switch ( rectype ) {
    case 1: // Linear with biomass ( y-1 )  BLM
      a = ( log_a );
      b = ( log_b );
      for( int i= 0; i< nyrs; i++ )
      {
        covars = Type(  0  );  // Initialize
        for ( int c= 0; c< ncov; c++ )
        {
          if( cov_type( c )==1 ){
            covars+=rs_parm( c )*( rs_cov( c,i )-Ration_scaled( i ) );
          }else{
            covars+=rs_parm( c )*rs_cov( c,i );
          }
        }
        logR_hat( i )= ( a+b*log( S_hat( i ) )+covars );
      } 
      k = 2+ncov;
      break;
      
    case 2: // Linear LM
      a = ( log_a );
      for( int i= 0; i< nyrs; i++ )
      {
        covars = Type(  0  );  // Initialize
        for ( int c= 0; c< ncov; c++ )
        {
          if( cov_type( c )==1 ){
            covars+=rs_parm( c )*( rs_cov( c,i )-Ration_scaled( i ) );
          }else{
            covars+=rs_parm( c )*rs_cov( c,i );
          }
        }
        logR_hat( i )= ( a+covars );
      } 
      k = 1+ncov;
      break;
      
    case 3: // Beverton Holt 
      a = exp( log_a );
      b = exp( log_b );
      R_hat  = 0.0;
      
      for( int i= 0; i< nyrs; i++ ){
        covars = Type(  0  );  // Initialize
        for ( int c= 0; c< ncov; c++ )
        {
          if( cov_type( c )==1 ){
            covars+=rs_parm( c )*( rs_cov( c,i )-Ration_scaled( i ) );
          }else{
            covars+=rs_parm( c )*rs_cov( c,i );
          }
        }

        R_hat( i )=( a* S_hat( i )/( Type(1.0) + b * S_hat( i ) ) ) + covars;
        
      }
      logR_hat=log( R_hat );
      k = 2+ncov;
      break;
      
    case 4: // Ricker
      a=exp( log_a );
      b=exp( log_b );
      for( int i= 0; i< nyrs; i++ )
      {
        covars = Type(  0  );  // Initialize
        for ( int c= 0; c< ncov; c++ )
        {
          if( cov_type( c ) == 1 ){
            covars+=rs_parm( c )*( rs_cov( c,i )-Ration_scaled( i ) );
          }else{
            covars+=rs_parm( c )*rs_cov( c,i );
          }
        }
        // logR_hat( i )= ( log( a*S_hat( i ) ) -b*S_hat( i )+covars ); // may 2017
        logR_hat( i )= a-b*S_hat( i )+covars+log( S_hat( i ) );  //a+b*S_hat+sum( covs*betas )+log( S_hat )  # mueter 2011 formulation
      }	
      k = 2+ncov;
      break;
    default:
      error( "Invalid 'rectype'" );
    }

    R_hat = exp( logR_hat );
     
  // ------------------------------------------------------------------------- //
  // 5. CALC OBJ FUN
  // ------------------------------------------------------------------------- //
  
    // 5.1 calculate obj function based on sigMethod
    
    // formulation based on: Porch CE, Lauretta MV (2016) On Making Statistical 
    // Inferences Regarding the Relationship  between Spawners and Recruits and the Irresolute
    // Case of Western Atlantic Bluefin Tuna (Thunnus thynnus). PLoS ONE 11(6): e0156767. 
    // doi:10.1371/journal.pone.0156767
    
    vector<Type> epsi    = log(Robs) - logR_hat;
    vector<Type> epsiSSB = log(SSB) - log(S_hat);
    Type sig  = exp(logsigma);
    Type sig_R = Type(0.0);
    Type sig_S = Type(0.0);
    Type nll = Type(0.0);
    switch ( sigMethod ) {
      case 1: // estimate sigma, random effects on SSB if tau >0
         sig  = exp(logsigma);
         sig_R = (1+tau)*sig;
         sig_S = tau*(sig);
         nll -= sum(dnorm(epsi,Type(0.0),sig_R,true));
         nll -= sum(dnorm(epsiSSB,Type(0.0),sig_S,true));
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
      case 2:// unbiased sigma (rather than biased estimate from MLE sensu Ludwig and Walters):
        sig = 0.0;
        for( int i= 0; i< nyrs; i++ )
         sig =  ( ((epsi[i]*epsi[i])/(1+tau)) + ((epsiSSB[i]*epsiSSB[i])/(tau)) );
         sig = (1/(nyrs-k))*sig;
         sig_R = (1+tau)*sig;
         sig_S = tau*(sig);
         nll -= sum(dnorm(epsi,Type(0.0),sig_R,true));
         nll -= sum(dnorm(epsiSSB,Type(0.0),sig_S,true));
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
      case 3:// as in 1 but with defined measurement error for rec (indep of random effects on S_hat)
        sig  = exp(logsigma);
         sig_S = tau*(sig);
         nll -= sum(dnorm(epsi,Type(0.0),sig + sdR,true));
         nll -= sum(dnorm(epsiSSB,Type(0.0),sig_S,true));
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
      case 4:// as in 1 but with defined measurement error for rec and SSB; tau = 0 defaults to no random effects
        sig  = exp(logsigma);
        sig_R =  ((1+tau)*sig);
        sig_S =  (tau*(sig));
        nll -= sum(dnorm(epsi,Type(0.0),sig_R + sdR,true));
        nll -= sum(dnorm(epsiSSB,Type(0.0),sig_S + sdS,true));
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
 // REPORT( SSB );
 // REPORT( Robs );
  REPORT( R_hat );
  REPORT( log_a );
  REPORT( log_b );
  REPORT (logsigma);
  REPORT( epsi );
  REPORT( rs_parm );
  REPORT( tau );
  REPORT( epsi_s );
  REPORT( epsiSSB );

  ADREPORT( S_hat );
  ADREPORT( R_hat );
  ADREPORT( log_a );
  ADREPORT( log_b );
  ADREPORT (logsigma);
  ADREPORT( epsi );
  ADREPORT( rs_parm );
  ADREPORT( epsi_s );
  ADREPORT( epsiSSB );  
  
  // ------------------------------------------------------------------------- //
  // 7. END
  // ------------------------------------------------------------------------- //
  
  return nll;
  
  
}
