#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
	DATA_VECTOR(S);
	DATA_VECTOR(R);
  DATA_VECTOR(Temp);
  int n = S.size();  // integer of y
  // Type ohno =R(100);

	PARAMETER(loga);
	PARAMETER(b);
  PARAMETER(c);
	PARAMETER(logSigma);
	vector<Type> logRhat(n);  // declare local integer
  vector<Type> Rhat(n);  // declare local integer

  Type neglogL = 0.0;
  Type a=exp(loga);

  logRhat=log(a*S)-b*S+c*Temp; //+ohno;
  Rhat=exp(logRhat);

  neglogL = -sum(dnorm(log(R), logRhat, exp(logSigma), true));
  REPORT(Rhat);
  ADREPORT(Rhat);
  ADREPORT(exp(loga));
  REPORT(loga);
  ADREPORT((b));
  REPORT(b);
  ADREPORT((c));
  REPORT(c);
  return neglogL;

}
