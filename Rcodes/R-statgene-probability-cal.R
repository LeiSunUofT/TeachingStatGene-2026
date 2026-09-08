# Calculcate genotype probabilities

/* for trans. prob., cond. prob. and other prob.*/

/* cp2( ) is the conditional probability of P(X|2 IBD ) 
   where X is the genotype data for the pair, 
   similar for cp1( ) and cp0( ), unordered genotype for a pair */
double cp2(int f1, int m1, int f2, int m2, double* q)
{      
  double temp;
 
  if ( f1 == m1 && m1 == f2 && f2 == m2  ) 
    temp = q[f1]*q[f1];
  else if ( (f1 == f2 && m1 == m2) || (f1 == m2 && m1 == f2) ) 
    temp = 2*q[f1]*q[m1];
  else 
    temp = 0;

  return temp;
}

double cp1(int f1, int m1, int f2, int m2, double* q)
{
  double temp;

  if ( f1 == m1 && m1 == f2 && f2 == m2 ) 
    temp = q[f1]*q[f1]*q[f1];
  else if ( (f1 == f2 && m1 == m2) || (f1 == m2 && m1 == f2) ) 
    temp = q[f1]*q[m1]*(q[f1] + q[m1]);
  else if ( (f1 == m1 && m1 == f2) || (f1 == m1 && m1 == m2) ) 
    temp = 2*q[f1]*q[f2]*q[m2];
  else if ( (f2 == m2 && m2 == f1) || (f2 == m2 && m2 == m1) ) 
    temp = 2*q[f1]*q[m1]*q[f2];
  else if ( f1 == f2 || f1 == m2 ) 
    temp = 2*q[m1]*q[f2]*q[m2];
  else if ( m1 == f2 || m1 == m2 ) 
    temp = 2*q[f1]*q[f2]*q[m2];
  else 
    temp = 0;

  return temp;
}

double cp0(int f1, int m1, int f2, int m2, double* q)
{
  double temp;

  temp = q[f1]*q[m1]*q[f2]*q[m2];
  if ( f1 == m1 && m1==f2 && f2 == m2 ) 
    temp = temp;
  else if ( f1 == m1 &&  f2 == m2 ) 
    temp = 2*temp;
  else if ( f1 == m1 || f2 == m2 || (f1 == f2 && m1 == m2) 
	  || (f1 == m2 && m1 == f2) )
    temp = 4*temp;
  else 
    temp = 8*temp;

  return temp;
}
