#ifndef COMMON_H
#define COMMON_H

#ifdef _WIN32
#define Subdirectory "TKS\\"
#else  // __linux
#define Subdirectory "TKS//"
#endif

#define PREP_old(Fpre, Fnum, R_W)             \
  {                                           \
    FILE *fp;                                 \
    {                                         \
      char PRT[100];                          \
      sprintf(PRT, "%s%dsA.dat", Fpre, Fnum); \
      if ((fp = fopen(PRT, R_W)) == NULL) {   \
        printf("Fale File=%s.\n", PRT);       \
        exit(0);                              \
      }                                       \
    }

#define PREP_check(Fpre, Fnum, R_W)           \
  {                                           \
    FILE *fp;                                 \
    {                                         \
      char PRT[100];                          \
      sprintf(PRT, "%s%dsA.dat", Fpre, Fnum); \
      if ((fp = fopen(PRT, R_W)) != NULL) {   \
        nowtime = Fnum;                       \
        yesfile++;                            \
        { fclose(fp); }                       \
      }                                       \
    }                                         \
  }

#define CLSP      \
  { fclose(fp); } \
  }

#define NR_END 1
#define FREE_ARG char *

double *double_vector(long nl, long nh) {
  double *v;
  v = (double *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));

  if (!v) {
    printf("it is just wrong to allocate");
    exit(0);
  };

  return v - nl + NR_END;
}

void free_double_vector(double *v, long nl, long nh) {
  free((FREE_ARG)(v + nl - NR_END));
}

void tridag(double a[], double b[], double c[], double r[], double u[],
            unsigned long n) {
  unsigned long j;
  double bet, *gam;
  gam = double_vector(1, n); /* One vector of workspace, gam is needed.*/

  if (b[1] == 0.0) printf("Error 1 in tridag");

  /*If this happens then you should rewrite your equations as a set of order
  N-1, with u2 trivially eliminated.*/
  u[1] = r[1] / (bet = b[1]);

  for (j = 2; j <= n; j++) /*Decomposition and forward substitution.*/
  {
    gam[j] = c[j - 1] / bet;
    bet = b[j] - a[j] * gam[j];

    if (bet == 0.0) printf("Error 2 in tridag"); /* Algorithm fails; see below*/

    u[j] = (r[j] - a[j] * u[j - 1]) / bet;
  }

  for (j = (n - 1); j >= 1; j--)
    u[j] -= gam[j + 1] * u[j + 1]; /*Backsubstitution.*/

  free_double_vector(gam, 1, n);
}

void cyclic(double a[], double b[], double c[], double alpha, double beta,
            double r[], double x[], unsigned long n) {
  unsigned long i;
  double fact, gamma, *bb, *u, *z;

  if (n <= 2) printf("n too small in cyclic");

  bb = double_vector(1, n);
  u = double_vector(1, n);
  z = double_vector(1, n);
  gamma = -b[1];        /* Avoid subtraction error in forming bb[1].*/
  bb[1] = b[1] - gamma; /*Set up the diagonal of the modied tridiagonalsystem*/
  bb[n] = b[n] - alpha * beta / gamma;

  for (i = 2; i < n; i++) bb[i] = b[i];

  tridag(a, bb, c, r, x, n); /* Solve A . x = r.*/
  u[1] = gamma;              /*Set up the vector u.*/
  u[n] = alpha;

  for (i = 2; i < n; i++) u[i] = 0.0;

  tridag(a, bb, c, u, z, n); /*Solve A . z = u.*/
  fact = (x[1] + beta * x[n] / gamma) /
         (1.0 + z[1] + beta * z[n] / gamma); /*Form v . x=(1 + v . z).*/

  for (i = 1; i <= n; i++)
    x[i] -= fact * z[i]; /*Now get the solution vector x.*/

  free_double_vector(z, 1, n);
  free_double_vector(u, 1, n);
  free_double_vector(bb, 1, n);
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

float ran2(long *idum) {
  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) {
    if (-(*idum) < 1)
      *idum = 1;
    else
      *idum = -(*idum);
    idum2 = (*idum);
    for (j = NTAB + 7; j >= 0; j--) {
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k * IQ1) - k * IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = (*idum) / IQ1;
  *idum = IA1 * (*idum - k * IQ1) - k * IR1;
  if (*idum < 0) *idum += IM1;
  k = idum2 / IQ2;
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
  if (idum2 < 0) idum2 += IM2;
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp = AM * iy) > RNMX)
    return RNMX;
  else
    return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

double gasdev(long *idum) {
  // Returns a normally distributed deviate with zero mean and unit variance,
  // using ran1(idum)  as the source of uniform deviates.
  float ran2(long *idum);
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
  if (iset == 0)  // We don't have an extra deviate handy, so
  {
    do {
      v1 = 2.0 * ran2(idum) - 1.0;  // pick two uniform numbers in the square
                                    // extendingfrom -1 to +1 in each direction,
      v2 = 2.0 * ran2(idum) - 1.0;
      rsq = v1 * v1 + v2 * v2;           // see if they are in the unit circle,
    } while (rsq >= 1.0 || rsq == 0.0);  // and if they are not, try again.
    fac = sqrt(-2.0 * log(rsq) / rsq);
    // Now make the Box-Muller transformation to get two normal deviates. Return
    // one and  save the other for next time.
    gset = v1 * fac;
    iset = 1;  // Set flag.
    return v2 * fac;
  } else {
    // We have an extra deviate handy,
    iset = 0;     // so unset the flag,
    return gset;  // and return it.
  }
}

#endif  // COMMON_H