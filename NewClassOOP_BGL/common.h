#ifndef COMMON_H
#define COMMON_H

#include "nr3.h"
#include "ran.h"
#include "tridag.h"

#ifdef _WIN32
#define Subdirectory "TKS\\"
#else // __linux
#define Subdirectory "TKS//"
#endif

#define PREP_old(Fpre, Fnum, R_W)             \
  {                                           \
    FILE *fp;                                 \
    {                                         \
      char PRT[100];                          \
      sprintf(PRT, "%s%dsA.dat", Fpre, Fnum); \
      if ((fp = fopen(PRT, R_W)) == NULL)     \
      {                                       \
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
      if ((fp = fopen(PRT, R_W)) != NULL)     \
      {                                       \
        nowtime = Fnum;                       \
        yesfile++;                            \
        {                                     \
          fclose(fp);                         \
        }                                     \
      }                                       \
    }                                         \
  }

#define CLSP    \
  {             \
    fclose(fp); \
  }             \
  }

const int TOTOTIME = 50000000;
const int INTERTM = 20000;
const int YN = 256;
const double PI = 3.141592653589793;
const double D_r = 1;
const long RSD[1] = {-1};

const double sig0 = 1.01;
const double DY = 1.5;
const double dt = 0.02;
const double dyy = 1.0 / (DY * DY);
extern char profile[300];
extern int nowtime;

#endif //COMMON_H
