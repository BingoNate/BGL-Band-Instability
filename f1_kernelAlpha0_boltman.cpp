#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "common.h"

const int YN = 256;
const double PI = 3.141592653589793;
const double D_r = 1;

char profile[300];
int TOTOTIME = 50000000;
int INTERTM = 20000;
int set_q0 = 0;
int nowtime = 0;
long RSD[1] = {-1};

double D_n;   // = 1;
double D_p;   // = 1;
double pk11;  // Second number represents different distribution.
double pk21;
double pk12;
double pk22;
double alpha0;
double gamma0;
double lambda0;
double D_0;
double rho_0;

double sig0 = 1.01;
double DY = 1.5;
double dt = 0.02;
double dyy = 1.0 / (DY * DY);

double rho[YN], p1[YN], q1[YN];
double rho_t[YN], p1_t[YN], q1_t[YN];

// According to Boltzmann Equations, the coefficients are changed.
void diff_nonlinear(double *p1_n, double *q1_n, double *rho_n) {
  double deltapq[YN];
  double deltap_q[YN];
  double muu[YN];

  int i;

  double coeff =
      1.0 / (PI * lambda0 * (pk21 - 1) -
             2 * alpha0 * sig0 * rho_0 * (8 / 15.0 * pk21 + 248 / 63.0));
  double Betaa = -8.0 * alpha0 * D_n * (pk21 + 1 / 15.0) * coeff;
  double xii = 128 * alpha0 * alpha0 * (13 / 9.0 - (1 + 6 * sqrt(2)) * pk11) *
               (pk21 + 1 / 15.0) * coeff / (35 * PI);
  double gammaa = 32 * gamma0 * alpha0 * D_n *
                  (13 / 9.0 - (1 + 6 * sqrt(2)) * pk11) * coeff / (35.0);
  double dd0 = 0.5;

  for (i = 0; i < YN; i++) {
    muu[i] = 16 * alpha0 / (3.0 * PI) * ((2 * sqrt(2) - 1) * pk11 - 7 / 5.0) *
                 rho[i] -
             lambda0 * (1 - pk11);
    deltap_q[i] = (p1[i] * p1[i] - q1[i] * q1[i]);
    deltapq[i] = (p1[i] * q1[i]);

    p1_n[i] = 0;
    q1_n[i] = 0;
    rho_n[i] = 0;
  }

  /*
   * Density Equation.
   */
  for (i = 0; i < YN; i++) {
    int Ai = i + 1;
    int Bi = i - 1;

    if (Ai >= YN) Ai = 0;
    if (Bi < 0) Bi = YN - 1;

    rho_n[i] = -0.5 * D_n * (p1[Ai] + p1[Bi] - 2 * p1[i]) * dyy;
  }

  /*
   * f_1^R Equation.
   */
  for (i = 0; i < YN; i++) {
    int Ai = i + 1;
    int Bi = i - 1;

    if (Ai >= YN) Ai = 0;
    if (Bi < 0) Bi = YN - 1;

    p1_n[i] = -0.25 * D_n * (rho[Ai] + rho[Bi] - 2 * rho[i]) * dyy;

    p1_n[i] -= 0.25 * Betaa * (deltap_q[Ai] + deltap_q[Bi] - 2 * deltap_q[i]) *
               dyy;  // Beta term.
    p1_n[i] += 0.125 * gammaa * q1[i] * (q1[Ai] + q1[Bi] - 2 * q1[i]) *
               dyy;  // Gamma term.

    p1_n[i] += muu[i] * p1[i];  // \mu term.

    p1_n[i] -= xii * (p1[i] * p1[i] + q1[i] * q1[i]) * p1[i];  // \xi term.
  }

  /*
   * f_1^I Equation.
   */
  for (i = 0; i < YN; i++) {
    int Ai = i + 1;
    int Bi = i - 1;

    if (Ai >= YN) Ai = 0;
    if (Bi < 0) Bi = YN - 1;

    q1_n[i] =
        -0.50 * Betaa * (deltapq[Ai] + deltapq[Bi] - 2 * deltapq[i]) * dyy;
    q1_n[i] -= 0.125 * gammaa * q1[i] * (p1[Ai] + p1[Bi] - 2 * p1[i]) * dyy;

    q1_n[i] += muu[i] * q1[i];
    q1_n[i] -= xii * (p1[i] * p1[i] + q1[i] * q1[i]) * q1[i];
  }
}

void implicit_method(double *p1_n, double *q1_n, double *rho_n) {
  double *aa, *bb, *cc, *rr, *ss;

  double *Di;
  Di = double_vector(1, YN);
  aa = double_vector(1, YN);

  bb = double_vector(1, YN);
  cc = double_vector(1, YN);

  rr = double_vector(1, YN);
  ss = double_vector(1, YN);

  double up_beta;
  double down_alpha;

  int i;
  double coeff =
      1.0 / (PI * lambda0 * (pk21 - 1) -
             2 * alpha0 * sig0 * rho_0 * (8 / 15.0 * pk21 + 248 / 63.0));
  double gammaa = 32 * gamma0 * alpha0 * D_n *
                  (13 / 9.0 - (1 + 6 * sqrt(2)) * pk11) * coeff / (35.0);
  /*
   * f_1^R Equation.
   */
  for (i = 0; i < YN; i++) {
    Di[i + 1] =
        0.5 * D_p +
        0.125 * gammaa * p1[i];  // Gamma terms introduce diffusion. D_p = 1.
  }

  for (i = 1; i <= YN; i++) {
    aa[i] = -Di[i] * dt * dyy;
    bb[i] = 1 + 2 * Di[i] * dt * dyy;
    cc[i] = -Di[i] * dt * dyy;
    rr[i] = p1[i - 1] + p1_n[i - 1] * dt;
  }

  down_alpha = -Di[YN] * dt * dyy;
  up_beta = -Di[1] * dt * dyy;

  cyclic(aa, bb, cc, down_alpha, up_beta, rr, ss, YN);

  for (i = 1; i <= YN; i++) {
    p1_t[i - 1] = ss[i];
  }

  /*
   * f_1^I Equation.
   */
  for (i = 0; i < YN; i++) {
    Di[i + 1] = 0.5 * D_p + 0.125 * gammaa * p1[i];
  }

  for (i = 1; i <= YN; i++) {
    aa[i] = -Di[i] * dt * dyy;
    bb[i] = 1 + 2 * Di[i] * dt * dyy;
    cc[i] = -Di[i] * dt * dyy;
    rr[i] = q1[i - 1] + q1_n[i - 1] * dt;
  }

  down_alpha = -Di[YN] * dt * dyy;
  up_beta = -Di[1] * dt * dyy;

  cyclic(aa, bb, cc, down_alpha, up_beta, rr, ss, YN);

  for (i = 1; i <= YN; i++) {
    q1_t[i - 1] = ss[i];
  }

  /*
   * \rho Equation.
   */
  for (i = 0; i < YN; i++) {
    Di[i + 1] = 0.5 * D_p;
  }

  for (i = 1; i <= YN; i++) {
    aa[i] = -Di[i] * dt * dyy;
    bb[i] = 1 + 2 * Di[i] * dt * dyy;
    cc[i] = -Di[i] * dt * dyy;
    rr[i] = rho[i - 1] + rho_n[i - 1] * dt;
  }

  down_alpha = -Di[YN] * dt * dyy;
  up_beta = -Di[1] * dt * dyy;

  cyclic(aa, bb, cc, down_alpha, up_beta, rr, ss, YN);

  for (i = 1; i <= YN; i++) {
    rho_t[i - 1] = ss[i];
  }

  ///////////////////////////////////////////////////////////////////////////////
  /*
   * Convergence.
   */
  double converge_rho = 0;
  double converge_p1 = 0;
  double converge_q1 = 0;

  for (i = 0; i < YN; i++) {
    converge_rho += (rho_t[i] - rho[i]) * (rho_t[i] - rho[i]);
    converge_p1 += (p1_t[i] - p1[i]) * (p1_t[i] - p1[i]);
    converge_q1 += (q1_t[i] - q1[i]) * (q1_t[i] - q1[i]);
  }

  if (set_q0 == 0) {
    double NN = 0, PNN = 0;

    for (i = 0; i < YN; i++) {
      NN += rho[i];

      PNN += fabs(p1[i]);
    }

    if (PNN > 1 * NN) exit(0);
  }

  if (nowtime % 2000 == 0)
    printf("cr,cp1,cq1=%d,%d,%e,%e,%e\n", set_q0, nowtime, converge_rho,
           converge_p1, converge_q1);

  if ((converge_rho < 1e-17) && (converge_p1 < 1e-17) &&
      (converge_q1 < 1e-17)) {
    if (set_q0 == 0) {
      set_q0 = 1;

      for (i = 0; i < YN; i++) q1_t[i] = p1_t[i] * 0.01 * gasdev(RSD);

      nowtime = 0;

      sprintf(profile, "%s", "decay");

    } else {
      if ((converge_rho < 1e-18) && (converge_p1 < 1e-18) &&
          (converge_q1 < 1e-20))

        exit(0);
    }
  }

  for (i = 0; i < YN; i++) {
    rho[i] = rho_t[i];
    p1[i] = p1_t[i];
    q1[i] = q1_t[i] * set_q0;
  }

  free_double_vector(aa, 1, YN);

  free_double_vector(bb, 1, YN);
  free_double_vector(cc, 1, YN);

  free_double_vector(rr, 1, YN);
  free_double_vector(ss, 1, YN);
  free_double_vector(Di, 1, YN);
}

void output(int nowtime) {
  if (nowtime % INTERTM != 0) return;

  int i;
  int j;

  char file_name[300];
  sprintf(file_name, "%s%s", Subdirectory, profile);

  PREP_old(file_name, nowtime, "w")

      for (i = 0; i < YN; i++) {
    fprintf(fp, "%d,%e,%e,%e\n", i, rho[i], p1[i], q1[i]);
  }

  CLSP
}

void evolution() {
  double p1_n[YN];
  double q1_n[YN];
  double rho_n[YN];

  for (; nowtime <= TOTOTIME; nowtime++) {
    diff_nonlinear(p1_n, q1_n, rho_n);

    implicit_method(p1_n, q1_n, rho_n);

    output(nowtime);
  }
}

void initialize() {
  int i;

  set_q0 = 0;

  for (i = 0; i < YN; i++) {
    rho[i] = rho_0 * sig0 + sin(2 * PI * i / YN + 3 * PI / 2) * rho_0 * sig0;

    //    p1[i] = (sin(2 * PI * i / YN + 3 * PI / 2) + 1.0) * 1.0;
    q1[i] = 0;

    //    p1_t[i] = p1[i];
    q1_t[i] = 0;
  }

  sprintf(profile, "%s", "profile");
}

int checkfiles() {
  int yesfile = 0;
  int i;
  char subfile[300];
  sprintf(subfile, "%s%s", Subdirectory, "profile");

  for (i = 0; i <= TOTOTIME; i += 1000) {
    PREP_check(subfile, i, "r");  // this obtains the file number and nowtime;
  }
  printf("checkfile::  nowtime=%d\n", nowtime);
  return yesfile;
}

void read_nowtime(int rtm) {
  double xy[2], uxy[2], vxy[2], uu;
  int a, rnn, flg;
  int nein;
  char subfile[300];

  sprintf(subfile, "%s%s", Subdirectory, "profile");

  PREP_old(subfile, rtm, "r");

  rnn = 0;
  do {
    flg = fscanf(fp, "%d,%lf,%lf,%lf", &a, &xy[0], &uxy[0], &uxy[1]);
    if (flg == 4) {
      rho[rnn] = xy[0];
      p1[rnn] = uxy[0];
      q1[rnn] = uxy[1];

      rho_t[rnn] = rho[rnn];
      p1_t[rnn] = p1[rnn];
      q1_t[rnn] = q1[rnn];

      rnn++;
    }

  } while (flg > 0);
  CLSP

      int n = rnn;  // update particle number

  printf("reading files time = %d,n = %d\n", rtm, n);
}

void read_from_file() {
  int findfile;

  findfile = checkfiles();

  if (findfile != 0) read_nowtime(nowtime);
}

void file_manu() {
  //	if(sid==0)
  {
    char subfile[300];
    sprintf(subfile, "mkdir %s", Subdirectory);  //  mkdir
    printf("%s\n", subfile);
    system(subfile);
  }
  //	MPI_Barrier(MPI_COMM_WORLD);
}

// Does not use pk21, pk22 since they are same as pk11 and pk12.
void read_parameters() {
  PREP_old("aaparameter", 0, "r")

      fscanf(fp, "%lf=pk11\n", &pk11);
  fscanf(fp, "%lf=pk21\n", &pk21);
  fscanf(fp, "%lf=pk12\n", &pk12);
  fscanf(fp, "%lf=pk22\n", &pk22);
  fscanf(fp, "%lf=alpha0\n", &alpha0);
  fscanf(fp, "%lf=gamma0\n", &gamma0);
  fscanf(fp, "%lf=lambda0\n", &lambda0);
  fscanf(fp, "%lf=D0\n", &D_0);
  fscanf(fp, "%lf=rho0\n", &rho_0);

  CLSP

      // Fixed transition density rho_t = 1.0. Then, pk11 is dependent on
      // collision rate \alpha_0.
      pk11 = (3 * PI * lambda0 + 7 * 16 * alpha0 * rho_0 / 5.0) /
             (16 * alpha0 * rho_0 * (2 * sqrt(2) - 1) + 3 * PI * lambda0);
	
//  alpha0 = 3 * PI * lambda0 * (1 - pk11) / (16 * rho_0 * ((2 * sqrt(2) - 1) * pk11 - 7 / 5.0));
  // Actually, the controlling parameters are pk12 and alpha_0;
  pk12 = pk11;
  pk22 = pk21;

  D_n = D_0;
  D_p = D_0;

  printf(
      "alpha0, lambda0, pk11, pk21, pk12, pk22, D_0, rho_0 = "
      "%f,%f,%f,%f,%f,%f,%f,%f\n",
      alpha0, lambda0, pk11, pk21, pk12, pk22, D_0, rho_0);
}

int main() {
  file_manu();
  read_parameters();

  initialize();

  nowtime = 0;

  read_from_file();

  evolution();
}
