#ifndef BOLTZMANN_ALPHA0_H
#define BOLTZMANN_ALPHA0_H

#include "common.h"
int nowtime = 0;
char profile[300];

class boltzmann_alpha0
{
private:
  int set_q0 = 0; // Used for convergence.
  double D_n;  // = 1;
  double D_p;  // = 1;
  double pk11; // Second number represents different distribution.
  double pk21;
  double pk12;
  double pk22;
  double alpha0;
  double gamma0;
  double lambda0;
  double D_0;
  double rho_0;

  //rho, f1R, f1I fields.
  NRvector<Doub> rho;  
  NRvector<Doub> p1;
  NRvector<Doub> q1;
  NRvector<Doub> rho_t;
  NRvector<Doub> p1_t;
  NRvector<Doub> q1_t;
  NRvector<Doub> rho_n;
  NRvector<Doub> p1_n;
  NRvector<Doub> q1_n;

public:
  boltzmann_alpha0(/* args */);
  ~boltzmann_alpha0();
  void read_parameters();
  void initialize();
  void differ_nonlinear();
  void implicit_method();
  void convergence();
  void output(int ntime);
  void evolution();
};

boltzmann_alpha0::boltzmann_alpha0(/* args */)
{
  set_q0 = 0;
  nowtime = 0;
  rho.resize(YN);
  p1.resize(YN);
  q1.resize(YN);
  rho_t.resize(YN);
  p1_t.resize(YN);
  q1_t.resize(YN);
  rho_n.resize(YN);
  p1_n.resize(YN);
  q1_n.resize(YN);
  read_parameters();
  initialize();
}

boltzmann_alpha0::~boltzmann_alpha0()
{
}

void boltzmann_alpha0::initialize()
{
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

  for (int i = 0; i < YN; ++i)
  {
    rho[i] = rho_0 * sig0 + sin(2 * PI * i / YN + 3 * PI / 2) * rho_0 * sig0;

    //    p1[i] = (sin(2 * PI * i / YN + 3 * PI / 2) + 1.0) * 1.0;
    q1[i] = 0;

    //    p1_t[i] = p1[i];
    q1_t[i] = 0;
  }

  sprintf(profile, "%s", "profile");
}

void boltzmann_alpha0::read_parameters()
{
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
}

void boltzmann_alpha0::differ_nonlinear()
{
  NRvector<Doub> deltapq(YN);
  NRvector<Doub> deltap_q(YN);
  NRvector<Doub> muu(YN);

  double coeff =
      1.0 / (PI * lambda0 * (pk21 - 1) -
             2 * alpha0 * sig0 * rho_0 * (8 / 15.0 * pk21 + 248 / 63.0));
  double Betaa = -8.0 * alpha0 * D_n * (pk21 + 1 / 15.0) * coeff;
  double xii = 128 * alpha0 * alpha0 * (13 / 9.0 - (1 + 6 * sqrt(2)) * pk11) *
               (pk21 + 1 / 15.0) * coeff / (35 * PI);
  double gammaa = 32 * gamma0 * alpha0 * D_n *
                  (13 / 9.0 - (1 + 6 * sqrt(2)) * pk11) * coeff / (35.0);

  int i;
  for (i = 0; i < YN; i++)
  {
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
  for (i = 0; i < YN; i++)
  {
    int Ai = i + 1;
    int Bi = i - 1;

    if (Ai >= YN)
      Ai = 0;
    if (Bi < 0)
      Bi = YN - 1;

    rho_n[i] = -0.5 * D_n * (p1[Ai] + p1[Bi] - 2 * p1[i]) * dyy;
  }

  /*
   * f_1^R Equation.
   */
  for (i = 0; i < YN; i++)
  {
    int Ai = i + 1;
    int Bi = i - 1;

    if (Ai >= YN)
      Ai = 0;
    if (Bi < 0)
      Bi = YN - 1;

    p1_n[i] = -0.25 * D_n * (rho[Ai] + rho[Bi] - 2 * rho[i]) * dyy;

    p1_n[i] -= 0.25 * Betaa * (deltap_q[Ai] + deltap_q[Bi] - 2 * deltap_q[i]) *
               dyy; // Beta term.
    p1_n[i] += 0.125 * gammaa * q1[i] * (q1[Ai] + q1[Bi] - 2 * q1[i]) *
               dyy; // Gamma term.

    p1_n[i] += muu[i] * p1[i]; // \mu term.

    p1_n[i] -= xii * (p1[i] * p1[i] + q1[i] * q1[i]) * p1[i]; // \xi term.
  }

  /*
   * f_1^I Equation.
   */
  for (i = 0; i < YN; i++)
  {
    int Ai = i + 1;
    int Bi = i - 1;

    if (Ai >= YN)
      Ai = 0;
    if (Bi < 0)
      Bi = YN - 1;

    q1_n[i] =
        -0.50 * Betaa * (deltapq[Ai] + deltapq[Bi] - 2 * deltapq[i]) * dyy;
    q1_n[i] -= 0.125 * gammaa * q1[i] * (p1[Ai] + p1[Bi] - 2 * p1[i]) * dyy;

    q1_n[i] += muu[i] * q1[i];
    q1_n[i] -= xii * (p1[i] * p1[i] + q1[i] * q1[i]) * q1[i];
  }
}

void boltzmann_alpha0::convergence()
{
  double converge_rho = 0;
  double converge_p1 = 0;
  double converge_q1 = 0;
  int i;
  for (i = 0; i < YN; ++i)
  {
    converge_rho += (rho_t[i] - rho[i]) * (rho_t[i] - rho[i]);
    converge_p1 += (p1_t[i] - p1[i]) * (p1_t[i] - p1[i]);
    converge_q1 += (q1_t[i] - q1[i]) * (q1_t[i] - q1[i]);
  }
  if (set_q0 == 0)
  {
    double NN = 0;
    double PNN = 0;
    for (i = 0; i < YN; ++i)
    {
      NN += rho[i];
      PNN += fabs(p1[i]);
    }
    if (PNN > NN)
      exit(0);
  }
  if (nowtime % 2000 == 0)
    printf("cr,cp1,cq1=%d,%d,%e,%e,%e\n", set_q0, nowtime, converge_rho,
           converge_p1, converge_q1);

  if ((converge_rho < 1e-17) && (converge_p1 < 1e-17) &&
      (converge_q1 < 1e-17))
  {
    if (set_q0 == 0)
    {
      set_q0 = 1;

      Normaldev gauss(0, 1, RSD[1]);
      for (i = 0; i < YN; i++)
      {
        q1_t[i] = p1_t[i] * 0.01 * gauss.dev();
      }

      nowtime = 0;

      sprintf(profile, "%s", "decay");
    }
    else
    {
      if ((converge_rho < 1e-18) && (converge_p1 < 1e-18) &&
          (converge_q1 < 1e-20))
        exit(0);
    }
  }
}

void boltzmann_alpha0::implicit_method()
{
  int i;
  double up_beta;
  double down_alpha;

  NRvector<Doub> Di(YN);
  NRvector<Doub> aa(YN);
  NRvector<Doub> bb(YN);
  NRvector<Doub> cc(YN);
  NRvector<Doub> rr(YN);
  NRvector<Doub> ss(YN);

  double coeff =
      1.0 / (PI * lambda0 * (pk21 - 1) -
             2 * alpha0 * sig0 * rho_0 * (8 / 15.0 * pk21 + 248 / 63.0));
  double gammaa = 32 * gamma0 * alpha0 * D_n *
                  (13 / 9.0 - (1 + 6 * sqrt(2)) * pk11) * coeff / (35.0);

  /*
   * f_1^R Equation.
   */
  for (i = 0; i < YN; ++i)
  {
    Di[i] = 0.5 * D_p + 0.125 * gammaa * p1[i]; // Gamma terms introduce diffusion. D_p = 1.
  }

  for (i = 0; i < YN; ++i)
  {
    aa[i] = -Di[i] * dt * dyy;
    bb[i] = 1 + 2 * Di[i] * dt * dyy;
    cc[i] = -Di[i] * dt * dyy;
    rr[i] = p1[i] + p1_n[i] * dt;
  }
  down_alpha = -Di[YN - 1] * dt * dyy;
  up_beta = -Di[0] * dt * dyy;

  cyclic(aa, bb, cc, down_alpha, up_beta, rr, ss);

  for (i = 0; i < YN; ++i)
  {
    p1_t[i] = ss[i];
  }

  /*
   * f_1^I Equation.
   */
  for (i = 0; i < YN; ++i)
  {
    Di[i] = 0.5 * D_p + 0.125 * gammaa * p1[i]; // Gamma terms introduce diffusion. D_p = 1.
  }

  for (i = 0; i < YN; ++i)
  {
    aa[i] = -Di[i] * dt * dyy;
    bb[i] = 1 + 2 * Di[i] * dt * dyy;
    cc[i] = -Di[i] * dt * dyy;
    rr[i] = q1[i] + q1_n[i] * dt;
  }
  down_alpha = -Di[YN - 1] * dt * dyy;
  up_beta = -Di[0] * dt * dyy;

  cyclic(aa, bb, cc, down_alpha, up_beta, rr, ss);

  for (i = 0; i < YN; ++i)
  {
    q1_t[i] = ss[i];
  }

  /*
   * \rho Equation.
   */
  for (i = 0; i < YN; ++i)
  {
    Di[i] = 0.5 * D_p;
  }

  for (i = 0; i < YN; ++i)
  {
    aa[i] = -Di[i] * dt * dyy;
    bb[i] = 1 + 2 * Di[i] * dt * dyy;
    cc[i] = -Di[i] * dt * dyy;
    rr[i] = rho[i] + rho_n[i] * dt;
  }
  down_alpha = -Di[YN - 1] * dt * dyy;
  up_beta = -Di[0] * dt * dyy;

  cyclic(aa, bb, cc, down_alpha, up_beta, rr, ss);

  for (i = 0; i < YN; ++i)
  {
    rho_t[i] = ss[i];
  }
  /////////////////////////////////////
  convergence();

  // Transfer current values.
  for (i = 0; i < YN; ++i)
  {
    rho[i] = rho_t[i];
    p1[i] = p1_t[i];
    q1[i] = q1_t[i] * set_q0;
  }
}

void boltzmann_alpha0::output(int ntime)
{
  if (ntime % INTERTM != 0)
    return;
  char filename[300];
  sprintf(filename, "%s%s", Subdirectory, profile);
  PREP_old(filename, ntime, "w") for (int i = 0; i < YN; ++i)
  {
    fprintf(fp, "%d,%e,%e,%e\n", i, rho[i], p1[i], q1[i]);
  }
  CLSP
}

void boltzmann_alpha0::evolution()
{
  for (; nowtime <= TOTOTIME; ++nowtime)
  {
    differ_nonlinear();
    implicit_method();
    output(nowtime);
  }
}

#endif //BOLTZMANN_ALPHA0_H
