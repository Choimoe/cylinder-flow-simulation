// if [ -x "$(command -v qsub)" ]; then ./q run.sh; else ./run.sh; fi
#include <iostream>
using namespace std;

#define _USE_MATH_DEFINES
#include <stdlib.h>

#include <CL/sycl.hpp>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <vector>

#include "Eigen/Dense"

using Eigen::MatrixXd;

const int M = 100;  // number of mesh along x
const int N = 100;  // number of mesh along y

typedef sycl::buffer<double, 2> mat_buf;
const auto siz = sycl::range<2>(M + 1, N + 1);
const auto AC_READ = sycl::access::mode::read;
const auto AC_WRITE = sycl::access::mode::write;

double d1 = 2;   // diameter of inner circle
double d2 = 40;  // diameter of outter circle

const int DIV = 1000;
const int SEC = 100;
const int EPOCH = DIV * SEC;

// const int DIV = 100;
// const int SEC = 1;
// const int EPOCH = 5000;

const double Re = 200;

const double dt = 1.0 / DIV;
const double err = 0.005;

// const double dt = 0.01;
// const double err = 0.01;

double u0 = 1;

MatrixXd X        (M + 1, N + 1);
MatrixXd Y        (M + 1, N + 1);
MatrixXd X_eta    (M + 1, N + 1);
MatrixXd X_xi     (M + 1, N + 1);
MatrixXd Y_eta    (M + 1, N + 1);
MatrixXd Y_xi     (M + 1, N + 1);
MatrixXd X_etaeta (M + 1, N + 1);
MatrixXd X_xixi   (M + 1, N + 1);
MatrixXd X_xieta  (M + 1, N + 1);
MatrixXd Y_etaeta (M + 1, N + 1);
MatrixXd Y_xixi   (M + 1, N + 1);
MatrixXd Y_xieta  (M + 1, N + 1);
MatrixXd J        (M + 1, N + 1);
MatrixXd Eta_x    (M + 1, N + 1);
MatrixXd Xi_x     (M + 1, N + 1);
MatrixXd Eta_y    (M + 1, N + 1);
MatrixXd Xi_y     (M + 1, N + 1);
MatrixXd Eta_xx   (M + 1, N + 1);
MatrixXd Xi_xx    (M + 1, N + 1);
MatrixXd Eta_xy   (M + 1, N + 1);
MatrixXd Eta_yy   (M + 1, N + 1);
MatrixXd Xi_yy    (M + 1, N + 1);
MatrixXd Xi_xy    (M + 1, N + 1);

const double d = 1.0;

MatrixXd U  (M + 1, N + 1);
MatrixXd V  (M + 1, N + 1);
MatrixXd W  (M + 1, N + 1);
MatrixXd PSI(M + 1, N + 1);

sycl::queue calc_queue;

void readPortfolio(MatrixXd& matrix, string b) {
  ifstream data(b);
  string lineOfData;

  if (data.is_open()) {
    int i = 0;
    while (data.good()) {
      char linebuff[4096];
      getline(data, lineOfData);
      strncpy(linebuff, lineOfData.c_str(), sizeof(linebuff) - 1);
      {
        int j = 0;
        double val;
        char* p_val;
        p_val = strtok(linebuff, " ,;");
        while (NULL != p_val) {
          val = atof(p_val);
          matrix(i, j) = val;
          j++;
          p_val = strtok(NULL, " ,;");
        }
      }
      i++;
    }
  } else {
    std::cout << "Unable to open file";
  }
}

void Jacobi() {
  int index = 0;
  X_xi      = MatrixXd::Zero(M + 1, N + 1);
  Y_xi      = MatrixXd::Zero(M + 1, N + 1);
  X_eta     = MatrixXd::Zero(M + 1, N + 1);
  Y_eta     = MatrixXd::Zero(M + 1, N + 1);
  X_xixi    = MatrixXd::Zero(M + 1, N + 1);
  X_etaeta  = MatrixXd::Zero(M + 1, N + 1);
  Y_xixi    = MatrixXd::Zero(M + 1, N + 1);
  Y_etaeta  = MatrixXd::Zero(M + 1, N + 1);
  J = MatrixXd::Zero(M + 1, N + 1);
  for (int i = 0; i < M; i++) {
    if (i > 0) {
      index = i - 1;
    } else {
      index = M - 1;
    }
    for (int j = 1; j < N; j++) {
      X_xi(i, j)      = (X(i + 1, j) - X(index, j)) / 2.0 / d;
      Y_xi(i, j)      = (Y(i + 1, j) - Y(index, j)) / 2.0 / d;
      X_eta(i, j)     = (X(i, j + 1) - X(i, j - 1)) / 2.0 / d;
      Y_eta(i, j)     = (Y(i, j + 1) - Y(i, j - 1)) / 2.0 / d;
      X_xixi(i, j)    = (X(i + 1, j) - 2 * X(i, j) + X(index, j)) / d / d;
      Y_xixi(i, j)    = (Y(i + 1, j) - 2 * Y(i, j) + Y(index, j)) / d / d;
      X_etaeta(i, j)  = (X(i, j + 1) - 2 * X(i, j) + X(i, j - 1)) / d / d;
      Y_etaeta(i, j)  = (Y(i, j + 1) - 2 * Y(i, j) + Y(i, j - 1)) / d / d;
      X_xieta(i, j)   = (X(i + 1, j + 1) - X(i + 1, j - 1) - X(index, j + 1) +
                         X(index, j - 1)) /
                        4.0 / d / d;
      Y_xieta(i, j)   = (Y(i + 1, j + 1) - Y(i + 1, j - 1) - Y(index, j + 1) +
                         Y(index, j - 1)) /
                        4.0 / d / d;
    }
  }

  for (int i = 0; i < M; i++) {
    index = i ? i - 1 : M - 1;
    int j = N;
    X_xi(i, j)      = (X(i + 1, j) - X(index, j)) / 2.0 / d;
    Y_xi(i, j)      = (Y(i + 1, j) - Y(index, j)) / 2.0 / d;
    X_eta(i, j)     = (X(i, j) - X(i, j - 1)) / d;
    Y_eta(i, j)     = (Y(i, j) - Y(i, j - 1)) / d;
    X_xixi(i, j)    = (X(i + 1, j) - 2 * X(i, j) + X(index, j)) / d / d;
    Y_xixi(i, j)    = (Y(i + 1, j) - 2 * Y(i, j) + Y(index, j)) / d / d;
    X_etaeta(i, j)  = X_etaeta(i, j - 1);
    Y_etaeta(i, j)  = Y_etaeta(i, j - 1);
    X_xieta(i, j)   =
          (X(i + 1, j) - X(i + 1, j - 1) - X(index, j) + X(index, j - 1)) / 2.0 /
          d / d;
    Y_xieta(i, j)   =
          (Y(i + 1, j) - Y(i + 1, j - 1) - Y(index, j) + Y(index, j - 1)) / 2.0 /
          d / d;
  }
  for (int i = 0; i < M; i++) {
    index = i ? i - 1 : M - 1;
    int j = 0;
    X_xi(i, j)      = (X(i + 1, j) - X(index, j)) / 2.0 / d;
    Y_xi(i, j)      = (Y(i + 1, j) - Y(index, j)) / 2.0 / d;
    X_eta(i, j)     = (X(i, j + 1) - X(i, j)) / d;
    Y_eta(i, j)     = (Y(i, j + 1) - Y(i, j)) / d;
    X_xixi(i, j)    = (X(i + 1, j) - 2 * X(i, j) + X(index, j)) / d / d;
    Y_xixi(i, j)    = (Y(i + 1, j) - 2 * Y(i, j) + Y(index, j)) / d / d;
    X_etaeta(i, j)  = X_etaeta(i, j + 1);
    Y_etaeta(i, j)  = Y_etaeta(i, j + 1);
    X_xieta(i, j)   =
          (X(i + 1, j + 1) - X(i + 1, j) - X(index, j + 1) + X(index, j)) / 2.0 /
          d / d;
    Y_xieta(i, j)   =
          (Y(i + 1, j + 1) - Y(i + 1, j) - Y(index, j + 1) + Y(index, j)) / 2.0 /
          d / d;
  }

  for (int i = 0; i < M; i++)
    for (int j = 0; j <= N; j++)
      J(i, j) = X_xi(i, j) * Y_eta(i, j) - X_eta(i, j) * Y_xi(i, j);
}

void deriviate() {
  Eigen::Vector3d A, B, C, D;
  Eigen::Matrix3d E;
  double temp;
  Xi_x    = MatrixXd::Zero(M + 1, N + 1);
  Eta_x   = MatrixXd::Zero(M + 1, N + 1);
  Xi_y    = MatrixXd::Zero(M + 1, N + 1);
  Eta_y   = MatrixXd::Zero(M + 1, N + 1);
  Xi_xx   = MatrixXd::Zero(M + 1, N + 1);
  Xi_xy   = MatrixXd::Zero(M + 1, N + 1);
  Xi_yy   = MatrixXd::Zero(M + 1, N + 1);
  Eta_xx  = MatrixXd::Zero(M + 1, N + 1);
  Eta_xy  = MatrixXd::Zero(M + 1, N + 1);
  Eta_yy  = MatrixXd::Zero(M + 1, N + 1);
  for (int i = 0; i < M; i++) {
    for (int j = 0; j <= N; j++) {
      Xi_x(i, j)  =  Y_eta(i, j)  / J(i, j);
      Eta_x(i, j) = -Y_xi(i, j)   / J(i, j);
      Xi_y(i, j)  = -X_eta(i, j)  / J(i, j);
      Eta_y(i, j) =  X_xi(i, j)   / J(i, j);
      A << Xi_x(i, j) * X_xixi(i, j)   + Xi_y(i, j) * Y_xixi(i, j),
           Xi_x(i, j) * X_xieta(i, j)  + Xi_y(i, j) * Y_xieta(i, j),
           Xi_x(i, j) * X_etaeta(i, j) + Xi_y(i, j) * Y_etaeta(i, j);
      E << pow(Y_eta(i, j), 2), -2 * Y_eta(i, j) * Y_xi(i, j),
          pow(Y_xi(i, j), 2), -X_eta(i, j) * Y_eta(i, j),
          (X_xi(i, j) * Y_eta(i, j) + X_eta(i, j) * Y_xi(i, j)),
          -X_xi(i, j) * Y_xi(i, j), pow(X_eta(i, j), 2),
          -2 * X_xi(i, j) * X_eta(i, j), pow(X_xi(i, j), 2);
      temp = pow(X_xi(i, j) * Y_eta(i, j) - X_eta(i, j) * Y_xi(i, j), 2);
      B = -E * A / temp;
      C << Eta_x(i, j) * X_xixi(i, j)   + Eta_y(i, j) * Y_xixi(i, j),
           Eta_x(i, j) * X_xieta(i, j)  + Eta_y(i, j) * Y_xieta(i, j),
           Eta_x(i, j) * X_etaeta(i, j) + Eta_y(i, j) * Y_etaeta(i, j);
      D = -E * C / temp;
      Xi_xx(i, j)   = B(0);
      Xi_xy(i, j)   = B(1);
      Xi_yy(i, j)   = B(2);
      Eta_xx(i, j)  = D(0);
      Eta_xy(i, j)  = D(1);
      Eta_yy(i, j)  = D(2);
    }
  }
}

MatrixXd ZERO = MatrixXd::Zero(M + 1, N + 1);

MatrixXd pre_A(M + 1, N + 1);
MatrixXd pre_B(M + 1, N + 1);
MatrixXd pre_R(M + 1, N + 1);
MatrixXd pre_O(M + 1, N + 1);
MatrixXd pre_E(M + 1, N + 1);

mat_buf buf_pre_A (ZERO.data(), siz);
mat_buf buf_pre_B (ZERO.data(), siz);
mat_buf buf_pre_R (ZERO.data(), siz);
mat_buf buf_pre_O (ZERO.data(), siz);
mat_buf buf_pre_E (ZERO.data(), siz);

mat_buf buf_Xi_x  (ZERO.data(), siz);
mat_buf buf_Xi_y  (ZERO.data(), siz);
mat_buf buf_Eta_x (ZERO.data(), siz);
mat_buf buf_Eta_y (ZERO.data(), siz);

mat_buf buf_X_xi  (ZERO.data(), siz);
mat_buf buf_X_eta (ZERO.data(), siz);
mat_buf buf_Y_xi  (ZERO.data(), siz);
mat_buf buf_Y_eta (ZERO.data(), siz);

void pre_process() {
  pre_A = MatrixXd::Zero(M + 1, N + 1);
  pre_B = MatrixXd::Zero(M + 1, N + 1);
  pre_R = MatrixXd::Zero(M + 1, N + 1);
  pre_O = MatrixXd::Zero(M + 1, N + 1);
  pre_E = MatrixXd::Zero(M + 1, N + 1);

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      pre_A(i, j) = pow(Xi_x(i, j), 2) + pow(Xi_y(i, j), 2);
      pre_B(i, j) = Xi_x(i, j) * Eta_x(i, j) + Xi_y(i, j) * Eta_y(i, j);
      pre_R(i, j) = pow(Eta_x(i, j), 2) + pow(Eta_y(i, j), 2);
      pre_O(i, j) = Xi_xx(i, j) + Xi_yy(i, j);
      pre_E(i, j) = Eta_xx(i, j) + Eta_yy(i, j);
    }
  }

  buf_pre_A = mat_buf(pre_A.data(), siz);
  buf_pre_B = mat_buf(pre_B.data(), siz);
  buf_pre_R = mat_buf(pre_R.data(), siz);
  buf_pre_O = mat_buf(pre_O.data(), siz);
  buf_pre_E = mat_buf(pre_E.data(), siz);

  buf_Xi_x  = mat_buf(Xi_x.data() , siz);
  buf_Xi_y  = mat_buf(Xi_y.data() , siz);
  buf_Eta_x = mat_buf(Eta_x.data(), siz);
  buf_Eta_y = mat_buf(Eta_y.data(), siz);

  buf_X_xi  = mat_buf(X_xi.data() , siz);
  buf_Y_xi  = mat_buf(Y_xi.data() , siz);
  buf_X_eta = mat_buf(X_eta.data(), siz);
  buf_Y_eta = mat_buf(Y_eta.data(), siz);
}

void push() {
  MatrixXd W_new(M + 1, N + 1);
  W_new = MatrixXd::Zero(M + 1, N + 1);

  mat_buf buf_U(U.data(), siz);
  mat_buf buf_V(V.data(), siz);
  mat_buf buf_W(W.data(), siz);
  mat_buf buf_new_W(W_new.data(), siz);

  auto tas = calc_queue.submit([&](sycl::handler& h) {
    auto pA     = buf_pre_A.get_access<AC_READ>(h);
    auto pB     = buf_pre_B.get_access<AC_READ>(h);
    auto pR     = buf_pre_R.get_access<AC_READ>(h);
    auto pO     = buf_pre_O.get_access<AC_READ>(h);
    auto pE     = buf_pre_E.get_access<AC_READ>(h);
    auto pXi_x  = buf_Xi_x.get_access <AC_READ>(h);
    auto pXi_y  = buf_Xi_y.get_access <AC_READ>(h);
    auto pEta_x = buf_Eta_x.get_access<AC_READ>(h);
    auto pEta_y = buf_Eta_y.get_access<AC_READ>(h);
    auto pU     = buf_U.get_access    <AC_READ>(h);
    auto pV     = buf_V.get_access    <AC_READ>(h);
    auto pW     = buf_W.get_access    <AC_READ>(h);

    auto pNew_W = buf_new_W.get_access<AC_WRITE>(h);

    h.parallel_for(sycl::range<2>(M, N), [=](sycl::id<2> idx) {
      int i = idx[0], index = i ? i - 1 : M - 1;
      int j = idx[1];
      if (j != 0) {
        double a = pA[j][i], b = pB[j][i], r = pR[j][i], o = pO[j][i],
               e = pE[j][i];
        double b_e = (a + o / 2.0) / Re, b_w = (a - o / 2.0) / Re,
               b_s = (r + e / 2.0) / Re, b_n = (r - e / 2.0) / Re;
        double b_w1 = (pU[j][i] * pXi_x[j][i] + pV[j][i] * pXi_y[j][i]) / 2.0,
               b_e1 = -b_w1,
               b_n1 = (pU[j][i] * pEta_x[j][i] + pV[j][i] * pEta_y[j][i]) / 2.0,
               b_s1 = -b_n1;
        pNew_W[j][i] =
            pW[j][i] +
            dt * (pW[j][i + 1] * (b_e + b_e1) + pW[j][index] * (b_w + b_w1) +
                  pW[j + 1][i] * (b_s + b_s1) + pW[j - 1][i] * (b_n + b_n1) -
                  2.0 / Re * (a + r) * pW[j][i] +
                  2.0 * b / Re *
                      (pW[j + 1][i + 1] - pW[j + 1][i - 1] - pW[j + 1][index] +
                       pW[j - 1][index]));
      }
    });
  });

  for (int j = 0; j < N; j++) W_new(M, j) = W_new(0, j);

  tas.wait();

  W = W_new;
}

void psi_iteration() {
  int cnt = 0;
  MatrixXd psi_new(M + 1, N + 1);
  psi_new = PSI;
  mat_buf buf_W(W.data(), siz);

  while (1) {
    mat_buf buf_PSI(PSI.data(), siz);
    mat_buf buf_new_PSI(psi_new.data(), siz);

    auto tas = calc_queue.submit([&](sycl::handler& h) {
      auto pA       = buf_pre_A.get_access  <AC_READ>(h);
      auto pB       = buf_pre_B.get_access  <AC_READ>(h);
      auto pR       = buf_pre_R.get_access  <AC_READ>(h);
      auto pO       = buf_pre_O.get_access  <AC_READ>(h);
      auto pE       = buf_pre_E.get_access  <AC_READ>(h);

      auto pW       = buf_W.get_access      <AC_READ>(h);
      auto pPSI     = buf_PSI.get_access    <AC_READ>(h);
      auto pNew_PSI = buf_new_PSI.get_access<AC_WRITE>(h);

      h.parallel_for(sycl::range<2>(M, N), [=](sycl::id<2> idx) {
        int i = idx[0], index = i ? i - 1 : M - 1;
        int j = idx[1];
        if (j != 0) {
          double a = pA[j][i], b = pB[j][i], r = pR[j][i], o = pO[j][i],
                 e = pE[j][i];

          pNew_PSI[j][i] =
              ((a - o / 2.0) * pPSI[j][index] + (a + o / 2.0) * pPSI[j][i + 1] +
               (r - e / 2.0) * pPSI[j - 1][i] + (r + e / 2.0) * pPSI[j + 1][i] +
               2 * b *
                   (pPSI[j + 1][i + 1] - pPSI[j - 1][i + 1] -
                    pPSI[j + 1][index] + pPSI[j - 1][index]) -
               pW[j][i]) /
              (2 * a + 2 * r);
        }
      });
    });

    tas.wait();

    for (int j = 1; j < N; j++) psi_new(M, j) = psi_new(0, j);

    double error =
        (PSI.block(0, 1, M + 1, N - 1) - psi_new.block(0, 1, M + 1, N - 1))
            .norm();

    PSI = psi_new;

    if (error < err) {
      return;
    }
  }
}

void velocity() {
  int index = 0;
  mat_buf buf_J   (J.data(), siz);
  mat_buf buf_PSI (PSI.data(), siz);

  mat_buf buf_U   (U.data(), siz);
  mat_buf buf_V   (V.data(), siz);

  auto tas = calc_queue.submit([&](sycl::handler& h) {
    auto pPSI   = buf_PSI.get_access  <AC_READ>(h);
    auto pJ     = buf_J.get_access    <AC_READ>(h);
    auto pX_xi  = buf_X_xi.get_access <AC_READ>(h);
    auto pY_xi  = buf_Y_xi.get_access <AC_READ>(h);
    auto pX_eta = buf_X_eta.get_access<AC_READ>(h);
    auto pY_eta = buf_Y_eta.get_access<AC_READ>(h);
    
    auto pU     = buf_U.get_access    <AC_WRITE>(h);
    auto pV     = buf_V.get_access    <AC_WRITE>(h);

    h.parallel_for(sycl::range<2>(M, N), [=](sycl::id<2> idx) {
      int i = idx[0], index = i ? i - 1 : M - 1;
      int j = idx[1];
      if (j != 0) {
        pU[j][i] = pX_xi[j][i] / pJ[j][i]  * (pPSI[j + 1][i] - pPSI[j - 1][i]) /
                       2.0 / d -
                   pX_eta[j][i] / pJ[j][i] * (pPSI[j][i + 1] - pPSI[j][index]) /
                       2.0 / d;
        pV[j][i] = pY_xi[j][i] / pJ[j][i]  * (pPSI[j + 1][i] - pPSI[j - 1][i]) /
                       2.0 / d -
                   pY_eta[j][i] / pJ[j][i] * (pPSI[j][i + 1] - pPSI[j][index]) /
                       2.0 / d;
      }
    });
  });

  tas.wait();
  U.row(M) = U.row(0);
  V.row(M) = V.row(0);
}

void clear() {
  double sum = 0;
  for (int i = 0; i < M; i++) sum += PSI(i, 1);

  double ave = sum / M;
  for (int i = 0; i < M; i++) PSI(i, 0) = ave;

  for (int i = 0; i < M; i++) {
    W(i, 0) = 2 * PSI(i, 1) * (pow(Eta_x(i, 0), 2) + pow(Eta_y(i, 0), 2));
  }
}

void init() {
  U = u0 * MatrixXd::Ones(M + 1, N + 1);
  U.row(M)  = MatrixXd::Zero(1, N + 1);
  U.col(0)  = MatrixXd::Zero(M + 1, 1);

  V         = MatrixXd::Zero(M + 1, N + 1);
  PSI       = MatrixXd::Zero(M + 1, N + 1);
  W         = MatrixXd::Zero(M + 1, N + 1);

  for (int i = 0; i < M; i++)
    for (int j = 1; j <= N; j++) PSI(i, j) = 0;

  clear();
}

void boundary() {
  int j = N;
  for (int i = 0; i < M; i++) {
    V(i, j)   = 0;
    W(i, j)   = 0;
    U(i, j)   = u0;
    PSI(i, j) = u0 * Y(i, j);
  }

  clear();
}

void save(int k) {
  if (k % DIV) return;

  char filename[30];
  sprintf(filename, "./output/PSI%d.dat", k / DIV);
  ofstream psiout(filename);
  psiout << PSI;
  psiout.close();
  sprintf(filename, "./output/U%d.dat", k / DIV);
  ofstream Uout(filename);
  Uout << U;
  Uout.close();
  sprintf(filename, "./output/V%d.dat", k / DIV);
  ofstream Vout(filename);
  Vout << V;
  Vout.close();
}

int main() {
  clock_t start, end;

  readPortfolio(X, "X.txt");
  readPortfolio(Y, "Y.txt");

  start = clock();

  Jacobi();
  deriviate();
  init();
  pre_process();

  PSI.row(M) = PSI.row(0);

  for (int k = 1; k <= EPOCH; k++) {
    // printf("[DEBUG] k = %d\n", k);
    push();
    psi_iteration();
    velocity();
    boundary();

    save(k);
  }

  end = clock();
  std::cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s"
            << std::endl;
  return 0;
}