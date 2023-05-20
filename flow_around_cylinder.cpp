#include <iostream>
using namespace std;
#include <Eigen/Dense>
#include <fstream>
#include <vector>
using namespace Eigen;

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <string.h>
#define random(x) (rand() % x)
double d1 = 2;   // diameter of inner circle
double d2 = 40;  // diameter of outter circle
#define M 100    // number of mesh along x
#define N 100    // number of mesh along y
double Re = 200;
double dt = 0.01;
double err;
double u0 = 1;

MatrixXd X(M + 1, N + 1);         // x coordinate
MatrixXd Y(M + 1, N + 1);         // y coordinate
MatrixXd X_eta(M + 1, N + 1);     // first order
MatrixXd X_xi(M + 1, N + 1);      // first order
MatrixXd Y_eta(M + 1, N + 1);     // first order
MatrixXd Y_xi(M + 1, N + 1);      // first order
MatrixXd X_etaeta(M + 1, N + 1);  // second order,eta xi same
MatrixXd X_xixi(M + 1, N + 1);    // second order,eta xi same
MatrixXd X_xieta(M + 1, N + 1);   // second order,eta xi same
MatrixXd Y_etaeta(M + 1, N + 1);  // second order,eta xi same
MatrixXd Y_xixi(M + 1, N + 1);    // second order,eta xi same
MatrixXd Y_xieta(M + 1, N + 1);   // second order,eta xi same
MatrixXd J(M + 1, N + 1);
MatrixXd Eta_x(M + 1, N + 1);   // first order
MatrixXd Xi_x(M + 1, N + 1);    // first order
MatrixXd Eta_y(M + 1, N + 1);   // first order
MatrixXd Xi_y(M + 1, N + 1);    // first order
MatrixXd Eta_xx(M + 1, N + 1);  // second order,eta xi same
MatrixXd Xi_xx(M + 1, N + 1);   // second order,eta xi same
MatrixXd Eta_xy(M + 1, N + 1);  // second order,eta xi same
MatrixXd Eta_yy(M + 1, N + 1);  // second order,eta xi same
MatrixXd Xi_yy(M + 1, N + 1);   // second order,eta xi same
MatrixXd Xi_xy(M + 1, N + 1);   // second order,eta xi same

double d = 1.0;

MatrixXd U(M + 1, N + 1);    // u velocity
MatrixXd V(M + 1, N + 1);    // v velocity
MatrixXd W(M + 1, N + 1);    // vorticity
MatrixXd PSI(M + 1, N + 1);  // flow function

void readPortfolio(MatrixXd& matrix, string b) {
  ifstream data(b);
  string lineOfData;

  if (data.is_open()) {
    int i = 0;
    while (data.good()) {
      char linebuff[4096];
      getline(data, lineOfData);
      strncpy_s(linebuff, lineOfData.c_str(), sizeof(linebuff) - 1);
      {
        int j = 0;
        double val;
        char *p_r = NULL, *p_val;
        p_val = strtok_s(linebuff, " ,;", &p_r);
        while (NULL != p_val) {
          val = atof(p_val);
          matrix(i, j) = val;
          j++;
          p_val = strtok_s(NULL, " ,;", &p_r);
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
  X_xi = MatrixXd::Zero(M + 1, N + 1);
  Y_xi = MatrixXd::Zero(M + 1, N + 1);
  X_eta = MatrixXd::Zero(M + 1, N + 1);
  Y_eta = MatrixXd::Zero(M + 1, N + 1);
  X_xixi = MatrixXd::Zero(M + 1, N + 1);
  X_etaeta = MatrixXd::Zero(M + 1, N + 1);
  Y_xixi = MatrixXd::Zero(M + 1, N + 1);
  Y_etaeta = MatrixXd::Zero(M + 1, N + 1);
  J = MatrixXd::Zero(M + 1, N + 1);
  for (int i = 0; i < M; i++) {
    if (i > 0) {
      index = i - 1;
    } else {
      index = M - 1;
    }
    for (int j = 1; j < N; j++) {
      X_xi(i, j) = (X(i + 1, j) - X(index, j)) / 2.0 / d;
      Y_xi(i, j) = (Y(i + 1, j) - Y(index, j)) / 2.0 / d;
      X_eta(i, j) = (X(i, j + 1) - X(i, j - 1)) / 2.0 / d;
      Y_eta(i, j) = (Y(i, j + 1) - Y(i, j - 1)) / 2.0 / d;
      X_xixi(i, j) = (X(i + 1, j) - 2 * X(i, j) + X(index, j)) / d / d;
      Y_xixi(i, j) = (Y(i + 1, j) - 2 * Y(i, j) + Y(index, j)) / d / d;
      X_etaeta(i, j) = (X(i, j + 1) - 2 * X(i, j) + X(i, j - 1)) / d / d;
      Y_etaeta(i, j) = (Y(i, j + 1) - 2 * Y(i, j) + Y(i, j - 1)) / d / d;
      X_xieta(i, j) = (X(i + 1, j + 1) - X(i + 1, j - 1) - X(index, j + 1) +
                       X(index, j - 1)) /
                      4.0 / d / d;
      Y_xieta(i, j) = (Y(i + 1, j + 1) - Y(i + 1, j - 1) - Y(index, j + 1) +
                       Y(index, j - 1)) /
                      4.0 / d / d;
    }
  }

  for (int i = 0; i < M; i++) {
    if (i > 0) {
      index = i - 1;
    } else {
      index = M - 1;
    }
    int j = N;
    X_xi(i, j) = (X(i + 1, j) - X(index, j)) / 2.0 / d;
    Y_xi(i, j) = (Y(i + 1, j) - Y(index, j)) / 2.0 / d;
    X_eta(i, j) = (X(i, j) - X(i, j - 1)) / d;
    Y_eta(i, j) = (Y(i, j) - Y(i, j - 1)) / d;
    X_xixi(i, j) = (X(i + 1, j) - 2 * X(i, j) + X(index, j)) / d / d;
    Y_xixi(i, j) = (Y(i + 1, j) - 2 * Y(i, j) + Y(index, j)) / d / d;
    X_etaeta(i, j) = X_etaeta(i, j - 1);
    Y_etaeta(i, j) = Y_etaeta(i, j - 1);
    X_xieta(i, j) =
        (X(i + 1, j) - X(i + 1, j - 1) - X(index, j) + X(index, j - 1)) / 2.0 /
        d / d;
    Y_xieta(i, j) =
        (Y(i + 1, j) - Y(i + 1, j - 1) - Y(index, j) + Y(index, j - 1)) / 2.0 /
        d / d;
  }
  for (int i = 0; i < M; i++) {
    if (i > 0) {
      index = i - 1;
    } else {
      index = M - 1;
    }
    int j = 0;
    X_xi(i, j) = (X(i + 1, j) - X(index, j)) / 2.0 / d;
    Y_xi(i, j) = (Y(i + 1, j) - Y(index, j)) / 2.0 / d;
    X_eta(i, j) = (X(i, j + 1) - X(i, j)) / d;
    Y_eta(i, j) = (Y(i, j + 1) - Y(i, j)) / d;
    X_xixi(i, j) = (X(i + 1, j) - 2 * X(i, j) + X(index, j)) / d / d;
    Y_xixi(i, j) = (Y(i + 1, j) - 2 * Y(i, j) + Y(index, j)) / d / d;
    X_etaeta(i, j) = X_etaeta(i, j + 1);
    Y_etaeta(i, j) = Y_etaeta(i, j + 1);
    X_xieta(i, j) =
        (X(i + 1, j + 1) - X(i + 1, j) - X(index, j + 1) + X(index, j)) / 2.0 /
        d / d;
    Y_xieta(i, j) =
        (Y(i + 1, j + 1) - Y(i + 1, j) - Y(index, j + 1) + Y(index, j)) / 2.0 /
        d / d;
  }
  for (int i = 0; i < M; i++) {
    for (int j = 0; j <= N; j++) {
      double temp = X_xi(i, j) * Y_eta(i, j) - X_eta(i, j) * Y_xi(i, j);
      J(i, j) = temp;
    }
  }
}

void deriviate() {
  int index = 0;
  Vector3d A, B, C, D;
  Matrix3d E;
  double temp;
  Xi_x = MatrixXd::Zero(M + 1, N + 1);
  Eta_x = MatrixXd::Zero(M + 1, N + 1);
  Xi_y = MatrixXd::Zero(M + 1, N + 1);
  Eta_y = MatrixXd::Zero(M + 1, N + 1);
  Xi_xx = MatrixXd::Zero(M + 1, N + 1);
  Xi_xy = MatrixXd::Zero(M + 1, N + 1);
  Xi_yy = MatrixXd::Zero(M + 1, N + 1);
  Eta_xx = MatrixXd::Zero(M + 1, N + 1);
  Eta_xy = MatrixXd::Zero(M + 1, N + 1);
  Eta_yy = MatrixXd::Zero(M + 1, N + 1);
  for (int i = 0; i < M; i++) {
    if (i > 0) {
      index = i - 1;
    } else {
      index = M - 1;
    }
    for (int j = 0; j <= N; j++) {
      Xi_x(i, j) = Y_eta(i, j) / J(i, j);
      Eta_x(i, j) = -Y_xi(i, j) / J(i, j);
      Xi_y(i, j) = -X_eta(i, j) / J(i, j);
      Eta_y(i, j) = X_xi(i, j) / J(i, j);
      A << Xi_x(i, j) * X_xixi(i, j) + Xi_y(i, j) * Y_xixi(i, j),
          Xi_x(i, j) * X_xieta(i, j) + Xi_y(i, j) * Y_xieta(i, j),
          Xi_x(i, j) * X_etaeta(i, j) + Xi_y(i, j) * Y_etaeta(i, j);
      E << pow(Y_eta(i, j), 2), -2 * Y_eta(i, j) * Y_xi(i, j),
          pow(Y_xi(i, j), 2), -X_eta(i, j) * Y_eta(i, j),
          (X_xi(i, j) * Y_eta(i, j) + X_eta(i, j) * Y_xi(i, j)),
          -X_xi(i, j) * Y_xi(i, j), pow(X_eta(i, j), 2),
          -2 * X_xi(i, j) * X_eta(i, j), pow(X_xi(i, j), 2);
      temp = pow(X_xi(i, j) * Y_eta(i, j) - X_eta(i, j) * Y_xi(i, j), 2);
      B = -E * A / temp;
      Xi_xx(i, j) = B(0);
      Xi_xy(i, j) = B(1);
      Xi_yy(i, j) = B(2);
      C << Eta_x(i, j) * X_xixi(i, j) + Eta_y(i, j) * Y_xixi(i, j),
          Eta_x(i, j) * X_xieta(i, j) + Eta_y(i, j) * Y_xieta(i, j),
          Eta_x(i, j) * X_etaeta(i, j) + Eta_y(i, j) * Y_etaeta(i, j);
      D = -E * C / temp;
      Eta_xx(i, j) = D(0);
      Eta_xy(i, j) = D(1);
      Eta_yy(i, j) = D(2);
    }
  }
}

void init() {
  U = u0 * MatrixXd::Ones(M + 1, N + 1);
  U.row(M) = MatrixXd::Zero(1, N + 1);
  U.col(0) = MatrixXd::Zero(M + 1, 1);

  V = MatrixXd::Zero(M + 1, N + 1);
  PSI = MatrixXd::Zero(M + 1, N + 1);
  W = MatrixXd::Zero(M + 1, N + 1);
  for (int i = 0; i < M; i++) {
    for (int j = 1; j <= N; j++) {
      PSI(i, j) = 0;
    }
  }
  double sum = 0;
  for (int i = 0; i < M; i++) {
    sum += PSI(i, 1);
  }
  double ave = sum / M;
  for (int i = 0; i < M; i++) {
    PSI(i, 0) = ave;
  }
  double a = 0;
  double b = 0;
  double r = 0;
  double o = 0;
  double e = 0;
  double temp = 0;
  int index = 0;
  for (int i = 0; i < M; i++) {
    W(i, 0) = 2 * PSI(i, 1) * (pow(Eta_x(i, 0), 2) + pow(Eta_y(i, 0), 2));
  }
}
void push() {
  int index = 0;
  MatrixXd W_new(M + 1, N + 1);
  W_new = MatrixXd::Zero(M + 1, N + 1);
  double a = 0;
  double b = 0;
  double r = 0;
  double o = 0;
  double e = 0;
  double temp = 0;
  double b_w = 0;
  double b_e = 0;
  double b_s = 0;
  double b_n = 0;
  double b_p = 0;
  double c_p = 0;
  double b_w1 = 0;
  double b_e1 = 0;
  double b_s1 = 0;
  double b_n1 = 0;
  for (int i = 0; i < M; i++) {
    if (i > 0) {
      index = i - 1;
    } else {
      index = M - 1;
    }
    for (int j = 1; j < N; j++) {
      a = pow(Xi_x(i, j), 2) + pow(Xi_y(i, j), 2);
      b = Xi_x(i, j) * Eta_x(i, j) + Xi_y(i, j) * Eta_y(i, j);
      r = pow(Eta_x(i, j), 2) + pow(Eta_y(i, j), 2);
      o = Xi_xx(i, j) + Xi_yy(i, j);
      e = Eta_xx(i, j) + Eta_yy(i, j);
      b_e = (a + o / 2.0) / Re;
      b_w = (a - o / 2.0) / Re;
      b_s = (r + e / 2.0) / Re;
      b_n = (r - e / 2.0) / Re;
      b_w1 = (U(i, j) * Xi_x(i, j) + V(i, j) * Xi_y(i, j)) / 2.0;
      b_e1 = -b_w1;
      b_n1 = (U(i, j) * Eta_x(i, j) + V(i, j) * Eta_y(i, j)) / 2.0;
      b_s1 = -b_n1;
      W_new(i, j) =
          W(i, j) +
          dt * (W(i + 1, j) * (b_e + b_e1) + W(index, j) * (b_w + b_w1) +
                W(i, j + 1) * (b_s + b_s1) + W(i, j - 1) * (b_n + b_n1) -
                2.0 / Re * (a + r) * W(i, j) +
                2.0 * b / Re *
                    (W(i + 1, j + 1) - W(i + 1, j - 1) - W(index, j + 1) +
                     W(index, j - 1)));
      W_new(M, j) = W_new(0, j);
    }
  }

  W = W_new;
}

void psi_iteration() {
  int index = 0;
  double b_w = 0;
  double b_e = 0;
  double b_s = 0;
  double b_n = 0;
  double b_p = 0;
  double c_p = 0;
  double a = 0;
  double b = 0;
  double r = 0;
  double o = 0;
  double e = 0;
  double temp = 0;
  MatrixXd psi_new(M + 1, N + 1);
  psi_new = PSI;
  int count = 0;

  while (1) {
    for (int i = 0; i < M; i++) {
      if (i > 0) {
        index = i - 1;
      } else {
        index = M - 1;
      }
      for (int j = 1; j < N; j++) {
        a = pow(Xi_x(i, j), 2) + pow(Xi_y(i, j), 2);
        b = Xi_x(i, j) * Eta_x(i, j) + Xi_y(i, j) * Eta_y(i, j);
        r = pow(Eta_x(i, j), 2) + pow(Eta_y(i, j), 2);
        o = Xi_xx(i, j) + Xi_yy(i, j);
        e = Eta_xx(i, j) + Eta_yy(i, j);
        psi_new(i, j) =
            ((a - o / 2.0) * PSI(index, j) + (a + o / 2.0) * PSI(i + 1, j) +
             (r - e / 2.0) * PSI(i, j - 1) + (r + e / 2.0) * PSI(i, j + 1) +
             2 * b *
                 (PSI(i + 1, j + 1) - PSI(i + 1, j - 1) - PSI(index, j + 1) +
                  PSI(index, j - 1)) -
             W(i, j)) /
            (2 * a + 2 * r);
        psi_new(M, j) = psi_new(0, j);
      }
    }
    double error =
        (PSI.block(0, 1, M + 1, N - 1) - psi_new.block(0, 1, M + 1, N - 1))
            .norm();
    PSI = psi_new;

    std::cout << "psi error" << error << endl;
    count++;
    if (error < err) {
      return;
    }
  }
}
void velocity() {
  int index = 0;
  for (int i = 0; i < M; i++) {
    if (i > 0) {
      index = i - 1;
    } else {
      index = M - 1;
    }
    for (int j = 1; j < N; j++) {
      U(i, j) =
          X_xi(i, j) / J(i, j) * (PSI(i, j + 1) - PSI(i, j - 1)) / 2.0 / d -
          X_eta(i, j) / J(i, j) * (PSI(i + 1, j) - PSI(index, j)) / 2.0 / d;
      V(i, j) =
          Y_xi(i, j) / J(i, j) * (PSI(i, j + 1) - PSI(i, j - 1)) / 2.0 / d -
          Y_eta(i, j) / J(i, j) * (PSI(i + 1, j) - PSI(index, j)) / 2.0 / d;
    }
  }
  U.row(M) = U.row(0);
  V.row(M) = V.row(0);
}
void boundary() {
  int index = 0;
  int j = N;
  for (int i = 0; i < M; i++) {
    if (i > 0) {
      index = i - 1;
    } else {
      index = M - 1;
    }
    U(i, j) = u0;
    V(i, j) = 0;
    PSI(i, j) = u0 * Y(i, j);
    W(i, j) = 0;
  }
  double sum = 0;
  for (int i = 0; i < M; i++) {
    sum += PSI(i, 1);
  }
  double ave = sum / M;
  for (int i = 0; i < M; i++) {
    PSI(i, 0) = ave;
  }
  double a = 0;
  double b = 0;
  double r = 0;
  double o = 0;
  double e = 0;
  double temp = 0;

  for (int i = 0; i < M; i++) {
    W(i, 0) = 2 * PSI(i, 1) * (pow(Eta_x(i, 0), 2) + pow(Eta_y(i, 0), 2));
  }
}

void init_breakpoint(int index1) {
  char filename[30];
  sprintf_s(filename, "W%d.txt", index1);
  readPortfolio(W, filename);
  sprintf_s(filename, "U%d.txt", index1);
  readPortfolio(U, filename);
  sprintf_s(filename, "V%d.txt", index1);
  readPortfolio(V, filename);
  sprintf_s(filename, "PSI%d.txt", index1);
  readPortfolio(PSI, filename);
}

int main() {
  readPortfolio(X, "X.txt");
  readPortfolio(Y, "Y.txt");
  Jacobi();
  deriviate();
  ofstream xout("./output/X.data");
  xout << X;
  xout.close();
  ofstream yout("./output/Y.data");
  yout << Y;
  yout.close();
  ofstream jout("./output/J.data");
  jout << J;
  jout.close();

  init();
  PSI.row(M) = PSI.row(0);
  ofstream Uout("./output/U.data");
  Uout << U;
  Uout.close();
  ofstream Vout("./output/V.data");
  Vout << V;
  Vout.close();
  ofstream wout("./output/W.data");
  wout << W;
  wout.close();
  ofstream pout("./output/PSI.data");
  pout << PSI;
  pout.close();
  char filename[30];
  err = 0.01;

  for (int k = 1; k <= 300; k++) {
    cout << "k" << k << endl;
    push();
    sprintf_s(filename, "./output/W%d.data", k);
    ofstream wout(filename);
    wout << W;
    wout.close();
    if (k == 1) {
      ofstream pout("./output/PSI.data");
      pout << PSI;
      pout.close();
    }
    psi_iteration();
    velocity();
    boundary();

    if (k % 10 == 0) {
      sprintf_s(filename, "./output/PSI%d.data", k);
      ofstream psiout(filename);
      psiout << PSI;
      psiout.close();
      sprintf_s(filename, "./output/U%d.data", k);
      ofstream Uout(filename);
      Uout << U;
      Uout.close();
      sprintf_s(filename, "./output/V%d.data", k);
      ofstream Vout(filename);
      Vout << V;
      Vout.close();
    }
  }

  std::cout << "====== end ======\n";
  return 0;
}