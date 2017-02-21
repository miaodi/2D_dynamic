#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <iomanip>
#include "spline.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/KroneckerProduct>

using namespace std;
using namespace spline;
using namespace Eigen;
typedef Eigen::SparseVector<double> SpVec;
typedef Eigen::SparseMatrix<double> SpMat;
const double Pi = 3.14159265358979323846264338327;

double
Deformation(double *knots_x, double *knots_y, double xi, double eta, MatrixXd DX, int m_x, int m_y, int p_x, int p_y);

void Geometry(double xi, double eta, double &pxpxi, double &pxpeta, double &pypxi, double &pypeta);

template<class T, class Fo>
T trapezoidal(T a, T b, Fo f, int n, T x, T y, T c0) {
    T h = (b - a) / n;
    T sum = f(x, y, b, a, c0) * .5;
    for (int i = 1; i < n; i++) sum += f(x, y, b, a + i * h, c0);
    sum += f(x, y, b, b, c0) * .5;
    return sum * h;
};

template<class T>
T f(T t) {
    if (t < 1 && t >= 0) {
        return 4 * (1 - pow(2 * t - 1, 2));
    } else {
        return 0;
    }
}

template<class T>
T Heaviside(T x) {
    if (x >= 0)
        return 1;
    return 0;
}

template<class T>
class integrand {
public:
    T operator()(T x, T y, T t, T tbar, T c0) {
        double a = Heaviside(c0 * (t - tbar) - sqrt(x * x + y * y)) / 2 / Pi / c0 /
                   sqrt(abs(pow(c0 * (t - tbar), 2) - x * x - y * y));
        double result = pow(c0, 2) * f(tbar) * a;
        return result;
    }
};

int main() {

    double length = 1;
    int order = 2;
    int p_x = order, p_y = order;
    double *knots_x, *knots_y;
    int elements_x = 18, elements_y = 18;
    int m_x = elements_x + 2 * p_x, m_y = elements_y + 2 * p_y;
    knots_x = new double[m_x + 1];
    for (int i = 0; i < p_x + 1; i++) {
        knots_x[i] = 0;
        knots_x[m_x - i] = 1;
    }
    for (int i = 0; i < elements_x; i++) {
        knots_x[i + p_x + 1] = (i + 1) * 1.0 / elements_x;
    }

    knots_y = new double[m_y + 1];
    for (int i = 0; i < p_y + 1; i++) {
        knots_y[i] = 0;
        knots_y[m_y - i] = 1;
    }
    for (int i = 0; i < elements_y; i++) {
        knots_y[i + p_y + 1] = (i + 1) * 1.0 / elements_y;
    }
    int dof_x = m_x - p_x, dof_y = m_y - p_y;
    int dof = dof_x * dof_y;

    double c0 = .1;
    double t = 0;
    double gamma = .5;
    double beta = 1.0 / 12;
    double h = length / elements_x;
    double tau = .25;
    double dt = tau * h / c0;
    double ratio_k = 6.0 * (9 * (3 - 40 * beta + 240 * pow(beta, 2)) * pow(c0 * dt, 4) + pow(h, 4)) /
                     (9 * (3 - 40 * beta + 240 * pow(beta, 2)) * pow(c0 * dt, 4) + 16 * pow(h, 4));
    double ratio_m = 6.0 *
                     (27 * pow(c0 * dt, 4) - 360 * beta * pow(c0 * dt, 4) + 2160 * pow(beta, 2) * pow(c0 * dt, 4) -
                      2 * pow(h, 4)) /
                     (27 * pow(c0 * dt, 4) - 360 * beta * pow(c0 * dt, 4) + 2160 * pow(beta, 2) * pow(c0 * dt, 4) -
                      32 * pow(h, 4));
    cout << ratio_m << " " << ratio_k << endl;
/*
    double gaussian[3];
    double weight[3];
    gaussian[0] = -pow(3 * (ratio_m - 2) / (ratio_m - 6), .5);
    gaussian[1] = 0;
    gaussian[2] = pow(3 * (ratio_m - 2) / (ratio_m - 6), .5);
    weight[0] = (ratio_m - 6) / (9 * (ratio_m - 2));
    weight[1] = 8.0 * (2 * ratio_m - 3) / (9 * (ratio_m - 2));;
    weight[2] = (ratio_m - 6) / (9 * (ratio_m - 2));
*/


    double *gaussian = x3;
    double *weight = w3;

    int gaussian_points = 3;
    MatrixXd M(dof, dof);
    MatrixXd K(dof, dof);
    VectorXd F(dof);
    M.setZero();
    K.setZero();
    F.setZero();
    for (int ii_x = 0; ii_x < elements_x; ii_x++) {
        double J_x = (knots_x[ii_x + p_x + 1] - knots_x[ii_x + p_x]) / 2;
        double Middle_x = (knots_x[ii_x + p_x + 1] + knots_x[ii_x + p_x]) / 2;
        int i_x = Findspan(m_x, p_x, knots_x, Middle_x);
        for (int ii_y = 0; ii_y < elements_y; ii_y++) {
            double J_y = (knots_y[ii_y + p_y + 1] - knots_y[ii_y + p_y]) / 2;
            double Middle_y = (knots_y[ii_y + p_y + 1] + knots_y[ii_y + p_y]) / 2;
            int i_y = Findspan(m_y, p_y, knots_y, Middle_y);
            MatrixXd M_ele;
            for (int jj_x = 0; jj_x < gaussian_points; jj_x++) {
                for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
                    double xi = Middle_x + J_x * gaussian[jj_x];
                    double eta = Middle_y + J_y * gaussian[jj_y];
                    double **ders_x, **ders_y;
                    DersBasisFuns(i_x, xi, p_x, knots_x, 1, ders_x);
                    DersBasisFuns(i_y, eta, p_y, knots_y, 1, ders_y);
                    VectorXd Nxi = VectorXd::Zero(p_x + 1), Nxi_xi = VectorXd::Zero(p_x + 1), Neta = VectorXd::Zero(
                            p_y + 1), Neta_eta = VectorXd::Zero(
                            p_y + 1);
                    for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
                        Nxi(kk_x) = ders_x[0][kk_x];
                        Nxi_xi(kk_x) = ders_x[1][kk_x];
                    }
                    for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
                        Neta(kk_y) = ders_y[0][kk_y];
                        Neta_eta(kk_y) = ders_y[1][kk_y];
                    }
                    for (int k = 0; k < 2; k++)
                        delete ders_x[k];
                    delete[] ders_x;
                    for (int k = 0; k < 2; k++)
                        delete ders_y[k];
                    delete[] ders_y;
                    VectorXd Nxi_xiNeta, NxiNeta_eta, NxiNeta;
                    Nxi_xiNeta = kroneckerProduct(Nxi_xi, Neta);
                    NxiNeta_eta = kroneckerProduct(Nxi, Neta_eta);
                    NxiNeta = kroneckerProduct(Nxi, Neta);
                    double pxpxi, pxpeta, pypxi, pypeta;
                    Geometry(xi, eta, pxpxi, pxpeta, pypxi, pypeta);
                    double Jacobian = pxpxi * pypeta - pxpeta * pypxi;
                    VectorXd Nx_xNy, NxNy_y;
                    Nx_xNy = 1.0 / Jacobian * (Nxi_xiNeta * pypeta - NxiNeta_eta * pypxi);
                    NxNy_y = 1.0 / Jacobian * (-Nxi_xiNeta * pxpeta + NxiNeta_eta * pxpxi);
                    M_ele = weight[jj_x] * weight[jj_y] * Jacobian * NxiNeta * NxiNeta.transpose() * J_x * J_y;
                    for (int iii_x = 0; iii_x < (p_x + 1) * (p_y + 1); iii_x++) {
                        for (int iii_y = 0; iii_y < (p_x + 1) * (p_y + 1); iii_y++) {
                            M(((i_x - p_x) * dof_y) + i_y - p_y + iii_x / (p_y + 1) * dof_y + iii_x % (p_y + 1),
                              ((i_x - p_x) * dof_y) + i_y - p_y + iii_y / (p_y + 1) * dof_y +
                              iii_y % (p_y + 1)) += M_ele(iii_x, iii_y);
                        }
                    }
                }
            }
        }
    }
/*
    gaussian[0] = -pow(3 * (ratio_k - 2) / (ratio_k - 6), .5);
    gaussian[1] = 0;
    gaussian[2] = pow(3 * (ratio_k - 2) / (ratio_k - 6), .5);
    weight[0] = (ratio_k - 6) / (9 * (ratio_k - 2));
    weight[1] = 8.0 * (2 * ratio_k - 3) / (9 * (ratio_k - 2));;
    weight[2] = (ratio_k - 6) / (9 * (ratio_k - 2));
*/
    for (int ii_x = 0; ii_x < elements_x; ii_x++) {
        double J_x = (knots_x[ii_x + p_x + 1] - knots_x[ii_x + p_x]) / 2;
        double Middle_x = (knots_x[ii_x + p_x + 1] + knots_x[ii_x + p_x]) / 2;
        int i_x = Findspan(m_x, p_x, knots_x, Middle_x);
        for (int ii_y = 0; ii_y < elements_y; ii_y++) {
            double J_y = (knots_y[ii_y + p_y + 1] - knots_y[ii_y + p_y]) / 2;
            double Middle_y = (knots_y[ii_y + p_y + 1] + knots_y[ii_y + p_y]) / 2;
            int i_y = Findspan(m_y, p_y, knots_y, Middle_y);
            MatrixXd M_ele;
            MatrixXd K_ele;
            for (int jj_x = 0; jj_x < gaussian_points; jj_x++) {
                for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
                    double xi = Middle_x + J_x * gaussian[jj_x];
                    double eta = Middle_y + J_y * gaussian[jj_y];
                    double **ders_x, **ders_y;
                    DersBasisFuns(i_x, xi, p_x, knots_x, 1, ders_x);
                    DersBasisFuns(i_y, eta, p_y, knots_y, 1, ders_y);
                    VectorXd Nxi = VectorXd::Zero(p_x + 1), Nxi_xi = VectorXd::Zero(p_x + 1), Neta = VectorXd::Zero(
                            p_y + 1), Neta_eta = VectorXd::Zero(
                            p_y + 1);
                    for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
                        Nxi(kk_x) = ders_x[0][kk_x];
                        Nxi_xi(kk_x) = ders_x[1][kk_x];
                    }
                    for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
                        Neta(kk_y) = ders_y[0][kk_y];
                        Neta_eta(kk_y) = ders_y[1][kk_y];
                    }
                    for (int k = 0; k < 2; k++)
                        delete ders_x[k];
                    delete[] ders_x;
                    for (int k = 0; k < 2; k++)
                        delete ders_y[k];
                    delete[] ders_y;
                    VectorXd Nxi_xiNeta, NxiNeta_eta, NxiNeta;
                    Nxi_xiNeta = kroneckerProduct(Nxi_xi, Neta);
                    NxiNeta_eta = kroneckerProduct(Nxi, Neta_eta);
                    NxiNeta = kroneckerProduct(Nxi, Neta);
                    double pxpxi, pxpeta, pypxi, pypeta;
                    Geometry(xi, eta, pxpxi, pxpeta, pypxi, pypeta);
                    double Jacobian = pxpxi * pypeta - pxpeta * pypxi;
                    VectorXd Nx_xNy, NxNy_y;
                    Nx_xNy = 1.0 / Jacobian * (Nxi_xiNeta * pypeta - NxiNeta_eta * pypxi);
                    NxNy_y = 1.0 / Jacobian * (-Nxi_xiNeta * pxpeta + NxiNeta_eta * pxpxi);
                    MatrixXd B = MatrixXd::Zero(2, (p_x + 1) * (p_y + 1));
                    for (int kk = 0; kk < (p_x + 1) * (p_y + 1); kk++) {
                        B(0, kk) = Nx_xNy(kk);
                        B(1, kk) = NxNy_y(kk);
                    }
                    K_ele = pow(c0, 2) * weight[jj_x] * weight[jj_y] * Jacobian * B.transpose() * B * J_x * J_y;
                    for (int iii_x = 0; iii_x < (p_x + 1) * (p_y + 1); iii_x++) {
                        for (int iii_y = 0; iii_y < (p_x + 1) * (p_y + 1); iii_y++) {
                            K(((i_x - p_x) * dof_y) + i_y - p_y + iii_x / (p_y + 1) * dof_y + iii_x % (p_y + 1),
                              ((i_x - p_x) * dof_y) + i_y - p_y + iii_y / (p_y + 1) * dof_y +
                              iii_y % (p_y + 1)) += K_ele(iii_x, iii_y);
                        }
                    }
                }
            }
        }
    }
    F(0) = pow(c0, 2) * f<double>(t) / 4;

    VectorXd d_current = VectorXd::Zero(dof);
    VectorXd v_current = VectorXd::Zero(dof);
    VectorXd a_current = VectorXd::Zero(dof);
    VectorXd d_pred = VectorXd::Zero(dof);
    VectorXd v_pred = VectorXd::Zero(dof);
    SpMat LHS = (M + beta * dt * dt * K).sparseView();
    ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;
    cg.compute(LHS);
    while (t < 8) {
        if (t + dt > 8) {
            dt = 8 - t;
        }
        t += dt;
        d_pred = d_current + dt * v_current + dt * dt / 2.0 * (1 - 2 * beta) * a_current;
        v_pred = v_current + (1 - gamma) * dt * a_current;
        F(0) = pow(c0, 2) * f<double>(t) / 4;
        VectorXd RHS = F - K * d_pred;
        a_current = cg.solve(RHS);
        d_current = d_pred + beta * dt * dt * a_current;
        v_current = v_pred + gamma * dt * a_current;
        /*
        int size = static_cast<int>(d_current.size());
        Matrix<Matrix<double, 1, 1>, Dynamic, 1> points(size);

        for (int i = 0; i < size; i++) {
            Matrix<double, 1, 1> point;
            point(0, 0) = d_current(i);
            points(i) = point;
        }
        NurbsCurve<double, 1> a(points, knotVector, 2);
        outf<<setprecision(10)<<"("<<t<<","<<a(.5)<<")"<<endl;
        outf1<<setprecision(10)<<"("<<t<<","<<trapezoidal(0.0,t, integrand<double>(),2000, .5)<<")"<<endl;
         */
        cout << t << " " << d_current(0) << " " << trapezoidal(0.0, t, integrand<double>(), 20000, 0.0, 0.000001, c0)
             << endl;
    }
    MatrixXd deform(dof_y, dof_x);
    for (int i = 0; i < dof_x; i++) {
        for (int j = 0; j < dof_y; j++) {
            deform(j, i) = d_current(j + i * dof_y);
        }
    }
    ofstream outf("20_regular.txt");
    ofstream outf1("exact.txt");
    outf << "x" << " " << "y" << " " << "z" << endl;
    for (int i = 0; i < 81; i++) {
        for (int j = 0; j < 81; j++) {
            outf << setprecision(10) << " " << i * 1.0 / 80 << " " << j * 1.0 / 80 << " "
                 << Deformation(knots_x, knots_y, i * 1.0 / 80, j * 1.0 / 80, deform, m_x,
                                m_y, p_x, p_y) << " " << endl;
        }
    }
    VectorXd exact1(201), exact2(201);
    for (int i = 0; i < 201; i++) {
        exact1(i) = trapezoidal(0.0, t - .01, integrand<double>(), 8000000, 0.0, i * 1.0 / 200, c0);
        exact2(i) = trapezoidal(0.0, t - .02, integrand<double>(), 8000000, 0.0, i * 1.0 / 200, c0);
    }
    for (int i = 0; i < 81; i++) {
        for (int j = 0; j < 81; j++) {
            outf1 << setprecision(10) << " " << i * 1.0 / 80 << " " << j * 1.0 / 80 << " "
                 << trapezoidal(0.0, t - .00001, integrand<double>(), 800000, i * 1.0 / 80, j * 1.0 / 80, c0) << " " << endl;
        }
    }
    return 0;
}

void Geometry(double xi, double eta, double &pxpxi, double &pxpeta, double &pypxi, double &pypeta) {
    double knot_x[] = {0, 0, 1, 1};
    double knot_y[] = {0, 0, 1, 1};
    MatrixXd B_x(2, 2);
    MatrixXd B_y(2, 2);
    MatrixXd weights(2, 2);
    weights << 1, 1, 1, 1;
    B_x << 0, 1, 0, 1;
    B_y << 0, 0, 1, 1;

    int p_x = 1, p_y = 1;
    int m_x = 3, m_y = 3;
    int dof_x = m_x - p_x, dof_y = m_y - p_y;
    int dof = dof_x * dof_y;
    MatrixXd B_xw, B_yw;
    B_xw = B_x.cwiseProduct(weights);
    B_yw = B_y.cwiseProduct(weights);
    int i_x = Findspan(m_x, p_x, knot_x, xi);
    int i_y = Findspan(m_y, p_y, knot_y, eta);
    double **ders_x, **ders_y;
    DersBasisFuns(i_x, xi, p_x, knot_x, 1, ders_x);
    DersBasisFuns(i_y, eta, p_y, knot_y, 1, ders_y);
    SpVec Nxi(dof_x), Nxi_xi(dof_x), Neta(dof_y), Neta_eta(dof_y);
    for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
        Nxi.insert(i_x - p_x + kk_x) = ders_x[0][kk_x];
        Nxi_xi.insert(i_x - p_x + kk_x) = ders_x[1][kk_x];
    }
    for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
        Neta.insert(i_y - p_y + kk_y) = ders_y[0][kk_y];
        Neta_eta.insert(i_y - p_y + kk_y) = ders_y[1][kk_y];
    }
    for (int k = 0; k < 2; k++)
        delete ders_x[k];
    delete[] ders_x;
    for (int k = 0; k < 2; k++)
        delete ders_y[k];
    delete[] ders_y;
    MatrixXd w, w_x, w_y;
    w = Neta.transpose() * weights * Nxi;
    w_x = Neta.transpose() * weights * Nxi_xi;
    w_y = Neta_eta.transpose() * weights * Nxi;
    MatrixXd pxpxi_temp, pxpeta_temp, pypxi_temp, pypeta_temp;
    pxpxi_temp = (Neta.transpose() * B_xw * Nxi_xi - Neta.transpose() * B_xw * Nxi * w_x(0, 0) / w(0, 0)) / w(0, 0);
    pypxi_temp = (Neta.transpose() * B_yw * Nxi_xi - Neta.transpose() * B_yw * Nxi * w_x(0, 0) / w(0, 0)) / w(0, 0);
    pxpeta_temp = (Neta_eta.transpose() * B_xw * Nxi - Neta.transpose() * B_xw * Nxi * w_y(0, 0) / w(0, 0)) / w(0, 0);
    pypeta_temp = (Neta_eta.transpose() * B_yw * Nxi - Neta.transpose() * B_yw * Nxi * w_y(0, 0) / w(0, 0)) / w(0, 0);
    pxpxi = pxpxi_temp(0, 0);
    pxpeta = pxpeta_temp(0, 0);
    pypxi = pypxi_temp(0, 0);
    pypeta = pypeta_temp(0, 0);
}

double
Deformation(double *knots_x, double *knots_y, double xi, double eta, MatrixXd DX, int m_x, int m_y, int p_x, int p_y) {
    int dof_x = m_x - p_x, dof_y = m_y - p_y;
    int dof = dof_x * dof_y;
    int i_x = Findspan(m_x, p_x, knots_x, xi);
    int i_y = Findspan(m_y, p_y, knots_y, eta);
    double **ders_x, **ders_y;
    DersBasisFuns(i_x, xi, p_x, knots_x, 0, ders_x);
    DersBasisFuns(i_y, eta, p_y, knots_y, 0, ders_y);
    SpVec Nxi(dof_x), Neta(dof_y);
    for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
        Nxi.insert(i_x - p_x + kk_x) = ders_x[0][kk_x];
    }
    for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
        Neta.insert(i_y - p_y + kk_y) = ders_y[0][kk_y];
    }
    for (int k = 0; k < 1; k++)
        delete ders_x[k];
    delete[] ders_x;
    for (int k = 0; k < 1; k++)
        delete ders_y[k];
    delete[] ders_y;
    MatrixXd x0 = Neta.transpose() * DX * Nxi;
    double position = x0(0, 0);
    return position;
}