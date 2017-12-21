package nystrom;

import static java.lang.Math.*;

import area.Area;
import systemsolver.Gauss;

public class newnestrom {

    private int n;
    private double[][] m;
    private double[] p;
    private double[] t;
    private Area area;
    private double alpha;

    public newnestrom(int n, Area area, double alpha) {
        this.n = n;
        m = new double[4 * n][4 * n];
        p = new double[4 * n];
        t = new double[2 * n];
        this.area = area;
        this.alpha = alpha;
        for (int i = 0; i < 2 * n; i++) {
            t[i] = i * PI / n;
        }
    }

    private void createSystem() {
        for (int i = 0; i < 2 * n; i++) {
            for (int j = 0; j < 2 * n; j++) {
                m[i][j] = K1(t[i], t[j]) * mod(area.a(t[i])) / (2 * n);
                if (i == j) {
                    m[i][j] -= 0.5/mod(area.aD(t[i]));
                }
            }
            for (int j = 2 * n; j < 4 * n; j++) {
                m[i][j] = K2(t[i], t[j - 2 * n]) * mod(area.b(t[i])) / (2 * n);
            }
            p[i] = f(t[i]);
        }
        for (int i = 2 * n; i < 4 * n; i++) {
            for (int j = 0; j < 2 * n; j++) {
                m[i][j] = (K3(t[i - 2 * n], t[j]) + alpha * K5(t[i - 2 * n], t[j]))* mod(area.a(t[i - 2 * n])) / (2 * n);
            }
            for (int j = 2 * n; j < 4 * n; j++) {
//                m[i][j] = (K4(t[i - 2 * n], t[j - 2 * n]) / (2 * n) + alpha * (K6(t[i - 2 * n], t[j - 2 * n]) / (2 * n) -
//                        K4(t[i - 2 * n], t[j - 2 * n]) * R(j - 2 * n, t[i - 2 * n])));
//                m[i][j] = (K4(t[i - 2 * n], t[j - 2 * n]) / (2 * n) + alpha * (K61(t[i - 2 * n], t[j - 2 * n]) / (2 * n) +
//                        K62(t[i - 2 * n], t[j - 2 * n]) * R(j - 2 * n, t[i - 2 * n])));
                m[i][j] = (K4(t[i - 2 * n], t[j - 2 * n]) / (2 * n) + alpha * (L2(t[i - 2 * n], t[j - 2 * n]) / (2 * n) +
                        L1(t[i - 2 * n], t[j - 2 * n]) * R(j - 2 * n, t[i - 2 * n])))* mod(area.b(t[i - 2 * n]));
                if (i == j) {
                    m[i][j] -= 0.5/mod(area.bD(t[i - 2 * n]));
                }
            }
            p[i] = g(t[i - 2 * n]);
        }
    }

    private double[] getPhi() {
        createSystem();
        Gauss gauss = new Gauss(m, p);
        return gauss.getX();
    }

    public double getSolution() {
        double[] phi = getPhi();
        double u = 0;
        for (int i = 0; i < 2 * n; i++) {
            u += phi[i] * K1Ex(t[i]) + phi[i + 2 * n] * K2Ex(t[i]);
        }
        u = u / (2 * n);
        return u;
    }

    public double getExactSolution() {
//        double[] u = new double[2 * n];
//        double[] y = yOut();
//        for (int i = 0; i < 2 * n; i++) {
//            double r1 = r(area.a(t[i]), y);
//            double r2 = r(area.b(t[i]), y);
//            u[i] = -log(r1) - log(r2);
//        }
//        return u;
        return -log(r(xIn(), yOut()))/(2*PI);
    }

    private double K1Ex(double t) {
        double[] x = xIn();
        double r = r(area.a(t), x);
        return ((x[0] - area.a1(t)) * area.a2D(t) - (x[1] - area.a2(t)) * area.a1D(t)) / (r * r * mod(area.aD(t)));
    }

    private double K2Ex(double t) {
        double[] x = xIn();
        return -log(r(x, area.b(t)));
    }

    private double K1(double t, double T) {
        double K;
        double mod = mod(area.aD(t));
        if (t != T) {
            double r = r(area.a(t), area.a(T));
            K = ((area.a1(t) - area.a1(T)) * area.a2D(T) - (area.a2(t) - area.a2(T)) * area.a1D(T)) / (r * r * mod);
        } else {
            K = (area.a1D(t) * area.a2DD(t) - area.a2D(t) * area.a1DD(t)) / (2 * mod * mod * mod);
        }
        return K;
    }

    private double K2(double t, double T) {
        return -log(r(area.a(t), area.b(T)));
    }

    private double K3(double t, double T) {
        double r = r(area.b(t), area.a(T));
        double K1;
        double K2;
        K1 = ((area.b1(t) - area.a1(T)) * area.b2D(t) - (area.b2(t) - area.a2(T)) * area.b1D(t)) / (mod(area.bD(t)));
        K2 = ((area.b1(t) - area.a1(T)) * area.a2D(T) - (area.b2(t) - area.a2(T)) * area.a1D(T)) / (mod(area.aD(T)));
        return (area.b2D(t) * area.a2D(T) + area.b1D(t) * area.a1D(T)) / ( r * r * mod(area.bD(t)) * mod(area.aD(T))) -
                2 * K1 * K2 / (pow(r, 4));
//        return  -K1 * K2;
    }


    ////
    private double K4(double t, double T) {
        double K;
        double mod = mod(area.bD(t));
        if (t != T) {
            double r = r(area.b(t), area.b(T));
            K = -((area.b1(t) - area.b1(T)) * area.b2D(t) - (area.b2(t) - area.b2(T)) * area.b1D(t)) / (r * r * mod);
        } else {
            K = (area.b2D(t) * area.b1DD(t) - area.b1D(t) * area.b2DD(t)) / (2 * mod * mod * mod);
        }
        return K;
    }

    private double K5(double t, double T) {
        double r = r(area.b(t), area.a(T));
        double mod = mod(area.aD(t));
        return ((area.b1(t) - area.a1(T)) * area.a2D(T) - (area.b2(t) - area.a2(T)) * area.a1D(T)) / (r * r * mod);
    }

    private double L2(double t, double T) {
        if (t != T) {
            double r = r(area.b(t), area.b(T));
            return log(4 * pow(sin((t - T) / 2), 2) / (r * r)) / 2;
        } else {
            return log(1 / pow(mod(area.bD(T)), 2)) / 2;
        }

    }

    private double L1(double t, double T) {
        return -0.5;
    }

    double R(int j, double t) {
        double sum = 0;
        for (int i = 1; i <= n - 1; i++) {
            sum = sum + cos(i * (t - this.t[j])) / i;
        }
        sum = sum + cos(n * (t - this.t[j])) / (2 * n);
        return -sum / n;
    }

    private double[] xIn() {
        return new double[]{2, 0};
    }

    private double[] yOut() {
        return new double[]{5, 5};
    }

    private double f(double t) {
        return -log(r(area.a(t), yOut()))/(2*PI);
//        return 1;
    }

    private double g(double t) {
        double[] y = yOut();
//        double r = r(area.b(t), y);
//        return -((area.b1(t) - y[0]) * area.b2D(t) - (area.b2(t) - y[1]) * area.b1D(t)) / (r * r * mod(area.bD(t)))
//                - alpha * log(r);
        double K;
        double mod = mod(area.aD(t));
//        if (t != T) {
        double r = r(area.a(t), y);
        K = ((area.b1(t) - y[0]) * area.b2D(t) - (area.b2(t) - y[1]) * area.b1D(t)) / (r * r * mod) - alpha * log(r);
//        } else {
//            K = (area.a1D(t) * area.a2DD(t) - area.a2D(t) * area.a1DD(t)) / (2 * mod * mod * mod);
//        }
        return K/(2*PI);
//        return 1;
    }

    private double mod(double[] x) {
        return sqrt(pow(x[0], 2) + pow(x[1], 2));
    }

    private double[] diff(double[] x, double[] y) {
        return new double[]{x[0] - y[0], x[1] - y[1]};
    }

    private double r(double[] x, double[] y) {
        return mod(diff(x, y));
    }

    public void print() {
        int n = p.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                System.out.print(m[i][j] + " ");
            }
            System.out.println(p[i]);
        }
        System.out.println();
    }

}
