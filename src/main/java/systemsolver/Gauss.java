package systemsolver;

public class Gauss {
    private int n;
    private double[][] a;
    private double[] b;

    public Gauss(double[][] a, double[] b) {
        this.a = a;
        this.b = b;
        n = b.length;
    }

    public double[] getX() {
        double[] x = new double[n];
        double[][] m = new double[n][n];
        double sum;
        for (int k = 0; k < n - 1; k++) {
            double max = a[k][k];
            int maxi = k;
            for (int i = k + 1; i < n; i++) {
                if (a[i][k] > max) {
                    max = a[i][k];
                    maxi = i;
                }
            }
            double temp;
            temp = b[k];
            b[k] = b[maxi];
            b[maxi] = temp;
            for (int j = k; j < n; j++) {
                temp = a[k][j];
                a[k][j] = a[maxi][j];
                a[maxi][j] = temp;
            }
            for (int i = k + 1; i < n; i++) {
                m[i][k] = -(a[i][k] / a[k][k]);
                b[i] = b[i] + m[i][k] * b[k];
                for (int j = k + 1; j < n; j++) {
                    a[i][j] = a[i][j] + m[i][k] * a[k][j];
                }
            }
        }
        x[n - 1] = b[n - 1] / a[n - 1][n - 1];
        for (int k = n - 2; k > -1; k--) {
            sum = 0;
            for (int j = k + 1; j < n; j++) {
                sum += a[k][j] * x[j];
            }
            x[k] = (b[k] - sum) / a[k][k];
        }
        return x;
    }

    public void print() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                System.out.print(this.a[i][j] + " ");
            }
            System.out.println(this.b[i]);
        }
        System.out.println();
    }
}
