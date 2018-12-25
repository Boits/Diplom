import org.knowm.xchart.QuickChart;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;

public class Solution {
    static double lambda = Math.sqrt(1 + 2 * Math.sqrt(7));
    static double b01 = lambda * lambda + 5;
    static double b02 = (lambda * lambda - 3) / lambda;
    static double b03 = (lambda * lambda + 3) / lambda;
    static double b04 = 2.0;
    static double z = 2.0;
    static double b1 = 1.0 / 2 * z;
    static double b2 = 1 - b1;
    static double eps = 0.00001;
    static double s = 1.0;

    static double P = 4.56 * Math.pow(10, -6);
    static double S = 400;
    static double m = 300;
    static double a = 5.95 * Math.pow(10, -5);
    static double alpha;

    static double norma(double q, double w, double e, double r) {
        return Math.sqrt(q * q + w * w + e * e + r * r);
    }
    //Функция опасности
    static double d1(double x1, double x2, double y1, double y2) {
        return b01 * (x1 - 1) + b02 * x2 + b03 * y1 + b04 * (y2 - 1);
    }
    //Правые части уравнений
    static double f1(double x1, double x2, double y1, double y2) {
        return x2 + y1;
    }
    static double f2(double x1, double x2, double y1, double y2) {
        return -x1 + y2;
    }
    static double f3(double x1, double x2, double y1, double y2) {
        return ((-3 / (Math.pow(norma(x1, x2, 0, 0), 3))) + 2) * x1 + y2 -(2.0 * P * S * Math.pow(Math.cos(alpha), 3)) / (a * m);
    }
    static double f4(double x1, double x2, double y1, double y2) {
        return ((-3 / (Math.pow(norma(x1, x2, 0, 0), 3))) - 1) * x2 - y1 -(2.0 * P * S * Math.pow(Math.cos(alpha), 2) * Math.sin(alpha)) / (a * m);
    }
    static double k11(double h, double x1, double x2, double y1, double y2) {
        return h * f1(x1, x2, y1, y2);
    }
    static double k21(double h, double x1, double x2, double y1, double y2) {
        return h * f2(x1, x2, y1, y2);
    }
    static double k31(double h, double x1, double x2, double y1, double y2) {
        return h * f3(x1, x2, y1, y2);
    }
    static double k41(double h, double x1, double x2, double y1, double y2) {
        return h * f4(x1, x2, y1, y2);
    }
    static double k12(double h, double x1, double x2, double y1, double y2) {
        return h * f1(x1 + z * k11(h, x1, x2, y1, y2), x2 + z * k21(h, x1, x2, y1, y2), y1 + z * k31(h, x1, x2, y1, y2), y2 + z * k41(h, x1, x2, y1, y2));
    }
    static double k22(double h, double x1, double x2, double y1, double y2) {
        return h * f2(x1 + z * k11(h, x1, x2, y1, y2), x2 + z * k21(h, x1, x2, y1, y2), y1 + z * k31(h, x1, x2, y1, y2), y2 + z * k41(h, x1, x2, y1, y2));
    }
    static double k32(double h, double x1, double x2, double y1, double y2) {
        return h * f3(x1 + z * k11(h, x1, x2, y1, y2), x2 + z * k21(h, x1, x2, y1, y2), y1 + z * k31(h, x1, x2, y1, y2), y2 + z * k41(h, x1, x2, y1, y2));
    }
    static double k42(double h, double x1, double x2, double y1, double y2) {
        return h * f4(x1 + z * k11(h, x1, x2, y1, y2), x2 + z * k21(h, x1, x2, y1, y2), y1 + z * k31(h, x1, x2, y1, y2), y2 + z * k41(h, x1, x2, y1, y2));
    }

    static double f(double alpha1) {
        return (2.0 * P * S * (b03 * (Math.pow(Math.cos(alpha1), 3)) + b04 * Math.pow(Math.cos(alpha1), 2) * Math.sin(alpha1)) / (a * m));
    }

    //Определение границ промежутка, в котором лежит альфа
    public static void Solution_1() {
        double alpha1[] = new double[500];
        double f[] = new double[500];
        int k = 0;
        double q = Math.PI / 2.0;
        for (double j = -q; j < q; j = j + 0.01) {
            alpha1[k] = j;
            k++;
        }
        k = 0;
        for (double j = -q; j < q; j = j + 0.01) {
            f[k] = f(j);
            k++;
        }
        //Поиск максимального значения альфа
        double f_max = -10;
        double alpha_max = -2;
        for (double j = -q; j < q; j = j + 0.00001) {
            if (f(j) > f_max) {
                f_max = f(j);
                alpha_max = j;
            }
        }
        System.out.println("f_max = " + f_max);
        System.out.println("alpha_max = " + alpha_max);
        //Поиск минимального значения альфа
        for (double alpha_ = -q; alpha_ < q; alpha_ = alpha_ + 0.000001) {
            alpha = alpha_;

            double pro_d1 = lambda * d1(1.01, 0, 0, 1) - f(alpha_);
            if ((Math.abs(pro_d1) < 0.000001)) {
                System.out.println("alpha0 = " + alpha_);
            }
        }
//График f(alpha) ========================================================
        double a[] = new double[k - 1];
        double b[] = new double[k - 1];

        for (int j = 0; j < k - 1; j++) {
            a[j] = alpha1[j];
        }
        for (int j = 0; j < k - 1; j++) {
            b[j] = f[j];
        }
        XYChart chart = QuickChart.getChart(" ", "alpha", "f(alpha)", "f(alpha)", a, b);
        new SwingWrapper(chart).displayChart();

        System.out.println("f(aplha0_1) = " + (f(-0.583)));
        System.out.println("f(aplha0_2) = " + (f(0.925)));

        chart.addSeries("alpha0_1", new double[]{-0.583}, new double[]{0.28});
        chart.addSeries("alpha0_2", new double[]{0.925}, new double[]{0.28});
        chart.addSeries("alpha_max", new double[]{0.168}, new double[]{0.792});
    }
    //Метод Рунге-Кутты
    public static void Solution_2() {
        //Выбор начального шага
        double w = Math.pow((1.0 / (1)), s + 1) + Math.pow(Math.abs(f1(0, 0, 0, 0)), s + 1);
        double h = Math.pow(eps / w, 1.0 / (s + 1));
        int n = (int) (Math.floor(2 * 1) / h);
        double[] y1 = new double[14000 * n];
        double[] y2 = new double[14000 * n];
        double[] y3 = new double[14000 * n];
        double[] y4 = new double[14000 * n];
        double[] p1 = new double[14000 * n];
        double[] p2 = new double[14000 * n];
        double[] p3 = new double[14000 * n];
        double[] p4 = new double[14000 * n];
        double[] t_graf1 = new double[14000 * n];
        y1[0] = 1.01;
        y2[0] = 0;
        y3[0] = 0;
        y4[0] = 1;
        for (double alpha_ = 0.168; alpha_ <= 1.2; alpha_ = alpha_ + 0.19) {
            alpha = alpha_;
            int i = 1, k = 1;
            for (double t = h; t <= 0.5; t = t + h) {
                t_graf1[k - 1] = t;
                double H = h;
                double y_1, y_2, y_3, y_4, y1a, y2a, y3a, y4a;
                y_1 = y1[i - 1] + (b2 * k11(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]) + b1 * k12(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]));
                y_2 = y2[i - 1] + (b2 * k21(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]) + b1 * k22(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]));
                y_3 = y3[i - 1] + (b2 * k31(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]) + b1 * k32(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]));
                y_4 = y4[i - 1] + (b2 * k41(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]) + b1 * k42(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]));

                y1a = y_1 + (b2 * k11(H / 2.0, y_1, y_2, y_3, y_4) + b1 * k12(H / 2.0, y_1, y_2, y_3, y_4));
                y2a = y_2 + (b2 * k21(H / 2.0, y_1, y_2, y_3, y_4) + b1 * k22(H / 2.0, y_1, y_2, y_3, y_4));
                y3a = y_3 + (b2 * k31(H / 2.0, y_1, y_2, y_3, y_4) + b1 * k32(H / 2.0, y_1, y_2, y_3, y_4));
                y4a = y_4 + (b2 * k41(H / 2.0, y_1, y_2, y_3, y_4) + b1 * k42(H / 2.0, y_1, y_2, y_3, y_4));

                y1[i] = y1[i - 1] + (b2 * k11(H, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]) + b1 * k12(H, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]));
                y2[i] = y2[i - 1] + (b2 * k21(H, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]) + b1 * k22(H, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]));
                y3[i] = y3[i - 1] + (b2 * k31(H, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]) + b1 * k32(H, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]));
                y4[i] = y4[i - 1] + (b2 * k41(H, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]) + b1 * k42(H, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]));
                //Оценка локальной погрешности
                p1[i] = (y1a - y1[i]) / (3.0);
                p2[i] = (y2a - y2[i]) / (3.0);
                p3[i] = (y3a - y3[i]) / (3.0);
                p4[i] = (y4a - y4[i]) / (3.0);
                //Автоматический выбор шага интегрирования
                if (norma(p1[i], p2[i], p3[i], p4[i]) > eps * 4.0) {
                    t = t - H / 2.0;
                    y1[i] = y1[i - 1] + (b2 * k11(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]) + b1 * k12(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]));
                    y2[i] = y2[i - 1] + (b2 * k21(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]) + b1 * k22(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]));
                    y3[i] = y3[i - 1] + (b2 * k31(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]) + b1 * k32(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]));
                    y4[i] = y4[i - 1] + (b2 * k41(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]) + b1 * k42(H / 2.0, y1[i - 1], y2[i - 1], y3[i - 1], y4[i - 1]));
                } else {
                    if ((norma(p1[i], p2[i], p3[i], p4[i]) > eps) && (norma(p1[i], p2[i], p3[i], p4[i]) <= eps * 4.0)) {
                        y1[i] = y1a;
                        y2[i] = y2a;
                        y3[i] = y3a;
                        y4[i] = y4a;
                        H = H / 2.0;
                    } else {
                        if ((norma(p1[i], p2[i], p3[i], p4[i]) >= eps / 8.0) && (norma(p1[i], p2[i], p3[i], p4[i]) <= eps)) {
                            y1[i] = y1[i];
                            y2[i] = y2[i];
                            y3[i] = y3[i];
                            y4[i] = y4[i];
                        } else {
                            H = 2.0 * H;
                        }
                    }
                }
                if ((t + H > 1) && (t < 1)) {
                    H = -t + 1;
                }
                h = H / 2;
                y1[i] = y1a + p1[i];
                y2[i] = y2a + p2[i];
                y3[i] = y3a + p3[i];
                y4[i] = y4a + p4[i];
                k = i;
                i++;
            }
            y1[k + 1] = y1[k] + p1[k];
            y2[k + 1] = y2[k] + p2[k];
//Графики ========================================================
            //Траектория движения КА
            double v1[] = new double[k-1];
            double v2[] = new double[k -1];
            for (int j = 0; j < k-1; j++) {
                v1[j] = y1[j];
                v2[j] = y2[j];
            }
            XYChart chart = QuickChart.getChart(" ", "x1", "x2", "alpha = " + (alpha_), v1, v2);
            new SwingWrapper(chart).displayChart();
            //Функция опасности
            double d11[] = new double[k-1];
            double[] t_graf = new double[k-1];
            for (int j = 0; j < k-1; j++) {
                d11[j] = d1((y1[j] + p1[j]),(y2[j] + p2[j]),(y3[j] + p3[j]),(y4[j] + p4[j]));
            }
            for(int j = 0; j < k-1; j++) {
                t_graf[j] = t_graf1[j];
            }
            XYChart chart1 = QuickChart.getChart(" ", "t", "d1", "alpha = " + (alpha_), t_graf, d11);
            new SwingWrapper(chart1).displayChart();
        }
    }
    public static void main(String[] args) {
        System.out.println("Поиск альфа: ");
        Solution_1();
        System.out.println("Траектории движения КА и графики функции опасности: ");
        Solution_2();
    }
}