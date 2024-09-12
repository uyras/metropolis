class My_gsl_vector
{
    gsl_vector* data;
public:



    My_gsl_vector(int n_)
        : data(gsl_vector_alloc(n_)) {}
    ~My_gsl_vector()
    {
        gsl_vector_free(data);
    }

    operator gsl_vector* () const { return data; }
};

class My_gsl_matrix
{
    gsl_matrix* data;

public:
    My_gsl_matrix(int n_, int l_)
        : data(gsl_matrix_alloc(n_, l_)) {}
    ~My_gsl_matrix()
    {
        gsl_matrix_free(data);
    }

    operator gsl_matrix* () const { return data; }
};

double chisq_compute(std::vector<double> x)
{
    double chisq = 0;
    for (int i = 0; i < 4; i++) chisq += x[i] * x[i];
    return chisq;
}

double f(double X, std::vector<double> A, int n, int order)
{
    double e_sum = 0;

    for (int i = 0; i < order + 1; i++) e_sum += A[i] * pow(X, i);

    return e_sum;

}

std::vector<double> gsl_polynomial_fit(std::vector<double>& data_x,
    std::vector<double>& data_y, const int n, const int order)
{
    My_gsl_vector y(n), c(order + 1);
    y.operator gsl_vector* (); c.operator gsl_vector* ();

    My_gsl_matrix X(n, order + 1), cov(order + 1, order + 1);
    X.operator gsl_matrix* (); cov.operator gsl_matrix* ();

    double chisq = chisq_compute(data_x);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < order + 1; j++) {
            gsl_matrix_set(X, i, j, pow(data_x[i], j));
        }
        gsl_vector_set(y, i, data_y[i]);
    }


    gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(n, order + 1);
    gsl_multifit_linear(X, y, c, cov, &chisq, work);
    gsl_multifit_linear_free(work);

    std::vector<double> vc;
    for (int i = 0; i < order + 1; i++) {
        vc.push_back(gsl_vector_get(c, i));
    }

    /*
    gsl_vector_free(y);
    gsl_vector_free(c);
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    */

    return vc;
}

std::vector<double> interpolation_temperatures_(std::vector<double> temperatures, std::vector<double> energies, int N)
{
    std::sort(temperatures.begin(), temperatures.end());
    std::sort(energies.begin(), energies.end());
   
    double p;
    int order = (N > 6) ? 6 : N - 1;

    std::vector<double> coefs = gsl_polynomial_fit(temperatures, energies, N, order);

    for (int i = 0; i < N - 1; i++)
    {

        double dT = (temperatures[i + 1] - temperatures[i]) / 2;
        p = std::pow(2.718282, (((double)energies[i + 1] - (double)energies[i]) * ((1 / (double)temperatures[i + 1])
            - 1 / (double)(temperatures[i]))));
        int j = i + 1;
        int n = 0;
        int flag = 0;


        while ((std::abs(p - 0.2) > 0.001))
        {
            n += 1;
            if (p > 0.2)
            {
                if (flag == -1) dT /= 2;
                if (flag == 1) dT *= 1.5;
                flag = 1;
                for (int k = j; k < N; k++) temperatures[k] += dT;
            }
            else
            {
                if (flag == 1) dT /= 2;
                if (flag == -1) dT *= 1.5;
                flag = -1;
                for (int k = j; k < N; k++) temperatures[k] -= dT;
            }
            if (n == 1000) break;
            for (int k = j; k < N; k++) energies[k] = f(temperatures[k], coefs, N, order);

            p = pow(2.718282, (((double)energies[i + 1] - (double)energies[i]) * ((1 / (double)temperatures[i + 1])
                - 1 / (double)(temperatures[i]))));
        }
    }

    std::sort(temperatures.begin(), temperatures.end());

    return temperatures;

}
