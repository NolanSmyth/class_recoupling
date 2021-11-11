#include "class.h"
#include "stdio.h"

int main()
{
    double g = 1.0;
    double r = 2.0e3;
    double dm = 1.0e-8;
    double w = g * g / (2 * M_PI * r);
    struct Model model = {.dm = dm, .w = w, .r = r, .g = g};
    const int ws_size = 1000;
    int size = 500;
    double log_tmin = -4.0;
    double log_tmax = 8.0;
    double step = (log_tmax - log_tmin) / ((double)size);
    double T;
    double error = 0.0;
    double val1 = 0.0;
    double val2 = 0.0;
    const char fname[] = "data.csv";
    const char mode[] = "w+";
    FILE *file = fopen("data.csv", "w+");
    for (int i = 0; i < size; i++)
    {
        T = pow(10.0, log_tmin + i * step);
        val1 = thermal_scattering_rate(T, &model, &error);
        fprintf(file, "%f, %.16e\n", T, val1);
    }
    fclose(file);
    return 0;
}