#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

int main(void)
{

    int data_size = 208; // this must be exactly the number of lines of the csv file
    FILE *fp;

    double fx[data_size];
    double fy[data_size];
    int i = 0;

    fp = fopen("interps/resonant_rate.csv", "r");

    if (fp == NULL)
    {
        printf("Error opening file");
        return 1;
    }

    while (fscanf(fp, "%lf", &fx[i]) != EOF)
    {
        fscanf(fp, ",");
        fscanf(fp, "%lf", &fy[i]);
        // printf("%e %e\n", fx[i], fy[i]);

        i++;
    }
    fclose(fp);

    double xi;
    double yi;

    {
        gsl_interp_accel *acc = gsl_interp_accel_alloc();
        gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, data_size);

        gsl_spline_init(spline, fx, fy, data_size);

        xi = 1000.1;
        yi = gsl_spline_eval(spline, xi, acc);
        printf("%g %g\n", xi, yi);

        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }
    return 0;
}