#include <stdlib.h>
#include <stdio.h>

int main(void)
{

    FILE *fp;

    float fx[210];
    float fy[210];
    int i = 0;

    fp = fopen("interps/resonant_rate.csv", "r");

    if (fp == NULL)
    {
        printf("Error opening file");
        return 1;
    }

    while (fscanf(fp, "%f", &fx[i]) != EOF)
    {
        fscanf(fp, ",");
        fscanf(fp, "%f", &fy[i]);
        // printf("%e %e\n", fx[i], fy[i]);

        i++;
    }
    fclose(fp);

    for (int j = 0; j < i; j++)
    {
        printf("%e %e\n", fx[j], fy[j]);
    }

    return 0;
}