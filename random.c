#include "graph.h"

/*Parisi-Rapuano [0,1) random number generator.*/

#define NormRANu (2.3283063671E-10F)
unsigned int  irr[256];
unsigned int  ir1;
unsigned char ind_ran, ig1, ig2, ig3;

double Random  (void)
{
    double r;

    ig1 = ind_ran - 24;
    ig2 = ind_ran - 55;
    ig3 = ind_ran - 61;
    irr[ind_ran] = irr[ig1] + irr[ig2];
    ir1 = (irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;

    return r;
}

void   ini_ran (int SEMILLA)
{
    int INI, FACTOR, SUM, i;

    srand(SEMILLA);

    INI    = SEMILLA;
    FACTOR = 67397;
    SUM    = 7364893;

    for (i = 0; i < 256; i++)
    {
        INI = (INI*FACTOR + SUM);
        irr[i] = INI;
    }

    ind_ran = ig1 = ig2 = ig3 = 0;
}

// copied from http://stackoverflow.com/a/10645091
double sampleNormal() {
        double u = Random() * 2 - 1;
        double v = Random() * 2 - 1;
        double r = u * u + v * v;
        if (r == 0 || r > 1) return sampleNormal();
        double c = sqrt(-2 * log(r) / r);
        return u * c;
}

int generateDegree (int m, double gamma, double norm_const)
{
        // min and max degree; we avoid k=0 because p~k^-GAMMA diverges
        double w, sum;
        int k;
        w = Random(); // random in [0,1)

        // Now we find x such that the sum from k=m to k=x of the area
        // of the probability distribution gives w. We are inverting the
        // quantile function.
        sum = 0;
        k   = m;
        while (sum < w) {
                sum += norm_const*pow(k,-gamma);
                k++;
        }

        if  ( k>sqrt(NODE_NR) || k<m) // one of the requirements is m<k<N^1/2
                return generateDegree(m,gamma,norm_const);
        return k;
}

