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
        if (r < 1e-10 || r > 1) return sampleNormal();  // I want r==0 but r is
                                                        // float so I use 1e-10
        double c = sqrt(-2 * log(r) / r);
        return u * c;
}

int generateDegree (int m, double gamma)
{
    // min and max degree; we avoid k=0 because p~k^-GAMMA diverges
    double norm, sum_norm, P, r;
    int k = 0;

    norm = sum_norm = P = r = 0;

    int k1, k_min, k_max;
    k_min = K_MIN;  // to generate connected net with prob 1,
                    // see PHYS.REVIEW E71,027103(2005)
    k_max = (int)(NODE_NR-1);
/*
    k = 2*NODE_NR;
    //while  ( k>sqrt(NODE_NR) || k<m) 
    {
        w = Random();
        k = (int) ( ((double)k_min-0.5)*pow(1-w,-1/(gamma-1)) + 0.5 );
    }
*/

    // normalization constant of p~k^-GAMMA
    sum_norm = 0;
    for (int k_idx = k_min; k_idx <= k_max ; k_idx++) {
        sum_norm += pow(k_idx,-GAMMA);
    }   
    norm = 1.0/sum_norm;

    k = -1; // this makes sure we get inside the loop, as m is posiitve
    while(  k<m || k>5*sqrt(NODE_NR) )
    {   
        P = 0;
        k1 =  k_min;
        for(int i_idx = k1; i_idx <= k_max; i_idx++) {
            P += norm*pow(i_idx,-gamma);
        }
        r = Random();
        while ( P > r )
        {
            //printf("bad loop4! P=%g 1-r=%g\n", P, 1-r);
            k1++;
            P = 0;
            for(int i_idx = k1; i_idx <= k_max; i_idx++) {
                P += norm*pow(i_idx,-gamma);
            }
        }
        k = k1-1;
        //printf("bad loop4! k=%d P=%g 1-r=%g\n", k, P, 1-r);
    }
    //printf("bad loop4! k=%d P=%g 1-r=%g\n", k, P, 1-r);
/*
    // Now we find x such that the sum from k=m to k=x of the area
    // of the probability distribution gives w. We are inverting the
    // quantile function.
    double sum, sum2, w;
    sum = 0;
    sum2= 0;
    k   = m;
    while (sum < w) 
    {
        sum += norm*pow(k,-gamma);
        sum2 = sum + norm*pow(k,-gamma);
        k++;
        if (sum2 > w) {
            if (w>sum && w < (sum2+sum)/2)  
            {
                k--;
                break;
            } 
            else
            {
                break;
            }
        }
    }
*/    

//    if  ( k>sqrt(NODE_NR) || k<m) // one of the requirements is m<k<N^1/2
//        return generateDegree(m,gamma,norm_const);

    return k;
}


// Returns a  to which to connect to proportionally to the degree
// WARNING: nodes must be length m, have it initialized to -1
void generate_node_BA (int m, int* nodes)
{
    int counter = 0;
    double w = 0;
    double sum = 0;
    double sum_norm = 0;

    for(int i = 0; i < m; i++) {
        nodes[i] = -1;
    }
    
    double *tags = calloc(NODE_NR, sizeof *tags);
    for(int i = 0; i < NODE_NR; i++) {
        tags[i] = degree[i];
    }

    for(int m_idx = 0; m_idx < m; m_idx++)
    {
        w = Random();

        sum_norm = 0;
        for(int i_idx = 0; i_idx < NODE_NR; i_idx++) {
            sum_norm += tags[i_idx];
        }

        counter = 0;
        sum = (double)tags[counter];
        while(sum < w*sum_norm){
            counter++;
            sum += (double)tags[counter];
        }
        //printf("sum: %lf counter: %d sum_norm: %d\n", sum, counter, sum_norm);

        nodes[m_idx] = counter;
        tags[counter] = 0;
    }

    for(int i = 0; i < m; i++)
    {
        if (nodes[i]<0 || nodes[i]>NODE_NR) {
            printf("error: something went wrong with BA node generator\n");
            exit(12);
        }
        
    }
    
    
    free(tags); tags = NULL;
}

