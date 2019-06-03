#include "graph.h"

//Parisi-Rapuano [0,1) random number generator.

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


void calcHist(char* filename)
{
    double min, max, d;
    int N_Intervalos = NODE_NR-1;
    int N_data = NODE_NR;
    double *Hist;
    Hist = (double*)calloc(N_Intervalos, sizeof *Hist);
    Histogram(degree, Hist, N_data, N_Intervalos, &d, &min, &max);
    #ifdef DEBUG
    printf("%s%g\n", "min=", min);
    printf("%s%g\n", "max=", max);
    printf("%s%g\n", "d=",   d);
    printf("%s\n", "Histogram calculated.");
    #endif

    //for(int i = 0; i < NODE_NR; i++) {
    //    printf("%lf\t", Hist[i]);
    //}   printf("\n");

    // Writing to file.
    FILE *file_hist;
    file_hist = fopen(filename, "w");
    if (file_hist != NULL) {
        #ifdef DEBUG
        printf("%s\n", "file_hist opened.");
        #endif
        fprintf(file_hist, "%s\t%s\n", "# Pasos","Frecuencia_relativa");
        for (int i = 0; i < N_Intervalos; i++) {
            #ifdef DEBUG
            printf("%s\n", "He escrito una linea en hist.txt!");
            #endif
            fprintf(file_hist, "%g %g\n", min+i*d, Hist[i]);
        }
    } else {
        printf("%s\n", "Something went wrong opening the histogram file.");
        exit(5);
    }
    fclose(file_hist);

}

void med_var (double *datos, int numero_de_datos,
              double *media, double *varianza)
{
    // Para las sumas.
    double sum_1, sum_2;
    sum_1 = sum_2 = 0;
    for (int i = 0; i < numero_de_datos; i++) {
            sum_1 += datos[i];
            sum_2 += datos[i]*datos[i];
    }
    *media     = sum_1 / numero_de_datos;
    *varianza  = sum_2 / numero_de_datos - (*media) * (*media);
}

void Histogram (int *data, double *Hist, int N_data, int N_Intervalos,
                double *d, double *min, double *max)
{
    // INFORMACION
    // *data        -- input,  dato
    // *Hist        -- output, histograma
    // N_data       -- input,  numero de datos
    // N_Intervalos -- input,  numero de divisiones/intervalos
    // *d           -- output, tamaño de cada intervalo
    // *min         -- output, valor minimo de los datos
    // *max         -- output, valor maximo de los datos

    // Calculo del minimo y maximo.
    // Inicializo los valores como me da la gana.
    double min_aux = +1e40; double max_aux = -1e40;

    for (int i = 0; i < N_Intervalos; ++i) {
        Hist[i] = 0;
    }

    for (int i = 0; i < N_data; i++)
    {
            if ( (double) data[i] < min_aux ) min_aux = data[i];
            if ( (double) data[i] > max_aux ) max_aux = data[i];
    }
    *max = max_aux;
    *min = min_aux;
    
    // Specifically for this type of networks:
    *min = 2;
    *max = NODE_NR;
    N_Intervalos = NODE_NR-2-1;

    // Calculo la longitud de cada intervalo una vez que
    // sabemos el maximo y minimo.
    *d = (*max-*min)/((double)N_Intervalos);
    if ( *d < 1e-10 ) { // i want d==0, but d is  afloar
        printf("%s\n", "La longitud del intervalo es nula.");
        exit(6);
    }

    #ifdef DEBUG
    printf("tamagnoINTERVALO: %g\n", *d);
    #endif

    //********** EL PASO MAS IMPORTANTE *************//
    for (int i = 0; i < N_data; i++){
        Hist[ (int)( ((double)(data[i])-*min)/ (*d)) ]++;
    }

    // Ahora normalizo.
    double factorNORMALIZANTE = (1.0/ ((*d) * N_data));

    #ifdef DEBUG
    printf("factorNORMALIZANTE: %g\n", factorNORMALIZANTE);
    #endif

    for (int i = 0; i < N_Intervalos; i++) {
        Hist[i]*=factorNORMALIZANTE;
    }

    #ifdef DEBUG
    printf("%s\n", "Histogram calculated.");
    #endif
}

int delta_kron (int i, int j)
{
    if (i==j) {
        return 1;
    } else
    {
        return 0;
    }
}


/*
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
*/

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
        if (r < 1e-6 || r > 1) return sampleNormal();  // I want r==0 but r is
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
    
    double *tags = (double*)calloc(NODE_NR, sizeof *tags);
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

}

int generate_node_FREQUENCY_GAP (double alpha, int *node_i, double *node_i_freq, 
                                int m, int *nodes, double *nodes_freq,
                            double t, double sigma, double *weff, int len_weff)
{
    int counter = 0;
    double w = 0;
    double sum = 0;
    double sum_norm = 0;

    for(int i = 0; i < m; i++) {
        nodes[i] = -1;
    }

    double tags[NODE_NR];
    int aux_name = 0;
    double tmp = 0;

/*
    double r_dom[NODE_NR];
    double weff_by_node[NODE_NR];
    double weff_by_domain[NODE_NR]; // initialized in weff_compt_efficient
    //for (size_t i = 0; i < NODE_NR; i++) {
    //    weff_by_domain[i] = FG_WEFF_LOWER_FREQUENCY;
    //}
    
    
    weff_compt_DOUBLY_efficient(weff_by_domain, weff_by_node, r_dom, NODE_NR, t, sigma);
    //for (int i = 0; i < NODE_NR; i++) {
    //    printf("%g ", weff_by_domain[i]);
    //}   printf("\n");
    
    //double rf = 0;

    double weff[NODE_NR];
    for (int i = 0; i < NODE_NR; i++) {
        //aux_name = GLOB_component_name[i];
        //weff[i] = weff_by_domain[aux_name];
        weff[i] = weff_by_node[i];
    }
*/

    double norm_omega_eff = 0;
    for (int i = 0; i < NODE_NR; i++) {
        norm_omega_eff += exp(FG_ALPHA*fabs(weff[i]));
    }

    double aux_rnd22 = Random();
    int iidx = 0;
    aux_name = GLOB_component_name[iidx];
    double aux_sum = exp(FG_ALPHA*fabs(weff[iidx]));
    //printf("glob sum: %g\n", GLOB_sum_omega_nat);
    while (aux_sum <= aux_rnd22*norm_omega_eff) {
        iidx++;
        aux_sum += exp(FG_ALPHA*fabs(weff[iidx]));
    }
    *node_i = iidx;
    *node_i_freq = weff[*node_i];
    //printf("first node: %d \n", iidx);

    double wi = weff[*node_i];
    double wf = 0;
    for(int i = 0; i < NODE_NR; i++) {
        wf = weff[i];
        //rf = r_dom[aux_name];
        
        
        //if (wf<(FG_WEFF_LOWER_FREQUENCY/10)) {
        //    printf("warning: something went wrong with FG generator\n");
        //    exit(123);
        //}
        tmp = diff_weff_weight(alpha, wi, wf, 0, 0);
        //printf("wi,wf=%g,%g\n",wi,wf);
        //printf("tmp:%g\n",tmp);
        tags[i] = tmp;
        //printf("%d tags[%d]=%lf\n",*node_i,i,tags[i]);
    }

    // Condition that neighbors have prob zero of being chosen
    // this is so that 
    for(int j = 0; j < degree[*node_i]; j++) {
        tags[ C[*node_i][j] ] = 0;
    }
    tags[*node_i] = 0; // *node_i also has 0 probability of being chosen
    //for (size_t i = 0; i < NODE_NR; i++){
    //    printf("%lf ", tags[i]);
    //}   printf("\n");
    
    sum_norm = 0;
    for(int i_idx = 0; i_idx < NODE_NR; i_idx++) {
        sum_norm += tags[i_idx];
    }
    //printf("sum_norm=%g\n",sum_norm);
    if (sum_norm < 1e-10)
    {
        #ifdef DEBUG
        printf("warning: sum_norm = %g something wrong!\n",sum_norm);
        printf("node_i %d:",*node_i);
        for(int j = 0; j < degree[*node_i]; j++) {
            printf("%d ", C[*node_i][j]  );
        }   printf("\n");
        #endif

        return 0;
    }

    double nr_compare = 0;
    for(int m_idx = 0; m_idx < m; m_idx++)
    {
        w = Random();

        counter = 0;
        sum = tags[counter];
        nr_compare =  w*sum_norm;
        while(sum < nr_compare){
            counter++;
            sum += tags[counter];
        }
        //printf("sum: %lf counter: %d sum_norm: %lf\n", sum, counter, sum_norm);

        nodes[m_idx] = counter;
        nodes_freq[m_idx] = weff[nodes[m_idx]];
        
        sum_norm = sum_norm-tags[counter];
        tags[counter] = 0; // we remove this node from the list
    }

    for(int i = 0; i < m; i++)
    {
        if (nodes[i]<0 || nodes[i]>NODE_NR) {
            printf("error: something went wrong with BA node generator\n");
            exit(12);
        }
        
    }

    return 1;
}

double diff_weff_weight(double alpha, double wi, double wj, double ri, double rj) {
    return exp(alpha*fabs(wi-wj));
}
