#include "graph.h"

void calcHist(char* filename)
{
    double min, max, d;
    int N_Intervalos = 20;//NODE_NR-1;
    int N_data = NODE_NR;
    double *Hist;
    Hist = calloc(N_Intervalos, sizeof *Hist);
    Histogram(GLOB_omega_nat, Hist, N_data, N_Intervalos, &d, &min, &max);
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
        exit(1);
    }
    fclose(file_hist); file_hist = NULL;
    free(Hist); Hist = NULL;
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

void Histogram (double *data, double *Hist, int N_data, int N_Intervalos,
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

    for (int i = 0; i < N_Intervalos; ++i)
        Hist[i] = 0;

    for (int i = 0; i < N_data; i++)
    {
            if ( (double) data[i] < min_aux ) min_aux = data[i];
            if ( (double) data[i] > max_aux ) max_aux = data[i];
    }
    *max = max_aux;
    *min = min_aux;
    
    // Specifically for this type of networks:
    //*min = 2;
    //*max = NODE_NR;
    //N_Intervalos = NODE_NR-2-1;

    // Calculo la longitud de cada intervalo una vez que
    // sabemos el maximo y minimo.
    *d = (*max-*min)/((double)N_Intervalos);
    if ( *d == 0 ) {
        printf("%s\n", "La longitud del intervalo es nula.");
        exit(1);
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
