/*
***Model of complex network synchronization and percolations
***Emanuel Cristian Boghiu, MSc Student 
***embogh@gmail.com OR 697566@unizar.es
*/

/*READ graph.h FOR CONTROLS*/
#include "graph.h"

int main(int args_number, char* args[])
{
    //ini_ran(time(NULL));  // Seed for the random generator.
    ini_ran(91028401);     // we want predictable results

    // We reserve memory for the global variables.
    init_C_memory(&C, NODE_NR, K_MAX);

    //////////////////////PERCOLATION//////////////////////////////////////////

    // saveAdj();
    
    int t_number = 10;
    double t_min = 0.0;
    double t_max = 0.99;
    double t_inc = (t_max-t_min)/(t_number); // sigma increments
    double t = t_min;
    double fractional_size_noavg[t_number][AVG_NUMBER];
    for(size_t i = 0; i < t_number; i++)
    {
        for(size_t j = 0; j < AVG_NUMBER; j++)
        {
            fractional_size_noavg[i][j] = 0;
        }
        
    }
/*
    double *fractional_size;
    fractional_size = calloc(t_number, sizeof *fractional_size);
    double *fractional_size_sigma;
    fractional_size_sigma = calloc(t_number, sizeof *fractional_size_sigma);
    //print_vec(&fractional_size, (int)t_number);
    double *edge_fraction;
    edge_fraction = calloc(t_number, sizeof *edge_fraction);
    int aux_int = 0;
    int idx = 0;

    for(size_t avg_idx = 0; avg_idx < AVG_NUMBER; avg_idx++)
    {
        init_C(&C, NODE_NR, K_MAX);
        idx = 0;
        for ( t=t_min; t<t_max; t += t_inc)
        {
            //init_C(&C, NODE_NR, K_MAX);
            initEXPL_product_rule(t);

            aux_int = int_max_vector(&component_size, NODE_NR);
            //printf("max comp: %d \n", aux_int);
            fractional_size_noavg[idx][avg_idx] = ((double)aux_int)/NODE_NR;
            printf("t=%g, fractional_size=%g\n", t, 
                        fractional_size_noavg[idx][avg_idx]);
            printf("Glob_nr_edges= %d\n", GLOB_nr_edges);
            //printf("component name unique elements: %d\n",
            //                unique_elements(component_name, NODE_NR));
            edge_fraction[idx] = t;
            idx++;    
        }
    }

    // calculate averages
    double *aux_array;
    aux_array = calloc(t_number, sizeof *aux_array);
    idx = 0;
    for ( t=t_min; t<t_max; t += t_inc) {
        for(size_t j = 0; j < AVG_NUMBER; j++) {
            aux_array[j] = fractional_size_noavg[idx][j];
        }
        med_var(aux_array, AVG_NUMBER, &(fractional_size[idx]),
                                        &(fractional_size_sigma[idx]));
        idx++;
    }

    write_C_to_file(); // to plot the graph

    // debug prints:
    //for(int i = 0; i < (int)t_number; i++) {
    //    printf("%g\t", fractional_size[i]);
    //}   printf("\n");
    
    //printf("component size: ");
    //for(int i = 0; i < NODE_NR; i++) {
    //    printf("%d\t", component_size[i]);
    //}   printf("\n");
    //printf("component size unique elemnts: %d\n",
    //                    unique_elements(component_size, NODE_NR));
    //printf("component name: ");
    //for(int i = 0; i < NODE_NR; i++) {
    //    printf("%d\t", component_name[i]);
    //}   printf("\n");
    //printf("component name unique elemnts: %d\n",
    //                    unique_elements(component_name, NODE_NR));
    //print_C();

    // We write to file percolation data.

    char filename1[128] = "frac_size_vs_t.txt";
    snprintf(filename1, sizeof(char) * 128, "N=%d_fracsize_vs_t.txt", NODE_NR);
    // we will store sigmas and phase coherence r
    FILE *f_out1 = fopen(filename1,"w");
    if (f_out1 == NULL) {
        printf("Could not open file.txt");
        exit(1);
    }
    for(size_t i = 0; i < (int)t_number; i++)
    {
        fprintf(f_out1, "%g %g %g\n", edge_fraction[i], fractional_size[i],
                                sqrt(fractional_size_sigma[i]/AVG_NUMBER);
    }
    fclose(f_out1);

*/
    ///////////////////////ERDOS-RENYI////////////////////////////////////////
    #ifdef ERDOS_RENYI
    initERmodel(ERDOS_RENYI_prob);         //erdos-renyi random model
    printf("Finished initializing ER network.\n");
    #endif
    #ifdef SCALE_FREE
    initScaleFree();    // following a power law degree distr., possible hubs
    printf("Finished initializing SF network.\n");
    #endif
    double deg_med, deg_var;
    deg_med = deg_var = 0;
    double *deg_aux;
    deg_aux = calloc(NODE_NR, sizeof *deg_aux);
    for(size_t i = 0; i < NODE_NR; i++) {
        deg_aux[i] = (double) (degree[i]);
    }
    med_var(deg_aux, NODE_NR, &deg_med, &deg_var);
    printf("Statistical properties. avg: %g var: %g\n", deg_med, deg_var);
    printf("K_C = %g\n", 1/sqrt(2*M_PI) * deg_med/deg_var);

    write_C_to_file(); // to plot the graph

    //for(int i = 0; i < NODE_NR; i++) {
    //    printf("%d\t", degree[i]);
    //}   printf("\n");

    #ifdef DEGREE_HISTOGRAM
    calcHist("hist.txt");
    #endif

    ///////////////////////OSCILLATOR/////////////////////////////////////////
/*
    printf("Starting oscillator test.\n");

    double h, *r_coh, r_med, r_var, sigma;
    int blind;

    // INITIAL CONDITIONS:
    // gives random initial values for all the phases, theta_i
    initThetas();
    // gives values distributed normally for the nat. frequancies, mean 0 var 1
    initOmegas();

    h       = DELTA_T;          // time increment
    t_max   = MAX_STEPS;        // number of time increments we measure
    blind   = IN_BETWEEN;       // number of times we dont measure
    #ifdef TERMALIZATION
    int termalization = TERMALIZATION;
    #endif

    r_coh = (double*)calloc(t_max,sizeof(*r_coh)); // r from phase coherence
    #ifndef TERMALIZATION
    FILE *theta_file = fopen("coh.txt","w");
    if (theta_file == NULL) {
        printf("%s\n", "Could not open coh.txt");
        exit(1);
    }
    #endif

    // this is for a varying filename
    #ifdef ERDOS_RENYI
    char filename2[128] = "SF_sigmaVSr.txt";
    snprintf(filename2, sizeof(char) * 128, "SF_sigmaVSr_ER_p=%g.txt", 
                                                ERDOS_RENYI_prob);
    #endif
    #ifdef SCALE_FREE
    char filename2[128] = "SF_sigmaVSr.txt";
    snprintf(filename2, sizeof(char) * 128, "SF_sigmaVSr_SF_gamma=%g.txt", 
                                                GAMMA);
    #endif
    // we will store sigmas and phase coherence r
    FILE *f_out2 = fopen(filename2,"w");
    if (f_out2 == NULL) {
        printf("Could not open sigmaVSr.txt");
        exit(1);
    }

    // If we wait to termalize we run through many sigmas, else we only measure
    // a single sigma and plot coh.txt to see how much we should wait until the
    // phase coherence reaches a stable value, if at all.
    #ifndef TERMALIZATION
    sigma = SIGMA_VAL; // no for loop, a SINGLE SIGMA VALUE
    #endif
    #ifdef TERMALIZATION
    double sigma_min = 0.00;
    double sigma_max = 1.00;
    double sigma_inc = (sigma_max-sigma_min)/10; // sigma increments
    #endif

    fprintf(f_out2, "# Node_nr = %d\n", NODE_NR);
    // Now we print the first with statistical data for the network.
    #ifdef ERDOS_RENYI
    fprintf(f_out2, "# ERDOS RENYI PROB=%g\n",ERDOS_RENYI_prob);
    #endif
    #ifdef SCALE_FREE
    fprintf(f_out2, "# SCALE FREE GAMMA=%g K_MIN=%d K_MAX=%d <k>=%g <k^2>=%d\n", 
                    GAMMA, K_MIN, NODE_NR, deg_med, deg_var);
    #endif
    fprintf(f_out2, "# K r sigma_r\n");

    for ( sigma = sigma_min; sigma < sigma_max; sigma += sigma_inc)
    {
        #ifdef TERMALIZATION // we wait for r to stabilize
        for (int i = 0; i < termalization; i++)
            for (int t_aux = 0; t_aux < blind; t_aux++)
                update(sigma, h);
        #endif

        for (int t = 0; t < t_max; t++)
        {
            for (int t_aux = 0; t_aux < blind; t_aux++) {
                update(sigma, h); // "blind" because we update without measuring
            }
            update(sigma, h);

            r_coh[t] = phase_coherence();

            #ifndef TERMALIZATION
            fprintf(theta_file, "%d\t%g\n", t, r_coh[t]);
            #endif
        }
        // Now we measure the system.
        med_var(r_coh, t_max, &r_med, &r_var);
        fprintf(f_out2, "%g\t%g\t%g\t%s%d\n", 
                        NODE_NR*sigma, r_med, sqrt(r_var/t_max));
            printf("%s%g\n", "One loop finished!\tsigma=", sigma);
    }
    fclose(f_out2);
    free(r_coh);

*/
    printf("\nYou've reached the end without dying! ;-)\n");
    return 0;
}
