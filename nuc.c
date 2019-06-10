#include "graph.h"

void oscillator_on()
{
      ///////////////////////ERDOS-RENYI////////////////////////////////////////
    #ifdef ERDOS_RENYI
    initERmodel(ERDOS_RENYI_prob);         //erdos-renyi random model
    printf("Finished initializing ER network.\n");
    #endif
    #ifdef SCALE_FREE
    initScaleFree();    // following a power law degree distr.,possible hubs
    printf("Finished initializing SF network.\n");
    #endif
    #ifdef READ_NETWORK_FROM_FILE
    init_C(&C, NODE_NR, K_MAX);
    read_edgelist_file_py (FILENAME_FROM_WHICH_TO_READ);
    #endif
    #ifdef BARABASI_ALBERT
    init_BA(BA_PARAM_M, NODE_NR);
    printf("Finished initializing BA network.\n");
    #endif
    double deg_med, deg_var;
    deg_med = deg_var = 0;
    double *deg_aux;
    deg_aux = (double*)calloc(NODE_NR, sizeof *deg_aux);
    if (deg_aux != NULL)
    {
        for(int i = 0; i < NODE_NR; i++) {
            deg_aux[i] = (double) (degree[i]);
        }
        med_var(deg_aux, NODE_NR, &deg_med, &deg_var);
        printf("Statistical properties. avg: %g var: %g\n", deg_med, deg_var);
        printf("K_C = %g\n", 1/sqrt(2*M_PI) * deg_med/deg_var);
        free(deg_aux);
    }
    else
    {
        free(deg_aux);
    }
    




    write_C_to_file(); // to plot the graph

    //for(int i = 0; i < NODE_NR; i++) {
    //    printf("%d\t", degree[i]);
    //}   printf("\n");

    #ifdef DEGREE_HISTOGRAM
    calcHist("hist.txt");
    #endif

    ///////////////////////OSCILLATOR/////////////////////////////////////////

    printf("Starting oscillator test.\n");


    double r_med, r_var;
    r_med = r_var = 0;

    // INITIAL CONDITIONS:
    // gives random initial values for all the phases, theta_i
    initThetas();
    // gives values distributed normally for the nat. frequancies, mean 0 var 1
    initOmegas();


    double h               = DELTA_T;   // time increment
    int nr_measurements    = MAX_STEPS; // number of time increments we measure
    int    blind           = IN_BETWEEN;  // number of times we dont measure
    #ifdef TERMALIZATION
    int termalization = TERMALIZATION;
    #endif

    double *r_coh;
    r_coh = (double*)calloc(nr_measurements, sizeof(*r_coh));
    #ifdef PRINT_EVOLUTION_OF_R
    char filename_coh[128] = ".";
    snprintf(filename_coh, sizeof(char) * 128, "coh_%s_%g.txt",
                                                "term", SIGMA_VAL);
    FILE *theta_file = fopen(filename_coh,"w");
    if (theta_file == NULL) {
        printf("%s\n", "Could not open coh.txt");
        exit(2);
    }
    #endif

    char filename2[128];
    // this is for a varying filename
    #ifdef ERDOS_RENYI
    snprintf(filename2, sizeof(char) * 128, "sync_sigmaVSr_ER_N=%d_p=%g.txt", 
                                                NODE_NR, ERDOS_RENYI_prob);
    #endif
    #ifdef SCALE_FREE
    snprintf(filename2, sizeof(char) * 128, "sync_sigmaVSr_SF_N=%d_gamma=%g.txt", 
                                                NODE_NR, GAMMA);
    #endif
    #ifdef READ_NETWORK_FROM_FILE
    snprintf(filename2, sizeof(char) * 128, "ES_sigmaVSr_file_N=%d_%s_%g.txt", 
                                            NODE_NR, NET_TYPE, NET_CHARACT);
    #endif
    #ifdef BARABASI_ALBERT
    snprintf(filename2, sizeof(char) * 128, "sync_sigmaVSr_file_N=%d_m=%d.txt", 
                                            NODE_NR, BA_PARAM_M);
    #endif


    FILE *f_out2 = fopen(filename2,"w");
    if (f_out2 == NULL) {
        printf("Could not open sigmaVSr.txt");
        exit(3);
    }
    fprintf(f_out2, "# Node_nr = %d\n", NODE_NR);
    // Now we print the first with statistical data for the network.
    #ifdef ERDOS_RENYI
        fprintf(f_out2,"# SCALE FREE p=%g K_MIN=%d K_MAX=%d <k>=%g <k^2>=%g \n", 
                    ERDOS_RENYI_prob, K_MIN, NODE_NR, deg_med, deg_var);
    #endif
    #ifdef SCALE_FREE
    fprintf(f_out2,"# SCALE FREE GAMMA=%g K_MIN=%d K_MAX=%d <k>=%g <k^2>=%g \n", 
                    GAMMA, K_MIN, NODE_NR, deg_med, deg_var);
    #endif
    fprintf(f_out2, "# K <r> sigma_r(sig not of average! divide by sqrt(n))\n");
    //printf("passed second ifdefs\n");

    // For Runge kutta I will define these auxiliary variables:

    
    double timp = 0;
    // If we wait to termalize we run through many sigmas, else we only measure
    // a single sigma and plot coh.txt to see how much we should wait until the
    // phase coherence reaches a stable value, if at all.
    double sigma = 0;
    #ifndef TERMALIZATION
    sigma = SIGMA_VAL; // no for loop, a SINGLE SIGMA VALUE
    int node_nr_aux_term = 3;
    #endif
    #ifdef TERMALIZATION
    
    double sigma_min = SIGMA_MIN;
    double sigma_max = SIGMA_MAX;
    double sigma_inc = (sigma_max-sigma_min)/NR_SIGMA;
    #endif
    #ifdef TERMALIZATION
    for ( sigma = sigma_min; sigma < sigma_max; sigma += sigma_inc)
    #endif
    {
        //for(int i = 0; i < NODE_NR; i++) {
        //    GLOB_theta[i] = M_PI*(-1 + Random()*2);
        //}
        

        #ifdef TERMALIZATION // we wait for r to stabilize
        for (int i = 0; i < termalization; i++) {
            for (int t_aux2 = 0; t_aux2 < blind; t_aux2++) {
                update_RK(timp, sigma, h);
                //update_EULER(sigma, h);
            }
        }
        #endif

        for (int t_idx = 0; t_idx < nr_measurements; t_idx++)
        {
            r_coh[t_idx] = phase_coherence();
            //r_coh[t_idx] = (GLOB_theta[node_nr_aux_term]);
            timp = t_idx*h*(1+blind);

            #ifdef PRINT_EVOLUTION_OF_R
            fprintf(theta_file, "%lf %lf %lf\n",
     timp + (sigma-sigma_min)/sigma_inc * nr_measurements*h*(1+blind),
            r_coh[t_idx], 0.0);//exp(-0.01*timp));
            #endif

            for (int t_aux = 0; t_aux < blind; t_aux++) {
                // "blind" because we update without measuring
                update_RK(timp, sigma, h);
                //update_EULER(sigma, h); 
            }

            update_RK(timp, sigma, h);
            //update_EULER(sigma, h);
        }

        
        med_var(r_coh, nr_measurements, &r_med, &r_var);
        fprintf(f_out2, "%g %g %g \n", 
                sigma, r_med, sqrt(r_var));
        printf("One loop finished! \t sigma=%g \t <r>=%g \t sigma_<r>=%g \n",
                    sigma, r_med, sqrt(r_var));
    }

    fclose(f_out2);     f_out2 = NULL;
    #ifdef PRINT_EVOLUTION_OF_R
    fclose(theta_file); theta_file = NULL ;
    #endif

}

void percolation_on()
{
       // saveAdj();
    int t_number = 100;
    double t_min = 0.0;
    double t_max = 1.0;
    double t_inc = (t_max-t_min)/(t_number); // sigma increments
    double t = t_min;
    double fractional_size_noavg[t_number][AVG_NUMBER];
    for(int i = 0; i < t_number; i++) {
        for(int j = 0; j < AVG_NUMBER; j++) {
            fractional_size_noavg[i][j] = 0;
        }
    }
    double sigma = SIGMA_VAL;


    double *fractional_size;
    fractional_size = (double*)calloc(t_number, sizeof *fractional_size);
    double *fractional_size_sigma ;
    fractional_size_sigma = (double*)calloc(t_number, sizeof *fractional_size_sigma);
    //print_vec(&fractional_size, (int)t_number);
    double *edge_fraction ;
    edge_fraction = (double*)calloc(t_number, sizeof *edge_fraction);
    int aux_int = 0;
    int idx = 0;

    for(int avg_idx = 0; avg_idx < AVG_NUMBER; avg_idx++)
    {
        init_C(&C, NODE_NR, K_MAX);
        idx = 0;
        for ( t=t_min; t<t_max; t += t_inc)
        {   //TODO: Make it a running average
            //init_C(&C, NODE_NR, K_MAX);
            initEXPL_product_rule(t, sigma);

            //aux_int = int_max_vector(&GLOB_component_size, NODE_NR);
            aux_int = (int) (GLOB_max_component_size);
            //printf("max comp: %d \n", aux_int);
            fractional_size_noavg[idx][avg_idx] = ((double)aux_int)/NODE_NR;
            //printf("t=%g, fractional_size=%g\n", t, 
            //            fractional_size_noavg[idx][avg_idx]);
            //printf("Glob_nr_edges= %d\n", GLOB_nr_edges);
            //printf("component name unique elements: %d\n",
            //                unique_elements(GLOB_component_name, NODE_NR));
            printf("GLOB_nr_edges=%d \t max_comp_size=%g unique_elem:%g\n",
                                            GLOB_nr_edges,
                                            GLOB_max_component_size,
                                            (GLOB_unique_elements_in_network));
            edge_fraction[idx] = t;
            idx++;    


        }
        if (AVG_NUMBER > 2)
        {
            printf("avgnr= %d\n", avg_idx);
        }
    }

    // calculate averages
    double *aux_array = NULL;
    aux_array = (double*)calloc(t_number, sizeof *aux_array);
    idx = 0;
    for ( t=t_min; t<t_max; t += t_inc) {
        for(int j = 0; j < AVG_NUMBER; j++) {
            aux_array[j] = fractional_size_noavg[idx][j];
        }
        med_var(aux_array, AVG_NUMBER, &(fractional_size[idx]),
                                        &(fractional_size_sigma[idx]));
        idx++;
    }
    free(aux_array); //aux_array = NULL;

    write_C_to_file(); // to plot the graph

    // debug prints:
    //for(int i = 0; i < (int)t_number; i++) {
    //    printf("%g\t", fractional_size[i]);
    //}   printf("\n");
    
    //printf("component size: ");
    //for(int i = 0; i < NODE_NR; i++) {
    //    printf("%d\t", GLOB_component_size[i]);
    //}   printf("\n");
    //printf("component size unique elemnts: %d\n",
    //                    unique_elements(GLOB_component_size, NODE_NR));
    //printf("component name: ");
    //for(int i = 0; i < NODE_NR; i++) {
    //    printf("%d\t", GLOB_component_name[i]);
    //}   printf("\n");
    //printf("component name unique elemnts: %d\n",
    //                    unique_elements(GLOB_component_name, NODE_NR));
    //print_C();

    // We write to file percolation data.

    char filename1[128] = "frac_size_vs_t.txt";
    snprintf(filename1, sizeof(char) * 128, "N=%d_fracsize_vs_t.txt", NODE_NR);
    // we will store sigmas and phase coherence r
    FILE *f_out1 = fopen(filename1,"w");
    if (f_out1 == NULL) {
        printf("Could not open file.txt");
        exit(4);
    }
    for(int i = 0; i < t_number; i++)
    {
        fprintf(f_out1, "%g\t%g\t%g\n", edge_fraction[i], fractional_size[i],
                                sqrt(fractional_size_sigma[i]/AVG_NUMBER));
    }
    fclose(f_out1);                
      
}

void epes_on ()
{
    //int t_blind  = 0;
    double t_min = 0.0;
    double t_max = 3.0;//1.0*(NODE_NR-1)/2;
    int t_number = t_max*NODE_NR;//10*(NODE_NR-1)/2.0;
    double t_inc = (t_max-t_min)/(t_number); // sigma increments
    double t     = t_min;
    double fractional_size_noavg[t_number][AVG_NUMBER];
    for(int i = 0; i < t_number; i++) {
        for(int j = 0; j < AVG_NUMBER; j++) {
            fractional_size_noavg[i][j] = 0;
        }
    }

    double *fractional_size;
    fractional_size = (double*)calloc(t_number, sizeof *fractional_size);
    double *fractional_size_sigma ;
    fractional_size_sigma = (double*)calloc(t_number, sizeof *fractional_size_sigma);
    //print_vec(&fractional_size, (int)t_number);
    double *edge_fraction ;
    edge_fraction = (double*)calloc(t_number, sizeof *edge_fraction);
    int aux_int = 0;
    int idx = 0;
    double clustering = 0;

    // For coherece
    double r_med, r_var;
    r_med = r_var = 0;

    
    double timp = 0;
    // If we wait to termalize we run through many sigmas, else we only measure
    // a single sigma and plot coh.txt to see how much we should wait until the
    // phase coherence reaches a stable value, if at all.
    double sigma = SIGMA_VAL;

    double h               = DELTA_T;   // time increment
    int nr_measurements    = MAX_STEPS; // number of time increments we measure
    int    blind           = IN_BETWEEN;  // number of times we dont measure
    #ifdef TERMALIZATION
    int termalization = TERMALIZATION;
    #endif

    double *r_coh;
    r_coh = (double*)calloc(nr_measurements, sizeof(*r_coh)); // r from phase coherence
    #ifndef TERMALIZATION
    char filename_coh[128] = ".";
    snprintf(filename_coh, sizeof(char) * 128, "coh_%s_%g.txt",
                                                "term", SIGMA_VAL);
    FILE *theta_file = fopen(filename_coh,"a");
    if (theta_file == NULL) {
        printf("%s\n", "Could not open coh.txt");
        exit(2);
    }
    #endif
    
    char filename2[128] = "aa";
    snprintf(filename2, sizeof(char) * 128, "EPES_N=%d_%s_sig=%g.txt", 
                                        NODE_NR, EPES_CHARACT, SIGMA_VAL);

    FILE *f_out2 = fopen(filename2,"w");
    if (f_out2 == NULL) {
        printf("Could not open %s.txt", filename2);
        exit(15);
    }
    fprintf(f_out2, "# Node_nr = %d\n", NODE_NR);
    // Now we print the first with statistical data for the network.
    fprintf(f_out2,"# K <r> sigma_r(sig not of average! divide by sqrt(n))\n");
    //printf("passed second ifdefs\n");

    // INITIAL CONDITIONS:
    // gives random initial values for all the phases, theta_i
    initThetas();
    // gives values distributed normally for the nat. frequancies, mean 0 var 1
    initOmegas();

    int avg_idx = 0;
    //for(int avg_idx = 0; avg_idx < AVG_NUMBER; avg_idx++)
    {
        init_C(&C, NODE_NR, K_MAX);
        idx = 0;
        for ( t=t_min; t<=t_max; t += t_inc)
        {   
            //init_C(&C, NODE_NR, K_MAX);
            initEXPL_product_rule(t, sigma);  

            //for(int i = 0; i < NODE_NR; i++) {
            //    GLOB_theta[i] = M_PI*(-1 + Random()*2);
            //}
            //for (int i = 0; i < NODE_NR; i++) {
            //    GLOB_omega_nat[i] = (double) degree[i];
            //            //0.5*(-1 + 2*Random());//sampleNormal();
            //}

            #ifdef TERMALIZATION // we wait for r to stabilize
            for (int i = 0; i < termalization; i++)
                for (int t_aux = 0; t_aux < blind; t_aux++)
                    update_RK(timp, sigma, h);
            #endif
            //printf("passed termalization\n");

            for (int t_idx = 0; t_idx < nr_measurements; t_idx++)
            {
                r_coh[t_idx] = phase_coherence();
                //r_coh[t_idx] = (GLOB_theta[node_nr_aux_term]);
                timp = t_idx*h*(1+blind);

                #ifndef TERMALIZATION
                fprintf(theta_file, "%lf %lf %lf\n",
         timp + (sigma-sigma_min)/sigma_inc * nr_measurements*h*(1+blind),
                r_coh[t_idx], //exp(-0.01*timp));
                            sin(timp));
                #endif

                for (int t_aux = 0; t_aux < blind; t_aux++) {
                    // "blind" because we update without measuring
                    update_RK(timp, sigma, h); 
                }

                update_RK(timp, sigma, h);
            }


            //aux_int = int_max_vector(&GLOB_component_size, NODE_NR);
            aux_int = (int) (GLOB_max_component_size);
            //printf("max comp: %d \n", aux_int);
            fractional_size_noavg[idx][avg_idx] = ((double)aux_int)/NODE_NR;
        
            //clustering = Clustering();
            clustering=0;
            med_var(r_coh, nr_measurements, &r_med, &r_var);
            fprintf(f_out2, "%g %g %g %g %g\n", 
    t, fractional_size_noavg[idx][avg_idx], r_med, sqrt(r_var), clustering);
    //printf("One loop finished! \t sigma=%g \t <r>=%g \t sigma_<r>=%g \n",
            //            sigma, r_med, sqrt(r_var));

            //printf("t=%g, fractional_size=%g\n", t, 
            //            fractional_size_noavg[idx][avg_idx]);
            //printf("Glob_nr_edges= %d\n", GLOB_nr_edges);
            //printf("component name unique elements: %d\n",
            //                unique_elements(GLOB_component_name, NODE_NR));
printf("t=%d/%d \t max_comp=%g \t unique_elem:%g \t <r>=%g \t sigma_<r>=%g clust=%g\n",
                                            GLOB_nr_edges, NODE_NR,
                                            GLOB_max_component_size,
                                            (GLOB_unique_elements_in_network),
                                            r_med,
                                            sqrt(r_var),
                                            clustering);
                        


            edge_fraction[idx] = t;
            idx++;
            
        }
        if (AVG_NUMBER >= 2)
        {
            printf("avgnr= %d\n", avg_idx);
        }
    }
    

    ////////////////////////////


                        
    
    ///////////////////////////
    fclose(f_out2);     f_out2 = NULL;
    #ifndef TERMALIZATION
    fclose(theta_file); theta_file = NULL ;
    #endif

    // calculate averages
    double *aux_array = 0;
    aux_array = (double*)calloc(t_number, sizeof *aux_array);
    idx = 0;
    for ( t=t_min; t<t_max; t += t_inc) {
        for(int j = 0; j < AVG_NUMBER; j++) {
            aux_array[j] = fractional_size_noavg[idx][j];
        }
        med_var(aux_array, AVG_NUMBER,  &(fractional_size[idx]),
                                        &(fractional_size_sigma[idx]));
        idx++;
    }
    free(aux_array);

    write_C_to_file(); // to plot the graph

    // debug prints:
    //for(int i = 0; i < (int)t_number; i++) {
    //    printf("%g\t", fractional_size[i]);
    //}   printf("\n");
    
    //printf("component size: ");
    //for(int i = 0; i < NODE_NR; i++) {
    //    printf("%d\t", GLOB_component_size[i]);
    //}   printf("\n");
    //printf("component size unique elemnts: %d\n",
    //                    unique_elements(GLOB_component_size, NODE_NR));
    //printf("component name: ");
    //for(int i = 0; i < NODE_NR; i++) {
    //    printf("%d\t", GLOB_component_name[i]);
    //}   printf("\n");
    //printf("component name unique elemnts: %d\n",
    //                    unique_elements(GLOB_component_name, NODE_NR));
    //print_C();

    /*
    // We write to file percolation data.
    char filename1[128] = "frac_size_vs_t.txt";
    snprintf(filename1, sizeof(char) * 128, "N=%d_fracsize_vs_t.txt", NODE_NR);
    // we will store sigmas and phase coherence r
    FILE *f_out1 = fopen(filename1,"w");
    if (f_out1 == NULL) {
        printf("Could not open file.txt");
        exit(4);
    }
    for(int i = 0; i < t_number; i++)
    {
        fprintf(f_out1, "%g\t%g\t%g\n", edge_fraction[i], fractional_size[i],
                                sqrt(fractional_size_sigma[i]/AVG_NUMBER));
    }
    fclose(f_out1);                 f_out1                  = NULL;
    */


}

void debug()
{
    add_edge(2,3);
    add_edge(3,3);
    add_edge(2,3);
    add_edge(0,1);
    add_edge(0,3);
    add_edge(6,7);
    //add_edge(3,5);

    print_linked_list();
/*
        for(size_t i = 0; i < NODE_NR; i++)
    {
        printf("comp: %d ||||| ", i);
        Node crawl;
        crawl = GLOB_dom->suc[i];
        while(crawl != NULL) {
            printf("%d ", crawl->id);
            crawl = crawl->next;
        }
        printf("size: %d\n", (int)GLOB_dom_size[i]);
        
    }*/
    printf("1\n");
    int m=3;
    int *nodes = (int*)calloc(m, sizeof *nodes);
    generate_node_BA(m, nodes);
    printf("generated nodes:\n");
    for(int i = 0; i < m; i++) {
        printf("%d ", nodes[i]);
    }
    printf("1\n");


    // INITIAL CONDITIONS:
    // gives random initial values for all the phases, theta_i
    initThetas();
    // gives values distributed normally for the nat. frequancies, mean 0 var 1
    initOmegas();

    double wi, wf;
    wi = wf = 0;
    printf("2\n");
    int node_i = 3;
    FILE *f1 = fopen("test_theor.txt","w");
    printf("3\n");
    double www = 0;
    for(int ii = 0; ii < NODE_NR; ii++)
    {
    //wi = weff_compt(GLOB_component_name[node_i], 0.0, SIGMA_VAL);
    //wf = weff_compt(GLOB_component_name[ii], 0.0, SIGMA_VAL);
    www = diff_weff_weight(FG_ALPHA, wi, wf, 0, 0);
    for (int j=0; j<degree[node_i]; j++)
    {
        if (C[node_i][j]==ii)
        {
            www = 0;
        }
    }

    fprintf(f1, "%d %g\n", ii,  www);
    }
    fclose(f1); f1=NULL;
    printf("4\n");

    FILE *f = fopen("test_FG.txt","w");
    printf("6\n");
    int fg_nodes[FG_M];
    for (size_t i = 0; i < FG_M; i++)
    {
        fg_nodes[i]=0;
    }
    
    for(long int i = 0; i < (long int)(1e3); i++)
    {
        //generate_node_FREQUENCY_GAP(FG_ALPHA, node_i, FG_M, fg_nodes,
        //                                                0, SIGMA_VAL);
        //printf("%d\n", fg_nodes[0]);
    fprintf(f, "%d\n", fg_nodes[0]);
    }
    printf("1\n");
    fclose(f); f=NULL;

    write_C_to_file();
}

void frequency_gap_on()
{
        //int t_blind  = 0;
    double t_min = FG_T_MIN;
    double t_max = FG_T_MAX;//1.0*(NODE_NR-1)/2;
    int t_number = FG_T_NUMBER;//t_max*NODE_NR;//10*(NODE_NR-1)/2.0;
    double t_inc = (t_max-t_min)/(t_number); // sigma increments
    double t     = t_min;
    double fractional_size_noavg[t_number][AVG_NUMBER];
    for(int i = 0; i < t_number; i++) {
        for(int j = 0; j < AVG_NUMBER; j++) {
            fractional_size_noavg[i][j] = 0;
        }
    }

    double fractional_size[t_number], fractional_size_sigma[t_number],
                                                edge_fraction[t_number];
    for (int i = 0; i < t_number; i++) {
        fractional_size[i] = fractional_size_sigma[i] = edge_fraction[i] = 0;
    }
    
    
    int aux_int = 0;
    int idx = 0;
    double clustering = 0;

    // For coherece
    double r_med, r_var;
    r_med = r_var = 0;
    double r_dpsidt, r_dpsidt_var;
    r_dpsidt = r_dpsidt_var = 0;    

    
    double timp = 0;
    // If we wait to termalize we run through many sigmas, else we only measure
    // a single sigma and plot coh.txt to see how much we should wait until the
    // phase coherence reaches a stable value, if at all.
    double sigma = SIGMA_VAL;
    double fg_alpha = FG_ALPHA;

    double h               = DELTA_T;   // time increment
    int nr_measurements    = MAX_STEPS; // number of time increments we measure
    int    blind           = IN_BETWEEN;  // number of times we dont measure
    #ifdef TERMALIZATION
    int termalization = TERMALIZATION;
    #endif

    double *r_coh, *r_psi;
    r_coh = (double*)calloc(nr_measurements, sizeof(*r_coh)); // r from phase coherence
    r_psi = (double*)calloc(nr_measurements, sizeof(*r_psi));
    #ifdef PRINT_EVOLUTION_OF_R
    char filename_coh[128] = ".";
    snprintf(filename_coh, sizeof(char) * 128, "data/coh_N=%d_a=%g_s=%g_m=%d.txt",
                                                NODE_NR,
                                                FG_ALPHA,
                                                SIGMA_VAL,
                                                FG_M);
    FILE *theta_file = fopen(filename_coh,"a");
    if (theta_file == NULL) {
        printf("%s\n", "Could not open coh.txt");
        exit(2);
    }
    #endif

////////
/*
    int col_nr = WEFF_MEMORY_DYNAMIC_MATRIX_INI_SIZE;
    C_dom =  (int**)malloc(NODE_NR * sizeof *C_dom);
    if (C_dom==NULL)
    {
        printf("warning: malloc gives bad results...\n");
        exit(111);
    }
    int *ax;
    for (int k = 0; k < NODE_NR; k++) {
        //ax = (int *)malloc(col_nr * sizeof *(C_dom[k]) );
        //C_dom[k] = ax;
        C_dom[k] = malloc(col_nr * sizeof *(C_dom[k]) );
        printf("apple's address = %p\n", C_dom[k]);
        if (C_dom[k]==NULL)
        {
            printf("warning: malloc gives bad results...\n");
            exit(111);
        }
    }
    for (int i = 0; i < NODE_NR; i++) {
        for (int j = 0; j < col_nr; j++) {
            C_dom[i][j] = -1;
        }
    }

    for (int i = 0; i < NODE_NR; i++)
    {
        C_dom_sizes[i] = col_nr;
    }
*/
////
    
    char filename2[128] = "aa";
    snprintf(filename2, sizeof(char) * 128, "data/FG_N=%d_m=%d_%d_a=%g_sig=%g%s.txt", 
             NODE_NR, FG_M, FG_ACLIOPTAS_K, FG_ALPHA, SIGMA_VAL, EXTRA_LABEL);

    FILE *f_out2 = fopen(filename2,"w");
    if (f_out2 == NULL) {
        printf("Could not open %s.txt", filename2);
        exit(15);
    }
    fprintf(f_out2, "# Relevant parameters to reproduce this network:\n");
    fprintf(f_out2, "# %s%d\t %s%d\n", "node_nr=", NODE_NR, "k_max=", K_MAX);
    fprintf(f_out2, "# %s%d\t %s%g\t %s%d\t %s%d\n",
                    "initial_seed=", GLOB_initial_seed,
                    "delta_t=", DELTA_T,
                    "max_steps=", MAX_STEPS,
                    "in_between=", IN_BETWEEN);
    #ifdef FREQUENCY_GAP
    fprintf(f_out2, "# %s%g\t %s%g\t %s%g\t %s%g\t %s%g\t %s%g\t %s%g\t %s%g\t %s%g\n",
                    "fg_m=", (double)FG_M,
                    "fg_alpha=", (double)FG_ALPHA,
                    "fg_t_min=", (double)FG_T_MIN,
                    "fg_t_max=", (double)FG_T_MAX,
                    "fg_t_number=", (double)FG_T_NUMBER,
                    "fg_achlioptas_k=", (double)FG_ACLIOPTAS_K,
                    "fg_weff_lower_frequency=", (double)FG_WEFF_LOWER_FREQUENCY,
                    "fg_weff_max_steps=", (double)FG_WEFF_MAX_STEPS,
                    "fg_weff_in_between=", (double)FG_WEFF_IN_BETWEEN);
    #endif

    // in the following file we will write edges added with their weff's
    char filename3[128] = "aa";
    snprintf(filename3, sizeof(char) * 128, 
                        "data/FG_N=%d_m=%d_%d_a=%g_sig=%g_EDGELIST%s.txt", 
            NODE_NR, FG_M, FG_ACLIOPTAS_K, FG_ALPHA, SIGMA_VAL, EXTRA_LABEL);

    FILE *f_out_edgelist = fopen(filename3,"w");
    if (f_out_edgelist == NULL) {
        printf("Could not open %s.txt", filename3);
        exit(16);
    } 

    // Now we print the first with statistical data for the network.
    fprintf(f_out2,"# K Giant_component <r> sigma_r*sqrt(NODE_NR) clustering\n");
    //printf("passed second ifdefs\n");

    // INITIAL CONDITIONS:
    // gives random initial values for all the phases, theta_i
    initThetas();
    // gives values distributed normally for the nat. frequancies, mean 0 var 1
    initOmegas();
    double aux_ang22=0;

    int avg_idx = 0;
    //for(int avg_idx = 0; avg_idx < AVG_NUMBER; avg_idx++)
    {
        init_C(&C, NODE_NR, K_MAX);
        idx = 0;
        #ifdef TERMALIZATION // we wait for r to stabilize
        for (int i = 0; i < termalization; i++)
            for (int t_aux = 0; t_aux < blind; t_aux++)
                update_RK(timp, sigma, h);
        #endif
        for ( t=t_min; t<=t_max; t += t_inc)
        {   
            increase_edges_FREQ_GAP(t, FG_M, fg_alpha, 0.0, sigma,
                                         f_out_edgelist, &r_med, &r_var);
            
            for(int i = 0; i < NODE_NR; i++) {
                //GLOB_theta[i] = M_PI*(-1 + Random()*2);
                aux_ang22 =  GLOB_theta[i];
                GLOB_theta[i] = atan2(sin(aux_ang22),cos(aux_ang22)); // reset to [0,2pi]
            }
            
            //for (int i = 0; i < NODE_NR; i++) {
            //    GLOB_omega_nat[i] = 2*(double) degree[i];
            //            //0.5*(-1 + 2*Random());//sampleNormal();
            //}
/* !!!!
            #ifdef TERMALIZATION // we wait for r to stabilize
            for (int i = 0; i < termalization; i++)
                for (int t_aux = 0; t_aux < blind; t_aux++)
                    update_RK(timp, sigma, h);
            #endif
            //printf("passed termalization\n");

            //r_coh[0] = phase_coherence();
            //r_psi[0] = psi_coherence();
            for (int t_idx = 0; t_idx < nr_measurements; t_idx++)
            {
                timp = t_idx*h*(1+blind);
                r_coh[t_idx] = phase_coherence();
                //r_psi[t_idx] = psi_coherence();
                //if (r_med > 0.5)
                {
                //    r_psi[t_idx] = (psi_coherence()-r_psi[t_idx-1])/(h*(1+blind));
                }
                //r_coh[t_idx] = (GLOB_theta[node_nr_aux_term]);

                #ifdef PRINT_EVOLUTION_OF_R
                fprintf(theta_file, "%lf %lf\n",
         timp + (t-t_min)/t_inc * nr_measurements*h*(1+blind),
                r_coh[t_idx]);
                #endif

                for (int t_aux = 0; t_aux < blind; t_aux++) {
                    // "blind" because we update without measuring
                    update_RK(timp, sigma, h); 
                }

                update_RK(timp, sigma, h);
            }
            //r_psi[0]=r_psi[1];
            //print_vec(r_psi, nr_measurements);
!!!
*/ 
            //aux_int = int_max_vector(&GLOB_component_size, NODE_NR);
            aux_int = (int) (GLOB_max_component_size);
            //printf("max comp: %d \n", aux_int);
            fractional_size_noavg[idx][avg_idx] = ((double)aux_int)/NODE_NR;
        
            //clustering = Clustering();
            clustering=0;

// *******            med_var(r_coh, nr_measurements, &r_med, &r_var);

            //med_var(r_psi, nr_measurements, &r_dpsidt, &r_dpsidt_var);
            //r_dpsidt = (r_psi[nr_measurements-1]-r_psi[0])/timp;


            fprintf(f_out2, "%g %g %g %g %g\n", 
    t, fractional_size_noavg[idx][avg_idx], r_med, sqrt(r_var), clustering);
    //printf("One loop finished! \t sigma=%g \t <r>=%g \t sigma_<r>=%g \n",
            //            sigma, r_med, sqrt(r_var));

            //printf("t=%g, fractional_size=%g\n", t, 
            //            fractional_size_noavg[idx][avg_idx]);
            //printf("Glob_nr_edges= %d\n", GLOB_nr_edges);
            //printf("component name unique elements: %d\n",
            //                unique_elements(GLOB_component_name, NODE_NR));
  if (GLOB_nr_edges%100==0)
  {  printf("t=%d/%d=%.3lf \t max_comp=%g\tunique_elem:%g\t r=%g\t N^0.5sig_r=%g\t psi=%g sig_psi=%g\n",
                                            GLOB_nr_edges, NODE_NR, t, 
                                            GLOB_max_component_size,
                                            (GLOB_unique_elements_in_network),
                                            r_med,
                                            sqrt(r_var),
                                            r_dpsidt,
                                            sqrt(r_dpsidt_var));
  }                     


            edge_fraction[idx] = t;
            idx++;
            
        }
        fclose(f_out_edgelist);

#ifdef HISTERESIS
        float aux11,aux22,aux33,aux44;//wase variabls
        int aux55;
        int line_nr_rev = count_lines_in_file(filename3);
        if (line_nr_rev != GLOB_nr_edges) {
            printf("warning: read edges not total nr of edges %d %d\n",
                                 line_nr_rev, GLOB_nr_edges);
        }
        
        int node_rev[line_nr_rev][2];
        FILE *f_rev = fopen(filename3,"r");
        if (f_rev != NULL) {
            int i=0;
            while ( fscanf(f_rev, "%d\t%d\t%d\t%g\t%g\t%g\t%g\n",
                                 &aux55,
                                 &node_rev[i][0], &node_rev[i][1],
                                 &aux11, &aux22,
                                 &aux33, &aux44) != EOF   ) 
            {
                i++;
            }
            fclose(f_rev);
        } else {
            printf("warning: could not open y0_ini filename!\n");
            exit(1);
        }


        int flagEOF = 0;
        int aux_nd1, aux_nd2;
        for ( t=t_max; t>=HISTERESIS_MIN; t -= t_inc)
        {
            while (  ( GLOB_nr_edges > (int)(t*((double)NODE_NR)) )   )  
            {
                aux_nd1 = GLOB_nr_edges-1;
                aux_nd2 = GLOB_nr_edges-1;
                //printf("removed edge: %d %d\n",
                //             node_rev[aux_nd1][0], node_rev[aux_nd2][1]);
                remove_edge(node_rev[aux_nd1][0],
                            node_rev[aux_nd2][1]);  
            }
            #ifdef TERMALIZATION // we wait for r to stabilize
            for (int i = 0; i < termalization; i++)
                for (int t_aux = 0; t_aux < blind; t_aux++)
                    update_RK(timp, sigma, h);
            #endif
            for (int t_idx = 0; t_idx < nr_measurements; t_idx++)
            {
                timp = t_idx*h*(1+blind);
                r_coh[t_idx] = phase_coherence();
                #ifdef PRINT_EVOLUTION_OF_R
                fprintf(theta_file, "%lf %lf\n",
         timp + (t-t_min)/t_inc * nr_measurements*h*(1+blind),
                r_coh[t_idx]);
                #endif
                for (int t_aux = 0; t_aux < blind; t_aux++) {
                    update_RK(timp, sigma, h); 
                }
                update_RK(timp, sigma, h);
            }
            med_var(r_coh, nr_measurements, &r_med, &r_var);
            fprintf(f_out2, "%g %g %g %g %g\n", t, 0., r_med, sqrt(r_var), clustering);
            printf("t=%d/%d=%.3lf \t max_comp=%g\tunique_elem:%g\t r=%g\t N^0.5sig_r=%g\t psi=%g sig_psi=%g\n",
                                            GLOB_nr_edges, NODE_NR, t, 
                                            0.,
                                            0.,
                                            r_med,
                                            sqrt(r_var),
                                            0.,
                                            0.);
        }
#endif

        if (AVG_NUMBER > 2)
        {
            printf("avgnr= %d\n", avg_idx);
        }
    }
    

    ////////////////////////////


                        
    
    ///////////////////////////
    //fclose(f_out2);
    //fclose(f_out_edgelist); 
    #ifndef TERMALIZATION
    fclose(theta_file);
    #endif

    // calculate averages
    double aux_array[t_number];
    for (int i = 0; i < t_number; i++) {
        aux_array[i] = 0;
    }
    idx = 0;
    for ( t=t_min; t<t_max; t += t_inc) {
        for(int j = 0; j < AVG_NUMBER; j++) {
            aux_array[j] = fractional_size_noavg[idx][j];
        }
        med_var(aux_array, AVG_NUMBER,  &(fractional_size[idx]),
                                        &(fractional_size_sigma[idx]));
        idx++;
    }

    write_C_to_file(); // to plot the graph

}

