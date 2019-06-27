#include "graph.h"

double GLOB_theta[NODE_NR];
void initThetas()
{
    for (int i = 0; i < NODE_NR; ++i) {
        GLOB_theta[i] = -M_PI + Random()*(2*M_PI);
        //GLOB_theta[i] = sampleNormal();
    }
}

// Generates natural frequencies about 0 according to a Gaussian distribution of
// variance 1.
double GLOB_omega_nat[NODE_NR];
void initOmegas()
{
    for (int i = 0; i < NODE_NR; i++) 
    {
        //GLOB_omega_nat[i] = sampleNormal();//(double)degree[i];//0.5*(-1 + 2*Random());//
        GLOB_omega_nat[i] = 0.5*(-1 + 2*Random());
        //GLOB_omega_nat[i] = 3+0.5*(-1 + 2*Random());
        //GLOB_omega_nat[i] = (double)degree[i];
    }

}

void update_EULER(double sigma, double h)
{ //TODO: RUNGE_KUTTA
    // we calculate the velocities of the system at this given moment
    double theta_dot[NODE_NR];
    for (int i  = 0; i < NODE_NR; i++) {
        theta_dot[i] = calculateTheta_dot_i(0, GLOB_theta, NODE_NR, sigma, i);
        //printf("theta_dot[i]=%g\n", theta_dot[i]);
    }
    for (int i = 0; i < NODE_NR; i++) { // simple Euler
        // theta(t+dt)=theta(t)+theta_dot(t)*h
        GLOB_theta[i] += theta_dot[i]*h;
    }
}

void update_RK (double t, double sigma, double h)
{
    // Here we calculate h\sum b_i k_i
    double k1[NODE_NR], k2[NODE_NR], k3[NODE_NR], k4[NODE_NR],
            theta_aux[NODE_NR], theta_aux2[NODE_NR];

    for(int i = 0; i < NODE_NR; i++) {
        k1[i] = calculateTheta_dot_i(t, GLOB_theta, NODE_NR, sigma, i);
        theta_aux[i]  = GLOB_theta[i] +  h/2 * k1[i];
    }
    for(int i = 0; i < NODE_NR; i++) {
        k2[i] = calculateTheta_dot_i(t, theta_aux, NODE_NR, sigma, i);
        theta_aux2[i] = GLOB_theta[i] + h/2 * k2[i];
    }
    for(int i = 0; i < NODE_NR; i++) {
        k3[i] = calculateTheta_dot_i(t, theta_aux2, NODE_NR, sigma, i);
        theta_aux[i] = GLOB_theta[i] + h * k3[i];
    }
    for(int i = 0; i < NODE_NR; i++) {
        k4[i] = calculateTheta_dot_i(t, theta_aux, NODE_NR, sigma, i);
    }

    for(int i = 0; i < NODE_NR; i++) {
        GLOB_theta[i] += h*(k1[i]/6 + k2[i]/3 + k3[i]/3 + k4[i]/6);
    }
}

double calculateThetaAverage()
{
    double sum = 0;
    for (int i = 0; i < NODE_NR; i++) {
        sum += GLOB_theta[i];
    }
    return sum/NODE_NR;
}

// 0 < r < 1
double phase_coherence()
{
    double rx, ry;
    rx = ry = 0;

    for (int i = 0; i < NODE_NR; i++) {
        rx += cos(GLOB_theta[i]);
        ry += sin(GLOB_theta[i]);
    }
    rx/=NODE_NR;
    ry/=NODE_NR;

    return sqrt(rx*rx+ry*ry);
}

// the following is inlined 

double psi_coherence(void)
{
    double Nrx, Nry;
    Nrx = Nry = 0;
    double phase=0;
    for (int i = 0; i < NODE_NR; i++) {
        phase=GLOB_theta[i];
        Nrx += cos(phase);
        Nry += sin(phase);
    }
    if  (Nrx < 1e-6){
        //printf("warning: ry too small for division?\n");
        return M_PI/2;
    }
    return atan2(Nry,Nrx);
}


double weff_compt_instant (int id_compt, double t, double sigma)
{
    int i=0;
    
    double weff = 0;
    int comp_size = GLOB_dom_size[id_compt];
    
#ifdef WEFF_MEMORY_LINKED_LIST
    Node crawl = GLOB_dom->suc[id_compt];
    while(crawl != NULL) {
        i = crawl->id;

        weff += calculateTheta_dot_i(t, GLOB_theta, NODE_NR, sigma, i);
        
        crawl = crawl->next;
    }
#endif
#ifdef WEFF_MEMORY_DYNAMIC_MATRIX
    
    for (i = 0; i < comp_size; i++)
    {
        weff += calculateTheta_dot_i(t, GLOB_theta, NODE_NR,
                                        sigma, C_dom[id_compt][i]);
    }
#endif
    



    weff /= comp_size;

    return weff;
}

double weff_compt (int id_compt, double t, double sigma)
{
    double timp = 0;
    double weff = 0.0;
    double w_avg = 0;
    for(int i = 0; i < FG_WEFF_MAX_STEPS; i++)
    {
        weff = weff_compt_instant (id_compt, t, sigma);
        w_avg = ((double)i)/(i+1) * w_avg + weff/(i+1);

        timp = i*DELTA_T*(1+FG_WEFF_IN_BETWEEN);
        for(int j = 0; j < FG_WEFF_IN_BETWEEN; j++) {
            update_RK(t+timp, sigma, DELTA_T); 
        }
        update_RK(t+timp, sigma, DELTA_T);
    }

    return w_avg;
}

// This function advances time globally certain amount of units. Then
// at each timestep it creates a moving average and stores it in weff_dom.
void weff_compt_efficient (double *weff_dom, int weff_dom_size,
                                double t, double sigma)
{
    
    double timp = 0;
    double weff = 0.0;
    int compont_name = 0;
    //double weff_by_domain[NODE_NR];
    for(int i = 0; i < weff_dom_size; i++)
    {
        weff_dom[i] = FG_WEFF_LOWER_FREQUENCY;
    }
    //We need to do the first step and change from default
    for(int i = 0; i < NODE_NR; i++)
    {
        // IF w<FR_WEFF.../10 then it means we have a new component whos weff
        // hasn't been calculated yet.
        compont_name = GLOB_component_name[i];
        if (weff_dom[compont_name]<(FG_WEFF_LOWER_FREQUENCY/10)) 
        {
            weff_dom[compont_name] =
                        weff_compt_instant(compont_name, t, sigma);
         //   printf("inst:%g i:%d compt:%d\n",weff_dom[compont_name],
         //                                   i,compont_name);
        }
    }

    for(int i = 0; i < FG_WEFF_MAX_STEPS; i++)
    {
        for(int idxx = 0; idxx < NODE_NR; idxx++)
        {
            // IF w>FR_WEFF.../10 then it means we it is a component which is 
            // part of the network and we update its average.
            compont_name = GLOB_component_name[idxx];
            if (weff_dom[compont_name]>(FG_WEFF_LOWER_FREQUENCY/10)) 
            {
                weff = weff_compt_instant (compont_name, t, sigma);
                weff_dom[compont_name] = 
                        ((double)i)/(i+1) * weff_dom[compont_name] + weff/(i+1);
            }
        }
        

        timp = i*DELTA_T*(1+FG_WEFF_IN_BETWEEN);
        for(int j = 0; j < FG_WEFF_IN_BETWEEN; j++) {
            update_RK(t+timp, sigma, DELTA_T); 
        }
        update_RK(t+timp, sigma, DELTA_T);
    }
}

void weff_compt_DOUBLY_efficient (double *weff_dom, double *weff_by_node,
                                    double *r_dom, int len,
                                    double *r_med, double *r_var,
                                    double t, double sigma)
{
    
    double timp = 0;

    int aux_deg=0;

    
    //double weff_by_domain[NODE_NR];
    for(int i = 0; i < NODE_NR; i++)
    {
        weff_dom[i] = FG_WEFF_LOWER_FREQUENCY;
        
//        r_dom[i] = -1;
    }
    //We need to do the first step and change from default

/*
    //double weff = 0.0;
    int compont_name = 0;
    for(int i = 0; i < NODE_NR; i++)
    {
        // IF w<FR_WEFF.../10 then it means we have a new component whos weff
        // hasn't been calculated yet.
        compont_name = GLOB_component_name[i];
        if (weff_dom[compont_name]<(FG_WEFF_LOWER_FREQUENCY/10)) 
        {
            weff_dom[compont_name] =
                        weff_compt_instant(compont_name, t, sigma);
         //   printf("inst:%g i:%d compt:%d\n",weff_dom[compont_name],
         //                                   i,compont_name);
        }
    }
*/

    
    FILE *theta_file = fopen("pattern.txt","a");

    // first hte time averages
    double r_coh[FG_WEFF_MAX_STEPS];
    double aux_time_average_theta_dot[NODE_NR];
    for (int i = 0; i < NODE_NR; i++) {
        aux_time_average_theta_dot[i]=0;
    }
    for(int i = 0; i < FG_WEFF_MAX_STEPS; i++)
    {
        for (int m = 0; m < NODE_NR; m++) {
            aux_time_average_theta_dot[m] = 
                ((double)i)/(i+1) * aux_time_average_theta_dot[m]
    + 1.0/(i+1.) * calculateTheta_dot_i(t, GLOB_theta, NODE_NR, sigma, m);
        }

/*
        for (int i_dom = 0; i_dom < NODE_NR; i_dom++)
        {
            if (C[i_dom][0] != 1)
            {
r_dom[i_dom] = ((double)i)/(i+1) * r_dom[i_dom]+ 1.0/(i+1) * phase_coherence_compt(i_dom);
            }
        }
*/
        r_coh[i] = phase_coherence();


        timp = i*DELTA_T*(1+FG_WEFF_IN_BETWEEN);
        for(int j = 0; j < FG_WEFF_IN_BETWEEN; j++) {
            update_RK(t+timp, sigma, DELTA_T); 
        }
        update_RK(t+timp, sigma, DELTA_T);
        if ((float)(GLOB_nr_edges)/NODE_NR > T_INI_MEAS && (float)(GLOB_nr_edges)/NODE_NR < T_FIN_MEAS)
        {
            fprintf(theta_file, "%d %f ",  GLOB_nr_edges, DELTA_T);
            for (int ii_aux = 0; ii_aux < NODE_NR; ii_aux++)
            {
                fprintf(theta_file, "%lf ", GLOB_theta[ii_aux]);
            }
            fprintf(theta_file, "\n");
        }
    }
    med_var(r_coh, FG_WEFF_MAX_STEPS, r_med, r_var);
    fclose(theta_file);
/*
    // after this we have the matrix of time averaged effective frequencies-
    // if space is a probelm use commented out version from above
    // now we need to do the "spatial" or component averaging 
printf("tdot avg: ");
for (int i = 0; i < NODE_NR; i++)
{
    printf("%g ", aux_time_average_theta_dot[i]);
}   printf("\n");
*/
    for (int i = 0; i < NODE_NR; i++) {
        weff_by_node[i] = aux_time_average_theta_dot[i];
    }


    double aux_size = 0;
    double aux_sum  = 0;
    for (int i = 0; i < NODE_NR; i++)
    {
        if  (C_dom[i][0]!=-1)
        {
            aux_sum = 0;
            aux_size = (double) GLOB_dom_size[i];
            for (int j = 0; j < aux_size; j++)
            {
                aux_sum += aux_time_average_theta_dot[ C_dom[i][j] ];
                //printf("%d ", C_dom[i][j]); 
            }//printf("\n");
            aux_sum /= aux_size;
            weff_dom[i]=aux_sum;
        }
        
        //printf("\n");
    }

/*
printf("wefdom: ");
for (int i = 0; i < NODE_NR; i++)
{
    printf("%g ",weff_dom[i]);
}   printf("\n");
*/
}

double phase_coherence_compt (int id_compt)
{

    int i=0;
    
    double rx, ry;
    int phase=0;

    //printf("node list: ");

    int comp_size = GLOB_dom_size[id_compt];
#ifdef WEFF_MEMORY_LINKED_LIST
    rx = ry = 0;
    phase = 0;
    Node crawl = GLOB_dom->suc[id_compt];
    while(crawl != NULL) {
        i = crawl->id;
        phase = GLOB_theta[i];
        rx += cos(phase);
        ry += sin(phase);
        //printf("%d ", i);
        crawl = crawl->next;
    } //printf("\n");
#endif
#ifdef WEFF_MEMORY_DYNAMIC_MATRIX
    rx = ry = 0;
    phase=0;
    for (i = 0; i < comp_size; i++)
    {
        phase = GLOB_theta[C_dom[id_compt][i]];
        rx += cos(phase);
        ry += sin(phase);
    }
#endif

    rx/=comp_size;
    ry/=comp_size;

    return sqrt(rx*rx+ry*ry);
}

double psi_coherence_compt (int id_compt)
{
    double Nrx, Nry;
    Nrx = Nry = 0;
    int phase=0;
    int i=0;
/*
    
    Node crawl = GLOB_dom->suc[id_compt];
    while(crawl != NULL) {
        i = crawl->id;
        Nrx += cos(GLOB_theta[i]);
        Nry += sin(GLOB_theta[i]);
        
        crawl = crawl->next;
    }
*/
    int comp_size = GLOB_dom_size[id_compt];
#ifdef WEFF_MEMORY_LINKED_LIST
    i=0;
    Nrx = Nry = 0;
    phase = 0;
    Node crawl = GLOB_dom->suc[id_compt];
    while(crawl != NULL) {
        i = crawl->id;
        phase = GLOB_theta[i];
        Nrx += cos(phase);
        Nry += sin(phase);
        //printf("%d ", i);
        crawl = crawl->next;
    } //printf("\n");
#endif
#ifdef WEFF_MEMORY_DYNAMIC_MATRIX
    Nrx = Nry = 0;
    phase=0;
    for (i = 0; i < comp_size; i++)
    {
        phase = GLOB_theta[ C_dom[id_compt][i] ];
        Nrx += cos(phase);
        Nry += sin(phase);
    }
#endif

    return atan2(Nry,Nrx);
}

//the following is inlined
double calculateTheta_dot_i(double t, double *phases, int phases_len,
                                            double sigma, int i)
{
    // \dot{theta_i}=\omega_i+|sigma\sum_{j=0}^{N} a_{ij}\sin{theta_j-theta_i}
    double sum = 0;//0*t*phases_len;
    //return cos(t);
    //double omega_nat_i = GLOB_omega_nat[i] ;
    
    #ifdef EPSILON_OSCILLATOR
    double coupling_all_to_all = EPSILON_OSCILLATOR;
    for(int j = 0; j < phases_len; j++) {
        sum += coupling_all_to_all * sin(phases[j]-phases[i]);
    }
    #endif

    //double coupling = sigma;
    
    
    //double weight = 1.0;
    //
    //if (degree[i] == 0) 
    //{
    //    weight = 1.0;
    //} 
    //else
    //{
    //    weight = 1.0 *1.0/pow(degree[i],0);
    //}
    
    sum = 0;
    double phase_i=phases[i];
    int deg=degree[i];
    for(int j = 0; j < deg; j++) {
        //coupling = sigma;
        //sum += sigma * sin( phases[ C[i][j] ] - phase[i]);
        sum += sin( phases[ C[i][j] ] - phase_i);
        //printf("%d-%d\n",i,C[i][j]);
    }
    
    return GLOB_omega_nat[i] + sigma*sum;
}


int initEXPL_product_rule (double t, double sigma)
{   
    int tot_nr_edges = 0;
    tot_nr_edges =  (int) (t* ((double)NODE_NR));
    //printf("tot_nr_edges=%d\n", tot_nr_edges);
    
    int rnd1, rnd2, rnd3, rnd4;
    //double size1 = 0;
    //double  size2, size3, size4;
    //size2 = size3 = size4 = 0;
    //rnd2 = rnd3 = rnd4 = 0;
    int bool_res = 0; // 0 false, 1 true
    //int name_aux1 = 0;
    //int name_aux2, name_aux3;
    //name_aux2 = name_aux3 = 0;

    
    while (GLOB_nr_edges < tot_nr_edges) 
{
    //GLOB_unique_components = unique_elements(GLOB_component_name, NODE_NR);
    GLOB_unique_components = (int) (GLOB_unique_elements_in_network);
    if (GLOB_unique_components < 0 || GLOB_unique_components > NODE_NR) {
        printf("warning: unique components wrong!\n");
        exit(11);
    }
    
    //printf("unique components: %d\n", GLOB_unique_components);
  //  if (GLOB_unique_components == 1) {
  //      return 0;
  //  } 
    /*  
    if (GLOB_unique_components == 2) {
        rnd1 = (int)(Random()*NODE_NR);
        rnd2 = 0;
        name_aux1 = GLOB_component_name[rnd1];
        while (GLOB_component_name[rnd2] == name_aux1) {
            rnd2++;
            if (rnd2 > NODE_NR) {
                printf("warning: rnd2>NODE_NR");
                return 0;
            }
            
        }
        add_edge(rnd1,rnd2);
        GLOB_nr_edges++;
        return 1;
    }
    */
    
    /*
    if (GLOB_unique_components == 3) {
        rnd1 = (int)(Random()*NODE_NR);
        rnd2 = 0;
        name_aux1 = GLOB_component_name[rnd1];
        while (GLOB_component_name[rnd2] == name_aux1) {
            rnd2++;
        }
        name_aux2 = GLOB_component_name[rnd2];
        rnd3 = 0;
        while ( (GLOB_component_name[rnd3] == name_aux1) ||
                (GLOB_component_name[rnd3] == name_aux2)    ) 
        {
            rnd3++;
        }
        //printf("rnds33: %d, %d, %d \n", rnd1, rnd2, rnd3);
        size1 = component_size[rnd1];
        size2 = component_size[rnd2];
        size3 = component_size[rnd3];
        if (size2 >= size3) {
            bool_res = add_edge(rnd1,rnd3);
        } else {
            bool_res = add_edge(rnd1,rnd2);
        }
        if (bool_res == 1) { // this means an edge has been added
            GLOB_nr_edges++;
            //printf("GLOB_nr_edges=%d\n", GLOB_nr_edges);
        }
        return 1;
    }*/

/*
        if (GLOB_nr_edges >= NODE_NR) {
            printf("IF STATEMENT:Glob_nr_edges= %d\n", GLOB_nr_edges);
            //exit(12);
        }
*/
        rnd1 = (int)(Random()*NODE_NR);
        rnd2 = (int)(Random()*NODE_NR);
        rnd3 = (int)(Random()*NODE_NR);
        rnd4 = (int)(Random()*NODE_NR);

/*
        name_aux1 = 0;
        rnd1 = (int)(Random()*NODE_NR);
        // TODO: remove thjese conditions
        rnd2 = (int)(Random()*NODE_NR);
        name_aux1 = GLOB_component_name[rnd1];
        while (GLOB_component_name[rnd2] == name_aux1) 
        {
            rnd2 = (int)(Random()*NODE_NR);
        } // now we have two nodes from different domains
        //printf("rnds1: %d, %d, %d \n", rnd1, rnd2, rnd3);
        rnd3 = (int)(Random()*NODE_NR);
        name_aux2 = GLOB_component_name[rnd2];
        while ( (GLOB_component_name[rnd3] == name_aux1) ||
                (GLOB_component_name[rnd3] == name_aux2)    ) 
        {
            rnd3 = (int)(Random()*NODE_NR);
        }


        //printf("rnds2: %d, %d, %d \n", rnd1, rnd2, rnd3);
        rnd4 = (int)(Random()*NODE_NR);
        name_aux3 = GLOB_component_name[rnd3];
        while (GLOB_component_name[rnd4] == name_aux3) 
        {
            rnd4 = (int)(Random()*NODE_NR);
        }
        //printf("rnds3: %d, %d, %d \n", rnd1, rnd2, rnd3);
*/

#ifdef EPES_MECH_Pure_perc
        size1 = (double) GLOB_component_size[rnd1];
        size2 = (double) GLOB_component_size[rnd2];
        size3 = (double) GLOB_component_size[rnd3];
        size4 = (double) GLOB_component_size[rnd4];
        bool_res = 0;
        if (size1*size2 >= size3*size4) 
        {
            bool_res = add_edge(rnd3,rnd4);
        } 
        else 
        {
            bool_res = add_edge(rnd1,rnd2);
        }
#endif
#ifdef EPES_MECH_compare_r
        size1 = phase_coherence_compt(GLOB_component_name[rnd1]);
        size2 = phase_coherence_compt(GLOB_component_name[rnd2]);
        size3 = phase_coherence_compt(GLOB_component_name[rnd3]);
        size4 = phase_coherence_compt(GLOB_component_name[rnd4]);
        bool_res = 0;
        if (size1*size2 >= size3*size4) 
        {
            bool_res = add_edge(rnd3,rnd4);
        } 
        else 
        {
            bool_res = add_edge(rnd1,rnd2);
        }
#endif
#ifdef EPES_MECH_Scale_by_dom_size
        size1 = phase_coherence_compt(GLOB_component_name[rnd1]);
        size2 = phase_coherence_compt(GLOB_component_name[rnd2]);
        size3 = phase_coherence_compt(GLOB_component_name[rnd3]);
        size4 = phase_coherence_compt(GLOB_component_name[rnd4]);
        size1 /= GLOB_component_size[rnd1];
        size2 /= GLOB_component_size[rnd2];
        size3 /= GLOB_component_size[rnd3];
        size4 /= GLOB_component_size[rnd4];
        bool_res = 0;
        if (size1*size2 >= size3*size4) 
        {
            bool_res = add_edge(rnd3,rnd4);
        } 
        else 
        {
            bool_res = add_edge(rnd1,rnd2);
        }
#endif
#ifdef EPES_MECH_weff
        size1 = weff_compt(GLOB_component_name[rnd1], 0, sigma);
        size2 = weff_compt(GLOB_component_name[rnd2], 0, sigma);
        size3 = weff_compt(GLOB_component_name[rnd3], 0, sigma);
        size4 = weff_compt(GLOB_component_name[rnd4], 0, sigma);
        size1 *= GLOB_component_size[rnd1];
        size2 *= GLOB_component_size[rnd2];
        size3 *= GLOB_component_size[rnd3];
        size4 *= GLOB_component_size[rnd4];
        bool_res = 0;
        if (size1*size2 >= size3*size4) 
        {
            bool_res = add_edge(rnd3,rnd4);
        } 
        else 
        {
            bool_res = add_edge(rnd1,rnd2);
        }
#endif
#ifdef EPES_MECH_iffs
        bool_res = 0;
        if (GLOB_component_name[rnd1] == GLOB_component_name[rnd2]) 
        {
            if (GLOB_component_name[rnd3] == GLOB_component_name[rnd4]) 
            {
                if (GLOB_component_size[rnd1] >= GLOB_component_size[rnd3] ) 
                {
                    bool_res = add_edge(rnd3,rnd4);
                } 
                else
                {
                    bool_res = add_edge(rnd1,rnd2);
                }
            }
            else
            {
                bool_res = add_edge(rnd1,rnd2);
            }
        }
        else
        {
            if (GLOB_component_name[rnd3] == GLOB_component_name[rnd4]) 
            {
                bool_res = add_edge(rnd3,rnd4);
            }
            else
            {
                size1 = (double) GLOB_component_size[rnd1];
                size2 = (double) GLOB_component_size[rnd2];
                size3 = (double) GLOB_component_size[rnd3];
                size4 = (double) GLOB_component_size[rnd4];
                if (size1*size2 >= size3*size4) 
                {
                    bool_res = add_edge(rnd3,rnd4);
                } 
                else 
                {
                    bool_res = add_edge(rnd1,rnd2);
                }
            }
        }
#endif
#ifdef EPES_MECH_selfloop
        name_aux1 = 0;
        double w = (double)name_aux1; // to get rid of 'unused' gcc flags
        rnd1 = (int)(Random()*NODE_NR);
        int rnd1_size = GLOB_component_size[rnd1];
        if ( (rnd1_size!=1 && rnd1_size!=2) && rnd1_size!=3 )
        {
            // the clustering condition is s.t. there exists teh possibility
            // of adding a node to the component, though the probability of
            // this not being the case is quite low for high enough 
            // clusters
            //if (localClustering(rnd1) < 1) 
            {
w = Random();
if (w > 10*GLOB_max_component_size/((double)NODE_NR)) 
{
    // we dont go outside our domain if this condition holds
    // we choose a random node within the domain and connect with it
    rnd2 = random_node_comp(rnd1);//GLOB_component_name[rnd1]);
    int count_aux=0;
    while ( exists_edge(rnd1,rnd2) )
    {
        rnd2 = random_node_comp(rnd1);
        count_aux++;
        if (count_aux > rnd1_size) {
            break;
        }
    }
    bool_res = add_edge(rnd1,rnd2);
}
else
{
    rnd2 = (int)(Random()*NODE_NR);
    name_aux1 = GLOB_component_name[rnd1];
    //while (GLOB_component_name[rnd2] == name_aux1) 
    {
        rnd2 = (int)(Random()*NODE_NR);
    } // now we have two nodes from different domains
    //printf("rnds1: %d, %d, %d \n", rnd1, rnd2, rnd3);
    rnd3 = (int)(Random()*NODE_NR);
    name_aux2 = GLOB_component_name[rnd2];
    //while ( (GLOB_component_name[rnd3] == name_aux1) ||
    //        (GLOB_component_name[rnd3] == name_aux2)    ) 
    {
        rnd3 = (int)(Random()*NODE_NR);
    }
    //printf("rnds2: %d, %d, %d \n", rnd1, rnd2, rnd3);
    rnd4 = (int)(Random()*NODE_NR);
    name_aux3 = GLOB_component_name[rnd3];
    //while (GLOB_component_name[rnd4] == name_aux3) 
    {
        rnd4 = (int)(Random()*NODE_NR);
    }
    size1 = (double) GLOB_component_size[rnd1];
    size2 = (double) GLOB_component_size[rnd2];
    size3 = (double) GLOB_component_size[rnd3];
    size4 = (double) GLOB_component_size[rnd4];
    bool_res = 0;
    if (size1*size2 >= size3*size4) 
    {
        bool_res = add_edge(rnd3,rnd4);
    } 
    else 
    {
        bool_res = add_edge(rnd1,rnd2);
    }
}
            } // clustering if
        } // component size if
        else
        {
            size1 = (double) GLOB_component_size[rnd1];
            size2 = (double) GLOB_component_size[rnd2];
            size3 = (double) GLOB_component_size[rnd3];
            size4 = (double) GLOB_component_size[rnd4];
            bool_res = 0;
            if (size1*size2 >= size3*size4) 
            {
                bool_res = add_edge(rnd3,rnd4);
            } 
            else 
            {
                bool_res = add_edge(rnd1,rnd2);
            }
        } // end size if
#endif





        if (bool_res == 1) { // this means an edge has been added
            GLOB_nr_edges++;

            //for (int i = 0; i < NODE_NR; i++) {
            //    GLOB_omega_nat[i] = (double) degree[i];//(double)GLOB_component_size[i];
            //            //0.5*(-1 + 2*Random());//sampleNormal();
            //}
        }

/*
        rnd1 = (int)(Random()*NODE_NR);
        // This takes edges at random without considering components.
        rnd2 = (int)(Random()*NODE_NR);
        while ( exists_edge(rnd1,rnd2) )    // it also checks for self loops,
                                            // where rnd1 == rnd2
        {
            rnd2 = (int)(Random()*NODE_NR);
        }
        
        rnd3 = (int)(Random()*NODE_NR);
        while (rnd3 == rnd2 || rnd3 == rnd1) 
        {
            rnd3 = (int)(Random()*NODE_NR); // we want one which is indep. from
                                            // the first edge
        }
        rnd4 = (int)Random()*NODE_NR;
        while ( exists_edge(rnd3, rnd4) ) 
        { // this also checks for rnd3==rnd4
            rnd4 = (int)(Random()*NODE_NR);
        }
        // I take the first three nodes to be each different from each other
        // and the last can be anything as long as its not the same as the
        // third. This is because the first and second form one edge, and the
        // third and fourth another. In order to have edges, 1st != 2nd,
        // 3rd!=4th, but one of the edges can of course share a node with the
        // other edge as long as its not two nodes, because then its the same 
        // edge.
        // printf("rnds: %d, %d, %d, %d \n", rnd1, rnd2, rnd3, rnd4);
        size1 = component_size[rnd1];
        size2 = component_size[rnd2];
        size3 = component_size[rnd3];
        size4 = component_size[rnd4];
        if (size1*size2 >= size3*size4) {
            bool_res = add_edge(rnd3,rnd4);
        } else {
            bool_res = add_edge(rnd1,rnd2);
        }
        if (bool_res == 1) { // this means an edge has been added
            GLOB_nr_edges++;
            printf("GLOB_nr_edges=%d \t max_comp_size=%g unique_elem:%g\n",
                                            GLOB_nr_edges,
                                            GLOB_max_component_size,
                                            (GLOB_unique_elements_in_network));
        }
*/
        
        //printf("nr_edges: %d \n", nr_edges);
    }
    return 1;
}

void increase_edges_FREQ_GAP (double t, int m, double alpha,
                                double tiempo, double sigma,
                                FILE *f_out, double *r_med, double *r_var)
{
    int tot_nr_edges = 0;
    tot_nr_edges =  (int) (t* ((double)NODE_NR));

    int rnd1 = 0;
    /*
    if (GLOB_nr_edges == 0) {
        int rnd2 = 0;
        rnd1 = (int)(Random()*NODE_NR); 
        rnd2 = (int)(Random()*NODE_NR); 
        while(rnd1 == rnd2){
            rnd2 = (int)(Random()*NODE_NR); 
        }
        add_edge(rnd1,rnd2);
        GLOB_nr_edges++;
    }
    */
#ifdef DEBUG
    if (m>NODE_NR || m<1) {
        printf("warning: m in frq.gap algorithm wrong\n");
        exit(45);
    }
#endif

    int bool_res = 0;
    int fg_nodes[m];
    //for(int i_idx = 0; i_idx < m; i_idx++) {
    //    fg_nodes[i_idx] = -1;
    //}
    

    double aux_rnd22=0;
    double aux_sum = 0;
    int iidx=0;
    int flag22=0;
    int flag32=1;
    int nd1, nd2, nd3, nd4;
    int size1, size2, size3, size4;
    int nd[FG_ACLIOPTAS_K][2];
    double nd_weff[FG_ACLIOPTAS_K][2];
    double achlioptas_min[FG_ACLIOPTAS_K];
    int achl_idx = 0;
    double achl_opt = 0;

    double rnd1_weff = 0;
    double weff_nd1, weff_nd2, weff_nd3, weff_nd4;
    double fg_nodes_weff[FG_M];

    double r_dom[NODE_NR];
    double weff_by_node[NODE_NR];
    double weff_by_domain[NODE_NR]; // initialized in weff_compt_efficient

    weff_compt_DOUBLY_efficient(weff_by_domain, weff_by_node,
                                r_dom, NODE_NR, r_med, r_var, t, sigma);
    
    double weff[NODE_NR];
    for (int i = 0; i < NODE_NR; i++) {
        weff[i] = weff_by_node[i];
    }
    

    while(GLOB_nr_edges <= tot_nr_edges)
    {
        GLOB_unique_components = GLOB_unique_elements_in_network;
        //if (GLOB_unique_components < 0 || GLOB_unique_components > NODE_NR) {
        //    printf("warning: unique components wrong!\n");
        //    exit(11);
        //}

        //rnd1 = (int)(Random()*NODE_NR);
        // esta suma se pod´ria evitar pero no quiero hacer una global solo para eso

        
        /*
        aux_rnd22 = Random();
        iidx = 0;
        aux_sum = exp(FG_ALPHA*fabs(GLOB_omega_nat[iidx]))/GLOB_sum_omega_nat;
        //printf("glob sum: %g\n", GLOB_sum_omega_nat);
        while (aux_sum < aux_rnd22) {
            iidx++;
            aux_sum += exp(FG_ALPHA*fabs(GLOB_omega_nat[iidx]))/GLOB_sum_omega_nat;
        }
        rnd1 = iidx;
        //printf("node: %d \n", rnd1);
        */
        #ifndef FG_ACLIOPTAS
        flag22 = generate_node_FREQUENCY_GAP(alpha,
&rnd1, &rnd1_weff, m, fg_nodes, fg_nodes_weff, tiempo, sigma, weff, NODE_NR);
        #endif
        #ifdef FG_ACLIOPTAS
if  (GLOB_unique_elements_in_network > 1) // once all connected, choice irrelevant
{
        // THIS ASUMES FG_M == 1 !!!!!!!!
        for (int i = 0; i < FG_ACLIOPTAS_K; i++) 
        {
            flag22 = generate_node_FREQUENCY_GAP(alpha,
&rnd1, &rnd1_weff, m, fg_nodes, fg_nodes_weff, tiempo, sigma, weff, NODE_NR);
            nd[i][0] = rnd1;
            nd[i][1] = fg_nodes[0];
            nd_weff[i][0] = rnd1_weff;
            nd_weff[i][1] = fg_nodes_weff[0];
        }
        /*
        printf("nodes: ");
        for (int i = 0; i < FG_ACLIOPTAS_K; i++)
        {
            printf("%d-%d ",nd[i][0],nd[i][1]);
        }   printf("\n");
        */

        achl_opt = +1e8;
        for (int i = 0; i < FG_ACLIOPTAS_K; i++) {
            achlioptas_min[i] = ( GLOB_component_size[ nd[i][0] ] *
                                  GLOB_component_size[ nd[i][1] ]   );
            if  ( (double)achlioptas_min[i] < achl_opt) {
                achl_opt = (double)achlioptas_min[i];
                achl_idx = i;
            }
        }
        /*
        printf("weights:");
        for (int i = 0; i < FG_ACLIOPTAS_K; i++)
        {
            printf("%g ", achlioptas_min[i]);
        }   printf("\n");
        printf("chosen: %d\n", achl_idx);
        */

        rnd1         = nd[achl_idx][0];
        rnd1_weff    = nd_weff[achl_idx][0];
        fg_nodes[0]  = nd[achl_idx][1];
        fg_nodes_weff[0] = nd_weff[achl_idx][1];

        // REMEMBER ALL THIS FOR FG_M=1
        //printf("chosen nodes: %d-%d\n", rnd1, fg_nodes[0]);
}
else
{
        flag22 = generate_node_FREQUENCY_GAP(alpha,
&rnd1, &rnd1_weff, m, fg_nodes, fg_nodes_weff, tiempo, sigma, weff, NODE_NR);
}
        #endif
        if (flag22!=0)
        {
            for(int i_idx = 0; i_idx < m; i_idx++) 
            {
                bool_res = add_edge(rnd1, fg_nodes[i_idx]);
                //printf("(%d,%d,%d)\n", rnd1, fg_nodes[i_idx], bool_res);
                if (bool_res == 1) 
                {
                    GLOB_nr_edges++;
                }
                flag32 = fprintf(f_out, "%d\t%d\t%d\t%g\t%g\t%g\t%g\n",
                                 GLOB_nr_edges, rnd1, fg_nodes[i_idx],
                                 rnd1_weff, fg_nodes_weff[i_idx],
                                 GLOB_omega_nat[rnd1],
                                 GLOB_omega_nat[fg_nodes[i_idx]]);
                if (flag32==0) {
                    printf("warning: fprintf failed for edgelist!\n");
                    exit(22);
                }
                
            }
            //printf("\n");
            // Now we thermalize to addapt to new edges.
#ifdef TERMALIZATION // we wait for r to stabilize
                for (int i = 0; i < TERMALIZATION; i++){
                    for (int t_aux = 0; t_aux < IN_BETWEEN; t_aux++) {
                        update_RK(tiempo, sigma, DELTA_T);
                    }
                }
#endif
        }
        flag22=0;
    }
    

}
