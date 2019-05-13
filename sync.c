#include "graph.h"

double *GLOB_theta;
void initThetas()
{
    GLOB_theta = malloc(NODE_NR*sizeof(*GLOB_theta));
    for (int i = 0; i < NODE_NR; ++i) {
        GLOB_theta[i] = -M_PI + Random()*(2*M_PI);
    }
}

// Generates natural frequencies about 0 according to a Gaussian distribution of
// variance 1.
double *GLOB_omega_nat;
void initOmegas()
{
    // omega_nat IS GLOBAL
    GLOB_omega_nat = malloc(NODE_NR*sizeof(*GLOB_omega_nat));
    for (int i = 0; i < NODE_NR; i++) 
    {
        //GLOB_omega_nat[i] = sampleNormal();//(double)degree[i];//0.5*(-1 + 2*Random());//
        GLOB_omega_nat[i] = 0.5*(-1 + 2*Random());
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
/*
double psi_coherence()
{
    double Nrx, Nry;
    Nrx = Nry = 0;
    for (int i = 0; i < NODE_NR; i++) {
        Nrx += cos(GLOB_theta[i]);
        Nry += sin(GLOB_theta[i]);
    }
    //if  (Nrx < 1e-6){
    //    printf("warning: ry too small for division?");
    //    return atan2(Nry,Nrx);
    //}
    return atan2(Nry,Nrx);
}
*/

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

void weff_compt_DOUBLY_efficient (double *weff_dom, double *r_dom, int weff_dom_size,
                                double t, double sigma)
{
    
    double timp = 0;

    int aux_deg=0;
    double aux_time_average_theta_dot[NODE_NR];
    for (int i = 0; i < NODE_NR; i++)
    {
        aux_time_average_theta_dot[i]=0;
    }
    
    //double weff_by_domain[NODE_NR];
    for(int i = 0; i < weff_dom_size; i++)
    {
        weff_dom[i] = FG_WEFF_LOWER_FREQUENCY;
        r_dom[i] = -1;
    }
    //We need to do the first step and change from default

    /* 
        double weff = 0.0;
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

    // first hte time averages

    for(int i = 0; i < FG_WEFF_MAX_STEPS; i++)
    {
        /*
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
        */

        for (int m = 0; m < NODE_NR; m++)
        {
            aux_deg=degree[m];
            for (int n = 0; n < aux_deg; n++)
            {
                // this means that add current theta_dot_ij to the running avg
aux_time_average_theta_dot[m] = ((double)i)/(i+1) * aux_time_average_theta_dot[m]
            + 1.0/(i+1) * calculateTheta_dot_i(t, GLOB_theta, NODE_NR, sigma, m);
            }
        }

        for (int i_dom = 0; i_dom < NODE_NR; i_dom++)
        {
            if (C[i_dom][0] != 1)
            {
r_dom[i_dom] = ((double)i)/(i+1) * r_dom[i_dom]+ 1.0/(i+1) * phase_coherence_compt(i_dom);
            }
        }
        
        
        

        timp = i*DELTA_T*(1+FG_WEFF_IN_BETWEEN);
        for(int j = 0; j < FG_WEFF_IN_BETWEEN; j++) {
            update_RK(t+timp, sigma, DELTA_T); 
        }
        update_RK(t+timp, sigma, DELTA_T);
    }


    // after this we have the matrix of time averaged effective frequencies-
    // if space is a probelm use commented out version from above
    // now we need to do the "spatial" or component averaging 
/*
printf("tdot avg: ");
for (int i = 0; i < NODE_NR; i++)
{
    printf("%g ",aux_time_average_theta_dot[i]);
}   printf("\n");
*/


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
/*
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
*/
