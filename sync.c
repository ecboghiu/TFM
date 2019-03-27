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
        GLOB_omega_nat[i] = sampleNormal();//(double)degree[i];//0.5*(-1 + 2*Random());//
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

double psi_coherence()
{
    double Nrx, Nry;
    Nrx = Nry = 0;

    for (int i = 0; i < NODE_NR; i++) {
        Nrx += cos(GLOB_theta[i]);
        Nry += sin(GLOB_theta[i]);
    }
    if  (Nrx < 1e-6){
        printf("warning: ry too small for division?");
        return atan2(Nry,Nrx);
    }
    return atan2(Nry,Nrx);
}

double phase_coherence_compt (int id_compt)
{

    int i=0;
    Node crawl = GLOB_dom->suc[id_compt];
    double rx, ry;
    rx = ry = 0;
    while(crawl != NULL) {
        i = crawl->id;
        rx += cos(GLOB_theta[i]);
        ry += sin(GLOB_theta[i]);
        
        crawl = crawl->next;
    }
    free(crawl); crawl = NULL;

    rx/=NODE_NR;
    ry/=NODE_NR;

    return sqrt(rx*rx+ry*ry);
}

double psi_coherence_compt (int id_compt)
{
    double Nrx, Nry;
    Nrx = Nry = 0;

    int i=0;
    Node crawl = GLOB_dom->suc[id_compt];
    while(crawl != NULL) {
        i = crawl->id;
        Nrx += cos(GLOB_theta[i]);
        Nry += sin(GLOB_theta[i]);
        
        crawl = crawl->next;
    }
    free(crawl); crawl = NULL;

    return atan2(Nry,Nrx);
}

double calculateTheta_dot_i(double t, double *phases, int phases_len,
                                            double sigma, int i)
{
    // \dot{theta_i}=\omega_i+|sigma\sum_{j=0}^{N} a_{ij}\sin{theta_j-theta_i}
    double sum = 0;
    //return cos(t);
    
    #ifdef EPSILON_OSCILLATOR
    double coupling_all_to_all = EPSILON_OSCILLATOR;
    for(int j = 0; j < phases_len; j++) {
        sum += coupling_all_to_all * sin(phases[j]-phases[i]);
    }
    #endif

    double coupling = 0;
    for(int j = 0; j < degree[i]; j++) {
        coupling = sigma;
        sum += coupling * sin( phases[ C[i][j] ] - phases[i] );
    }
    
    return GLOB_omega_nat[i] + sum;
}
