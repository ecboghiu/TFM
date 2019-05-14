#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> // For initializing the random generator.
#include <string.h> // For atoi(),atof() mainly, amongst other things

#define DEBUG

/*******CONTROL: type of graph in case of oscillator dynamc*******/
//#define SMALL_WORLD
//#define ERDOS_RENYI
#define ERDOS_RENYI_prob 0.006
//#define SCALE_FREE
//#define READ_NETWORK_FROM_FILE
#ifdef READ_NETWORK_FROM_FILE
    #define FILENAME_FROM_WHICH_TO_READ "nx_edgelist_BA.txt"
    #define NET_TYPE "growth"
    #define NET_CHARACT 0.0
#endif

//#define BARABASI_ALBERT
#ifdef BARABASI_ALBERT
    #define BA_PARAM_M 3
#endif
#define EPES
#ifdef EPES
    #define EPES_CHARACT "tribe"
#endif

//#define INITIAL_SEED 1557414566
#define INITIAL_SEED 1557414222


// Define if you want a histogram of the degrees of the graph
//#define DEGREE_HISTOGRAM

// Number of nodes in the graph.
#define NODE_NR 10
#define K_MAX 100
#define K_MIN 2
#define AVG_NUMBER 1

// Value not chosen arbitrarily, but so that theta_dot*h~1e-4,ie,
// sufficently. small
#define DELTA_T 1e-1
// How many times we measure.
#define MAX_STEPS 3000
// Number of updates in between measures.
#define IN_BETWEEN 0
#define SIGMA_MIN 0.0
#define SIGMA_MAX 0.2
#define NR_SIGMA 10

// only when termalization is undefined, when termalization is
// defined we make a look going through many sigma values, not one
#define SIGMA_VAL 1.5
//#define PRINT_EVOLUTION_OF_R

// If its not defined we dont wait to termalize
// we need to wait around 4s
#define TERMALIZATION 100

//#define WEFF_MEMORY_LINKED_LIST 1   // THIS IS OBSOLTE NOW
#define WEFF_MEMORY_DYNAMIC_MATRIX
#ifdef WEFF_MEMORY_DYNAMIC_MATRIX
    #define WEFF_MEMORY_DYNAMIC_MATRIX_INI_SIZE NODE_NR/2
#endif

//#define EPSILON_OSCILLATOR 1e-3

//#define OSCILLATOR_ON
//#define PERCOLATION_ON
//#define SYNC_AND_PERC_ON
//#ifdef SYNC_AND_PERC_ON
    //#define EPES_MECH_Pure_perc
    //#define EPES_MECH_compare_r
    //#define EPES_MECH_Scale_by_dom_size
    //#define EPES_MECH_weff
    //#define EPES_MECH_iffs
    #define EPES_MECH_selfloop
//#endif
#define FREQUENCY_GAP
//#ifdef FREQUENCY_GAP
    #define FG_M 1
    #define FG_ALPHA -1.0
    #define FG_T_MIN 0.0
    #define FG_T_MAX 1.2
    #define FG_T_NUMBER (FG_T_MAX-FG_T_MIN)*NODE_NR
    #define FG_WEFF_LOWER_FREQUENCY -1e8
    #define FG_WEFF_MAX_STEPS MAX_STEPS
    #define FG_WEFF_IN_BETWEEN IN_BETWEEN
//#endif


// For the scale-free probability distribution.
// Input this into wolfram alpha if you want ot know the gamma for a certain
// mean degree:
// 6 = 1/(\int x^(-y)dx for x=2..100) *\int x*x^(-y) dx for x=2..100
#define GAMMA 2.0842189
//2.27837758 1.69285028

#define M_PI 3.14159265358979323846



// Where we store our connection matrix. C[i][j] gives node i's j-th neighbor.
// there are degree[i] neighbors. Everything after the last number has value -1.
// All nodes are tagged with natural numbers from 0 to NODE_NR-1.
extern int **C;
extern int **C_dom;
extern int C_dom_sizes[NODE_NR];

extern int GLOB_initial_seed; // this is so that we store the seed we used
                         // and then we can access it to write it to file

// Global array where degrees are stored.
extern int *degree;

// GLOB_component_name[i] gives node i's component id.
extern int *GLOB_component_name;

// GLOB_component_name[i] gives node i's component size (neighbors of neighbors
// of neighbors etc.)
extern int *GLOB_component_size;
extern int GLOB_dom_array_lengths[NODE_NR];

// Gives the total number of edges in the graph.
extern int GLOB_nr_edges;

// Gives the number of different components in the graph. Serves as a measure
// of fragmentation.
extern int GLOB_unique_components;

// Gives the size of the biggest component of the graph. Serves as a measure of
// fragmentation.
extern double GLOB_max_component_size;
extern double GLOB_unique_elements_in_network;

// GLOB_dom_size[id_compt] gives the size of component of name id_compt
extern double GLOB_dom_size[NODE_NR];

// GLOB_theta[i] gives node i's current phase
extern double *GLOB_theta;
// GLOB_omega_nat[j] gives node j's natural frequency \omega_j
extern double *GLOB_omega_nat;

#define NormRANu (2.3283063671E-10F)
extern unsigned int  irr[256];
extern unsigned int  ir1;
extern unsigned char ind_ran, ig1, ig2, ig3;

/*******STRUCT DEFINITIONS*******/
struct _Node {
        int id;
        struct _Node *next;
        //struct _Node *next_domain; // ignore this, only useful for percolation
}; typedef struct _Node *Node;

//struct _Connected {
//        struct _Node *main_head;
//}; typedef struct _Connected *Graph_domains;

struct _Graph {
    Node *suc;
}; typedef struct _Graph *Graph;

struct _Graph_array {
    int **suc;
}; typedef struct _Graph_array *Graph_array;


//Graph_array GLOB_dom_array;
Graph GLOB_dom;

/*******FUNCTION PROTOTYPES*******/
void update_EULER   (double sigma, double h);

void update_RK      (double t, double sigma, double h);

int  initEXPL_product_rule      (double t, double sigma);
int  init_FREQ_GAP              (double t, double sigma);
void increase_edges_FREQ_GAP    (double t, int m, double alpha,
                                double tiempo, double sigma);

void oscillator_on      (void);
void percolation_on     (void);
void epes_on            (void);
void frequency_gap_on   (void);

// initializing various things
void initOmegas     (void);
void initThetas     (void);
void init_C         (int ***data_ptr, int dim_x, int dim_y);
void init_C_memory  (int ***data_ptr, int dim_x, int dim_y);
void init_glob_vect_memory (void);
void initDom        (void);
void clear_C_memory (int ***data_ptr, int dim_x, int dim_y);

void initERmodel    (double prob);
void initScaleFree  (void);
void init_BA        (int m, int N);

int  add_edge        (int i, int j);
int  remove_edge     (int I, int J);
int  exists_edge     (int i, int j);
void insertNode      (Node *I, int i);
void removeNode      (Node *I, int j);
int  delta_kron      (int i, int j);

int  read            (int i, int j); // works like a[i][j] where a is
                                    // adjacency matris
void readAdj                (int l); // reads from file
void saveAdj                (void);
void printAdj               (void);
void print_vec              (int *vec, int size);
int  int_max_vector         (int **vector, int size);
void copy_vector            (int *target, int *source,
                             int i_from, int len);
void print_C                (void);
void print_linked_list      (void);
void write_C_to_file        (void);
void saveAdjGephi           (int C_MAT[][K_MAX]);
void read_edgelist_file_py  (const char* filename);


// Graph observables, measurables for debugging and other
void    debug                   (void);
void    calculateDegree         (void);

//double  calculateTheta_dot_i    (double t, double *phases, int phases_len,
//                                double sigma, int i);

inline 
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
    //
    //double weight = 1.0;
    //if (degree[i] == 0)  {
    //    weight = 1.0;
    //}  else {
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



//double  diff_weff_weight        (double alpha, double wi, double wj);
inline
double diff_weff_weight(double alpha, double wi, double wj, double ri, double rj)
{
    return exp(alpha*rj*fabs(wi-wj));
}

double  calculateThetaAverage   (void);
double  weff_compt              (int id_compt, double t, double sigma);
void    weff_compt_efficient    (double *weff_dom, int weff_dom_size,
                                 double t, double sigma);
void    weff_compt_DOUBLY_efficient    (double *weff_dom, double *r_dom,
                                        int weff_dom_size,
                                        double t, double sigma);
double  weff_compt_instant      (int id_compt, double t, double sigma);
double  phase_coherence         (void);

//double  psi_coherence           (void);
inline
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
    //if  (Nrx < 1e-6){
    //    printf("warning: ry too small for division?");
    //    return atan2(Nry,Nrx);
    //}
    return atan2(Nry,Nrx);
}

double  phase_coherence_compt   (int id_compt);
double  psi_coherence_compt     (int id_compt);
int     debug_count_nodes       (void);
int     number_of_domains       (void);
void    domain_size_array       (int *size_array, int i, int* domain, int *size);
int     domain_of_node_i_size   (int i);
int     max_domain_size         (void);
int     join_domains            (int i, int j);
int     random_node_comp        (int id_compt);
int     number_of_edges         (char* filename);
int     unique_elements         (int arr[], int len);
// Clustering
double localClustering (int i);
double Clustering      (void);

// Functions for generating random numbers.
//double  Random                      (void); // Random number in [0,1).
inline
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

void    ini_ran                     (int SEMILLA); // initializes the generator
int     generateDegree              (int m, double gamma);
void    generate_node_BA            (int m, int* nodes);
void    generate_node_FREQUENCY_GAP (double alpha, int node_i, int m,
                                    int *nodes, double t, double sigma);

// copied from http://stackoverflow.com/a/10645091
double sampleNormal     (void); 

// Statistics functions.
void   Histogram        (int *data, double *Hist, int N_data, int N_Intervalos,
                        double *d, double *min, double *max);
void   med_var          (double *datos, int numero_de_datos,
                        double *media, double *varianza);
void   calcHist         (char *filename);

// Some other generic functions.
int  wheel_sum          (signed int i, signed int aux, int wheel_len);
int  in_pair_list       (int list[][2], int len, int i, int j);
int  exists_link        (int i, int j); // return 1 if there is link 
                                        // from j to i
void make_edge_list     (void);

// Functions for lists
void append_to_list     (signed int *list, signed int elem);
int  list_len           (signed int *list);
void remove_from_list   (signed int *list, signed int loc, int length);
void shuffle            (int *array, size_t n);


#endif
