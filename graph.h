#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> // For initializing the random generator.
#include <string.h> // For atoi(),atof() mainly, amongst other things

#define DEBUG

/*******CONTROL: type of graph*******/
//#define SMALL_WORLD
//#define ERDOS_RENYI
#define ERDOS_RENYI_prob 0.006
//#define SCALE_FREE
#define READ_NETWORK_FROM_FILE
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


// Define if you want a histogram of the degrees of the graph
//#define DEGREE_HISTOGRAM

// Number of nodes in the graph.
#define NODE_NR 1000
#define K_MAX 100
#define K_MIN 2
#define AVG_NUMBER 1

// Value not chosen arbitrarily, but so that theta_dot*h~1e-4,ie,
// sufficently. small
#define DELTA_T 5e-2
// How many times we measure.
#define MAX_STEPS 1000
// Number of updates in between measures.
#define IN_BETWEEN 2
#define SIGMA_MIN 0.0
#define SIGMA_MAX 0.2
#define NR_SIGMA 10

// If its not defined we dont wait to termalize
// we need to wait around 4s
#define TERMALIZATION 40

//#define EPSILON_OSCILLATOR 1e-3

#define OSCILLATOR_ON
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
//#define FREQUENCY_GAP
//#ifdef FREQUENCY_GAP
    #define FG_M 1
    #define FG_ALPHA 0.5
    #define FG_T_MIN 0.0
    #define FG_T_MAX 1.0
    #define FG_T_NUMBER 100
    #define FG_WEFF_LOWER_FREQUENCY -1e6
    #define FG_WEFF_MAX_STEPS MAX_STEPS
    #define FG_WEFF_IN_BETWEEN IN_BETWEEN
//#endif



// only when termalization is undefined, when termalization is
// defined we make a look going through many sigma values, not one
#define SIGMA_VAL 1.0
#define PRINT_EVOLUTION_OF_R

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

// Global array where degrees are stored.
extern int *degree;

// GLOB_component_name[i] gives node i's component id.
extern int *GLOB_component_name;

// GLOB_component_name[i] gives node i's component size (neighbors of neighbors
// of neighbors etc.)
extern int *GLOB_component_size;

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
extern double *GLOB_dom_size;

// GLOB_theta[i] gives node i's current phase
extern double *GLOB_theta;
// GLOB_omega_nat[j] gives node j's natural frequency \omega_j
extern double *GLOB_omega_nat;

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
void print_vec              (int **vec, int size);
int  int_max_vector         (int **vector, int size);
void print_C                (void);
void print_linked_list      (void);
void write_C_to_file        (void);
void saveAdjGephi           (int C_MAT[][K_MAX]);
void read_edgelist_file_py  (const char* filename);


// Graph observables, measurables for debugging and other
void    debug                   (void);
void    calculateDegree         (void);
double  calculateTheta_dot_i    (double t, double *phases, int phases_len,
                                double sigma, int i);
double diff_weff_weight         (double alpha, double wi, double wj);
double  calculateThetaAverage   (void);
double  weff_compt              (int id_compt, double t, double sigma);
double  phase_coherence         (void);
double  psi_coherence           (void);
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
double  Random                      (void); // Random number in [0,1).
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
