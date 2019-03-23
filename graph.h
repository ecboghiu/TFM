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
#define SCALE_FREE
//#define READ_NETWORK_FROM_FILE    //WARNING: as of time moment this
                                    // doesn't work well

// Define if you want a histogram of the degrees of the graph
//#define DEGREE_HISTOGRAM

// Number of nodes in the graph.
#define NODE_NR 1000
#define K_MAX 1000
#define K_MIN 2
#define AVG_NUMBER 1

//#define OSCILLATOR_ON
//#define PERCOLATION_ON


// For the scale-free probability distribution.
// Input this into wolfram alpha if you want ot know the gamma for a certain
// mean degree:
// 6 = 1/(\int x^(-y)dx for x=2..100) *\int x*x^(-y) dx for x=2..100
#define GAMMA 1.69285028
//2.27837758 1.69285028

#define M_PI 3.14159265358979323846

// Value not chosen arbitrarily, but so that theta_dot*h~1e-4,ie,
// sufficently. small
#define DELTA_T 1e-2
// How many times we measure.
#define MAX_STEPS 10000
// Number of updates in between measures.
#define IN_BETWEEN 10

// If its not defined we dont wait to termalize
// we need to wait around 4s
#define TERMALIZATION 40

//#define EPSILON_OSCILLATOR 1e-3

// only when termalization is undefined, when termalization is
// defined we make a look going through many sigma values, not one
#define SIGMA_VAL 0.4


// Where we store our connection matrix
extern int **C;

// Global array where degrees are stored.
extern int *degree;

// Where I will store the component it's connected to.
extern int *GLOB_component_name;
extern int *GLOB_component_size;
extern int GLOB_nr_edges;
extern int GLOB_unique_components;
extern double GLOB_max_component_size;
extern double GLOB_unique_elements_in_network;

// We will store the natural frequencies of all the oscillators.
extern double *GLOB_theta;
extern double *GLOB_omega_nat;

/*******FUNCTION PROTOTYPES*******/
void update_EULER   (double sigma, double h);
void update_RK      (double t, double sigma, double h);

int initEXPL_product_rule (double t);

void initOmegas     ();
void initThetas     ();

void init_C         (int ***data_ptr, int dim_x, int dim_y);
void init_C_memory  (int ***data_ptr, int dim_x, int dim_y);
void clear_C_memory (int ***data_ptr, int dim_x, int dim_y);
void initERmodel    (double prob);
void initScaleFree  ();

int add_edge        (int i, int j);
int remove_edge     (int I, int J);
int exists_edge     (int i, int j);
int unique_elements (int arr[], int len);

int  read           (int i, int j); // works like a[i][j] where a is
                                    // adjacency matris
void readAdj        (int l); // reads from file
void saveAdj        ();
void printAdj       ();
void print_vec      (int **vec, int size);
int int_max_vector  (int **vector, int size);
void print_C        ();
void write_C_to_file();
void saveAdjGephi   (); // writes edges to file in a format usable by Gephi
int number_of_edges (char* filename);

// Graph observables, measurables.
void   calculateDegree       (); // calculates array with degrees of all nodes
double calculateTheta_dot_i  (double t, double *phases, int phases_len,
                                double sigma, int i);
double calculateThetaAverage ();
double phase_coherence       ();
int debug_count_nodes        ();
int number_of_domains        ();
void domain_size_array       (int *size_array, int i, int* domain, int *size);
int domain_of_node_i_size    (int i);
int max_domain_size          ();

// Clustering
double localClustering  (int i);
double Clustering       ();

// Functions for generating random numbers.
double Random           (); // Random number in [0,1).
void   ini_ran          (int SEMILLA); // initializes the generator
int    generateDegree   (int m, double gamma, double norm_const);
double sampleNormal     (); // copied from http://stackoverflow.com/a/10645091

// Statistics functions.
void   Histogram        (int *data, double *Hist, int N_data, int N_Intervalos,
                         double *d, double *min, double *max);
void   med_var          (double *datos, int numero_de_datos,
                         double *media, double *varianza);
void   calcHist         (char *filename);

// Some other generic functions.
int wheel_sum           (signed int i, signed int aux, int wheel_len);
int in_pair_list        (int list[][2], int len, int i, int j);
int  exists_link        (int i, int j); // return 1 if there is link 
                                        // from j to i
void make_edge_list     ();

// Functions for lists
void append_to_list     (signed int *list, signed int elem);
int list_len            (signed int *list);
void remove_from_list   (signed int *list, signed int loc, int length);

#endif
