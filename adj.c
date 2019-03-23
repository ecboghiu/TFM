#include "graph.h"

double GLOB_max_component_size;
double GLOB_unique_elements_in_network;
int *degree, *GLOB_component_size, *GLOB_component_name;
int GLOB_nr_edges, GLOB_unique_components;
int **C;

int add_edge (int i, int j) {
    if (i==j) {
        //printf("warning: no self-loops!\n");
        return 0;
    }

    if (exists_edge(i,j)) {
        //printf("warning: there exists an edge already!\n");
        return 0;
    }
    

    C[ i ][ degree[i] ] = j;
    C[ j ][ degree[j] ] = i;
    degree[i]++;
    degree[j]++;
    
    int new_name = GLOB_component_name[j];
    int old_name = GLOB_component_name[i];
    int new_size = 0;
    if (new_name == old_name) // then they are from the same component
    {  
        new_size = GLOB_component_size[i]; // the size DOESN'T CHANGE
    }
    else
    {
        new_size = GLOB_component_size[i] + GLOB_component_size[j]; // diff comps sum
        (GLOB_unique_elements_in_network)--;
        //printf("names: %d %d\n", new_name, old_name);
    }

    if (new_size > GLOB_max_component_size) {
        GLOB_max_component_size = new_size;
        //printf("max components: %d of %d\n", GLOB_max_component_size, NODE_NR);
        if (GLOB_max_component_size > NODE_NR) {
            printf("warning: too big a component");
            exit(9);
        }
    }
    
    
    for(int i_aux = 0; i_aux < NODE_NR; i_aux++)
    {
        if (GLOB_component_name[i_aux] == old_name) {
            GLOB_component_name[i_aux]  = new_name;
            GLOB_component_size[i_aux]  = new_size;
        }
        if (GLOB_component_name[i_aux] == new_name) {
            GLOB_component_size[i_aux]  = new_size;
        }
    }



    return 1;
}

int remove_edge(int I, int J) {
    for(int j = 0; j < K_MAX; j++) {
        if (C[I][j] == J) {
            for(int i = j; i < K_MAX-j-1; i++) {
                C[I][i] = C[I][i+1];
            }
            degree[I]--;
            break;
        }
    }
    for(int j = 0; j < K_MAX; j++) {
        if (C[J][j] == I) {
            for(int i = j; i < K_MAX-j-1; i++) {
                C[J][i] = C[J][i+1];
            }
            degree[J]--;
            return 1;
        }
    }
    return 0; // false if no edge to remove
}

int exists_edge(int i, int j)
{
    if (i==j) {
        //printf("warning: no self-loops!\n");
        return 1;
    }

    int k_max = degree[i];
    for(int k = 0; k < k_max; k++)
    {
        if ( C[i][k] == j) {
            //printf("warning: edge already exists\n");
            return 1;
        }
    }
    return 0;
}


void init_C_memory(int ***data_ptr, int dim_x, int dim_y)
{
    int i,j,k;
    
    int **data;
    data = (int **) malloc(sizeof(int *) * dim_x);
    for (k = 0; k < dim_x; k++) {
        data[k] = (int *) malloc(sizeof(int) * dim_y);
    }
    for (i = 0; i < dim_x; i++) {
        for (j = 0; j < dim_y; j++) {
            data[i][j] = -1; // -1 means no connection
        }
    }
    *data_ptr = data;

    
    GLOB_max_component_size = -1;
    GLOB_unique_elements_in_network = NODE_NR;
    degree         = calloc(NODE_NR, sizeof *degree);
    GLOB_component_size = malloc(NODE_NR * sizeof *GLOB_component_size);
    for(int i = 0; i < NODE_NR; i++){
        GLOB_component_size[i] = 1; // all nodes are their own cluster
    }
    GLOB_component_name = malloc(NODE_NR* sizeof *GLOB_component_name);
    for(int i = 0; i < NODE_NR; i++){
        GLOB_component_name[i] = i; // all nodes are their own cluster
    }

    GLOB_nr_edges = 0;
    GLOB_unique_components = NODE_NR;
}

void clear_C_memory(int ***data_ptr, int dim_x, int dim_y)
{
    for (int i = 0; i < dim_x; i++) {
        free((*data_ptr)[i]);
        (*data_ptr)[i] = NULL;
    }
    free(*data_ptr); *data_ptr = NULL;
}

// TODO: MAKE A FUNC TO ASSIGN MEMORY; MAKE ANOTHER TO CLEAR VALUES!!!!!!
// Taken from:: https://stackoverflow.com/questions/11463455/c-programming-initialize-2d-array-dynamically
void init_C(int ***data_ptr, int dim_x, int dim_y) 
{
    for (int i = 0; i < dim_x; i++) {
        for (int j = 0; j < dim_y; j++) {
            (*data_ptr)[i][j] = -1; // -1 means no connection
        }
    }

    for(int i = 0; i < NODE_NR; i++)
    {
        degree[i] = 0;
        GLOB_component_size[i] = 1;
        GLOB_component_name[i] = i;
    }

    GLOB_max_component_size = -1;
    GLOB_unique_elements_in_network = NODE_NR;

    GLOB_nr_edges = 0;
    GLOB_unique_components = NODE_NR;
}

void initERmodel(double prob)
{
    for (int i = 0; i < NODE_NR-1; i++) {
        for (int j = i+1; j < NODE_NR; j++) { // so that we dont repeat edges
            if (Random()<prob) {
                add_edge(i,j);
            }
        }
    }
}

/* This should be run after initAdj to ensure a null initial graph. */
void initScaleFree ()
{
    int i, j, sum, len, k_min, k_max;
    double norm_const, sum_norm;
/*
    k_min = K_MIN;  // to generate connected net with prob 1,
                    // see PHYS.REVIEW E71,027103(2005)
    k_max = NODE_NR-1;

    // normalization constant of p~k^-GAMMA
    sum_norm = 0;
    for (int k = k_min; k <= k_max ; k++) {
        sum_norm += pow(k,-GAMMA);
    }
    norm_const = 1.0/sum_norm;
    #ifdef DEBUG
    printf("%s norm const: %g\n", "Check1.", norm_const);
    #endif
*/

    printf("Beginning generating degree distribution.\n");
    int *aux_deg;
    aux_deg = calloc(NODE_NR, sizeof *aux_deg);
    if (aux_deg == NULL)
    {
        printf("warning: null memory allocation\n");
        exit(10);
    }
    sum = 1; // so we start the first loop
    while (sum%2 != 0)  // requirement of the algorithm: \sum k_i even and less
                        // edges than than the maximum possible in
                        // an undirected graph
    {
         printf("Beginning generating degree distribution.\n");
        for (i = 0; i < NODE_NR; i++) {
            aux_deg[i]= generateDegree(k_min,GAMMA,norm_const);
            #ifdef DEBUG
            printf("just generated a degree: %d\n", aux_deg[i]);
            #endif
        }
        sum = 0;
        for (i = 0; i < NODE_NR; i++) {
            sum += aux_deg[i];
        }
    }
    /* We will use a list of ints whose elements are nodes. The element with value
        -1 marks the END OF THE LIST. The functions are made with this in mind.
        The length of the list is \sum k_i. We repeat the node according to the degree
        number present. We then pair nodes at random and if we pair them and form a
        link we remove the two nodes from the list.
    */

    printf("Finished generating degree distribution.\n");
    signed int *list = calloc((sum+1),sizeof(*list)); // sum+1 because we need
    list[0] = -1;                                     // to store -1
    for (signed int i = 0; i < NODE_NR; i++) {
        for (j = 0; j < aux_deg[i]; j++) {
            append_to_list(list, i);
        }
    }
    len = list_len(list);
    #ifdef DEBUG
    printf("%s len: %d sum: %d\n", "Check1.", len, sum);
    #endif
    printf("Finished declaring list.\n");

    
    printf("Entering nasty loop.\n");
    while (len!=0)
    { //TODO: SOLVE LIST GOES TO LENGTH 1, DOESN?T MAKE SENSe
        i = (int) (Random()*(len));
        j = (int) (Random()*(len));
        while (list[i]==-1 || list[j]==-1){
            printf("warning: you got list[i]==-1 which shouldnt happen!\n");
            i = (int) (Random()*(len));
            j = (int) (Random()*(len));
        }
        
        if ( !( exists_edge(list[i],list[j]) ) )
        {
            add_edge(list[i], list[j]);

            //#ifdef DEBUG
            //printf("%s%d\t%s%d%s%d\t (%d,%d)\n","initScaleFree: list_len --",
            //            len, " i=", i, ",j=", j, list[i], list[j]);
            //#endif
            //len = list_len(list);
            #ifdef DEBUG
            printf("%s len: %d sum: %d\n", "Check1.", len, sum);
            #endif

            if (j>i) {
                remove_from_list(list, j, sum+1);
                remove_from_list(list, i, sum+1);
            }
            if (i>j) {
                remove_from_list(list, i, sum+1);
                remove_from_list(list, j, sum+1);
            }
            len = list_len(list);
        }
    }
    printf("Finished nasty loop.\n");

    // Check if correct
    for(int i = 0; i < NODE_NR; i++)
    {
        if (degree[i] != aux_deg[i]) {
            printf("warning: node %d does not have the degree it should!\n",
                            i);
            exit(7);
        }
        
    }
    free(aux_deg); aux_deg = NULL;
    free(list); list = NULL;
}


// This saves an edges table suitable for the open source program
// called Gephi.
void saveAdjGephi (int C[][K_MAX])
{
    int i, j;
    i = j = 0;
    FILE *net;
    net = fopen("net_edges_Gephi.txt","w");
    if (net != NULL) {
        fprintf(net, "%s\n", "Source,Target,Type");
        for (i = 0; i < NODE_NR; i++) {
            for (j = 0; j < K_MAX; j++) {
                if (C[i][j] != 0) {
                    fprintf(net, "%d,%d,Undirected\n", i, j);
                }
            }
        }
        printf("%s\n", "Successfully saves in Gephi format!");
    }
    else
    {
        printf("%s\n", "Error opening file, saveAdjGephi");
    }
}


void calculateDegree()
{
    /*
    Node app = malloc(sizeof(*app));
    free(degree);
    degree = malloc(NODE_NR*sizeof(*degree));
    for (int i = 0; i < NODE_NR; i++) {
        degree[i] = 0;
        app = adj->suc[i];
        while (app!=NULL) {
            degree[i]++;
            app = app->next;
        }
    }
    free(app);
    */
}

// Taken from: https://stackoverflow.com/questions/28556641/
// a-program-that-counts-unique-elements-in-an-un-ordered-array
int unique_elements(int arr[], int len) {

    int counted[len], j, n, count, flag;

    counted[0] = arr[0]; 

    count = 1;/*one element is counted*/

        for(j=0; j <= len-1; ++j) {
        flag = 1;;
        /*the counted array will always have 'count' elements*/
        for(n=0; n < count; ++n) {
            if(arr[j] == counted[n]) {
                flag = 0;
            }
        }
        if(flag == 1) {
            ++count;
            counted[count-1] = arr[j];
        }
    }
    return count;
}

void print_vec (int **vec, int size)
{
    for(int i = 0; i < size; i++)
    {
        printf("%d    ", (int)(*vec)[i]);
    }
    printf("\n");
    
}


void print_C ()
{
    int i,j;
    for(i=0;i<NODE_NR;i++) {
        for(j=0;j<K_MAX;j++) {
            printf("%d\t", C[i][j]);
        }
        printf("\t\t %d \t %d \t %d", degree[i], GLOB_component_name[i],
                                                    GLOB_component_size[i]);
        printf("\n");
    }
    printf("\n");
}

void write_C_to_file()
{
    FILE *f = fopen("adj_C.txt","w");
    for(int i = 0; i < NODE_NR; i++)
    {
        fprintf(f, "%d ", i);
        for(int j = 0; j < degree[i]; j++)
        {
            //if (C[i][j] > i) // this is if you don't want repeated edges 
            {
                fprintf(f, "%d ", C[i][j]);
            }
        }
        fprintf(f, "\n");
    }
    
}

int int_max_vector (int **vector, int size) {
    int aux = 0;
    for(int i = 0; i < size; i++) {
        if ((*vector)[i] >= aux) {
            aux  = (*vector)[i];
        }
    }
    return aux;
}
