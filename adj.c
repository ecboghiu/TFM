#include "graph.h"

double GLOB_max_component_size;
double GLOB_unique_elements_in_network;
int *degree, *GLOB_component_size, *GLOB_component_name;
int GLOB_nr_edges, GLOB_unique_components;
int **C;
int **C_dom;
double GLOB_dom_size[NODE_NR];
int GLOB_dom_array_lengths[NODE_NR];
int C_dom_sizes[NODE_NR];

void init_C_memory(int ***data_ptr, int dim_x, int dim_y)
{
    int i,j,k;
    
    int **data;
    data = (int**)malloc(sizeof(*data) * dim_x);
    if (data==NULL)
    {
        printf("warning: malloc gives bad results...\n");
        exit(111);
    }
    
    int *ax2;
    for (k = 0; k < dim_x; k++) {
        ax2 = (int*)malloc(sizeof (*(data[k])) * dim_y);
        data[k] =  ax2;
        if (data[k]==NULL)
        {
            printf("warning: malloc gives bad results...\n");
            exit(111);
        }
    }
    for (i = 0; i < dim_x; i++) {
        for (j = 0; j < dim_y; j++) {
            data[i][j] = -1; // -1 means no connection
        }
    }
    *data_ptr = data;

#ifdef WEFF_MEMORY_DYNAMIC_MATRIX
    int col_nr = WEFF_MEMORY_DYNAMIC_MATRIX_INI_SIZE;
    C_dom =  malloc(NODE_NR * sizeof *C_dom);
    if (C_dom==NULL)
    {
        printf("warning: malloc gives bad results...\n");
        exit(111);
    }

    for (k = 0; k < dim_x; k++) {
        //ax = (int *)malloc(col_nr * sizeof *(C_dom[k]) );
        //C_dom[k] = ax;
        C_dom[k] = malloc(col_nr * sizeof *(C_dom[k]) );
                
        //printf("apple's address = %p\n", (void*)(C_dom[k]));
        if (C_dom[k]==NULL)
        {
            printf("warning: malloc gives bad results...\n");
            exit(111);
        }
    }
    for (i = 0; i < dim_x; i++) {
        for (j = 0; j < col_nr; j++) {
            C_dom[i][j] = -1;
            if(j==0)
            {
                C_dom[i][j] = i;
            }
        }
    }

    for (i = 0; i < dim_x; i++)
    {
        C_dom_sizes[i] = col_nr;
    }
#endif
}

void init_glob_vect_memory(void)
{    
    int i=0;
    GLOB_unique_elements_in_network = NODE_NR;

    degree = (int*)malloc(NODE_NR * sizeof *degree);
    if (degree==NULL)
        {
            printf("warning: malloc gives bad results...\n");
            exit(111);
        }
    else
    {
        for (i = 0; i < NODE_NR; i++)
        {
            degree[i]=0;
        }
        
    }
    
    

    GLOB_component_size = (int*)malloc(NODE_NR * sizeof *GLOB_component_size);
    if (GLOB_component_size==NULL)
    {
        printf("warning: malloc gives bad results...\n");
        exit(111);
    } 
    for(i = 0; i < NODE_NR; i++){
        GLOB_component_size[i] = 1; // all nodes are their own cluster
    }
    GLOB_component_name = (int*)malloc(NODE_NR* sizeof *GLOB_component_name);
    if (GLOB_component_name==NULL)
        {
            printf("warning: malloc gives bad results...\n");
            exit(111);
        }
    for(i = 0; i < NODE_NR; i++){
        GLOB_component_name[i] = i; // all nodes are their own cluster
    }
    
}

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

    if (new_name != old_name) {
        join_domains(old_name,new_name);
/*        
printf("add edge: %d %d\n", i, j);
printf("weff matrix:\n");
for (int i_idx = 0; i_idx < NODE_NR; i_idx++)
{
    printf("%d:%d: ",i_idx,C_dom_sizes[i_idx]);
    for (int j_idx = 0; j_idx < C_dom_sizes[i_idx]; j_idx++)
    {
        printf("%d ", C_dom[i_idx][j_idx]);
    }
    printf("\n");
    
}

print_linked_list();
*/
    }
    

    int new_size = 0;
    if (new_name == old_name) // then they are from the same component
    {  
        new_size = GLOB_component_size[i]; // the size DOESN'T CHANGE
    }
    else
    {
         // diff comps sum
        new_size = GLOB_component_size[i] + GLOB_component_size[j];
        //printf("new_size: %d, %d %d", new_size, i, j);
        (GLOB_unique_elements_in_network)--;
        //printf("names: %d %d %d\n", new_name, old_name, new_size);
    }

    if (new_size > GLOB_max_component_size) {
        GLOB_max_component_size = (double)new_size;
        //printf("max components: %f of %d\n", GLOB_max_component_size, NODE_NR);
        if (GLOB_max_component_size > NODE_NR) 
        {
            printf("warning: too big a component\n");
            //exit(9);
        }
    }
    
    
    
    for(int i_aux = 0; i_aux < NODE_NR; i_aux++)
    {
        if (GLOB_component_name[i_aux] == old_name) {
            //GLOB_dom_size[GLOB_component_name[i_aux]] = new_size;
            GLOB_component_name[i_aux]  = new_name;
            GLOB_component_size[i_aux]  = new_size;
        }
        if (GLOB_component_name[i_aux] == new_name) {
            //GLOB_dom_size[GLOB_component_name[i_aux]] = new_size;
            GLOB_component_size[i_aux]  = new_size;
        }
    }

    // we join i to j
    
    if (old_name != new_name) {
        GLOB_dom_size[old_name]=0;
        GLOB_dom_size[new_name]=new_size;
    }
    else
    {
        GLOB_dom_size[old_name]=new_size;
        GLOB_dom_size[new_name]=new_size;
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

    for(int k = 0; k < degree[i]; k++)
    {
        if ( C[i][k] == j) {
            //printf("warning: edge already exists\n");
            return 1;
        }
    }
    return 0;
}

void read_edgelist_file_py (const char* filename)
{
    // Note, you should have previously initialized C.
    // Take care with how many nodes the network has!

    FILE *f_in;
    f_in = fopen(filename, "r");
    if (f_in == NULL) {
        printf("Could not open file for edgelist!\n");
        exit(11);
    } 
    else {
        int i,j;
        i = j = 0;
        #ifdef READ_NETWORK_FROM_FILE
        printf("Beginning adding edges from %s\n", FILENAME_FROM_WHICH_TO_READ);
        #else
        printf("Beginning adding edges from file.\n");
        #endif

        while ( fscanf(f_in, "%d %d {}\n", &i, &j) != EOF   ) {
            add_edge(i,j);
        }
        printf("Last added edge: %d %d\n", i,j);
        fclose(f_in); f_in = NULL;
    }
}




void initDom()
{
    GLOB_dom      =    (Graph) malloc(sizeof(*GLOB_dom));
    GLOB_dom->suc = (Node*)calloc(NODE_NR,sizeof(*(GLOB_dom->suc)));
    for(int i = 0; i < NODE_NR; i++)
    {
        insertNode(&(GLOB_dom->suc[i]),i);
    }
    
    /*
    GLOB_dom_array = malloc(sizeof(*GLOB_dom_array));
    GLOB_dom_array->suc = calloc(NODE_NR,sizeof(*(GLOB_dom_array->suc)));
    for (size_t i = 0; i < NODE_NR; i++)
    {
        GLOB_dom_array->suc[i] = calloc(1,sizeof *GLOB_dom_array->suc[i]);
    }
    for (size_t i = 0; i < NODE_NR; i++)
    {
        for (size_t j = 0; j < 10; j++)
        {
            GLOB_dom_array->suc[i][j] = -1;
        }
    }
    for (size_t i = 0; i < NODE_NR; i++)
    {
            GLOB_dom_array->suc[i][0] = i;
    }

    for (size_t i = 0; i < NODE_NR; i++)
    {
        GLOB_dom_array_lengths[i] = 1;
    }
    */
  

    //GLOB_dom_size = malloc(NODE_NR *sizeof *GLOB_dom_size);
    for(int i = 0; i < NODE_NR; i++)
    {
        GLOB_dom_size[i] = 1;
    }
    
}

void clear_C_memory(int ***data_ptr, int dim_x, int dim_y)
{
    //int a = dim_y;
    for (int i = 0; i < dim_x; i++) {
        free((*data_ptr)[i]);
        (*data_ptr)[i] = NULL;
    }
    free(*data_ptr);
}

// TODO: MAKE A FUNC TO ASSIGN MEMORY; MAKE ANOTHER TO CLEAR VALUES!!!!!!
// Taken from:: https://stackoverflow.com/questions/11463455/
// c-programming-initialize-2d-array-dynamically
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

    GLOB_max_component_size = 1;
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
    int i, j, sum, len;
    i = j = sum = len = 0;

    int k_min = K_MIN;  // to generate connected net with prob 1,
                    // see PHYS.REVIEW E71,027103(2005)
    //int k_max = NODE_NR-1;


    printf("Beginning generating degree distribution.\n");
    int *aux_deg;
    aux_deg = (int*)calloc(NODE_NR, sizeof *aux_deg);
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
        for (i = 0; i < NODE_NR; i++) {
            aux_deg[i]= generateDegree(k_min,GAMMA);
            #ifdef DEBUG
            //printf("just generated a degree: %d\n", aux_deg[i]);
            #endif
        }
        sum = 0;
        for (i = 0; i < NODE_NR; i++) {
            sum += aux_deg[i];
        }
    }
    /* We will use a list of ints whose elements are nodes. The element with
    value -1 marks the END OF THE LIST. The functions are made with this
    in mind. The length of the list is \sum k_i. We repeat the node according
    to the degree number present. We then pair nodes at random and if we
    pair them and form a link we remove the two nodes from the list.
    */

    printf("Finished generating degree distribution.\n");
    signed int *list = (int* )calloc((sum+1),sizeof(*list)); // sum+1 because we 
    list[0] = -1;                                     // need to store -1
    for (signed int i_idx = 0; i_idx < NODE_NR; i_idx++) {
        for (j = 0; j < aux_deg[i_idx]; j++) {
            append_to_list(list, i_idx);
        }
    }
    len = list_len(list);
    #ifdef DEBUG
    printf("%s len: %d sum: %d\n", "Check1.", len, sum);
    #endif
    printf("Finished declaring list.\n");

    printf("Shuffling list.\n");
    shuffle(list, len);
    printf("Finished shuffling.\n");

    printf("Entering nasty loop.\n");
    while (len!=0)
    { //TODO: SOLVE LIST GOES TO LENGTH 1, DOESN?T MAKE SENSe
        i = (int) (Random()*(len));
        while(list[i]==-1) {
            printf("warning: you got list[i]==-1 which shouldnt happen!\n");
            i = (int) (Random()*(len));
        }
        j = (int) (Random()*(len));
        while(list[j]==-1) {
            printf("warning: you got list[i]==-1 which shouldnt happen!\n");
            j = (int) (Random()*(len));
        }

        while(exists_edge(list[i],list[j]))
        {
            //i = (int) (Random()*(len));
            //j = (int) (Random()*(len));
            j = (j+1)%len;
        }
/*
        i=0; j=0;
        if (len < 80) {
            for(int i_1 = 0; i_1 < len; i_1++) {
                for(int j_1 = 0; j_1 < len; j_1++) {
                    if (! (exists_edge(list[i_1],list[j_1]))) {
                        i=i_1;
                        j=j_1;
                        break;
                    } else
                    {
                        printf("exists link\n");
                    }
                    
                    
                }
                if (! (exists_edge(list[i_1],list[j]))) {
                        i=i_1;
                        break;
                }
            }
        }
*/

        //printf("i,j: %d %d\n", i, j);
        //if ( !( exists_edge(list[i],list[j]) ) )
        {
            add_edge(list[i], list[j]);

            #ifdef DEBUG
            //printf("%s%d\t%s%d%s%d\t (%d,%d)\n","initScaleFree: list_len --",
            //            len, " i=", i, ",j=", j, list[i], list[j]);
            #endif

            if (j>i) {
                remove_from_list(list, j, sum+1);
                remove_from_list(list, i, sum+1);
            }
            if (i>=j) {
                remove_from_list(list, i, sum+1);
                remove_from_list(list, j, sum+1);
            }
            #ifdef DEBUG
            printf("%s len: %d sum: %d\n", "Check1.", len, sum);
            #endif
            len = list_len(list);
            if (len < 60) {
                for(int i_4 = 0; i_4 < len; i_4++)
                {
                    printf("%d ", list[i_4]);
                }   printf("\n");
            }
            

            
        }
    }
    printf("Finished nasty loop.\n");

    // Check if correct
    for(int i_5 = 0; i_5 < NODE_NR; i_5++)
    {
        if (degree[i_5] != aux_deg[i_5]) {
            printf("warning: node %d does not have the degree it should!\n",
                            i_5);
            exit(7);
        }
        
    }
    free(aux_deg);
    free(list); 
}


/*
This function genereated a Barabasi Albert network.
We assume that N is greater than 2.
*/
void init_BA (int m, int N)
{
    if (N<2) {
        printf("warning: N<2 in Ba algorithm\n");
        exit(45);
    }

    if (m>N || m<1) {
        printf("warning: m in Ba algorithm wrong\n");
        exit(45);
    }
    

    add_edge(0,1); // We need to start with some edges

    int *ba_nodes = (int*)calloc(m, sizeof *ba_nodes);
    for(int node_id = 2; node_id < N; node_id++)
    {
        generate_node_BA(m, ba_nodes);
        for(int i = 0; i < m; i++) {
            add_edge(node_id, ba_nodes[i]);
        }
        
    }
    free(ba_nodes);
}

// This saves an edges table suitable for the open source program
// called Gephi.
void saveAdjGephi (int C_MAT[][K_MAX])
{
    int i, j;
    i = j = 0;
    FILE *net;
    net = fopen("net_edges_Gephi.txt","w");
    if (net != NULL) {
        fprintf(net, "%s\n", "Source,Target,Type");
        for (i = 0; i < NODE_NR; i++) {
            for (j = 0; j < K_MAX; j++) {
                if (C_MAT[i][j] != 0) {
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

void print_vec (int *vec, int size)
{
    for(int i = 0; i < size; i++)
    {
        printf("%d    ", (int)vec[i]);
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

// Inserting a node is the same as inserting a new link.
void insertNode (Node *I, int i)
{
    // The new node is the new head. *I is a pointer to the old
    // head of the linked list.
    // We create a new aux pointer and make it hold the location
    // of the old head struct
    Node aux;
    aux = *I; // s is now a pointer to a struct of type Node
    // Now we make *I point to new blank memory and give it id
    // value i and make it point to the old head: head insertion.
    (*I) = (Node)malloc( sizeof *(*I) );
    (*I)->id=i;
    (*I)->next=aux;
}

void removeNode (Node *I, int j)
{
    Node crawl = *I;
    Node crawled_before = *I;
    if (j == ((*I)->id)) { // SPECIAL CASE: if we have to remove the head
        crawl = (*I);
        if ((*I)->next == NULL) { // There is only the head in the list.
            (*I) = NULL;
        } else { // There is more than one element in the linked list.
            (*I) = (*I)->next;
        }
        free(crawl);
    } else {
        while (crawl != NULL) {
            if (crawl->id == j) {
                (crawled_before->next) = (crawl->next); // we skip node j
                free(crawl); // we free the memory;
                return;
            }
            crawled_before = crawl;
            crawl = crawl->next;
        }
    }
}

int join_domains(int dom_i, int dom_j)
{
    
    if (dom_i == dom_j) {
        //printf("warning: no self domain joining allowed");
        return 0;
    }

#ifdef WEFF_MEMORY_DYNAMIC_MATRIX
// first without realloc, assume size is NODE_NR
//int idx_old = GLOB_dom_size[dom_j];
//int length_i = GLOB_dom_size[dom_i];

int size_i_old = GLOB_dom_size[dom_i];
int size_j_new = GLOB_dom_size[dom_j];


/*
if (length_i+idx_old >= C_dom_sizes[dom_j])
{
    int new_size = 2*(length_i+idx_old);
    int *aux_ptr = NULL;
    aux_ptr = (int *)realloc(C_dom[dom_j], new_size * sizeof(*aux_ptr)  );
    if  (aux_ptr==NULL)
    {
        free(aux_ptr);
        printf("warning: error using realloc!\n");
        exit(222);
    }
    else
    {
        C_dom[dom_j] = aux_ptr;
        
        for (int i = idx_old ; i < new_size; i++)
        {
            C_dom[dom_j][i]=-1;
        }
        C_dom_sizes[dom_j] = new_size;
    }
}
*/


if(C_dom[dom_j][0] == -1)
{
    printf("warning: an inexisting component should not have been selected!\n");
    exit(111);
}

int old_capacity = C_dom_sizes[dom_j];
int new_size = size_i_old + size_j_new;
int new_capacity = 2*new_size;
// https://stackoverflow.com/questions/3850749/
// does-realloc-overwrite-old-contents
if ( new_size >= old_capacity)
{
    //printf("sizenew, capacity: %d %d\n", new_size, old_capacity);
    int *oldpointer = C_dom[dom_j];
    int *newpointer = realloc(oldpointer, new_capacity*sizeof(int));
    if (newpointer == NULL) {
        printf("warning: problems with pointer\n");
        exit(11);
    } else {
        /* everything ok                                                                 */
        /* `newpointer` now points to a new memory block with the contents of oldpointer */
        /* `oldpointer` points to an invalid address                                     */
        oldpointer = newpointer;
        C_dom[dom_j] = newpointer;
        
        /* oldpointer points to the correct address                                */
        /* the contents at oldpointer have been copied while realloc did its thing */
        /* if the new size is smaller than the old size, some data was lost        */
        
    }
    C_dom_sizes[dom_j]=new_capacity;
}

for (int i = 0; i < size_i_old; i++)
{
    C_dom[dom_j][size_j_new + i] = C_dom[dom_i][i];
}
int *old_pointer2 = C_dom[dom_i];
int *new_pointer2 = realloc(old_pointer2, 1*sizeof(int));
if (new_pointer2 == NULL)
{
    printf("warning: problems with pointer\n");
    exit(11);
}
else
{
    old_pointer2 = new_pointer2;
    C_dom[dom_i] = new_pointer2;
    C_dom[dom_i][0] = -1;
    C_dom_sizes[dom_i] = 1;
}



for (int i = new_size; i < C_dom_sizes[dom_j]; i++)
{
    C_dom[dom_j][i] = -1;
}



/*
int actual_size = C_dom_sizes[dom_j];
if (size_j_new+size_i_old > actual_size)
{
    realloc(C_dom[dom_j], 100*sizeof(int));
    C_dom_sizes[dom_j]=100;
}
*/

/*
if //(size_j_new+size_i_old > C_dom_sizes[dom_j])
{
    printf("step1 domi %d domj %d\n", dom_i, dom_j);
    int size_lb = size_j_new+size_i_old;
    int aux_vect[NODE_NR];
    //copy_vector(aux_vect, C_dom[dom_i], 0, C_dom_sizes[dom_i]);
    for (int i = 0; i < C_dom_sizes[dom_i]; i++) {
        aux_vect[i] = C_dom[dom_i][i];
    }
    //for (int i = 0; i < NODE_NR; i++)
    {
        printf("apple's address = %p %p dom = %d\n", (void*)C_dom[dom_i], (void*)(C_dom_original)[dom_i], dom_i);
    }
    
    print_vec(aux_vect, C_dom_sizes[dom_i]);

    //free((C_dom_original)[0]);

    // Now we will try to realoc dom_j.
    if (C_dom[dom_j]!=NULL) {
        printf("step2\n");
        //int *aux_ptr2 = (void*)C_dom[dom_i];
        //free(aux_ptr2);
        printf("step3\n");
    }
    
    printf("step4\n");
    //C_dom[dom_i] =  (int*) malloc(2*size_lb * sizeof(int));
    int* tmp =  realloc(C_dom[dom_j], 2*size_lb * sizeof(int));
    if (tmp) {
        C_dom[dom_j] = tmp;
    } else
    {
        printf("warning: malloc failed! \n");
        exit(1111);
    }
    printf("step4\n");

    // Now we copy from i onto newly resized j
    for (int i = 0; i < size_i_old; i++) {
        C_dom[dom_j][size_j_new + i] = aux_vect[i];
        C_dom[dom_i][i] = -1;
    }
    for (int i = 0; i < size_lb; i++)
    {
        C_dom[dom_j][size_lb + i] = -1; 
    }
    
    C_dom_sizes[dom_i] = 1;
    C_dom_sizes[dom_j] = 2*size_lb;
}

*/
/*
if(C_dom[dom_i][0] != -1) //this should probably be uncommented
{
    //printf("dom_i dom_j: %d %d\n", dom_i, dom_j);
    //printf("idx_old length_i: %d %d\n", idx_old, length_i);
    for (int i = 0; i < length_i; i++)
    {
        C_dom[dom_j][idx_old+i]=C_dom[dom_i][i];
        
    }
    C_dom[dom_i][0]=-1;
    //C_dom[dom_i] = realloc(C_dom[dom_i], sizeof(int));
}
*/
#endif

#ifdef WEFF_MEMORY_LINKED_LIST
    Node crawl_j_end = GLOB_dom->suc[dom_j];
    if (crawl_j_end != NULL) {
        while( (crawl_j_end->next) != NULL) {
            crawl_j_end = crawl_j_end->next;
        }
        crawl_j_end->next = GLOB_dom->suc[dom_i];
        GLOB_dom->suc[dom_i] = NULL;
        return 1;
    }
    else
    {
        if (GLOB_dom->suc[dom_i] != NULL) {
            GLOB_dom->suc[dom_j] = GLOB_dom->suc[dom_i];
            GLOB_dom->suc[dom_i] = NULL;
            return 1;
        }
        else
        {
            return 1; // both are null so we are joining inexisting domains
        }
    }
#endif

return 1;


    /*
    int size_j, size_i;
    size_i = GLOB_dom_size[dom_i];
    size_j = GLOB_dom_size[dom_j];

    int max_size = GLOB_dom_array_lengths[dom_j];
    if (size_j+size_i > max_size)
    {
        GLOB_dom_array->suc[dom_j] = 
                    realloc(GLOB_dom_array->suc[dom_j], 2*max_size);
        GLOB_dom_array_lengths[dom_j] = 2*max_size;
    }
    
    
    for (int i = 0; i < GLOB_dom_size[dom_i]; i++)
    {
        GLOB_dom_array->suc[dom_j][size_j+i]=GLOB_dom_array->suc[dom_i][i];
        GLOB_dom_array->suc[dom_i][i] = -1;
    }
    

    return 1;
*/


}

// Returns 1 if there is a link from j to i.
int read (int i, int j)
{
    // This should work as reading from an adjacency matrix though slower.
    if (i==j) {
        return 0;
    }
    
    int deg = degree[i];
    for(int k = 0; k < deg; k++)
    {
        if (C[i][k] == j) {
            return 1;
        }
    }
    return 0;
}

double localClustering (int i)
{
    int nr_neigh = degree[i];
    double max_edges_in_G_i = 1.0/2 * nr_neigh * (nr_neigh-1);
    if (nr_neigh == 0) {
        return 0;
    }
    if (nr_neigh == 1) {
        return 0;
    }
   
   /* double sum = 0;
    for (int j = 0; j < NODE_NR; j++) {
        for (int m = 0; m < NODE_NR; m++) {
            sum += read(i,j)*read(j,m)*read(m,i);
        }
    }
    return (1.0/2)*sum/max_edges_in_G_i;
*/
    double nr_edges = 0;
    int node_j = -1;
    int node_k = -1;
    
    for(int j = 0; j < nr_neigh; j++)
    {
        node_j = C[i][j];
        for(int k = 0; k < nr_neigh; k++) {
            node_k = C[i][k];
            //printf("node j,k: %d %d\n", node_j, node_k);
            if ( read(node_j,node_k) ) {
                nr_edges++;
                //printf("nr edges: %g\n", nr_edges);
            }
        }
    }
    nr_edges /= 2.0;
    nr_edges /= max_edges_in_G_i;
    //printf("nr_edges: %g", nr_edges);

    return nr_edges;

}

// Average of the local clustering over the whole graph.
double Clustering ()
{
    double sum = 0;
    for (int i = 0; i < NODE_NR; i++) {
        sum += localClustering(i);
    }
    //printf("sum: %g\n", sum);
    return sum/NODE_NR;
}

int random_node_comp (int node)
{
    int rnd_idx = (int) ( Random() * GLOB_component_size[node] );
    //printf("rnd_idx/tot=%d %d\n",rnd_idx,GLOB_component_size[node]);
    Node crawl = GLOB_dom->suc[GLOB_component_name[node]];
    int i=0;
    int counter=0;
    while(crawl != NULL) {
        i = crawl->id;
        if (counter == rnd_idx) {
            break;
        }
        counter++;
        crawl = crawl->next;
        //printf("coiunter=%d %d %d %d\n",counter, rnd_idx, i, GLOB_component_size[node] );
    }
    if (crawl == NULL) 
    {
        // this means for some reason it didn't stop when it should
        return -1;
    }

    return i;
}