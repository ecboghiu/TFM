/*
***Model of complex network synchronization
***Prácticas BIFI Julio 2016
***Emanuel Cr. Boghiu
***embogh@gmail.com
*/

/*IMPORTANT:*/
/*READ graph.h FOR CONTROLS*/

#include "graph.h"

int main(int args_number, char* args[])
{
    printf("Initializing program.\n");

    ini_ran(time(NULL));  // Seed for the random generator.
    //ini_ran(662323);     // we want predictable results

    init_C_memory(&C, NODE_NR, K_MAX);
    initDom();
/*
    join_domains(2,3);
    join_domains(0,3);
    join_domains(0,3);
    join_domains(0,2);
    join_domains(3,0);

    add_edge(2,3);
    add_edge(3,3);
    add_edge(2,3);
    add_edge(0,1);
    add_edge(0,3);

    print_linked_list();

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
        
    }

        int m=3;
    int *nodes = calloc(m, sizeof *nodes);
    generate_node_BA(m, nodes);
    printf("generated nodes:\n");
    for(int i = 0; i < m; i++) {
        printf("%d", nodes[i]);
    }
*/




////////////////////////// PERCOLATION ////////////////////////////////////////
#ifdef PERCOLATION_ON
percolation_on();
#endif
////////////////////////// OSCILLATOR /////////////////////////////////////////
#ifdef OSCILLATOR_ON
oscillator_on();
#endif
////////////////////////// EPES ///////////////////////////////////////////////
#ifdef SYNC_AND_PERC_ON
epes_on();
#endif // endif SYNC_AND_PERC_ON
//////////////////////////////////////////////////////////////////////////////

    // Freeing other global arrays.
    free(GLOB_component_name); GLOB_component_name = NULL;
    free(GLOB_component_size); GLOB_component_size = NULL;

    clear_C_memory(&C, NODE_NR, K_MAX);

    printf("\nYou've reached the end without dying! ;-)\n");
    return 0;
}
