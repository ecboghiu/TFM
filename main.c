/*
***Emanuel Cr. Boghiu
***embogh@gmail.com
*/

/*IMPORTANT:*/
/*READ graph.h FOR CONTROLS*/

#include "graph.h"

//int main(int args_number, char* args[])
int main(void)
{
    printf("Initializing program.\n");

    //ini_ran(time(NULL));  // Seed for the random generator.
    ini_ran(662323);     // we want predictable results

    init_C_memory(&C, NODE_NR, K_MAX);
    initDom();

    //debug();

////////////////////////// PERCOLATION ////////////////////////////////////////
#ifdef PERCOLATION_ON
//percolation_on();
#endif
////////////////////////// OSCILLATOR /////////////////////////////////////////
#ifdef OSCILLATOR_ON
    oscillator_on();
#endif
////////////////////////// EPES ///////////////////////////////////////////////
#ifdef SYNC_AND_PERC_ON
    epes_on();
#endif
////////////////////////// FREQUENCY GAP///////////////////////////////////////
#ifdef FREQUENCY_GAP
    frequency_gap_on();
#endif
////////////////////////////////////////////////////////////////////////////////

    // Freeing other global arrays.
    //free(GLOB_component_name); GLOB_component_name = NULL;
    //free(GLOB_component_size); GLOB_component_size = NULL;

    //clear_C_memory(&C, NODE_NR, K_MAX);

    printf("\nYou've reached the end without dying! ;-)\n");
    return 0;
}
