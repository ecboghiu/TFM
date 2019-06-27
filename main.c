/*
***Emanuel Cr. Boghiu
***embogh@gmail.com
*/

/*IMPORTANT:*/
/*READ graph.h FOR CONTROLS*/

#include "graph.h"

int GLOB_initial_seed;

//int main(int args_number, char* args[])
int main(void)
{
    printf("\nInitializing program.\n\n");

#ifdef INITIAL_SEED
    GLOB_initial_seed = INITIAL_SEED;
#else
    GLOB_initial_seed = time(NULL);
#endif
    ini_ran(GLOB_initial_seed);     
    
    init_C_memory(&C, NODE_NR, K_MAX);
    init_glob_vect_memory();
    initDom();

    //debug();

int tiempo_ini = time(NULL);

FILE *theta_file = fopen("pattern.txt","w");
fclose(theta_file);

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

int tiempo_fin = time(NULL);

    printf("\nYou've reached the end without dying! ;-) (runtime: %g min)\n",
            (tiempo_fin-tiempo_ini)/60.0);
    return 0;
}
