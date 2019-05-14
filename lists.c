#include "graph.h"

// This list ends in -1.
// WARNING: list should not be longer than number of nodes to avoid going over
// assigned space.
void append_to_list (signed int * list, signed int elem)
{
    int i = 0;
    while (list[i]!=-1) {
        i++;
    }
    list[i]   = elem;
    list[i+1] = -1;
}

int list_len (signed int * list)
{
    int end = 0;
    while (list[end]!=-1) {
        end++;
    }
    return end;
}


// !!! CARE NOT TO REMOVE THE -1 AT THE END. EVERYTHING WILLL EXPLOOO0000000OODE!
void remove_from_list (signed int * list, signed int loc, int length)
{
    /*
    int end = 0;
    while (list[end]!=-1) {
        end++;
    } // now end points to the end of the list
    if (loc > end) {
        printf("You want to remove something out of the list!\n");
    } else if ( loc == end) {
        #ifdef DEBUG
        printf("YOU SHOULD NOT REMOVE THE END!\n");
        #endif
        return;
    }
    */
    int i = loc;
    while (list[i] != -1 && i<length) {
        list[i] = list[i+1];
        i++;
    }
}

// taken form https://stackoverflow.com/questions/6127503/shuffle-array-in-c
void shuffle(int *array, size_t n)
{
    if (n > 1) 
    {   
        int t = 0;
        size_t i, j;
        for (i = 0; i < n - 1; i++) 
        {
            j = i + rand() / (RAND_MAX / (n - i) + 1);
            t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}

// Funcion para imprimir en pantalla la red con el fin de poder debugear mas
// facilmente. No coge parametros porque usa las variables globales.
void print_linked_list()
{
    Node aux_inside = NULL;
    printf("\n");
    for(int i = 0; i < NODE_NR; i++)     
    {
        aux_inside = GLOB_dom->suc[i];
        printf("fila %d [ --", i);
        while ((aux_inside) != NULL) {
            printf(" %d -", aux_inside->id);
            aux_inside = aux_inside->next;
        }
        if (aux_inside == NULL) {
            printf("- NULL -");
        }
        printf("- ]");
        printf("\n");
    }
}

void copy_vector (int *target, int *source, int i_from, int len)
{
    for (int i = i_from; i < len; i++) {
        target[i] = source[i];
    }    
}