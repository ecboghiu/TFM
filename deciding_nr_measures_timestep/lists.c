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