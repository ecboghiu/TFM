#include "graph.h"

int initEXPL_product_rule (double t, double sigma)
{   
    int tot_nr_edges = 0;
    tot_nr_edges =  (int) (t* ((double)NODE_NR));
    //printf("tot_nr_edges=%d\n", tot_nr_edges);
    
    int rnd1, rnd2, rnd3, rnd4;
    double size1, size2, size3, size4;
    size1 = size2 = size3 = size4 = 0;
    rnd1 = rnd2 = rnd3 = rnd4 = 0;
    int bool_res = 0; // 0 false, 1 true
    int name_aux1, name_aux2, name_aux3;
    name_aux1 = name_aux2 = name_aux3 = 0;

    
    while (GLOB_nr_edges < tot_nr_edges) 
{
    //GLOB_unique_components = unique_elements(GLOB_component_name, NODE_NR);
    GLOB_unique_components = (int) (GLOB_unique_elements_in_network);
    if (GLOB_unique_components < 0 || GLOB_unique_components > NODE_NR) {
        printf("warning: unique components wrong!\n");
        exit(11);
    }
    
    //printf("unique components: %d\n", GLOB_unique_components);
  //  if (GLOB_unique_components == 1) {
  //      return 0;
  //  } 
    /*  
    if (GLOB_unique_components == 2) {
        rnd1 = (int)(Random()*NODE_NR);
        rnd2 = 0;
        name_aux1 = GLOB_component_name[rnd1];
        while (GLOB_component_name[rnd2] == name_aux1) {
            rnd2++;
            if (rnd2 > NODE_NR) {
                printf("warning: rnd2>NODE_NR");
                return 0;
            }
            
        }
        add_edge(rnd1,rnd2);
        GLOB_nr_edges++;
        return 1;
    }
    */
    
    /*
    if (GLOB_unique_components == 3) {
        rnd1 = (int)(Random()*NODE_NR);
        rnd2 = 0;
        name_aux1 = GLOB_component_name[rnd1];
        while (GLOB_component_name[rnd2] == name_aux1) {
            rnd2++;
        }
        name_aux2 = GLOB_component_name[rnd2];
        rnd3 = 0;
        while ( (GLOB_component_name[rnd3] == name_aux1) ||
                (GLOB_component_name[rnd3] == name_aux2)    ) 
        {
            rnd3++;
        }
        //printf("rnds33: %d, %d, %d \n", rnd1, rnd2, rnd3);
        size1 = component_size[rnd1];
        size2 = component_size[rnd2];
        size3 = component_size[rnd3];
        if (size2 >= size3) {
            bool_res = add_edge(rnd1,rnd3);
        } else {
            bool_res = add_edge(rnd1,rnd2);
        }
        if (bool_res == 1) { // this means an edge has been added
            GLOB_nr_edges++;
            //printf("GLOB_nr_edges=%d\n", GLOB_nr_edges);
        }
        return 1;
    }*/

/*
        if (GLOB_nr_edges >= NODE_NR) {
            printf("IF STATEMENT:Glob_nr_edges= %d\n", GLOB_nr_edges);
            //exit(12);
        }
*/
        rnd1 = (int)(Random()*NODE_NR);
        rnd2 = (int)(Random()*NODE_NR);
        rnd3 = (int)(Random()*NODE_NR);
        rnd4 = (int)(Random()*NODE_NR);

/*
        name_aux1 = 0;
        rnd1 = (int)(Random()*NODE_NR);
        // TODO: remove thjese conditions
        rnd2 = (int)(Random()*NODE_NR);
        name_aux1 = GLOB_component_name[rnd1];
        while (GLOB_component_name[rnd2] == name_aux1) 
        {
            rnd2 = (int)(Random()*NODE_NR);
        } // now we have two nodes from different domains
        //printf("rnds1: %d, %d, %d \n", rnd1, rnd2, rnd3);
        rnd3 = (int)(Random()*NODE_NR);
        name_aux2 = GLOB_component_name[rnd2];
        while ( (GLOB_component_name[rnd3] == name_aux1) ||
                (GLOB_component_name[rnd3] == name_aux2)    ) 
        {
            rnd3 = (int)(Random()*NODE_NR);
        }


        //printf("rnds2: %d, %d, %d \n", rnd1, rnd2, rnd3);
        rnd4 = (int)(Random()*NODE_NR);
        name_aux3 = GLOB_component_name[rnd3];
        while (GLOB_component_name[rnd4] == name_aux3) 
        {
            rnd4 = (int)(Random()*NODE_NR);
        }
        //printf("rnds3: %d, %d, %d \n", rnd1, rnd2, rnd3);
*/

#ifdef EPES_MECH_Pure_perc
        size1 = (double) GLOB_component_size[rnd1];
        size2 = (double) GLOB_component_size[rnd2];
        size3 = (double) GLOB_component_size[rnd3];
        size4 = (double) GLOB_component_size[rnd4];
        bool_res = 0;
        if (size1*size2 >= size3*size4) 
        {
            bool_res = add_edge(rnd3,rnd4);
        } 
        else 
        {
            bool_res = add_edge(rnd1,rnd2);
        }
#endif
#ifdef EPES_MECH_compare_r
        size1 = phase_coherence_compt(GLOB_component_name[rnd1]);
        size2 = phase_coherence_compt(GLOB_component_name[rnd2]);
        size3 = phase_coherence_compt(GLOB_component_name[rnd3]);
        size4 = phase_coherence_compt(GLOB_component_name[rnd4]);
        bool_res = 0;
        if (size1*size2 >= size3*size4) 
        {
            bool_res = add_edge(rnd3,rnd4);
        } 
        else 
        {
            bool_res = add_edge(rnd1,rnd2);
        }
#endif
#ifdef EPES_MECH_Scale_by_dom_size
        size1 = phase_coherence_compt(GLOB_component_name[rnd1]);
        size2 = phase_coherence_compt(GLOB_component_name[rnd2]);
        size3 = phase_coherence_compt(GLOB_component_name[rnd3]);
        size4 = phase_coherence_compt(GLOB_component_name[rnd4]);
        size1 /= GLOB_component_size[rnd1];
        size2 /= GLOB_component_size[rnd2];
        size3 /= GLOB_component_size[rnd3];
        size4 /= GLOB_component_size[rnd4];
        bool_res = 0;
        if (size1*size2 >= size3*size4) 
        {
            bool_res = add_edge(rnd3,rnd4);
        } 
        else 
        {
            bool_res = add_edge(rnd1,rnd2);
        }
#endif
#ifdef EPES_MECH_weff
        size1 = weff_compt(GLOB_component_name[rnd1], 0, sigma);
        size2 = weff_compt(GLOB_component_name[rnd2], 0, sigma);
        size3 = weff_compt(GLOB_component_name[rnd3], 0, sigma);
        size4 = weff_compt(GLOB_component_name[rnd4], 0, sigma);
        size1 *= GLOB_component_size[rnd1];
        size2 *= GLOB_component_size[rnd2];
        size3 *= GLOB_component_size[rnd3];
        size4 *= GLOB_component_size[rnd4];
        bool_res = 0;
        if (size1*size2 >= size3*size4) 
        {
            bool_res = add_edge(rnd3,rnd4);
        } 
        else 
        {
            bool_res = add_edge(rnd1,rnd2);
        }
#endif
#ifdef EPES_MECH_iffs
        bool_res = 0;
        if (GLOB_component_name[rnd1] == GLOB_component_name[rnd2]) 
        {
            if (GLOB_component_name[rnd3] == GLOB_component_name[rnd4]) 
            {
                if (GLOB_component_size[rnd1] >= GLOB_component_size[rnd3] ) 
                {
                    bool_res = add_edge(rnd3,rnd4);
                } 
                else
                {
                    bool_res = add_edge(rnd1,rnd2);
                }
            }
            else
            {
                bool_res = add_edge(rnd1,rnd2);
            }
        }
        else
        {
            if (GLOB_component_name[rnd3] == GLOB_component_name[rnd4]) 
            {
                bool_res = add_edge(rnd3,rnd4);
            }
            else
            {
                size1 = (double) GLOB_component_size[rnd1];
                size2 = (double) GLOB_component_size[rnd2];
                size3 = (double) GLOB_component_size[rnd3];
                size4 = (double) GLOB_component_size[rnd4];
                if (size1*size2 >= size3*size4) 
                {
                    bool_res = add_edge(rnd3,rnd4);
                } 
                else 
                {
                    bool_res = add_edge(rnd1,rnd2);
                }
            }
        }
#endif
#ifdef EPES_MECH_selfloop
        name_aux1 = 0;
        double w = (double)name_aux1; // to get rid of 'unused' gcc flags
        rnd1 = (int)(Random()*NODE_NR);
        int rnd1_size = GLOB_component_size[rnd1];
        if ( (rnd1_size!=1 && rnd1_size!=2) && rnd1_size!=3 )
        {
            // the clustering condition is s.t. there exists teh possibility
            // of adding a node to the component, though the probability of
            // this not being the case is quite low for high enough 
            // clusters
            //if (localClustering(rnd1) < 1) 
            {
w = Random();
if (w > 10*GLOB_max_component_size/((double)NODE_NR)) 
{
    // we dont go outside our domain if this condition holds
    // we choose a random node within the domain and connect with it
    rnd2 = random_node_comp(rnd1);//GLOB_component_name[rnd1]);
    int count_aux=0;
    while ( exists_edge(rnd1,rnd2) )
    {
        rnd2 = random_node_comp(rnd1);
        count_aux++;
        if (count_aux > rnd1_size) {
            break;
        }
    }
    bool_res = add_edge(rnd1,rnd2);
}
else
{
    rnd2 = (int)(Random()*NODE_NR);
    name_aux1 = GLOB_component_name[rnd1];
    //while (GLOB_component_name[rnd2] == name_aux1) 
    {
        rnd2 = (int)(Random()*NODE_NR);
    } // now we have two nodes from different domains
    //printf("rnds1: %d, %d, %d \n", rnd1, rnd2, rnd3);
    rnd3 = (int)(Random()*NODE_NR);
    name_aux2 = GLOB_component_name[rnd2];
    //while ( (GLOB_component_name[rnd3] == name_aux1) ||
    //        (GLOB_component_name[rnd3] == name_aux2)    ) 
    {
        rnd3 = (int)(Random()*NODE_NR);
    }
    //printf("rnds2: %d, %d, %d \n", rnd1, rnd2, rnd3);
    rnd4 = (int)(Random()*NODE_NR);
    name_aux3 = GLOB_component_name[rnd3];
    //while (GLOB_component_name[rnd4] == name_aux3) 
    {
        rnd4 = (int)(Random()*NODE_NR);
    }
    size1 = (double) GLOB_component_size[rnd1];
    size2 = (double) GLOB_component_size[rnd2];
    size3 = (double) GLOB_component_size[rnd3];
    size4 = (double) GLOB_component_size[rnd4];
    bool_res = 0;
    if (size1*size2 >= size3*size4) 
    {
        bool_res = add_edge(rnd3,rnd4);
    } 
    else 
    {
        bool_res = add_edge(rnd1,rnd2);
    }
}
            } // clustering if
        } // component size if
        else
        {
            size1 = (double) GLOB_component_size[rnd1];
            size2 = (double) GLOB_component_size[rnd2];
            size3 = (double) GLOB_component_size[rnd3];
            size4 = (double) GLOB_component_size[rnd4];
            bool_res = 0;
            if (size1*size2 >= size3*size4) 
            {
                bool_res = add_edge(rnd3,rnd4);
            } 
            else 
            {
                bool_res = add_edge(rnd1,rnd2);
            }
        } // end size if
#endif





        if (bool_res == 1) { // this means an edge has been added
            GLOB_nr_edges++;

            //for (int i = 0; i < NODE_NR; i++) {
            //    GLOB_omega_nat[i] = (double) degree[i];//(double)GLOB_component_size[i];
            //            //0.5*(-1 + 2*Random());//sampleNormal();
            //}
        }

/*
        rnd1 = (int)(Random()*NODE_NR);
        // This takes edges at random without considering components.
        rnd2 = (int)(Random()*NODE_NR);
        while ( exists_edge(rnd1,rnd2) )    // it also checks for self loops,
                                            // where rnd1 == rnd2
        {
            rnd2 = (int)(Random()*NODE_NR);
        }
        
        rnd3 = (int)(Random()*NODE_NR);
        while (rnd3 == rnd2 || rnd3 == rnd1) 
        {
            rnd3 = (int)(Random()*NODE_NR); // we want one which is indep. from
                                            // the first edge
        }
        rnd4 = (int)Random()*NODE_NR;
        while ( exists_edge(rnd3, rnd4) ) 
        { // this also checks for rnd3==rnd4
            rnd4 = (int)(Random()*NODE_NR);
        }
        // I take the first three nodes to be each different from each other
        // and the last can be anything as long as its not the same as the
        // third. This is because the first and second form one edge, and the
        // third and fourth another. In order to have edges, 1st != 2nd,
        // 3rd!=4th, but one of the edges can of course share a node with the
        // other edge as long as its not two nodes, because then its the same 
        // edge.
        // printf("rnds: %d, %d, %d, %d \n", rnd1, rnd2, rnd3, rnd4);
        size1 = component_size[rnd1];
        size2 = component_size[rnd2];
        size3 = component_size[rnd3];
        size4 = component_size[rnd4];
        if (size1*size2 >= size3*size4) {
            bool_res = add_edge(rnd3,rnd4);
        } else {
            bool_res = add_edge(rnd1,rnd2);
        }
        if (bool_res == 1) { // this means an edge has been added
            GLOB_nr_edges++;
            printf("GLOB_nr_edges=%d \t max_comp_size=%g unique_elem:%g\n",
                                            GLOB_nr_edges,
                                            GLOB_max_component_size,
                                            (GLOB_unique_elements_in_network));
        }
*/
        
        //printf("nr_edges: %d \n", nr_edges);
    }
    return 1;
}

void increase_edges_FREQ_GAP (double t, int m, double alpha,
                                double tiempo, double sigma)
{
    int tot_nr_edges = 0;
    tot_nr_edges =  (int) (t* ((double)NODE_NR));

    int rnd1 = 0;
    if (GLOB_nr_edges == 0) {
        int rnd2 = 0;
        rnd1 = (int)(Random()*NODE_NR); 
        rnd2 = (int)(Random()*NODE_NR); 
        while(rnd1 == rnd2){
            rnd2 = (int)(Random()*NODE_NR); 
        }
        add_edge(rnd1,rnd2);
        GLOB_nr_edges++;
    }

    if (m>NODE_NR || m<1) {
        printf("warning: m in frq.gap algorithm wrong\n");
        exit(45);
    }

    int bool_res = 0;
    int fg_nodes[m];
    //for(int i_idx = 0; i_idx < m; i_idx++) {
    //    fg_nodes[i_idx] = -1;
    //}
    
    while(GLOB_nr_edges <= tot_nr_edges)
    {
        GLOB_unique_components = GLOB_unique_elements_in_network;
        //if (GLOB_unique_components < 0 || GLOB_unique_components > NODE_NR) {
        //    printf("warning: unique components wrong!\n");
        //    exit(11);
        //}

        rnd1 = (int)(Random()*NODE_NR);
        generate_node_FREQUENCY_GAP(alpha, rnd1, m, fg_nodes, tiempo, sigma);
        for(int i_idx = 0; i_idx < m; i_idx++) 
        {
            bool_res = add_edge(rnd1, fg_nodes[i_idx]);
            //printf("(%d,%d,%d)\n", rnd1, fg_nodes[i_idx], bool_res);
            if (bool_res == 1) 
            {
                GLOB_nr_edges++;
            }
        }
        //printf("\n");
        // Now we thermalize to addapt to new edges.
#ifdef TERMALIZATION // we wait for r to stabilize
            for (int i = 0; i < TERMALIZATION; i++)
                for (int t_aux = 0; t_aux < IN_BETWEEN; t_aux++)
                    update_RK(tiempo, sigma, DELTA_T);
#endif
    }
    

}