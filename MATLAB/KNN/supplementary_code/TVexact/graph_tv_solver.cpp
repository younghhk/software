#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>


#include "graph.h"

extern "C" {
  void graph_tv(double *,int,int,int*,int*,double,float);
  void graph_tv_weighted(double *,int,int,int*,int*,double*,double,float);
}



//solve min 0.5 \|y-x\|^2 + \lambda \|WD x\|_1
//where D is the incidence matrix, W is the diagonal matrix of edge weights
void graph_tv_weighted(double *Y,// value of nodes
        int n, //number of nodes
        int m, // number of edges
        int* edges1, 
        int* edges2,
//an array of edges of size m. There is an edge edges1[i] -> edges2[i]
        double* edge_weights,
        double lambda,
        float erreur)
{
    
    int i;
    double dval,dr,drd;
    
    dr=lambda;
    
    Graph::node_id * nodes, no;

    if (n <=2 || m <= 2)
        { fprintf(stderr,"error: bad size\n"); exit(0); }
      
    nodes = (Graph::node_id *) malloc(n*sizeof(Graph::node_id)); 
    Graph *BKG = new Graph();
    for (i=0;i<n;i++) {
        no=nodes[i]=BKG->add_node();
        BKG->set_tweights(no,Y[i],0.);
    }
    for (i=0;i<m;i++) {
        dr=lambda*edge_weights[i];
        BKG->add_edge(nodes[edges1[i]-1],nodes[edges2[i]-1],dr,dr);
    }
    BKG->dyadicparametricTV(erreur);   
    for (i=0;i<n;i++) Y[i]=BKG->what_value(nodes[i]);
         delete BKG; free(nodes);
}

//solve min 0.5 \|y-x\|^2 + \lambda \|D x\|_1
//where D is the incidence matrix
void graph_tv(double *Y,// value of nodes
        int n, //number of nodes
        int m, // number of edges
        int* edges1, 
        int* edges2,
//an array of edges of size m. There is an edge edges1[i] -> edges2[i]
        double lambda,
        float erreur)
{
        int i;
    double dval,dr,drd;
    
    dr=lambda;
    
    Graph::node_id * nodes, no;

    if (n <=2 || m <= 2)
        { fprintf(stderr,"error: bad size\n"); exit(0); }
      
    nodes = (Graph::node_id *) malloc(n*sizeof(Graph::node_id)); 
    Graph *BKG = new Graph();
    for (i=0;i<n;i++) {
        no=nodes[i]=BKG->add_node();
        BKG->set_tweights(no,Y[i],0.);
    }
    for (i=0;i<m;i++) {
        BKG->add_edge(nodes[edges1[i]-1],nodes[edges2[i]-1],dr,dr);
    }
    BKG->dyadicparametricTV(erreur);   
    for (i=0;i<n;i++) Y[i]=BKG->what_value(nodes[i]);
         delete BKG; free(nodes);

}