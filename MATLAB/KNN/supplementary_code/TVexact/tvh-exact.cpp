#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "graph.h"

extern "C" {
  void TV4(double *,int,int,double,float);
  void TV8(double *,int,int,double,float);
}
// TV-minimization 
// this can be linked with C code (tested only on linux) provided the
// flag -lstdc++ is given to the linker
// TV4 = with nearest neighbours interaction
// TV8 = with also next nearest
// TV16 = with 16 neighbours


void TV4(double *G, // solution is returned in the same array
	 int sx,int sy,// size of image
	 double lambda, // weight on the total variation
	 float erreur)  // erreur pour le calcul du min

{

  int i,j,k,sxy=sx*sy;
  double dval, dr,drd;

  Graph::node_id * nodes, no;

  if (sx <=2 || sy <= 2)
    { fprintf(stderr,"error: bad size\n"); exit(0); }

  nodes = (Graph::node_id *) malloc(sxy*sizeof(Graph::node_id));

  lambda *= 0.78539816339744830962 ; // pi/4
  // renormalization to ensure that the TV
  // is "as close as possible" to the isotropic TV [TV4 is anyway very far]

  dr=lambda;

  Graph *BKG = new Graph();

  for (i=0;i<sxy;i++) {
    no=nodes[i]=BKG->add_node();
    BKG->set_tweights(no,G[i],0.);
  }
  
  for (i=1;i<sx;i++) BKG->add_edge(nodes[i-1],nodes[i],dr,dr);
  for (j=1,k=sx;j<sy;j++) {
    BKG->add_edge(nodes[k-sx],nodes[k],dr,dr);
    for (i=1,k++;i<sx;i++,k++) {
      BKG->add_edge(nodes[k-1],nodes[k],dr,dr);
      BKG->add_edge(nodes[k-sx],nodes[k],dr,dr);
    }
  }
  BKG->dyadicparametricTV(erreur);
  //if (BKG->error()) { fprintf(stderr,"error in maxflow\n"); exit(0); }
  
  for (k=0;k<sxy;k++) G[k]=BKG->what_value(nodes[k]);
  delete BKG; free(nodes);
}

void TV8(double *G, // solution is returned in the same array
	 int sx,int sy, // size of image
	 double lambda, // weight on the total variation
	 float erreur)  // erreur pour le calcul du min

{

  int i,j,k,sxy=sx*sy;
  double dval, dr,drd;

  Graph::node_id * nodes, no;

  if (sx <=2 || sy <= 2)
    { fprintf(stderr,"error: bad size\n"); exit(0); }
  // levels ranges from 0 to numlevel (numlevel+1 values)

  nodes = (Graph::node_id *) malloc(sxy*sizeof(Graph::node_id));

  lambda *= 0.39269908169872415481 ; // pi/8
  // pi/8: renormalization to ensure that the TV
  // is "as close as possible" to the isotropic TV

  dr=lambda;
  drd=lambda*0.70710678118654752440;

  Graph *BKG = new Graph();

  for (i=0;i<sxy;i++) {
    no=nodes[i]=BKG->add_node();
    BKG->set_tweights(no,G[i],0.);
  }
  for (i=1;i<sx;i++) BKG->add_edge(nodes[i-1],nodes[i],dr,dr);
  for (j=1,k=sx;j<sy;j++) {
    BKG->add_edge(nodes[k-sx],nodes[k],dr,dr);
    for (i=1,k++;i<sx;i++,k++) {
      BKG->add_edge(nodes[k-1],nodes[k],dr,dr);
      BKG->add_edge(nodes[k-sx],nodes[k],dr,dr);
      BKG->add_edge(nodes[k-1],nodes[k-sx],drd,drd);
      BKG->add_edge(nodes[k-sx-1],nodes[k],drd,drd);
    }
  }
  BKG->dyadicparametricTV(erreur);
  //if (BKG->error()) { fprintf(stderr,"error in maxflow\n"); exit(0); }
  
  for (k=0;k<sxy;k++) G[k]=BKG->what_value(nodes[k]);
  delete BKG; free(nodes);
}
