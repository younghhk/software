#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "graph.h"

extern "C" {
  void TV16(double *,int,int,double,float);
}
// TV-minimization 
// this can be linked with C code (tested only on linux) provided the
// flag -lstdc++ is given to the linker
// TV4 = with nearest neighbours interaction
// TV8 = with also next nearest
// TV16 = with 16 neighbours

void TV16(double *G, // solution is returned in the same array
	  int sx,int sy, // size of image
	  double lambda, // weight on the total variation
	  float erreur)  // erreur pour le calcul du min
{

  int i,j,k,sxy=sx*sy;
  Graph::captype dr,drd,dr12;

  Graph::node_id * nodes, no;

  if (sx <=2 || sy <= 2)
    { fprintf(stderr,"error: bad size\n"); exit(0); }

  nodes = (Graph::node_id *) malloc(sxy*sizeof(Graph::node_id));

  lambda *= 0.39269908169872415481/2.; ;
  // pi/16: renormalization to ensure that the TV
  // is "closer" to the isotropic TV

  dr=lambda;
  drd=lambda*0.70710678118654752440;
  dr12=lambda*.447213595499957939281834;

  Graph *BKG = new Graph();

  for (i=0;i<sxy;i++) {
    no=nodes[i]=BKG->add_node();
    BKG->set_tweights(no,G[i],0.);
  }
  k=1;
  for (i=1;i<sx;i++,k++) BKG->add_edge(nodes[k-1],nodes[k],dr,dr);
  BKG->add_edge(nodes[k-sx],nodes[k],dr,dr);
  k++;
  BKG->add_edge(nodes[k-1],nodes[k],dr,dr);
  BKG->add_edge(nodes[k-sx],nodes[k],dr,dr);
  BKG->add_edge(nodes[k-1],nodes[k-sx],drd,drd);
  BKG->add_edge(nodes[k],nodes[k-1-sx],drd,drd);
  for (i=2,k++;i<sx;i++,k++) {
    BKG->add_edge(nodes[k-1],nodes[k],dr,dr);
    BKG->add_edge(nodes[k-sx],nodes[k],dr,dr);
    BKG->add_edge(nodes[k-1],nodes[k-sx],drd,drd);
    BKG->add_edge(nodes[k],nodes[k-1-sx],drd,drd);
    BKG->add_edge(nodes[k],nodes[k-2-sx],dr12,dr12);
    BKG->add_edge(nodes[k-sx],nodes[k-2],dr12,dr12);
  }
  for (j=2;j<sy;j++) {
    BKG->add_edge(nodes[k-sx],nodes[k],dr,dr);
    k++;
    BKG->add_edge(nodes[k-1],nodes[k],dr,dr);
    BKG->add_edge(nodes[k-sx],nodes[k],dr,dr);
    BKG->add_edge(nodes[k-1],nodes[k-sx],drd,drd);
    BKG->add_edge(nodes[k],nodes[k-1-sx],drd,drd);
    BKG->add_edge(nodes[k],nodes[k-2*sx-1],dr12,dr12);
    BKG->add_edge(nodes[k-2*sx],nodes[k-1],dr12,dr12);
    for (i=2,k++;i<sx;i++,k++) {
      BKG->add_edge(nodes[k-1],nodes[k],dr,dr);
      BKG->add_edge(nodes[k-sx],nodes[k],dr,dr);
      BKG->add_edge(nodes[k-1],nodes[k-sx],drd,drd);
      BKG->add_edge(nodes[k],nodes[k-1-sx],drd,drd);      
      BKG->add_edge(nodes[k],nodes[k-2*sx-1],dr12,dr12);
      BKG->add_edge(nodes[k-2*sx],nodes[k-1],dr12,dr12);
      BKG->add_edge(nodes[k],nodes[k-2-sx],dr12,dr12);
      BKG->add_edge(nodes[k-sx],nodes[k-2],dr12,dr12);
    }
  }
  BKG->dyadicparametricTV(erreur);
  //if (BKG->error()) { fprintf(stderr,"error in maxflow\n"); exit(0); }
  
  for (k=0;k<sxy;k++) G[k]=BKG->what_value(nodes[k]);
  delete BKG;
  free(nodes);
}
