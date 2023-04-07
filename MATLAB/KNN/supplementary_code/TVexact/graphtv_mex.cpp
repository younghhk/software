/* -*- mode: c++; c-basic-offset: 3 -*-
 *
 * TV-minimization 
 * this can be linked with C code (tested only on linux) provided the
 * flag -lstdc++ is given to the linker
 * TV4 = with nearest neighbours interaction
 * TV8 = with also next nearest
 * TV16 = with 16 neighbours
 * 
 * Based on the code by 
 *	A. Chambolle and J. Darbon: On total variation
 * 	minimization and surface evolution using parametric maximum flows,
 *	preprint (2008).
 * Their code implements Dorit Hochbaum's algorithm:
 *    D. S. Hochbaum: An efficient algorithm for image segmentation,
 *    Markov random fields and related problems. J. ACM, 48(4):686--701,
 *    2001.     
 * 
 *
 * 
 * 
 * 
 */

/**
 * @authors
 *  Jalal Fadili
 * @date 2008-04-22
 *
 */

/*
 * @file tvmin_mex.cpp
 * @brief Fast TV-minimization
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <iostream>
#include "mex.h"
#include "matrix.h"


//using namespace std ;

extern "C" void graph_tv_weighted(double *,int,int,int*,int*,double*,double,float);
extern "C" void graph_tv(double *,int,int,int*,int*,double,float);




//Need nodes, edges1,edges2, lambda and (edgeweights)

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  double * tmp1, *tmp2, *weights, *tmp3;
  int * edges1,*edges2;
  int m, n,i;
  double lambda;
  int numdeep;

  
  if(nrhs==4)//unweighted version
  {
      tmp1	= mxGetPr(prhs[0]);
      n 	= mxGetM(prhs[0]);
      if (sizeof(int)==sizeof(int32_t)){
          
         edges1	= (int * )mxGetData(prhs[1]); //need to make sure int32 is used.
         edges2	= (int * )mxGetData(prhs[2]);//assume the data is correct
      }
      else
      {
         mexErrMsgTxt("Interger size not right");
      }
      
      m 	= mxGetM(prhs[1]);  
      lambda  = *(mxGetPr(prhs[3]));
      tmp2 = (double*)malloc(n*sizeof(double));
      memcpy(tmp2,tmp1,n*sizeof(double));
      graph_tv(tmp2,n,m,edges1,edges2,lambda,0);
  }
 else if (nrhs==5)
 {//weighted version
      tmp1	= mxGetPr(prhs[0]);
      n 	= mxGetM(prhs[0]); 
      if (sizeof(int)==sizeof(int32_t)){
         edges1	= (int* )mxGetData(prhs[1]); //need to make sure int32 is used.
         edges2	= (int* )mxGetData(prhs[2]);//assume the data is correct
      }
      else
      {
         mexErrMsgTxt("Interger size not right");
      }
      m 	= mxGetM(prhs[1]);  
      lambda  = *(mxGetPr(prhs[3]));
      weights = mxGetPr(prhs[4]);
      tmp3= (double*)malloc(m*sizeof(double));
      memcpy(tmp3,weights,m*sizeof(double));
      tmp2 = (double*)malloc(n*sizeof(double));
      memcpy(tmp2,tmp1,n*sizeof(double));
      for(i=0;i<m;i++)
      graph_tv_weighted(tmp2,n,m,edges1,edges2,tmp3,lambda,0);
      free(tmp3);
    }


  
   plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
   tmp1 = mxGetPr(plhs[0]);
  memcpy(tmp1,tmp2,n*sizeof(double));
  
  free(tmp2);
  return;
}
