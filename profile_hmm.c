#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>


#define EXTERN 
#include "profile_hmm.h"
#undef EXTERN

#define EXTERN extern
#include "utils.h"
#undef EXTERN

double **Xfmatrix;
double **Yfmatrix;
double **Mfmatrix;
double **Xbmatrix;
double **Ybmatrix;
double **Mbmatrix;
#ifdef LONG
double **XLfmatrix;
double **YLfmatrix;
double **XLbmatrix;
double **YLbmatrix;
#endif
#ifdef NCSTATE
double *XNfmatrix;
double *YNfmatrix;
double *XNbmatrix;
double *YNbmatrix;
double *XCfmatrix;
double *YCfmatrix;
double *XCbmatrix;
double *YCbmatrix;
#endif



void init_profile_HMM(int len){
  int i,j;
  double f;
  // Forward Algorithm Probability matrix
  
  Xfmatrix = (double **)(malloc((len+3)*sizeof(double *)));
  Yfmatrix = (double **)(malloc((len+3)*sizeof(double *)));
  Mfmatrix = (double **)(malloc((len+3)*sizeof(double *)));
  Xfmatrix++;
  Yfmatrix++;
  Mfmatrix++;

#ifdef LONG
  XLfmatrix = (double **)(malloc((len+3)*sizeof(double *)));
  YLfmatrix = (double **)(malloc((len+3)*sizeof(double *)));
  XLfmatrix++;
  YLfmatrix++;
#endif
#ifdef NCSTATE
  XNfmatrix = (double *)(malloc((len+3)*sizeof(double)));
  YNfmatrix = (double *)(malloc((len+3)*sizeof(double)));
  XCfmatrix = (double *)(malloc((len+3)*sizeof(double)));
  YCfmatrix = (double *)(malloc((len+3)*sizeof(double)));
  XNfmatrix++;
  YNfmatrix++;
  XCfmatrix++;
  YCfmatrix++;
#endif

  for(i=-1;i<=len+1;i++){
    Xfmatrix[i] = (double *)(malloc((len+3)*sizeof(double)));
    Yfmatrix[i] = (double *)(malloc((len+3)*sizeof(double)));
    Mfmatrix[i] = (double *)(malloc((len+3)*sizeof(double)));
    Xfmatrix[i]++;
    Yfmatrix[i]++;
    Mfmatrix[i]++;
#ifdef LONG
    XLfmatrix[i] = (double *)(malloc((len+3)*sizeof(double)));
    YLfmatrix[i] = (double *)(malloc((len+3)*sizeof(double)));
    XLfmatrix[i]++;
    YLfmatrix[i]++;
#endif
  }
  
  


  
  // Backward Algorithm Probability Matrix

  Xbmatrix = (double **)(malloc((len+3)*sizeof(double *)));
  Ybmatrix = (double **)(malloc((len+3)*sizeof(double *)));
  Mbmatrix = (double **)(malloc((len+3)*sizeof(double *)));
#ifdef LONG
  XLbmatrix = (double **)(malloc((len+3)*sizeof(double *)));
  YLbmatrix = (double **)(malloc((len+3)*sizeof(double *)));
#endif
#ifdef NCSTATE
  XNbmatrix = (double *)(malloc((len+3)*sizeof(double)));
  YNbmatrix = (double *)(malloc((len+3)*sizeof(double)));
  XCbmatrix = (double *)(malloc((len+3)*sizeof(double)));
  YCbmatrix = (double *)(malloc((len+3)*sizeof(double)));
#endif

  for(i=-1;i<=len+1;i++){
    Xbmatrix[i] = (double *)(malloc((len+3)*sizeof(double)));
    Ybmatrix[i] = (double *)(malloc((len+3)*sizeof(double)));
    Mbmatrix[i] = (double *)(malloc((len+3)*sizeof(double)));
#ifdef LONG
    XLbmatrix[i] = (double *)(malloc((len+3)*sizeof(double)));
    YLbmatrix[i] = (double *)(malloc((len+3)*sizeof(double)));
#endif
  }
  // Normalize Parameters
    
  for(i=0;i<N_PEPT;i++){
    for(j=i+1;j<N_PEPT;j++){
      emitPairsDefault[i][j] = emitPairsDefault[j][i];
    }
  }
  
  
  f = 0.0;
  for(i=0;i<N_PEPT;i++){
    for(j=0;j< N_PEPT;j++){
      f = f+emitPairsDefault[i][j];
    }
  }
  
  

  
  
  for(i=0;i<N_PEPT;i++){
    for(j=0;j<=i;j++){
      emitPairsDefault[i][j] /= f;
    }
  }
  
  

  f = 0.0;
  
  
  for(i=0;i<N_PEPT;i++){
    f +=  emitSingleDefault[i];
  }
  
  
  for(i=0;i<N_PEPT;i++){
    emitSingleDefault[i] =  emitSingleDefault[i]/f;
  }
  
  // Take LOGS
#ifndef LONG
  for(i=0;i<3;i++){
    initDistribDefault[i] = log(initDistribDefault[i]);
  }
#else
  for(i=0;i<5;i++){
    initDistribDefault[i] = log(initDistribDefault[i]);
  }
#endif
  

  for(i=0;i<N_PEPT;i++){
    emitSingleDefault[i] = log(emitSingleDefault[i]);
    for(j=0;j<N_PEPT;j++){
      emitPairsDefault[i][j] = log(emitPairsDefault[i][j]);
    }
  }

#ifndef NCSTATE
#ifndef LONG
  MtoXY = log(gapOpenDefault);
  selfM = log(1 - 2*gapOpenDefault);
  selfXY = log(gapExtendDefault);
  XYtoM = log(1 - gapExtendDefault);
#else
  MtoXY = log(gapOpenDefault);
  selfM = log(1 - 2*gapOpenDefault-2*gapOpenDefaultL);
  selfXY = log(gapExtendDefault);
  XYtoM = log(1 - gapExtendDefault);
  MtoXYL = log(gapOpenDefaultL);
  selfXYL = log(gapExtendDefaultL);
  XYLtoM = log(1 - gapExtendDefaultL);
#endif  
#else
#ifndef LONG
  MtoXY = log(gapOpenDefault);
  selfM = log(1 - 2*gapOpenDefault-2*MtoC);
  selfXY = log(gapExtendDefault);
  XYtoM = log(1 - gapExtendDefault);
#else
  MtoXY = log(gapOpenDefault);
  selfM = log(1 - 2*gapOpenDefault-2*gapOpenDefaultL-2*MtoC);
  selfXY = log(gapExtendDefault);
  XYtoM = log(1 - gapExtendDefault);
  MtoXYL = log(gapOpenDefaultL);
  selfXYL = log(gapExtendDefaultL);
  XYLtoM = log(1 - gapExtendDefaultL);
#endif  
  selfN = log(1-NtoM);
  selfC = log(1-NtoM);
  NtoM = log(NtoM);
  MtoC = log(MtoC);
#endif
  
  fprintf(stderr,"Initialized HMM\n");
  fprintf(stderr,"selfM %f MtoXY %f MtoXYL %f\n",selfM,MtoXY,MtoXYL);
  fprintf(stderr,"selfXY %f XYtoM %f\nselfXYL %f XYLtoM %f\n",selfXY,XYtoM,
	  selfXYL,XYLtoM);
  fprintf(stderr,"selfN %f selfC %f NtoM %f MtoC %f\n",selfN,selfC,NtoM,MtoC);
  
}


void init



