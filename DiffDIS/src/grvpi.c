#include "../include/grvpi.h"
#include <math.h>
#include <string.h>

  void grvpiho_(double* ZX, double* ZQ, double* ZUV, double* ZDV, double* ZUDB,
              double* ZSB, double* ZCB, double* ZBB, double* ZGL);
  void GRVpi_HO(double x, double QQ, double xF[]) {
    //--- xF[q] = x*(q + \bar q)
    double ZX,ZQ,ZUV,ZDV,SEA,ZSB,ZCB,ZBB,ZGL;
    if(x > 0.99999) {
      memset(xF,0,6*sizeof(double));
      return;
    }
    ZX = x;
    ZQ = sqrt(QQ);
    grvpiho_(&ZX,&ZQ,&ZUV,&ZDV,&SEA,&ZSB,&ZCB,&ZBB,&ZGL);
    xF[0] = ZGL;
    //xF[1] = ZDV+SEA;
    xF[1] = ZDV;
    xF[2] = SEA;
    xF[3] = ZSB;
    xF[4] = ZCB;
    xF[5] = ZBB;
  }

  void grvpilo_(double* ZX, double* ZQ, double* ZUV, double* ZDV, double* ZUDB,
              double* ZSB, double* ZCB, double* ZBB, double* ZGL);
  void GRVpi_LO(double x, double QQ, double xF[]) {
    double ZX,ZQ,ZUV,ZDV,SEA,ZSB,ZCB,ZBB,ZGL;
    if(x > 0.99999) {
      memset(xF,0,6*sizeof(double));
      return;
    }
    ZX = x;
    ZQ = sqrt(QQ);
    grvpilo_(&ZX,&ZQ,&ZUV,&ZDV,&SEA,&ZSB,&ZCB,&ZBB,&ZGL);
    xF[0] = ZGL;
    //xF[1] = ZDV+SEA;
    xF[1] = ZDV;
    xF[2] = SEA;
    xF[3] = ZSB;
    xF[4] = ZCB;
    xF[5] = ZBB;
  }
