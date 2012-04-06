#ifndef LINPACK_H_
#define LINPACK_H_

#include <Blas.h>

extern void DCHDC(double*, int, int, double*, int*, int, int*);
extern void DCHDD(double*, int, int, double*, double*, int, int, double*, double*, double*, double*, int*);
extern void DCHEX(double*, int, int, int, int, double*, int, int, double*, double*, int);
extern void DCHUD(double*, int, int, double*, double*, int, int, double*, double*, double*, double*);
extern void DGBCO(double*, int, int, int, int, int*, double*, double*);
extern void DGBDI(double*, int, int, int, int, int*, double*);
extern void DGBFA(double*, int, int, int, int, int*, int*);
extern void DGBSL(double*, int, int, int, int, int*, double*, int);
extern void DGECO(double*, int, int, int*, double*, double*);
extern void DGEDI(double*, int, int, int*, double*, double*, int);
extern void DGEFA(double*, int, int, int*, int*);
extern void DGESL(double*, int, int, int*, double*, int);
extern void DGTSL(int, double*, double*, double*, double*, int*);
extern void DPBCO(double*, int, int, int, double*, double*, int*);
extern void DPBDI(double*, int, int, int, double*);
extern void DPBFA(double*, int, int, int, int*);
extern void DPBSL(double*, int, int, int, double*);
extern void DPOCO(double*, int, int, double*, double*, int*);
extern void DPODI(double*, int, int, double*, int);
extern void DPOFA(double*, int, int, int*);
extern void DPOSL(double*, int, int, double*);
extern void DPPCO(double*, int, double*, double*, int*);
extern void DPPDI(double*, int, double*, int);
extern void DPPFA(double*, int, int*);
extern void DPPSL(double*, int, double*);
extern void DPTSL(int, double*, double*, double*);
extern void DQRDC(double*, int, int, int, double*, int*, double*, int);
extern void DQRSL(double*, int, int, int, double*, double*, double*, double*, double*, double*, double*, int, int*);
extern void DSICO(double*, int, int, int*, double*, double*);
extern void DSIDI(double*, int, int, int*, double*, int*, double*, int);
extern void DSIFA(double*, int, int, int*, int*);
extern void DSISL(double*, int, int, int*, double*);
extern void DSPCO(double*, int, int*, double*, double*);
extern void DSPDI(double*, int, int*, double*, int*, double*, int);
extern void DSPFA(double*, int, int*, int*);
extern void DSPSL(double*, int, int*, double*);
extern void DSVDC(double*, int, int, int, double*, double*, double*, int, double*, int, double*, int, int*);
extern void DTRCO(double*, int, int, double*, double*, int);
extern void DTRDI(double*, int, int, double*, int, int*);
extern void DTRSL(double*, int, int, double*, int, int*);

extern int dchdc_(double*, int*, int*, double*, int*, int*, int*);
extern int dchdd_(double*, int*, int*, double*, double*, int*, int*, double*, double*, double*, double*, int*);
extern int dchex_(double*, int*, int*, int*, int*, double*, int*, int*, double*, double*, int*);
extern int dchud_(double*, int*, int*, double*, double*, int*, int*, double*, double*, double*, double*);
extern int dgbco_(double*, int*, int*, int*, int*, int*, double*, double*);
extern int dgbdi_(double*, int*, int*, int*, int*, int*, double*);
extern int dgbfa_(double*, int*, int*, int*, int*, int*, int*);
extern int dgbsl_(double*, int*, int*, int*, int*, int*, double*, int*);
extern int dgeco_(double*, int*, int*, int*, double*, double*);
extern int dgedi_(double*, int*, int*, int*, double*, double*, int*);
extern int dgefa_(double*, int*, int*, int*, int*);
extern int dgesl_(double*, int*, int*, int*, double*, int*);
extern int dgtsl_(int*, double*, double*, double*, double*, int*);
extern int dpbco_(double*, int*, int*, int*, double*, double*, int*);
extern int dpbdi_(double*, int*, int*, int*, double*);
extern int dpbfa_(double*, int*, int*, int*, int*);
extern int dpbsl_(double*, int*, int*, int*, double*);
extern int dpoco_(double*, int*, int*, double*, double*, int*);
extern int dpodi_(double*, int*, int*, double*, int*);
extern int dpofa_(double*, int*, int*, int*);
extern int dposl_(double*, int*, int*, double*);
extern int dppco_(double*, int*, double*, double*, int*);
extern int dppdi_(double*, int*, double*, int*);
extern int dppfa_(double*, int*, int*);
extern int dppsl_(double*, int*, double*);
extern int dptsl_(int*, double*, double*, double*);
extern int dqrdc_(double*, int*, int*, int*, double*, int*, double*, int*);
extern int dqrsl_(double*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, int*, int*);
extern int dsico_(double*, int*, int*, int*, double*, double*);
extern int dsidi_(double*, int*, int*, int*, double*, int*, double*, int*);
extern int dsifa_(double*, int*, int*, int*, int*);
extern int dsisl_(double*, int*, int*, int*, double*);
extern int dspco_(double*, int*, int*, double*, double*);
extern int dspdi_(double*, int*, int*, double*, int*, double*, int*);
extern int dspfa_(double*, int*, int*, int*);
extern int dspsl_(double*, int*, int*, double*);
extern int dsvdc_(double*, int*, int*, int*, double*, double*, double*, int*, double*, int*, double*, int*, int*);
extern int dtrco_(double*, int*, int*, double*, double*, int*);
extern int dtrdi_(double*, int*, int*, double*, int*, int*);
extern int dtrsl_(double*, int*, int*, double*, int*, int*);

#endif
