#include "mpdefs.h"

/* mp.c ou mp.s */

GEN     gerepile(long l, long p, GEN q);
GEN     mptrunc(GEN x),mpent(GEN x),shifts(long x, long y),shifti(GEN x, long n);
GEN     addsi(long x, GEN y),addsr(long x, GEN y),addii(GEN x, GEN y),addir(GEN x, GEN y),addrr(GEN x, GEN y), addss(long x, long y);
GEN     mulss(long x, long y),mulsi(long x, GEN y),mulsr(long x, GEN y),mulii(GEN x, GEN y),mulir(GEN x, GEN y),mulrr(GEN x, GEN y);
GEN     divsi(long x, GEN y),divis(GEN y, long x),divsr(long x, GEN y),divrs(GEN x, long y),divir(GEN x, GEN y);
GEN     divri(GEN x, GEN y),divrr(GEN x, GEN y),convi(GEN x),confrac(GEN x);
GEN     modss(long x, long y),resss(long x, long y),modsi(long x, GEN y),modii(GEN x, GEN y);
GEN     dvmdii(GEN x, GEN y, GEN *z);
long    vals(long x);
int     cmpss(long x, long y),cmpsi(long x, GEN y),cmpii(GEN x, GEN y),cmprr(GEN x, GEN y);
void    affir(GEN x, GEN y),affrr(GEN x, GEN y);
void    diviiz(GEN x, GEN y, GEN z),cgiv(GEN x);
void    mpdivz(GEN x, GEN y, GEN z);
void    modiiz(GEN x, GEN y, GEN z);

#ifdef __cplusplus
  #include "INLINE.h"
#endif /* defined(__cplusplus) */

#include "mpinline.h"
