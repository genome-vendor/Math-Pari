#include "genpari.h"

#ifdef __cplusplus
char *pariversion = "                    (C++ 68020 version)\n";
#else
char *pariversion = "                      (68020 version)\n";
#endif

uLong overflow,hiremainder;

#ifdef __cplusplus
#define IL inline
#else
#define IL
#endif

IL long
addll(uLong x, uLong y)
{
  uLong z;

  z=x+y;overflow=(z<x)?1:0; return z;
}

IL long
addllx(uLong x, uLong y)
{
  uLong z;

  z=x+y+overflow;overflow=(z<x)||((z==x)&&overflow)?1:0; return z;
}

IL long
subll(uLong x, uLong y)
{
  uLong z;

  z=x-y;overflow=(z>x)?1:0; return z;
}

IL long
subllx(uLong x, uLong y)
{
  uLong z;

  z=x-y-overflow;overflow=(z>x)||((z==x)&&overflow)?1:0; return z;
}

IL long
shiftl(uLong x, uLong y)
{
  hiremainder=x>>(32-y);return (x<<y);
}

IL long
shiftlr(uLong x, uLong y)
{
  hiremainder=x<<(32-y);return (x>>y);
}

int
bfffo(uLong x)
{
  int sc;
  static int tabshi[16]={4,3,2,2,1,1,1,1,0,0,0,0,0,0,0,0};

  if(x&(0xffff0000)) sc=0;else {sc=16;x<<=16;}
  if(!(x&(0xff000000))) {sc+=8;x<<=8;}
  if(x&(0xf0000000)) x>>=28;else {sc+=4;x>>=24;}
  sc+=tabshi[x];return sc;
}

IL long
mulll(uLong x, uLong y)
{
  uLong xlo,xhi,ylo,yhi;
  uLong z;

  xlo=x&65535;xhi=x>>16;ylo=y&65535;yhi=y>>16;
  z=addll(xlo*yhi,xhi*ylo);
  hiremainder=(overflow)?xhi*yhi+65536+(z>>16):xhi*yhi+(z>>16);
  z=addll(xlo*ylo,(z<<16));hiremainder+=overflow;
  return z;
}

IL long
addmul(uLong x, uLong y)
{
  uLong xlo,xhi,ylo,yhi;
  uLong z,z2;

  xlo=x&65535;xhi=x>>16;ylo=y&65535;yhi=y>>16;
  z=addll(xlo*yhi,xhi*ylo);
  z2=(overflow)?xhi*yhi+65536+(z>>16):xhi*yhi+(z>>16);
  z=addll(xlo*ylo,(z<<16));z2+=overflow;
  z=addll(z,hiremainder);hiremainder=z2+overflow;
  return z;
}

long
divll(uLong x, uLong y)
{
#define HIBIT 0x80000000
#define HIMASK 0xffff0000
#define LOMASK 0xffff
#define HIWORD(a) (a >> 16)
/* si le compilateur est bugge, il faut mettre (a >> 16) & LOMASK) */
#define LOWORD(a) (a & LOMASK)
#define GLUE(hi, lo) ((hi << 16) + lo)
#define SPLIT(a, b, c) b = HIWORD(a); c = LOWORD(a)

  uLong v1, v2, u3, u4, q1, q2, aux, aux1, aux2;
  int k;
    
  for(k = 0; !(y & HIBIT); k++)
  {
    hiremainder <<= 1;
    if (x & HIBIT) hiremainder++;
    x <<= 1;
    y <<= 1;
  }
        
  SPLIT(y, v1, v2);
  SPLIT(x, u3, u4);
    
  q1 = hiremainder / v1; if (q1 & HIMASK) q1 = LOMASK;
  hiremainder -= q1 * v1;
  aux = v2 * q1;
  again:
  SPLIT(aux, aux1, aux2);
  if (aux2 > u3) aux1++;
  if (aux1 > hiremainder) {q1--; hiremainder += v1; aux -= v2; goto again;}
  u3 -= aux2;
  hiremainder -= aux1;
  hiremainder <<= 16; hiremainder += u3 & LOMASK;
    
  q2 = hiremainder / v1; if (q2 & HIMASK) q2 = LOMASK;
  hiremainder -= q2 * v1;
  aux = v2 * q2;
  again2:
  SPLIT(aux, aux1, aux2);
  if (aux2 > u4) aux1++;
  if (aux1 > hiremainder) {q2--; hiremainder += v1; aux -= v2; goto again2;}
  u4 -= aux2;
  hiremainder -= aux1;
  hiremainder <<= 16; hiremainder += u4 & LOMASK;
  hiremainder >>= k;
  return GLUE(q1, q2);
}

/* The following program is a C translation of the gerepile
program of mp.s and is intended solely for debugging. In normal
use, it is never called by the 68020 version of PARI/GP */

GEN
gerepilc(GEN l, GEN p, GEN q)
{
  long av,declg,tl;
  GEN ll,pp,l1,l2,l3;

  declg=(long)l-(long)p;if(declg<=0) return q;
  for(ll=l,pp=p;pp>(GEN)avma;) *--ll= *--pp;
  av=(long)ll;
  while((ll<l)||((ll==l)&&(long)q))
  {
    l2=ll+lontyp[tl=typ(ll)];
    if(tl==10) {l3=ll+lgef(ll);ll+=lg(ll);if(l3>ll) l3=l2;}
    else {ll+=lg(ll);l3=ll;} 
    for(;l2<l3;l2++) 
    {
      l1=(GEN)(*l2);
      if((l1<l)&&(l1>=(GEN)avma))
      {
	if(l1<p) *l2+=declg;
	else
	  if(ll<=l) err(gerper);
      }
    }
  }
  if((!((long)q))||((q<p)&&(q>=(GEN)avma)))
  {
    avma=av;return q+(declg>>2);
  }
  else {avma=av;return q;}
}

double
rtodbl(GEN x)
{
  double y,ma;
  long ex,s=signe(x),co=32;
  
  if((!s)||((ex=expo(x))< -1023)) return 0.0;
  if(ex>=0x3ff) err(rtodber);
  ma=(uLong)x[2]+(uLong)x[3]/exp2l(co);
  y=exp2l(ex+1-co);
  return (s<0) ? -y*ma : y*ma;
}

GEN
dbltor(double x)
{
  GEN z;
  long co=32;
  uLong y;
  long ex;
  int n;
  
  if(x==0) {z=cgetr(3);z[2]=0;z[1]=0x800000-308;return z;}
  z=cgetr(4);if(x<0) {setsigne(z,-1);x= -x;} else setsigne(z,1);
  ex=(long)flog2(x);setexpo(z,ex);x=x*exp2l(-ex+co-1);y=(uLong)x;x-=y;z[2]=y;
  y=(long)(exp2l(co)*x);z[3]=y;

  if(!(z[2]&(1<<31)))
  {
    n=bfffo(z[2]);z[2]=shiftl(z[2],n);z[3]=shiftl(z[3],n);
    z[2]+=hiremainder;setexpo(z,expo(z)-n);
  }
  return z;
}

double
gtodouble(GEN x)
{
  GEN     x1;
  static long reel4[4]={0x2010004,0,0,0};
  
  if (typ(x)==2) x1=x;else gaffect(x,x1=(GEN)reel4);
  return rtodbl(x1);
}

#ifdef macintosh

#include <Events.h>
     
long
timer(void)
{
  static long oldticks;
  long ticks = TickCount();
  long delay = ticks - oldticks;
  oldticks = ticks;
  return 50 * delay / 3;
}

long
timer2(void)
{
  static long oldticks;
  long ticks = TickCount();
  long delay = ticks - oldticks;
  oldticks = ticks;
  return 50 * delay / 3;
}

#else

#include <sys/time.h>
#include <sys/resource.h>

long
timer(void)
{
  static long oldmusec;
  static long oldsec;
  long delay;
  struct rusage r;
  struct timeval t;
  getrusage(0,&r);t=r.ru_utime;
  delay = 1000 * (t.tv_sec - oldsec) + (t.tv_usec - oldmusec) / 1000;
  oldmusec = t.tv_usec;
  oldsec = t.tv_sec;
  return delay;
}

long
timer2(void)
{
  static long oldmusec;
  static long oldsec;
  long delay;
  struct rusage r;
  struct timeval t;
  getrusage(0,&r);t=r.ru_utime;
  delay = 1000 * (t.tv_sec - oldsec) + (t.tv_usec - oldmusec) / 1000;
  oldmusec = t.tv_usec;
  oldsec = t.tv_sec;
  return delay;
}

#endif



