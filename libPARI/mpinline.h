/* This file contains declarations/definitions of functions that _can_
   be inlined. */

#ifdef INCLUDE_INLINE

INLINE GEN
cgeti(long x)
{
  uLong p1;
  GEN z;
  
  p1=avma-(((uLong)x)<<TWOPOTBYTES_IN_LONG);
  if(p1<bot) err(errpile);
  avma=p1;z=(GEN)p1;z[0]=evaltyp(1)+evalpere(1)+evallg(x);
  return z;
}

INLINE GEN
cgetr(long x)
{
  uLong p1;
  GEN z;
  
  p1=avma-(((uLong)x)<<TWOPOTBYTES_IN_LONG);
  if(p1<bot) err(errpile);
  avma=p1;z=(GEN)p1;z[0]=evaltyp(2)+evalpere(1)+evallg(x);
  return z;
}

INLINE GEN
icopy(GEN x)
{
  GEN y;
  long lx=lgef(x),i;
  
  y=cgeti(lx);
  for(i=1;i<lx;i++) y[i]=x[i];
  return y;
}

INLINE GEN
rcopy(GEN x)
{
  GEN y;
  long lx=lg(x),i;
  
  y=cgetr(lx);
  for(i=1;i<lx;i++) y[i]=x[i];
  return y;
}

INLINE GEN
cgetg(long x, long y)
{
  uLong p1;
  GEN z;
  
  p1=avma-(((uLong)x)<<TWOPOTBYTES_IN_LONG);
  if(p1<bot) err(errpile);
  avma=p1;z=(GEN)p1;z[0]=evalpere(1)+evaltyp(y)+evallg(x);
  return z;
}

INLINE GEN
negi(GEN x)
{
  long s=signe(x);
  GEN y;
  
  if(!s) return gzero;
  y=icopy(x);setsigne(y,-s);
  return y;
}

INLINE GEN
negr(GEN x)
{
  GEN y;
  
  y=rcopy(x);setsigne(y,-signe(x));
  return y;
}


INLINE GEN
absi(GEN x)
{
  GEN y;
  long s=signe(x);
  
  if(!s) return gzero;
  y=icopy(x);setsigne(y,1);return y;
}

INLINE GEN
absr(GEN x)
{
  GEN y;
  long s=signe(x);
  
  y=rcopy(x);
  if(s) setsigne(y,1);
  return y;
}

INLINE int
expi(GEN x)
{
  long lx=lgef(x);
  
  return lx==2 ? -HIGHEXPOBIT : ((lx-2)<<TWOPOTBITS_IN_LONG)-bfffo(x[2])-1;
}

INLINE GEN
stoi(long x)
{
  GEN y;
  
  if(!x) return gzero;
  y=cgeti(3);
  if(x>0) {y[1]=evalsigne(1)+evallgef(3);y[2]=x;}
  else {y[1]=evalsigne(-1)+evallgef(3);y[2]= -x;}
  return y;
}

INLINE long
itos(GEN x)
{
  long s=signe(x),p2;
  uLong p1;
  
  if(!s) return 0;
  if(lgef(x)>3) err(affer2);
  p1=x[2];if(p1>=(uLong)HIGHBIT) err(affer2);
  p2=(s>0)?p1:(-((long)p1));return p2;
}

INLINE void
affii(GEN x, GEN y)
{
  long lx=lgef(x),i;
  
  if(x==y) return;
  if(lg(y)<lx) err(affer3);
  for(i=1;i<lx;i++) y[i]=x[i];
}

INLINE void
mpaff(GEN x, GEN y)
{
  long tx=typ(x),ty=typ(y);
  if(tx==1)
  {
    if(ty==1) affii(x,y);else affir(x,y);
  }
  else
  {
    if(ty==1) affri(x,y);else affrr(x,y);
  }
}

INLINE void
affsi(long s, GEN x)
{
  long lx;
  
  if(!s) {x[1]=2;return;}
  lx=lg(x);if(lx<3) err(affer1);
  if(s>0) {x[1]=evalsigne(1)+evallgef(3);x[2]=s;}
  else {x[1]=evalsigne(-1)+evallgef(3);x[2]= -s;}
}

INLINE void
affsr(long s, GEN x)
{
  long l,i,d;
  
  if(!s)
  {
    l=(lg(x)-2)<<TWOPOTBITS_IN_LONG;x[1]=HIGHEXPOBIT-l;x[2]=0;
  }
  else
  {
    d=1;if(s<0) {d= -1;s= -s;}
    l=bfffo(s);x[1]=evalexpo((BITS_IN_LONG-1)-l)+evalsigne(d);
    x[2]=(s<<l);for(i=3;i<lg(x);i++) x[i]=0;
  }
}

INLINE GEN
shiftr(GEN x, long n)
{
  long l;
  GEN y;
  
  y=rcopy(x);l=expo(x)+n;
  if(l>=HIGHEXPOBIT||l<-HIGHEXPOBIT) err(shier2);
  setexpo(y,l);return y;
}

INLINE int
cmpir(GEN x, GEN y)
{
  long av=avma;
  int p;
  GEN z;
  
  if(!signe(x)) return -signe(y);
  z=cgetr(lg(y));affir(x,z);
  p=cmprr(z,y);avma=av;return p;
}

INLINE int
mpcmp(GEN x, GEN y)
{
  if(typ(x)==1) return (typ(y)==1) ? cmpii(x,y) : cmpir(x,y);
  return (typ(y)==1) ? -cmpir(y,x) : cmprr(x,y);
}

INLINE int
cmpsr(long x, GEN y)
{
  int p;
  long av;
  GEN z;
  
  if(!x) return -signe(y);
  av=avma;z=cgetr(3);affsr(x,z);
  p=cmprr(z,y);avma=av;return p;
}	


INLINE GEN
mpadd(GEN x, GEN y)
{
  if(typ(x)==1) return (typ(y)==1) ? addii(x,y) : addir(x,y);
  return (typ(y)==1) ? addir(y,x) : addrr(x,y);
}

INLINE void
addssz(long x, long y, GEN z)
{
  long av=avma;
  GEN p1;
  
  if(typ(z)==1) gops2ssz(addss,x,y,z);
  else
  {
    p1=cgetr(lg(z));affsr(x,p1);p1=addrs(p1,y);
    affrr(p1,z);avma=av;
  }
}

INLINE GEN
subii(GEN x, GEN y)
{
  long s=signe(y);
  GEN z;
  
  if(x==y) return gzero;
  setsigne(y,-s);z=addii(x,y);setsigne(y,s);
  return z;
}

INLINE GEN
subrr(GEN x, GEN y)
{
  long s=signe(y);
  GEN z;
  
  if(x==y)
  {
    z=cgetr(3);z[2]=0;z[1]=HIGHEXPOBIT-(lg(x)<<TWOPOTBITS_IN_LONG);return z;
  }
  setsigne(y,-s);z=addrr(x,y);setsigne(y,s);return z;
}

INLINE GEN
subir(GEN x, GEN y)
{
  long s=signe(y);
  GEN z;
  
  setsigne(y,-s);z=addir(x,y);setsigne(y,s);return z;
}

INLINE GEN
subri(GEN x, GEN y)
{
  long s=signe(y);
  GEN z;
  
  setsigne(y,-s);z=addir(y,x);setsigne(y,s);return z;
}

INLINE GEN
mpsub(GEN x, GEN y)
{
  if(typ(x)==1) return (typ(y)==1) ? subii(x,y) : subir(x,y);
  return (typ(y)==1) ? subri(x,y) : subrr(x,y);
}

INLINE GEN
subsi(long x, GEN y)
{
  long s=signe(y);
  GEN z;
  
  setsigne(y,-s);z=addsi(x,y);setsigne(y,s);return z;
}

INLINE GEN
subsr(long x, GEN y)
{
  long s=signe(y);
  GEN z;
  
  setsigne(y,-s);z=addsr(x,y);setsigne(y,s);return z;
}

INLINE GEN
subss(long x, long y) {return addss(-y,x);}

INLINE void
subssz(long x, long y, GEN z)
{
  long av=avma;
  GEN p1;
  
  if(typ(z)==1) gops2ssz(addss,x,-y,z);
  else
  {
    p1=cgetr(lg(z));affsr(x,p1);p1=addrs(p1,-y);
    affrr(p1,z);avma=av;
  }
}

INLINE GEN
mpmul(GEN x, GEN y)
{
  if(typ(x)==1) return (typ(y)==1) ? mulii(x,y) : mulir(x,y);
  return (typ(y)==1) ? mulir(y,x) : mulrr(x,y);
}

INLINE void
mulssz(long x, long y, GEN z)
{
  long av=avma;
  GEN p1;
  
  if(typ(z)==1) gops2ssz(mulss,x,y,z);
  else
  {
    p1=cgetr(lg(z));affsr(x,p1);p1=mulsr(y,p1);
    mpaff(p1,z);avma=av;
  }
}

INLINE void
mulsii(long x, GEN y, GEN z)
{
  long av=avma;
  GEN p1;
  
  p1=mulsi(x,y);affii(p1,z);avma=av;
}

INLINE void
addsii(long x, GEN y, GEN z)
{
  long av=avma;
  GEN p1;
  
  p1=addsi(x,y);affii(p1,z);avma=av;
}

INLINE long
divisii(GEN x, long y, GEN z)
{
  long av=avma,k;
  GEN p1;
  
  p1=divis(x,y);affii(p1,z);avma=av;
  k=hiremainder;return k;
}

INLINE long
vali(GEN x)
{
  long i,lx=lgef(x);
  
  if(!signe(x)) return -1;
  for(i=lx-1;(i>=2)&&(!x[i]);i--);
  return ((lx-1-i)<<TWOPOTBITS_IN_LONG)+vals(x[i]);
}

INLINE GEN
divss(long x, long y)
{
  long p1;
  
  if(!y) err(diver1);
  hiremainder=0;p1=divll((uLong)labs(x),(uLong)labs(y));
  if(y<0) {hiremainder= -((long)hiremainder);p1= -p1;}
  if(x<0) p1= -p1;
  return stoi(p1);
}

INLINE GEN
mpdiv(GEN x, GEN y)
{
  if(typ(x)==1) return (typ(y)==1) ? divii(x,y) : divir(x,y);
  return (typ(y)==1) ? divri(x,y) : divrr(x,y);
}

INLINE GEN
dvmdss(long x, long y, GEN *z)
{
  GEN p1;

  p1=divss(x,y);*z=stoi(hiremainder);return p1;
}

INLINE void
dvmdssz(long x, long y, GEN z, GEN t)
{
  long av=avma;
  GEN p1;

  p1=divss(x,y);affsi(hiremainder,t);mpaff(p1,z);avma=av;
}

INLINE GEN
dvmdsi(long x, GEN y, GEN *z)
{
  GEN p1;
  p1=divsi(x,y);*z=stoi(hiremainder);return p1;
}

INLINE void
dvmdsiz(long x, GEN y, GEN z, GEN t)
{
  long av=avma;
  GEN p1;
  
  p1=divsi(x,y);affsi(hiremainder,t);mpaff(p1,z);avma=av;
}

INLINE GEN
dvmdis(GEN x, long y, GEN *z)
{
  GEN p1;
  p1=divis(x,y);*z=stoi(hiremainder);
  return p1;
}

INLINE void
dvmdisz(GEN x, long y, GEN z, GEN t)
{
  long av=avma;
  GEN p1;
  
  p1=divis(x,y);affsi(hiremainder,t);mpaff(p1,z);avma=av;
}

INLINE void
dvmdiiz(GEN x, GEN y, GEN z, GEN t)
{
  long av=avma;
  GEN p1,p2;

  p1=dvmdii(x,y,&p2);mpaff(p1,z);mpaff(p2,t);avma=av;
}

INLINE GEN
ressi(long x, GEN y)
{
  divsi(x,y);return stoi(hiremainder);
}

INLINE GEN
modis(GEN x, long y)
{
  divis(x,y);if(!hiremainder) return gzero;
  return (signe(x)>0) ? stoi(hiremainder) : stoi(labs(y)+hiremainder);
}

INLINE GEN
resis(GEN x, long y)
{
  divis(x,y);return stoi(hiremainder);
}
     
INLINE void
divisz(GEN x, long y, GEN z)
{
  long av=avma;
  GEN p1;
  
  if(typ(z)==1) gops2gsz(divis,x,y,z);
  else {p1=cgetr(lg(z));affir(x,p1);p1=divrs(p1,y);affrr(p1,z);avma=av;}
}

INLINE void
divsiz(long x, GEN y, GEN z)
{
  long av=avma,lz;
  GEN p1,p2;
  
  if(typ(z)==1) gops2sgz(divsi,x,y,z);
  else
  {
    lz=lg(z);p1=cgetr(lz);p2=cgetr(lz);affsr(x,p1);affir(y,p2);
    p1=divrr(p1,p2);affrr(p1,z);avma=av;
  }
}

INLINE void
divssz(long x, long y, GEN z)
{
  long av=avma;
  GEN p1;
  
  if(typ(z)==1) gops2ssz(divss,x,y,z);
  else {p1=cgetr(lg(z));affsr(x,p1);p1=divrs(p1,y);affrr(p1,z);avma=av;}
}

INLINE void
divrrz(GEN x, GEN y, GEN z)
{
  long av=avma;
  GEN p1;

  p1=divrr(x,y);mpaff(p1,z);avma=av;
}

INLINE void
resiiz(GEN x, GEN y, GEN z)
{
  long av=avma;
  GEN p1;

  p1=resii(x,y);affii(p1,z);avma=av;
}

INLINE int
mpdivis(GEN x, GEN y, GEN z)
{
  long av=avma;
  GEN p1,p2;

  p1=dvmdii(x,y,&p2);
  if(signe(p2)) {avma=av;return 0;}
  affii(p1,z);avma=av;return 1;
}

INLINE int
divise(GEN x, GEN y)
{
  long av=avma;
  GEN p1;

  p1=dvmdii(x,y,(GEN *)-1);avma=av;return signe(p1) ? 0 : 1;
}

INLINE double
gtodouble(GEN x)
{
  GEN x1;
  long t=typ(x);
  static long reel4[4]={evaltyp(2)+evalpere(1)+evallg(4),0,0,0};
  
  if (t==2) x1=x;
  else gaffect(x,x1=(GEN)reel4);
  return rtodbl(x1);
}

#ifdef __HAS_NO_ASM__

INLINE long
addll(uLong x, uLong y)
{
  uLong z;

  z=x+y;overflow=(z<x)?1:0; return (long)z;
}

INLINE long
addllx(uLong x, uLong y)
{
  uLong z;

  z=x+y+overflow;overflow=(z<x)||((z==x)&&overflow)?1:0; return (long)z;
}

INLINE long
subll(uLong x, uLong y)
{
  uLong z;

  z=x-y;overflow=(z>x)?1:0; return (long)z;
}

INLINE long
subllx(uLong x, uLong y)
{
  uLong z;

  z=x-y-overflow;overflow=(z>x)||((z==x)&&overflow)?1:0; return z;
}

INLINE long
shiftl(uLong x, uLong y)
{
  hiremainder=x>>(BITS_IN_LONG-y);return (x<<y);
}

INLINE long
shiftlr(uLong x, uLong y)
{
  hiremainder=x<<(BITS_IN_LONG-y);return (x>>y);
}

/* Version Peter Montgomery */

INLINE long
mulll(uLong x, uLong y)
{
  uLong xlo,xhi,ylo,yhi;
  uLong xylo,xymid,xyhi;
  uLong xymidhi, xymidlo;
/*
        Assume (for presentation) that BITS_IN_LONG = 32.
        Then 0 <= xhi, xlo, yhi, ylo <= 2^16 - 1.  Hence

-2^31  + 2^16 <= (xhi - 2^15)*(ylo - 2^15) + (xlo - 2^15)*(yhi - 2^15) <= 2^31.

        If xhi*ylo + xlo*yhi = 2^32*overflow + xymid, then

-2^32 + 2^16 <= 2^32*overflow + xymid - 2^15*(xhi + ylo + xlo + yhi) <= 0.

2^16*overflow <= (xhi+xlo+yhi+ylo)/2 - xymid/2^16 <= 2^16*overflow + 2^16-1

        This inequality was derived using exact (rational) arithmetic;
        it remains valid when we truncate the two middle terms.
*/

  xlo = LOWWORD(x); xhi = HIGHWORD(x);
  ylo = LOWWORD(y); yhi = HIGHWORD(y);

  xylo = xlo*ylo; xyhi = xhi*yhi;
  xymid = (xhi + xlo)*(yhi + ylo) - (xyhi + xylo);

  xymidhi = HIGHWORD(xymid);
  xymidlo = xymid << BITS_IN_HALFULONG;

  xylo += xymidlo;
  hiremainder = xyhi + xymidhi + (xylo < xymidlo)
     + (((((xhi + xlo) + (yhi + ylo)) >> 1) - xymidhi) & HIGHMASK);

  return xylo;
}

INLINE long
addmul(uLong x, uLong y)
{
  uLong xlo,xhi,ylo,yhi;
  uLong xylo,xymid,xyhi;
  uLong xymidhi, xymidlo;

  xlo = LOWWORD(x); xhi = HIGHWORD(x);
  ylo = LOWWORD(y); yhi = HIGHWORD(y);

  xylo = xlo*ylo; xyhi = xhi*yhi;
  xymid = (xhi + xlo)*(yhi + ylo) - (xyhi + xylo);

  xylo += hiremainder; xyhi += (xylo < hiremainder);

  xymidhi = HIGHWORD(xymid);
  xymidlo = xymid << BITS_IN_HALFULONG;

  xylo += xymidlo;
  hiremainder = xyhi + xymidhi + (xylo < xymidlo)
     + (((((xhi + xlo) + (yhi + ylo)) >> 1) - xymidhi) & HIGHMASK);

  return xylo;
}

#endif

#endif /* defined(INCLUDE_INLINE) */

#ifndef __cplusplus

GEN cgeti(long x),cgetr(long x),stoi(long x);
GEN cgetg(long x, long y),negi(GEN x),negr(GEN x),absi(GEN x),absr(GEN x);
GEN shiftr(GEN x, long n),mpadd(GEN x, GEN y),subii(GEN x, GEN y);
GEN subrr(GEN x, GEN y),subir(GEN x, GEN y),subri(GEN x, GEN y);
GEN mpsub(GEN x, GEN y),subsi(long x, GEN y),subsr(long x, GEN y);
GEN subss(long x, long y),mpmul(GEN x, GEN y);
GEN divss(long x, long y),mpdiv(GEN x, GEN y),dvmdss(long x, long y, GEN *z);
GEN ressi(long x, GEN y),modis(GEN x, long y),resis(GEN x, long y);
double  gtodouble(GEN x);
long itos(GEN x),divisii(GEN x, long y, GEN z),vali(GEN x);
int expi(GEN x),cmpir(GEN x, GEN y),mpcmp(GEN x, GEN y);
int cmpsr(long x, GEN y),mpdivis(GEN x, GEN y, GEN z),divise(GEN x, GEN y);
void affii(GEN x, GEN y),mpaff(GEN x, GEN y),affsi(long s, GEN x);
void affsr(long s, GEN x),addssz(long x, long y, GEN z);
void mulssz(long x, long y, GEN z),mulsii(long x, GEN y, GEN z);
void addsii(long x, GEN y, GEN z),dvmdssz(long x, long y, GEN z, GEN t);
void dvmdsiz(long x, GEN y, GEN z, GEN t),dvmdisz(GEN x, long y, GEN z, GEN t);
void dvmdiiz(GEN x, GEN y, GEN z, GEN t),divisz(GEN x, long y, GEN z);
void divsiz(long x, GEN y, GEN z),divssz(long x, long y, GEN z);
void divrrz(GEN x, GEN y, GEN z),resiiz(GEN x, GEN y, GEN z);

GEN     icopy(GEN x), rcopy(GEN x);
GEN     dvmdsi(long x, GEN y, GEN *z),dvmdis(GEN x, long y, GEN *z);


#endif  /* !defined __cplusplus */

/* We should use the following definitions if there is assembler code
 * or it is not c++ (since for c these functions are extern anyway). */

#if !defined(__cplusplus) || !defined(__HAS_NO_ASM__)

#ifdef __cplusplus
extern "C" {
#endif

extern long addll(uLong x, uLong y);
extern long addllx(uLong x, uLong y);
extern long subll(uLong x, uLong y);
extern long subllx(uLong x, uLong y);
extern long shiftl(uLong x, uLong y);
extern long shiftlr(uLong x, uLong y);
extern long mulll(uLong x, uLong y);
extern long addmul(uLong x, uLong y);
extern long divll(uLong x, uLong y);
extern int bfffo(uLong x);

#ifdef __cplusplus
}
#endif

#endif /* !defined(__cplusplus) || !defined(__HAS_NO_ASM__) */
