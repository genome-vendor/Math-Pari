/********************************************************************/
/********************************************************************/
/**                                                                **/
/**                  BIBLIOTHEQUE  MATHEMATIQUE                    **/
/**                     (deuxieme partie)                          **/
/**                                                                **/
/**                     Copyright Babe Cool                        **/
/**                                                                **/
/********************************************************************/
/********************************************************************/

#include "genpari.h"
#include "rect.h"

static GEN rectmstring(long ne, char *s, int l);

/********************************************************************/
/********************************************************************/
/**                                                                **/
/**                        ITERATIONS                              **/
/**                                                                **/
/********************************************************************/
/********************************************************************/

GEN
forpari(entree *ep, GEN a, GEN b, char *ch)
{
  GEN p1;
  long av=avma,tx=typ(a),la,lx,i,e;
  if(tx>2) err(forer1);
  if(tx==1)
  {
    la=lx=lgef(a);if(lx<3) lx=3;
    if(!gcmp0(b)) {e=(gexpo(b)>>TWOPOTBITS_IN_LONG)+3;if(e>lx) lx=e;}
    p1=newbloc(lx);for(i=0;i<la;i++) p1[i]=a[i];
    setlg(p1,lx);p1[-1]=(long)ep->value;ep->value=(void*)p1;
  }
  else newvalue(ep,a);
  p1=(GEN)ep->value;
  while (gcmp(p1,b)<=0)
  {
    lisseq(ch);gaddgsz(p1,1,p1);avma=av;
  }
  killvalue(ep);return gnil;
}

GEN
forstep(entree *ep, GEN a, GEN b, GEN s, char *ch)
{
  GEN p1;
  long av=avma,tx=typ(a),lx,la,s1,i,e;
  s1=gsigne(s);if(!s1) err(forer2);
  if(tx>2) err(forer1);
  if(tx==1)
  {
    la=lx=lgef(a);if(lx<3) lx=3;
    if(!gcmp0(b)) {e=(gexpo(b)>>TWOPOTBITS_IN_LONG)+3;if(e>lx) lx=e;}
    if(typ(s)==2) {affir(a,p1=cgetr(lg(s)));newvalue(ep,p1);avma=av;}
    else
    {
      p1=newbloc(lx);for(i=0;i<la;i++) p1[i]=a[i];
      setlg(p1,lx);p1[-1]=(long)ep->value;ep->value=(void*)p1;
    }
  }
  else newvalue(ep,a);
  p1=(GEN)ep->value;
  if(s1>0)
    while (gcmp(p1,b)<= 0)
    {
      lisseq(ch);gaddz(p1,s,p1);avma = av;
    }
  else
    while (gcmp(p1,b)>= 0)
    {
      lisseq(ch);gaddz(p1,s,p1);avma = av;
    }
  killvalue(ep);return gnil;
}

GEN
forprime(entree *ep, GEN a, GEN b, char *ch)
{
  GEN p1;
  long av=avma,prime=0;
  byteptr p=diffptr;
  
  newvalue(ep,gun); p1 = (GEN)ep->value;
  while((*p)&&gcmpgs(a,prime)>0) prime += *p++;
  if(!*p) err(recprimer);
  affsi(prime,p1);
  while(gcmp(p1,b)<=0)
  {
    if(!*p) err(recprimer);
    lisseq(ch);addsiz(*p++,p1,p1);avma=av;
  }
  killvalue(ep);return gnil;
}

GEN
fordiv(entree *ep, GEN a, char *ch)
{
  long i,l,av2,av=avma;
  GEN p1,t=divisors(a);
  l=lg(t);
  newvalue(ep,a); p1=(GEN)ep->value;
  av2=avma;
  for(i=1;i<l;i++)
  {
    gaffect((GEN)t[i],p1);
    lisseq(ch);
    avma=av2;
  }
  avma=av;
  killvalue(ep);return gnil;
}

GEN
forvec(entree *ep, GEN x, char *ch)
{
  GEN p1,y,a,b,c,pinit;
  long av=avma,tx=typ(x),ty,tz,n,e,la,lx,i,j;
  
  if((tx!=17)&&(tx!=18)) err(talker,"not a vector in forvec");
  n=lg(x)-1;if(!n) return gnil;
  y=newbloc(n+1);y[0]=evaltyp(17)+evalpere(1)+evallg(n+1);
  for(i=1;i<=n;i++)
  {
    c=(GEN)x[i];ty=typ(c);
    if(((ty!=17)&&(ty!=18))||(lg(c)!=3))
      err(talker,"not a vector of two-component vectors in forvec");
    a=(GEN)c[1];b=(GEN)c[2];tz=typ(a);
    if(gcmp(a,b)>0) {avma=av;return gnil;}
    if(tz==1)
    {
      la=lx=lgef(a);if(lx<3) lx=3;
      if(!gcmp0(b)) {e=(gexpo(b)>>TWOPOTBITS_IN_LONG)+3;if(e>lx) lx=e;}
      p1=newbloc(lx);for(j=0;j<la;j++) p1[j]=a[j];
      setlg(p1,lx);
    }
    else p1=gclone(a);
    y[i]=(long)p1;
  }
  y[-1]=(long)ep->value;ep->value=(void*)y;pinit=gclone(y);
  do
  {
    lisseq(ch);i=n+1;
    do
    {
      i--;
      if(i<n) gaffect((GEN)pinit[i+1],(GEN)y[i+1]);
      if(i) gaddgsz((GEN)y[i],1,(GEN)y[i]);
    }
    while(i&&(gcmp((GEN)y[i],(GEN)(((GEN)x[i])[2]))>0));
    avma=av;
  }
  while(i);
  killbloc(pinit);for(i=1;i<=n;i++) killbloc((GEN)y[i]);
  killvalue(ep);return gnil;
}

/********************************************************************/
/********************************************************************/
/**                                                                **/
/**       SOMMES,PRODUITS,VECTEURS,MATRICES ET RECURRENCES         **/
/**                                                                **/
/********************************************************************/
/********************************************************************/

GEN
somme(entree *ep, GEN x, GEN a, GEN b, char *ch)
{
  GEN   y,z,p1;
  long  av = avma,tetpil,limite=(av+bot)/2;
  
  newvalue(ep,gun); p1 = (GEN)ep->value; gaffect(a, p1);
  y = x;
  tetpil = avma;
  if(gcmp(a,b)>0)
  {
    gaddgsz(p1,-1,p1);
    while(gcmp(p1,b)>0)
    {
      if(avma<limite) {tetpil=avma; y=gerepile(av,tetpil,gcopy(y));}
      z=gneg(lisexpr(ch));
      tetpil=avma;y=gadd(y,z);
      gaddgsz(p1,-1,p1);
    }
  }
  else
    do
    {
      if(avma<limite) {tetpil=avma; y=gerepile(av,tetpil,gcopy(y));}
      z=lisexpr(ch);
      tetpil=avma;y=gadd(y,z);
      gaddgsz(p1,1,p1);
    }
    while(gcmp(p1,b)<=0);
  killvalue(ep);
  return gerepile(av,tetpil,y);
}

GEN
suminf(entree *ep, GEN a, char *ch, long prec)
{
  GEN   y,z,p1;
  long  av=avma,tetpil,limite=(bot+avma)/2,fl=0;
  
  newvalue(ep,gun); p1 = (GEN)ep->value; gaffect(a, p1);
  y=gun;
  do
  {
    if (avma < limite) {tetpil = avma; y=gerepile(av,tetpil,gcopy(y));}
    z=lisexpr(ch); y=gadd(y,z);
    gaddgsz(p1,1,p1);
    if((!gcmp0(z))&&(gexpo(z)>gexpo(y)-((prec-2)<<TWOPOTBITS_IN_LONG)-5)) fl=0;
    else fl++;
  }
  while(fl<3);
  killvalue(ep);
  tetpil=avma; return gerepile(av,tetpil,gsubgs(y,1));
}

GEN
produit(entree *ep, GEN x, GEN a, GEN b, char *ch)
{
  GEN   y,z,p1;
  long  av=avma,tetpil,limite=(av+bot)/2;
  
  newvalue(ep,gun); p1 = (GEN)ep->value; gaffect(a,p1);
  y = x;
  tetpil = avma;
  if(gcmp(a,b)>0)
  {
    gaddgsz(p1,-1,p1);
    while(gcmp(p1,b)>0)
    {
      if(avma<limite) {tetpil=avma; y=gerepile(av,tetpil,gcopy(y));}
      z=lisexpr(ch);
      tetpil=avma;y=gdiv(y,z);
      gaddgsz(p1,-1,p1);
    }
  }
  else
    do
    {
      if(avma<limite) {tetpil=avma; y=gerepile(av,tetpil,gcopy(y));}
      z=lisexpr(ch);
      tetpil=avma;y=gmul(y,z);
      gaddgsz(p1,1,p1);
    }
    while(gcmp(p1,b)<=0);
  killvalue(ep);
  return gerepile(av,tetpil,y);
}

GEN
prodinf(entree *ep, GEN a, char *ch, long prec)
{
  GEN   y,z,p1;
  long  av=avma,tetpil,limite=(av+bot)/2,fl=0;
  
  newvalue(ep, gun); p1=(GEN)ep->value; gaffect(a,p1);
  affsr(1,y=cgetr(prec));
  do
  {
    if(avma<limite) {tetpil=avma; y=gerepile(av,tetpil,gcopy(y));}
    z=lisexpr(ch); y=gmul(y,z);
    if(DEBUGLEVEL>4) {outerr(p1);outerr(z);outerr(y);}
    gaddgsz(p1,1,p1);
    if(gexpo(gsubgs(z,1))>-((prec-2)<<TWOPOTBITS_IN_LONG)-5) fl=0;else fl++;
  }
  while(fl<3);
  killvalue(ep);
  tetpil = avma; return gerepile(av,tetpil,gcopy(y));
}

GEN
prodinf1(entree *ep, GEN a, char *ch, long prec)
{
  GEN   y,z,p1,p2;
  long  av=avma,tetpil,limite=(av+bot)/2,fl=0;
  
  newvalue(ep,gun); p1=(GEN)ep->value; gaffect(a,p1);
  affsr(1,y=cgetr(prec));
  do
  {
    if(avma<limite) {tetpil=avma; y=gerepile(av,tetpil,gcopy(y));}
    p2=lisexpr(ch);z=gaddsg(1,p2);
    y=gmul(y,z);
    gaddgsz(p1,1,p1);
    if((!gcmp0(z))&&(gexpo(p2)>-((prec-2)<<TWOPOTBITS_IN_LONG)-5)) fl=0;else fl++;
  }
  while(fl<3);
  killvalue(ep);
  tetpil = avma; return gerepile(av,tetpil,gcopy(y));
}

GEN
prodeuler(entree *ep, GEN a, GEN b, char *ch, long prec)
{
  GEN y,z,p1;
  long av=avma,tetpil,prime=0;
  byteptr p=diffptr;
  
  newvalue(ep,gun); p1 = (GEN)ep->value;
  while((*p)&&gcmpgs(a,prime)>0) prime += *p++;
  if(!*p) err(recprimer);
  affsi(prime,p1);affsr(1,y=cgetr(prec));
  if(gcmp(p1,b)>0)
  {
    if(!gcmp1(gsub(a,b))) err(recer1);
    tetpil=avma;
    y=gcopy(y);
  }
  else
    do
    {
      if(!*p) err(recprimer);
      z=lisexpr(ch);
      tetpil=avma; y=gmul(y,z);
      addsiz(*p++,p1,p1);
    }
    while(gcmp(p1,b)<=0);
  killvalue(ep);
  return gerepile(av,tetpil,y);
}

GEN
vecteur(entree *ep, GEN nmax, char *ch)
{
  GEN   y,p1,t;
  long  i,m;
  
  if((typ(nmax)!=1) || (signe(nmax)<0)) err(vecer1);
  if(gcmpgs(nmax,MAXLEN)>=0)
    err(talker,"too many components in vector");
  m=itos(nmax);
  y=cgetg(m+1 ,17);
  newvalue(ep, gun); p1 = (GEN)ep->value;
  for(i=1;i<=m;i++)
  {
    affsi(i,p1);
    t=lisexpr(ch);
    y[i] = isonstack(t) ? (long)t : lcopy(t);
  }
  killvalue(ep);
  return y;
}

GEN
vvecteur(entree *ep, GEN nmax, char *ch)
{
  GEN   y=vecteur(ep,nmax,ch);
  settyp(y,18);
  return y;
}

GEN
matrice(entree *ep1, entree *ep2, GEN nlig, GEN ncol, char *ch)
{
  GEN   y,z,t,p1,p2;
  long  i,j,m,n;
  
  if((typ(nlig)!=1) || (signe(nlig)<=0)) err(mater1);
  if((typ(ncol)!=1) || (signe(ncol)<=0)) err(mater1);
  if(gcmpgs(nlig,MAXLEN)>=0)
    err(talker,"too many rows in matrix");
  if(gcmpgs(ncol,MAXLEN)>=0)
    err(talker,"too many columns in matrix");
  n=itos(nlig); m=itos(ncol);
  newvalue(ep1, gun); p1 = (GEN)ep1->value;
  newvalue(ep2, gun); p2 = (GEN)ep2->value;
  y=cgetg(m+1 ,19);
  for(i=1;i<=m;i++)
  {
    affsi(i,p2);
    z=cgetg(n+1 ,18);
    y[i]=(long)z;
    for(j=1;j<=n;j++)
    {
      affsi(j,p1);
      t=lisexpr(ch);
      z[j] = isonstack(t) ? (long)t : lcopy(t);
    }
  }
  killvalue(ep1); killvalue(ep2); 
  return y;
}

GEN
divsomme(entree *ep, GEN num, char *ch)
{
  GEN   y,z,p1;
  long  av=avma,tetpil,d,n,d2;
  
  newvalue(ep, gun); p1 = (GEN)ep->value;
  n=itos(num);/* provisoire */
  tetpil=av=avma;y=gzero;
  for(d = d2 = 1; d2 < n; d++, d2 += d+d-1)
    if (!(n%d))
    {
      affsi(d,p1);y=gadd(y, lisexpr(ch));
      affsi(n/d,p1); z = lisexpr(ch);
      tetpil=avma; y=gadd(y,z);
    }
  if (d2 == n)
  {
    affsi(d,p1); z = lisexpr(ch);
    tetpil=avma; y=gadd(y,z);
  }
  killvalue(ep);
  return gerepile(av,tetpil,y);
}

/********************************************************************/
/********************************************************************/
/**                                                                **/
/**                   CALCUL D'UNE INTEGRALE                       **/
/**                   ( Methode de ROMBERG )                       **/
/**                                                                **/
/********************************************************************/
/********************************************************************/

GEN
qromb(entree *ep, GEN a, GEN b, char *ch, long prec)
     
#define JMAX 25
#define JMAXP JMAX+3
#define KLOC 5
     
{
  GEN     ss,dss,s,h,q1,p1,p2,p3,p4,p,qlint,del,sz,x,sum;
  long    av,tetpil,j,j1,j2,lim,l,it,tnm,sig;
  
  av=avma;l=prec;
  if(typ(a)!=2) { p=cgetr(prec);gaffect(a,p);a=p;}
  if(typ(b)!=2) { p=cgetr(prec);gaffect(b,p);b=p;}
  qlint=subrr(b,a);sig=signe(qlint);
  if(!sig) {avma=av;return gzero;}
  newvalue(ep,cgetr(prec)); q1=(GEN)ep->value;
  if(sig<0) {setsigne(qlint,1);s=a;a=b;b=s;}
  p3=cgetg(KLOC+1,17);p4=cgetg(KLOC+1,17);
  s=cgetg(JMAXP,17);h=cgetg(JMAXP,17);
  affsr(1,(GEN)(h[1]=lgetr(l)));
  gaffect(a,q1); p1=lisexpr(ch); if (!isonstack(p1)) p1=gcopy(p1);
  gaffect(b,q1); p2=lisexpr(ch);
  s[1]=lmul2n(gmul(qlint,gadd(p1,p2)),-1);it=1;sz=gmul(gzero,(GEN)s[1]);
  s[2]=s[1];h[2]=lshiftr((GEN)h[1],-2);
  for(j=2;j<=JMAX;j++)
  {
    tnm=it;del=divrs(qlint,tnm);x=addrr(a,shiftr(del,-1));
    for(sum=sz,j1=1;j1<=it;j1++,x=addrr(x,del))
    {
      affrr(x,q1);p1=lisexpr(ch);sum=gadd(sum,p1);
    }
    it*=2;
    s[j]=lmul2n(gadd((GEN)s[j-1],gmul(sum,del)),-1);
    if(j>=KLOC)
    {
      for(j1=1;j1<=KLOC;j1++)
      {
	p3[j1]=s[j1+j-KLOC];p4[j1]=h[j1+j-KLOC];
      }
      tetpil=avma;ss=polint(p4,p3,gzero,&dss);
      j1=gexpo(ss);j2=gexpo(dss);lim=((prec-2)<<TWOPOTBITS_IN_LONG)-j-5;
      if((((j1-j2)>lim))||((j1< -lim)&&(j2<j1-1)))
      {
	if(gcmp0(gimag(ss))) ss=greal(ss);
	tetpil=avma;
	killvalue(ep);
	return gerepile(av,tetpil,gmulsg(sig,ss));
      }
    }
    s[j+1]=s[j];h[j+1]=lshiftr((GEN)h[j],-2);
  }
  err(intger2);return gnil;
}         

GEN
qromo(entree *ep, GEN a, GEN b, char *ch, long prec)
     
#undef JMAX
#define JMAX 16
#define JMAXP JMAX+3
#define KLOC 5
     
{
  GEN     ss,dss,s,h,q1,sz,p1,p3,p4,p,qlint,del,ddel,x,sum;
  long    av,tetpil,j,j1,j2,lim,l,it,tnm,sig;
  
  av=avma;l=prec;
  if(typ(a)!=2) { p=cgetr(prec);gaffect(a,p);a=p;}
  if(typ(b)!=2) { p=cgetr(prec);gaffect(b,p);b=p;}
  qlint=subrr(b,a);sig=signe(qlint);
  if(!sig) {avma=av;return gzero;}
  if(sig<0) {setsigne(qlint,1);s=a;a=b;b=s;}
  newvalue(ep,cgetr(prec)); q1=(GEN)ep->value;
  p3=cgetg(KLOC+1,17);p4=cgetg(KLOC+1,17);
  s=cgetg(JMAXP,17);h=cgetg(JMAXP,17);
  affsr(1,(GEN)(h[1]=lgetr(l)));affrr(shiftr(addrr(a,b),-1),q1);p1=lisexpr(ch);
  s[1]=lmul(qlint,p1);it=1;sz=gmul(gzero,(GEN)s[1]);
  s[2]=s[1];h[2]=ldivrs((GEN)h[1],9);
  for(j=2;j<=JMAX;j++)
  {
    tnm=it;del=divrs(qlint,3*tnm);ddel=shiftr(del,1);
    x=addrr(a,shiftr(del,-1));sum=sz;
    for(j1=1;j1<=it;j1++)
    {
      affrr(x,q1);p1=lisexpr(ch);sum=gadd(sum,p1);x=addrr(x,ddel);
      affrr(x,q1);p1=lisexpr(ch);sum=gadd(sum,p1);x=addrr(x,del);
    }
    it*=3;
    s[j]=ladd(gdivgs((GEN)s[j-1],3),gmul(sum,del));
    if(j>=KLOC)
    {
      for(j1=1;j1<=KLOC;j1++)
      {
	p3[j1]=s[j1+j-KLOC];p4[j1]=h[j1+j-KLOC];
      }
      tetpil=avma;ss=polint(p4,p3,gzero,&dss);
      j1=gexpo(ss);j2=gexpo(dss);lim=((prec-2)<<TWOPOTBITS_IN_LONG)-(3*j/2)-5;
      if((((j1-j2)>lim))||((j1< -lim)&&(j2<j1-1)))
      {
	if(gcmp0(gimag(ss))) ss=greal(ss);
	tetpil=avma; killvalue(ep);
	return gerepile(av,tetpil,gmulsg(sig,ss));
      }
    }
    s[j+1]=s[j];h[j+1]=ldivrs((GEN)h[j],9);
  }
  err(intger2);return gnil;
}         

GEN
qromi(entree *ep, GEN a, GEN b, char *ch, long prec)
                
#undef JMAX
#define JMAX 16
#define JMAXP JMAX+3
#define KLOC 5
     
{
  GEN     ss,dss,s,h,q1,sz,p1,p3,p4,p,qlint,del,ddel,x,sum;
  long    av,tetpil,j,j1,j2,lim,l,it,tnm,sig;
  
  av=avma;l=prec;
  p=cgetr(prec);gaffect(ginv(a),p);a=p;
  p=cgetr(prec);gaffect(ginv(b),p);b=p;
  qlint=subrr(b,a);sig= -signe(qlint);
  if(!sig) {avma=av;return gzero;}
  if(sig>0) {setsigne(qlint,1);s=a;a=b;b=s;}
  newvalue(ep,cgetr(prec)); q1=(GEN)ep->value;
  p3=cgetg(KLOC+1,17);p4=cgetg(KLOC+1,17);
  s=cgetg(JMAXP,17);h=cgetg(JMAXP,17);
  affsr(1,(GEN)(h[1]=lgetr(l)));x=divsr(2,addrr(a,b));
  affrr(x,q1);
  p1=gmul(lisexpr(ch),mulrr(x,x));s[1]=lmul(qlint,p1);it=1;
  sz=gmul(gzero,(GEN)s[1]);s[2]=s[1];h[2]=ldivrs((GEN)h[1],9);
  for(j=2;j<=JMAX;j++)
  {
    tnm=it;del=divrs(qlint,3*tnm);ddel=shiftr(del,1);
    x=addrr(a,shiftr(del,-1));sum=sz;
    for(j1=1;j1<=it;j1++)
    {
      divsrz(1,x,q1);p1=gmul(lisexpr(ch),mulrr(q1,q1));
      sum=gadd(sum,p1);x=addrr(x,ddel);
      divsrz(1,x,q1);p1=gmul(lisexpr(ch),mulrr(q1,q1));
      sum=gadd(sum,p1);x=addrr(x,del);
    }
    it*=3;
    s[j]=ladd(gdivgs((GEN)s[j-1],3),gmul(sum,del));
    if(j>=KLOC)
    {
      for(j1=1;j1<=KLOC;j1++)
      {
	p3[j1]=s[j1+j-KLOC];p4[j1]=h[j1+j-KLOC];
      }
      tetpil=avma;ss=polint(p4,p3,gzero,&dss);
      j1=gexpo(ss);j2=gexpo(dss);lim=((prec-2)<<TWOPOTBITS_IN_LONG)-(3*j/2)-5;
      if((((j1-j2)>lim))||((j1< -lim)&&(j2<j1-1)))
      {
	if(gcmp0(gimag(ss))) ss=greal(ss);
	tetpil=avma; killvalue(ep);
	return gerepile(av,tetpil,gmulsg(sig,ss));
      }
    }
    s[j+1]=s[j];h[j+1]=ldivrs((GEN)h[j],9);
  }
  err(intger2);return gnil;
}         

GEN
rombint(entree *ep, GEN a, GEN b, char *ch, long prec)
{
  GEN     s,p1,p2,p3;
  static  long p4[3]={evaltyp(1)+evalpere(MAXPERE)+evallg(3),evalsigne(-1)+evallgef(3),1};
  long    l,av,tetpil;
  
  l=gcmp(b,a);
  if(!l) return gzero;
  if(l<0) {s=a;a=b;b=s;}
  av=avma;
  if(gcmpgs(b,100)>=0)
  {
    if(gcmpgs(a,1)>=0) return qromi(ep,a,b,ch,prec);
    p1=qromi(ep,gun,b,ch,prec);
    if(gcmpgs(a,-100)>=0)
    {
      p2=qromo(ep,a,gun,ch,prec);tetpil=avma;
      return gerepile(av,tetpil,gadd(p1,p2));
    }
    p2=qromo(ep,p4,gun,ch,prec);p3=gadd(p2,qromi(ep,a,p4,ch,prec));
    tetpil=avma;return gerepile(av,tetpil,gadd(p1,p3));
  }
  if(gcmpgs(a,-100)>=0) return qromo(ep,a,b,ch,prec);
  if(gcmpgs(b,-1)>=0)
  {
    p1=qromi(ep,a,p4,ch,prec);p2=qromo(ep,p4,b,ch,prec);tetpil=avma;
    return gerepile(av,tetpil,gadd(p1,p2));
  }
  return qromi(ep,a,b,ch,prec);
}         

/********************************************************************/
/********************************************************************/
/**                                                                **/
/**                    SOMMATION DE SERIES                         **/
/**                                                                **/
/********************************************************************/
/********************************************************************/

void
eulsum(GEN *sum, GEN term, long jterm, GEN *tab, long *dsum, long prec)
{
  long j;
  static long nterm;
  GEN  tmp,dum,p1;
  static GEN unreel;
  
  if(jterm==1)
  {
    nterm=1;tab[1]=term;*sum=gmul2n(term,-1);affsr(1,unreel=cgetr(prec));
  }
  else
  {
    tmp=tab[1];tab[1]=gmul(unreel,term);
    for(j=1;j<nterm;j++)
    {
      dum=tab[j+1];tab[j+1]=gmul2n(gadd(tab[j],tmp),-1);tmp=dum;
    }
    tab[nterm+1]=gmul2n(gadd(tab[nterm],tmp),-1);
    if(gcmp(gabs(tab[nterm+1],prec),gabs(tab[nterm],prec))<=0)
      p1=gmul2n(tab[++nterm],-1);
    else p1=tab[nterm+1];
    *sum=gadd(*sum,p1);*dsum=gexpo(p1);
  }
}

GEN
sumalt2(entree *ep, GEN a, char *ch, long prec)

#define JMIT 10000
{
  long av,tetpil,j,nterm,jterm;
  GEN  p1,sum,sum0,q1,tmp,dum,*tab,unreel;
  
  newvalue(ep,gun); q1=(GEN)ep->value; gaffect(a,q1);
  tab=(GEN *)newbloc(JMIT);
  av=avma;
  p1=lisexpr(ch);affsr(1,unreel=cgetr(prec));
  nterm=1;tab[1]=p1;sum=gmul2n(p1,-1);
  for(jterm=2;jterm<=JMIT;jterm++)
  {
    tmp=tab[1];a=gaddsg(1,a);gaffect(a,q1);tab[1]=gmul(unreel,lisexpr(ch));
    for(j=1;j<nterm;j++)
    {
      dum=tab[j+1];tab[j+1]=gmul2n(gadd(tab[j],tmp),-1);tmp=dum;
    }
    tab[nterm+1]=gmul2n(gadd(tab[nterm],tmp),-1);
    if(gcmp(gabs(tab[nterm+1],prec),gabs(tab[nterm],prec))<=0)
      p1=gmul2n(tab[++nterm],-1);
    else p1=tab[nterm+1];
    sum0=sum;sum=gadd(sum,p1);
    if((gcmp0(p1)||(gexpo(p1)<-((prec-2)<<TWOPOTBITS_IN_LONG)+5)||(gegal(sum,sum0)))&&(jterm>=10))
    {
      tetpil=avma;
      killbloc((GEN)tab);
      killvalue(ep);
      return gerepile(av,tetpil,gcopy(sum));
    }
  }
  err(eulsumer1);return gnil;
}

GEN
sumalt1(entree *ep, GEN a, char *ch, long prec)
{
  long av,tetpil,n,N;
  GEN  q1,unreel,comb,w,u,sum;
  
  newvalue(ep,gun);q1=(GEN)ep->value;gaffect(a,q1);
  av=avma;affsr(1,unreel=cgetr(prec));sum=gmul(unreel,lisexpr(ch));
  N=(((prec-2)<<TWOPOTBITS_IN_LONG)+7)/3;
  for(n=1;n<N;n++)
  {
    a=gaddsg(1,a);gaffect(a,q1);sum=gadd(sum,lisexpr(ch));
  }
  comb=gneg(gun);w=gmul2n(gun,N);u=gmul2n(unreel,-N);
  for(n=0;n<N;n++)
  {
    a=gaddsg(1,a);gaffect(a,q1);
    w=addii(w,comb);sum=gadd(sum,gmul(gmul(lisexpr(ch),w),u));
    comb=divis(mulsi(N-n,comb),n+1);
  }
  tetpil=avma;killvalue(ep);
  return gerepile(av,tetpil,gcopy(sum));
}

GEN
sumalt4(entree *ep, GEN a, char *ch, long prec)
{
  long av,tetpil,k,N;
  GEN  q1,unreel,s,t,az,c,x;
  
  newvalue(ep,gun);q1=(GEN)ep->value;gaffect(a,q1);
  av=avma;affsr(1,unreel=cgetr(prec));
  az=gun;c=gun;s=gzero;t=gzero;
  N=(long)(0.4*(((prec-2)<<TWOPOTBITS_IN_LONG)+7));
  for(k=0;k<N;k++)
  {
    x=gmul(unreel,lisexpr(ch));s=k?gadd(s,x):gmul(x,unreel);
    t=gadd(t,gmul(x,c));
    az=divii(mulii(mulss(N-k,N+k),shifti(az,1)),mulss(k+1,k+k+1));
    c=addii(az,c);
    if(k<N-1) {a=gaddsg(1,a);gaffect(a,q1);}
  }
  x=gdiv(t,c);tetpil=avma;killvalue(ep);
  return gerepile(av,tetpil,gsub(s,x));
}

GEN
sumalt(entree *ep, GEN a, char *ch, long prec)
{
  long av,tetpil,k,N;
  GEN  q1,unreel,s,az,c,x,e1,d;
  
  newvalue(ep,gun);q1=(GEN)ep->value;gaffect(a,q1);
  av=avma;affsr(1,unreel=cgetr(prec));
  e1=addsr(3,gsqrt(stoi(8),prec));
  N=(long)(0.4*(((prec-2)<<TWOPOTBITS_IN_LONG)+7));
  d=gpuigs(e1,N);d=shiftr(addrr(d,divsr(1,d)),-1);
  az=negi(gun);c=d;s=gzero;
  for(k=0;k<N;k++)
  {
    x=lisexpr(ch);c=gadd(az,c);
    if(DEBUGLEVEL>4) outerr(c);
    s=gadd(s,gmul(x,c));
    az=divii(mulii(mulss(N-k,N+k),shifti(az,1)),mulss(k+1,k+k+1));
    if(k<N-1) {a=gaddsg(1,a);gaffect(a,q1);}
  }
  tetpil=avma;killvalue(ep);
  return gerepile(av,tetpil,gdiv(s,d));
}

GEN
sumalt3(entree *ep, GEN a, char *ch, long prec)
{
  long av,tetpil,k,N;
  GEN  q1,unreel,s,az,c,x,e1,d;
  
  newvalue(ep,gun);q1=(GEN)ep->value;gaffect(a,q1);
  av=avma;affsr(1,unreel=cgetr(prec));
  e1=addsr(3,gsqrt(stoi(8),prec));
  N=(long)(0.4*(((prec-2)<<TWOPOTBITS_IN_LONG)+7));
  d=gpuigs(e1,N);d=shiftr(addrr(d,divsr(1,d)),-1);
  az=negr(unreel);c=d;s=gzero;
  for(k=0;k<N;k++)
  {
    x=lisexpr(ch);c=gadd(az,c);
    s=gadd(s,gmul(x,c));
    az=divri(mulir(mulss(N-k,N+k),shiftr(az,1)),mulss(k+1,k+k+1));
    if(k<N-1) {a=gaddsg(1,a);gaffect(a,q1);}
  }
  tetpil=avma;killvalue(ep);
  return gerepile(av,tetpil,gdiv(s,d));
}

GEN
sumpos(entree *ep, GEN a, char *ch, long prec)
{
  long av,tetpil,k,kk,N;
  GEN  p1,r,q1,unreel,s,az,c,x,e1,d,stock;
  
  av=avma;a=gsubgs(a,1);affsr(1,unreel=cgetr(prec));
  e1=addsr(3,gsqrt(stoi(8),prec));
  N=(long)(0.4*(((prec-2)<<TWOPOTBITS_IN_LONG)+7));
  stock=cgetg(N+1,17);for(k=1;k<=N;k++) stock[k]=(long)gnil;
  d=gpuigs(e1,N);d=shiftr(addrr(d,divsr(1,d)),-1);
  az=negi(gun);c=d;s=gzero;
  p1=cgeti(prec+1);p1[1]=evalsigne(1)+evallgef(1+prec);
  newvalue(ep,p1);q1=(GEN)ep->value;
  for(k=0;k<N;k++)
  {
    if(odd(k)&&(stock[k]!=(long)gnil)) x=(GEN)stock[k];
    else
    {
      x=gzero;r=stoi(2*k+2);kk=0;
      do
      {
	gaddz(r,a,q1);p1=gmul2n(gmul(unreel,lisexpr(ch)),kk);
	x=gadd(x,p1);r=shifti(r,1);kk++;
      }
      while((kk<=1)||(gexpo(p1)>=(-((prec-2)<<TWOPOTBITS_IN_LONG)+5)));
      if(2*k<N) stock[2*k+1]=(long)x;
      r=stoi(k+1);gaddz(r,a,q1);
      x=gadd(gmul(unreel,lisexpr(ch)),gmul2n(x,1));
    }
    c=gadd(az,c);s=odd(k)?gsub(s,gmul(x,c)):gadd(s,gmul(x,c));
    if(DEBUGLEVEL>4) outerr(c);
    az=divii(mulii(mulss(N-k,N+k),shifti(az,1)),mulss(k+1,k+k+1));
  }
  tetpil=avma;killvalue(ep);
  return gerepile(av,tetpil,gdiv(s,d));
}

GEN
sumposold(entree *ep, GEN a, char *ch, long prec)
{
  long av,tetpil,k,jterm,dsum;
  GEN  sum,term,q1,p1,*tab,unreel,r;
  
  tab=(GEN *)newbloc(JMIT);
  av=avma;a=gsubgs(a,1);
  affsr(1,unreel=cgetr(prec));
  p1=cgeti(prec+1);p1[1]=evalsigne(1)+evallgef(1+prec);
  newvalue(ep,p1);q1=(GEN)ep->value;sum=gzero;
  for(jterm=1;jterm<=JMIT;jterm++)
  {
    term=gzero;r=stoi(jterm);k=0;
    do
    {
      gaddz(r,a,q1);p1=gmul2n(gmul(unreel,lisexpr(ch)),k);
      term=gadd(term,p1);r=shifti(r,1);k++;
    }
    while(gexpo(p1)>=(-((prec-2)<<TWOPOTBITS_IN_LONG)+5));
    if(!odd(jterm)) term=gneg(term);
    eulsum(&sum,term,jterm,tab,&dsum,prec);
    if(dsum< -((prec-2)<<TWOPOTBITS_IN_LONG)+5)
    {
      tetpil=avma;killbloc((GEN)tab);killvalue(ep);
      return gerepile(av,tetpil,gcopy(sum));
    }
  }
  killbloc((GEN)tab);killvalue(ep);
  err(eulsumer1);return gnil;
}

/********************************************************************/
/********************************************************************/
/**                                                                **/
/**                        TRACE GROSSIER                          **/
/**                                                                **/
/********************************************************************/
/********************************************************************/


GEN
plot(entree *ep, GEN a, GEN b, char *ch)
     
#define ISCR 60
#define JSCR 21
#define BLANK ' '
#define ZERO '-'
#define YY '|'
#define XX '-'
#define FF 'x'
     
{
  long av,av2,jz,j,i,sig;
  GEN p1,p2,ysml,ybig,x,diff,dyj,dx,y[ISCR+1];
  char scr[ISCR+1][JSCR+1], thestring[100];
  
  sig=gcmp(b,a);if(!sig) return gnil;
  av=avma; pariputc('\n');
  if(sig<0) {x=a;a=b;b=x;}
  newvalue(ep,cgetr(5)); x=(GEN)ep->value;
  for(i=1;i<=ISCR;i++) y[i]=cgetr(3);
  gaffect(a,x);
  dx=gdivgs(gsub(b,a),(ISCR-1));ysml=gzero;ybig=gzero;
  for(j=1;j<=JSCR;j++) scr[1][j]=scr[ISCR][j]=YY;
  for(i=2;i<ISCR;i++)
  {
    scr[i][1]=scr[i][JSCR]=XX;
    for(j=2;j<JSCR;j++) scr[i][j]=BLANK;
  }
  av2=avma;
  for(i=1;i<=ISCR;i++)
  {
    gaffect(lisexpr(ch),y[i]);
    if(gcmp(y[i],ysml)<0) ysml=y[i];
    if(gcmp(y[i],ybig)>0) ybig=y[i];
    gaddz(x,dx,x);avma=av2;
  }
  diff=gsub(ybig,ysml);
  if(gcmp0(diff)) {ybig=gaddsg(1,ybig);diff=gun;}
  dyj=gdivsg(JSCR-1,diff);jz=1-gtolong(gmul(ysml,dyj));
  av2=avma;
  for(i=1;i<=ISCR;i++)
  {
    scr[i][jz]=ZERO;j=1+gtolong(gmul(gsub(y[i],ysml),dyj));
    scr[i][j]=FF;avma=av2;
  }
  p1=cgetr(4);gaffect(ybig,p1);
  sprintf(thestring, " %10.3f ",rtodbl(p1)); pariputs(thestring);
  for(i=1;i<=ISCR;i++)  pariputc(scr[i][JSCR]); pariputc('\n');
  for(j=(JSCR-1);j>=2;j--)
  {
    pariputs("            ");
    for(i=1;i<=ISCR;i++) pariputc(scr[i][j]);
    pariputc('\n');
  }
  p1=cgetr(4);gaffect(ysml,p1);
  sprintf(thestring, " %10.3f ",rtodbl(p1)); pariputs(thestring);
  for(i=1;i<=ISCR;i++)  pariputc(scr[i][1]);
  pariputc('\n');
  p1=cgetr(4);p2=cgetr(4);gaffect(a,p1);gaffect(b,p2);
  sprintf(thestring, "%8s %10.3f %44s %10.3f\n"," ",rtodbl(p1)," ",rtodbl(p2));
  pariputs(thestring);
  killvalue(ep);
  avma=av; pariputc('\n'); return gnil;
}

/********************************************************************/
/********************************************************************/
/**                                                                **/
/**                  RECHERCHE DE ZEROS REELS                      **/
/**                                                                **/
/********************************************************************/
/********************************************************************/

GEN
zbrent(entree *ep, GEN a, GEN b, char *ch, long prec)
{
  long av=avma,tetpil,sig,iter,itmax;
  GEN q1,c,d,e,tol,toli,min1,min2,fa,fb,fc,p,q,r,s,xm;
  
  if(typ(a)!=2) {p=cgetr(prec);gaffect(a,p);a=p;}
  if(typ(b)!=2) {p=cgetr(prec);gaffect(b,p);b=p;}
  sig=cmprr(b,a);if(!sig) {avma = av; return gzero;}
  if(sig<0) {c=a;a=b;b=c;}
  newvalue(ep,a); fa=lisexpr(ch);
  changevalue(ep,b); q1=(GEN)ep->value; fb=lisexpr(ch);
  if(gsigne(fa)*gsigne(fb)>0) err(zbrenter1);
  itmax=(prec<<TWOPOTBITS_IN_LONG)+50;affsr(1,tol=cgetr(3));tol=shiftr(tol,-(((prec-2)<<TWOPOTBITS_IN_LONG)-5));
  fc=fb;
  for(iter=1;iter<=itmax;iter++)
  {
    if(gsigne(fb)*gsigne(fc)>0)
    {
      c=a;fc=fa;d=subrr(b,a);e=d;
    }
    if(gcmp(gabs(fc,0),gabs(fb,0))<0)
    {
      a=b;b=c;c=a;fa=fb;fb=fc;fc=fa;
    }
    toli=mulrr(tol,absr(b));
    xm=shiftr(subrr(c,b),-1);
    if((cmprr(absr(xm),toli)<=0)||gcmp0(fb))
    {
      tetpil=avma;
      killvalue(ep);
      return gerepile(av,tetpil,gcopy(b));
    }
    if((cmprr(absr(e),toli)>=0)&&(gcmp(gabs(fa,0),gabs(fb,0))>0))
    {
      s=gdiv(fb,fa);
      if(cmprr(a,b)==0)
      {
	p=gmul2n(gmul(xm,s),1);q=gsubsg(1,s);
      }
      else
      {
	q=gdiv(fa,fc);r=gdiv(fb,fc);
	p=gmul2n(gmul(gsub(q,r),gmul(xm,q)),1);
	p=gmul(s,gsub(p,gmul(gsub(b,a),gsubgs(r,1))));
	q=gmul(gmul(gsubgs(q,1),gsubgs(r,1)),gsubgs(s,1));
      }
      if(gsigne(p)>0) q=gneg(q);
      p=gabs(p,0);
      min1=gsub(gmulsg(3,gmul(xm,q)),gabs(gmul(q,toli),0));
      min2=gabs(gmul(e,q),0);
      if(gcmp(gmul2n(p,1),gmin(min1,min2))<0) { e=d;d=gdiv(p,q);}
      else {d=xm;e=d;}
    }
    else {d=xm;e=d;}
    a=b;fa=fb;
    if(gcmp(gabs(d,0),toli)>0) b=addrr(b,d);
    else
    {
      if(gsigne(xm)>0) b=addrr(b,absr(toli));
      else b=subrr(b,absr(toli));
    }
    gaffect(b,q1); ;fb=lisexpr(ch);
  }
  err(zbrenter2);return gnil;
}

GEN
lllgen(GEN x)
{
  long lx=lg(x), tx=typ(x),i,j,av,av1;
  GEN g;

  if(tx!=19) err(lller1);
  av=avma;
  g=cgetg(lx,19);
  for(j=1;j<lx;j++) g[j]=lgetg(lx,18);
  for(i=1;i<lx;i++)
    for(j=1;j<=i;j++) coeff(g,i,j)=coeff(g,j,i)=(long)gscal((GEN)x[i],(GEN)x[j]);
  av1=avma;return gerepile(av,av1,lllgramallgen(g,2));
}

GEN
lllkerimgen(GEN x)
{
  long lx=lg(x), tx=typ(x),i,j,av,av1;
  GEN g;

  if(tx!=19) err(lller1);
  av=avma;
  g=cgetg(lx,19);
  for(j=1;j<lx;j++) g[j]=lgetg(lx,18);
  for(i=1;i<lx;i++)
    for(j=1;j<=i;j++) coeff(g,i,j)=coeff(g,j,i)=(long)gscal((GEN)x[i],(GEN)x[j]);
  av1=avma;return gerepile(av,av1,lllgramallgen(g,0));
}

GEN
lllgramgen(GEN x)
{
  return lllgramallgen(x,2);
}

GEN
lllgramkerimgen(GEN x)
{
  return lllgramallgen(x,0);
}

int
pslg(GEN x)
{
  if(gcmp0(x)) return 2;
  return (typ(x)<10)?3:lgef(x);
}

GEN
lllgramallgen(GEN x, long all)
{
  long av=avma,tetpil,lx=lg(x),tx=typ(x),i,j,k,l,n,lim,dec;
  GEN u,B,lam,q,cq,h,la,bb,p1,p2,p3,p4,y,fl;
  int ps1,ps2,flc;

  if(tx!=19) err(lllger1);
  n=lx-1;if(n<=1) return idmat(n);
  if(lg((GEN)x[1])!=lx) err(lllger2);
  av=avma;lim=(avma+bot)>>1;
  B=cgetg(lx+1,18);
  fl=cgetg(lx,17);
  B[1]=un;lam=cgetg(lx,19);
  for(j=1;j<lx;j++) lam[j]=lgetg(lx,18);
  for(i=1;i<lx;i++) 
    for(j=1;j<=i;j++)
    {
      if((j<i)&&(!signe((GEN)fl[j]))) coeff(lam,i,j)=coeff(lam,j,i)=zero;
      else
      {
	u=(GEN)coeff(x,i,j);
	if(typ(u)>10) err(lllger4);
	for(k=1;k<j;k++)
	  if(signe((GEN)fl[k])) u=gdiv(gsub(gmul((GEN)B[k+1],u),gmul(gcoeff(lam,i,k),gcoeff(lam,j,k))),(GEN)B[k]);
	if(j<i) {coeff(lam,i,j)=(long)u;coeff(lam,j,i)=zero;}
	else 
	{
	  if(!gcmp0(u)) {B[i+1]=(long)u;coeff(lam,i,i)=fl[i]=un;}
	  else {B[i+1]=B[i];coeff(lam,i,i)=fl[i]=zero;}
	}
      }
    }
  k=2;h=idmat(n);flc=0;
  for(;;)
  {
    u=(GEN)coeff(lam,k,k-1);
    if(pslg(u)>=pslg((GEN)B[k]))
    {
      q=gdeuc(u,(GEN)B[k]);cq=gdivsg(1,content(q));q=gmul(q,cq);flc=1;
      h[k]=lsub(gmul(cq,(GEN)h[k]),gmul(q,(GEN)h[k-1]));
      coeff(lam,k,k-1)=lsub(gmul(cq,gcoeff(lam,k,k-1)),gmul(q,(GEN)B[k]));
      for(i=1;i<k-1;i++) 
	coeff(lam,k,i)=lsub(gmul(cq,gcoeff(lam,k,i)),gmul(q,gcoeff(lam,k-1,i)));
    }
    ps1=pslg(gmul((GEN)B[k],(GEN)B[k]));
    ps2=pslg(gadd(p3=gmul((GEN)B[k-1],(GEN)B[k+1]),p4=gmul(la=gcoeff(lam,k,k-1),gcoeff(lam,k,k-1))));
    if(signe((GEN)fl[k-1])&&((ps1>ps2)||((ps1==ps2)&&flc)||(!signe((GEN)fl[k]))))
    {
      p1=(GEN)h[k-1];h[k-1]=h[k];h[k]=(long)p1;
      if((ps1==ps2)&&flc) flc=0;else flc=1;
      for(j=1;j<=k-2;j++) 
      {
	p1=gcoeff(lam,k-1,j);coeff(lam,k-1,j)=coeff(lam,k,j);
	coeff(lam,k,j)=(long)p1;
      }
      if(signe((GEN)fl[k]))
      {
	for(i=k+1;i<=n;i++)
	{
	  bb=(GEN)coeff(lam,i,k);
	  coeff(lam,i,k)=ldiv(gsub(gmul((GEN)B[k+1],gcoeff(lam,i,k-1)),gmul(la,bb)),(GEN)B[k]);
	  coeff(lam,i,k-1)=ldiv(gadd(gmul(la,gcoeff(lam,i,k-1)),gmul((GEN)B[k-1],bb)),(GEN)B[k]);
	}
	B[k]=ldiv(gadd(p3,p4),(GEN)B[k]);
      }
      else
      {
	if(!gcmp0(la))
	{
	  p2=(GEN)B[k];p1=gdiv(p4,p2);
	  for(i=k+1;i<lx;i++)
	    coeff(lam,i,k-1)=ldiv(gmul(la,gcoeff(lam,i,k-1)),p2);
	  for(j=k+1;j<lx-1;j++)
	    for(i=j+1;i<lx;i++)
	      coeff(lam,i,j)=ldiv(gmul(p1,gcoeff(lam,i,j)),p2);
	  B[k+1]=B[k]=(long)p1;
	  for(i=k+2;i<=lx;i++)
	    B[i]=ldiv(gmul(p1,(GEN)B[i]),p2);
	}
	else
	{
	  coeff(lam,k,k-1)=zero;
	  for(i=k+1;i<lx;i++)
	  {coeff(lam,i,k)=coeff(lam,i,k-1);coeff(lam,i,k-1)=zero;}
	  B[k]=B[k-1];fl[k]=un;fl[k-1]=zero;
	}
      }
      if(k>2) k--;
    }
    else
    {
      for(l=k-2;l>=1;l--)
      {
	u=(GEN)coeff(lam,k,l);
	if(pslg(u)>=pslg((GEN)B[l+1]))
	{
	  q=gdeuc(u,(GEN)B[l+1]);cq=gdivsg(1,content(q));q=gmul(q,cq);flc=1;
	  h[k]=lsub(gmul(cq,(GEN)h[k]),gmul(q,(GEN)h[l]));
	  coeff(lam,k,l)=lsub(gmul(cq,gcoeff(lam,k,l)),gmul(q,(GEN)B[l+1]));
	  for(i=1;i<l;i++) coeff(lam,k,i)=lsub(gmul(cq,gcoeff(lam,k,i)),gmul(q,gcoeff(lam,l,i)));
	}
      }
      k++;
      if(k>n) 
      {
	for(k=1;(k<=n)&&(!signe((GEN)fl[k]));k++);
	tetpil=avma;
	if(!all)
	{
	  y=cgetg(3,17);
	  p2=cgetg(k,19);for(i=1;i<k;i++) p2[i]=lcopy((GEN)h[i]);
	  y[1]=(long)p2;p2=cgetg(n-k+2,19);y[2]=(long)p2;
	  for(i=k;i<=n;i++) p2[i-k+1]=lcopy((GEN)h[i]);
	}
	else
	{
	  if(all==1)
	  {
	    y=cgetg(k,19);for(i=1;i<k;i++) y[i]=lcopy((GEN)h[i]);
	  }
	  else
	  {
	    y=cgetg(n-k+2,19);
	    for(i=k;i<=n;i++) y[i-k+1]=lcopy((GEN)h[i]);
	  }
	}
	return gerepile(av,tetpil,y);
      }
    }
    if(avma<lim)
    {
      tetpil=avma;
      B=gcopy(B);h=gcopy(h);lam=gcopy(lam);fl=gcopy(fl);
      dec=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;
      B+=dec;h+=dec;lam+=dec;fl+=dec;
    }
  }
}

/********************************************************************/
/********************************************************************/
/**                                                                **/
/**                   FONCTIONS DE RECTPLOT                        **/
/**                                                                **/
/********************************************************************/
/********************************************************************/

GEN
initrect(long ne, long x, long y)
{
  Rect *e;
  RectObj *z;

  if(!GOODRECT(ne)) err(rploter2);
  if((x<=1)||(y<=1)) err(rploter1);
  e=rectgraph[ne];
  if(!(z = (RectObj*)malloc(sizeof(RectObj)))) err(memer);
  RoNext(z) = 0;
  RoType(z) = ROt_NULL;
  RHead(e)=RTail(e)=z;
  RXsize(e)=x;RYsize(e)=y;
  affsr(0,RXcursor(e)); affsr(0,RYcursor(e));
  affsr(1,RXscale(e)); affsr(0,RXshift(e));
  affsr(1,RYscale(e)); affsr(0,RYshift(e));
  return gnil;
}

GEN
rectscale(long ne, GEN x1, GEN x2, GEN y1, GEN y2)
{
  long av;
  Rect *e;
  GEN p1,p2,p3;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];av=avma;
  p2=gadd(RXshift(e),gmul(RXscale(e),RXcursor(e)));p3=gadd(RYshift(e),gmul(RYscale(e),RYcursor(e)));
  gaffect(p1=gdivsg(RXsize(e),gsub(x2,x1)),RXscale(e));gaffect(gneg(gmul(x1,p1)),RXshift(e));
  gaffect(p1=gdivsg(RYsize(e),gsub(y1,y2)),RYscale(e));gaffect(gneg(gmul(y2,p1)),RYshift(e));
  gaffect(gdiv(gsub(p2,RXshift(e)),RXscale(e)),RXcursor(e));gaffect(gdiv(gsub(p3,RYshift(e)),RYscale(e)),RYcursor(e));
  avma=av;return gnil;
}

GEN
rectcursor(long ne)
{
  GEN z;
  Rect *e;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  z=cgetg(3,17);z[1]=lcopy(RXcursor(e));z[2]=lcopy(RYcursor(e));
  return z;
}

GEN
rectmove(long ne, GEN x, GEN y) /* code = 0 */
{
  long av;
  Rect *e;
  RectObj *z;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  if(!(z = (RectObj*)malloc(sizeof(RectObj1P)))) err(memer);
  av=avma;gaffect(x,RXcursor(e));gaffect(y,RYcursor(e));
  RoNext(z) = 0; RoType(z) = ROt_MV;
  RoMVx(z) = gtolong(gadd(gmul(x,RXscale(e)),RXshift(e)));
  RoMVy(z) = gtolong(gadd(gmul(y,RYscale(e)),RYshift(e)));
  avma=av;
  if(!RHead(e)) RHead(e)=RTail(e)=z;
  else {RoNext(RTail(e))=z;RTail(e)=z;}
  return gnil;
}

GEN
rectrmove(long ne, GEN x, GEN y) /* code = 0 */
{
  long av;
  Rect *e;
  RectObj *z;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  if(!(z = (RectObj*)malloc(sizeof(RectObj1P)))) err(memer);
  av=avma;gaffect(gadd(x,RXcursor(e)),RXcursor(e));gaffect(gadd(y,RYcursor(e)),RYcursor(e));
  RoNext(z) = 0; RoType(z) = ROt_MV;
  RoMVx(z) = gtolong(gadd(gmul(RXcursor(e),RXscale(e)),RXshift(e)));
  RoMVy(z) = gtolong(gadd(gmul(RYcursor(e),RYscale(e)),RYshift(e)));
  avma=av;
  if(!RHead(e)) RHead(e)=RTail(e)=z;
  else {RoNext(RTail(e))=z;RTail(e)=z;}
  return gnil;
}

GEN
rectpoint(long ne, GEN x, GEN y) /* code = 1 */
{
  long av;
  Rect *e;
  RectObj *z;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  if(!(z = (RectObj*)malloc(sizeof(RectObj1P)))) err(memer);
  av=avma;gaffect(x,RXcursor(e));gaffect(y,RYcursor(e));
  RoNext(z)=0;
  RoPTx(z) = gtolong(gadd(gmul(x,RXscale(e)),RXshift(e)));
  RoPTy(z) = gtolong(gadd(gmul(y,RYscale(e)),RYshift(e)));
  avma=av;
  RoType(z) = ((RoPTx(z)<0)||(RoPTy(z)<0)||(RoPTx(z)>RXsize(e))
	       ||(RoPTy(z)>RYsize(e)))
    ? ROt_MV
      : ROt_PT;
  if(!RHead(e)) RHead(e)=RTail(e)=z;
  else {RoNext(RTail(e))=z;RTail(e)=z;}
  return gnil;
}

GEN
rectrpoint(long ne, GEN x, GEN y) /* code = 1 */
{
  long av;
  Rect *e;
  RectObj *z;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  if(!(z = (RectObj*)malloc(sizeof(RectObj1P)))) err(memer);
  av=avma;gaffect(gadd(x,RXcursor(e)),RXcursor(e));gaffect(gadd(y,RYcursor(e)),RYcursor(e));
  RoNext(z)=0;
  RoPTx(z) = gtolong(gadd(gmul(RXcursor(e),RXscale(e)),RXshift(e)));
  RoPTy(z) = gtolong(gadd(gmul(RYcursor(e),RYscale(e)),RYshift(e)));
  avma=av;
  RoType(z) = ((RoPTx(z)<0)||(RoPTy(z)<0)||(RoPTx(z)>RXsize(e))
	       ||(RoPTy(z)>RYsize(e)))
    ? ROt_MV
      : ROt_PT;
  if(!RHead(e)) RHead(e)=RTail(e)=z;
  else {RoNext(RTail(e))=z;RTail(e)=z;}
  return gnil;
}

GEN
rectline(long ne, GEN gx2, GEN gy2) /* code = 2 */
{
  long dx,dy,dxy,xmin,xmax,ymin,ymax,x1,y1,x2,y2,av;
  Rect *e;
  RectObj *z;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  if(!(z = (RectObj*)malloc(sizeof(RectObj2P)))) err(memer);
  av=avma;
  x1=gtolong(gadd(gmul(RXcursor(e),RXscale(e)),RXshift(e)));y1=gtolong(gadd(gmul(RYcursor(e),RYscale(e)),RYshift(e)));
  gaffect(gx2,RXcursor(e));gaffect(gy2,RYcursor(e));
  x2=gtolong(gadd(gmul(gx2,RXscale(e)),RXshift(e)));y2=gtolong(gadd(gmul(gy2,RYscale(e)),RYshift(e)));
  avma=av;
  xmin=max(min(x1,x2),0);xmax=min(max(x1,x2),RXsize(e));
  ymin=max(min(y1,y2),0);ymax=min(max(y1,y2),RYsize(e));
  dxy=x1*y2-y1*x2;dx=x2-x1;dy=y2-y1;
  if(dy)
  {
    if(dx*dy<0) {xmin=max(xmin,(dxy+RYsize(e)*dx)/dy);xmax=min(xmax,dxy/dy);}
    else {xmin=max(xmin,dxy/dy);xmax=min(xmax,(dxy+RYsize(e)*dx)/dy);}
  }
  if(dx)
  {
    if(dx*dy<0) {ymin=max(ymin,(RXsize(e)*dy-dxy)/dx);ymax=min(ymax,-dxy/dx);}
    else {ymin=max(ymin,-dxy/dx);ymax=min(ymax,(RXsize(e)*dy-dxy)/dx);}
  }
  RoNext(z)=0;
  RoLNx1(z) = xmin; RoLNx2(z) = xmax;
  if(dx*dy<0) {RoLNy1(z) = ymax; RoLNy2(z) = ymin;}
  else {RoLNy1(z) = ymin; RoLNy2(z) = ymax;}
  RoType(z) = ((xmin>xmax)||(ymin>ymax))
    ? ROt_MV
      : ROt_LN;
  if(!RHead(e)) RHead(e)=RTail(e)=z;
  else {RoNext(RTail(e))=z;RTail(e)=z;}
  return gnil;
}

GEN
rectrline(long ne, GEN gx2, GEN gy2) /* code = 2 */
{
  long x1,y1,x2,y2,dx,dy,dxy,xmin,xmax,ymin,ymax,av;
  Rect *e;
  RectObj *z;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  if(!(z = (RectObj*)malloc(sizeof(RectObj2P)))) err(memer);
  av=avma;
  x1=gtolong(gadd(gmul(RXcursor(e),RXscale(e)),RXshift(e)));y1=gtolong(gadd(gmul(RYcursor(e),RYscale(e)),RYshift(e)));
  gaffect(gadd(gx2,RXcursor(e)),RXcursor(e));gaffect(gadd(gy2,RYcursor(e)),RYcursor(e));
  x2=gtolong(gadd(gmul(RXcursor(e),RXscale(e)),RXshift(e)));y2=gtolong(gadd(gmul(RYcursor(e),RYscale(e)),RYshift(e)));
  avma=av;
  xmin=max(min(x1,x2),0);xmax=min(max(x1,x2),RXsize(e));
  ymin=max(min(y1,y2),0);ymax=min(max(y1,y2),RYsize(e));
  dxy=x1*y2-y1*x2;dx=x2-x1;dy=y2-y1;
  if(dy)
  {
    if(dx*dy<0) {xmin=max(xmin,(dxy+RYsize(e)*dx)/dy);xmax=min(xmax,dxy/dy);}
    else {xmin=max(xmin,dxy/dy);xmax=min(xmax,(dxy+RYsize(e)*dx)/dy);}
  }
  if(dx)
  {
    if(dx*dy<0) {ymin=max(ymin,(RXsize(e)*dy-dxy)/dx);ymax=min(ymax,-dxy/dx);}
    else {ymin=max(ymin,-dxy/dx);ymax=min(ymax,(RXsize(e)*dy-dxy)/dx);}
  }
  RoNext(z)=0;
  RoLNx1(z) = xmin; RoLNx2(z) = xmax;
  if(dx*dy<0) {RoLNy1(z) = ymax; RoLNy2(z) = ymin;}
  else {RoLNy1(z) = ymin; RoLNy2(z) = ymax;}
  RoType(z) = ((xmin>xmax)||(ymin>ymax))
    ? ROt_MV
      : ROt_LN;
  if(!RHead(e)) RHead(e)=RTail(e)=z;
  else {RoNext(RTail(e))=z;RTail(e)=z;}
  return gnil;
}


GEN
rectbox(long ne, GEN gx2, GEN gy2) /* code = 3 */
{
  long av;
  Rect *e;
  long x1,y1,x2,y2,xmin,ymin,xmax,ymax;
  RectObj *z;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  if(!(z = (RectObj*)malloc(sizeof(RectObj2P)))) err(memer);
  av=avma;
  x1=gtolong(gadd(gmul(RXcursor(e),RXscale(e)),RXshift(e)));y1=gtolong(gadd(gmul(RYcursor(e),RYscale(e)),RYshift(e)));
  x2=gtolong(gadd(gmul(gx2,RXscale(e)),RXshift(e)));y2=gtolong(gadd(gmul(gy2,RYscale(e)),RYshift(e)));
  avma=av;
  xmin=max(min(x1,x2),0);xmax=min(max(x1,x2),RXsize(e));
  ymin=max(min(y1,y2),0);ymax=min(max(y1,y2),RYsize(e));
  RoNext(z)=0; RoType(z) = ROt_BX;
  RoBXx1(z) = xmin; RoBXy1(z) = ymin;
  RoBXx2(z) = xmax; RoBXy2(z) = ymax;
  if(!RHead(e)) RHead(e)=RTail(e)=z;
  else {RoNext(RTail(e))=z;RTail(e)=z;}
  return gnil;
}

GEN
rectrbox(long ne, GEN gx2, GEN gy2) /* code = 3 */
{
  long av;
  Rect *e;
  RectObj *z;
  long x1,y1,x2,y2,xmin,ymin,xmax,ymax;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  if(!(z = (RectObj*)malloc(sizeof(RectObj2P)))) err(memer);
  av=avma;
  x1=gtolong(gadd(gmul(RXcursor(e),RXscale(e)),RXshift(e)));y1=gtolong(gadd(gmul(RYcursor(e),RYscale(e)),RYshift(e)));
  x2=gtolong(gadd(gmul(gadd(gx2,RXcursor(e)),RXscale(e)),RXshift(e)));
  y2=gtolong(gadd(gmul(gadd(gy2,RYcursor(e)),RYscale(e)),RYshift(e)));
  avma=av;
  xmin=max(min(x1,x2),0);xmax=min(max(x1,x2),RXsize(e));
  ymin=max(min(y1,y2),0);ymax=min(max(y1,y2),RYsize(e));
  RoNext(z) = 0; RoType(z) = ROt_BX;
  RoBXx1(z) = xmin;RoBXy1(z) = ymin;
  RoBXx2(z) = xmax;RoBXy2(z) = ymax;
  if(!RHead(e)) RHead(e)=RTail(e)=z;
  else {RoNext(RTail(e))=z;RTail(e)=z;}
  return gnil;
}

GEN
killrect(long ne)
{
  RectObj *p1,*p2;
  Rect *e;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  p1=RHead(e);
  RHead(e) = RTail(e) = NULL;
  RXsize(e) = RYsize(e)=0;
  affsr(0,RXcursor(e));affsr(0,RYcursor(e));affsr(1,RXscale(e));affsr(0,RXshift(e));affsr(1,RYscale(e));affsr(0,RYshift(e));
  while(p1) 
  {
    if((RoType(p1)==ROt_MP)||(RoType(p1)==ROt_ML)) {
      free((long *)RoMPxs(p1));
      free((long *)RoMPys(p1));
    }
    if((RoType(p1)==ROt_ST)) free((long *)RoSTs(p1));
    p2=RoNext(p1);free(p1);p1=p2;
  }
  return gnil;
}

GEN
rectpoints(long ne, GEN listx, GEN listy) /* code = 4 */
{
  long lx,*px,*py,*ptx,*pty,x,y,i,cp,av;
  Rect *e;
  RectObj *z;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  if((typ(listx)<17)||(typ(listx)>18)||(typ(listy)<17)||(typ(listy)>18))
    err(ploter4);
  lx=lg(listx);if(lg(listy)!=lx) err(ploter5);
  lx--;if(!lx) return gnil;
  if(!(px=(long*)malloc(lx<<TWOPOTBYTES_IN_LONG))) err(memer);
  if(!(py=(long*)malloc(lx<<TWOPOTBYTES_IN_LONG))) err(memer);
  cp=0;
  for(i=0;i<lx;i++) 
  {
    av=avma;
    x=px[i]=gtolong(gadd(gmul(RXscale(e),(GEN)listx[i+1]),RXshift(e)));
    y=py[i]=gtolong(gadd(gmul(RYscale(e),(GEN)listy[i+1]),RYshift(e)));
    if((x>=0)&&(y>=0)&&(x<=RXsize(e))&&(y<=RYsize(e))) cp++;
    avma=av;
  }
  if(!cp) {free(px);free(py);return gnil;}
  if(!(ptx=(long*)malloc(cp<<TWOPOTBYTES_IN_LONG))) err(memer);
  if(!(pty=(long*)malloc(cp<<TWOPOTBYTES_IN_LONG))) err(memer);
  cp=0;
  for(i=0;i<lx;i++) 
  {
    x=px[i];y=py[i];
    if((x>=0)&&(y>=0)&&(x<=RXsize(e))&&(y<=RYsize(e))) {ptx[cp]=x;pty[cp]=y;}
    cp++;
  }
  free(px);free(py);
  if(!(z = (RectObj*)malloc(sizeof(RectObjMP)))) err(memer);
  RoNext(z)=0; RoType(z) = ROt_MP; 
  RoMPcnt(z) = cp; RoMPxs(z) = ptx; RoMPys(z) = pty;
  if(!RHead(e)) RHead(e)=RTail(e)=z;
  else {RoNext(RTail(e))=z;RTail(e)=z;}
  return gnil;
}

GEN
rectlines(long ne, GEN listx, GEN listy) /* code = 5 */
{
  long lx,*ptx,*pty,i,av;
  Rect *e;
  RectObj *z;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  if((typ(listx)<17)||(typ(listx)>18)||(typ(listy)<17)||(typ(listy)>18))
    err(ploter4);
  lx=lg(listx);if(lg(listy)!=lx) err(ploter5);
  lx--;if(!lx) return gnil;
  if(!(ptx=(long*)malloc(lx<<TWOPOTBYTES_IN_LONG))) err(memer);
  if(!(pty=(long*)malloc(lx<<TWOPOTBYTES_IN_LONG))) err(memer);
  for(i=0;i<lx;i++) 
  {
    av=avma;
    ptx[i]=gtolong(gadd(gmul(RXscale(e),(GEN)listx[i+1]),RXshift(e)));
    pty[i]=gtolong(gadd(gmul(RYscale(e),(GEN)listy[i+1]),RYshift(e)));
    avma=av;
  }
  if(!(z=(RectObj*)malloc(sizeof(RectObj2P)))) err(memer);
  RoNext(z)=0;
  RoType(z)=ROt_ML;
  RoMLcnt(z)=lx;
  RoMLxs(z)=ptx;
  RoMLys(z)=pty;
  if(!RHead(e)) RHead(e)=RTail(e)=z;
  else {
    RoNext(RTail(e))=z;
    RTail(e)=z;
  }
  return gnil;
}

GEN
rectstring(long ne, GEN x) /* code = 6 */
{
  long i,a,lx,av;
  char *c;
  GEN p1;
  Rect *e;
  RectObj *z;

  if(typ(x)<17)
  {
    if(!(c=(char*)malloc(20))) err(memer);
    sprintf(c,"%9.3f",gtodouble(x));
    lx = strlen(c) + 2;
  }
  else
  {
    lx=lg(x);if(!(c=(char*)malloc(lx))) err(memer);
    for(i=1;i<lx;i++)
    {
      p1=(GEN)x[i];if(typ(p1)!=1) err(rploter6);
      a=itos(p1);if((a<0)||(a>255)) err(rploter6);
      c[i-1]=(char)a;
    }
    c[lx-1]=(char)0;
  }
  return rectmstring(ne,c,lx - 2);
}

/* Put malloced string: */

static GEN
rectmstring(long ne, char *s, int l) /* code = 6 */
{
  long i,a,lx,av;
  GEN p1;
  Rect *e;
  RectObj *z;
  char *s1;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  if(!(z = (RectObj*)malloc(sizeof(RectObjST)))) err(memer);
  RoNext(z)=0; RoType(z) = ROt_ST;
  RoSTl(z) = l; RoSTs(z) = s;
  av=avma;
  RoSTx(z) = gtolong(gadd(gmul(RXscale(e),RXcursor(e)),RXshift(e)));
  RoSTy(z) = gtolong(gadd(gmul(RYscale(e),RYcursor(e)),RYshift(e)));
  avma=av;
  if(!RHead(e)) RHead(e)=RTail(e)=z;
  else {RoNext(RTail(e))=z;RTail(e)=z;}
  return gnil;
}

GEN
rectcstring(long ne, char *s) /* code = 6 */
{
  long i,l;
  char *s1;

  l = strlen(s);
  s1 = (char *)malloc(l + 1);
  strcpy(s1,s);
  return rectmstring(ne, s1, l);
}

GEN
rectpointtype(long ne, long type) /* code = 0 */
{
  long av;
  Rect *e;
  RectObj *z;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  if(!(z = (RectObj*)malloc(sizeof(RectObjPN)))) err(memer);
  RoNext(z) = 0; RoType(z) = ROt_PTT;
  RoPTTpen(z) = type;
  if(!RHead(e)) RHead(e)=RTail(e)=z;
  else {RoNext(RTail(e))=z;RTail(e)=z;}
  return gnil;
}

GEN
rectlinetype(long ne, long type) 
{
  long av;
  Rect *e;
  RectObj *z;

  if(!GOODRECT(ne)) err(rploter2);
  e=rectgraph[ne];
  if(!(z = (RectObj*)malloc(sizeof(RectObjPN)))) err(memer);
  RoNext(z) = 0; RoType(z) = ROt_LNT;
  RoLNTpen(z) = type;
  if(!RHead(e)) RHead(e)=RTail(e)=z;
  else {RoNext(RTail(e))=z;RTail(e)=z;}
  return gnil;
}

GEN
rectcopy(long source, long dest, long xoff, long yoff)
{
  Rect *s = rectgraph[source], *d = rectgraph[dest];
  RectObj *p1 = RHead(s);
  RectObj *tail = RTail(d);
  RectObj *next;
  int i;

  while (p1) {
      switch(RoType(p1))
      {
      case ROt_PT:
	if(!(next = (RectObj*)malloc(sizeof(RectObj1P)))) err(memer);
	memcpy(next,p1,sizeof(RectObj1P));
	RoPTx(next) += xoff;
	RoPTy(next) += yoff;
	RoNext(tail) = next;
	tail = next;
	break;
      case ROt_LN:
      case ROt_BX:
	if(!(next = (RectObj*)malloc(sizeof(RectObj2P)))) err(memer);
	memcpy(next,p1,sizeof(RectObj2P));
	RoLNx1(next) += xoff;
	RoLNy1(next) += yoff;
	RoLNx2(next) += xoff;
	RoLNy2(next) += yoff;
	RoNext(tail) = next;
	tail = next;
	break;
      case ROt_MP:
      case ROt_ML:
	if(!(next = (RectObj*)malloc(sizeof(RectObjMP)))) err(memer);
	memcpy(next,p1,sizeof(RectObjMP));
	if(!(RoMPxs(next) = (long*)malloc(sizeof(long)*RoMPcnt(next)))) 
	  err(memer);
	if(!(RoMPys(next) = (long*)malloc(sizeof(long)*RoMPcnt(next)))) 
	  err(memer);
	memcpy(RoMPxs(next),RoMPxs(p1),sizeof(long)*RoMPcnt(next));
	memcpy(RoMPys(next),RoMPys(p1),sizeof(long)*RoMPcnt(next));
	for (i=0;i<RoMPcnt(next);i++) {
	  RoMPxs(next)[i] += xoff;
	  RoMPys(next)[i] += yoff;
	}
	RoNext(tail) = next;
	tail = next;
	break;
      case ROt_ST:
	if(!(next = (RectObj*)malloc(sizeof(RectObjMP)))) err(memer);
	memcpy(next,p1,sizeof(RectObjMP));
	if(!(RoSTs(next) = (char*)malloc(RoSTl(p1)+1))) err(memer);
	memcpy(RoSTs(next),RoSTs(p1),RoSTl(p1)+1);
	RoSTx(next) += xoff;
	RoSTy(next) += yoff;
	RoNext(tail) = next;
	tail = next;
	break;
      case ROt_PTT:
      case ROt_LNT:
	if(!(next = (RectObj*)malloc(sizeof(RectObjPN)))) err(memer);
	memcpy(next,p1,sizeof(RectObjPN));
	RoNext(tail) = next;
	tail = next;
	break;
      default: break;
      }
    p1=RoNext(p1);
  }
  RoNext(tail) = NULL;
  RTail(d) = tail;
  return gnil;
}

/*************************************************************************/
/*                                                                       */
/*                                                                       */
/*                         POSTSCRIPT OUTPUT                             */
/*                                                                       */
/*                                                                       */
/*************************************************************************/

typedef struct spoint {
  int x,y;} SPoint; 
typedef struct ssegment {
  int x1,y1,x2,y2;} SSegment;
typedef struct srectangle {
  int x,y,width,height;} SRectangle;

void ps_point(FILE *psfile, int x, int y),ps_line(FILE *psfile, int x1, int y1, int x2, int y2),ps_rect(FILE *psfile, int x1, int y1, int x2, int y2),ps_string(FILE *psfile, int x, int y, char *c);

#undef ISCR
#undef JSCR
#define ISCR 1120 /* 1400 en haute resolution */     
#define JSCR 800  /* 1120 en haute resolution */     
#define DECI 95 
#define DECJ 30 

GEN
postploth(entree *ep, GEN a, GEN b, char *ch)
{
  long av,av2,jz,j1,i,sig,is,is2,js,js2;
  GEN p1,ysml,ybig,x,diff,dyj,dx,y[ISCR+1];
  char c1[20];
  FILE *psfile;

  is=ISCR-DECI;js=JSCR-DECJ;is2=is-DECI;js2=js-DECJ;
  sig=gcmp(b,a); if(!sig) return gnil;
  psfile = fopen("pari.ps", "a");
  if (!psfile) err(poster1);
  fprintf(psfile,"%%!\n50 50 translate\n/Times-Roman findfont 16 scalefont setfont\n0.65 0.65 scale\n");
  ps_line(psfile,DECI,DECJ,DECI,js);
  ps_line(psfile,DECI,DECJ,is,DECJ);
  ps_line(psfile,is,DECJ,is,js);
  ps_line(psfile,DECI,js,is,js);

  av=avma;
  if(sig<0) {x=a;a=b;b=x;}
  for(i=1;i<=is2;i++) y[i]=cgetr(3);
  newvalue(ep,cgetr(3)); x=(GEN)ep->value; gaffect(a,x);
  dx=gdivgs(gsub(b,a),is2-1);ysml=gzero;ybig=gzero;
  av2=avma;
  for(i=1;i<=is2;i++)
  {
    gaffect(lisexpr(ch),y[i]);
    if(gcmp(y[i],ysml)<0) ysml=y[i];
    if(gcmp(y[i],ybig)>0) ybig=y[i];
    gaddz(x,dx,x);avma=av2;
  }
  diff=gsub(ybig,ysml);
  if(gcmp0(diff)) {ybig=gaddsg(1,ybig);diff=gun;}
  dyj=gdivsg(js2-1,diff);jz=js+itos(ground(gmul(ysml,dyj)));
  ps_line(psfile,DECI,jz,is,jz);
  if(gsigne(a)*gsigne(b)<0)
  {
    jz=1-itos(ground(gdiv(a,dx)))+DECI;
    ps_line(psfile,jz,DECJ,jz,js);
  }
  av2=avma;
  for(i=1;i<=is2;i++)
  {
    j1=js-itos(ground(gmul(gsub((GEN)y[i],ysml),dyj)));
    if(i==1) fprintf(psfile,"%ld %d moveto\n",j1,DECI);
    else fprintf(psfile,"%ld %ld lineto\n",j1,i-1+DECI);
    avma=av2;
  }
  p1=cgetr(4);gaffect(ysml,p1);sprintf(c1," %9.3f ",rtodbl(p1));
  ps_string(psfile,5,js,c1);gaffect(ybig,p1);
  sprintf(c1," %9.3f ",rtodbl(p1));
  ps_string(psfile,5,DECJ,c1);gaffect(a,p1);
  sprintf(c1," %9.3f ",rtodbl(p1));
  ps_string(psfile,DECI-36,js+20,c1);gaffect(b,p1);
  sprintf(c1," %9.3f ",rtodbl(p1));
  ps_string(psfile,is-36,js+20,c1);
  fprintf(psfile,"stroke showpage\n");fclose(psfile);
  avma = av;killvalue(ep);return gnil;
}

GEN
postploth2(entree *ep, GEN a, GEN b, char *ch)
{
  long av,av2,jz,iz,k1,j1,i,sig,is,is2,js,js2;
  GEN p1,ysml,ybig,xsml,xbig,diffx,diffy,dxj,t,dyj,dt,y[ISCR+1],x[ISCR+1];
  char c1[20];
  FILE *psfile;

  is=ISCR-DECI;js=JSCR-DECJ;is2=is-DECI;js2=js-DECJ;
  sig=gcmp(b,a); if(!sig) return gnil;
  psfile = fopen("pari.ps", "a");
  if (!psfile) err(poster1);
  fprintf(psfile,"%%!\n50 50 translate\n/Times-Roman findfont 16 scalefont setfont\n0.65 0.65 scale\n");
  ps_line(psfile,DECI,DECJ,DECI,js);
  ps_line(psfile,DECI,DECJ,is,DECJ);
  ps_line(psfile,is,DECJ,is,js);
  ps_line(psfile,DECI,js,is,js);

  av=avma;
  if(sig<0) {p1=a;a=b;b=p1;}
  for(i=1;i<=is2;i++) {x[i]=cgetr(3);y[i]=cgetr(3);}
  newvalue(ep,cgetr(3)); t=(GEN)ep->value; gaffect(a,t);
  dt=gdivgs(gsub(b,a),is2-1);ysml=ybig=xsml=xbig=gzero;
  av2=avma;
  for(i=1;i<=is2;i++)
  {
    p1=lisexpr(ch);gaffect((GEN)p1[1],x[i]);gaffect((GEN)p1[2],y[i]);
    if(gcmp(y[i],ysml)<0) ysml=y[i];
    if(gcmp(y[i],ybig)>0) ybig=y[i];
    if(gcmp(x[i],xsml)<0) xsml=x[i];
    if(gcmp(x[i],xbig)>0) xbig=x[i];
    gaddz(t,dt,t);avma=av2;
  }
  diffy=gsub(ybig,ysml);
  if(gcmp0(diffy)) {ybig=gaddsg(1,ybig);diffy=gun;}
  diffx=gsub(xbig,xsml);
  if(gcmp0(diffx)) {xbig=gaddsg(1,xbig);diffx=gun;}
  dyj=gdivsg(js2-1,diffy);jz=js+itos(ground(gmul(ysml,dyj)));
  dxj=gdivsg(is2-1,diffx);iz=DECI-itos(ground(gmul(xsml,dxj)));
  if(gsigne(ysml)*gsigne(ybig)<0)
    ps_line(psfile,DECI,jz,is,jz);
  if(gsigne(xsml)*gsigne(xbig)<0)
    ps_line(psfile,iz,DECJ,iz,js);
  av2=avma;
  for(i=1;i<=is2;i++)
  {
    j1=js-itos(ground(gmul(gsub(y[i],ysml),dyj)));
    k1=DECI+itos(ground(gmul(gsub(x[i],xsml),dxj)));
    if(i==1) fprintf(psfile,"%ld %ld moveto\n",j1,k1);
    else fprintf(psfile,"%ld %ld lineto\n",j1,k1);
    avma=av2;
  }
  p1=cgetr(4);gaffect(ysml,p1);sprintf(c1," %9.3f ",rtodbl(p1));
  ps_string(psfile,5,js,c1);gaffect(ybig,p1);
  sprintf(c1," %9.3f ",rtodbl(p1));
  ps_string(psfile,5,DECJ,c1);gaffect(xsml,p1);
  sprintf(c1," %9.3f ",rtodbl(p1));
  ps_string(psfile,DECI-36,js+20,c1);gaffect(xbig,p1);
  sprintf(c1," %9.3f ",rtodbl(p1));
  ps_string(psfile,is-36,js+20,c1);
  fprintf(psfile,"stroke showpage\n");fclose(psfile);
  avma = av;killvalue(ep);return gnil;
}

GEN
postplothraw(GEN listx, GEN listy)
{
  long av = avma,av2,i,lx,is,js,is2,js2;
  char c1[20];
  SPoint *points;
  GEN p1,xsml,xbig,ysml,ybig,dx,dy,scal,scaly;
  FILE *psfile;

  if((typ(listx)<17)||(typ(listx)>18)||(typ(listy)<17)||(typ(listy)>18))
    err(ploter4);
  lx=lg(listx);
  if(lg(listy)!=lx) err(ploter5);
  if(lx==1) return gnil;
  points = (SPoint*)malloc(lx*sizeof(SPoint));
  if(!points) err(ploter6);
  is=ISCR-DECI;js=JSCR-DECJ;is2=is-DECI;js2=js-DECJ;
  av=avma;xsml=xbig=(GEN)listx[1];ysml=ybig=(GEN)listy[1];
  for(i = 0; i < lx-1; i++)
  {
    p1=(GEN)listx[i+1];
    if(gcmp(p1,xsml)<0) xsml=p1;if(gcmp(p1,xbig)>0) xbig=p1;
    p1=(GEN)listy[i+1];
    if(gcmp(p1,ysml)<0) ysml=p1;if(gcmp(p1,ybig)>0) ybig=p1;
  }
  dx=gsub(xbig,xsml);dy=gsub(ybig,ysml);
  if(gcmp0(dx))
  {
    if(gcmp0(dy)) 
    {
      scal=gun;dx=gsubsg(is2>>1,xsml);
      dy=gsubsg(js2>>1,ysml);
    }
    else
    {
      scal=gdivsg(js2,dy);
      dx=gneg(gmul(scal,xsml));dy=gneg(gmul(scal,ysml));
    }
  }
  else
  {
    scal=gdivsg(is2,dx);
    if(!gcmp0(dy))
    {
      scaly=gdivsg(js2,dy);if(gcmp(scaly,scal)<0) scal=scaly;
    }
    dx=gneg(gmul(scal,xsml));dy=gneg(gmul(scal,ysml));
  }
  for(i = 0; i < lx-1; i++)
  {
    av2=avma;
    points[i].x = DECI + itos(ground(gadd(gmul((GEN)listx[i+1],scal),dx)));
    points[i].y = js - itos(ground(gadd(gmul((GEN)listy[i+1],scal),dy)));
    avma=av2;
  }
  xsml=gneg(gdiv(dx,scal));xbig=gdiv(gsubsg(is2,dx),scal);
  ysml=gneg(gdiv(dy,scal));ybig=gdiv(gsubsg(js2,dy),scal);

  psfile = fopen("pari.ps", "a");
  if (!psfile) err(poster1);
  fprintf(psfile,"%%!\n50 50 translate\n/Times-Roman findfont 16 scalefont setfont\n0.65 0.65 scale\n");
  ps_line(psfile,DECI,DECJ,DECI,js);
  ps_line(psfile,DECI,DECJ,is,DECJ);
  ps_line(psfile,is,DECJ,is,js);
  ps_line(psfile,DECI,js,is,js);
  for(i = 0; i < lx-1; i++)
    ps_point(psfile,points[i].x,points[i].y);
  free(points);avma = av;
  p1=cgetr(4);gaffect(ysml,p1);sprintf(c1," %9.3f ",rtodbl(p1));
  ps_string(psfile,5,js,c1);gaffect(ybig,p1);
  sprintf(c1," %9.3f ",rtodbl(p1));
  ps_string(psfile,5,DECJ,c1);gaffect(xsml,p1);
  sprintf(c1," %9.3f ",rtodbl(p1));
  ps_string(psfile,DECI-36,js+20,c1);gaffect(xbig,p1);
  sprintf(c1," %9.3f ",rtodbl(p1));
  ps_string(psfile,is-36,js+20,c1);
  fprintf(psfile,"stroke showpage\n");fclose(psfile);
  return gnil;
}

GEN
postdraw(GEN list)
{
  long *p1,*ptx,*pty,*numpoints,*numtexts,*xtexts,*ytexts;
  Rect *e;
  long n,i,j,x0,y0,av=avma;
  long a,b,c,d,nd[10],ne;
  char **texts;
  FILE *psfile;

  SPoint *points, **lines, *SLine;
  SSegment *segments; 
  SRectangle *rectangles, SRec;

  if(typ(list)!=17) err(rploter3);
  n=lg(list)-1;if(n%3) err(rploter4);
  n=n/3;if(!n) return gnil;
  nd[0]=nd[1]=nd[2]=nd[3]=nd[4]=nd[5]=nd[6]=0;
  for(i=0;i<n;i++)
  {
    if(typ((GEN)list[3*i+1])!=1) err(rploter5);
    ne=itos((GEN)list[3*i+1]);if(!GOODRECT(ne)) err(rploter2);
    e=rectgraph[ne];
    p1=(long*)RHead(e);while((long)p1) 
    {
      if(p1[1]!=4) nd[p1[1]]++;
      else nd[1]+=p1[2];
      p1=(long*)p1[0];
    }
  }
  points=(SPoint*)malloc(nd[1]*sizeof(SPoint));
  segments=(SSegment*)malloc(nd[2]*sizeof(SSegment));
  rectangles=(SRectangle*)malloc(nd[3]*sizeof(SRectangle));
  lines=(SPoint**)malloc(nd[5]*sizeof(SPoint*));
  numpoints=(long*)malloc(nd[5]*sizeof(long));
  texts=(char**)malloc(nd[6]*sizeof(char*));
  numtexts=(long*)malloc(nd[6]*sizeof(long));
  xtexts=(long*)malloc(nd[6]*sizeof(long));
  ytexts=(long*)malloc(nd[6]*sizeof(long));
  nd[1]=nd[2]=nd[3]=nd[5]=nd[6]=0;
  for(i=0;i<n;i++)
  {
    e=rectgraph[itos((GEN)list[3*i+1])];x0=list[3*i+2];y0=list[3*i+3];
    if((typ((GEN)x0)!=1)||(typ((GEN)y0)!=1)) err(rploter5);
    x0=itos((GEN)x0);y0=itos((GEN)y0);
    p1=(long*)RHead(e);
    while((long)p1)
    {
      switch(p1[1])
      {
	case 1: 
	  points[nd[1]].x=p1[2]+x0;
	  points[nd[1]].y=p1[3]+y0;
	  nd[1]++;break;
	case 2:
	  segments[nd[2]].x1=p1[2]+x0;
	  segments[nd[2]].y1=p1[3]+y0;
	  segments[nd[2]].x2=p1[4]+x0;
	  segments[nd[2]].y2=p1[5]+y0;
	  nd[2]++;break;
	case 3:
	  a=rectangles[nd[3]].x=p1[2]+x0;
	  b=rectangles[nd[3]].y=p1[3]+y0;
	  rectangles[nd[3]].width=p1[4]+x0-a;
	  rectangles[nd[3]].height=p1[5]+y0-b;
	  nd[3]++;break;
	case 4:
	  ptx=(long*)p1[3];pty=(long*)p1[4];
	  for(j=0;j<p1[2];j++)
	  {
	    points[nd[1]+j].x=ptx[j]+x0;
	    points[nd[1]+j].y=pty[j]+y0;
	  }
	  nd[1]+=p1[2];break;
	case 5:
	  ptx=(long*)p1[3];pty=(long*)p1[4];
	  numpoints[nd[5]]=p1[2];
	  lines[nd[5]]=(SPoint*)malloc(p1[2]*sizeof(SPoint));
	  for(j=0;j<p1[2];j++)
	  {
	    lines[nd[5]][j].x=ptx[j]+x0;
	    lines[nd[5]][j].y=pty[j]+y0;
	  }
	  nd[5]++;break;
	case 6: 
	  texts[nd[6]]=(char*)p1[3];numtexts[nd[6]]=p1[2];
	  xtexts[nd[6]]=p1[4]+x0;ytexts[nd[6]]=p1[5]+y0;
	  nd[6]++;break;
	default: break;
      }
      p1=(long*)p1[0];
    }
  }
  psfile = fopen("pari.ps", "a");
  if (!psfile) err(poster1);
  fprintf(psfile,"%%!\n50 50 translate\n/Times-Roman findfont 16 scalefont setfont\n0.65 0.65 scale\n");
  for(i=0;i<nd[1];i++) ps_point(psfile,points[i].x,points[i].y);
  for(i=0;i<nd[2];i++) ps_line(psfile,segments[i].x1,segments[i].y1,segments[i].x2,segments[i].y2);
  for(i=0;i<nd[3];i++) 
  {
    SRec=rectangles[i];a=SRec.x;b=SRec.y;c=a+SRec.width;
    d=b+SRec.height;ps_rect(psfile,a,b,c,d);
  }
  for(i=0;i<nd[5];i++) 
  {
    SLine=lines[i];
    for(j=0;j<numpoints[i];j++)
    {
      if(!j) fprintf(psfile,"%d %d moveto\n",SLine[0].y,SLine[0].x);
      else fprintf(psfile,"%d %d lineto\n",SLine[j].y,SLine[j].x);
    }
  }
  for(i=0;i<nd[6];i++) 
    ps_string(psfile,xtexts[i],ytexts[i],texts[i]);
  fprintf(psfile,"stroke showpage\n");fclose(psfile);
  free(points);free(segments);free(rectangles);
  free(numpoints);for(i=0;i<nd[5];i++) free(lines[i]);
  free(lines);free(texts);free(numtexts);free(xtexts);free(ytexts);
  avma = av;return gnil;
}

void
ps_point(FILE *psfile, int x, int y)
{
  fprintf(psfile,"%d %d moveto\n0 2 rlineto 2 0 rlineto 0 -2 rlineto closepath fill\n",y,x);
  return;
}

void
ps_line(FILE *psfile, int x1, int y1, int x2, int y2)
{
  fprintf(psfile,"%d %d moveto\n%d %d lineto\n",y1,x1,y2,x2);
  return;
}

void
ps_rect(FILE *psfile, int x1, int y1, int x2, int y2)
{
  fprintf(psfile,"%d %d moveto\n%d %d lineto\n%d %d lineto\n%d %d lineto\nclosepath\n",y1,x1,y1,x2,y2,x2,y2,x1);
  return;
}

void
ps_string(FILE *psfile, int x, int y, char *c)
{
  fprintf(psfile,"%d %d moveto 90 rotate\n(",y,x);
  fputs(c,psfile);fprintf(psfile,") show -90 rotate\n");
  return;
}

/***********************************************************************/
/*                                                                     */
/*                                                                     */
/*                          SET OPERATIONS                             */
/*                                                                     */
/*                                                                     */
/***********************************************************************/

GEN
gtoset(GEN x)
{
  long tx=typ(x),av=avma,tetpil,c,i,ly;
  GEN y,z;

  if((tx<17)||(tx==19)) y=gtovec(x);else y=x;
  ly=lg(y);if(ly==1) {avma=av;return cgetg(1,17);}
  y=sort(y);z=cgetg(ly,17);
  c=1;z[1]=y[1];
  for(i=2;i<ly;i++) 
    if(!gegal((GEN)y[i],(GEN)z[c])) z[++c]=y[i];
  tetpil=avma;y=cgetg(c+1,17);for(i=1;i<=c;i++) y[i]=lcopy((GEN)z[i]);
  return gerepile(av,tetpil,y);
}

long
isvecset(GEN x)
{
  long lx,i;

  if(typ(x)!=17) return 0;
  lx=lg(x);
  for(i=1;i<lx-1;i++) if(gcmp((GEN)x[i+1],(GEN)x[i])<=0) return 0;
  return 1;
}

long setsearch(GEN x, GEN y)
/* looks if y belongs to the set x and returns the index if yes, 0 if no */
{
  long lx,j,li,ri,fl;

  if(!isvecset(x)) err(talker,"not a set in setsearch");
  lx=lg(x);li=1;ri=lx-1;
  while(ri>=li)
  {
    j=(ri+li)>>1;
    fl=gcmp((GEN)x[j],y);
    if(!fl) return j;
    if(fl<0) li=j+1;else ri=j-1;
  }
  return 0;
}

GEN
setunion(GEN x, GEN y)
{
  long av=avma,tetpil;
  GEN z;

  if((!isvecset(x))||(!isvecset(y))) err(talker,"not a set in setunion");
  z=concat(x,y);tetpil=avma;return gerepile(av,tetpil,gtoset(z));
}

GEN
setintersect(GEN x, GEN y)
{
  long av=avma,tetpil,i,lx,c;
  GEN z,p1;

  if((!isvecset(x))||(!isvecset(y))) err(talker,"not a set in setintersect");  
  lx=lg(x);z=cgetg(lx,17);c=0;
  for(i=1;i<lx;i++) {p1=(GEN)x[i];if(setsearch(y,p1)) z[++c]=(long)p1;}
  tetpil=avma;p1=cgetg(c+1,17);for(i=1;i<=c;i++) p1[i]=lcopy((GEN)z[i]);
  return gerepile(av,tetpil,p1);
}

GEN
setminus(GEN x, GEN y)
{
  long av=avma,tetpil,i,lx,c;
  GEN z,p1;

  if((!isvecset(x))||(!isvecset(y))) err(talker,"not a set in setintersect");  
  lx=lg(x);z=cgetg(lx,17);c=0;
  for(i=1;i<lx;i++) {p1=(GEN)x[i];if(!setsearch(y,p1)) z[++c]=(long)p1;}
  tetpil=avma;p1=cgetg(c+1,17);for(i=1;i<=c;i++) p1[i]=lcopy((GEN)z[i]);
  return gerepile(av,tetpil,p1);
}


/***********************************************************************/
/*                                                                     */
/*                                                                     */
/*               OPERATIONS ON DIRICHLET SERIES                        */
/*                                                                     */
/*                                                                     */
/***********************************************************************/

/* Addition, subtraction and scalar multiplication of Dirichlet series
   are done on the corresponding vectors */

long
dirval(GEN x)
{
/* a usage interne, pas de verifs */
  long i;
  for(i=1;(i<lg(x))&&gcmp0((GEN)x[i]);i++);
  return i;
}

GEN
dirmul(GEN x, GEN y)
{
  long lx,ly,lz,dx,dy,av,tetpil,i,j;
  GEN z,p1;
  
  if((typ(x)!=17)||(typ(y)!=17)) err(talker,"not a dirseries in dirmul");
  av=avma;dx=dirval(x);dy=dirval(y);lx=lg(x);ly=lg(y);
  if(ly-dy<lx-dx) {z=y;y=x;x=z;lz=ly;ly=lx;lx=lz;lz=dy;dy=dx;dx=lz;}
  lz=min(lx*dy,ly*dx);
  z=cgetg(lz,17);for(i=1;i<lz;i++) z[i]=zero;
  for(j=dx;j<lx;j++)
  {
    p1=(GEN)x[j];
    if(!gcmp0(p1))
    {
      if(gcmp1(p1))
	for(i=j*dy;i<lz;i+=j) z[i]=ladd((GEN)z[i],(GEN)y[i/j]);
      else
      {
	if(gcmp_1(p1))
	  for(i=j*dy;i<lz;i+=j) z[i]=lsub((GEN)z[i],(GEN)y[i/j]);
	else
	  for(i=j*dy;i<lz;i+=j) z[i]=ladd((GEN)z[i],gmul(p1,(GEN)y[i/j]));
      }
    }
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(z));
}

GEN
dirdiv(GEN x, GEN y)
{
  long lx,ly,lz,dx,dy,av,tetpil,i,j;
  GEN z,p1;
  
  if((typ(x)!=17)||(typ(y)!=17)) err(talker,"not a dirseries in dirmul");
  av=avma;dx=dirval(x);dy=dirval(y);lx=lg(x);ly=lg(y);
  if(dy!=1) err(talker,"not an invertible dirseries in dirdiv");
  lz=min(lx,ly*dx);p1=(GEN)y[1];
  if(!gcmp1(p1)) {y=gdiv(y,p1);x=gdiv(x,p1);}
  else x=gcopy(x);
  z=cgetg(lz,17);for(i=1;i<dx;i++) z[i]=zero;
  for(j=dx;j<lz;j++)
  {
    p1=(GEN)x[j];z[j]=(long)p1;
    if(!gcmp0(p1))
    {
      if(gcmp1(p1))
	for(i=j+j;i<lz;i+=j) x[i]=lsub((GEN)x[i],(GEN)y[i/j]);
      else
      {
	if(gcmp_1(p1))
	  for(i=j+j;i<lz;i+=j) x[i]=ladd((GEN)x[i],(GEN)y[i/j]);
	else
	  for(i=j+j;i<lz;i+=j) x[i]=lsub((GEN)x[i],gmul(p1,(GEN)y[i/j]));
      }
    }
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(z));
}

GEN
dirzetak(GEN nf, GEN b)
{
  long av=avma,tetpil,av1,lim,i,j,lx,bs,c,q,p2s;
  byteptr d=diffptr;
  GEN p,p1,p2,z,z2,pol;

  if((typ(nf)!=17)||(lg(nf)<10)) err(idealer1);
  lim=(avma+bot)>>1;pol=(GEN)nf[1];
  if(typ(b)!=1) err(talker,"not an integer type in dirzeta");
  if(signe(b)<=0) b=gun;
  if(gcmpgs(b,MAXHALFULONG)>=0) bs=MAXHALFULONG;
  else bs=itos(b)+1;
  p=cgeti(3);p[1]=evalsigne(1)+evallgef(3);p[2]=2;c= *d++;
  av1=avma;
  z=cgetg(bs,17);z[1]=un;for(i=2;i<bs;i++) z[i]=zero;
  do
  {
    if(signe(modis((GEN)nf[4],p[2])))
    {p1=(GEN)simplefactmod(pol,p)[1];lx=lg(p1);}
    else
    {
      p2=primedec(nf,p);lx=lg(p2);p1=cgetg(lx,18);
      for(i=1;i<lx;i++) p1[i]=((GEN)p2[i])[4];
    }
    for(j=1;j<lx;j++)
    {
      p2=gpuigs(p,itos((GEN)p1[j]));
      if(gcmpgs(p2,bs)<0)
      {
	z2=gcopy(z);q=p2s=itos(p2);
	while(q<bs)
	{
	  for(i=q;i<bs;i+=q) z2[i]=ladd((GEN)z2[i],(GEN)z[i/q]);
	  q*=p2s;
	}
	z=z2;
      }
    }
    if(avma<lim) {tetpil=avma;z=gerepile(av1,tetpil,gcopy(z));}
    c= *d++;p[2]+=c;
  }
  while(c&&(p[2]<bs));
  tetpil=avma;return gerepile(av,tetpil,gcopy(z));
}
