/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                       BASE D'ENTIERS                            */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

# include "genpari.h"

#define TRUE 1
#define FALSE 0

GEN pradical(GEN nf,GEN p),pol_min(GEN alpha,GEN nf,GEN p,GEN algebre,GEN algebre1);
GEN eval_pol(GEN nf,GEN pol,GEN alpha,GEN p,GEN algebre,GEN algebre1);
GEN lens(GEN nf,GEN p,GEN a),two_elt(GEN nf,GEN p,GEN ideal);
GEN element_mulid(GEN nf, GEN x, long i);
GEN element_muli(GEN nf, GEN x, GEN y);
GEN nfker(GEN nf, GEN R, GEN prhall);
GEN nfsuppl(GEN nf, GEN x, long n, GEN prhall);
GEN nfidealdet1(GEN nf, GEN a, GEN b);
long idealvalint(GEN nf, GEN x, GEN vp);

GEN maxord(GEN p,GEN f,long mf);
GEN nilord(GEN p,GEN fx,long mf,GEN gx);
GEN Decomp(GEN p,GEN f,long mf,GEN theta,GEN chi,GEN nu);
GEN respm(GEN f1,GEN f2,GEN pm);
GEN nbasis(GEN ibas,GEN pd);
GEN eltppm(GEN f,GEN pd,GEN theta,GEN k);
GEN testb(GEN p,GEN fa,long Da,GEN theta,long Dt);
GEN testc(GEN p, GEN fa, long c, GEN alph2, long Ma, GEN thet2, long Mt);
GEN testd(GEN p,GEN fa,long c,long Da,GEN alph2,long Ma,GEN theta);
GEN csrch(GEN p,GEN fa,GEN gamma);
GEN dedek(GEN f, long mf, GEN p,GEN g);
GEN dbasis(GEN p, GEN f, long mf, GEN alpha, GEN U);
long cbezout(long a,long b,long *u,long *v);
long clcm(long a,long b);

#define coef1(a,i,j)      (*((long*)(*(a+(j)+1))+(i)+1))
#define gcoef1(a,i,j)     (GEN)coef1(a,i,j)

GEN rquot(GEN x, GEN y),ordmax(GEN f, GEN p, GEN e, GEN *ptdelta),rtran(GEN v, GEN w, GEN q),mtran(GEN v, GEN w, GEN q, GEN m),matinv(GEN x, GEN d);

void rowred(GEN a, long rlim, GEN rmod);

/*******************************************************************
                             ROUND 2

  Entree:     x polynome unitaire a coefficients dans Z de deg n
	    definissant un corps de nombres K=Q(theta);
              code 0, 1 ou (long)p selon que l'on veut base, smallbase
            ou factoredbase.
	      y pointeur sur un GEN destine a recevoir
	    le discriminant du corps K.
  Sortie:    retourne 1) un vecteur (horizontal) a n composantes, 
            de polynomes a coeff dans Q (de deg 0,1...n-1)
	    constituant une base de l'anneau des entiers de K.
	        2) le discriminant de K (dans *y).
	    Rem: le denominateur commun des coef. est dans da.
*******************************************************************/

GEN
allbase(GEN x, long code, GEN *y)
{
  GEN p,a,at,bt,b,da,db,q;
  long av=avma,tetpil,n,h,j,je,k,r,s,t,pro,v;

  if(typ(x)!=10) err(allbaser1);
  n=lgef(x)-3;if(n<=0) err(allbaser1);
  v=varn(x);*y=discsr(x);
  if(DEBUGLEVEL) timer2();
  switch(code)
  {
    case 0: p=auxdecomp(absi(*y),1);h=lg((GEN)p[1])-1;break; /* base */
    case 1: p=auxdecomp(absi(*y),0);h=lg((GEN)p[1])-1;break; /* smallbase */
    default: p=(GEN)code;
      if((typ(p)!=19)||(lg(p)!=3)) err(factoreder1); /* factoredbase */
      h=lg((GEN)p[1])-1;
      q=gun;for(je=1;je<=h;je++) q=gmul(q,gpui(gcoeff(p,je,1),gcoeff(p,je,2),0));
      if(gcmp(absi(q),absi(*y))) err(factoreder2);
  }
  if(DEBUGLEVEL>0) {fprintferr("temps factpol: ");fprintferr("%ld\n",timer2());flusherr();}
  a=idmat(n);da=gun;
  for(je=1;je<=h;je++)
  {
    if(gcmpgs(gcoeff(p,je,2),1)>0)
    {
      b=ordmax(x,gcoeff(p,je,1),gcoeff(p,je,2),&db);
      a=gmul(db,a);b=gmul(da,b);
      da=mulii(db,da);db=da;
      at=gtrans(a);bt=gtrans(b);
      for(r=n-1;r>=0;r--)
	for(s=r;s>=0;s--)
	  while(signe(gcoef1(bt,s,r)))
	  {
	    q=rquot(gcoef1(at,s,s),gcoef1(bt,s,r));
	    at[s+1]=(long)rtran((GEN)at[s+1],(GEN)bt[r+1],q);
	    for(t=s-1;t>=0;t--)
	    {
	      q=rquot(gcoef1(at,t,s),gcoef1(at,t,t));
	      at[s+1]=(long)rtran((GEN)at[s+1],(GEN)at[t+1],q);
	    }
	    pro=at[s+1];at[s+1]=bt[r+1];bt[r+1]=pro;
	  }
      for(j=n-1;j>=0;j--)
      {
	for(k=0;k<j;k++)
	{
	  while(signe(gcoef1(at,j,k)))
	  {
	    q=rquot(gcoef1(at,j,j),gcoef1(at,j,k));
	    at[j+1]=(long)rtran((GEN)at[j+1],(GEN)at[k+1],q);
	    pro=at[j+1];at[j+1]=at[k+1];at[k+1]=pro;
	  }
	}
	if(signe(gcoef1(at,j,j))<0)
	  for(k=0;k<=j;k++) coef1(at,k,j)=lnegi(gcoef1(at,k,j));
	for(k=j+1;k<n;k++)
	{
	  q=rquot(gcoef1(at,j,k),gcoef1(at,j,j));
	  at[k+1]=(long)rtran((GEN)at[k+1],(GEN)at[j+1],q);
	}
      }
      for(j=1;j<n;j++)
	if(!cmpii(gcoef1(at,j,j),gcoef1(at,j-1,j-1)))
	{
	  coef1(at,0,j)=zero;
	  for(k=1;k<=j;k++)
	    coef1(at,k,j)=coef1(at,k-1,j-1);
	}
      a=gtrans(at);
    }
  }
  for(j=1;j<=n;j++)
  {
    *y=divii(mulii(gcoeff(a,j,j),*y),da);
    *y=divii(mulii(gcoeff(a,j,j),*y),da);
  }
  tetpil=avma;*y=gcopy(*y);at=cgetg(n+1,17);
  for(j=1;j<=n;j++)
  {
    q=cgetg(j+2,10);q[1]=evalsigne(1)+evallgef(2+j)+evalvarn(v);at[j]=(long)q;
    for(k=2;k<=j+1;k++) q[k]=ldiv(gcoeff(a,j,k-1),da);
  }
  pro=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;at+=pro;(*y)+=pro;
  return at;
}

GEN
base(GEN x, GEN *y)
{
  return allbase4(x,0,y,(GEN *)0);
}

GEN
base2(GEN x, GEN *y)
{
  return allbase(x,0,y);
}

GEN
smallbase(GEN x, GEN *y)
{
  return allbase4(x,1,y,(GEN *)0);
}

GEN
factoredbase(GEN x, GEN p, GEN *y)
{
  return allbase4(x,(long)p,y,(GEN *)0);
}

GEN
discf(GEN x)
{
  GEN y;
  long av,tetpil;

  av=avma;allbase4(x,0,&y,(GEN *)0);tetpil=avma;
  return gerepile(av,tetpil,gcopy(y));
}

GEN
discf2(GEN x)
{
  GEN y;
  long av,tetpil;

  av=avma;allbase(x,0,&y);tetpil=avma;
  return gerepile(av,tetpil,gcopy(y));
}

GEN
smalldiscf(GEN x)
{
  GEN y;
  long av,tetpil;

  av=avma;allbase4(x,1,&y,(GEN *)0);tetpil=avma;
  return gerepile(av,tetpil,gcopy(y));
}

GEN
factoreddiscf(GEN x, GEN p)
{
  GEN y;
  long av,tetpil;

  av=avma;allbase4(x,(long)p,&y,(GEN *)0);tetpil=avma;
  return gerepile(av,tetpil,gcopy(y));
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*   Quotient et Reste normalises   ( -1/2 < r = x-q*y <= 1/2 )    */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

GEN
rquot(GEN x, GEN y)
{
  GEN u,v,w,p;
  long av,av1;

  av=avma;
  u=absi(y);v=shifti(x,1);w=shifti(y,1);
  if ( cmpii(u,v)>0) p=subii(v,u);
  else p=addsi(-1,addii(u,v));
  av1=avma;
  return gerepile(av,av1,divii(p,w));
}
 
GEN
rrmdr(GEN x, GEN y)
{
  GEN p;
  long av,av1;

  av=avma;
  p=mulii(rquot(x,y),y);
  av1=avma;
  return gerepile(av,av1,subii(x,p));
}

GEN
rinv(GEN x, GEN y)
{
  GEN a,c,q,r,t;
  long av,av1;

  av=avma;
  a=gun;c=gzero;
  while( signe(y))
  {
    q=rquot(x,y);
    r=subii(a,mulii(q,c));a=c;c=r;
    t=subii(x,mulii(q,y));x=y;y=t;
  }
  av1=avma;
  if (signe(x)<0) a=negi(a);
  if (signe(c)) { av1=avma; a=rrmdr(a,c); }
  return gerepile(av,av1,a);
}

GEN
rgcd(GEN x, GEN y)
{
  GEN z;
  long av,av1;

  av=avma;
  while(signe(y))
  {
    z=rrmdr(x,y);x=y;y=z;
  }
  av1=avma;
  return gerepile(av,av1,absi(x));
}

GEN
rlcm(GEN x, GEN y)
{
  GEN d,z;
  long av,av1;

  av=avma;
  z=mulii(x,y);d=rgcd(x,y);
  av1=avma;
  return gerepile(av,av1,divii(z,d));
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*           Matrice compagnon du polynome unitaire x              */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

GEN
companion(GEN x)
{
  long    i,j,l;
  GEN     y;
  
  l=lgef(x)-2;y=cgetg(l,19);
  for(i=1;i<l;i++) y[i]=lgetg(l,18);
  for(i=0;i<l-2;i++)
    for(j=0;j<l-1;j++) coef1(y,i,j)=((i+1)==j) ? un : zero;
  for(j=0;j<l-1;j++) coef1(y,l-2,j)=lneg((GEN)x[j+2]);
  return y;
}



GEN
ordmax(GEN f, GEN p, GEN e, GEN *ptdelta)
{
  GEN m,hh,pp,dd,ppdd,index,q,r,s,b,c,t,jp,v,delta;
  GEN cf[100],w[100],a;
  long h,i,j,k,sp,epsilon,n=lgef(f)-3,av=avma,tetpil,av3,dec;

  a=cgetg(n*n+1,19);
  for(j=1;j<=n*n;j++)
  {
    a[j]=lgetg(n+1,18);
    for(k=1;k<=n;k++) coeff(a,k,j)=zero;
  }
  v=cgetg(n+1,18);
  cf[0]=idmat(n);
  cf[1]=companion(f);
  for(j=2;j<n;j++) cf[j]=gmul(cf[1],cf[j-1]);
  delta=gun; epsilon=itos(e);
  m=idmat(n);

  do
  {
    pp=mulii(p,p);
    dd=mulii(delta,delta);
    ppdd=mulii(dd,pp);
    b=matinv(m,delta);
    for(i=0;i<n;i++)
    {
      t=gscalsmat(0,n); /* t <--- matrice nulle d'ordre n */
      for(h=0;h<n;h++)
	for(j=0;j<n;j++)
	  for(k=0;k<n;k++)
	    coef1(t,j,k)=(long)rrmdr(addii(gcoef1(t,j,k),mulii(gcoef1(m,i,h),gcoef1(cf[h],j,k))),ppdd);
      c=gmul(t,b);
      w[i]=gmul(m,c);
      for(j=0;j<n;j++)
	for(k=0;k<n;k++)
	  coef1(w[i],j,k)=(long)rrmdr(divii(gcoef1(w[i],j,k),dd),pp);
    }
    if(cmpis(p,n)>0)
    {
      for(i=0;i<n;i++)
	for(j=0;j<n;j++)
	{
	  coeff(t,i+1,j+1)=zero;
	  for(k=0;k<n;k++)
	    for(h=0;h<n;h++)
	    {
	      r=modii(gcoef1(w[i],k,h),p);
	      s=modii(gcoef1(w[j],h,k),p);
	      coef1(t,i,j)=lmodii(addii(gcoef1(t,i,j),mulii(r,s)),p);
	    }
	}
    }
    else
    {
      for(j=0;j<n;j++)
      {
	for(i=0;i<n;i++)
	  coef1(b,i,j)=(i==0)? un : zero;
/* ici la boucle en k calcule la puissance p mod p de w[j] */
	sp=itos(p);
	for(k=0;k<sp;k++)
	{
	  for(i=0;i<n;i++)
	  {
	    v[i+1]=zero;
	    for(h=0;h<n;h++)
	      v[i+1]=lmodii(addii((GEN)v[i+1],mulii(gcoef1(b,h,j),gcoef1(w[j],h,i))),p);
	  }
	  for(i=0;i<n;i++) coef1(b,i,j)=v[i+1];
	}
      }
      q=p;t=b;
      while(cmpis(q,n)<0)
      {
	q=mulii(q,p);
	t=gmul(b,t);
      }
    }
    for(i=0;i<n;i++)
      for(j=0;j<n;j++)
      {
	coef1(a,j,i)=(i==j)? (long)p : zero;
	coef1(a,j,n+i)=lmodii(gcoef1(t,i,j),p);
      }
    rowred(a,2*n-1,pp);
    for(i=0;i<n;i++)
      for(j=0;j<n;j++)
	coef1(b,j,i)=coef1(a,j,i);
    jp=matinv(b,p);
    for(k=0;k<n;k++)
    {
      t=gmul(jp,w[k]);
      t=gmul(t,b);
      for(i=0;i<n;i++)
	for(j=0;j<n;j++)
	  coef1(t,i,j)=ldivii(gcoef1(t,i,j),p);
      h=0;
      for(i=0;i<n;i++)
	for(j=0;j<n;j++)
	{
	  coef1(a,k,h)=coef1(t,i,j);
	  h++;
	}
    }
    rowred(a,n*n-1,pp);
    index=gun;
    for(i=0;i<n;i++)
      index=mulii(index,gcoef1(a,i,i));
    if (cmpsi(1,index))
    {
      delta=mulii(index,delta);
      for(i=0;i<n;i++)
	for(j=0;j<n;j++)
	  coef1(c,i,j)=coef1(a,i,j);
      b=matinv(c,index);
      m=gmul(b,m);
      hh=delta;
      for(i=0;i<n;i++)
	for(j=0;j<n;j++)
	  hh=rgcd(gcoef1(m,i,j),hh);
      if(cmpis(hh,1)>1)
      {
	delta=divii(delta,hh);
	for(i=0;i<n;i++)
	  for(j=0;j<n;j++)
	    coef1(m,i,j)=ldivii(gcoef1(m,i,j),hh);
      }
      q=index;
      while(!signe(modii(q,p)))
      {
	q=divii(q,p);
	epsilon=epsilon-2;
      }
    }
  }
  while(!gcmp1(index) && (epsilon>=2));
  tetpil=avma;delta=gcopy(delta);m=gcopy(m);
  av3=avma;dec=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;
  *ptdelta=adecaler(delta,tetpil,av3)?delta+dec:delta;
  return m+dec;
}

void
rowred(GEN a, long rlim, GEN rmod)
{
  long j,k,n,pro;
  GEN q;

  n=lg((GEN)a[1])-1;
  for(j=0;j<n;j++)
  {
    for(k=j+1;k<=rlim;k++)
      while (signe(gcoef1(a,j,k)))
      {
	q=rquot(gcoef1(a,j,j),gcoef1(a,j,k));
	a[j+1]=(long)mtran((GEN)a[j+1],(GEN)a[k+1],q,rmod);
	pro=a[j+1];a[j+1]=a[k+1];a[k+1]=pro;
      }
    if (signe(gcoef1(a,j,j))<0)
      for(k=j;k<n;k++) coef1(a,k,j)=lnegi(gcoef1(a,k,j));
    for(k=0;k<j;k++)
    {
      q=rquot(gcoef1(a,j,k),gcoef1(a,j,j));
      a[k+1]=(long)mtran((GEN)a[k+1],(GEN)a[j+1],q,rmod);
    }
  }
}

GEN
rtran(GEN v, GEN w, GEN q)
{
  long av,tetpil;
  GEN p1;

  if (signe(q))
  {
    av=avma;p1=gmul(q,w);tetpil=avma;
    return gerepile(av,tetpil,gsub(v,p1));
  }
  else return v;
}

GEN
mtran(GEN v, GEN w, GEN q, GEN m)
{
  long k;
  
  if (signe(q))
  {
    for(k=0;k<lg(v)-1;k++)
    {
      v[k+1]=(long)rrmdr(subii((GEN)v[k+1],modii(mulii(q,(GEN)w[k+1]),m)),m);
    }
  }
  return v;
}


GEN matinv(GEN x, GEN d)
             
/*=======================================================================
    Calcule d/x  ou  d est entier et x matrice triangulaire inferieure
  entiere dont les coeff diagonaux divisent
  d ( resultat entier).
========================================================================*/
{
  long n,i,j,k,av,av1;
  GEN y,h;

  av=avma;
  y=idmat(n=lg(x)-1);
  for(i=1;i<=n;i++)
    coeff(y,i,i)=ldivii(d,gcoeff(x,i,i));
  for(i=2;i<=n;i++)
    for(j=i-1;j;j--)
    {
      for(h=gzero,k=j+1;k<=i;k++)
	h=gadd(h,mulii(gcoeff(y,i,k),gcoeff(x,k,j)));
      coeff(y,i,j)=ldivii(negi(h),gcoeff(x,j,j));
    }
  av1=avma;
  return gerepile(av,av1,gcopy(y));
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                    BASE D'ENTIERS (ROUND 4)                     */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

int
fnz(GEN x,long j)
{
  long i;
  i=1;while(!signe((GEN)x[i])) i++;
  return i==j;
}

GEN
allbase4(GEN f,long code, GEN *y, GEN *ptw)
{
/* retourne la base, dans y le discf et dans ptw la factorisation (peut
 etre partielle) de discf */
   
  GEN w,w1,w2,a,da,b,db,bas,q,bdiag,ab,centre;
  long v,n,mf,h,templevel,lfa;
  long l,i,j,k,av=avma,tetpil,pro,first;
  
  if(typ(f)!=10) err(allbaser1);
  n=lgef(f)-3;if(n<=0) err(allbaser1);
  v=varn(f);
  *y=discsr(f);
  if(DEBUGLEVEL) {timer2();templevel=DEBUGLEVEL;DEBUGLEVEL=5;}
  switch(code)
  {
    case 0: w=auxdecomp(absi(*y),1);h=lg((GEN)w[1])-1;break; /* base */
    case 1: w=auxdecomp(absi(*y),0);h=lg((GEN)w[1])-1;break; /* smallbase */
    default: w=(GEN)code;
      if((typ(w)!=19)||(lg(w)!=3)) err(factoreder1); /* factoredbase */
      h=lg((GEN)w[1])-1;
      q=gun;for(i=1;i<=h;i++) q=gmul(q,gpui((GEN)coeff(w,i,1),(GEN)coeff(w,i,2),0));
      if(gcmp(absi(q),absi((GEN)*y))) err(factoreder2);
  }
  if(DEBUGLEVEL) 
  {
    DEBUGLEVEL=templevel;
    fprintferr("temps factorisation disc: ");fprintferr("%ld\n",timer2());flusherr();
  }
  a=idmat(n);da=gun;first=TRUE;
  for(i=1;i<=h;i++)
  {
    mf=itos((GEN)coeff(w,i,2));
    if(mf>1)
    { 
      if(DEBUGLEVEL)
      {
	fprintferr("On traite le cas p^k = ");
	bruterr((GEN)coeff(w,i,1),'g',-1);
	fprintferr("^%ld\n",mf);
      }
      b=(GEN)maxord((GEN)coeff(w,i,1),f,mf);
      bdiag=cgetg(n+1,17);for(j=1;j<=n;j++) bdiag[j]=coeff(b,j,j);
      db=denom(bdiag); 
      if (!(gcmp1(db)))    /* la matrice est identite   */
      {
	da=gmul(da,db);
	if (first!=TRUE)
	{
	  b=gmul(da,b);a=gmul(db,a);
	  for(j=1;(j<=n)&&(fnz((GEN)a[j],j)&&fnz((GEN)b[j],j));j++);
	  k=j-1;ab=cgetg(2*n-k+1,19);
	  for(j=1;j<=k;j++)
	  {
	    ab[j]=a[j];
	    coeff(ab,j,j)=(long)mppgcd(gcoeff(a,j,j),gcoeff(b,j,j));
	  }
	  for(;j<=n;j++) ab[j]=a[j];
	  for(;j<=2*n-k;j++) ab[j]=b[j+k-n];
	  a=hnf(ab);
	}
	else {a=gmul(b,db);first=FALSE;}
      }
      if(DEBUGLEVEL>=3)
      {
	fprintferr("Le resultat pour ce nombre p est : \n ");
	outerr(b);fprintferr("\n");
      }
    }
  } 
  for(j=1;j<=n;j++)
  {
    *y=divii(mulii((GEN)coeff(a,j,j),*y),da);
    *y=divii(mulii((GEN)coeff(a,j,j),*y),da);
  }
  if(ptw)
  {
    w1=(GEN)w[1];w2=(GEN)w[2];lfa=0;
    for(j=1;j<=h;j++)
    {
      k=ggval(*y,(GEN)w1[j]);
      w2[j]=lstoi(k);if(k) lfa++;
    }
  }
  tetpil=avma;
  *y=gcopy(*y);bas=cgetg(n+1,17);
  for(j=n-1;j>0;j--)
    if(cmpis((GEN)coeff(a,j,j),2)==1)
    {
      centre=shifti((GEN)coeff(a,j,j),-1);
      for(k=j+1;k<=n;k++)
	if(cmpii((GEN)coeff(a,j,k),centre)==1)
	  for(l=1;l<=j;l++)
	    coeff(a,l,k)=lsubii((GEN)coeff(a,l,k),(GEN)coeff(a,l,j));
    }

  for(k=1;k<=n;k++)
  {
    q=cgetg(k+2,10);q[1]=evalsigne(1)+evallgef(2+k)+evalvarn(v);bas[k]=(long)q;
    for(j=2;j<=k+1;j++) q[j]=ldiv((GEN)coeff(a,j-1,k),da);
  }
  if(ptw)
  {
    w=cgetg(3,19);w[1]=lgetg(lfa+1,18);w[2]=lgetg(lfa+1,18);
    for(l=0,j=1;j<=h;j++)
      if(signe((GEN)w2[j]))
      {l++;coeff(w,l,1)=lcopy((GEN)w1[j]);coeff(w,l,2)=lcopy((GEN)w2[j]);}
  }
  pro=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;
  bas+=pro;(*y)+=pro;if(ptw) *ptw=w+pro;
  return bas;
}

    /*     p-maximal order of Af        */
    /*     p^m does not divide Df       */

GEN
maxord(GEN p,GEN f,long mf)
{
  GEN w,g,h,res,fmp;
  long j,r,v=varn(f),n=lgef(f)-3,av=avma,tetpil,flw;

  if((flw=(cmpsi(n,p)<0)))
  {
    fmp=gmul(gmodulcp(gun,p),f);
    g=gdeuc(fmp,polgcd(fmp,deriv(fmp,v)));
  }
  else 
  {
    w=factmod(f,p);r=lg((GEN)w[1])-1;g=gun;
    for(j=1;j<=r;j++) g=gmul((GEN)coeff(w,j,1),g);
  }
  res=dedek(f, mf,p,g);
  if(itos((GEN)res[1]))
  {
    tetpil=avma;
    return gerepile(av,tetpil,dbasis(p,f,mf,(GEN)polx[v],(GEN)res[2]));
  }
  else
  {
    if(flw) {w=factmod(f,p);r=lg((GEN)w[1])-1;}
    h=bestnu(w);
    if (r==1)
    { 
      res=nilord(p,f,mf,h);
      tetpil=avma;return gerepile(av,tetpil,gcopy(res));
    }
    else
    { 
      tetpil=avma;return gerepile(av,tetpil,Decomp(p,f,mf,polx[v],f,h));
    }
  }
}

GEN
dedek(GEN f, long mf, GEN p,GEN g)
/* Return res[1] = 1 : if Z[alpha] is maximal or 2*dU >= m-1 else return 0 */
/* Return res[2] = U if res[2] == 1 else res[2] = f */
{
  
  long av=avma,tetpil,dk;
  GEN k,h,unmodp,res;
  
  res=cgetg(3,17);
  if(DEBUGLEVEL>=3)
  {
    fprintferr(" On est dans Dedekind ");
    if(DEBUGLEVEL>=4)
    {
      fprintferr(" avec les parametres \n" );
      fprintferr(" p=");bruterr(p,'g',-1);
      fprintferr(",  f=");bruterr(f,'g',-1);
    }
    fprintferr("\n");
  } 
  unmodp=gmodulcp(gun,p);
  g=gmul(g,unmodp);
 
  h=gdivent(gmul(f,unmodp),g);
  k=gdiv(gsub(lift(f),gmul(lift(g),lift(h))),p);
  k=ggcd(gmul(k,unmodp),ggcd(g,h));
  dk=lgef(k)-3;
  if(DEBUGLEVEL>=4)
    fprintferr(" Le pgcd est de degre %ld \n",dk );
  res[1]= ((dk==0)||(2*dk >= mf-1))?un:zero;
  if (dk!=0)
    res[2]= (long)lift(gdiv(gmul(f,unmodp),k));
  else
    res[2]= lcopy(f);
  tetpil=avma;
  return gerepile(av,tetpil,gcopy(res));
}

GEN
dbasis(GEN p, GEN f, long mf, GEN alpha, GEN U)
{
  long av=avma,tetpil,n=lgef(f)-3,m;
  long dU=lgef(U)-3,c,i,dh;
  GEN unmodpdd,b,p1,ha,pd;

  if(DEBUGLEVEL>=3)
  {
    fprintferr(" On est dans Dedekind Basis ");
    if(DEBUGLEVEL>=4)
    {
      fprintferr(" avec les parametres \n" );
      fprintferr(" p=");bruterr(p,'g',-1);
      fprintferr(",  f=");bruterr(f,'g',-1);
      fprintferr(",  alpha=");bruterr(alpha,'g',-1);
    }
    fprintferr("\n");
  }     
  m = n - dU;
  pd = gpuigs(p, mf/2 );
  unmodpdd = gmodulcp(gun,gmul(pd,p));
  
  b=cgetg(n+1,19);                   /* Z[a] + U/p Z[a] is maximal */
  ha = pd;
  
  p1=cgetg(n+1,18);b[1]=(long)p1;
  p1[1]=(long)pd;for(i=2;i<=n;i++) p1[i]=zero;
  for(c=2;c<=n;c++)
  {
    p1=cgetg(n+1,18);b[c]=(long)p1;
    if( c == n-m+1)
      ha = lift(gmul(gdiv(gmul(pd,eleval(f, U, alpha)),p),unmodpdd));
    else
      ha = lift(gmul(gmod(gmul(ha,alpha),f),unmodpdd));
    dh = lgef(ha) - 3 ;
    for(i=1;i<=dh+1;i++) p1[i]=ha[1+i];
    for(i=dh+2;i<=n;i++) p1[i]=zero;
  }
  if(DEBUGLEVEL>=4)
  {
    fprintferr(" On construit un nouvel ordre  \n" );
    if(DEBUGLEVEL>=5) outerr(b);
    fprintferr(" On fait sa HNF \n");
  }
  b=gdiv(hnfmodid(b,pd),pd);
  if(DEBUGLEVEL>=4)
  {
    fprintferr(" Sa HNF est finie \n");
    if(DEBUGLEVEL>=5) outerr(b);
  }

  tetpil=avma;
  return gerepile(av,tetpil,gcopy(b));
}

GEN
Decomp(GEN p,GEN f,long mf,GEN theta,GEN chi,GEN nu)
{
  long n1,n2,j,i,v2,v1,v=varn(f),av=avma,tetpil;
  GEN unmodpdr,unmodp,unmodpdrp,unmodpkpdr,unmodpmr;
  GEN pk,ph,pmr,pdr;
  GEN b1,b2,b3,a2,a1,e,f1,f2;
  GEN ib1,ib2,ibas,h;

  if(DEBUGLEVEL>=3)
  {
    fprintferr(" On entre dans Decomp ");
    if(DEBUGLEVEL>=4)
    {
      fprintferr(" avec les parametres suivants \n ");
      fprintferr(" p=");bruterr(p,'g',-1);
      fprintferr(",  f=");bruterr(f,'g',-1);
      fprintferr(",  exposant=%ld ",mf);
    }
    fprintferr("\n");
  }
  
  unmodp=gmodulcp(gun,p);

  pdr=(GEN)respm(f,deriv(f,v),gpuigs(p,mf));
  pmr=mulii(pdr,mulii(pdr,p));

  unmodpmr=gmodulcp(gun,pmr);
  unmodpdr=gmodulcp(gun,pdr);
  unmodpdrp=gmodulcp(gun,mulii(pdr,p));

  b1=gmul(chi,unmodp);  a2=gzero;
  b2=unmodp;            a1=gun;
  b3=gmul(nu,unmodp);
  
  while ( lgef(b3)>3 )
  {
    b1=gdivent(b1,b3);
    b2=gmul(b2,b3);
    b3=lift(gbezout(b2,b1,&a1,&a2));
  }
                           
  e=(GEN)eleval(f,lift(gmul(a1,b2)),theta);
  e=gdiv(lift(gmul(gmul(pdr,e),unmodpdrp)),pdr); 

  pk=p;
  ph=mulii(pdr,pmr); 

 /*    E(t)- e(t) belongs to p^k Op, which is contained in p^(k-df)*Zp[xi]  */

  while (cmpii(pk,ph)==-1)
  {
    e=gmod(gmul(e,gmul(e,gsubsg(3,gmulsg(2,e)))),f);
    pk=gmul(pk,pk);
    unmodpkpdr=gmodulcp(gun,mulii(pk,pdr ));
    e=gdiv(lift(gmul(gmul(pdr,e),unmodpkpdr)),pdr);
  }  
  
  f1=(GEN)gcdpm(f,gmul(pdr,gsubsg(1,e)),mulii(pmr,pdr));f1=lift(gmul(gmod(f1,f),unmodpmr));
  f2=gdivent(f,f1);f2=lift(gmul(gmod(f2,f),unmodpmr));

  n1=lgef(f1)-3;v1=ggval(discsr(f1),p); b1=(GEN)maxord(p,f1,v1);

  ib1=cgetg(n1+1,17);
  for(i=1;i<=n1;i++)
  {
    h=gzero;
    for(j=1;j<=i;j++)
      h=gadd(h,gmul((GEN)coeff(b1,j,i),gpuigs(polx[v],j-1)));
    ib1[i]=(long)h;
  }  
  
  n2=lgef(f2)-3; v2=ggval(discsr(f2),p); b2=(GEN)maxord(p,f2,v2);

  ib2=cgetg(n2+1,17);
  for(i=1;i<=n2;i++)
  {
    h=gzero;
    for(j=1;j<=i;j++)
      h=gadd(h,gmul((GEN)coeff(b2,j,i),gpuigs(polx[v],j-1)));
    ib2[i]=(long)h;
  }
  
  ibas=cgetg(n1+n2+1,17);


  for(j=1;j<=n1;j++)
    ibas[j]=(long)lift(gmul(gmod(gmul(gmul(pdr,(GEN)ib1[j]),e),f),unmodpdr));
  for(j=n1+1;j<=n1+n2;j++)
    ibas[j]=(long)lift(gmul(gmod(gmul(gsubsg(1,e),gmul(pdr,(GEN)ib2[j-n1])),f),unmodpdr));
  tetpil=avma;
  return gerepile(av,tetpil,nbasis(ibas,pdr));
 
}

GEN
nilord(GEN p,GEN fx,long mf,GEN gx)
{
  long La,Ma,first=TRUE,v=varn(fx),av=avma,tetpil; 
  GEN alpha,chi,nu,eta,w,phi;
  GEN res,pm,Dchi,unmodp,unmodpm;
  

  if(DEBUGLEVEL>=3)
  {
    fprintferr(" On entre dans Nilord ");
    if(DEBUGLEVEL>=4)
    {
      fprintferr(" avec les parametres suivants \n ");
      fprintferr(" p=");bruterr(p,'g',-1);
      fprintferr(",  fx=");bruterr(fx,'g',-1);
      fprintferr(",  exposant=%ld,  gx= ",mf);bruterr(gx,'g',-1);
    }
    fprintferr("\n");
  } 
  
  pm=gpuigs(p,mf+1);

  alpha=polx[v];     chi=fx;     nu=gx;       Dchi=gpuigs(p,mf);

  unmodpm=gmodulcp(gun,pm);
  unmodp=gmodulcp(gun,p);
 
  res=cgetg(4,17);

  while (TRUE)
  {
    if (gcmp0(Dchi))
      alpha=gadd(alpha,gmul(p,polx[v]));
    else
    {
      if (first!=TRUE) res=dedek(chi, mf, p,nu);
      else {res[1]=zero;res[2]=(long)alpha;first=FALSE;}
      if (itos((GEN)res[1])==1) 
      { 
	tetpil=avma;
	return gerepile(av,tetpil,dbasis(p,fx,mf,alpha,(GEN)res[2]));
      }
      else
      { 
	if (gcmp(vstar(p,chi),gzero)==1)
	{
	  alpha=gadd(alpha,gun);
	  chi=gsubst(chi,v,gsub(polx[v],gun));
	  nu=lift(gmul(gsubst(nu,v,gsub(polx[v],gun)),unmodp));
	}
	w=(GEN)setup(p,chi,polx[v],nu);
	eta=(GEN)w[2];
	La=itos((GEN)w[3]);
	Ma=itos((GEN)w[4]);
	if (La>1)
	  alpha=gadd(alpha,eleval(fx,eta,alpha));
	else
	{
	  w=(GEN)bsrch(p,chi,ggval(Dchi,p),eta,Ma);
	  phi=(GEN)eleval(fx,(GEN)w[2],alpha);
	  if (gcmp1((GEN)w[1]))
	  {
	    tetpil=avma;
	    return gerepile(av,tetpil,Decomp(p,fx,mf,phi,(GEN)w[3],(GEN)w[4]));
	  }
	  else alpha=phi;
	}
      }
    }
    w=(GEN)factcp(p,fx,alpha);
    chi=(GEN)w[1];nu=(GEN)w[2];
    if(cmpis((GEN)w[4],1)==1)
    {
      tetpil=avma;
      return gerepile(av,tetpil,Decomp(p,fx,mf,alpha,chi,nu));
    }
    Dchi=lift(gmul(discsr(lift(gmul(chi,unmodpm))),unmodpm));
    if (gcmp0(Dchi))
      Dchi=discsr(chi);
     
  }
}


/***********************************************************************/
/****             returns                                           ****/
/****       [1,theta,chi,nu]  if theta non-primary                  ****/
/****       [2,phi, * , * ]   if D_phi > D_alpha or M_phi > M_alpha ****/
/***********************************************************************/

GEN
bsrch(GEN p,GEN fa,long ka,GEN eta,long Ma)
{
  long n=lgef(fa)-3,Da=lgef(eta)-3;
  long c,r,field,j,MaVb,deg,av=avma,tetpil;
  GEN pc,pcc,unmodpcc,Vb;
  GEN beta,b,gamma,delta,pik,w;
  
  pc=respm(fa,deriv(fa,varn(fa)),gpuigs(p,ka));
  c=ggval(pc,p);
  pcc=gmul(pc,pc);
  unmodpcc=gmodulcp(gun,pcc);
  
  r=1+(long)ceil(c/(double)(Da)+gtodouble(gdivsg(c*n-2,mulsi(Da,subis(p,1)))));
  
  b=cgetg(5,17);
  
  beta=gdiv(lift(gpuigs(gmodulcp(eta,fa),Ma)),p);
  
  while(TRUE)
  { 
    beta=gdiv(lift(gmul(gmul(pc,beta),unmodpcc)),pc);
    w=testd(p,fa,c,Da,eta,Ma,beta);
    if(cmpis((GEN)w[1],3)==-1) 
    { tetpil=avma;
    return gerepile(av,tetpil,gcopy((GEN)w));
    } 
      
    Vb=vstar(p,(GEN)w[3]);
      
    MaVb=itos(gmulsg(Ma,Vb));
      
    pik=lift(gpuigs(gmodulcp(eta,fa),MaVb));
      
    gamma=gmod(gmul(beta,(GEN)(vecbezout(pik,fa))[1]),fa);
    gamma=gdiv(lift(gmul(gmul(pc,gamma),unmodpcc)),pc);
    w=(GEN)testd(p,fa,c,Da,eta,Ma,gamma);
    if (cmpis((GEN)w[1],3)==-1) 
    {
      tetpil=avma;
      return gerepile(av,tetpil,gcopy((GEN)w));
    } 

    delta=eltppm(fa,pc,gamma,gpuigs(p,r*Da));
    delta=gdiv(lift(gmul(gmul(pc,delta),unmodpcc)),pc);
    w=(GEN)testd(p,fa,c,Da,eta,Ma,delta);
    if (cmpis((GEN)w[1],3)==-1)
    {
      tetpil=avma;
      return gerepile(av,tetpil,gcopy((GEN)w));
    } 
      
    field=TRUE;
    deg=lgef(delta)-3;
    for(j=0;j<=deg;j++)
      if (!(gcmp0((GEN)delta[j+2])))
	if (ggval((GEN)delta[j+2],p) < 0)  field=FALSE;
    if (field) 
      beta=gsub(beta,gmod(gmul(pik,delta),fa));
    else
    { 
      tetpil=avma;
      return gerepile(av,tetpil,csrch(p,fa,gamma));
    } 
  }
}

/***********************************************************************/
/****    returns                                                    ****/
/****    [1,phi,chi,nu]      if theta non-primary                   ****/
/****    [2,phi,chi,nu]      if D_phi > D_aplha or M_phi > M_alpha  ****/
/****    [3,phi,chi,nu]      otherwise                              ****/
/***********************************************************************/

GEN
testd(GEN p,GEN fa,long c,long Da,GEN alph2,long Ma,GEN theta)
{
  long Mt,Dt,av=avma,tetpil;
  GEN chit,nut,thet2,b,w;
  
  b=cgetg(5,17);
  
  
  w=factcp(p,fa,theta);
  chit=(GEN)w[1];
  nut=(GEN)w[2];
  Dt=itos((GEN)w[3]);

  if (cmpis((GEN)w[4],1)==1)
  {
    b[1]=un;
    b[2]=(long)theta;
    b[3]=(long)chit;
    b[4]=(long)nut;
    tetpil=avma;
    return gerepile(av,tetpil,gcopy((GEN)b));
  } 

  if (Da< clcm(Da,Dt)) 
  { 
    tetpil=avma;
    return gerepile(av,tetpil,testb(p,fa,Da,theta,Dt));
  }
  
  w=setup(p,fa,theta,nut);
  thet2=(GEN)w[2];
  Mt=itos((GEN)w[4]);
  
  if (Ma < clcm(Ma,Mt))
  {
    tetpil=avma;
    return gerepile(av,tetpil,testc(p,fa,c,alph2,Ma,thet2,Mt));
  }
  else
  {   
    b[1]=(long)stoi(3);
    b[2]=(long)theta;
    b[3]=(long)chit;
    b[4]=(long)nut; 
    tetpil=avma;
    return gerepile(av,tetpil,gcopy((GEN)b));
  }
}


/***********************************************************************/
/*****    Returns [1,phi,chi,nu] if phi non-primary                *****/
/*****            [2,phi,chi,nu] if D_phi = lcm (D_alpha, D_theta) *****/
/***********************************************************************/

GEN
testc(GEN p, GEN fa, long c, GEN alph2, long Ma, GEN thet2, long Mt)

{
  GEN b,pc,ppc,c1,c2,c3,psi,unmodppc,phi,w;
  long g,r,s,t,v=varn(fa),av=avma,tetpil;

  b=cgetg(5,17);
  pc=gpuigs(p,c);
  ppc=mulii(pc,p);
  unmodppc=gmodulcp(gun,ppc);

  g=cbezout(Ma,Mt,&r,&s);
  t=0;
  while (r<0)
  {
    r=r+Mt;
    t++;
  }
  while (s<0)
  {
    s=s+Ma;
    t++;
  }
  c1=lift(gpuigs(gmodulcp(alph2,fa),s));
  c2=lift(gpuigs(gmodulcp(thet2,fa),r));
  c3=gdiv(gmod(gmul(c1,c2),fa),gpuigs(p,t));
  psi=gdiv(lift(gmul(gmul(pc,c3),unmodppc)),pc);
  phi=gadd(polx[v],psi);

  w=factcp(p,fa,phi);
  if(cmpis((GEN)w[4],1)==1)
  {
    b[1]=un;
    b[2]=(long)phi;
    b[3]=w[1];
    b[4]=w[2];
    tetpil=avma;
    return gerepile(av,tetpil,gcopy((GEN)b));
  }
  else
  {   
    b[1]=deux;
    b[2]=(long)phi;
    b[3]=w[1];
    b[4]=w[2]; 
    tetpil=avma;
    return gerepile(av,tetpil,gcopy((GEN)b));
  }
}


/************************************************************************/
/*****   Returns [1,phi,chi,nu] if phi non-primary                  *****/
/*****           [2,phi,chi,nu] if D_phi = lcm (D_alpha, D_theta)   *****/
/************************************************************************/


GEN
testb(GEN p,GEN fa,long Da,GEN theta,long Dt)
{
  long Dat,t,j,vf=varn(fa),av=avma,tetpil;
  GEN b,w,r,v;
  GEN phi,h;
  
  
  Dat=clcm(Da,Dt);
  b=cgetg(5,17);
  t=0;
  
  while (TRUE)
  {
    t++;
    v=stoi(t);
    h=gzero;
    j=0;
    while (!(gcmp0(v)))
    {
      r=gmod(v,p);
      v=gdivent(v,p);
      h=gadd(h,gmul(r,gpuigs(polx[vf],j)));
      j++;
    }
    phi=gadd(theta,gmod(h,fa));
    w=factcp(p,fa,phi);
    if (cmpis((GEN)w[4],1)==1)
    {
      b[1]=un;
      b[2]=(long)phi;
      b[3]=w[1];
      b[4]=w[2];
      tetpil=avma;
      return gerepile(av,tetpil,gcopy((GEN)b));
    }
    if (cmpis((GEN)w[3],Dat)==0)
    {
      b[1]=deux;
      b[2]=(long)phi;
      b[3]=w[1];
      b[4]=w[2];
      tetpil=avma;
      return gerepile(av,tetpil,gcopy((GEN)b));
    }
  }
}

/***********************************************************************/
/*****     Factorize characteristic polynomial of beta mod p       *****/
/***********************************************************************/

GEN
factcp(GEN p,GEN f,GEN beta)
{
  GEN chi,nu,b;
  long v,av=avma,tetpil;
  
  v=varn(f);
  chi=lift(caradj0(gmodulcp(beta,f),v));
  nu=lift(factmod(chi,p));
  
  b=cgetg(5,17);
  
  b[1]=(long)chi;
  b[2]=coeff(nu,1,1);
  b[3]=lstoi(lgef((GEN)b[2])-3);
  b[4]=lstoi(lg((GEN)nu[1])-1);
  tetpil=avma;
  return gerepile(av,tetpil,gcopy((GEN)b));
}

/***********************************************************************/
/*****************          minimum extension valuation     ************/
/***********************************************************************/

GEN
vstar(GEN p,GEN h)
{
  long m,first,j,av=avma,tetpil;
  GEN w,v;
  
  m=lgef(h)-3;
  first=TRUE;
  v=gzero;
  for(j=1;j<=m;j++)
    if (!(gcmp0((GEN)h[m-j+2])))
    {
      w=gdiv(stoi(ggval((GEN)h[m-j+2],p)),stoi(j));
      if (first) 
	v=w;
      else 
	if(gcmp(w,v)==-1) v=w;
      first=FALSE;
    }
  tetpil=avma;
  return gerepile(av,tetpil,gcopy(v));
}

/************************************************************************/
/*** Returns [theta_1,theta_2,L_theta,M_theta] with theta non-primary ***/
/***            [1]      [2]     [3]     [4]                         ****/
/************************************************************************/

GEN
setup(GEN p,GEN f,GEN theta,GEN nut)
{
  GEN b,t1,t2,v1;
  long Lt,Mt,r,s,c,v,av=avma,tetpil;
  
  v=varn(f);
  b=cgetg(5,17);
  
  t1=eleval(f,nut,theta);
  v1=vstar(p,lift(caradj0(gmodulcp(t1,f),v)));
  
  if (typ(v1)==1) 
  {     
    Lt=itos(v1);
    Mt=1;
  }
  else
  {
    Lt=itos((GEN)v1[1]);
    Mt=itos((GEN)v1[2]);
  }
  
  c=cbezout(Lt,-Mt,&r,&s);
  
  while(r<=0)
  {
    r=r+Mt;
    s=s+Lt;
  }
  t2=gdiv(lift(gpuigs(gmodulcp(t1,f),r)),gpuigs(p,s));
  
  b[1]=(long)t1;
  b[2]=(long)t2;
  b[3]=lstoi(Lt);
  b[4]=lstoi(Mt);
  
  tetpil=avma;
  return gerepile(av,tetpil,gcopy((GEN)b));
}

/***********************************************************************/
/**************          evaluate g(a)                ******************/
/***********************************************************************/

GEN
eleval(GEN f,GEN h,GEN a)
{
  long n,k,v=varn(f),av=avma,tetpil;
  GEN g,y;
  
  g=gmul(h,polun[v]);
  n=lgef(g)-3;
  y=gzero;
  for(k=n;k>=0;k--)
    y=gmod(gadd(gmul(y,a),(GEN)g[k+2]),f);
  tetpil=avma;
  return gerepile(av,tetpil,gcopy(y));
}

/************************************************************************/
/************ Returns [theta,chi,nu ] with theta non-primary ************/
/************************************************************************/

GEN
csrch(GEN p,GEN fa,GEN gamma)
{
  GEN b,h,theta,w,v,r;
  long t,j,vf=varn(fa),av=avma,tetpil;
  
  b=cgetg(5,17);
  
  t=0;
  while (TRUE)
  {
    t++;
    v=stoi(t);
    h=gzero;
    j=0;
    while (!(gcmp0(v)))
    {
      r=gmod(v,p);
      v=gdivent(v,p);
      h=gadd(h,gmul(r,gpuigs(polx[vf],j)));
      j++;
    }
    theta=gadd(gamma,gmod(h,fa));
    w=factcp(p,fa,theta);
    if (cmpis((GEN)w[4],1)==1)
    {
      b[1]=un;
      b[2]=(long)theta;
      b[3]=w[1];
      b[4]=w[2];
      tetpil=avma;
      return gerepile(av,tetpil,gcopy((GEN)b));
    }
  }
}

/***********************************************************************/
/*****************          Modular power of an elment        **********/
/***********************************************************************/

GEN
eltppm(GEN f,GEN pd,GEN theta,GEN k)
{
  GEN pdd,phi,psi,unmodpdd,q;
  long r,av=avma,tetpil;
  
  pdd=gmul(pd,pd);
  unmodpdd=gmodulcp(gun,pdd);
  phi=pd;
  psi=gmul(pd,theta);
  q=k;
  
  while (cmpis(q,0)!=0)
  { 
    r=q[lgef(q)-1]&1;
    if (r !=0)
    { 
      phi=gmod(gdiv(gmul(phi,psi),pd),f);
      phi=lift(gmul(phi,unmodpdd));
    }
    q=gshift(q,-1);
    if (cmpis(q,0) != 0)
    {
      psi=gmod(gdiv(gmul(psi,psi),pd),f);
      psi=lift(gmul(psi,unmodpdd));
    }
  }
  tetpil=avma;
  return gerepile(av,tetpil,gdiv(phi,pd));
}

/**********************************************************************/
/******       polynomial gcd mod p^m (assumes f1 monic)        ********/
/**********************************************************************/

GEN
gcdpm(GEN f1,GEN f2,GEN pm)
{
  long n,c,deg,k,j,v=varn(f1),av=avma,tetpil;
  GEN a,h,unmodpm,b;
  
  
  n=lgef(f1)-3;
  unmodpm=gmodulcp(gun,pm);
  a=cgetg(n+1,19);
  
  h=lift(gmul(gmod(f2,f1),unmodpm));
  for(k=1;k<=n;k++)
    a[k]=lgetg(n+1,18);
  for(j=1;j<=n;j++)
  { deg=lgef(h)-3;
  for(k=1;k<=deg+1;k++)
    coeff(a,k,j)=h[k+1];
  for(k=deg+2;k<=n;k++)
    coeff(a,k,j)=zero;
      
  if (j<n) h=lift(gmul(gmod(gmul(polx[v],h),f1),unmodpm));
  }

  a=hnfmodid(a,pm);c=0;
  
  for(j=n;j>=1;j--)
    if (!(gcmp0(lift(gmul((GEN)coeff(a,j,j),unmodpm)))))
      c=j;
  b=gmul(gzero,polun[v]);
  for(k=1;k<=c;k++)
    b=gadd(b,gmul(gdiv((GEN)coeff(a,k,c),(GEN)coeff(a,c,c)),gpuigs(polx[v],k-1)));
  tetpil=avma;
  return gerepile(av,tetpil,gcopy((GEN)b));
}

/***********************************************************************/
/*****      reduced resultant mod p^m (assumes f1 monic)        ********/
/***********************************************************************/

GEN
respm(GEN f1,GEN f2,GEN pm)
{
  long av=avma,tetpil;
  GEN a1,a2,pc,g,unmodpm;
 
  unmodpm=gmodulcp(gun,pm);
  g=gbezout(f1,f2,&a1,&a2);

  a1=lift(gmul(gmul(pm,a1),unmodpm));
  a2=lift(gmul(gmul(pm,a2),unmodpm));

  pc=ggcd(pm,content(a1));
  pc=ggcd(pc,content(a2));

  tetpil=avma;

  return gerepile(av,tetpil,gdiv(pm,pc));
}

/***********************************************************************/
/***********  Normalized integral basis                   **************/
/***********************************************************************/

GEN
nbasis(GEN ibas,GEN pd)
{
  long n,j,k,m,av=avma,tetpil;
  GEN a,unmodpd;
  
  unmodpd=gmodulcp(gun,pd);
  n=lg(ibas)-1;
  
  a=cgetg(n+1,19);
  m=lgef((GEN)ibas[1])-2;
  for(k=1;k<=n;k++)
  {
    m=lgef((GEN)ibas[k])-2;
    a[k]=lgetg(n+1,18);
    for(j=1;j<=m;j++)
      coeff(a,j,k)=(long)((GEN)ibas[k])[j+1];
    for(j=m+1;j<=n;j++)
      coeff(a,j,k)=zero;
  }
  a=hnfmodid(a,pd);
  tetpil=avma;
  return gerepile(av,tetpil,gdiv(a,pd));
}

/***********************************************************************/
/**************        Pick best divisor of chi           **************/
/***********************************************************************/

GEN
bestnu(GEN w)
{
  long r,j,av=avma,tetpil;
  GEN g,h;
  
  r=lg((GEN)w[1])-1;
  g=polun[0];
  
  for(j=1;j<=r;j++)
  {
    h=(GEN)coeff(w,j,1);
    if (lgef(h)>lgef(g)) g=h;
  }
  tetpil=avma;
  return gerepile(av,tetpil,lift(g));
}

/**********************************************************************/
/******                        bezout etendu                    *******/
/******                  Return d=pgcd(a,b)  and  &u,&v         *******/
/**********************************************************************/

long
cbezout(long a,long b,long *u,long *v)
{
  long d,v3,v1,q,r,t;
  
  (*u)=1;
  d=a;
  v1=0;
  v3=b;
  
  while (1)
  {
    if (v3==0) 
    { 
      (*v)=(d-a*(*u))/b;
      return d;
    }
    q=d/v3;
    r=d%v3;
    t=(*u)-v1*q;
    (*u)=v1;
    d=labs(v3);
    v1=t;
    v3=r;
  }
}

long
clcm(long a,long b)
{
  long d,r,v1,q;
  
  d=a;
  r=b;

  while (1)
  {
    if (r==0) 
    { 
      return (a*b)/d;
    }
    v1=r;
    q=d/r;
    r=d%r;
    d=labs(v1);
  }
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                 ALGORITHME DE BUCHMANN-LENSTRA                  */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



GEN
pradical(GEN nf, GEN p)
              
/* Calcule une F_p base du p-radical de Z_K,i.e une
 matrice N fois r (r=dimension sur F_p du p-radical) */

{
  long av=avma,tetpil,j,k,N=lgef((GEN)nf[1])-3;
  GEN p1,unmodp,zmodp,x,m;

  p1=p;while(cmpis(p1,N)<0) p1=mulii(p1,p);
  m=cgetg(N+1,19);
  unmodp=gmodulcp(gun,p);zmodp=gmul(gzero,unmodp);
  for(k=1;k<=N;k++)
  {
    x=cgetg(N+1,18);
    for(j=1;j<=N;j++) x[j]=(long)zmodp;
    x[k]=(long)unmodp;
    m[k]=(long)element_pow(nf,x,p1);
  }
  tetpil=avma;return gerepile(av,tetpil,gmul(unmodp,ker(m)));
}
      
GEN
pol_min(GEN alpha, GEN nf, GEN p, GEN algebre, GEN algebre1)
                                     
/* Calcule le polynome minimal de alpha dans algebre (polynome a
 coefficients dans F_p) */

{
  long av=avma,tetpil,i,j,N=lgef((GEN)nf[1])-3,k=N-lg(algebre1)+1;
  GEN puiss,puiss2,noyau,unmodp,zmodp,vecteur;

  unmodp=gmodulcp(gun,p);zmodp=gmul(unmodp,gzero);
  vecteur=gmul(unmodp,(GEN)(idmat(N)[k+1]));
  puiss=cgetg(N+2,19);puiss[1]=(long)vecteur;
  for(i=1;i<=N;i++)
    puiss[i+1]=(long)(inverseimage(algebre,element_pow(nf,alpha,stoi(i))));
  puiss2=cgetg(N+2,19);
  for(i=1;i<=N+1;i++)
  {
    puiss2[i]=lgetg(N-k+1,18);
    for(j=1;j<=N-k;j++)
      coeff(puiss2,j,i)=coeff(puiss,k+j,i);
  }
  noyau=gmul(ker(puiss2),unmodp);
  tetpil=avma;return gerepile(av,tetpil,gtopolyrev((GEN)noyau[1],0));
}

GEN
eval_pol(GEN nf, GEN pol, GEN alpha, GEN p, GEN algebre, GEN algebre1)
                                         
/* Evalue le polynome pol en alpha,element de nf */

{
  long av=avma,tetpil,i;
  long N=lgef((GEN)nf[1])-3,k=N-lg(algebre1)+1,lx=lgef(pol)-3;
  GEN res,valeur,unmodp,zmodp;

  res=cgetg(N+1,18);
  unmodp=gmodulcp(gun,p);zmodp=gmul(unmodp,gzero);
  for(i=2;i<=N;i++) res[i]=(long)zmodp;
  res[1]=lmul(unmodp,(GEN)pol[lx+2]);
  for(i=lx+1;i>=2;i--)
  {
    res=element_mul(nf,alpha,res);res[1]=ladd((GEN)res[1],(GEN)pol[i]);
  }
  res=inverseimage(algebre,res);
  valeur=cgetg(N-k+1,18);
  for(i=1;i<=N-k;i++) valeur[i]=res[k+i];
  tetpil=avma;return gerepile(av,tetpil,gmul(algebre1,valeur));
}

GEN
kerlens2(GEN x, GEN pgen)
{
  long i,j,k,t,nbc,nbl,av,av1;
  GEN a,c,l,d,y,q;

  av=avma;
  a=gmul(x,gmodulcp(gun,pgen));
  nbl=nbc=lg(x)-1;
  c=cgetg(nbl+1,17);l=cgetg(nbc+1,17);
  d=cgetg(nbc+1,17);
  for(i=1;i<=nbl;i++) c[i]=0;
  k=1;t=1;
  while((t<=nbl)&&(k<=nbc))
  {
    for(j=1;j<k;j++)
      for(i=1;i<=nbl;i++)
	if(i!=l[j])
	  coeff(a,i,k)=lsub(gmul((GEN)d[j],gcoeff(a,i,k)),gmul(gcoeff(a,l[j],k),gcoeff(a,i,j)));
    t=1;while((t<=nbl)&&((c[t])||gcmp0(gcoeff(a,t,k)))) t++;
    if (t<=nbl) {d[k]=coeff(a,t,k);c[t]=k;l[k++]=t;}
  }
  if(k>nbc) err(kerlenser);
  y=cgetg(nbc+1,18);
  y[1]=(k>1)?(long)coeff(a,l[1],k):un;
  for(q=gun,j=2;j<k;j++)
  {
    q=gmul(q,(GEN)d[j-1]);
    y[j]=lmul(gcoeff(a,l[j],k),q);
  }
  if(k>1) y[k]=lneg(gmul(q,(GEN)d[k-1]));
  for(j=k+1;j<=nbc;j++) y[j]=zero;
  av1=avma;
  return gerepile(av,av1,gcopy(lift(y))); 
}

GEN
kerlens(GEN x, GEN pgen)
{
  long i,j,k,t,nbc,nbl,p,q,*c,*l,*d,**a;
  GEN y;

  if(cmpis(pgen,(MAXHALFULONG>>1))>0)
    return kerlens2(x,pgen);
/* ici p<=(MAXHALFULONG>>1) ==> simple precision (long de C) */
  p=itos(pgen);
  nbl=nbc=lg(x)-1;
  a=(long**)newbloc(nbc+1);
  for(j=1;j<=nbc;j++)
  {
    c=a[j]=newbloc(nbl+1);
    for(i=1;i<=nbl;i++) c[i]=itos(modis(gcoeff(x,i,j),p));
  }
  c=newbloc(nbl+1);
  l=newbloc(nbc+1);
  d=newbloc(nbc+1);
  for(i=1;i<=nbl;i++) c[i]=0;
  k=1;t=1;
  while((t<=nbl)&&(k<=nbc))
  {
    for(j=1;j<k;j++)
      for(i=1;i<=nbl;i++)
	if(i!=l[j]) a[k][i]=(d[j]*a[k][i]-a[j][i]*a[k][l[j]]) % p;
    t=1;while((t<=nbl)&&((c[t])||(!a[k][t]))) t++;
    if (t<=nbl) {d[k]=a[k][t];c[t]=k;l[k++]=t;}
  }
  if(k>nbc) err(kerlenser);
  y=cgetg(nbc+1,18);
  t=(k>1) ? a[k][l[1]]:1;
  y[1]=(t>0)? lstoi(t):lstoi(t+p);
  for(q=1,j=2;j<k;j++)
  {
    q=(q*d[j-1])%p;
    t=(a[k][l[j]]*q)%p;
    y[j]=(t>0)? lstoi(t):lstoi(t+p);
  }
  if(k>1)
  {
    t=(q*d[k-1])%p;
    y[k]=(t>0)? lstoi(p-t):lstoi(-t);
  }
  for(j=k+1;j<=nbc;j++) y[j]=zero;
  killbloc(c);killbloc(l);killbloc(d);
  for(j=1;j<=nbc;j++) killbloc(a[j]);killbloc((GEN)a);
  return y;
}
  
GEN
lens(GEN nf, GEN p, GEN a)
                
/* Calcule la constante de lenstra de l'ideal p.Z_K+a.Z_K ou a est un
vecteur sur la base d'entiers */

{
  long av=avma,tetpil,N=lgef((GEN)nf[1])-3,j;
  GEN mat;

  mat=cgetg(N+1,19);for(j=1;j<=N;j++) mat[j]=(long)element_mulid(nf,a,j);
  tetpil=avma;return gerepile(av,tetpil,kerlens(mat,p));
}

GEN
two_elt(GEN nf, GEN p, GEN ideal)
                    
/* Recoit un ideal (mod p) et calcule une representation a deux
 elements (ideal non egal a Z_K) */

{
  long av=avma,av1,tetpil,N=lgef((GEN)nf[1])-3,m,r,i,j,fl;
  GEN beta,alpha,lambda,norme,pf;

  m=lg(ideal)-1;
  if(!m)
  {alpha=cgetg(N+1,18);for(i=1;i<=N;i++) alpha[i]=zero;return alpha;}
  beta=gmodulcp(gmul((GEN)nf[7],lift(ideal)),(GEN)nf[1]);
  pf=gpuigs(p,N-m);fl=r=1;
  for(i=1;(i<=m)&&fl;i++)
  {
    alpha=(GEN)beta[i];norme=gnorm(alpha);
    if(signe(modii(divii(norme,pf),p))) fl=0;
    else
    {
      alpha=gadd(alpha,p);norme=gnorm(alpha);
      if(signe(modii(divii(norme,pf),p))) fl=0;
    }
  }
  if(fl)
  {
    lambda=cgeti(m+1);
    av1=avma;
    for(i=1;i<=m;i++) lambda[i]=r;
    do
    {
      avma=av1;
      alpha=gmodulcp(gzero,(GEN)nf[1]);
      for(i=1;i<=m;i++) alpha=gadd(alpha,gmulsg(lambda[i],(GEN)beta[i]));
      norme=gnorm(alpha);
      if(signe(modii(divii(norme,pf),p))) fl=0;
      else
      {
	alpha=gadd(alpha,p);norme=gnorm(alpha);
	if(signe(modii(divii(norme,pf),p))) fl=0;
      }
      if(fl)
      {
	for(j=m;(lambda[j]+r)==0;j--);
	lambda[j]--;
	for(i=j+1;i<=m;i++) lambda[i]=r;
	for(j=1;(j<m)&&(!lambda[j]);j++);
	if(!lambda[j])
	{
	  r++;for(i=1;i<=m;i++) lambda[i]=r;
	  if(cmpis(p,(r<<1))<0) err(talker,"bug in two_elt");
	}
      }
    }
    while(fl);
  }
  alpha=lift(alpha);beta=cgetg(N+1,18);
  for(i=1;i<=N;i++) beta[i]=(long)truecoeff(alpha,i-1);
  alpha=gmul((GEN)nf[8],beta);
  alpha=gmul(gmodulcp(gun,p),alpha);
  alpha=centerlift(alpha);
  if(!signe(modii(divii(subres(gmul((GEN)nf[7],alpha),(GEN)nf[1]),pf),p)))
    alpha[1]=(long)gadd((GEN)alpha[1],p);
  tetpil=avma;return gerepile(av,tetpil,gcopy(alpha));
}


GEN
primedec(GEN nf, GEN p)
              
/* Recoit un corps de nombres nf et un premier p,ressort une liste 
 des ideaux premiers au dessus de p dans le format vu plus haut, dans
l'ordre croissant des degres residuels */

{
  long av=avma,tetpil,i,j,k,v,kbar,l3,np,c,i1,i2,indice,N,lp;
  GEN f,ff,list,list2,ip,elementh,hensemble;
  GEN algebre,algebre1,b,b2,mat1,mat2;
  GEN alpha,beta,beta1,p1,p2,p3,unmodp,zmodp,vecteur,pol,f1,polg;
  GEN T,pidmat,ppuin;

  if(DEBUGLEVEL>=3) timer2();
  nf=checknf(nf);  
  T=(GEN)nf[1];N=lgef(T)-3;
  if(signe(modii((GEN)nf[4],p)))
  {
    f=centerlift(ff=factmod(T,p));np=lg((GEN)f[1]);
    if(DEBUGLEVEL>=6) {fprintferr("temps factmod: ");fprintferr("%ld\n",timer2());flusherr();}
    list=cgetg(np,17);
    for(i=1;i<np;i++)
    {
      p1=(GEN)(list[i]=lgetg(6,17));
      p1[1]=(long)p;p3=gcoeff(f,i,1);l3=lgef(p3)-1;
      p2=cgetg(N+1,18);
      if(l3==(N+2)) 
      {
	p1[2]=(long)p2;p1[3]=un;p1[4]=lstoi(N);
	p3=cgetg(N+1,18);p1[5]=(long)p3;p3[1]=un;
	p2[1]=(long)p;for(j=2;j<=N;j++) p3[j]=p2[j]=zero;
      }
      else
      {
	v=ggval(subres(p3,T),p)/(l3-2);if(v>1) p3[2]=ladd((GEN)p3[2],p);
	for(j=1;j<l3;j++) p2[j]=p3[j+1];
	for(j=l3;j<=N;j++) p2[j]=zero;
	p1[2]=lmul((GEN)nf[8],p2);
	p1[3]=(long)coeff(f,i,2);
	p1[4]=lstoi(l3-2);
	p3=gdiv(T,gcoeff(ff,i,1));l3=lgef(p3)-1;
	for(j=1;j<l3;j++) p2[j]=p3[j+1];
	for(j=l3;j<=N;j++) p2[j]=zero;
	p1[5]=(long)centerlift(gmul((GEN)nf[8],p2));
      }
    }
    if(DEBUGLEVEL>=6) {fprintferr("temps lens etc...: ");fprintferr("%ld\n",timer2());flusherr();}      
    p1=stoi(4);tetpil=avma;return gerepile(av,tetpil,vecsort(list,p1));
  }
  else 
  {
    unmodp=gmodulcp(gun,p);zmodp=gmodulcp(gzero,p);
    list=cgetg(N+1,17);for(i=1;i<=N;i++) list[i]=lgetg(6,17);
    f=centerlift(ff=factmod(T,p));indice=0;
    if(DEBUGLEVEL>=6) {fprintferr("temps factmod: ");fprintferr("%ld\n",timer2());flusherr();}
    f1=(GEN)ff[1];np=lg(f1);polg=(GEN)f1[1];
    for(i=2;i<np;i++) polg=gmul(polg,(GEN)f1[i]);
    polg=gmul(unmodp,gdiv(gsub(gmul(lift(polg),lift(gdiv(T,polg))),T),p));
    for(i=1;i<np;i++)
    {
      if((gcmp1(gcoeff(f,i,2)))||(!gdivise(polg,(GEN)f1[i])))
      {
	indice++;p1=(GEN)list[indice];
	p1[1]=(long)p;p3=gcoeff(f,i,1);l3=lgef(p3)-1;
	p2=cgetg(N+1,18);
	if(l3==(N+2)) 
	{
	  p1[2]=(long)p2;p1[3]=un;p1[4]=lstoi(N);
	  p3=cgetg(N+1,18);p1[5]=(long)p3;p3[1]=un;
	  p2[1]=(long)p;for(j=2;j<=N;j++) p3[j]=p2[j]=zero;
	  p3=gmul(unmodp,polun[varn(T)]);
	}
	else
	{
	  v=ggval(subres(p3,T),p)/(l3-2);
	  if(v>1) p3[2]=ladd((GEN)p3[2],p);
	  for(j=1;j<l3;j++) p2[j]=p3[j+1];
	  for(j=l3;j<=N;j++) p2[j]=zero;
	  p1[2]=lmul((GEN)nf[8],p2);
	  p1[3]=(long)coeff(f,i,2);
	  p1[4]=lstoi(l3-2);
	  p3=gdiv(T,gcoeff(ff,i,1));l3=lgef(p3)-1;
	  for(j=1;j<l3;j++) p2[j]=p3[j+1];
	  for(j=l3;j<=N;j++) p2[j]=zero;
	  p1[5]=(long)centerlift(gmul((GEN)nf[8],p2));
	}
	beta=(indice==1)?p3:gdiv(beta,gcoeff(ff,i,1));
      }
    }
    if(DEBUGLEVEL>=3) {fprintferr("temps %ld facteurs non ramifies : ",indice);fprintferr("%ld\n",timer2());flusherr();}      
    ip=pradical(nf,p);
    if(DEBUGLEVEL>=3) {fprintferr("temps pradical: ");fprintferr("%ld\n",timer2());flusherr();}
    if(indice)
    {
      if(typ(beta)!=10) err(talker,"bugbeta in primedec");
      p2=cgetg(N+1,18);l3=lgef(beta)-1;
      for(j=1;j<l3;j++) p2[j]=beta[j+1];
      for(j=l3;j<=N;j++) p2[j]=zero;
      beta=gmul((GEN)nf[8],p2);beta1=lift(beta);
      lp=lg(ip)-1;p1=cgetg(lp+lp+N+1,19);
      for(i=1;i<=N;i++) p1[i]=(long)element_mulid(nf,beta,i);
      for(;i<=N+lp;i++) 
      {
	p2=lift((GEN)ip[i-N]);p1[i]=ldiv(element_mul(nf,p2,beta1),p);
	p1[i+lp]=(long)p2;
      }
      ip=image(gmul(unmodp,p1));
    }
    pidmat=cgetg(N+1,18);vecteur=cgetg(N+1,18);
    for(i1=2;i1<=N;i1++) {pidmat[i1]=zero;vecteur[i1]=(long)zmodp;}
    pidmat[1]=(long)p;vecteur[1]=(long)unmodp;
    ppuin=gpuigs(p,N);c=0;
    hensemble=cgetg(N+1,17); 
    if(lg(ip)<N+1) {c=1;hensemble[1]=(long)ip;}
    while(c)
    {  
      elementh=(GEN)(hensemble[c]);k=lg(elementh)-1;kbar=N-k;
      algebre=gmul(unmodp,suppl(concat(elementh,vecteur)));
      algebre1=cgetg(kbar+1,19);
      for(i1=1;i1<=kbar;i1++) algebre1[i1]=algebre[i1+k];
      b=cgetg(kbar+1,19);
      for(i1=1;i1<=kbar;i1++)
	b[i1]=lsub(element_pow(nf,(GEN)algebre1[i1],p),(GEN)algebre1[i1]);
      b2=inverseimage(algebre,b);
      mat1=cgetg(kbar+1,19);
      for(i1=1;i1<=kbar;i1++) mat1[i1]=lgetg(kbar+1,18);
      for(i1=1;i1<=kbar;i1++)
	for(i2=1;i2<=kbar;i2++)
	  coeff(mat1,i2,i1)=coeff(b2,k+i2,i1);
      mat2=cgetg(k+N+1,19);
      for(i1=1;i1<=k;i1++) mat2[i1]=elementh[i1];
      mat1=gmul(unmodp,ker(mat1));
      if(lg(mat1)>2)
      {
	alpha=gmul(algebre1,(GEN)mat1[2]);
	pol=pol_min(alpha,nf,p,algebre,algebre1);
	setvarn(pol,0);
	p1=(GEN)factmod(pol,p)[1];
	for(i=1;i<lg(p1);i++)
	{
	  beta1=eval_pol(nf,(GEN)p1[i],alpha,p,algebre,algebre1);
	  for(i1=1;i1<=N;i1++)
	    mat2[k+i1]=(long)element_mulid(nf,beta1,i1);
	  hensemble[c]=(long)(image(mat2));c++;
	}
	c--;
	if(DEBUGLEVEL>=3) {fprintferr("temps hensemble[%ld]: ",c);fprintferr("%ld\n",timer2());flusherr();}
      }
      else
      {
	indice++;p1=(GEN)list[indice];
	p1[1]=(long)p;p1[4]=lstoi(kbar);
	p1[2]=(long)two_elt(nf,p,elementh);
	p1[5]=(long)lens(nf,p,(GEN)p1[2]);
	p1[3]=lstoi(element_val2(nf,pidmat,ppuin,p1));
	c--;
	if(DEBUGLEVEL>=3) {fprintferr("temps hensemble[%ld]: ",c);fprintferr("%ld\n",timer2());flusherr();}
      }
    }
    list2=cgetg(indice+1,17);
    for(i=1;i<=indice;i++) list2[i]=list[i];
    p1=stoi(4);tetpil=avma;
    return gerepile(av,tetpil,vecsort(list2,p1));
  }
}

long
idealval(GEN nf, GEN ix, GEN vp)
                  
/* recoit un ideal ix et un ideal premier vp dans le format
donne par primedec et calcule la valuation de ix en vp */

{
  long N,v,vd,w,av=avma,i,j,bo;
  GEN mat,x,d,bp,p,p1,r,denx;

  nf=checknf(nf);  
  if((typ(vp)!=17)||(lg(vp)!=6)) err(idealer3);
  if((typ(ix)<=10)||(typ(ix)==18)) return element_val(nf,ix,vp);
  N=lgef((GEN)nf[1])-3;p=(GEN)vp[1];
  if((typ(ix)==17)&&(lg(ix)==3)) x=(GEN)ix[1]; else x=ix;
  if(typ(x)!=19) err(idealer2);
  denx=denom(x);
  if(!gcmp1(denx)) x=gmul(denx,x);
  if(lg((GEN)x[1])!=(N+1)) err(idealer4);
  if(lg(x)!=(N+1)) x=idealmul(nf,x,idmat(N));
  for(d=gun,i=1;i<=N;i++) d=mulii(d,(GEN)coeff(x,i,i));
  v=ggval(d,p);vd=ggval(denx,p);
  if(!v) return -vd*itos((GEN)vp[3]);
  bo=0;w=0;bp=(GEN)vp[5];
  do
  {
    if(w) {for(i=1;i<=N;i++) mat[i]=(long)element_muli(nf,(GEN)mat[i],bp);}
    else 
    {
      mat=cgetg(N+1,19);
      for(i=1;i<=N;i++) mat[i]=(long)element_mulh(nf,i,N,(GEN)x[i],bp);
    }
    if(divise(gcoeff(mat,N,N),p))
    {
      for(j=1;j<=N;j++)
	for(i=1;i<=N;i++)
	{
	  p1=dvmdii(gcoeff(mat,i,j),p,&r);
	  if(signe(r)) goto labeliv; else coeff(mat,i,j)=(long)p1;
	}
      w++;
    }
    else bo=1;
  }
  while((bo==0)&&(w<v));
  labeliv:
  avma=av;return w-vd*itos((GEN)vp[3]);
}

/*******************************************************************
                           ROUND 2 relatif

  Entree:   nf = corps de base K dans le format initalg.
            x polynome unitaire a coefficients dans Z_K de deg n
	    definissant une extension relative L=K(theta);
	    La variable de x doit etre de numero strictement
	    inferieur a celle de nf[1].
  Sortie:   retourne une pseudo-base [A,I] de Z_L, ou A est une matrice
            nxn a coefficients dans nf[1] sous forme HNF et I un vecteur
	    d'ideaux a n composantes
            
*******************************************************************/

GEN rnfordmax(GEN nf, GEN pol, GEN pr, GEN unnf, GEN zeronf, GEN id, GEN psid, GEN powbasis);
GEN rnfjoinmodules(GEN nf, GEN x, GEN y);

void
checkbnf(GEN bnf)
{
  if(typ(bnf)!=17) err(idealer1);
  if(lg(bnf)!=9)
  {
    if((lg(bnf)==10)&&(typ((GEN)bnf[1])==10))
      err(talker,"please apply buchinit first in rnf function");
    else err(idealer1);
  }
}

GEN
checknf(GEN nf)
{
  if(typ(nf)!=17) err(idealer1);
  if(lg(nf)==10) return nf;
  if(lg(nf)==9) return checknf((GEN)nf[7]);
  err(idealer1);return gnil;
}
    
GEN
rnfround2all(GEN nf, GEN pol, long all)
{
  long av=avma,tetpil,i,j,n,N,nbidp,ipr,vpol,cpt;
  GEN p1,p2,p3,p4,polnf,fact,list,ep,unnf,zeronf,id,A,I,W,pseudo,y,discpol,psid,powbasis,d,D;

  nf=checknf(nf);
  polnf=(GEN)nf[1];vpol=varn(pol);
  if((typ(pol)!=10)||(vpol>=varn(polnf)))
    err(talker,"incorrect polynomial in relativeround2");
  N=lgef(polnf)-3;n=lgef(pol)-3;discpol=discsr(pol);
  fact=idealfactor(nf,discpol);list=(GEN)fact[1];ep=(GEN)fact[2];
  nbidp=lg(list)-1;
  if(DEBUGLEVEL>1)
  {
    fprintferr(" Ideaux a considerer :\n");
    for(i=1;i<=nbidp;i++)
    {
      j=itos((GEN)ep[i]);
      if(j>1) {bruterr((GEN)list[i],'g',-1);fprintferr("^%ld\n",j);}
    }
    flusherr();
  }
  unnf=cgetg(N+1,18);for(i=2;i<=N;i++) unnf[i]=zero;unnf[1]=un;
  zeronf=cgetg(N+1,18);for(i=1;i<=N;i++) zeronf[i]=zero;
  id=idmat(N);A=cgetg(n+1,19);
  for(j=1;j<=n;j++)
  {
    p1=cgetg(n+1,18);A[j]=(long)p1;
    for(i=1;i<=n;i++) p1[i]=(i==j)?(long)unnf:(long)zeronf;
  }
  I=cgetg(n+1,17);for(i=1;i<=n;i++) I[i]=(long)id;
  pseudo=cgetg(3,17);pseudo[1]=(long)A;pseudo[2]=(long)I;
  psid=gcopy(pseudo);
  powbasis=cgetg(n+1,17);powbasis[1]=(long)polun[vpol];
  for(i=2;i<=n;i++) powbasis[i]=lmul((GEN)powbasis[i-1],polx[vpol]);
  cpt=0;
  for(ipr=1;ipr<=nbidp;ipr++)
  {
    if(!gcmp1((GEN)ep[ipr]))
    {
      y=rnfordmax(nf,pol,(GEN)list[ipr],unnf,zeronf,id,psid,powbasis);
      if(cpt) pseudo=rnfjoinmodules(nf,pseudo,y);
      else {cpt++;pseudo=y;}
    }
  }
  W=gmodulcp(gmul(powbasis,basistoalg(nf,(GEN)pseudo[1])),pol);
  I=(GEN)pseudo[2];
  p2=cgetg(n+1,19);for(j=1;j<=n;j++) p2[j]=lgetg(n+1,18);
  for(j=1;j<=n;j++) for(i=j;i<=n;i++)
  {
    coeff(p2,i,j)=(long)trace(gmul((GEN)W[i],(GEN)W[j]));
    if(i!=j) coeff(p2,j,i)=coeff(p2,i,j);
  }
  d=algtobasis(nf,det(p2));
  for(i=1;(i<=n)&&vecegal((GEN)I[i],id);i++);
  if(i>n) D=id;
  else
  {
    D=(GEN)I[i];
    for(i++;i<=n;i++)
      if(!vecegal((GEN)I[i],id)) D=idealmul(nf,D,(GEN)I[i]);
    D=idealmul(nf,D,D);
  }
  p2=gtomat(d);
  p3=auxdecomp(content(d),0);
  p4=gun;
  for(i=1;i<lg((GEN)p3[1]);i++)
    p4=gmul(p4,gpuigs(gcoeff(p3,i,1),(itos(gcoeff(p3,i,2)))>>1));
  p4=gsqr(p4);
  tetpil=avma;
  if(all)
  {
    p1=cgetg(5,17);p1[1]=lcopy((GEN)pseudo[1]);
    p1[2]=lcopy(I);p1[3]=(long)idealmul(nf,D,p2);p1[4]=(long)gdiv(d,p4);
  }
  else
  {
    p1=cgetg(3,17);p1[1]=(long)idealmul(nf,D,p2);p1[2]=(long)gdiv(d,p4);
  }
  tetpil=avma;return gerepile(av,tetpil,p1);
}

GEN
rnfpseudobasis(GEN nf, GEN pol)
{
  return rnfround2all(nf,pol,1);
}

GEN
rnfdiscf(GEN nf, GEN pol)
{
  return rnfround2all(nf,pol,0);
}

GEN
nfreducemodpr(GEN nf, GEN x, GEN prhall)
{
 /* a usage interne, pas de gestion de pile */
  
  long N=lg(x)-1,i,flx,v;
  GEN p,prh,den;

  flx=1;for(i=1;(i<=N)&&flx;i++) flx=(typ((GEN)x[i])!=3);
  if(!flx) x=lift(x);prh=(GEN)prhall[1];
  p=gcoeff(prh,1,1);
  if(gcmp1(p)) err(talker,"bug in reducemodpr");
  den=denom(x);
  if(!gcmp1(den))
  {
    v=ggval(den,p);
    if(v) x=element_mul(nf,x,element_pow(nf,(GEN)prhall[2],stoi(v)));
  }
  x=gmod(x,p);
  for(i=N;i>=1;i--)
    if(gcmp1(gcoeff(prh,i,i))) x=gsub(x,gmul((GEN)x[i],(GEN)prh[i]));
  return gmul(gmodulcp(gun,p),x);
}

GEN
rnfelement_mulmod(GEN nf, GEN multab, GEN zeronf, GEN unnf, GEN x, GEN y, GEN prhall)
{
  /* a usage interne, pas de gestion de pile : x et y sont des vecteurs dont
     les coefficients sont les composantes sur nf[7] ; avec reduction mod pr sauf
     si prhall=gzero */

  long i,j,k,n;
  GEN p1,p2,z,s;

  n=lg(x)-1;x=lift(x);y=lift(y);z=cgetg(n+1,18);
  for(k=1;k<=n;k++)
  {
    s=zeronf;
    for(i=1;i<=n;i++)
    {
      for(j=1;j<=n;j++)
      {
	p2=gcoeff(multab,k,(i-1)*n+j);
	if(!gcmp0(p2))
	{
	  p1=element_mul(nf,(GEN)x[i],(GEN)y[j]);
	  if(vecegal(p2,unnf)) s=gadd(s,p1);
	  else s=gadd(s,element_mul(nf,p1,p2));
	}
      }
    }
    if(!gcmp0(prhall)) z[k]=(long)nfreducemodpr(nf,s,prhall);
    else z[k]=(long)s;
  }
  return z;
}

GEN
rnfelement_sqrmod(GEN nf, GEN multab, GEN zeronf, GEN unnf, GEN x, GEN prhall)
  /* a usage interne, pas de gestion de pile : x est un vecteur dont
     les coefficients sont les composantes sur nf[7] */
              
{
  long i,j,k,n;
  GEN p1,p2,z,s;

  n=lg(x)-1;x=lift(x);
  z=cgetg(n+1,18);
  for(k=1;k<=n;k++)
  {
    s=zeronf;
    for(i=1;i<=n;i++)
    {
      if(!gcmp0(p2=gcoeff(multab,k,(i-1)*n+i)))
      {
	p1=element_sqr(nf,(GEN)x[i]);
	if(vecegal(p2,unnf)) s=gadd(s,p1);
	else s=gadd(s,element_mul(nf,p1,p2));
      }
    }
    for(i=1;i<=n;i++) for(j=i+1;j<=n;j++)
    {
      if(!gcmp0(p2=gcoeff(multab,k,(i-1)*n+j)))
      {
	p1=gmul2n(element_mul(nf,(GEN)x[i],(GEN)x[j]),1);
	if(vecegal(p2,unnf)) s=gadd(s,p1);
	else s=gadd(s,element_mul(nf,p1,p2));
      }
    }
    if(!gcmp0(prhall)) z[k]=(long)nfreducemodpr(nf,s,prhall);
    else z[k]=(long)s;
  }
  return z;
}

GEN
rnfelement_powmod(GEN nf, GEN multab, GEN zeronf, GEN unnf, GEN x, GEN k, GEN prhall)
                
/* Calcule x^k mod pr dans l'extension . */
 

{
  long i,f,n,av=avma,tetpil;
  GEN k1,y,z;

  n=lg(x)-1;k1=k;z=x;f=1;y=cgetg(n+1,18);
  for(i=2;i<=n;i++) y[i]=(long)zeronf;y[1]=(long)unnf;
  while(f)
  {
    if(mpodd(k1)) y=rnfelement_mulmod(nf,multab,zeronf,unnf,z,y,prhall);
    k1=shifti(k1,-1);f=signe(k1);
    if(f) z=rnfelement_sqrmod(nf,multab,zeronf,unnf,z,prhall);
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(y));
}

GEN
rnfordmax(GEN nf, GEN pol, GEN pr, GEN unnf, GEN zeronf, GEN id, GEN psid, GEN powbasis)
{
  long av=avma,tetpil,av1,av2,lim,dec,i,j,k,n,N,notfinished=1,v1,v2,vpol,m,fl;
  GEN polnf,p,q,q1,prh,prhall,A,Aa,Aaa,A1,den,I,R,p1,p2,p3,multab,Aainv,alphalist;
  GEN pip,baseIp,baseOp,alpha,matprod,alphainv,matC,matG,matV,vecpro,matH;
  GEN neworder,matId,H,Hid,alphalistinv,epr,betae;
  
  long cmpt;

  polnf=(GEN)nf[1];N=lgef(polnf)-3;n=lgef(pol)-3;vpol=varn(pol);
  p=(GEN)pr[1];q=gpui(p,(GEN)pr[4],0);pip=(GEN)pr[2];
  q1=q;while(cmpis(q1,n)<0) q1=mulii(q1,q);
  prh=idealmulprime(nf,id,pr);
  epr=(GEN)pr[3];betae=gdiv(element_pow(nf,(GEN)pr[5],epr),gpui(p,addis(epr,-1),0));
  p1=cgetg(2,19);p1[1]=(long)betae;
  p1=idealadd(nf,gmul(p,id),idealmul(nf,p1,id));
  prhall=cgetg(3,17);prhall[1]=(long)prh;
  prhall[2]=idealaddone(nf,pr,p1)[2];
/* ceci contient un alpha congru a 1 mod pr et a 0 mod q^{e_q} pour tous
   les autres ideaux premiers q au dessus de p */
  A=(GEN)psid[1];I=(GEN)psid[2];matId=(GEN)psid[1];

  cmpt=0;
  if(DEBUGLEVEL>1)
  {fprintferr("\n Ideal traite : ");outerr(pr);flusherr();}

  lim=(bot+avma)>>1;av1=avma;
  while(notfinished)
  {
    if(cmpt&&(avma<lim))
    {
      tetpil=avma;A=gcopy(A);I=gcopy(I);av2=avma;
      dec=lpile(av1,tetpil,0)>>TWOPOTBYTES_IN_LONG;
      if(adecaler(A,tetpil,av2)) A+=dec;
      if(adecaler(I,tetpil,av2)) I+=dec;
    }
    cmpt++;
    if(DEBUGLEVEL>1)
    {fprintferr("\n\n     %ld eme passe \n",cmpt);flusherr();}
    
    alphalist=cgetg(n+1,17);alphalistinv=cgetg(n+1,17);
    for(i=1;i<=n;i++)
    {
      if(vecegal((GEN)I[i],id)) alphalist[i]=alphalistinv[i]=(long)unnf;
      else
      {
	den=denom((GEN)I[i]);p1=gcmp1(den)? (GEN)I[i]:gmul(den,(GEN)I[i]);
	p1=gdiv(ideal_two_elt(nf,p1),den);
	v1=(gcmp0((GEN)p1[1]))?EXP220:element_val(nf,p2=gmul((GEN)p1[1],unnf),pr);
	v2=(gcmp0((GEN)p1[2]))?EXP220:element_val(nf,(GEN)p1[2],pr);
	if(v1>v2) p2=(GEN)p1[2];
	alphalist[i]=(long)p2;alphalistinv[i]=(long)element_inv(nf,p2);
      }
    }
    A1=cgetg(n+1,19);
    for(j=1;j<=n;j++)
    {
      p1=cgetg(n+1,18);A1[j]=(long)p1;
      for(i=1;i<=n;i++) p1[i]=(long)element_mul(nf,gcoeff(A,i,j),(GEN)alphalist[j]);
    }
    Aa=basistoalg(nf,A1);Aainv=ginv(Aa);
    Aaa=gmodulcp(gmul(powbasis,Aa),pol);
    multab=cgetg(n*n+1,19);for(j=1;j<=n*n;j++) multab[j]=lgetg(n+1,18);
    for(i=1;i<=n;i++) for(j=i;j<=n;j++)
    {
      p1=gmul((GEN)Aaa[i],(GEN)Aaa[j]);p2=cgetg(n+1,18);
      if((typ(p1)==9)&&(varn((GEN)p1[1])==vpol)) p1=(GEN)p1[2];
      if((gcmp0(p1))||(typ(p1)<=9)||((typ(p1)==10)&&(varn(p1)>vpol)))
      {p2[1]=(long)p1;for(k=2;k<=n;k++) p2[k]=(long)gmodulcp(gzero,polnf);}
      else for(k=1;k<=n;k++) p2[k]=(long)truecoeff(p1,k-1);
      p3=algtobasis(nf,gmul(Aainv,p2));
      for(k=1;k<=n;k++)
      {coeff(multab,k,(i-1)*n+j)=(long)p3[k];coeff(multab,k,(j-1)*n+i)=(long)p3[k];}
    }
    R=cgetg(n+1,19);
    for(j=1;j<=n;j++)
      R[j]=(long)rnfelement_powmod(nf,multab,zeronf,unnf,(GEN)matId[j],q1,prhall);
    baseIp=nfker(nf,R,prhall);
    baseOp=nfsuppl(nf,baseIp,n,prhall);
    alpha=cgetg(n+1,19);
    for(j=1;j<lg(baseIp);j++) alpha[j]=(long)lift((GEN)baseOp[j]);
    for(;j<=n;j++)
    {
      p1=cgetg(n+1,18);alpha[j]=(long)p1;
      for(i=1;i<=n;i++) p1[i]=(long)element_mul(nf,pip,lift(gcoeff(baseOp,i,j)));
    }
    matprod=cgetg(n+1,19);
    for(j=1;j<=n;j++)
    {
      p1=cgetg(n+1,18);matprod[j]=(long)p1;
      for(i=1;i<=n;i++)
	p1[i]=(long)rnfelement_mulmod(nf,multab,zeronf,unnf,(GEN)matId[j],(GEN)alpha[i],gzero);
    }
    p1=basistoalg(nf,alpha);alphainv=ginv(p1);
    matC=cgetg(n+1,19);
    for(j=1;j<=n;j++)
    {
      p1=cgetg(n*n+1,18);matC[j]=(long)p1;
      for(i=1;i<=n;i++)
      {
	p2=gmul(alphainv,basistoalg(nf,gcoeff(matprod,i,j)));
	for(k=1;k<=n;k++)
	  p1[(i-1)*n+k]=(long)nfreducemodpr(nf,algtobasis(nf,(GEN)p2[k]),prhall);
      }
    }
    matG=nfker(nf,matC,prhall);m=lg(matG)-1;
    matV=cgetg(n+m+1,19);
    for(j=1;j<=m;j++) matV[j]=(long)lift((GEN)matG[j]);
    for(j=1;j<=n;j++)
    {
      p1=cgetg(n+1,18);matV[j+m]=(long)p1;
      for(k=1;k<=n;k++) matV[j+m]=(long)matId[j];
    }
    vecpro=cgetg(3,17);vecpro[1]=(long)matV;
    p1=cgetg(n+m+1,17);vecpro[2]=(long)p1;
    for(i=1;i<=m;i++) p1[i]=(long)idealinv(nf,prh);
    for(i=m+1;i<=n+m;i++)
      p1[i]=(long)idealmul(nf,(GEN)I[i-m],(GEN)alphalistinv[i-m]);
    matH=nfhermite(nf,vecpro);
    p1=algtobasis(nf,gmul(basistoalg(nf,A1),basistoalg(nf,(GEN)matH[1])));
    p2=(GEN)matH[2];
    H=cgetg(n+1,19);
    for(j=1;j<=n;j++)
    {
      p3=cgetg(n+1,18);H[j]=(long)p3;
      for(i=1;i<=n;i++) p3[i]=(long)element_mul(nf,gcoeff(p1,i,j),(GEN)alphalistinv[j]);
    }
    Hid=cgetg(n+1,17);
    for(j=1;j<=n;j++) Hid[j]=(long)idealmul(nf,(GEN)p2[j],(GEN)alphalist[j]);
    if(DEBUGLEVEL>1)
    {
      fprintferr(" Nouvel ordre :\n");outerr((GEN)H);
      outerr((GEN)Hid);flusherr();
    }
    fl=i=1;while(fl&&(i<=n)){if(!vecegal((GEN)I[i],(GEN)Hid[i])) fl=0;i++;}
    if(fl) notfinished=0;
    A=H;I=Hid;
  }
  neworder=cgetg(3,17);neworder[1]=(long)A;neworder[2]=(long)I;
  tetpil=avma;return gerepile(av,tetpil,gcopy(neworder));
}
      
GEN
rnfjoinmodules(GEN nf, GEN x, GEN y)
{
/* given MODULES x and y by their pseudo-bases in HNF, gives a
   pseudo-basis of the module generated by x and y. A usage interne,
   pas de verifications, mais gestion de pile. */

  long av=avma,tetpil,j,lx,ly;
  GEN p1,p2,z,Hx,Hy,Ix,Iy;

  Hx=(GEN)x[1];Ix=(GEN)x[2];Hy=(GEN)y[1];Iy=(GEN)y[2];
  lx=lg(Hx);ly=lg(Hy);
  z=cgetg(3,17);
  p1=cgetg(lx+ly-1,19);z[1]=(long)p1;
  p2=cgetg(lx+ly-1,17);z[2]=(long)p2;
  for(j=1;j<lx;j++) {p1[j]=Hx[j];p2[j]=Ix[j];}
  for(;j<lx+ly-1;j++) {p1[j]=Hy[j-lx+1];p2[j]=Iy[j-lx+1];}
  tetpil=avma;return gerepile(av,tetpil,nfhermite(nf,z));
}

GEN
rnfsimplifybasis(GEN bnf, GEN order)
/* given bnf as output by buchinit and a pseudo-basis of an order
   in HNF [A,I] (or [A,I,D,d] it does not matter), tries to simplify the
   HNF as much as possible. The resulting matrix will be upper triangular
   but the diagonal coefficients will not be equal to 1. The ideals
   are guaranteed to be integral and primitive. */
{
  long av=avma,tetpil,j,N,n;
  GEN p1,id,Az,Iz,nf,A,I;

  checkbnf(bnf);
  if((typ(order)!=17)||(lg(order)<3))
    err(talker,"not a pseudo-basis in nfsimplifybasis");
  A=(GEN)order[1];I=(GEN)order[2];n=lg(A)-1;nf=(GEN)bnf[7];
  N=lgef((GEN)nf[1])-3;id=idmat(N);Iz=cgetg(n+1,17);Az=cgetg(n+1,19);
  for(j=1;j<=n;j++)
  {
    if(vecegal((GEN)I[j],id)) {Iz[j]=(long)id;Az[j]=A[j];}
    else
    {
      p1=content((GEN)I[j]);
      if(!gcmp1(p1))
      {
	Iz[j]=(long)gdiv((GEN)I[j],p1);Az[j]=lmul((GEN)A[j],p1);
      }
      else Az[j]=A[j];
      if(!vecegal((GEN)Iz[j],id))
      {
	p1=isprincipalgen(bnf,(GEN)Iz[j]);
	if(gcmp0((GEN)p1[1]))
	{
	  p1=(GEN)p1[2];Iz[j]=(long)id;
	  Az[j]=(long)element_mulvec(nf,p1,(GEN)Az[j]);
	}
      }
    }
  }
  tetpil=avma;p1=cgetg(lg(order),17);p1[1]=lcopy(Az);p1[2]=lcopy(Iz);
  for(j=3;j<lg(order);j++) p1[j]=lcopy((GEN)order[j]);
  return gerepile(av,tetpil,p1);
}

GEN
rnfsteinitz(GEN nf, GEN order)
/* given a pseudo-basis of an order in HNF [A,I] (or [A,I,D,d] it does
   not matter), gives an nxn matrix (not in HNF) of a pseudo-basis and
   an ideal vector [id,id,...,id,I] such that
   order=nf[7]^(n-1)xI. Since it uses the approximation theorem,
   can be long. */
{
  long av=avma,tetpil,N,j,n;
  GEN id,A,I,p1,p2,a,b;

  nf=checknf(nf);
  N=lgef((GEN)nf[1])-3;id=idmat(N);
  if(typ(order)==10) order=rnfpseudobasis(nf,order);
  if((typ(order)!=17)||(lg(order)<3))
    err(talker,"not a pseudo-matrix in rnfsteinitz");
  A=gcopy((GEN)order[1]);I=gcopy((GEN)order[2]);n=lg(A)-1;
  for(j=1;j<=n-1;j++)
  {
    a=(GEN)I[j];
    if(!vecegal(a,id))
    {
      b=(GEN)I[j+1];
      if(vecegal(b,id))
      {
	p1=(GEN)A[j];A[j]=A[j+1];A[j+1]=lneg(p1);
	I[j]=(long)b;I[j+1]=(long)a;
      }
      else
      {
	p2=nfidealdet1(nf,a,b);
	p1=gadd(element_mulvec(nf,(GEN)p2[1],(GEN)A[j]),element_mulvec(nf,(GEN)p2[2],(GEN)A[j+1]));
	A[j+1]=(long)gadd(element_mulvec(nf,(GEN)p2[3],(GEN)A[j]),element_mulvec(nf,(GEN)p2[4],(GEN)A[j+1]));
	A[j]=(long)p1;
	I[j]=(long)id;I[j+1]=(long)idealmul(nf,a,b);
	p1=content((GEN)I[j+1]);
	if(!gcmp1(p1))
	{I[j+1]=(long)gdiv((GEN)I[j+1],p1);A[j+1]=lmul(p1,(GEN)A[j+1]);}
      }
    }
  }
  tetpil=avma;p1=cgetg(lg(order),17);
  p1[1]=lcopy(A);p1[2]=lcopy(I);
  for(j=3;j<lg(order);j++) p1[j]=lcopy((GEN)order[j]);
  return gerepile(av,tetpil,p1);
}

GEN
rnfbasis(GEN bnf, GEN order)
/* Given bnf as output by buchinit and either an order as output by
   rnfpseudobasis or a polynomial, and outputs a basis if it is free,
   an n+1-generating set if it is not */
{
  long av=avma,tetpil,j,N,n;
  GEN nf,A,I,classe,p1,p2,id;

  checkbnf(bnf);
  nf=(GEN)bnf[7];N=lgef((GEN)nf[1])-3;id=idmat(N);
  if(typ(order)==10) order=rnfpseudobasis(nf,order);
  if((typ(order)!=17)||(lg(order)<3))
    err(talker,"not a pseudo-matrix in rnfbasis");
  A=(GEN)order[1];I=(GEN)order[2];n=lg(A)-1;
  for(j=1;(j<=(n-1))&&vecegal((GEN)I[j],id);j++);
  if(j<n) order=rnfsteinitz(nf,order);
  A=(GEN)order[1];I=(GEN)order[2];classe=(GEN)I[n];
  p1=isprincipalgen(bnf,classe);
  if(gcmp0((GEN)p1[1]))
  {
    tetpil=avma;p2=cgetg(n+1,19);
    for(j=1;j<=n-1;j++) p2[j]=lcopy((GEN)A[j]);
    p2[n]=(long)element_mulvec(nf,(GEN)p1[2],(GEN)A[n]);
  }
  else
  {
    p1=ideal_two_elt(nf,classe);
    tetpil=avma;p2=cgetg(n+2,19);
    for(j=1;j<=n-1;j++) p2[j]=lcopy((GEN)A[j]);
    p2[n]=lmul((GEN)p1[1],(GEN)A[n]);
    p2[n+1]=(long)element_mulvec(nf,(GEN)p1[2],(GEN)A[n]);
  }
  return gerepile(av,tetpil,p2);
}

GEN
rnfhermitebasis(GEN bnf, GEN order)
/* Given bnf as output by buchinit and either an order as output by
   rnfpseudobasis or a polynomial, and outputs a basis (not pseudo)
   in Hermite Normal Form if it exists, zero if not */
{
  long av=avma,tetpil,j,N,n;
  GEN nf,A,I,p1,id;

  checkbnf(bnf);
  nf=(GEN)bnf[7];N=lgef((GEN)nf[1])-3;id=idmat(N);
  if(typ(order)==10)
  {
    order=rnfpseudobasis(nf,order);
    A=(GEN)order[1];
  }
  else
  {
    if((typ(order)!=17)||(lg(order)<3))
      err(talker,"not a pseudo-matrix in rnfbasis");
    A=gcopy((GEN)order[1]);
  }
  I=(GEN)order[2];n=lg(A)-1;
  for(j=1;j<=n;j++)
  {
    if(!vecegal((GEN)I[j],id))
    {
      p1=isprincipalgen(bnf,(GEN)I[j]);
      if(gcmp0((GEN)p1[1]))
	A[j]=(long)element_mulvec(nf,(GEN)p1[2],(GEN)A[j]);
      else {avma=av;return gzero;}
    }
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(A));
}

long
rnfisfree(GEN bnf, GEN order)
{
  long av=avma,n,N,j;
  GEN nf,p1,id,I;
  
  checkbnf(bnf);
  if(gcmp1((GEN)((GEN)((GEN)bnf[8])[1])[1])) return 1;
  nf=(GEN)bnf[7];N=lgef((GEN)nf[1])-3;id=idmat(N);
  if(typ(order)==10) order=rnfpseudobasis(nf,order);
  if((typ(order)!=17)||(lg(order)<3))
    err(talker,"not a pseudo-matrix in rnfisfree");
  I=(GEN)order[2];n=lg(I)-1;
  for(j=1;(j<=n)&&vecegal((GEN)I[j],id);j++);
  if(j>n) {avma=av;return 1;}
  p1=(GEN)I[j];
  for(j++;j<=n;j++) if(!vecegal((GEN)I[j],id)) p1=idealmul(nf,p1,(GEN)I[j]);
  j=gcmp0(isprincipal(bnf,p1));avma=av;
  return j;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~									~*/
/*~			      (Z_K/I)^*					~*/
/*~									~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#define element_mulmodideal(nf,x,y,ideal) (nfreducemodideal(element_mul(nf,x,y),ideal))
#define element_sqrmodideal(nf,x,ideal) (nfreducemodideal(element_sqr(nf,x),ideal))
#define element_divmodideal(nf,x,y,ideal) (nfreducemodideal(element_mul(nf,x,element_div(nf,(GEN)idealaddone(nf,y,ideal)[1],y)),ideal))

GEN
zsigne(GEN nf,GEN alpha,GEN arch,long nba)
{
/* retourne le vecteur des signatures de alpha, a coefficients modulo 2 */
  long i,i1,va;
  GEN vecsign,zeromod2,unmod2,rac;

  if(typ(alpha)==18) alpha=gmul((GEN)nf[7],alpha);
  if(typ(alpha)==9) alpha=(GEN)alpha[2];
  vecsign=cgetg(nba+1,18);
  zeromod2=gmodulcp(gzero,gdeux);unmod2=gmodulcp(gun,gdeux);
  rac=(GEN)nf[6];va=varn((GEN)nf[1]);
  for(i1=0,i=1;i<lg(arch);i++)
    if(signe((GEN)arch[i]))
      vecsign[++i1]=(signe(gsubst(alpha,va,(GEN)rac[i]))>0)?(long)zeromod2:(long)unmod2;
  return vecsign;
}

GEN
nfreducemodideal(GEN x,GEN ideal)
/* a usage interne...Reduction de x modulo l'ideal */
{
  long N,i;
  GEN y,p1,q;
    
  N=lg(x)-1;y=gcopy(x);
  for(i=N;i>=1;i--)
  {
    p1=gcoeff(ideal,i,i);q=gdivround((GEN)y[i],p1);
    if(signe(q)) y=gsub(y,gmul(q,(GEN)ideal[i]));
  }
  return y;
}

GEN
reducemodmatrix(GEN x, GEN y)
/* a usage interne...Reduction de la matrice x modulo la matrice y */
{
  long i,k;
  GEN z;

  y=hnf(y);
  k=lg(x)-1;z=cgetg(k+1,19);
  for(i=1;i<=k;i++) z[i]=(long)nfreducemodideal((GEN)x[i],y);
  return z;
}

GEN
nfreducemodidele(GEN nf,GEN x,GEN arch,GEN structarch,GEN g)
{
/* un element g etant donne, calcule un element congru a g modulo l'ideal
   x et de meme signature aux places de arch */
  
  long i,nba;
  GEN p1,p2,generator;
  
  p1=nfreducemodideal(g,x);
  nba=lg((GEN)structarch[1])-1;generator=(GEN)structarch[2];
  p2=lift(gadd(zsigne(nf,p1,arch,nba),zsigne(nf,g,arch,nba)));
  for(i=1;i<=nba;i++)
    if(signe((GEN)p2[i])) p1=element_mul(nf,p1,(GEN)generator[i]);
  if(gcmp(gnorml2(p1),gnorml2(g))>0) return g;else return p1;
}

GEN
element_powmodideal(GEN nf,GEN x,GEN k,GEN ideal)
{
  long i,f,N;
  GEN k1,y,z;

  N=lgef((GEN)nf[1])-3;k1=k;z=x;f=1;y=cgetg(N+1,18);
  y[1]=un;for(i=2;i<=N;i++) y[i]=zero;
  while(f)
  {
    if(mpodd(k1)) y=element_mulmodideal(nf,z,y,ideal);
    k1=shifti(k1,-1);f=signe(k1);
    if(f) z=element_sqrmodideal(nf,z,ideal);
  }
  return y;
}

GEN
zidealij(GEN x, GEN y)
/* etant donnes deux ideaux entiers x et y en HNF tels que x|y|x^2, trouve le
   quotient (1+x)/(1+y) [h,[clh],[gen],ux^-1]. A usage interne, pas de
   verifs. */
{
  long av=avma,tetpil,j,N,c;
  GEN p1,p2,p3,p4,p5,u,v,d,clh,z;

  p1=invmulmat(x,y);p2=smith2(p1);u=(GEN)p2[1];v=(GEN)p2[2];
  p3=gmul(x,reducemodmatrix(ginv(u),p1));N=lg(p3)-1;
  p4=gtrans(gmul(u,ginv(x)));
  for(j=1;j<=N;j++) coeff(p3,1,j)=(long)addsi(1,gcoeff(p3,1,j));
  d=gmul(u,gmul(p1,v));c=0;clh=gun;
  for(j=1;j<=N;j++) if(!gcmp1(gcoeff(d,j,j))) c++;
  p5=cgetg(c+1,19);for(j=1;j<=c;j++) p5[j]=p4[j];
  tetpil=avma;z=cgetg(4,17);p1=cgetg(c+1,17);p2=cgetg(c+1,17);z[1]=(long)p1;z[2]=(long)p2;
  for(j=1;j<=c;j++) {p1[j]=lcopy(gcoeff(d,j,j));p2[j]=lcopy((GEN)p3[j]);}
  z[3]=(long)gtrans(p5);return gerepile(av,tetpil,z);
}

GEN
zarchstar(GEN nf,GEN x,GEN arch,long nba)
{
/* un ideal x etant donne, calcule des elements de 1+x ayant comme signature
   0 sauf en 1 unique place ou arch=1 */
  
  long av=avma,av1,tetpil,N,i,j,fl,r,rr,limr,k,kk,rankinit,ranknew;
  GEN p1,p2,p3,y,bas,v2,v3,vecsign,genarch,alpha,lambda;

  if(!nba)
  {
    y=cgetg(3,17);p1=cgetg(1,17);y[1]=(long)p1;
    p2=cgetg(1,17);y[2]=(long)p2;return y;
  }
  x=gmul(x,lllintpartial(x));
  N=lgef((GEN)nf[1])-3;
  v2=cgetg(1,19);
  r=1;rr=3;bas=gmul((GEN)nf[7],x);rankinit=0;j=0;fl=1;
  lambda=cgeti(N+1);genarch=cgetg(nba+1,17);
  /* dans cette boucle on va chercher des elements pour que la matrice
     des signatures modulo 2 soit inversible */
  while(fl)
  {
    limr=(itos(gmin(stoi(BIGINT),gpuigs(stoi(rr),N)))-1)>>1;
    for(k=rr;(k<=limr)&&fl;k++)
    {
      kk=k;for(i=1;i<=N;i++) {lambda[i]=(kk+r)%rr-r;kk/=rr;}
      av1=avma;
      alpha=gun;for(i=1;i<=N;i++) alpha=gadd(alpha,gmulsg(lambda[i],(GEN)bas[i]));
      vecsign=zsigne(nf,alpha,arch,nba);
      if(lg(v2)>1) v3=concat(v2,gtomat(vecsign));
      else v3=gtomat(vecsign);ranknew=rank(v3);
      if(ranknew>rankinit)
      {
	v2=v3;rankinit++;fl=(rankinit<nba);
	genarch[++j]=(long)algtobasis(nf,alpha);
      }
      else avma=av1;
    }
    if(fl) {r++;rr+=2;}
  }
  if(rankinit<nba) err(talker,"bug in zidealstarinit");
  v3=lift(ginv(v2));p2=cgetg(nba+1,17);
  p3=cgetg(N+1,18);p3[1]=un;for(i=2;i<=N;i++) p3[i]=zero;
  for(j=1;j<=nba;j++)
  {
    p1=gcopy(p3);
    for(i=1;i<=nba;i++)
      if(signe(gcoeff(v3,i,j))) p1=element_mul(nf,p1,(GEN)genarch[i]);
    p2[j]=(long)p1;
  }
  tetpil=avma;y=cgetg(3,17);p1=cgetg(nba+1,17);y[1]=(long)p1;
  for(i=1;i<=nba;i++) p1[i]=(long)gdeux;y[2]=lcopy(p2);
  return gerepile(av,tetpil,y);
}

GEN
zconvert(GEN nf,GEN uv,GEN x,GEN arch,GEN structarch,GEN g)
{
/* un element g generateur d'un p^k divisant x etant donne, calcule
   un element congru a g modulo p^k et a 1 modulo x/p^k et de plus
   positif aux places de arch */
  
  long i,nba;
  GEN p1,p2,generator;
  
  p1=nfreducemodideal(gadd((GEN)uv[1],element_mul(nf,g,(GEN)uv[2])),x);
  nba=lg((GEN)structarch[1])-1;generator=(GEN)structarch[2];
  p2=lift(zsigne(nf,p1,arch,nba));
  for(i=1;i<=nba;i++)
    if(signe((GEN)p2[i])) p1=element_mul(nf,p1,(GEN)generator[i]);
  return p1;
}

GEN
zprimestar(GEN nf,GEN pr,GEN ep,GEN x,GEN arch,GEN structarch)
{
/* Calcule les donnees necessaires correspondant a pr^ep divisant x */
  
  long av=avma,tetpil,N,f,fl2,j,m,n,i,fl,psim,k,e,nbp,nba,a,b;
  GEN prh,p,pef,pefm1,fa,list,v,prhall,p1,p2,p3,p4,prk,pefm2,pra,prb;
  GEN list_of_struct,newgen,g0,uv,q;

  prh=idealhermite(nf,pr);N=lg(prh)-1;
  f=itos((GEN)pr[4]);k=itos(ep);p=(GEN)pr[1];pef=gpuigs(p,f);
  pefm1=gaddgs(pef,-1);fa=factor(pefm1);list=(GEN)fa[1];nbp=lg(list)-1;
  v=cgetg(N+1,18);psim=itos(p);nba=lg((GEN)structarch[1])-1;
  if(f==1) {v[1]=gener(p)[2];for(j=2;j<=N;j++) v[j]=zero;}
  else
  {
    prhall=cgetg(3,17);prhall[1]=(long)prh;prhall[2]=zero;
    fl=1;
    for(n=psim;fl;n++)
    {
      m=n;
      for(i=1;i<=N;i++)
      {if(!gcmp1(gcoeff(prh,i,i))) {v[i]=lstoi(m%psim);m/=psim;}else v[i]=zero;}
      fl2=1;
      for(j=1;(j<=nbp)&&fl2;j++)
      {
	p1=lift(element_powmodpr(nf,v,divii(pefm1,(GEN)list[j]),prhall));
	p1[1]=(long)addsi(-1,(GEN)p1[1]);
	if(gcmp0(p1)||idealval(nf,p1,pr)) fl2=0;
      }
      if(fl2) fl=0;
    }
  }
  prk=idealpow(nf,pr,ep);pefm2=gaddgs(pef,-2);e=1;
  while(e<k)
  {
    p1=element_powmodideal(nf,v,pefm2,prk);p2=element_mulmodideal(nf,v,p1,prk);
    p2[1]=(long)addsi(-1,(GEN)p2[1]);
    v=gsub(v,element_divmodideal(nf,p2,gmul(pefm1,p1),prk));e<<=1;
  }
  q=idealpow(nf,pr,ep);
  uv=idealaddone(nf,q,idealdivexact(nf,x,q));
/*  g0=zconvert(nf,uv,x,arch,structarch,v); */
  g0=nfreducemodideal(gadd((GEN)uv[1],element_mul(nf,v,(GEN)uv[2])),x);  
  list_of_struct=cgetg(2,17);
  p1=cgetg(5,17);list_of_struct[1]=(long)p1;
  p2=cgetg(2,17);p1[1]=(long)p2;p2[1]=(long)pefm1;
  p2=cgetg(2,17);p1[2]=(long)p2;p2[1]=(long)g0;
  p1[3]=un;p2=cgetg(2,17);p1[4]=(long)p2;p2[1]=(long)zsigne(nf,g0,arch,nba);
  if(gcmp1(ep)){tetpil=avma;return gerepile(av,tetpil,gcopy(list_of_struct));}
  else
  {
    e=itos(ep);a=1;b=2;pra=idealhermite(nf,pr);prb=idealmul(nf,pra,pr);
    while(a<b)
    {
      p1=zidealij(pra,prb);newgen=(GEN)p1[2];p3=cgetg(lg(newgen),17);
      for(i=1;i<lg(newgen);i++)
/*	newgen[i]=(long)zconvert(nf,uv,x,arch,structarch,(GEN)newgen[i]); */
      {
	newgen[i]=(long)nfreducemodideal(gadd((GEN)uv[1],element_mul(nf,(GEN)newgen[i],(GEN)uv[2])),x);
	p3[i]=(long)zsigne(nf,(GEN)newgen[i],arch,nba);
      }
      p4=cgetg(5,17);for(i=1;i<=3;i++) p4[i]=p1[i];p4[4]=(long)p3;
      p2=cgetg(2,17);p2[1]=(long)p4;list_of_struct=concat(list_of_struct,p2);
      a=b;b=min(e,b<<1);
      if(a!=b)
      {
	pra=prb;
	if(b==(a>>1)) prb=idealmul(nf,pra,pra);
	else prb=idealmul(nf,pra,idealpow(nf,pr,stoi(e-a)));
      }
    }
    tetpil=avma;return gerepile(av,tetpil,gcopy(list_of_struct));
  }
}

GEN
compute_prhall(GEN nf,GEN pr)
{
/* calcule le prhall associe a pr necessaire pour faire des reductions modulo
   pr meme pour des non entiers algebriques */
  
  long N;
  GEN prh,epr,betae,p1,prhall,id;

  N=lgef((GEN)nf[1])-3;id=idmat(N);
  prh=idealmul(nf,id,pr);epr=(GEN)pr[3];
  betae=gdiv(element_pow(nf,(GEN)pr[5],epr),gpui((GEN)pr[1],addis(epr,-1),0));
  p1=cgetg(2,19);p1[1]=(long)betae;
  p1=idealadd(nf,gmul((GEN)pr[1],id),idealmul(nf,p1,id));
  prhall=cgetg(3,17);prhall[1]=(long)prh;
  prhall[2]=idealaddone(nf,pr,p1)[2];
  return prhall;
}

long
nfsearch(GEN smalltable,GEN x)
/* rend l'indice de x dans la table smalltable, sinon 0 */
{
  long l,u,i,s;

  u=lg(smalltable)-1;l=1;
  while(u>=l)
  {
    i=(l+u)>>1;s=lexcmp(x,(GEN)smalltable[i]);
    if(!s) return i;
    if(s<0) u=i-1;else l=i+1;
  }
  return 0;
}

GEN
nfshanks(GEN nf,GEN x,GEN g0,GEN pr,GEN prhall)
/* rend le plus petit entier positif n tel que g0^n=x modulo pr */
{
  long av=avma,tetpil,i,N,lbaby,fl,k;
  GEN pf1,pfqr,p1,smalltable,p2,giant,smalltable2,permtable,v,g0inv;

  pf1=addsi(-1,gpui((GEN)pr[1],(GEN)pr[4],0));pfqr=racine(pf1);
  if(cmpis(pfqr,65535)>=0) err(talker,"module is too large in nfshanks");
  lbaby=itos(pfqr);N=lgef((GEN)nf[1])-3;
  g0inv=element_invmodpr(nf,g0,prhall);
  p1=x;smalltable=cgetg(lbaby+2,17);fl=0;
  for(i=0;(i<=lbaby)&&(!fl);i++)
  {
    p2=lift(p1);smalltable[i+1]=(long)p2;p2[1]=(long)addsi(-1,(GEN)p2[1]);fl=gcmp0(p2);
    if(!fl) fl=element_val(nf,p2,pr);p2[1]=(long)addsi(1,(GEN)p2[1]);
    if((i<lbaby)&&(!fl)) p1=element_mulmodpr(nf,p1,g0inv,prhall);
  }
  if(fl) {avma=av;return stoi(i-1);}
  permtable=indexlexsort(smalltable);
  smalltable2=cgetg(lbaby+2,17);
  for(i=1;i<=lbaby+1;i++) smalltable2[i]=smalltable[itos((GEN)permtable[i])];
  smalltable=smalltable2;
  giant=element_divmodpr(nf,x,p1,prhall);p1=giant;k=1;
  while(1)
  {
    if((i=nfsearch(smalltable,lift(p1))))
    {
      v=addii(mulsi(lbaby,stoi(k)),(GEN)permtable[i]);tetpil=avma;
      return gerepile(av,tetpil,addsi(-1,v));
    }
    p1=element_mulmodpr(nf,p1,giant,prhall);k++;
  }
}

GEN
zinternallog(GEN nf,GEN set_of_list,long nbgen,GEN x,GEN arch,GEN fa,GEN a)
/* A usage interne : retourne la decomposition de a sur les nbgen generateurs
   successifs contenus dans set_of_list */
{
  long av=avma,tetpil,nbp,cp,k,e,j,i,nba;
  GEN list,ep,y,ainit,list_of_struct,primetop,g0,pr,p,prhall,structarch,pk1,pk;
  GEN prk,p1,p2,b,p3,p4,p6,psigne;

  list=(GEN)fa[1];ep=(GEN)fa[2];nbp=lg(list)-1;cp=0;ainit=gcopy(a);
  structarch=(GEN)set_of_list[nbp+1];nba=lg((GEN)structarch[1])-1;
  psigne=zsigne(nf,ainit,arch,nba);
  y=cgetg(nbgen-nba+1,18);
  for(k=1;k<=nbp;k++)
  {
    a=ainit;
    list_of_struct=(GEN)set_of_list[k];primetop=(GEN)list_of_struct[1];
    g0=(GEN)((GEN)primetop[2])[1];
    pr=(GEN)list[k];p=(GEN)pr[1];e=itos((GEN)ep[k]);prhall=compute_prhall(nf,pr);
    y[++cp]=(long)nfshanks(nf,a,g0,pr,prhall);
    if(mpodd((GEN)y[cp])) psigne=gadd(psigne,(GEN)((GEN)primetop[4])[1]);
    pk1=gpuigs(p,e-1);pk=gmul(pk1,p);prk=idealpow(nf,pr,(GEN)ep[k]);
    p2=element_powmodideal(nf,g0,subii((GEN)primetop[1],(GEN)y[cp]),prk);
    a=element_mulmodideal(nf,a,p2,prk);
    for(j=2;j<lg(list_of_struct);j++)
    {
      p1=(GEN)list_of_struct[j];
      b=gcopy(a);b[1]=(long)addsi(-1,(GEN)b[1]);
      p2=gmul((GEN)p1[3],b);p3=(GEN)p1[1];p4=(GEN)p1[2];
      if(lg(p2)!=lg(p3)) err(talker,"bug in zinternallog");
      for(i=1;i<lg(p3);i++)
      {
	y[++cp]=(long)negi(p6=modii(negi((GEN)p2[i]),(GEN)p3[i]));
	if(mpodd((GEN)y[cp])) psigne=gadd(psigne,(GEN)((GEN)p1[4])[i]);
	a=element_mul(nf,element_pow(nf,(GEN)p4[i],p6),a);
      }
    }
  }
  p1=lift(psigne);tetpil=avma;return gerepile(av,tetpil,concat(y,p1));
}

GEN
zidealstarinit(GEN nf, GEN ideal)
/* Calcule [[ideal,arch],[h,[cyc],[gen]],idealfact,[liste],U] */
{
  long av=avma,tetpil,i,j,k,nba,nbp,N,fl,c,s,R1,nbgen,cp,jj;
  GEN p1,p2,p3,p3plus,p3moins,p4,y,h,clh,met,u1,basecl,generator,mot;
  GEN fa,fa2,list,ep,x,arch,allgenerator,list_of_struct,structarch,u1u2,u2;
  
  nf=checknf(nf);N=lgef((GEN)nf[1])-3;R1=itos((GEN)((GEN)nf[2])[1]);
  if((typ(ideal)==17)&&(lg(ideal)==3))
  {
    fl=1;x=(GEN)ideal[1];arch=(GEN)ideal[2];nba=0;
    if((typ(arch)!=17)&&(typ(arch)!=18)&&(lg(arch)!=1+R1))
      err(talker,"incorrect archimedean component in zidealstarinit");
    for(i=1;i<=R1;i++) if(signe((GEN)arch[i])) nba++;
  }
  else
  {
    fl=0;x=ideal;arch=cgetg(R1+1,17);for(i=1;i<=R1;i++) arch[i]=zero;nba=0;
    ideal=cgetg(3,17);ideal[1]=(long)x;ideal[2]=(long)arch;
  }
  x=idealhermite(nf,x);
  if(!gcmp1(denom(x))) err(talker,"zidealstarinit needs an integral ideal");
  y=cgetg(6,17);y[1]=(long)ideal;
  fa=idealfactor(nf,x);list=(GEN)fa[1];ep=(GEN)fa[2];nbp=lg(list)-1;
  y[3]=(long)fa;fa2=cgetg(nbp+2,17);y[4]=(long)fa2;
  structarch=zarchstar(nf,x,arch,nba);allgenerator=cgetg(1,17);
  for(i=1;i<=nbp;i++)
  {
    p1=zprimestar(nf,(GEN)list[i],(GEN)ep[i],x,arch,structarch);
    fa2[i]=(long)p1;
    for(j=1;j<lg(p1);j++)
      allgenerator=concat(allgenerator,(GEN)((GEN)p1[j])[2]);
  }
  fa2[nbp+1]=(long)structarch;
  allgenerator=concat(allgenerator,(GEN)structarch[2]);
  nbgen=lg(allgenerator)-1;
  h=cgetg(nbgen+1,19);cp=0;
  for(k=1;k<lg(fa2)-1;k++)
  {
    list_of_struct=(GEN)fa2[k];
    for(j=1;j<lg(list_of_struct);j++)
    {
      p1=(GEN)list_of_struct[j];p2=(GEN)p1[1];p3=(GEN)p1[2];
      for(jj=1;jj<lg((GEN)p1[1]);jj++)
      {
	h[++cp]=(long)gneg(zinternallog(nf,fa2,nbgen,x,arch,fa,element_pow(nf,(GEN)p3[jj],(GEN)p2[jj])));
	coeff(h,cp,cp)=((GEN)p1[1])[jj];
      }
    }
  }
  p1=cgetg(nbgen+1,18);for(i=1;i<=nbgen;i++) p1[i]=zero;
  for(j=1;j<=nba;j++) {h[++cp]=lcopy(p1);coeff(h,cp,cp)=deux;}
  if(cp!=nbgen) err(talker,"bug in zidealstarinit");
  if(nbgen)
  {
    u1u2=smith2(h);u1=(GEN)u1u2[1];u2=(GEN)u1u2[2];
    met=gmul(u1,gmul(h,u2));y[5]=(long)u1;u1=reducemodmatrix(ginv(u1),h);
    clh=gun;c=0;
    for(i=1;i<=nbgen;i++)
      if(!gcmp1(gcoeff(met,i,i))) {clh=mulii(clh,gcoeff(met,i,i));c++;}
  }
  else {clh=gun;met=cgetg(1,19);u1=cgetg(1,19);y[5]=(long)u1;c=0;}
  basecl=cgetg(c+1,17);
  p4=cgetg(N+1,18);p4[1]=un;for(i=2;i<=N;i++) p4[i]=zero;
  for(j=1;j<=c;j++)
  {
    p3plus=p4;p3moins=p4;
    for(i=1;i<=nbgen;i++)
    {
      p1=gcoeff(u1,i,j);s=signe(p1);
      if(s)
      {
	if(s>0) p3plus=element_mul(nf,p3plus,(GEN)element_pow(nf,(GEN)allgenerator[i],p1));
	else p3moins=element_mul(nf,p3moins,(GEN)element_pow(nf,(GEN)allgenerator[i],negi(p1)));
      }
    }
    p4=element_div(nf,(GEN)idealaddone(nf,idealhermite(nf,p3moins),x)[1],p3moins);
/* p4 est l'inverse de p3moins modulo x */
    p3=element_mulmodideal(nf,p3plus,p4,x);
    generator=(GEN)structarch[2];
    p2=lift(gadd(gadd(zsigne(nf,p3,arch,nba),zsigne(nf,p3plus,arch,nba)),zsigne(nf,p3moins,arch,nba)));
    for(i=1;i<=nba;i++)
      if(signe((GEN)p2[i])) p3=element_mul(nf,p3,(GEN)generator[i]);
/* on a corrige pour que p3 ait la meme signature que p3plus/p3moins */    
    basecl[j]=(long)p3;
  }
  p1=cgetg(4,17);p1[1]=(long)clh;
  mot=cgetg(c+1,17);for(i=1;i<=c;i++) mot[i]=coeff(met,i,i);
  p1[2]=(long)mot;p1[3]=(long)basecl;
  y[2]=(long)p1;
  tetpil=avma;return gerepile(av,tetpil,gcopy(y));
}
  
GEN
smithspec(GEN nf, GEN fa, GEN fa2, GEN x)
{
  long av=avma,tetpil,sizeh,nbp,i,j,ind,c,jj;
  GEN list,ep,h,generators,uv,p1,p2,p3,p4,q,moduli,a,d,z,u,v,U,hid,hjd;
    
  list=(GEN)fa[1];ep=(GEN)fa[2];nbp=lg(list)-1;
  sizeh=0;for(i=1;i<=nbp;i++) sizeh+=(lg((GEN)((GEN)fa2[i])[2])-1);
  h=cgetg(sizeh+1,17);U=idmat(sizeh);
  generators=cgetg(sizeh+1,17);moduli=cgetg(sizeh+1,17);
  ind=0;
  for(i=1;i<=nbp;i++)
  {
    p1=(GEN)fa2[i];p2=(GEN)p1[2];p3=(GEN)p1[3];
    p4=idealpow(nf,(GEN)list[i],(GEN)ep[i]);
    for(j=1;j<lg(p2);j++)
    {
      ind++;h[ind]=p2[j];moduli[ind]=(long)p4;
      generators[ind]=p3[j];
    }
  }
  for(i=1;i<=sizeh;i++)
  {
    q=(GEN)moduli[i];a=(GEN)generators[i];
    uv=idealaddone(nf,q,idealdivexact(nf,x,q));
    generators[i]=(long)nfreducemodideal(gadd((GEN)uv[1],element_mul(nf,a,(GEN)uv[2])),x);
  }
  for(i=sizeh;i>=2;i--)
  {
    for(j=i-1;j>=1;j--)
    {
      if(!divise((GEN)h[j],(GEN)h[i]))
      {
	d=bezout((GEN)h[i],(GEN)h[j],&u,&v);
	hid=divii((GEN)h[i],d);hjd=gneg(divii((GEN)h[j],d));
	if(signe(v)<0) {v=addii(v,hid);u=addii(u,hjd);}
	h[j]=(long)mulii((GEN)h[i],q=divii((GEN)h[j],d));h[i]=(long)d;
	generators[j]=(long)nfreducemodideal(element_mul(nf,(GEN)generators[j],element_div(nf,(GEN)idealaddone(nf,(GEN)generators[i],x)[1],(GEN)generators[i])),x);
	for(jj=1;jj<=sizeh;jj++)
	  coeff(U,i,jj)=ladd(gcoeff(U,i,jj),gcoeff(U,j,jj));
	generators[i]=(long)element_mulmodideal(nf,(GEN)generators[i],element_powmodideal(nf,(GEN)generators[j],mulii(v,q),x),x);
	p2=gmul(hjd,v);
	for(jj=1;jj<=sizeh;jj++)
	  coeff(U,j,jj)=ladd(gcoeff(U,j,jj),gmul(p2,gcoeff(U,i,jj)));
      }
    }
  }
  q=gun;for(i=1;(i<=sizeh)&&(!gcmp1((GEN)h[i]));i++) q=mulii(q,(GEN)h[i]);
  c=i-1;tetpil=avma;
  z=cgetg(5,17);z[1]=lcopy(q);
  p1=cgetg(c+1,17);z[2]=(long)p1;
  for(i=1;i<=c;i++) p1[i]=lcopy((GEN)h[i]);
  p1=cgetg(c+1,17);z[3]=(long)p1;
  for(i=1;i<=c;i++) p1[i]=lcopy((GEN)generators[i]);
  z[4]=lcopy(U);
  return gerepile(av,tetpil,z);
}

GEN
zunramifiedprimestar(GEN nf, GEN pr, GEN ep)
{
  long av=avma,tetpil,N,i,j,k,fl,fl2,f,psim,n,m,nbp,e;
  GEN y,p1,p2,prk,prh,h,p,pef,pefm1,pefm2,pkm1,generators,list,v,fa,prhall;

  prh=idealhermite(nf,pr);N=lg(prh)-1;
  f=itos((GEN)pr[4]);k=itos(ep);p=(GEN)pr[1];pef=gpuigs(p,f);
  pefm1=gaddgs(pef,-1);fa=factor(pefm1);list=(GEN)fa[1];nbp=lg(list)-1;
  v=cgetg(N+1,18);psim=itos(p);
  if(f==1) {v[1]=gener(p)[2];for(j=2;j<=N;j++) v[j]=zero;}
  else
  {
    prhall=cgetg(3,17);prhall[1]=(long)prh;prhall[2]=zero;
    fl=1;
    for(n=psim;fl;n++)
    {
      m=n;
      for(i=1;i<=N;i++)
      {if(!gcmp1(gcoeff(prh,i,i))) {v[i]=lstoi(m%psim);m/=psim;}else v[i]=zero;}
      fl2=1;
      for(j=1;(j<=nbp)&&fl2;j++)
      {
	p1=lift(element_powmodpr(nf,v,divii(pefm1,(GEN)list[j]),prhall));
	p1[1]=(long)addsi(-1,(GEN)p1[1]);
	if(gcmp0(p1)||idealval(nf,p1,pr)) fl2=0;
      }
      if(fl2) fl=0;
    }
  }
  prk=idealpow(nf,pr,ep);pefm2=gaddgs(pef,-2);e=1;
  while(e<k)
  {
    p1=element_powmodideal(nf,v,pefm2,prk);p2=element_mulmodideal(nf,v,p1,prk);
    p2[1]=(long)addsi(-1,(GEN)p2[1]);
    v=gsub(v,element_divmodideal(nf,p2,gmul(pefm1,p1),prk));e<<=1;
  }
  if(k==1)
  {
    h=cgetg(2,17);h[1]=(long)pefm1;generators=cgetg(2,17);
    generators[1]=(long)v;
  }
  else
  {
    if((psim!=2)||((k==2)&&(f>1)))
    {
      pkm1=gpuigs(p,k-1);h=cgetg(f+2,17);
      h[1]=(long)pefm1;for(i=2;i<=f+1;i++) h[i]=(long)pkm1;
      generators=cgetg(f+2,17);generators[1]=(long)v;
      v=cgetg(N+1,18);for(j=1;j<=N;j++) v[j]=zero;
      for(j=1,i=1;i<=N;i++)
      {
	if(!gcmp1(gcoeff(prh,i,i)))
	{
	  v[i]=un;p1=(i==1)?gcopy((GEN)pr[2]):element_mul(nf,(GEN)pr[2],v);
	  v[i]=zero;
	  if(psim!=2)
	  {p1[1]=(long)addsi(1,(GEN)p1[1]);generators[++j]=(long)p1;}
	  else
	  {
	    p1[1]=(long)addsi(2,(GEN)p1[1]);
	    if((!gcmp0(p1))&&(element_val(nf,p1,pr)==1))
	    {
	      p1[1]=(long)addsi(-1,(GEN)p1[1]);
	      generators[++j]=(long)p1;
	    }
	  }
	}
      }
      if(psim==2) {v[1]=lneg(gun);generators[++j]=lcopy(v);}
    }
    else /* ici p=2 */
    {
      if(k==2) /* ici p=2 et f=1 */
      {
	h=cgetg(2,17);h[1]=(long)gdeux;generators=cgetg(2,17);
	p1=gcopy((GEN)pr[2]);p1[1]=(long)addsi(1,(GEN)p1[1]);
	generators[1]=(long)p1;
      }
      else /* ici p=2 et k>=3 */
      {
	pkm1=gpuigs(p,k-1);
	if(f>1)
	{
	  h=cgetg(f+3,17);
	  h[1]=(long)pefm1;for(i=2;i<=f;i++) h[i]=(long)pkm1;
	  h[f+1]=(long)divii(pkm1,p);h[f+2]=(long)gdeux;
	  generators=cgetg(f+3,17);generators[1]=(long)v;
	  v=cgetg(N+1,18);for(j=1;j<=N;j++) v[j]=zero;
	  for(j=1,i=1;i<=N;i++)
	  {
	    if(!gcmp1(gcoeff(prh,i,i)))
	    {
	      v[i]=un;p1=(i==1)?gcopy((GEN)pr[2]):element_mul(nf,(GEN)pr[2],v);
	      v[i]=zero;p1[1]=(long)addsi(2,(GEN)p1[1]);
	      if((!gcmp0(p1))&&(element_val(nf,p1,pr)==1))
	      {
		p1[1]=(long)addsi(-1,(GEN)p1[1]);
		generators[++j]=(long)p1;
	      }
	    }
	  }
	  v[1]=lstoi(5);generators[f+1]=lcopy(v);v[1]=lneg(gun);
	  generators[f+2]=lcopy(v);
	}
	else
	{
	  h=cgetg(3,17);h[1]=(long)divii(pkm1,p);h[2]=(long)gdeux;
	  generators=cgetg(3,17);v=cgetg(N+1,18);for(j=2;j<=N;j++) v[j]=zero;
	  v[1]=lstoi(5);generators[1]=lcopy(v);v[1]=lneg(gun);
	  generators[2]=lcopy(v);
	}
      }
    }
    p1=gpuigs(pef,k-1);
  }
  tetpil=avma;
  y=cgetg(4,17);y[1]=(k==1)?lcopy(pefm1):lmul(p1,pefm1);
  y[2]=lcopy(h);y[3]=lcopy(generators);
  return gerepile(av,tetpil,y);
}

GEN
zprimestarold(GEN nf,GEN pr,GEN ep)
{
  if(!gcmp1((GEN)pr[3]))
  {err(impl,"zprimestar for ramified modules");return gnil;}
  else return zunramifiedprimestar(nf,pr,ep);
}

GEN compute_class_number(GEN mit,GEN *met,GEN *u1);

GEN
zidealstarinitold(GEN nf, GEN ideal)
/* Calcule [[ideal,arch],[h,[cyc],[gen]],idealfact,[liste de [h,[cyc],[gen]] brute],
   [U,v2^-1,U']] si arch est inclus
   [ideal,[h,[cyc],[gen]],idealfact,[liste de [h,[cyc],[gen]] brute],U] si
   arch est omis
   */
{
  long av=avma,tetpil,i,i1,j,k,nba,nbp,N,va,ngen,sizeh,r,rr,rankinit,ranknew;
  long fl,limr,kk,av1,c,s,R1;
  GEN p1,p2,p3,y,rac,zeromod2,unmod2,cyclic,cyclicplus,generator,genplus;
  GEN v2,v3,v,bas,lambda,vecsign,h,clh,met,u1,u1init,basecl,mot,alpha;
  GEN fa,fa2,list,ep,BIGU,x,arch;
  
  nf=checknf(nf);
  if((typ(ideal)==17)&&(lg(ideal)==3))
  {
    fl=1;x=(GEN)ideal[1];arch=(GEN)ideal[2];
    if((typ(arch)!=17)&&(typ(arch)!=18)&&(lg(arch)!=1+itos((GEN)((GEN)nf[2])[1])))
      err(talker,"incorrect archimedean component in zidealstarinit");
  }
  else {fl=0;x=ideal;}
  x=idealhermite(nf,x);
  if(!gcmp1(denom(x))) err(talker,"zidealstarinit needs an integral ideal");
  y=cgetg(6,17);y[1]=(long)ideal;
  fa=idealfactor(nf,x);list=(GEN)fa[1];ep=(GEN)fa[2];nbp=lg(list)-1;
  y[3]=(long)fa;fa2=cgetg(nbp+1,17);y[4]=(long)fa2;
  for(i=1;i<=nbp;i++)  fa2[i]=(long)zprimestarold(nf,(GEN)list[i],(GEN)ep[i]);
  p1=smithspec(nf,fa,fa2,x); /* on s'est occupe des places non-archimediennes */
  if(!fl)
  { /* si fl=0, pas de partie archimedienne donc on s'en va */
    p2=cgetg(4,17);y[2]=(long)p2;
    for(i=1;i<=3;i++) p2[i]=p1[i];
    y[5]=p1[4];tetpil=avma;return gerepile(av,tetpil,gcopy(y));
  }
  else
  { 
    BIGU=cgetg(4,17);y[5]=(long)BIGU;BIGU[1]=p1[4];
    va=varn((GEN)nf[1]);N=lgef((GEN)nf[1])-3;R1=itos((GEN)((GEN)nf[2])[1]);
    rac=(GEN)nf[6];
    nba=0;for(i=1;i<lg(arch);i++) if(signe((GEN)arch[i])) nba++;
    if(!nba)
    { /* en fait pas de partie archimedienne. On s'en va, mais le format de
	 sortie est un peu different */
      BIGU[2]=lgetg(1,19);BIGU[3]=(long)idmat(lg((GEN)p1[2])-1);
      p2=cgetg(4,17);y[2]=(long)p2;
      for(i=1;i<=3;i++) p2[i]=p1[i];
      tetpil=avma;return gerepile(av,tetpil,gcopy(y));
    }
    zeromod2=gmodulcp(gzero,gdeux);unmod2=gmodulcp(gun,gdeux);
    cyclic=(GEN)p1[2];ngen=lg(cyclic)-1;generator=(GEN)p1[3];
	/* cyclic et generators contiennent la structure et les generateurs
	   de ce qu'il faut sans partie archimedienne, de longueur ngen */
    sizeh=ngen+nba;
    genplus=cgetg(sizeh+1,18);for(i=1;i<=ngen;i++) genplus[i]=generator[i];
    cyclicplus=cgetg(sizeh+1,18);for(i=1;i<=ngen;i++) cyclicplus[i]=cyclic[i];
    for(i=ngen+1;i<=sizeh;i++) cyclicplus[i]=deux;
    cyclic=cyclicplus;
    v2=cgetg(1,19);
    r=1;rr=3;bas=gmul((GEN)nf[7],x);rankinit=0;j=ngen;fl=1;
    lambda=cgeti(N+1);vecsign=cgetg(nba+1,18);
	/* dans cette boucle on va chercher des elements pour que la matrice
	   des signatures modulo 2 soit inversible */
    while(fl)
    {
      limr=(itos(gmin(stoi(BIGINT),gpuigs(stoi(rr),N)))-1)>>1;
      for(k=rr;(k<=limr)&&fl;k++)
      {
	kk=k;for(i=1;i<=N;i++) {lambda[i]=(kk+r)%rr-r;kk/=rr;}
	av1=avma;
	alpha=gun;for(i=1;i<=N;i++) alpha=gadd(alpha,gmulsg(lambda[i],(GEN)bas[i]));
	for(i1=0,i=1;i<=R1;i++)
	  if(signe((GEN)arch[i]))
	    vecsign[++i1]=(signe(gsubst(alpha,va,(GEN)rac[i]))>0)?(long)zeromod2:(long)unmod2;
	if(lg(v2)>1) v3=concat(v2,gtomat(vecsign));
	else v3=gtomat(vecsign);ranknew=rank(v3);
	if(ranknew>rankinit)
	{
	  v2=v3;rankinit++;fl=(rankinit<nba);j++;
	  genplus[j]=(long)algtobasis(nf,alpha);
	}
	else avma=av1;
      }
      if(fl) {r++;rr+=2;}
    }
    if(rankinit<nba) err(talker,"bug in zidealstarinit");
    generator=genplus;
    v=ginv(v2);
    h=cgetg(sizeh+1,19);
	/* h va etre la matrice carree complete de taille sizeh=ngen+nba */
    for(j=1;j<=ngen;j++)
    {
      p1=cgetg(sizeh+1,18);h[j]=(long)p1;
      if(!mpodd((GEN)cyclic[j])) for(i=ngen+1;i<=sizeh;i++) p1[i]=zero;
      else
      {
	for(i1=0,i=1;i<=R1;i++)
	  if(signe((GEN)arch[i]))
	    vecsign[++i1]=(signe(gsubst(gmul((GEN)nf[7],(GEN)generator[j]),va,(GEN)rac[i]))>0)?(long)zeromod2:(long)unmod2;
	p3=gmul(v,vecsign);
	for(i=ngen+1;i<=sizeh;i++) p1[i]=(long)lift((GEN)p3[i-ngen]);
      }
      for(i=1;i<=ngen;i++) p1[i]=(i==j)?cyclic[i]:zero;
    }
    for(j=ngen+1;j<=sizeh;j++)
    {
      p1=cgetg(sizeh+1,18);h[j]=(long)p1;
      for(i=1;i<=sizeh;i++) p1[i]=(i==j)?deux:zero;
    }
  }
      /* on fait un smith, comme le programme est deja ecrit ailleurs
	 inutile de le reecrire. u1 est la matrice INVERSE de la matrice
	 de transformation des LIGNES dans smith, et reduite par
	 reducemodmatrix */
  clh=compute_class_number(h,&met,&u1);
  u1init=u1;u1=reducemodmatrix(u1,h);
  c=0;for(i=1;i<=sizeh;i++) if(!gcmp1(gcoeff(met,i,i))) c++;
  basecl=cgetg(c+1,17);
  for(j=1;j<=c;j++)
  {
    p1=gcoeff(u1,1,j);
    if(signe(p1)<0)
    {p2=(GEN)cyclic[1];if(mpodd(p2)) p2=shifti(p2,1);p1=modii(p1,p2);}
	/* Comme on veut des entiers, on ne veut elever qu'a des puissances
	   positives ou nulles. Il est loisible de remplacer a^-k par
	   a^(c-k) si c est pair ou par a^(2c-k) si c est impair, puisque
	   cela ne change pas la signature, ou c est l'ordre de a modulo l'ideal */
    p3=(GEN)element_pow(nf,(GEN)generator[1],p1);
    for(i=2;i<=sizeh;i++)
    {
      p1=gcoeff(u1,i,j);s=signe(p1);
      if(s)
      {
	if(s<0)
	{p2=(GEN)cyclic[i];if(mpodd(p2)) p2=shifti(p2,1);p1=modii(p1,p2);}
	p3=element_mul(nf,p3,(GEN)element_pow(nf,(GEN)generator[i],p1));
      }
    }
    basecl[j]=(long)p3;
  }
  p1=cgetg(4,17);p1[1]=(long)clh;
  mot=cgetg(c+1,17);for(i=1;i<=c;i++) mot[i]=coeff(met,i,i);
  p1[2]=(long)mot;p1[3]=(long)basecl;
  y[2]=(long)p1;BIGU[2]=(long)v;BIGU[3]=(long)ginv(u1init);
      /* il faut recalculer le u1 initial (ou compute_class_number pourrait
	 le rendre aussi */
  tetpil=avma;return gerepile(av,tetpil,gcopy(y));
}
    
GEN
zidealstarold(GEN nf, GEN x)
{
  long av=avma,tetpil;
  GEN y;
  
  y=zidealstarinitold(nf,x);
  tetpil=avma;return gerepile(av,tetpil,gcopy((GEN)y[2]));
}

GEN
zidealstar(GEN nf, GEN x)
{
  long av=avma,tetpil;
  GEN y;
  
  y=zidealstarinit(nf,x);
  tetpil=avma;return gerepile(av,tetpil,gcopy((GEN)y[2]));
}

GEN
nflogpadic(GEN nf,GEN x,GEN pr,GEN prk,long k)
/* a usage interne : l'ideal pr est non ramifie et x-1 est dans pr */
{
  long av=avma,tetpil,i,v,n,lk,a;
  double lp;
  GEN p,y,s,yn,beta,pk,pkp,pa,uu,vv,p1,m;

  p=(GEN)pr[1];
  y=gneg(x);y[1]=laddsi(1,(GEN)y[1]);
  if(gcmp0(y))
  {
    yn=gadd(gun,ggrandocp(p,k));tetpil=avma;
    return gerepile(av,tetpil,gmul(y,yn));
  }
  s=gzero;v=element_val(nf,y,pr);n=max((k+v-1)/v,2);yn=y;
  lp=log(gtodouble(p));lk=(long)ceil(log(2.0*k/v)/lp);
  pkp=idealmul(nf,prk,idealpow(nf,pr,stoi(lk)));pk=gpuigs(p,k+lk);
  beta=(GEN)idealaddone(nf,pkp,idealdivexact(nf,idealhermite(nf,p),idealpow(nf,pr,(GEN)pr[3])))[2];
  for(i=1;(i<=n)||(i*v<k+log((double)i)/lp);i++)
  {
    a=ggval(stoi(i),p);m=divsi(i,pa=gpuigs(p,a));bezout(m,pk,&uu,&vv);
    p1=gdiv(gmul(element_mul(nf,element_pow(nf,beta,stoi(a)),yn),uu),pa);
    s=nfreducemodideal(gadd(s,p1),prk);
    yn=element_mulmodideal(nf,yn,y,pkp);
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(s));
}

GEN
zideallog(GEN nf, GEN x, GEN bigideal)
{
  long av=avma,tetpil,nbp,sizeh,i,ind,fl,flid,v,f,k,j,N,c,lim,av1,va,R1,nba,i1;
  long psim,flagen,fl2,flbug,fl22,fbis;
  GEN fa,fa2,list,ep,y,fap,generateurs,h,pr,prk,prh,prhall,g0,p1,p2,p3,m,pk,p;
  GEN unpadic,pk1,b,x1,den,dinv,ideal,ideal1,zeromod2,unmod2,BIGU,cyclic,arch,rac;
  GEN epr,betae,id;

  nf=checknf(nf);
  if((typ(bigideal)!=17)||(lg(bigideal)!=6))
    err(talker,"not a big ideal in zideallog");
  N=lgef((GEN)nf[1])-3;R1=itos((GEN)((GEN)nf[2])[1]);
  ideal=(GEN)bigideal[1];
  if((typ(ideal)==17)&&(lg(ideal)==3))
  {flid=1;ideal1=(GEN)ideal[1];arch=(GEN)ideal[2];}
  else {flid=0;ideal1=ideal;}
  fa=(GEN)bigideal[3];fa2=(GEN)bigideal[4];list=(GEN)fa[1];ep=(GEN)fa[2];
  nbp=lg(list)-1;id=idmat(N);
  sizeh=0;for(i=1;i<=nbp;i++) sizeh+=(lg((GEN)((GEN)fa2[i])[2])-1);
  if(gcmp0(x)) err(talker,"element is not coprime to ideal in zideallog");
  switch(typ(x))
  {
    case 1: case 4: case 5:
      p1=cgetg(N+1,18);p1[1]=(long)x;for(j=2;j<=N;j++) p1[j]=zero;
      x=p1;break;
    case 9: case 10: x=algtobasis(nf,x);break;
    case 18: break;
    default: err(talker,"not an element in zideallog");
  }
  den=denom(x);
  if(!gcmp1(den))
  {
	/* ici x n'est pas entier, mais on le suppose tel que v_p(x)=0 pour
	   tout p divisant l'ideal. Il faut ecrire x=a/b avec a et b entiers
	   premiers a l'ideal */
    k=1;
    for(i=1;i<=nbp;i++)
    {
      pr=(GEN)list[i];v=ggval(den,(GEN)pr[1])/itos((GEN)ep[i])+1;
      if(v>k) k=v;
    }
    p2=gmul(den,id);p3=idealpow(nf,ideal1,stoi(k));
    p1=idealadd(nf,p2,p3);
    dinv=idealinv(nf,p1);
    p2=idealmullll(nf,p2,dinv);p3=idealmullll(nf,p3,dinv);
    x1=(GEN)idealaddone(nf,p2,p3)[1];
	/* on a trouve x=a/b comme on veut, avec a=x1*x, b=x1 */
    p1=zideallog(nf,element_mul(nf,x1,x),bigideal);
    p2=zideallog(nf,x1,bigideal);
    p1=gsub(p1,p2);
    p2=(GEN)((GEN)bigideal[2])[2];c=lg(p2)-1;
    tetpil=avma;y=cgetg(c+1,18);
    for(i=1;i<=c;i++) y[i]=(long)modii((GEN)p1[i],(GEN)p2[i]);
    return gerepile(av,tetpil,y);
  }
  y=cgetg(sizeh+1,18);ind=0;
      /* Attention: pour l'instant sizeh est la taille de la matrice traitee
	 par smithspec dans zidealstar, c'est a dire ne prenant pas en compte
	 les places archimediennes */
  for(i=1;i<=nbp;i++)
  {
    fap=(GEN)fa2[i];generateurs=(GEN)fap[3];h=(GEN)fap[2];
    pr=(GEN)list[i];p=(GEN)pr[1];psim=itos(p);
    if(gcmp0(x)||element_val(nf,x,pr))
      err(talker,"element is not coprime to ideal in zideallog");
    prk=idealpow(nf,pr,(GEN)ep[i]);f=itos((GEN)pr[4]);
    flagen=((psim>=3)||(f>1));
    if(flagen)
    {
      prh=idealhermite(nf,pr);
      epr=(GEN)pr[3];
      betae=gdiv(element_pow(nf,(GEN)pr[5],epr),gpui((GEN)pr[1],addis(epr,-1),0));
      p1=cgetg(2,19);p1[1]=(long)betae;
      p1=idealadd(nf,gmul((GEN)pr[1],id),idealmul(nf,p1,id));
      prhall=cgetg(3,17);prhall[1]=(long)prh;
      prhall[2]=idealaddone(nf,pr,p1)[2];
      g0=nfreducemodpr(nf,(GEN)generateurs[1],prhall);
      lim=(avma+bot)>>1;av1=avma;
      p1=cgetg(N+1,18);for(j=2;j<=N;j++) p1[j]=zero;
      p1[1]=un;fl=0;
/* Ceci est une recherche la plus bete possible du logarithme discret
   dans Z_K/\goth p. En cas de necessite, c'est tres facile a mettre
   un baby-step giant step

    for(v=0;!fl;v++)
    {
      p2=lift(gsub(p1,x));fl=gcmp0(p2);
      if(!fl) fl=element_val(nf,p2,pr);
      if(!fl) p1=element_mulmodpr(nf,p1,g0,prhall);
      if(avma<lim) {tetpil=avma;p1=gerepile(av1,tetpil,gcopy(p1));}
      if(DEBUGLEVEL>=2) if(!v%1000) {fprintferr("%ld ",v);flusherr();}
    }
    y[++ind]=lstoi(v-1);
*/
      y[++ind]=(long)nfshanks(nf,x,g0,pr,prhall);
    }
    k=itos((GEN)ep[i]);fl22=(psim==2)&&(k==2);fbis=f-fl22;m=cgetg(fbis+1,19);
    if(k>1)
    {
	  /* calcul des autres exposants par log p-adiques. Il faut verifier
	     que c'est correct. Sinon on doit faire ca modulo p^2 et relever
	     par Newton multivariable */
      for(j=1;j<=fbis;j++)
	m[j]=(long)nflogpadic(nf,(GEN)generateurs[j+flagen],pr,prk,k);
      if(flagen)
      {
	g0=(GEN)generateurs[1];pk1=gpuigs(p,k-1);pk=gmul(pk1,p);
	p2=element_powmodideal(nf,g0,addii(negi(addsi(1,(GEN)y[ind])),gpuigs(p,f)),prk);
	p2=element_mulmodideal(nf,x,p2,prk);
      }
      else p2=x;
      p3=nflogpadic(nf,p2,pr,prk,k);unpadic=gadd(gun,ggrandocp(p,k));
      m=gmul(unpadic,m);
      p1=indexrank(m);
      if(lg((GEN)p1[2])!=(fbis+1)) err(talker,"bug in zideallog");
      m=matextract(m,(GEN)p1[1],(GEN)p1[2]);
      p3=extract(p3,(GEN)p1[1]);
      p3=gmul(unpadic,p3);
      b=gauss(m,p3);
      fl2=((psim==2)&&(k>=3));
      flbug=1;for(j=1;(j<fbis)&&flbug;j++) flbug=(padicprec((GEN)b[j],p)<(k-1));
      if(flbug&&fbis) flbug=(padicprec((GEN)b[fbis],p)<(k-1-fl2));
      if(flbug)	err(impl,"zideallog in special case");
      p1=cgetg(fbis+1,18);
      for(j=1;j<fbis;j++) p1[j]=(long)lift(gmul(gmodulcp(gun,pk1),(GEN)b[j]));
      if(fbis)
	p1[fbis]=(long)lift(gmul(gmodulcp(gun,fl2?shifti(pk1,-1):pk1),(GEN)b[j]));
      for(j=1;j<=fbis;j++) y[++ind]=p1[j];
      if(psim==2)
      {
	p3=gzero;
	for(j=1;j<=fbis;j++)
	  p3=gadd(p3,gmul((GEN)p1[j],gsub((GEN)generateurs[j+flagen],(GEN)id[1])));
	p3=gsub(gadd(p3,(GEN)id[1]),p2);
	y[++ind]=(gcmp0(p3)||(element_val(nf,p3,pr)>=2)) ? zero : un;
      }
    }
  }
  cyclic=(GEN)((GEN)bigideal[2])[2];
  if(!flid)
  { /* ici pas de partie archimedienne, on termine */
    c=lg(cyclic)-1;
    p1=gmul((GEN)bigideal[5],y); /* on applique le smith a y */
    tetpil=avma;y=cgetg(c+1,18);
    for(i=1;i<=c;i++) y[i]=(long)modii((GEN)p1[i],(GEN)cyclic[i]);
    return gerepile(av,tetpil,y);
  }
  else
  {
    BIGU=(GEN)bigideal[5];p1=gmul((GEN)BIGU[1],y); /* on applique le smith a y */
    nba=lg((GEN)BIGU[2])-1;c=lg((GEN)BIGU[3])-1;
    if(!nba)
    { /* ici en fait pas de partie archimedienne non plus, on termine */
      tetpil=avma;y=cgetg(c+1,18); 
      for(i=1;i<=c;i++) y[i]=(long)modii((GEN)p1[i],(GEN)cyclic[i]);
      return gerepile(av,tetpil,y);
    }
    else
    { /* ici on a une partie archimedienne */
      y=cgetg(c+1,18);for(i=1;i<=c-nba;i++) y[i]=p1[i];
      p2=cgetg(nba+1,18);p3=gmul((GEN)nf[7],x);rac=(GEN)nf[6];
      va=varn((GEN)nf[1]);
      zeromod2=gmodulcp(gzero,gdeux);unmod2=gmodulcp(gun,gdeux);
      for(i1=0,i=1;i<=R1;i++)
	if(signe((GEN)arch[i]))
	  p2[++i1]=(signe(gsubst(p3,va,(GEN)rac[i]))>0)?(long)zeromod2:(long)unmod2;
	  /* Ceci est la matrice des signatures de notre element x */
      p2=gmul((GEN)BIGU[2],p2);
	  /* idem, mais exprime en fonction de la matrice inversible modulo 2 */
      for(i=c-nba+1;i<=c;i++) y[i]=(long)lift((GEN)p2[i-c+nba]);
      p1=gmul((GEN)BIGU[3],y); /* on applique le deuxieme smith a y */
      c=lg(cyclic)-1;
      tetpil=avma;y=cgetg(c+1,18);
      for(i=1;i<=c;i++) y[i]=(long)modii((GEN)p1[i],(GEN)cyclic[i]);
      return gerepile(av,tetpil,y);
    }
  }
}

long cgcd(long, long);

GEN
subcyclo(GEN p, GEN d)
/* finds an equation for the d-th degree subfield of Q(zeta_p), where p
   must be prime or a power of a prime. */
{
  long av=avma,tetpil,j,k,e,necprec,q,phip,dsim;
  GEN p1,a,om,pol;

  dsim=itos(d);
  necprec=max(5,((dsim+2*gexpo(p))>>TWOPOTBITS_IN_LONG)+2);
  a=gener(p);p1=cgetg(3,6);p1[1]=zero;
  p1[2]=(long)divri(gmul2n(mppi(necprec),1),p);
  om=gexp(p1,necprec);
  phip=itos(phi(p));
  if(gcmpgs(p,phip+1)) err(impl,"subcyclo for prime powers");
  if(phip%dsim) err(talker,"degree does not divide phi(p) in subcyclo");
  q=phip/dsim;
  pol=gun;
  for(k=0;k<dsim;k++)
  {
    p1=gzero;
    for(j=0;j<q;j++)
    {
      e=dsim*j+k;p1=gadd(p1,gpui(om,(GEN)gpuigs(a,e)[2],necprec));
    }
    pol=gmul(pol,gsub(polx[0],p1));
  }
  pol=greal(pol);tetpil=avma;
  return gerepile(av,tetpil,ground(pol));
}
