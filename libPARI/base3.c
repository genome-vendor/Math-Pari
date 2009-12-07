/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                OPERATIONS DANS LES CORPS DE NOMBRES             */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

# include "genpari.h"
GEN nfreducemat(GEN nf, GEN a, GEN b, GEN c, GEN d);
GEN nfgauss(GEN nf, GEN a, GEN b, GEN prhall);
GEN fasthnf(GEN x,GEN detmat);
GEN nfbezout(GEN nf, GEN a, GEN b, GEN ida, GEN idb, GEN *u, GEN *v, GEN *w, GEN *di);
/*******************************************************************

Operations sur les elements et ideaux de corps de nombres.
Un element sera represente par un vecteur colonne dans la base d'entiers nf[7].
Un ideal premier est represente par [p,a,e,f,b] ou l'ideal est p.Z_K+a.Z_K, ou
a est un element de Z_K, e l'indice de ramification, f le degre residuel et b
l'element de Z_K "constante de Lenstra".
Un ideal sera represente par un couple [M,V] ou M est une HNF de l'ideal dans
la base d'entiers, et V un vecteur ligne a r1+r2 composantes complexes
representant la "partie archimedienne" de l'ideal (considere alors comme
idele). Par exemple, si l'ideal est principal
engendre par a, V contiendra le vecteur des log complexes des r1+r2 premiers
conjugues de a (ceci depend bien sur de a et pas seulement de l'ideal).
Les programmes marchent aussi si seulement M est fourni. 

********************************************************************/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                  OPERATIONS SUR LES ELEMENTS                    */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

int
isnfscalar(GEN x)
{
  long lx=lg(x),i;
  
  for(i=2;(i<lx)&&gcmp0((GEN)x[i]);i++);
  return (i==lx);
}

GEN
element_mul(GEN nf, GEN x, GEN y)
                
/* Recoit deux vecteurs de longueur N representant x et y dans la base 
 d'entiers (peut etre modulo p) et ressort leur produit sous forme d'un
 vecteur. */

{
  long av=avma,tetpil,i,j,k,N;
  GEN s,v,c;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if(isnfscalar(x)) return gmul((GEN)x[1],y);
  if(isnfscalar(y)) return gmul((GEN)y[1],x);
  v=cgetg(N+1,18);
  for(k=1;k<=N;k++)
  {
    s=gzero;
    for(i=1;i<=N;i++)
    {
      c=gcoeff((GEN)nf[9],k,(i-1)*N+i);
      if(signe(c))
      {
	if(gcmp1(c)) s=gadd(s,gmul((GEN)x[i],(GEN)y[i]));
	else s=gadd(s,gmul(gmul((GEN)x[i],(GEN)y[i]),c));
      }
      for(j=i+1;j<=N;j++)
      {
	c=gcoeff((GEN)nf[9],k,(i-1)*N+j);
	if(signe(c))
	{
	  if(gcmp1(c)) s=gadd(s,gadd(gmul((GEN)x[i],(GEN)y[j]),gmul((GEN)x[j],(GEN)y[i])));
	  else s=gadd(s,gmul(gadd(gmul((GEN)x[i],(GEN)y[j]),gmul((GEN)x[j],(GEN)y[i])),c));
	}
      }
    }
    v[k]=(long)s;
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(v));
}

GEN
element_inv(GEN nf, GEN x)
                
/* Recoit un vecteurs de longueur N representant x dans la base 
 d'entiers et ressort son inverse sous forme d'un vecteur */
{
  long av=avma,tetpil,flx,i,N;
  GEN p1,p,z,unmod;
  
  N=lgef((GEN)nf[1])-3;
  if(isnfscalar(x))
  {
    z=cgetg(N+1,18);z[1]=(long)ginv((GEN)x[1]);
    for(i=2;i<=N;i++) z[i]=lcopy((GEN)x[i]);
    return z;
  }
  flx=1;for(i=1;(i<=N)&&flx;i++)
  {
    flx=(typ((GEN)x[i])!=3);if(!flx) p=(GEN)((GEN)x[i])[1];
  }
  if(!flx) x=lift(x);
  p1=lift(gdiv(gmodulcp(polun[varn((GEN)nf[1])],(GEN)nf[1]),gmodulcp(gmul((GEN)nf[7],x),(GEN)nf[1])));
  z=cgetg(N+1,18);
  for(i=1;i<=N;i++) z[i]=(long)truecoeff(p1,i-1);
  if(flx)
  {tetpil=avma;return gerepile(av,tetpil,gmul((GEN)nf[8],z));}
  else
  {
    p1=gmul((GEN)nf[8],z);unmod=gmodulcp(gun,p);tetpil=avma;
    return gerepile(av,tetpil,gmul(unmod,p1));
  }
}

GEN
element_div(GEN nf, GEN x, GEN y)
                
/* Recoit deux vecteurs de longueur N representant x et y dans la base 
 d'entiers et ressort leur quotient sous forme d'un vecteur */
{
  long av=avma,tetpil,flx,fly,i,N;
  GEN p1,p,z,unmod;

  nf=checknf(nf);  
  N=lgef((GEN)nf[1])-3;
  if(isnfscalar(y)) return gdiv(x,(GEN)y[1]);
  if(isnfscalar(x))
  {
    p1=element_inv(nf,y);tetpil=avma;
    return gerepile(av,tetpil,gmul((GEN)x[1],p1));
  }
  flx=1;for(i=1;(i<=N)&&flx;i++)
  {
    flx=(typ((GEN)x[i])!=3);if(!flx) p=(GEN)((GEN)x[i])[1];
  }
  if(!flx) x=lift(x);
  fly=1;for(i=1;(i<=N)&&fly;i++)
  {
    fly=(typ((GEN)y[i])!=3);if(!fly) p=(GEN)((GEN)y[i])[1];
  }
  if(!fly) y=lift(y);

  p1=lift(gdiv(gmodulcp(gmul((GEN)nf[7],x),(GEN)nf[1]),gmodulcp(gmul((GEN)nf[7],y),(GEN)nf[1])));
  z=cgetg(N+1,18);
  for(i=1;i<=N;i++) z[i]=(long)truecoeff(p1,i-1);
  if(flx&&fly)
  {tetpil=avma;return gerepile(av,tetpil,gmul((GEN)nf[8],z));}
  else
  {
    p1=gmul((GEN)nf[8],z);unmod=gmodulcp(gun,p);tetpil=avma;
    return gerepile(av,tetpil,gmul(unmod,p1));
  }
}

GEN
element_muli(GEN nf, GEN x, GEN y)
                
/* Recoit deux vecteurs de longueur N representant x et y dans la base 
 d'entiers (uniquement coeff entiers) et ressort leur produit sous forme
 d'un vecteur */

{
  long av=avma,tetpil,i,j,k,N=lgef((GEN)nf[1])-3;
  GEN s,v,c;

  v=cgetg(N+1,18);
  for(k=1;k<=N;k++)
  {
    s=gzero;
    for(i=1;i<=N;i++)
    {
      c=(GEN)coeff((GEN)nf[9],k,(i-1)*N+i);
      if(signe(c))
      {
	if(gcmp1(c)) s=addii(s,mulii((GEN)x[i],(GEN)y[i]));
	else s=addii(s,mulii(mulii((GEN)x[i],(GEN)y[i]),c));
      }
      for(j=i+1;j<=N;j++)
      {
	c=(GEN)coeff((GEN)nf[9],k,(i-1)*N+j);
	if(signe(c))
	{
	  if(gcmp1(c)) s=addii(s,addii(mulii((GEN)x[i],(GEN)y[j]),mulii((GEN)x[j],(GEN)y[i])));
	  else s=addii(s,mulii(addii(mulii((GEN)x[i],(GEN)y[j]),mulii((GEN)x[j],(GEN)y[i])),c));
	}
      }
    }
    v[k]=(long)s;
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(v));
}

GEN
element_mulh(GEN nf, long limi, long limj, GEN x, GEN y)

/* Recoit deux vecteurs de longueur N representant x et y dans la base 
 d'entiers (uniquement coeff entiers) et ressort leur produit sous forme d'un
 vecteur */

{
  long av=avma,tetpil,i,j,k,N=lgef((GEN)nf[1])-3;
  GEN s,v,c;

  if(limi<limj) {i=limi;limi=limj;limj=i;s=x;x=y;y=s;}
  v=cgetg(N+1,18);
  for(k=1;k<=N;k++)
  {
    s=gzero;
    for(i=1;i<=limj;i++)
    {
      c=(GEN)coeff((GEN)nf[9],k,(i-1)*N+i);
      if(signe(c))
      {
	if(gcmp1(c)) s=addii(s,mulii((GEN)x[i],(GEN)y[i]));
	else s=addii(s,mulii(mulii((GEN)x[i],(GEN)y[i]),c));
      }
      for(j=i+1;j<=limj;j++)
      {
	c=(GEN)coeff((GEN)nf[9],k,(i-1)*N+j);
	if(signe(c))
	{
	  if(gcmp1(c)) s=addii(s,addii(mulii((GEN)x[i],(GEN)y[j]),mulii((GEN)x[j],(GEN)y[i])));
	  else s=addii(s,mulii(addii(mulii((GEN)x[i],(GEN)y[j]),mulii((GEN)x[j],(GEN)y[i])),c));
	}
      }
    }
    for(i=limj+1;i<=limi;i++)
    {
      for(j=1;j<=limj;j++)
      {
	c=(GEN)coeff((GEN)nf[9],k,(i-1)*N+j);
	if(signe(c))
	{
	  if(gcmp1(c)) s=addii(s,mulii((GEN)x[i],(GEN)y[j]));
	  else s=addii(s,mulii(mulii((GEN)x[i],(GEN)y[j]),c));
	}
      }
    }
    v[k]=(long)s;
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(v));
}

GEN
element_sqr(GEN nf, GEN x)
              
/* Recoit un vecteurs de longueur N representant x dans la base 
 d'entiers et ressort son carre sous forme d'un
 vecteur */

{
  long av=avma,tetpil,i,j,k,N=lgef((GEN)nf[1])-3;
  GEN s,v,c;

  if(isnfscalar(x))
  {
    s=cgetg(N+1,18);s[1]=(long)gsqr((GEN)x[1]);
    for(i=2;i<=N;i++) s[i]=lcopy((GEN)x[i]);
    return s;
  }
  v=cgetg(N+1,18);
  for(k=1;k<=N;k++)
  {
    s=gzero;
    for(i=1;i<=N;i++)
    {
      c=(GEN)coeff((GEN)nf[9],k,(i-1)*N+i);
      if(signe(c))
      {
	if(gcmp1(c)) s=gadd(s,gmul((GEN)x[i],(GEN)x[i]));
	else s=gadd(s,gmul(gmul((GEN)x[i],(GEN)x[i]),c));
      }
      for(j=i+1;j<=N;j++)
      {
	c=(GEN)coeff((GEN)nf[9],k,(i-1)*N+j);
	if(signe(c))
	{
	  if(gcmp1(c)) s=gadd(s,gmul2n(gmul((GEN)x[i],(GEN)x[j]),1));
	  else s=gadd(s,gmul(gmul((GEN)x[i],(GEN)x[j]),shifti(c,1)));
	}
      }
    }
    v[k]=(long)s;
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(v));
}

GEN
element_pow(GEN nf, GEN x, GEN k)
                
/* Calcule x^k dans le corps de nombres nf */

{
  long i,f,s,av=avma,tetpil,N;
  GEN k1,y,z;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if(typ(k)!=1) err(talker,"not an integer exponent in nfpow");
  if(isnfscalar(x))
  {
    z=cgetg(N+1,18);z[1]=(long)gpui((GEN)x[1],k,0);
    for(i=2;i<=N;i++) z[i]=lcopy((GEN)x[i]);
    return z;
  }
  s=signe(k);k1=(s>=0)?k:negi(k);z=x;f=1;y=cgetg(N+1,18);
  for(i=2;i<=N;i++) y[i]=zero;
  y[1]=un;
  while(f)
  {
    if(mpodd(k1)) y=element_mul(nf,z,y);
    k1=shifti(k1,-1);f=signe(k1);
    if(f) z=element_sqr(nf,z);
  }
  tetpil=avma;return gerepile(av,tetpil,(s>=0)?gcopy(y):element_inv(nf,y));
}

GEN
element_mulid(GEN nf, GEN x, long i)
                
/* Recoit un vecteur x de longueur N representant x dans la base d'entiers
(peut etre modulo p) et ressort x.w_i sous forme d'un
 vecteur */

{
  long av=avma,tetpil,j,k,N=lgef((GEN)nf[1])-3;
  GEN s,v,c,p1;

  v=cgetg(N+1,18);
  for(k=1;k<=N;k++)
  {
    s=gzero;
    for(j=1;j<=N;j++)
    {
      c=gcoeff((GEN)nf[9],k,(i-1)*N+j);p1=(GEN)x[j];
      if((signe(c))&&(!gcmp0(p1)))
	s=gcmp1(c)?gadd(s,p1):gadd(s,gmul(p1,c));
    }
    v[k]=(long)s;
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(v));
}


long
element_val(GEN nf, GEN x, GEN vp)
{
  long av=avma,N,w,bo,i,vd,v,e=itos((GEN)vp[3]);
  GEN denx,p1,p=(GEN)vp[1],bp,d,r;

  nf=checknf(nf);
  if((typ(vp)!=17)||(lg(vp)!=6)) err(idealer3);
  N=lgef((GEN)nf[1])-3;
  switch(typ(x))
  {
    case 1: case 4: case 5:
      return ggval(x,p)*e;
    case 10: d=gabs(subres(x,(GEN)nf[1]),0);
      x=(GEN)principalideal(nf,x)[1];break;
    case 9: d=gabs(subres((GEN)x[2],(GEN)nf[1]),0);
      x=(GEN)principalideal(nf,x)[1];break;
    case 18:
      if(lg(x)==N+1)
      {d=gabs(subres(gmul((GEN)nf[7],x),(GEN)nf[1]),0);break;}
    default: err(talker,"incorrect type in element_val");
  }
  if(isnfscalar(x)) return ggval((GEN)x[1],p)*e;
  denx=denom(x);
  if(!gcmp1(denx)) {x=gmul(denx,x);vd=ggval(denx,p);} else vd=0;
  v=ggval(d,p)+N*vd;
  if(!v) return -vd*e;
  bo=0;w=0;bp=(GEN)vp[5];
  do
  {
    x=element_muli(nf,x,bp);
    if(divise((GEN)x[N],p))
    {
      for(i=1;i<=N;i++)
      {
	p1=dvmdii((GEN)x[i],p,&r);
	if(signe(r)) goto labelelv; else x[i]=(long)p1;
      }
      w++;
    }
    else bo=1;
  }
  while((bo==0)&&(w<v));
  labelelv:
  avma=av;return w-vd*e;
}

long
element_val2(GEN nf, GEN x, GEN d, GEN vp)
/* a usage interne, pas de verifs. */
{
  long av=avma,N=lgef((GEN)nf[1])-3,w,bo,i,v;
  GEN p1,p=(GEN)vp[1],bp,r;

  v=ggval(d,p);
  if(!v) return 0;
  bo=0;w=0;bp=(GEN)vp[5];
  do
  {
    x=element_muli(nf,x,bp);
    if(divise((GEN)x[N],p))
    {
      for(i=1;i<=N;i++)
      {
	p1=dvmdii((GEN)x[i],p,&r);
	if(signe(r)) goto labelelv2; else x[i]=(long)p1;
      }
      w++;
    }
    else bo=1;
  }
  while((bo==0)&&(w<v));
  labelelv2:
  avma=av;return w;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                   OPERATIONS SUR LES IDEAUX                     */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

GEN
ideal_two_elt(GEN nf, GEN ix)
/* Etant donne un ideal ix, ressort un vecteur [a,alpha] a deux
   composantes tel que a soit rationnel et ix=aZ_K+alpha Z_K, alpha etant
   un vecteur colonne sur la base d'entiers. On peut avoir a=0 ou
   alpha=0, mais on ne cherche pas a determiner si ix est principal.
   */
{
  long av,av1,tetpil,tx=typ(ix),r,i,j,fl,N;
  GEN p1,z,alpha,beta,lambda,norme,pnorm,ideal,den,idz;
  
  nf=checknf(nf);  
  N=lgef((GEN)nf[1])-3;
  if((tx==17)&&(lg(ix)==3)) {ix=(GEN)ix[1];tx=typ(ix);}
  av=avma;z=cgetg(3,17);
  if(tx<=10)
  {
    switch(tx)
    {
      case 1: case 4: case 5: z[1]=lcopy(ix);
	p1=cgetg(N+1,18);z[2]=(long)p1;for(i=1;i<=N;i++) p1[i]=zero;
	return z;
      case 9: if(!gegal((GEN)nf[1],(GEN)ix[1])) 
	err(talker,"incompatible number fields in ideal_two_elt");
	ix=(GEN)ix[2];   /* fall through */
      case 10: z[1]=zero;av=avma;p1=cgetg(N+1,18);
	for(i=1;i<=N;i++) p1[i]=(long)truecoeff(ix,i-1);
	tetpil=avma;z[2]=lpile(av,tetpil,gmul((GEN)nf[8],p1));
	return z;
      default: err(talker,"incorrect type in ideal_two_elt");
	return gzero;
    }
  }
  else
  {
    switch(tx)
    {
      case 17:
	if(lg(ix)==6)
	{z[1]=lcopy((GEN)ix[1]);z[2]=lcopy((GEN)ix[2]);return z;}
	else err(talker,"incorrect type in ideal_two_elt");break;
      case 18:
	if(lg(ix)==N+1)
	{z[1]=zero;z[2]=lcopy(ix);return z;}
	else err(talker,"incorrect type in ideal_two_elt");break;
      case 19:
	if(lg(ix)==1)
	{
	  z[1]=zero;p1=cgetg(N+1,18);z[2]=(long)p1;
	  for(i=1;i<=N;i++) p1[i]=zero;
	  return z;
	}
	if(lg((GEN)ix[1])!=N+1)
	  err(talker,"incorrect type in ideal_two_elt");
	if(lg(ix)==2)
	{z[1]=zero;z[2]=lcopy((GEN)ix[1]);return z;}
	ideal=idealhermite(nf,ix);den=denom(ideal);
	ideal=gmul(ideal,den);idz=gcoeff(ideal,1,1);
	pnorm=idealnorm(nf,ideal);
	beta=gmodulcp(gmul((GEN)nf[7],ideal),(GEN)nf[1]);
	fl=r=1;
	for(i=2;(i<=N)&&fl;i++)
	{
	  alpha=(GEN)beta[i];norme=gnorm(alpha);
	  if(!cmpii(mppgcd(divii(norme,pnorm),pnorm),gun)) fl=0;
	  else
	  {
	    alpha=gadd(alpha,idz);norme=gnorm(alpha);
	    if(!cmpii(mppgcd(divii(norme,pnorm),pnorm),gun)) fl=0;
	  }
	}
	if(fl)
	{
	  lambda=cgeti(N+1);
	  av1=avma;
	  for(i=1;i<=N;i++) lambda[i]=r;
	  do
	  {
	    avma=av1;
	    alpha=gmodulcp(gzero,(GEN)nf[1]);
	    for(i=1;i<=N;i++)
	      alpha=gadd(alpha,gmulsg(lambda[i],(GEN)beta[i]));
	    norme=gnorm(alpha);
	    if(!cmpii(mppgcd(divii(norme,pnorm),pnorm),gun)) fl=0;
	    else
	    {
	      alpha=gadd(alpha,idz);norme=gnorm(alpha);
	      if(!cmpii(mppgcd(divii(norme,pnorm),pnorm),gun)) fl=0;
	    }
	    if(fl)
	    {
	      for(j=N;(lambda[j]+r)==0;j--);
	      lambda[j]--;
	      for(i=j+1;i<=N;i++) lambda[i]=r;
	      for(j=1;(j<N)&&(!lambda[j]);j++);
	      if(!lambda[j])
	      {
		r++;for(i=1;i<=N;i++) lambda[i]=r;
		if(cmpis(idz,(r<<1))<0)
		  err(talker,"ideal_two_elt fails");
	      }
	    }
	  }
	  while(fl);
	}
	alpha=lift(alpha);beta=cgetg(N+1,18);
	for(i=1;i<=N;i++) beta[i]=(long)truecoeff(alpha,i-1);
	alpha=gmul((GEN)nf[8],beta);
	alpha=gmul(gmodulcp(gun,idz),alpha);
	z[1]=ldiv(idz,den);
	z[2]=ldiv(centerlift(alpha),den);
	tetpil=avma;return gerepile(av,tetpil,gcopy(z));
      default: err(talker,"incorrect type in ideal_two_elt");
    }
  }
  return gnil;
}

GEN
idealfactor(GEN nf, GEN x)
{
  long av=avma,tetpil,i,j,k,lf,lff,N,ls,v,vd;
  GEN d,f,f1,f2,ff,ff1,ff2,y1,y2,y,p1,p2,denx;

  nf=checknf(nf);  
  if((typ(x)<=10)||(typ(x)==18)) x=principalideal(nf,x);
  if((typ(x)==17)&&(lg(x)==3)) x=(GEN)x[1];
  if(typ(x)!=19) err(idealer2);
  N=lgef((GEN)nf[1])-3;if(lg(x)!=(N+1)) x=idealmul(nf,x,idmat(N));
  denx=denom(x);if(!gcmp1(denx)) x=gmul(denx,x);
  for(d=gun,i=1;i<=N;i++) d=mulii(d,gcoeff(x,i,i));
  f=factor(absi(d));f1=(GEN)f[1];f2=(GEN)f[2];lf=lg(f1);
  if(!gcmp1(denx))
  {ff=factor(denx);ff1=(GEN)ff[1];ff2=(GEN)ff[2];lff=lg(ff1);}
  else lff=1;
  y1=cgetg((lf+lff-2)*N+1,18);y2=cgeti((lf+lff-2)*N+1);k=0;
  for(i=1;i<lf;i++)
  {
    p1=primedec(nf,(GEN)f1[i]);ls=itos((GEN)f2[i]);
    vd=ggval(denx,(GEN)f1[i]);
    for(j=1;j<lg(p1);j++)
    {
      p2=(GEN)p1[j];
      if(ls)
      {
	v=idealval(nf,x,p2);
	ls-=(v*itos((GEN)p2[4]));
	v-=vd*itos((GEN)p2[3]);
      }
      else v=-vd*itos((GEN)p2[3]);
      if(v) {y1[++k]=(long)p2;y2[k]=v;}
    }
  }
  if(!gcmp1(denx))
  {
    for(i=1;i<lff;i++)
    {
      if(!divise(d,(GEN)ff1[i]))
      {
	p1=primedec(nf,(GEN)ff1[i]);
	for(j=1;j<lg(p1);j++)
	{
	  p2=(GEN)p1[j];y1[++k]=(long)p2;y2[k]=-itos((GEN)ff2[i])*itos((GEN)p2[3]);
	}
      }
    }
  }
  tetpil=avma;
  y=cgetg(3,19);p1=cgetg(k+1,18);p2=cgetg(k+1,18);y[1]=(long)p1;y[2]=(long)p2;
  for(i=1;i<=k;i++) {p1[i]=lcopy((GEN)y1[i]);p2[i]=lstoi(y2[i]);}
  return gerepile(av,tetpil,y);
}

GEN
idealadd(GEN nf, GEN x, GEN y)
{
  long av=avma,tetpil,N;
  GEN z,p1,denx,deny,denz;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if((typ(x)<=18)||((typ(x)==19)&&(lg(x)!=N+1))) x=idealhermite(nf,x);
  if((typ(y)<=18)||((typ(y)==19)&&(lg(y)!=N+1))) y=idealhermite(nf,y);
  denx=denom(x);deny=denom(y);denz=gmul(denx,deny);
  if(!gcmp1(denz)) {x=gmul(x,denz);y=gmul(y,denz);}
  p1=ggcd(detint(x),detint(y));
  if(gcmp1(p1))
  {
    if(gcmp1(denz)) {avma=av;return idmat(N);}
    else
    {p1=ginv(denz);tetpil=avma;return gerepile(av,tetpil,gscalmat(p1,N));}
  }
  else
  {
    z=concat(x,y);tetpil=avma;z=hnfmod(z,p1);
    if(gcmp1(denz)) return gerepile(av,tetpil,z);
    else {tetpil=avma;return gerepile(av,tetpil,gdiv(z,denz));}
  }
}

int
ishnfall(GEN x)
{
  long i,j;
  for(i=2;i<lg(x);i++)
    for(j=1;j<i;j++)
      if(!gcmp0(gcoeff(x,i,j))) return 0;
  return 1;
}

GEN hermiteidealmul(GEN nf, GEN ix, GEN iy);

GEN
idealhermite(GEN nf, GEN x)
{
  long av=avma,tetpil,N;
  GEN z,dp,id;
  
  nf=checknf(nf);N=lgef((GEN)nf[1])-3;  
  if((typ(x)==17)&&(lg(x)==3)) x=(GEN)x[1];
  id=idmat(N);
  if((typ(x)<=10)||(typ(x)==18))
  {
    x=principalideal(nf,x);tetpil=avma;
    return gerepile(av,tetpil,hermiteidealmul(nf,x,id));
  }
  if((typ(x)==17)&&(lg(x)==6)) {avma=av;return idealmulprime(nf,gun,x);}
  if(typ(x)!=19) err(idealer2);
  if((lg(x)==1)||gcmp0(x)) {avma=av;return gscalmat(gzero,N);}
  if(lg((GEN)x[1])!=N+1) err(idealer2);
  if(lg(x)>N+1)
  {
    z=denom(x);if(!gcmp1(z)) x=gmul(z,x);
    dp=detint(x);tetpil=avma;x=hnfmod(x,dp);
    if(!gcmp1(z)) {tetpil=avma;return gerepile(av,tetpil,gdiv(x,z));}
    else return gerepile(av,tetpil,x);
  }
  if(lg(x)<=N)
  {tetpil=avma;return gerepile(av,tetpil,hermiteidealmul(nf,x,id));}
  if(ishnfall(x)) {tetpil=avma;return gerepile(av,tetpil,gcopy(x));}
  z=denom(x);
  if(!gcmp1(z))
  {
    x=gmul(z,x);x=hnfmod(x,detint(x));tetpil=avma;
    return gerepile(av,tetpil,gdiv(x,z));
  }
  else
  {
    dp=detint(x);tetpil=avma;
    return gerepile(av,tetpil,hnfmod(x,dp));
  }
}

GEN
idealhermite2(GEN nf, GEN a, GEN b)
/* HNF of aZ_K+bZ_K */

{
  long av=avma,tetpil,N;
  GEN p1,id;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if((typ(a)<=10)||(typ(a)==18)) a=principalideal(nf,a);
  if((typ(a)!=19)||(lg(a)!=2))
    err(talker,"incorrect type in idealhermite2");
  if((typ(b)<=10)||(typ(b)==18)) b=principalideal(nf,b);
  if((typ(b)!=19)||(lg(b)!=2))
    err(talker,"incorrect type in idealhermite2");
  p1=concat(a,b);id=idmat(N);tetpil=avma;
  return gerepile(av,tetpil,idealmul(nf,id,p1));
}

int
ishnf1(GEN x)
{
  long i,lx;
  
  if(lg(x)==1) return 0;
  x=(GEN)x[1];lx=lg(x);
  for(i=2;(i<lx)&&gcmp0((GEN)x[i]);i++);
  return (i>=lx);
}

GEN
idealaddone(GEN nf, GEN x, GEN y)
{
  long av=avma,tetpil,N,i,j,flx=0,fly=0,flz=0;
  GEN v,v1,v2,v3,v4,p1,p3,id,u,unnf;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  id=idmat(N);
  if(vecegal(x,id))
  {
    avma=av;v=cgetg(3,17);p1=cgetg(N+1,18);v[1]=(long)p1;
    p1[1]=un;for(i=2;i<=N;i++) p1[i]=zero;
    p1=cgetg(N+1,18);v[2]=(long)p1;for(i=1;i<=N;i++) p1[i]=zero;
    return v;
  }
  if(vecegal(y,id))
  {
    avma=av;v=cgetg(3,17);p1=cgetg(N+1,18);v[2]=(long)p1;
    p1[1]=un;for(i=2;i<=N;i++) p1[i]=zero;
    p1=cgetg(N+1,18);v[1]=(long)p1;for(i=1;i<=N;i++) p1[i]=zero;
    return v;
  }
  if(DEBUGLEVEL>2)
  {
    fprintferr(" entree dans idealaddone() :\n");
    fprintferr(" x = ");outerr(x);fprintferr(" y = ");outerr(y);
    flusherr();
  }
  if((typ(x)!=19)||(lg(x)!=lg((GEN)x[1])))
  {x=idealhermite(nf,x);flx=1;}
  if((typ(y)!=19)||(lg(y)!=lg((GEN)y[1])))  
  {y=idealhermite(nf,y);fly=1;}
  if((flx||ishnf1(x))&&(fly||ishnf1(y)))
  {
    if(gcmp1(bezout(gcoeff(x,1,1),gcoeff(y,1,1),&u,&v)))
    {
      p1=gmul(u,(GEN)x[1]);flz=1;
      unnf=cgetg(N+1,18);unnf[1]=un;
      for(i=2;i<=N;i++) unnf[i]=zero;
    }
  }
  if(!flz)
  {
    v=hnfperm(concat(x,y));v1=(GEN)v[1];v2=(GEN)v[2];v3=(GEN)v[3];j=0;
    for(i=1;i<=N;i++)
    {
      if(!gcmp1(gcoeff(v1,i,i)))
	err(talker,"ideals don't sum to Z_K in idealaddone");
      if(gcmp1((GEN)v3[i])) j=i;
    }
    v4=(GEN)v2[N+j];p1=cgetg(N+1,18);
    for(i=1;i<=N;i++) p1[i]=v4[i];
    p1=gmul(x,p1);unnf=(GEN)v1[1];
  }
  p3=idealmullll(nf,x,y);
  tetpil=avma;v=cgetg(3,17);v[1]=(long)element_reduce(nf,p1,p3);
  v[2]=lsub(unnf,(GEN)v[1]);
  if(DEBUGLEVEL>2)
  {fprintferr(" sortie de idealaddone : v = ");outerr(v);flusherr();}
  return gerepile(av,tetpil,v);
}

GEN
ideleaddone(GEN nf,GEN x,GEN idele)
{
  long av=avma,tetpil,i,R1,va;
  GEN n,k,u,v,p1,p2,p3,ideal,rac,arch;

  nf=checknf(nf);
  if((typ(idele)==17)&&(lg(idele)==2))
  {
    R1=itos((GEN)((GEN)nf[2])[1]);
    ideal=(GEN)idele[1];arch=(GEN)idele[2];
    if((typ(arch)!=17)&&(typ(arch)!=18)&&(lg(arch)!=R1+1))
      err(talker,"incorrect idele in ideleaddone");
    p1=idealaddone(nf,x,ideal);u=gmul((GEN)nf[7],(GEN)p1[1]);
    rac=(GEN)nf[6];va=varn((GEN)nf[1]);n=gcoeff(idealhermite(nf,ideal),1,1);
    k=gzero;
    for(i=1;i<=R1;i++)
    {
      if(signe((GEN)arch[i]))
      {p2=gsubst(u,va,(GEN)rac[i]);if(signe(p2)<0) k=gmax(k,gceil(gdiv(gneg(p2),n)));}
    }
    p3=gmul(k,n);
    u=(GEN)p1[1];v=(GEN)p1[2];u[1]=ladd((GEN)u[1],p3);v[1]=lsub((GEN)v[1],p3);
    tetpil=avma;return gerepile(av,tetpil,gcopy(p1));
  }
  else return idealaddone(nf,x,idele);
}

GEN
idealaddmultone(GEN nf, GEN list)
{
  long av=avma,tetpil,N,i,i1,j,k;
  GEN z,v,v1,v2,v3,p1;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if(DEBUGLEVEL>1)
  {
    fprintferr(" entree dans idealaddmultone() :\n");
    fprintferr(" list = ");outerr(list);flusherr();
  }
  if((typ(nf)!=17)||(lg(nf)<10)) err(idealer1);
  if((typ(list)!=17)&&(typ(list)!=18))
    err(talker,"not a list in idealaddmultone");
  k=lg(list);
  if(k==1) err(talker,"ideals don't sum to Z_K in idealaddmultone");
  p1=(GEN)list[1];
  if((typ(p1)!=19)||(lg(p1)!=lg((GEN)p1[1])))
    list[1]=(long)idealhermite(nf,p1);
  z=(GEN)list[1];
  for(i=2;i<k;i++)
  {
    p1=(GEN)list[i];
    if((typ(p1)!=19)||(lg(p1)!=lg((GEN)p1[1])))
      list[i]=(long)idealhermite(nf,p1);
    z=concat(z,(GEN)list[i]);
  }
  v=hnfperm(z);v1=(GEN)v[1];v2=(GEN)v[2];v3=(GEN)v[3];j=0;
  for(i=1;i<=N;i++)
  {
    if(!gcmp1(gcoeff(v1,i,i)))
      err(talker,"ideals don't sum to Z_K in idealaddmultone");
    if(gcmp1((GEN)v3[i])) j=i;
  }
  v=(GEN)v2[(k-2)*N+j];
  z=cgetg(k,17);
  for(i=1;i<k;i++)
  {
    p1=cgetg(N+1,18);z[i]=(long)p1;
    for(i1=1;i1<=N;i1++) p1[i1]=v[(i-1)*N+i1];
  }
  tetpil=avma;v=cgetg(k,typ(list));
  for(i=1;i<k;i++) v[i]=lmul((GEN)list[i],(GEN)z[i]);
  if(DEBUGLEVEL>1)
  {fprintferr(" sortie de idealaddmultone v = ");outerr(v);flusherr();}
  return gerepile(av,tetpil,v);
}  

	
#define isideletype(x) ((typ(x)==17)&&(lg(x)==3))
#define principalid(nf,x,y) (((typ(y)==17)&&(lg(y)==3))?principalidele(nf,x):principalideal(nf,x))

GEN
idealmul(GEN nf, GEN ix, GEN iy)
                  
/*
  recoit deux ideaux ix et iy comme ci-dessus avec ou sans leur composante
  archimedienne et ressort leur produit sans le reduire */

{
  long av=avma,f,tetpil,rx,ry,tx=typ(ix),ty=typ(iy),N,i,j;
  GEN x,y,dx,dy,dz,y1,m,dp;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if((tx<=10)||(tx==18)) {ix=principalid(nf,ix,iy);tx=typ(ix);}
  if((ty<=10)||(ty==18)) {iy=principalid(nf,iy,ix);ty=typ(iy);}
  if((tx==17)&&(lg(ix)==6)) 
  {
    if((ty==17)&&(lg(iy)==6))
    {
      x=idealmulprime(nf,idmat(N),ix);
      tetpil=avma;
      return gerepile(av,tetpil,idealmulprime(nf,x,iy));
    }
    else return idealmulprime(nf,iy,ix);
  }
  if((ty==17)&&(lg(iy)==6)) return idealmulprime(nf,ix,iy);
  f=1;
  if(tx==17) {if(lg(ix)==3) x=(GEN)ix[1];else err(idealer2);}
  else {f=0;x=ix;}
  if(ty==17) {if(lg(iy)==3) y=(GEN)iy[1];else err(idealer2);}
  else {f=0;y=iy;}
  if((typ(x)!=19)||(typ(y)!=19)) err(idealer2);
  rx=lg(x)-1;ry=lg(y)-1;
  if((!rx)||(!ry)) {tetpil=avma;y=gscalmat(gzero,N);}
  else
  {
    dx=denom(x);dy=denom(y);x=gmul(dx,x);y=gmul(dy,y);dz=mulii(dx,dy);
    if((rx<=2)||(ry<=2))
    {
      m=cgetg(rx*ry+1,19);
      for(i=1;i<=rx;i++)
	for(j=1;j<=ry;j++)
	  m[(i-1)*ry+j]=(long)element_muli(nf,(GEN)x[i],(GEN)y[j]);
      dp=detint(m);tetpil=avma;y=hnfmod(m,dp);
    }
    else
    {
      x=idealhermite(nf,x);y=idealhermite(nf,y);
      tetpil=avma;y=idealmulh(nf,x,y);
    }
    if(!gcmp1(dz)) {tetpil=avma;y=gdiv(y,dz);}
  }
  if(f) 
  {
    y1=cgetg(3,17);y1[1]=(long)y;
    y1[2]=ladd((GEN)ix[2],(GEN)iy[2]);
    return gerepile(av,tetpil,y1);
  }
  else return gerepile(av,tetpil,y);
}

GEN
hermiteidealmul(GEN nf, GEN ix, GEN iy)
                  
/*
  recoit deux ideaux ix et iy comme ci-dessus avec ou sans leur composante
  archimedienne et ressort leur produit sans le reduire. Version
  sans ideal_two_elt, indispensable dans idealhermite. */

{
  long av=avma,f,tetpil,rx,ry,i,j,tx=typ(ix),ty=typ(iy),N;
  GEN m,x,y,dp,dx,dy,dz,y1;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if((tx<=10)||(tx==18)) {ix=principalid(nf,ix,iy);tx=typ(ix);}
  if((ty<=10)||(ty==18)) {iy=principalid(nf,iy,ix);ty=typ(iy);}
  if((tx==17)&&(lg(ix)==6)) 
  {
    if((ty==17)&&(lg(iy)==6))
    {
      x=idealmulprime(nf,idmat(N),ix);
      tetpil=avma;
      return gerepile(av,tetpil,idealmulprime(nf,x,iy));
    }
    else return idealmulprime(nf,iy,ix);
  }
  if((ty==17)&&(lg(iy)==6)) return idealmulprime(nf,ix,iy);
  f=1;
  if(tx==17) {if(lg(ix)==3) x=(GEN)ix[1];else err(idealer2);}
  else {f=0;x=ix;}
  if(ty==17) {if(lg(iy)==3) y=(GEN)iy[1];else err(idealer2);}
  else {f=0;y=iy;}
  if((typ(x)!=19)||(typ(y)!=19)) err(idealer2);
  rx=lg(x)-1;ry=lg(y)-1;
  if((!rx)||(!ry)) {tetpil=avma;y=gscalmat(gzero,N);}
  else
  {
    dx=denom(x);dy=denom(y);x=gmul(dx,x);y=gmul(dy,y);
    dz=mulii(dx,dy);
    m=cgetg(rx*ry+1,19);
    for(i=1;i<=rx;i++)
      for(j=1;j<=ry;j++)
	m[(i-1)*ry+j]=(long)element_muli(nf,(GEN)x[i],(GEN)y[j]);
    if((rx==lg((GEN)x[1])-1)&&(ry==(lg((GEN)y[1])-1)))
    {
      dp=mulii(det(x),det(y));tetpil=avma;
      y=hnfmod(m,dp);
    }
    else {dp=detint(m);tetpil=avma;y=hnfmod(m,dp);}
    if(!gcmp1(dz)) {tetpil=avma;y=gdiv(y,dz);}
  }
  if(f) 
  {
    y1=cgetg(3,17);y1[1]=(long)y;
    y1[2]=ladd((GEN)ix[2],(GEN)iy[2]);
    return gerepile(av,tetpil,y1);
  }
  else return gerepile(av,tetpil,y);
}

GEN
idealmulred(GEN nf, GEN ix, GEN iy, long prec)
{
  long av=avma,tetpil;
  GEN p1;

  p1=idealmul(nf,ix,iy);tetpil=avma;
  return gerepile(av,tetpil,ideallllred(nf,p1,gzero,prec));
}

GEN idealmulspec(GEN nf, GEN x, GEN spec);

GEN
idealmulh(GEN nf, GEN ix, GEN iy)
                  
/* recoit deux ideaux ix et iy comme ci-dessus avec ou sans leur composante
archimedienne et ressort leur produit sans le reduire. On suppose les ideaux
sous forme HNF et de meme taille. A usage interne donc aucune verification. */

{
  long av=avma,f,tetpil,N,i;
  GEN m,x,y,dy;
  
  f=1;
  if(typ(ix)==17) x=(GEN)ix[1];else {f=0;x=ix;}
  if(typ(iy)==17) y=(GEN)iy[1];else {f=0;y=iy;}
  N=lg(x)-1;
  dy=gcoeff(y,1,1);for(i=2;i<=N;i++) dy=mulii(dy,gcoeff(y,i,i));
  y=ideal_two_elt(nf,y);
  m=cgetg(4,17);m[1]=y[1];m[2]=y[2];m[3]=(long)dy;
  tetpil=avma;
  if(f)
  {
    y=cgetg(3,17);y[1]=(long)idealmulspec(nf,x,m);
    y[2]=ladd((GEN)ix[2],(GEN)iy[2]);
  }
  else y=idealmulspec(nf,x,m);
  return gerepile(av,tetpil,y);
}


GEN
hermiteidealmulh(GEN nf, GEN ix, GEN iy)
                  
/* recoit deux ideaux ix et iy comme ci-dessus avec ou sans leur composante
   archimedienne et ressort leur produit sans le reduire. On suppose les ideaux
   sous forme HNF et de meme taille. A usage interne donc aucune verification.
   Ancienne version sans ideal_two_elt. */

{
  long av=avma,f,tetpil,N,i,j;
  GEN m,x,y,dx,dy,dz;
  
  f=1;
  if(typ(ix)==17) x=(GEN)ix[1];else {f=0;x=ix;}
  if(typ(iy)==17) y=(GEN)iy[1];else {f=0;y=iy;}
  N=lg(x)-1;m=cgetg(N*N+1,19);dx=gcoeff(x,1,1);dy=gcoeff(y,1,1);
  for(i=2;i<=N;i++) {dx=mulii(dx,gcoeff(x,i,i));dy=mulii(dy,gcoeff(y,i,i));}
  dz=mulii(dx,dy);
  for(i=1;i<=N;i++)
    for(j=1;j<=N;j++)
      m[(N-i)*N+j]=(long)element_mulh(nf,i,j,(GEN)x[i],(GEN)y[j]);
  tetpil=avma;
  if(f) {y=cgetg(3,17);y[1]=(long)fasthnf(m,dz);y[2]=ladd((GEN)ix[2],(GEN)iy[2]);}
  else y=fasthnf(m,dz);
  return gerepile(av,tetpil,y);
}

GEN
idealmulprime(GEN nf, GEN ix, GEN vp)
                  
/*
  recoit un ideal ix et un ideal premier vp en format
  primedec et ressort leur produit.
  Remarque importante: ce programme peut etre utilise pour tout ideal
  vp=[p,a,e,f,b] non necessairement premier, a condition que vp=pZ_K+aZ_K,
  que p soit entier et que la norme de vp soit p^f; e et b ne sont pas
  utilises.
*/

{
  long av=avma,tetpil,i,f,N;
  GEN m,x,y,dx,denx,p1;

  if((typ(nf)!=17)||(lg(nf)<10)) err(idealer1);
  N=lgef((GEN)nf[1])-3;
  if((typ(vp)!=17)||(lg(vp)!=6)) err(idealer3);
  if((typ(ix)<=10)||(typ(ix)==18)) ix=principalideal(nf,ix);
  if((typ(ix)==17)&&(lg(ix)==3)) {f=1;x=(GEN)ix[1];} else {f=0;x=ix;}
  if(typ(x)!=19) err(idealer2);
  if(N!=(lg((GEN)x[1])-1)) err(idealer4);
  denx=denom(x);if(!gcmp1(denx)) x=gmul(denx,x);
  if(lg(x)!=(N+1)) x=idealmul(nf,x,idmat(N));
  dx=gpui((GEN)vp[1],(GEN)vp[4],0);for(i=1;i<=N;i++) dx=mulii(dx,gcoeff(x,i,i));
  m=cgetg((N<<1)+1,19);
  for(i=1;i<=N;i++) m[i]=(long)element_muli(nf,(GEN)vp[2],(GEN)x[i]);
  for(i=N+1;i<=(N<<1);i++) m[i]=lmul((GEN)vp[1],(GEN)x[i-N]);
  tetpil=avma;p1=fasthnf(m,dx);
  if(gcmp1(denx))
  {
    if(f) {y=cgetg(3,17);y[1]=(long)p1;y[2]=lcopy((GEN)ix[2]);}
    else y=p1;
  }
  else
  {
    tetpil=avma;
    if(f) {y=cgetg(3,17);y[1]=ldiv(p1,denx);y[2]=lcopy((GEN)ix[2]);}
    else y=gdiv(p1,denx);
  }
  return gerepile(av,tetpil,y);
}

GEN
idealmulspec(GEN nf, GEN x, GEN spec)
/* a usage interne, aucune verif. ix est un ideal entier sans partie 
archimedienne en HNF, spec=[a,alpha,n] representant un ideal aZ_K+alpha Z_K
de norme n. Calcule le produit */
{
  long av=avma,tetpil,i,N=lgef((GEN)nf[1])-3;
  GEN dx,m;

  dx=(GEN)spec[3];for(i=1;i<=N;i++) dx=mulii(dx,gcoeff(x,i,i));
  m=cgetg((N<<1)+1,19);
  for(i=1;i<=N;i++) m[i]=(long)element_muli(nf,(GEN)spec[2],(GEN)x[i]);
  for(i=N+1;i<=(N<<1);i++) m[i]=lmul((GEN)spec[1],(GEN)x[i-N]);
  tetpil=avma;return gerepile(av,tetpil,fasthnf(m,dx));
}

long
idealvalint(GEN nf, GEN x, GEN vp)
                  
/* recoit un ideal entier x et un ideal premier vp dans le format
donne par primedec et calcule la valuation de ix en vp. A usage interne, pas de verifs.
*/

{
  long N,v,w,av=avma,i,j,bo;
  GEN mat,d,bp,p=(GEN)vp[1],p1,r;

  N=lg((GEN)x[1])-1;for(d=gun,i=1;i<=N;i++) d=mulii(d,(GEN)coeff(x,i,i));
  v=ggval(d,p);if(!v) return 0;
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
	  if(signe(r)) goto labelivint; else coeff(mat,i,j)=(long)p1;
	}
      w++;
    }
    else bo=1;
  }
  while((bo==0)&&(w<v));
  labelivint:
  avma=av;return w;
}


/****************************************/
/* Calcul de la differente d'un corps K */
/****************************************/

GEN
differente(GEN nf, GEN premiers)
            
/* Calcule la differente de nf */

{
  long av=avma,tetpil,i,j,vi,ei,v,nb_p,N,vpc,a;
  GEN ideal,mat_diff,liste_id,p1,p2,pcon,pr,pol;

  pol=(GEN)nf[1];N=lgef(pol)-3;
  if(DEBUGLEVEL) fprintferr("calcul de la differente\n");
  if(gcmp1((GEN)nf[4]))
  {
    p1=gmodulcp(deriv(pol,varn(pol)),pol);p2=idmat(N);
    tetpil=avma;return gerepile(av,tetpil,idealmul(nf,p1,p2));
  }
  ideal=gmul((GEN)nf[3],ginv((GEN)((GEN)nf[5])[4]));
  ideal=gdiv(ideal,pcon=content(ideal));
  if(DEBUGLEVEL) {fprintferr("temps D*delta^-1: ");fprintferr("%ld\n",timer2());flusherr();}
  ideal=hnfmodid(ideal,divii((GEN)nf[3],pcon));
  if(DEBUGLEVEL) {fprintferr("temps hnf(D*delta^-1): ");fprintferr("%ld\n",timer2());flusherr();}
  mat_diff=idmat(lgef((GEN)nf[1])-3);
  if(!premiers) premiers=factor(gabs((GEN)nf[3],0));
  nb_p=lg((GEN)premiers[1]);
  if(DEBUGLEVEL) {fprintferr("temps factor(D): ");fprintferr("%ld\n",timer2());flusherr();}
  for(i=1;i<nb_p;i++)
  {
    liste_id=primedec(nf,pr=gcoeff(premiers,i,1));vi=itos(gcoeff(premiers,i,2));
    vpc=ggval(pcon,pr);
    for(j=1;j<lg(liste_id);j++)
    {
      p1=(GEN)liste_id[j];ei=itos((GEN)p1[3]);
      if(ei>1)
      {
	if(DEBUGLEVEL) 
	{
	  fprintferr("traitement de ");outerr(p1);
	  flusherr();
	}
	if(signe(ressi(ei,pr))) v=ei-1;
	else
	{
	  v=ei*(vi-vpc)-idealval(nf,ideal,p1);
	}
	a=1+(v-1)/ei;p2=cgetg(4,17);p2[1]=(long)gpuigs(pr,a);
	p2[2]=(long)element_pow(nf,(GEN)p1[2],stoi(v));
	p2[3]=(long)gpui(pr,mulsi(v,(GEN)p1[4]),0);
	mat_diff=idealmulspec(nf,mat_diff,p2);
      }
    }
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(mat_diff));
}

GEN
idealnorm(GEN nf, GEN x)
/* Return the norm of an ideal, as a GEN */
{
  const long ltop = avma;
  long lbot,N;
  GEN x1, x2;
  const long tx = typ(x);

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;

  if (tx == 19)
  {
    x1 = x;
    if (lg(x) != N + 1) x1 = idealmul(nf, idmat(N), x);
                /* Perhaps from principal ideal */
  }
  else if (tx == 17 && lg(x) == 6)
    x1 = idealmulprime(nf, idmat(N), x);   /* Prime ideal */
  else err(idealer2);
  x2 = det(x1);
  lbot = avma;
  return gerepile(ltop, lbot, gabs(x2,0));
}

/****************************************/
/* Calcul de l'inverse d'un ideal       */
/****************************************/

GEN
idealinv(GEN nf, GEN ix)

/* Calcule le dual de mat_id pour la forme trace */

{
  long av=avma,tetpil,f,tx=typ(ix),N;
  GEN mat_dual,di,x,y,p1;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if((tx==18)&&(lg(ix)==N+1))
  {ix=gmodulcp(gmul((GEN)nf[7],ix),(GEN)nf[1]);tx=9;}
  if(tx==10) {ix=gmodulcp(ix,(GEN)nf[1]);tx=9;}
  if(tx<=9)
  {
    x=principalideal(nf,ginv(ix));y=idmat(N);
    tetpil=avma;return gerepile(av,tetpil,idealmul(nf,x,y));
  }
  if(tx==17)
  {
    if(lg(ix)==6)
    {
      p1=cgetg(6,17);p1[1]=ix[1];p1[2]=ix[5];
      p1[3]=p1[5]=zero;p1[4]=(long)subsi(N,(GEN)ix[4]);
      p1=idealmulprime(nf,idmat(N),p1);tetpil=avma;
      return gerepile(av,tetpil,gdiv(p1,(GEN)ix[1]));
    }
    if(lg(ix)!=3) err(idealer2);
    else {f=1;x=(GEN)ix[1];}
  }
  else {f=0;x=ix;}
  if(typ(x)!=19) err(idealer2);
  {					 /* THIS CODE NEW P.M & M.H. */
    GEN detx = det(x);
    GEN denx = denom(x);
    GEN xadjust = gcmp1(denx) ? detx : gmul(detx, gpuigs(denx, N));

    if (gcmp0(detx)) err(talker, "Cannot invert zero ideal");
    x = idealmul(nf, x, (GEN)((GEN)nf[5])[7]);
                         /* Multiply by fieldd * (inverse of different) */

    mat_dual = invmulmat(x, gmul(xadjust, (GEN)((GEN)nf[5])[6]));
    mat_dual = gdiv(gtrans(mat_dual), xadjust);

        /* nf[5][4] is a dense symmetric matrix.  We computed
           nf[5][6] = fieldd * ginv(nf[5][4]) in initalg.
           x is upper triangular (HNF), and easily inverted.
           The factor fieldd cancels while solving the linear equations.
        */

    di = denom(mat_dual);
    mat_dual = gmul(mat_dual, di);
    mat_dual = hnfmod(mat_dual, gdiv(gpuigs(di, N), detx));
  }
  tetpil=avma;
  if(f) {y=cgetg(3,17);y[1]=ldiv(mat_dual,di);y[2]=lneg((GEN)ix[2]);}
  else y=gdiv(mat_dual,di);
  return gerepile(av,tetpil,y);
} 

GEN
oldidealinv(GEN nf, GEN ix)
               
/* Calcule le dual de mat_id pour la forme trace */

{
  long av=avma,tetpil,f,tx=typ(ix),N;
  GEN mat_dual,di,x,y,p1;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if((tx==18)&&(lg(ix)==N+1))
  {ix=gmodulcp(gmul((GEN)nf[7],ix),(GEN)nf[1]);tx=9;}
  if(tx==10) {ix=gmodulcp(ix,(GEN)nf[1]);tx=9;}
  if(tx<=9)
  {
    x=principalideal(nf,ginv(ix));y=idmat(lgef((GEN)nf[1])-3);
    tetpil=avma;return gerepile(av,tetpil,idealmul(nf,x,y));
  }
  if(tx==17)
  {
    if(lg(ix)==6) 
    {
      p1=cgetg(6,17);p1[1]=ix[1];p1[2]=ix[5];N=lgef((GEN)nf[1])-3;
      p1[3]=p1[5]=zero;p1[4]=(long)subsi(N,(GEN)ix[4]);
      p1=idealmulprime(nf,idmat(N),p1);tetpil=avma;
      return gerepile(av,tetpil,gdiv(p1,(GEN)ix[1]));
    }
    if(lg(ix)!=3) err(idealer2);
    else {f=1;x=(GEN)ix[1];}
  }
  else {f=0;x=ix;}
  if(typ(x)!=19) err(idealer2);
  mat_dual=ginv(gmul(gtrans(x),(GEN)((GEN)nf[5])[4]));
  di=denom(mat_dual);mat_dual=gmul(di,mat_dual);
  mat_dual=idealmul(nf,(GEN)((GEN)nf[5])[5],mat_dual);
  tetpil=avma;
  if(f) {y=cgetg(3,17);y[1]=ldiv(mat_dual,di);y[2]=lneg((GEN)ix[2]);}
  else y=gdiv(mat_dual,di);
  return gerepile(av,tetpil,y);
}

/* Eleve un ideal premier vp a la puissance n ou n est dans Z */

GEN
idealpowprime(GEN nf, GEN vp, GEN n, long prec)
{
  long N,RU,av=avma,tetpil,s,i,m,ns;
  unsigned long j;
  GEN x,p1,p2,p3,vpinv;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if((typ(vp)!=17)||(lg(vp)!=6)) err(idealer3);
  p1=(GEN)nf[2];RU=itos((GEN)p1[1])+itos((GEN)p1[2]);
  x=cgetg(3,17);x[1]=(long)idmat(N);
  p1=cgetg(RU+1,17);
  for(i=1;i<=RU;i++)
  {
    p2=cgetg(3,6);p1[i]=(long)p2;
    p3=cgetr(prec);p2[1]=(long)p3;
    p3=cgetr(prec);p2[2]=(long)p3;
    affsr(0,(GEN)p2[1]);affsr(0,(GEN)p2[2]);
  }
  x[2]=(long)p1;
  s=signe(n);if(!s) return x;
  if(s<0) 
  {
    n=negi(n);vpinv=cgetg(6,17);vpinv[1]=vp[1];vpinv[2]=vp[5];
    vpinv[3]=vp[3];vpinv[4]=(long)subsi(N,(GEN)vp[4]);vpinv[5]=vp[2];
    vp=vpinv;
  }
  if(gcmpgs(n,16)<0)
  {
    ns=n[2];
    for(j=1;j<ns;j++) x=ideallllred(nf,idealmulprime(nf,x,vp),gzero,prec);
    tetpil=avma;x=ideallllred(nf,idealmulprime(nf,x,vp),gzero,prec);
  }
  else
  {
    m=n[lgef(n)-1];j=HIGHBIT;
    while((m&j)==0) 
      j>>=1;
    x=idealmulprime(nf,x,vp);j>>=1;
    if(gcmp1(n)) {tetpil=avma;x=ideallllred(nf,x,gzero,prec);}
    for(;j;j>>=1)
    {
      x=idealmulh(nf,x,x);
      if(m&j) x=idealmulprime(nf,x,vp);
      tetpil=avma;x=ideallllred(nf,x,gzero,prec);
    }
    for (i=lgef(n)-2;i>=2;i--)
    {
      for (m=n[i],j=HIGHBIT;j;j>>=1)
      {
	x=idealmulh(nf,x,x);
	if (m&j) x=idealmulprime(nf,x,vp);
	tetpil=avma;x=ideallllred(nf,x,gzero,prec);
      }
    }
  }
  return gerepile(av,tetpil,x);
}

/* Eleve un ideal ix a la puissance n ou n est dans Z */

GEN
idealpow(GEN nf, GEN ix, GEN n)
{
  long N,av=avma,tetpil,s,i,j,m,f;
  GEN iy,iz,y,x,denx,denz;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;iy=idmat(N);
  if((typ(ix)==18)&&(lg(ix)==N+1))
    ix=gmodulcp(gmul((GEN)nf[7],ix),(GEN)nf[1]);
  if(typ(ix)==10) ix=gmodulcp(ix,(GEN)nf[1]);
  if(typ(ix)<=9)
  {
    x=principalideal(nf,gpui(ix,n,0));tetpil=avma;
    return gerepile(av,tetpil,idealmul(nf,x,iy));
  }
  if((typ(ix)==17)&&(lg(ix)==6)) ix=idealmulprime(nf,iy,ix);
  if((typ(ix)==17)&&(lg(ix)==3)) {f=1;x=(GEN)ix[1];} else {f=0;x=ix;}
  if(typ(x)!=19) err(idealer2);
  if((N!=lg(x)-1)||(N!=lg((GEN)x[1])-1)) err(idealer4);
  s=signe(n);
  if(!s) 
  {
    if(f) {y=cgetg(3,17);y[1]=(long)iy;y[2]=lmul(gzero,(GEN)ix[2]);}
    else y=iy;
    return y;
  }
  if(s<0) n=negi(n);
  denx=denom(x);iz=gcmp1(denx)?x:gmul(x,denx);
  for (i=lgef(n)-1;i>2;i--)
  {
    for (m=n[i],j=0;j<BITS_IN_LONG;j++,m>>=1)
    {
      if (m&1) iy=idealmulh(nf,iy,iz);
      iz=idealmulh(nf,iz,iz);
    }
  }
  for (m=n[2];m>1;m>>=1)
  {
    if (m&1) iy=idealmulh(nf,iy,iz);
    iz=idealmulh(nf,iz,iz);
  }
  tetpil=avma;iy=idealmulh(nf,iy,iz);
  if(s<0) {n=negi(n);tetpil=avma;iy=idealinv(nf,iy);}
  if(!gcmp1(denx))
  {denz=gpui(denx,negi(n),0);tetpil=avma;iy=gmul(denz,iy);}
  if(f) {y=cgetg(3,17);y[1]=(long)iy;y[2]=lmul(n,(GEN)ix[2]);}
  else y=iy;
  return gerepile(av,tetpil,y);
}

GEN
idealpows(GEN nf, GEN ideal, long iexp)
    /* Return ideal**iexp in number field nf.  iexp is a C integer. */
{
  const long ltop = avma;
  GEN gexponent = stoi(iexp);
  long lbot = avma;
  return gerepile(ltop, lbot, idealpow(nf, ideal, gexponent));
} 

GEN
idealpowred(GEN nf, GEN ix, GEN n, long prec)
{
  long N,av=avma,tetpil,s,i,j,m,f;
  GEN iy,iz,y,x;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;iy=idmat(N);
  if((typ(ix)<=10)||(typ(ix)==18))
    ix=idealmul(nf,principalideal(nf,ix),iy);
  if((typ(ix)==17)&&(lg(ix)==6)) return idealpowprime(nf,ix,n,prec);
  if((typ(ix)==17)&&(lg(ix)==3)) {f=1;x=(GEN)ix[1];} else {f=0;x=ix;}
  if(typ(x)!=19) err(idealer2);
  if((N!=lg(x)-1)||(N!=lg((GEN)x[1])-1)) err(idealer4);
  s=signe(n);
  if(!s) 
  {
    if(f) {y=cgetg(3,17);y[1]=(long)iy;y[2]=lmul(gzero,(GEN)ix[2]);}
    else y=iy;
    return y;
  }
  if(s<0) n=negi(n);
  iz=x;
  for (i=lgef(n)-1;i>2;i--)
  {
    for (m=n[i],j=0;j<BITS_IN_LONG;j++,m>>=1)
    {
      if (m&1) iy=ideallllred(nf,idealmulh(nf,iy,iz),gzero,prec);
      iz=ideallllred(nf,idealmulh(nf,iz,iz),gzero,prec);
    }
  }
  for (m=n[2];m>1;m>>=1)
  {
    if (m&1) iy=ideallllred(nf,idealmulh(nf,iy,iz),gzero,prec);
    iz=ideallllred(nf,idealmulh(nf,iz,iz),gzero,prec);
  }
  iy=idealmulh(nf,iy,iz);if(s<0) {n=negi(n);iy=idealinv(nf,iy);}
  tetpil=avma;iy=ideallllred(nf,iy,gzero,prec);
  if(f) {y=cgetg(3,17);y[1]=(long)iy;y[2]=lmul(n,(GEN)ix[2]);}
  else y=iy;
  return gerepile(av,tetpil,y);
}
  
long
isideal(GEN nf,GEN x)  
{
  long N,av,f,i,j,k,tx=typ(x);
  GEN p1,minv,be;

  nf=checknf(nf);  
  if(tx<=9)
    return ((tx==1)||(tx==4)||(tx==5)||(tx==10)||((tx==9)&&gegal((GEN)nf[1],(GEN)x[1])));
  if((typ(x)==17)&&(lg(x)==3)) x=(GEN)x[1];
  if(typ(x)!=19) err(idealer2);
  N=lgef((GEN)nf[1])-3;if(lg((GEN)x[1])!=(N+1)) err(idealer4);
  av=avma;be=idmat(N);
  if(lg(x)!=(N+1)) x=idealmul(nf,x,idmat(N));
  x=gdiv(x,content(x));minv=ginv(x);f=1;
  for(i=1;(i<=N)&&f;i++)
    for(j=1;(j<=N)&&f;j++)
    {
      p1=gmul(minv,element_muli(nf,(GEN)x[i],(GEN)be[j]));
      for(k=1;(k<=N)&&f;k++) if(typ((GEN)p1[k])!=1) f=0;
    }
  avma=av;return f;
}

GEN
idealdiv(GEN nf, GEN x, GEN y)
{
  long av=avma,tetpil;
  GEN z;

  z=idealinv(nf,y);tetpil=avma;return gerepile(av,tetpil,idealmul(nf,x,z));
}

GEN
idealdivexact(GEN nf, GEN x, GEN y)
/*
	This routine computes the quotient x/y of two ideals
	in the number field nf.  It assumes that the quotient is
	an integral ideal.

	The idea is to find an ideal z dividing y
	such that gcd(N(x)/N(z), N(z)) = 1. Then

		x + (N(x)/N(z))    x
		--------------- = -----
		y + (N(y)/N(z))    y

	When x and y are integral ideals, this
	identity can be checked by looking at the exponent of
	a prime ideal p on both sides of the equation.

	Specifically, if a prime ideal p divides N(z),
	then it divides neither N(x)/N(z) nor N(y)/N(z)
	(since N(x)/N(z) is the product of the integers N(x/y) and N(y/z)).
	Both the numerator and the denominator on the left will
	be coprime to p.  So will x/y, since x/y is assumed integral
	and its norm N(x/y) is coprime to p

	If instead p does not divide N(z), then the power of p
	dividing N(x)/N(z) is the same as the power of p
	dividing N(x), which is at least as large as the
	power of p dividing x.  Likewise for N(y)/N(z).
	So the power of p dividing the left side equals the
	power of dividing the right side.

			Peter Montgomery

			July, 1994.

*/
{
  const long ltop = avma;
  long lbot,N;
  GEN cy = content(y);
  GEN x1, y1, detx1, dety1, detq, gcancel, gtemp;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;

  if (gcmp0(cy)) err(talker, "Cannot invert zero ideal");

  x1 = gdiv(x, cy);
  y1 = gdiv(y, cy);
  detx1 = idealnorm(nf, x1);
  dety1 = idealnorm(nf, y1);
  if (gcmp0(detx1))
  {
    avma = ltop;return gcopy(x);	/* If numerator is zero ideal */
  }
  detq = gdiv(detx1, dety1);
  if (!gcmp1(denom(x1)) || typ(detq) != 1)
    err(talker, "Quotient not integral in idealdivexact");
/*
	Find a norm gcancel such that
		(1) gcancel divides dety1;
		(2) gcd(detx1/gcancel, gcancel) = 1.
*/
  gcancel = dety1;
  do
  {
    gtemp = ggcd(gcancel, gdiv(detx1, gcancel));
    gcancel = gdiv(gcancel, gtemp);
  }
  while (!gcmp1(gtemp));

/*
		Replace x1/y1 by the ideal quotient

			x1 + (detx1/gcancel)
                        --------------------
                        y1 + (dety1/gcancel)
*/

  x1 = idealadd(nf, x1, gscalmat(gdiv(detx1, gcancel), N));

  if (gegal(gcancel, dety1)) /* If reducing y1 to unit ideal */
  { 
    gtemp = idmat(N);
			/* Would like to return x1, but idealadd
			   does not generate reduced HNF if original x1
			   was principal rather than HNF */
    lbot = avma;
    return gerepile(ltop, lbot, idealmul(nf, x1, gtemp));
  }
  else
  {
    y1 = idealadd(nf, y1, gscalmat(gdiv(dety1, gcancel), N));
    lbot = avma;
    return gerepile(ltop, lbot, idealdiv(nf, x1, y1));
  }
}

GEN
idealintersect(GEN nf, GEN x, GEN y)
{
  long av=avma,tetpil,lz,i,j,N;
  GEN z,p1,p2,dx,dy,dz;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if((typ(x)<=18)||((typ(x)==19)&&(lg(x)!=N+1))) x=idealhermite(nf,x);
  if((typ(y)<=18)||((typ(y)==19)&&(lg(y)!=N+1))) y=idealhermite(nf,y);
  dx=denom(x);dy=denom(y);dz=gmul(dx,dy);
  if(!gcmp1(dx)) y=gmul(y,dx);if(!gcmp1(dy)) x=gmul(x,dy);
  z=kerint(concat(x,y));lz=lg(z);p1=cgetg(lz,19);
  for(j=1;j<lz;j++)
  {
    p2=cgetg(N+1,18);p1[j]=(long)p2;
    for(i=1;i<=N;i++) p2[i]=coeff(z,i,j);
  }
  p2=gmul(x,p1);
  if(gcmp1(dz)) {tetpil=avma;return gerepile(av,tetpil,hnf(p2));}
  else {p2=hnf(p2);tetpil=avma;return gerepile(av,tetpil,gdiv(p2,dz));}
}

GEN
principalideal(GEN nf, GEN a)  
{
  long av,tetpil,N,i;
  GEN y,z;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  switch(typ(a))
  {
    case 1: case 4: case 5: z=cgetg(2,19);y=cgetg(N+1,18);z[1]=(long)y;
      y[1]=lcopy(a);for(i=2;i<=N;i++) y[i]=zero;
      return z;
    case 9: if(!gegal((GEN)nf[1],(GEN)a[1])) 
      err(talker,"incompatible number fields in principalideal");
      a=(GEN)a[2];   /* fall through */
    case 10: av=avma;z=cgetg(N+1,18);
      for(i=1;i<=N;i++) z[i]=(long)truecoeff(a,i-1);
      tetpil=avma;y=gmul((GEN)nf[8],z);z=cgetg(2,19);z[1]=(long)y;
      return gerepile(av,tetpil,z);
    case 18:
      if(lg(a)==N+1) {z=cgetg(2,19);z[1]=lcopy(a);return z;}
    default: err(talker,"incorrect type in principalideal");
      return gnil;
  }
}

GEN
principalidele(GEN nf, GEN a)  
{
  long av,tetpil,N,RU,R1,R2,i;
  GEN y,z,res,arc,p1;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  y=(GEN)nf[2];R1=itos((GEN)y[1]);R2=itos((GEN)y[2]);
  RU=R1+R2;
  switch(typ(a))
  {
    case 1: case 4: case 5: z=cgetg(2,19);y=cgetg(N+1,18);z[1]=(long)y;
      y[1]=lcopy(a);for(i=2;i<=N;i++) y[i]=zero;
      arc=cgetg(RU+1,18);for(i=1;i<=RU;i++) arc[i]=zero;
      res=cgetg(3,17);res[1]=(long)z;res[2]=(long)arc;return res;
    case 9: if(!gegal((GEN)nf[1],(GEN)a[1])) 
      err(talker,"incompatible number fields in principalidele");
      a=(GEN)a[2];   /* fall through */
    case 10: av=avma;z=cgetg(N+1,18);
      for(i=1;i<=N;i++) z[i]=(long)truecoeff(a,i-1);
      y=gmul((GEN)nf[8],z);p1=gmul((GEN)((GEN)nf[5])[1],y);arc=cgetg(RU+1,18);
      for(i=1;i<=R1;i++) arc[i]=(long)glog((GEN)p1[i],MEDDEFAULTPREC);
      for(i=R1+1;i<=RU;i++) arc[i]=lmul2n(glog((GEN)p1[i],MEDDEFAULTPREC),1);
      tetpil=avma;z=cgetg(2,19);z[1]=lcopy(y);
      res=cgetg(3,17);res[1]=(long)z;res[2]=lcopy(arc);
      return gerepile(av,tetpil,res);
    case 18:
      if(lg(a)==N+1)
      {
	av=avma;y=gmul((GEN)nf[7],a);tetpil=avma;
	return gerepile(av,tetpil,principalidele(nf,y));
      }
    default: err(talker,"incorrect type in principalidele");
      return gnil;
  }
}

GEN
ideallllredall(GEN nf, GEN ix, GEN vdir, long prec, long precint)
{
  long N,av=avma,tetpil,i,j,f,r1,r2,ru;
  GEN T,p1,p2,p3,y,alpha,beta,x,x2,v,z,detmat;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if((!gcmp0(vdir))&&(typ(vdir)!=17)) err(idealer5);
  if((typ(ix)<=10)||(typ(ix)==18)) ix=principalidele(nf,ix);
  if((typ(ix)==17)&&(lg(ix)==6)) {f=0;ix=idealmulprime(nf,idmat(N),ix);}
  if((typ(ix)==17)&&(lg(ix)==3)) {f=1;x=(GEN)ix[1];} else {f=0;x=ix;}
  if(typ(x)!=19) err(idealer2);
  if(DEBUGLEVEL>=6)
    fprintferr("entree dans ideallllred, time = %ld\n",timer2());flusherr();
  if(!gcmp1(gcoeff(x,N,N))) 
  {
    p1=content(x);if(!gcmp1(p1)) x=gdiv(x,p1);
    if(DEBUGLEVEL>=6)
      fprintferr("temps content : %ld\n",timer2());flusherr();
  }
  p1=(GEN)nf[2];r2=itos((GEN)p1[2]);ru=N-r2;r1=ru-r2;p1=(GEN)nf[5];
  if(!gcmp0(vdir))
  {
    if(lg(vdir)!=(ru+1)) err(idealer5);
    p3=(GEN)p1[2];p2=cgetg(ru+1,19);
    for(j=1;j<=ru;j++) 
    {
      if(!gcmp0((GEN)vdir[j]))
      {
	if(typ((GEN)vdir[j])==1) 
	  p2[j]=lmul2n((GEN)p3[j],itos((GEN)vdir[j])<<1);
	else
	  p2[j]=lmul((GEN)p3[j],gpui(stoi(4),(GEN)vdir[j],0));
      }
      else p2[j]=p3[j];
    }
    p1=greal(gmul(p2,(GEN)p1[1]));
  }
  else p1=(GEN)p1[3];
  if(DEBUGLEVEL>=6)
    fprintferr("temps initialisations : %ld\n",timer2());flusherr();
  y=gmul(x,(GEN)(lllgram(gmul(gtrans(x),gmul(p1,x)),2*precint-2)[1]));
  if(DEBUGLEVEL>=6)
    fprintferr("temps lllgram : %ld\n",timer2());flusherr();
  for(i=2;(i<=N)&&gcmp0((GEN)y[i]);i++);
  if(i>N) 
  {
    tetpil=avma;
    if(f) {y=cgetg(3,17);y[1]=lcopy(x);y[2]=lcopy((GEN)ix[2]);}
    else y=gcopy(x);
    return gerepile(av,tetpil,y);
  }
  T=(GEN)nf[1];alpha=gmodulcp(gmul((GEN)nf[7],y),T);
  beta=lift(gdiv(gnorm(alpha),alpha));
  if(DEBUGLEVEL>=6)
    fprintferr("temps alpha/beta : %ld\n",timer2());flusherr();
  z=gmul((GEN)((GEN)(nf[5]))[1],y);
  p1=cgetg(N+1,18);for(i=1;i<=N;i++) p1[i]=(long)truecoeff(beta,i-1);
  p1=gmul((GEN)nf[8],p1);p2=cgetg(N+1,19);
  for(i=1;i<=N;i++) p2[i]=(long)element_muli(nf,p1,(GEN)x[i]);
  p1=content(p2);if(!gcmp1(p1)) p2=gdiv(p2,p1);
  if(DEBUGLEVEL>=6)
    fprintferr("temps nouvel ideal : %ld\n",timer2());flusherr();
  detmat=detint(p2);
  if(f)
  {
    x2=(GEN)ix[2];v=cgetg(ru+1,17);
    for(i=1;i<=r1;i++) v[i]=(long)glog((GEN)z[i],prec);
    for(i=r1+1;i<=ru;i++) v[i]=lmul2n(glog((GEN)z[i],prec),1);
    tetpil=avma;y=cgetg(3,17);y[1]=(long)hnfmod(p2,detmat);p1=cgetg(ru+1,17);
    y[2]=(long)p1;for(i=1;i<=ru;i++) p1[i]=lsub((GEN)x2[i],(GEN)v[i]);
  }
  else {tetpil=avma;y=hnfmod(p2,detmat);}
  if(DEBUGLEVEL>=6)
    fprintferr("temps hermite final : %ld\n",timer2());flusherr();
  return gerepile(av,tetpil,y);
}

GEN
ideallllred(GEN nf, GEN ix, GEN vdir, long prec)
{
  return ideallllredall(nf,ix,vdir,prec,prec);
}

GEN
ideallllredpart1(GEN nf, GEN x, GEN vdir, long flprem, long prec)
/* a usage interne, pas de verifs et pas de gestion de pile */
{
  long N=lgef((GEN)nf[1])-3,i,j,r2,ru,N2,av,tetpil;
  GEN p1,p2,p3,p4,y,alpha,z,xpro;

  av=avma;
  if(!gcmp1(gcoeff(x,N,N))) {p1=content(x);if(!gcmp1(p1)) x=gdiv(x,p1);}
  p1=(GEN)nf[2];r2=itos((GEN)p1[2]);ru=N-r2;p1=(GEN)nf[5];
  if(!gcmp0(vdir))
  {
    p3=(GEN)p1[2];p2=cgetg(ru+1,19);
    for(j=1;j<=ru;j++) 
    {
      if(!gcmp0((GEN)vdir[j]))
      {
	if(typ((GEN)vdir[j])==1) 
	  p2[j]=lmul2n((GEN)p3[j],itos((GEN)vdir[j])<<1);
	else
	  p2[j]=lmul((GEN)p3[j],gpui(stoi(4),(GEN)vdir[j],0));
      }
      else p2[j]=p3[j];
    }
    p1=greal(gmul(p2,(GEN)p1[1]));
  }
  else p1=(GEN)p1[3];
  if(!flprem)
  {
	/* N2=max(2,(N+1)>>1); */ 
    N2=N;xpro=cgetg(N2+1,19);for(i=1;i<=N2;i++) xpro[i]=x[i];
  }
  else xpro=x;
  p4=lllgram(gmul(gtrans(xpro),gmul(p1,xpro)),prec+1);
  y=gmul(xpro,(GEN)p4[1]);
  for(i=2;(i<=N)&&gcmp0((GEN)y[i]);i++);
  if(i>N) y=gmul(xpro,(GEN)p4[2]);
  alpha=gmul((GEN)nf[7],y);p1=subres(alpha,(GEN)nf[1]);
  z=cgetg(4,17);z[1]=(long)x;z[2]=(long)y;z[3]=(long)gabs(p1,0);
  tetpil=avma;return gerepile(av,tetpil,gcopy(z));
}

GEN
ideallllredpart1spec(GEN nf, GEN x, GEN matt2, long flprem, long prec)
/* a usage interne, pas de verifs et pas de gestion de pile */
{
  long N=lgef((GEN)nf[1])-3,i,N2,av,tetpil;
  GEN p1,p4,y,alpha,z,xpro;

  av=avma;
  if(!gcmp1(gcoeff(x,N,N))) {p1=content(x);if(!gcmp1(p1)) x=gdiv(x,p1);}
  if(!flprem) 
  {
	/* N2=max(2,(N+1)>>1); */
    N2=N;xpro=cgetg(N2+1,19);for(i=1;i<=N2;i++) xpro[i]=x[i];
  }
  else xpro=x;
  p4=lllgram(gmul(gtrans(xpro),gmul(matt2,xpro)),prec+1);
  y=gmul(xpro,(GEN)p4[1]);
  for(i=2;(i<=N)&&gcmp0((GEN)y[i]);i++);
  if(i>N) y=gmul(xpro,(GEN)p4[2]);
  alpha=gmul((GEN)nf[7],y);p1=subres(alpha,(GEN)nf[1]);
  z=cgetg(4,17);z[1]=(long)x;z[2]=(long)(y);z[3]=(long)gabs(p1,0);
  tetpil=avma;return gerepile(av,tetpil,gcopy(z));
}

GEN
ideallllredpart2(GEN nf, GEN arch, GEN gamma, long prec)
{
  long N=lgef((GEN)nf[1])-3,av=avma,tetpil,i,r1,r2,ru;
  GEN p1,v,z;
  
  z=gmul((GEN)((GEN)(nf[5]))[1],gamma);
  p1=(GEN)nf[2];r2=itos((GEN)p1[2]);ru=N-r2;r1=ru-r2;
  v=cgetg(ru+1,17);
  for(i=1;i<=r1;i++) v[i]=(long)glog((GEN)z[i],prec);
  for(i=r1+1;i<=ru;i++) v[i]=lmul2n(glog((GEN)z[i],prec),1);
  tetpil=avma;p1=cgetg(ru+1,17);
  for(i=1;i<=ru;i++) p1[i]=lsub((GEN)arch[i],(GEN)v[i]);
  return gerepile(av,tetpil,p1);
}

GEN
minideal(GEN nf, GEN ix, GEN vdir, long prec)
{
  long N,av=avma,tetpil,i,j,f,r1,r2,ru;
  GEN p1,p2,p3,y,x,v,z;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if((!gcmp0(vdir))&&(typ(vdir)!=17)) err(idealer5);
  if((typ(ix)==17)&&(lg(ix)==3)) {f=1;x=(GEN)ix[1];} else {f=0;x=ix;}
  if(typ(x)!=19) err(idealer2);
  p1=(GEN)nf[2];r2=itos((GEN)p1[2]);ru=N-r2;r1=ru-r2;p1=(GEN)nf[5];
  if(!gcmp0(vdir))
  {
    p3=(GEN)p1[2];p2=cgetg(ru+1,19);
    for(j=1;j<=ru;j++) p2[j]=lmul2n((GEN)p3[j],itos((GEN)vdir[j])<<1);
    p1=greal(gmul(p2,(GEN)p1[1]));
  }
  else p1=(GEN)p1[3];
  y=gmul(x,(GEN)(lllgram(gmul(gtrans(x),gmul(p1,x)),prec)[1]));
  z=gmul((GEN)((GEN)(nf[5]))[1],y);
  tetpil=avma;p2=cgetg(3,17);p2[1]=(long)gtomat(y);v=cgetg(ru+1,17);
  for(i=1;i<=r1;i++) v[i]=(long)glog((GEN)z[i],prec);
  for(i=r1+1;i<=ru;i++) v[i]=lmul2n(glog((GEN)z[i],prec),1);
  p2[2]=(long)v;
  return gerepile(av,tetpil,p2);
}

GEN
idealapprall(GEN nf, GEN x, long fl)
/* Given a fractional ideal x (if fl=0) or a prime ideal factorization
   with possibly zero or negative exponents (if fl=1), gives a b such that
   v_p(b)=v_p(x) for all prime ideals p dividing x (or in the ideal factorization)
   and v_p(b)>=0 for all other p, using the (standard) proof given in GTM 138.
   Certainly not the most efficient, but sure */
{
  long av=avma,tetpil,i,j,k,l,N,r,r2;
  GEN fact,fact2,list,ep,ep1,ep2,y,z,v,p1,p2,p3,p4,s,pr,alpha,beta,den,u1;

  if(DEBUGLEVEL>2)
  {
    fprintferr(" entree dans idealapprall() :\n");
    fprintferr(" x = ");outerr(x);flusherr();
  }
  if(fl)
  {
    nf=checknf(nf);N=lgef((GEN)nf[1])-3;
    if((typ(x)!=19)||(lg(x)!=3))
      err(talker,"not a prime ideal factorization in idealapprall");
    fact=x;list=(GEN)fact[1];ep=(GEN)fact[2];r=lg(list);
    if(r==1)
    {
      avma=av;z=cgetg(N+1,18);z[1]=un;for(i=2;i<=N;i++) z[i]=zero;
      return z;
    }
    else
    {
      for(i=1;(i<r)&&(signe((GEN)ep[i])>=0);i++);
      if(i<r)
      {
	ep1=cgetg(r,18);
	for(i=1;i<r;i++) ep1[i]=(long)gmax(negi((GEN)ep[i]),gzero);
	fact[2]=(long)ep1;beta=idealapprall(nf,fact,1);
	fact2=idealfactor(nf,gtomat(beta));
	p1=(GEN)fact2[1];r2=lg(p1);ep2=(GEN)fact2[2];
	l=r+r2-1;z=cgetg(l,18);for(i=1;i<r;i++) z[i]=list[i];
	ep1=cgetg(l,18);
	for(i=1;i<r;i++) ep1[i]=(long)gmax((GEN)ep[i],gzero);
	j=r-1;
	for(i=1;i<r2;i++)
	{
	  p3=(GEN)p1[i];
	  for(k=1;(k<r)&&(!vecegal((GEN)list[k],p3));k++);
	  if(k==r){j++;z[j]=(long)p3;ep1[j]=ep2[i];}
	}
	p2=cgetg(j+1,17);for(i=1;i<=j;i++) p2[i]=z[i];
	p3=cgetg(j+1,17);for(i=1;i<=j;i++) p3[i]=ep1[i];
	fact=cgetg(3,19);fact[1]=(long)p2;fact[2]=(long)p3;
	alpha=idealapprall(nf,fact,1);
	if(DEBUGLEVEL>2)
	{
	  fprintferr(" alpha = ");outerr(alpha);
	  fprintferr(" beta = ");outerr(beta);flusherr();
	}
	tetpil=avma;return gerepile(av,tetpil,element_div(nf,alpha,beta));
      }	
      y=idmat(N);
      for(i=1;i<r;i++)
      {
	pr=(GEN)list[i];
	if(signe((GEN)ep[i]))
	{
	  if(cmpis((GEN)pr[4],N))
	  {
	    p1=element_pow(nf,(GEN)pr[2],p4=addii(gun,(GEN)ep[i]));
	    p2=cgetg(3,19);p3=cgetg(N+1,18);p2[1]=(long)p3;
	    p3[1]=(long)gpui((GEN)pr[1],p4,0);for(j=2;j<=N;j++) p3[j]=zero;
	    p2[2]=(long)p1;
	    y=idealmul(nf,y,p2);
	  }
	  else y=gmul(gpui((GEN)pr[1],addii(gun,(GEN)ep[i]),0),y);
	}
	else y=idealmulprime(nf,y,pr);
      }
    }
  }
  else
  {
    x=idealhermite(nf,x);N=lgef((GEN)nf[1])-3;
    den=denom(x);if(!gcmp1(den)) x=gmul(den,x);
    fact=idealfactor(nf,x);
    list=(GEN)fact[1];ep=(GEN)fact[2];r=lg(list);
    if(!gcmp1(den))
    {
      fact2=idealfactor(nf,den);
      p1=(GEN)fact2[1];r2=lg(p1);
      l=r+r2-1;z=cgetg(l,18);for(i=1;i<r;i++) z[i]=list[i];
      ep1=cgetg(l,18);
      for(i=1;i<r;i++) ep1[i]=ep[i];
      j=r-1;
      for(i=1;i<r2;i++)
      {
	p3=(GEN)p1[i];
	for(k=1;(k<r)&&(!vecegal((GEN)list[k],p3));k++);
	if(k==r){j++;z[j]=(long)p3;ep1[j]=zero;}
      }
      p2=cgetg(j+1,17);for(i=1;i<=j;i++) p2[i]=z[i];
      p3=cgetg(j+1,17);for(i=1;i<=j;i++) p3[i]=ep1[i];
      fact=cgetg(3,19);fact[1]=(long)p2;fact[2]=(long)p3;
      alpha=idealapprall(nf,fact,1);
      if(DEBUGLEVEL>2){fprintferr(" alpha = ");outerr(alpha);flusherr();}
      tetpil=avma;return gerepile(av,tetpil,gdiv(alpha,den));
    }
    y=x;for(i=1;i<r;i++) y=idealmul(nf,y,(GEN)list[i]);
  }
  z=cgetg(r,17);
  for(i=1;i<r;i++)
  {
    pr=(GEN)list[i];
    if(cmpis((GEN)pr[4],N))
    {
      p1=element_pow(nf,(GEN)pr[5],p4=addii(gun,(GEN)ep[i]));
      p2=cgetg(3,19);p3=cgetg(N+1,18);p2[1]=(long)p3;
      p3[1]=(long)gpui((GEN)pr[1],p4,0);for(j=2;j<=N;j++) p3[j]=zero;
      p2[2]=(long)p1;
      z[i]=ldiv(idealmul(nf,y,p2),(GEN)p3[1]);
    }
    else z[i]=(long)gdiv(y,gpui((GEN)pr[1],addii(gun,(GEN)ep[i]),0));
  }
  v=idealaddmultone(nf,z);
  s=cgetg(N+1,18);for(i=1;i<=N;i++) s[i]=zero;
  for(i=1;i<r;i++)
  {
    pr=(GEN)list[i];
    if(signe((GEN)ep[i]))
      s=gadd(s,element_mul(nf,(GEN)v[i],element_pow(nf,(GEN)pr[2],(GEN)ep[i])));
    else s=gadd(s,(GEN)v[i]);
  }
  s=gmod(s,gcoeff(y,1,1));
  y=gmul(y,lllint(y));
  z=cgetg(N+2,19);for(i=1;i<=N;i++) z[i]=y[i];z[N+1]=(long)s;
  u1=(GEN)ker(z)[1];u1=gmul(u1,denom(u1));p2=(GEN)u1[N+1];
  p1=cgetg(N+1,18);for(i=1;i<=N;i++) p1[i]=lround(gdiv((GEN)u1[i],p2));
  p3=gmul(y,p1);p3=gadd(s,p3);
  if(DEBUGLEVEL>2)
  {fprintferr(" sortie de idealapprall p3 = ");outerr(p3);flusherr();}
  tetpil=avma;
  return gerepile(av,tetpil,gcopy(p3));
}

GEN
idealchinese(GEN nf, GEN x, GEN y)
/* Given a prime ideal factorization x with possibly zero or negative exponents,
   and a vector y of elements of nf, gives a b such that
   v_p(b-y_p)>=v_p(x) for all prime ideals p in the ideal factorization
   and v_p(b)>=0 for all other p, using the (standard) proof given in GTM 138.
   Certainly not the most efficient, but sure */
{
  long av=avma,tetpil,i,j,k,l,N,r,r2;
  GEN fact,fact2,list,ep,ep1,ep2,z,t,v,p1,p2,p3,p4,s,pr,den,u1;

  if(DEBUGLEVEL>2)
  {
    fprintferr(" entree dans idealchinese() :\n");
    fprintferr(" x = ");outerr(x);
    fprintferr(" y = ");outerr(y);flusherr();
  }
  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if((typ(x)!=19)||(lg(x)!=3))
    err(talker,"not a prime ideal factorization in idealchinese");
  fact=x;list=(GEN)fact[1];ep=(GEN)fact[2];r=lg(list);
  if((typ(y)<17)||(typ(y)>18)||(lg(y)!=r))
    err(talker,"not a suitable vector of elements in idealchinese");
  if(r==1){avma=av;z=cgetg(N+1,18);z[1]=un;for(i=2;i<=N;i++) z[i]=zero;return z;}
  den=denom(y);
  if(!gcmp1(den))
  {
    fact2=idealfactor(nf,den);
    p1=(GEN)fact2[1];r2=lg(p1);ep2=(GEN)fact2[2];
    l=r+r2-1;z=cgetg(l,18);for(i=1;i<r;i++) z[i]=list[i];
    ep1=cgetg(l,18);
    for(i=1;i<r;i++) ep1[i]=ep[i];
    j=r-1;
    for(i=1;i<r2;i++)
    {
      p3=(GEN)p1[i];
      for(k=1;(k<r)&&(!vecegal((GEN)list[k],p3));k++);
      if(k==r){j++;z[j]=(long)p3;ep1[j]=ep2[i];}
      else ep1[k]=(long)gadd((GEN)ep1[k],(GEN)ep2[i]);
    }
    r=j+1;
    list=cgetg(r,17);for(i=1;i<r;i++) list[i]=z[i];
    ep=cgetg(r,17);for(i=1;i<r;i++) ep[i]=ep1[i];
  }
  for(i=1;i<r;i++) ep[i]=(long)gmax((GEN)ep[i],gzero);
  t=idmat(N);
  for(i=1;i<r;i++)
  {
    pr=(GEN)list[i];
    if(signe((GEN)ep[i]))
    {
      if(cmpis((GEN)pr[4],N))
      {
	p1=element_pow(nf,(GEN)pr[2],p4=(GEN)ep[i]);
	p2=cgetg(3,19);p3=cgetg(N+1,18);p2[1]=(long)p3;
	p3[1]=(long)gpui((GEN)pr[1],p4,0);for(j=2;j<=N;j++) p3[j]=zero;
	p2[2]=(long)p1;
	t=idealmul(nf,t,p2);
      }
      else t=gmul(gpui((GEN)pr[1],(GEN)ep[i],0),t);
    }
  }
  z=cgetg(r,17);
  for(i=1;i<r;i++)
  {
    pr=(GEN)list[i];
    if(cmpis((GEN)pr[4],N))
    {
      p1=element_pow(nf,(GEN)pr[5],p4=(GEN)ep[i]);
      p2=cgetg(3,19);p3=cgetg(N+1,18);p2[1]=(long)p3;
      p3[1]=(long)gpui((GEN)pr[1],p4,0);for(j=2;j<=N;j++) p3[j]=zero;
      p2[2]=(long)p1;
      z[i]=ldiv(idealmul(nf,t,p2),(GEN)p3[1]);
    }
    else z[i]=ldiv(t,gpui((GEN)pr[1],(GEN)ep[i],0));
  }
  v=idealaddmultone(nf,z);
  s=cgetg(N+1,18);for(i=1;i<=N;i++) s[i]=zero;
  for(i=1;i<r;i++)
  {
    pr=(GEN)list[i];
    s=gadd(s,element_mul(nf,(GEN)v[i],(GEN)y[i]));
  }
  s=gmod(s,gcoeff(t,1,1));
  t=gmul(t,lllint(t));
  z=cgetg(N+2,19);for(i=1;i<=N;i++) z[i]=t[i];z[N+1]=(long)s;
  u1=(GEN)ker(z)[1];u1=gmul(u1,denom(u1));p2=(GEN)u1[N+1];
  p1=cgetg(N+1,18);for(i=1;i<=N;i++) p1[i]=lround(gdiv((GEN)u1[i],p2));
  p3=gmul(t,p1);tetpil=avma;p3=gadd(s,p3);
  if(DEBUGLEVEL>2)
  {fprintferr(" sortie de idealchinese() : p3 = ");outerr(p3);flusherr();}
  return gerepile(av,tetpil,p3);
}

GEN
idealappr(GEN nf, GEN x)
{
  return idealapprall(nf,x,0);
}

GEN

idealapprfact(GEN nf, GEN x)
{
  return idealapprall(nf,x,1);
}

GEN
ideal_two_elt2(GEN nf, GEN x, GEN a)
/* Given an integral ideal x and a in x, gives a b such that
   x=aZ_K+bZ_K using a different algorihm than ideal_two_elt */
{
  long av=avma,tetpil,i,N,r;
  GEN p1,p2,b,li,list,z;

  x=idealhermite(nf,x);N=lgef((GEN)nf[1])-3;
  if((typ(a)<=10)||(typ(a)==18)) a=principalideal(nf,a);
  if((typ(a)!=19)||(lg(a)!=2))
    err(talker,"incorrect type in ideal_two_elt2");
  if(gcmp0(x))
  {
    if(!gcmp0(a)) err(talker,"element not in ideal in ideal_tow_elt2");
    tetpil=avma;return gerepile(av,tetpil,gcopy(a));
  }
  if(!gcmp1(denom(gauss(x,a))))
    err(talker,"element does not belong to ideal in ideal_two_elt2");
  p1=idealfactor(nf,a);list=(GEN)p1[1];r=lg(list);
  z=cgetg(3,19);p2=cgetg(r,18);z[1]=(long)list;z[2]=(long)p2;
  for(i=1;i<r;i++) p2[i]=lstoi(idealval(nf,x,(GEN)list[i]));
  li=gcoeff(x,1,1);b=idealapprall(nf,z,1);
  b=gmul(gmodulcp(gun,li),b);tetpil=avma;
  return gerepile(av,tetpil,centerlift(b));
}

GEN
idealcoprime(GEN nf, GEN x, GEN y)
/* Given 2 integral ideals x and y in a number field nf gives a beta belonging to nf
   such that beta.x is an integral ideal coprime to y */
{
  long av=avma,tetpil,i,r;
  GEN fact,p1,p2,ep;

  if(DEBUGLEVEL>2)
  {
    fprintferr(" entree dans idealcoprime() :\n");
    fprintferr(" x = ");outerr(x);
    fprintferr(" y = ");outerr(y);flusherr();
  }
  fact=idealfactor(nf,y);p1=(GEN)fact[1];r=lg(p1);
  ep=cgetg(r,18);
  for(i=1;i<r;i++) ep[i]=lstoi(-idealval(nf,x,(GEN)p1[i]));
  fact[2]=(long)ep;tetpil=avma;p2=idealapprall(nf,fact,1);
  if(DEBUGLEVEL>2)
  {fprintferr(" sortie de idealcoprime() : p2 = ");outerr(p2);flusherr();}
  return gerepile(av,tetpil,p2);
}

GEN
idealcoprimeinv(GEN nf, GEN x, GEN y)
/* Given 2 integral ideals x and y in a number field nf gives a beta belonging to nf
   such that beta.x^(-1) is an integral ideal coprime to y */
{
  long av=avma,tetpil,i,r;
  GEN fact,p1,p2,ep;

  if(DEBUGLEVEL>2)
  {
    fprintferr(" entree dans idealcoprimeinv() :\n");
    fprintferr(" x = ");outerr(x);
    fprintferr(" y = ");outerr(y);flusherr();
  }
  fact=idealfactor(nf,idealmul(nf,x,y));p1=(GEN)fact[1];r=lg(p1);
  ep=cgetg(r,18);
  for(i=1;i<r;i++) ep[i]=lstoi(idealval(nf,x,(GEN)p1[i]));
  fact[2]=(long)ep;tetpil=avma;p2=idealapprall(nf,fact,1);
  if(DEBUGLEVEL>2)
  {fprintferr(" sortie de idealcoprimeinv() : p2 = ");outerr(p2);flusherr();}
  return gerepile(av,tetpil,p2);
}

long
isinvector(GEN v, GEN x, long n)
/* returns the first index i<=n such that x=v[i] if it exits, 0 otherwise */
{
  long i;

  for(i=1;i<=n;i++) if(gegal((GEN)v[i],x)) return i;
  return 0;
}
  
GEN
idealcoprimeinvabc(GEN nf, GEN x, GEN a, GEN b, GEN c)
/* Given an integral ideal x and three algebraic integers a, b and c in a number
   field nf gives a beta belonging to nf such that beta.x^(-1) is an integral
   ideal coprime to abc.Z_k */
{
  long av=avma,tetpil,i,j,r,ra,rb,rc;
  GEN facta,factb,factc,fact,factall,p1,p2,ep;

  if(DEBUGLEVEL>2)
  {
    fprintferr(" entree dans idealcoprimeinvabc() :\n");
    fprintferr(" x = ");outerr(x);fprintferr(" a = ");outerr(a);
    fprintferr(" b = ");outerr(b);fprintferr(" c = ");outerr(c);
    flusherr();
  }
  facta=(GEN)idealfactor(nf,a)[1];factb=(GEN)idealfactor(nf,b)[1];
  factc=(GEN)idealfactor(nf,c)[1];ra=lg(facta);rb=lg(factb);rc=lg(factc);
  factall=cgetg(ra+rb+rc-2,18);
  for(i=1;i<ra;i++) factall[i]=facta[i];j=ra-1;
  for(i=1;i<rb;i++) if(!isinvector(factall,(GEN)factb[i],j)) factall[++j]=factb[i];
  for(i=1;i<rc;i++) if(!isinvector(factall,(GEN)factc[i],j)) factall[++j]=factc[i];
  r=j+1;fact=cgetg(3,19);p1=cgetg(r,18);ep=cgetg(r,18);
  for(i=1;i<r;i++) p1[i]=factall[i];
  for(i=1;i<r;i++) ep[i]=lstoi(idealval(nf,x,(GEN)p1[i]));
  fact[1]=(long)p1;fact[2]=(long)ep;tetpil=avma;p2=idealapprall(nf,fact,1);
  if(DEBUGLEVEL>2)
  {
    fprintferr(" sortie de idealcoprimeinvabc() : p2 = ");outerr(p2);
    flusherr();
  }
  return gerepile(av,tetpil,p2);
}


  
GEN
newidealcoprimeinvabc(GEN nf, GEN x, GEN a, GEN b, GEN c)
/* Given an integral ideal x and three algebraic integers a, b and c in a number
   field nf gives a beta belonging to nf such that beta.x^(-1) is an integral
   ideal coprime to abc.Z_k */
{
  long av=avma,tetpil,i,j,r,r2,ra,rb,rc,N=lgef((GEN)nf[1])-3;
  GEN facta,factb,factc,fact,factall,list,pr,p1,p2,p3,p4,abc,ep;

  fact=idealfactor(nf,x);list=(GEN)fact[1];ep=(GEN)fact[2];r=lg(list);
  p1=idealapprall(nf,fact,1);
  facta=(GEN)idealfactor(nf,a)[1];factb=(GEN)idealfactor(nf,b)[1];
  factc=(GEN)idealfactor(nf,c)[1];ra=lg(facta);rb=lg(factb);rc=lg(factc);
  factall=cgetg(ra+rb+rc-2,18);
  for(i=1;i<ra;i++) factall[i]=facta[i];j=ra-1;
  for(i=1;i<rb;i++) if(!isinvector(factall,(GEN)factb[i],j)) factall[++j]=factb[i];
  for(i=1;i<rc;i++) if(!isinvector(factall,(GEN)factc[i],j)) factall[++j]=factc[i];
  r2=j+1;
  p2=x;p3=idmat(N);
  for(i=1;i<r2;i++)
  {
    pr=(GEN)factall[i];
    if(isinvector(list,pr,r-1)) p2=idealmulprime(nf,p2,pr);
    else p3=idealmulprime(nf,p3,pr);
  }
  p4=idealaddone(nf,p2,p3);
  p1=gadd(element_mul(nf,p1,(GEN)p4[2]),(GEN)p4[1]);
  abc=element_mul(nf,element_mul(nf,a,b),c);
  tetpil=avma;return gerepile(av,tetpil,nfmod(nf,p1,abc));
}

GEN
findX(GEN nf, GEN a, GEN b, GEN J, GEN M)
/* A usage interne (pas de verification ni de gestion pile) :
   Solve the equation ((b+aX).Z_k/((a,b).J),M)=Z_k. */
{
  long N,i,k,r,v;
  GEN p1,p2,abJ,fact,list,ve,ep,int0,int1,int2,pr;
    
  if(DEBUGLEVEL>2)
  {
    fprintferr(" entree dans findX() :\n");
    fprintferr(" a = ");outerr(a);fprintferr(" b = ");outerr(b);
    fprintferr(" J = ");outerr(J);fprintferr(" M = ");outerr(M);
    flusherr();
  }
  N=lgef((GEN)nf[1])-3;
  p1=cgetg(3,19);p1[1]=(long)a;p1[2]=(long)b;
  if(N==2) p1=idealmul(nf,p1,idmat(2));
  abJ=idealmul(nf,p1,J);
  fact=idealfactor(nf,M);list=(GEN)fact[1];r=lg(list);
  ve=cgetg(r,17);ep=cgetg(r,17);
  int0=cgetg(N+1,18);for(i=1;i<=N;i++) int0[i]=zero;
  int1=cgetg(N+1,18);int1[1]=un;for(i=2;i<=N;i++) int1[i]=zero;
  int2=cgetg(N+1,18);int2[1]=deux;for(i=2;i<=N;i++) int2[i]=zero;
  for(i=1;i<r;i++)
  {
    pr=(GEN)list[i];v=element_val(nf,a,pr);
    if(v){ep[i]=un;if(element_val(nf,b,pr)<=v) ve[i]=(long)int0;else ve[i]=(long)int1;}
    else
    {
      v=idealval(nf,abJ,pr);
      p1=element_div(nf,(GEN)idealaddone(nf,gtomat(a),pr)[1],a);
      ep[i]=lstoi(v+1);k=1;
      while(k<=v){p1=element_mul(nf,p1,gsub(int2,element_mul(nf,a,p1)));k<<=1;}
      p1=element_mul(nf,p1,gsub(element_pow(nf,(GEN)pr[2],stoi(v)),b));
      ve[i]=(long)gmod(p1,gpuigs((GEN)pr[1],v+1));
    }
  }
  fact[2]=(long)ep;p2=idealchinese(nf,fact,ve);
  if(DEBUGLEVEL>2)
  {fprintferr(" sortie de findX() : p2 = ");outerr(p2);flusherr();}
  return p2;
}  

GEN
threetotwo1(GEN nf, GEN a, GEN b, GEN c)
/* Given 3 algebraic integers a,b,c in a number field nf given by their components
   on the integral basis, gives a three-component vector [d,e,U] whose 
   first two components are algebraic integers d,e and the third a unimodular
   3x3-matrix U such that [a,b,c]*U=[0,d,e] */
{
  long av=avma,tetpil,i,N,r;
  GEN y,p1,p2,p3,id,M,X,Y,J,e,e1,b1,c1,u,v,U,int0,Z,pk,list,pr;
  
  if(DEBUGLEVEL>2)
  {
    fprintferr(" On entre dans threetotwo1() : \n");
    fprintferr(" a = ");outerr(a);fprintferr(" b = ");outerr(b);
    fprintferr(" c = ");outerr(c);flusherr();
  }
  if(gcmp0(a))
  {
    y=cgetg(4,17);y[1]=lcopy(b);y[2]=lcopy(c);y[3]=(long)idmat(3);
    return y;
  }
  if(gcmp0(b))
  {
    y=cgetg(4,17);y[1]=lcopy(a);y[2]=lcopy(c);
    p1=idmat(3);y[3]=(long)p1;p2=(GEN)p1[1];p1[1]=p1[2];p1[2]=(long)p2;
    return y;
   }
  if(gcmp0(c))
  {
    y=cgetg(4,17);y[1]=lcopy(a);y[2]=lcopy(b);
    p1=idmat(3);y[3]=(long)p1;p2=(GEN)p1[1];p1[1]=p1[3];p1[3]=p1[2];
    p1[2]=(long)p2;return y;
   }
  N=lgef((GEN)nf[1])-3;id=idmat(N);
  p1=cgetg(4,19);p1[1]=(long)a;p1[2]=(long)b;p1[3]=(long)c;p1=idealmul(nf,p1,id);
  if(DEBUGLEVEL>2)
  {fprintferr(" ideal a.Z_k+b.Z_k+c.Z_k = ");outerr(p1);flusherr();}
  if(gcmp1(gcoeff(p1,1,1)))
  {
/* cas ou a.Z_K+b.Z_K+c.Z_K=Z_K ; beaucoup plus simple et frequent, donc
   on traite separement */
    list=(GEN)idealfactor(nf,c)[1];r=lg(list);
    p2=id;p3=id;
    for(i=1;i<r;i++)
    {
      pr=(GEN)list[i];
      if(element_val(nf,b,pr)) p3=idealmulprime(nf,p3,pr);
      else if(!element_val(nf,a,pr)) p2=idealmulprime(nf,p2,pr);
    }
    X=(GEN)idealaddone(nf,p2,p3)[1];
    b1=gadd(b,element_mul(nf,a,X));c1=c;e=(GEN)id[1];
    if(DEBUGLEVEL>2)
    {
      fprintferr(" b1 = ");outerr(b1);
      fprintferr(" c1 = ");outerr(c1);flusherr();
    }
    p1=idealhermite(nf,gtomat(b1));p2=idealhermite(nf,gtomat(c1));
    p3=idealaddone(nf,p1,p2);
    u=element_div(nf,(GEN)p3[1],b1);v=element_div(nf,(GEN)p3[2],c1);
    if(DEBUGLEVEL>2)
    {
      fprintferr(" u = ");outerr(u);fprintferr(" v = ");
      outerr(v);flusherr();
    }
    U=cgetg(4,19);
    p1=cgetg(4,18);p2=cgetg(4,18);p3=cgetg(4,18);
    U[1]=(long)p1;U[2]=(long)p2;U[3]=(long)p3;
    p1[1]=(long)element_mul(nf,X,c1);
    p1[2]=(long)c1;p1[3]=lneg(b1);
    int0=cgetg(N+1,18);for(i=1;i<=N;i++) int0[i]=zero;
    p2[1]=id[1];p2[2]=p2[3]=(long)int0;
    Z=element_mul(nf,X,u);
    pk=nfreducemat(nf,c1,(GEN)p1[3],u,v);
    p3[1]=(long)int0;p3[2]=lsub(u,element_mul(nf,pk,c1));
    p3[3]=(long)gadd(v,element_mul(nf,pk,b1));
    e=gadd(e,element_mul(nf,a,gsub(element_mul(nf,pk,(GEN)p1[1]),Z)));
    tetpil=avma;
    y=cgetg(4,17);y[1]=lcopy(a);y[2]=lcopy(e);y[3]=lcopy(U);
    return gerepile(av,tetpil,y);
  }
  p1=idealinv(nf,p1);
  p1=gmul(p1,e=denom(p1));
  M=idealmul(nf,id,gtomat(element_mul(nf,element_mul(nf,a,b),c)));
  J=idealmul(nf,p1,gtomat(e1=idealcoprime(nf,p1,M)));
  e=gmul(e,e1);
  if(DEBUGLEVEL>2)
  {
    fprintferr(" ideal J = ");outerr(J);
    fprintferr(" e = ");outerr(e);
    flusherr();
  }
  p1=cgetg(3,19);p1[1]=(long)a;p1[2]=(long)b;M=idealmul(nf,p1,J);
  if(DEBUGLEVEL>2)
  {fprintferr(" ideal M=(a.Z_k+b.Z_k).J = ");outerr(M);flusherr();}
  X=findX(nf,a,b,J,M);
  if(DEBUGLEVEL>2){fprintferr(" X = ");outerr(X);flusherr();}
  p1=gadd(b,element_mul(nf,a,X));
  p2=cgetg(3,19);p2[1]=(long)element_mul(nf,a,p1);p2[2]=(long)element_mul(nf,c,p1);
  if(N==2) p2=idealmul(nf,p2,id);
  p3=cgetg(3,19);p3[1]=(long)a;p3[2]=(long)b;p3=idealmul(nf,p3,id);
  if(DEBUGLEVEL>2)
  {fprintferr(" ideal a.Z_k+b.Z_k = ");outerr(p3);flusherr();}
  Y=findX(nf,a,c,J,idealdiv(nf,p2,p3));
  if(DEBUGLEVEL>2){fprintferr(" Y = ");outerr(Y);flusherr();}
  b1=element_div(nf,p1,e);
  if(DEBUGLEVEL>2){fprintferr(" b1 = ");outerr(b1);flusherr();}
  p2=gadd(c,element_mul(nf,a,Y));
  c1=element_div(nf,p2,e);
  if(DEBUGLEVEL>2){fprintferr(" c1 = ");outerr(c1);flusherr();}
  p1=idealhermite(nf,gtomat(b1));p2=idealhermite(nf,gtomat(c1));
  p3=idealaddone(nf,p1,p2);
  u=element_div(nf,(GEN)p3[1],b1);v=element_div(nf,(GEN)p3[2],c1);
  if(DEBUGLEVEL>2)
  {
    fprintferr(" u = ");outerr(u);fprintferr(" v = ");outerr(v);
    flusherr();
  }
  U=cgetg(4,19);
  p1=cgetg(4,18);p2=cgetg(4,18);p3=cgetg(4,18);
  U[1]=(long)p1;U[2]=(long)p2;U[3]=(long)p3;
  p1[1]=lsub(element_mul(nf,X,c1),element_mul(nf,Y,b1));
  p1[2]=(long)c1;p1[3]=lneg(b1);
  int0=cgetg(N+1,18);for(i=1;i<=N;i++) int0[i]=zero;
  p2[1]=id[1];p2[2]=p2[3]=(long)int0;
  Z=gadd(element_mul(nf,X,u),element_mul(nf,Y,v));
  pk=nfreducemat(nf,c1,(GEN)p1[3],u,v);
  p3[1]=(long)int0;p3[2]=lsub(u,element_mul(nf,pk,c1));
  p3[3]=(long)gadd(v,element_mul(nf,pk,b1));
  e=gadd(e,element_mul(nf,a,gsub(element_mul(nf,pk,(GEN)p1[1]),Z)));
  tetpil=avma;
  y=cgetg(4,17);y[1]=lcopy(a);y[2]=lcopy(e);y[3]=lcopy(U);
  if(DEBUGLEVEL>2)
  {fprintferr(" sortie de threetotwo1() : y = ");outerr(y);flusherr();}
  return gerepile(av,tetpil,y);
}

GEN
threetotwo2(GEN nf, GEN a, GEN b, GEN c)
/* Given 3 algebraic integers a,b,c in a number field nf given by their components
   on the integral basis, gives a three-component vector [d,e,U] whose
   first two components are algebraic integers d,e and the third a unimodular
   3x3-matrix U such that [a,b,c]*U=[0,d,e] */
{
  long av=avma,tetpil,i,N;
  GEN y,p1,p2,p3,id,M,X,Y,J,e,b1,c1,u,v,U,int0,Z,pk;
  
  if(DEBUGLEVEL>2)
  {
    fprintferr(" On entre dans threetotwo2() : \n");
    fprintferr(" a = ");outerr(a);fprintferr(" b = ");outerr(b);
    fprintferr(" c = ");outerr(c);flusherr();
  }
  if(gcmp0(a))
  {
    y=cgetg(4,17);y[1]=lcopy(b);y[2]=lcopy(c);y[3]=(long)idmat(3);
    return y;
  }
  if(gcmp0(b))
  {
    y=cgetg(4,17);y[1]=lcopy(a);y[2]=lcopy(c);
    y[3]=(long)lisexpr("[0,1,0;1,0,0;0,0,1]");return y;
  }
  if(gcmp0(c))
  {
    y=cgetg(4,17);y[1]=lcopy(a);y[2]=lcopy(b);
    y[3]=(long)lisexpr("[0,1,0;0,0,1;1,0,0]");return y;
  }
  N=lgef((GEN)nf[1])-3;id=idmat(N);
  p1=cgetg(4,19);p1[1]=(long)a;p1[2]=(long)b;p1[3]=(long)c;p1=idealmul(nf,p1,id);
  if(DEBUGLEVEL>2)
  {fprintferr(" ideal a.Z_k+b.Z_k+c.Z_k = ");outerr(p1);flusherr();}
  J=idealdiv(nf,e=idealcoprimeinvabc(nf,p1,a,b,c),p1);
  if(DEBUGLEVEL>2)
  {
    fprintferr(" ideal J = ");outerr(J);
    fprintferr(" e = ");outerr(e);
    flusherr();
  }
  p1=cgetg(3,19);p1[1]=(long)a;p1[2]=(long)b;M=idealmul(nf,p1,J);
  if(DEBUGLEVEL>2)
  {fprintferr(" ideal M=(a.Z_k+b.Z_k).J = ");outerr(M);flusherr();}
  X=findX(nf,a,b,J,M);
  if(DEBUGLEVEL>2){fprintferr(" X = ");outerr(X);flusherr();}
  p1=gadd(b,element_mul(nf,a,X));
  p2=cgetg(3,19);p2[1]=(long)element_mul(nf,a,p1);p2[2]=(long)element_mul(nf,c,p1);
  if(N==2) p2=idealmul(nf,p2,id);
  p3=cgetg(3,19);p3[1]=(long)a;p3[2]=(long)b;p3=idealmul(nf,p3,id);
  if(DEBUGLEVEL>2)
  {fprintferr(" ideal a.Z_k+b.Z_k = ");outerr(p3);flusherr();}
  Y=findX(nf,a,c,J,idealdiv(nf,p2,p3));
  if(DEBUGLEVEL>2){fprintferr(" Y = ");outerr(Y);flusherr();}
  b1=element_div(nf,p1,e);
  if(DEBUGLEVEL>2){fprintferr(" b1 = ");outerr(b1);flusherr();}
  p2=gadd(c,element_mul(nf,a,Y));
  c1=element_div(nf,p2,e);
  if(DEBUGLEVEL>2){fprintferr(" c1 = ");outerr(c1);flusherr();}
  p1=idealhermite(nf,gtomat(b1));p2=idealhermite(nf,gtomat(c1));
  p3=idealaddone(nf,p1,p2);
  u=element_div(nf,(GEN)p3[1],b1);v=element_div(nf,(GEN)p3[2],c1);
  if(DEBUGLEVEL>2)
  {
    fprintferr(" u = ");outerr(u);fprintferr(" v = ");outerr(v);
    flusherr();
  }
  U=cgetg(4,19);
  p1=cgetg(4,18);p2=cgetg(4,18);p3=cgetg(4,18);
  U[1]=(long)p1;U[2]=(long)p2;U[3]=(long)p3;
  p1[1]=lsub(element_mul(nf,X,c1),element_mul(nf,Y,b1));
  p1[2]=(long)c1;p1[3]=lneg(b1);
  int0=cgetg(N+1,18);for(i=1;i<=N;i++) int0[i]=zero;
  p2[1]=id[1];p2[2]=p2[3]=(long)int0;
  Z=gadd(element_mul(nf,X,u),element_mul(nf,Y,v));
  pk=nfreducemat(nf,c1,(GEN)p1[3],u,v);
  p3[1]=(long)int0;p3[2]=lsub(u,element_mul(nf,pk,c1));
  p3[3]=(long)gadd(v,element_mul(nf,pk,b1));
  e=gadd(e,element_mul(nf,a,gsub(element_mul(nf,pk,(GEN)p1[1]),Z)));
  tetpil=avma;
  y=cgetg(4,17);y[1]=lcopy(a);y[2]=lcopy(e);y[3]=lcopy(U);
  if(DEBUGLEVEL>2)
  {fprintferr(" sortie de threetotwo2() : y = ");outerr(y);flusherr();}
  return gerepile(av,tetpil,y);
}

GEN
basistoalg(GEN nf, GEN x)
{
  long tx=typ(x),lx=lg(x),av=avma,tetpil,i;
  GEN p1,z;

  nf=checknf(nf);  
  switch(tx)
  {
    case 17: case 18:
      for(i=1;(i<lx)&&(typ((GEN)x[i])<17);i++);
      if(i==lx)
      {
	if(tx==17) x=gtrans(x);
	tetpil=avma;p1=gmul((GEN)nf[7],x);
	return gerepile(av,tetpil,gmodulcp(p1,(GEN)nf[1]));
      }
    case 19:
      z=cgetg(lx,tx);for(i=1;i<lx;i++) z[i]=(long)basistoalg(nf,(GEN)x[i]);
      return z;
    default: z=cgetg(3,19);z[1]=lcopy((GEN)nf[1]);
      z[2]=lmul(x,polun[varn((GEN)nf[1])]);return z;
  }
}

GEN
algtobasis(GEN nf, GEN x)
{
  long tx=typ(x),lx=lg(x),av=avma,tetpil,i,v,N;
  GEN y,z;
  
  nf=checknf(nf);  
  switch(tx)
  {
    case 17: case 18: case 19:
      z=cgetg(lx,tx);for(i=1;i<lx;i++) z[i]=(long)algtobasis(nf,(GEN)x[i]);
      return z;
    case 9:
      y=(GEN)x[1];v=varn(y);setvarn(y,varn((GEN)nf[1]));
      if(!gegal((GEN)nf[1],y))
      {
	setvarn(y,v);err(talker,"not the same number field in algtobasis");
      }
      x=(GEN)x[2];
    case 10: N=lgef((GEN)nf[1])-3;
      if((tx==10)&&((lgef(x)-3)>=N)) x=gmod(x,(GEN)nf[1]);
      z=cgetg(N+1,18);for(i=1;i<=N;i++) z[i]=(long)truecoeff(x,i-1);
      tetpil=avma;return gerepile(av,tetpil,gmul((GEN)nf[8],z));
    default: N=lgef((GEN)nf[1])-3;
      z=cgetg(N+1,18);for(i=2;i<=N;i++) z[i]=zero;z[1]=lcopy(x);
      return z;
  }
}

GEN
nfdiveuc(GEN nf, GEN a, GEN b)
/* Given a and b in nf, gives an algebraic integer y in nf such that a-b.y
   is "small" */
{
  long av=avma,tetpil;
  GEN p1;
  
  p1=element_div(nf,a,b);tetpil=avma;
  return gerepile(av,tetpil,ground(p1));
}

GEN
nfmod(GEN nf, GEN a, GEN b)
/* Given a and b in nf, gives a "small" algebraic integer r in nf
   of the form a-b.y is "small" */
{
  long av=avma,tetpil;
  GEN p1;
  
  p1=gneg(element_mul(nf,b,ground(element_div(nf,a,b))));tetpil=avma;
  return gerepile(av,tetpil,gadd(a,p1));
}

GEN
nfdivres(GEN nf, GEN a, GEN b)
/* Given a and b in nf, gives a two-component vector [y,r] in nf such
   that r=a-b.y is "small" */
{
  long av=avma,tetpil;
  GEN p1,y,z;
  
  p1=gneg(element_mul(nf,b,y=ground(element_div(nf,a,b))));tetpil=avma;
  z=cgetg(3,17);z[1]=lcopy(y);z[2]=(long)gadd(a,p1);
  return gerepile(av,tetpil,z);
}



GEN
nfreducemat(GEN nf, GEN a, GEN b, GEN c, GEN d)
/* A usage interne : pas de verif ni gestion de pile. Given a, b, c, d
   in nf, gives an algebraic integer y in nf such that [c,d]-y.[a,b]
   is "small" */
{
  long av=avma,tetpil,i,i1,i2,j,j1,j2,k,N;
  GEN p1,p2,X,M,y,mult,s;

  mult=(GEN)nf[9];N=lgef((GEN)nf[1])-3;X=cgetg(N+1,18);
  for(j=1;j<=N;j++)
  {
    s=gzero;
    for(i=1;i<=N;i++) for(k=1;k<=N;k++)
    {
      p1=gcoeff(mult,k,j+(i-1)*N);
      if(!gcmp0(p1))
	s=gadd(s,gmul(p1,gadd(gmul((GEN)a[i],(GEN)c[k]),gmul((GEN)b[i],(GEN)d[k]))));
    }
    X[j]=(long)s;
  }
  M=cgetg(N+1,19);
  for(j2=1;j2<=N;j2++)
  {
    p1=cgetg(N+1,18);M[j2]=(long)p1;
    for(j1=1;j1<=N;j1++)
    {
      s=gzero;
      for(i1=1;i1<=N;i1++) for(i2=1;i2<=N;i2++) for(k=1;k<=N;k++)
      {
	p2=gmul(gcoeff(mult,k,j1+(i1-1)*N),gcoeff(mult,k,j2+(i2-1)*N));
	if(!gcmp0(p2))
        s=gadd(s,gmul(p2,gadd(gmul((GEN)a[i1],(GEN)a[i2]),gmul((GEN)b[i1],(GEN)b[i2]))));
      }
      p1[j1]=(long)s;
    }
  }
  y=gauss(M,X);tetpil=avma;
  return gerepile(av,tetpil,ground(y));
}

GEN
threetotwo(GEN nf, GEN a, GEN b, GEN c)
/* Given 3 algebraic integers a,b,c in a number field nf given by their components
   on the integral basis, gives a three-component vector [d,e,U] whose
   first two components are algebraic integers d,e and the third a unimodular
   3x3-matrix U such that [a,b,c]*U=[0,d,e] Naive method which may not
   work, but fast and small coefficients. */
{
  long av=avma,tetpil,i,N;
  GEN pol,p1,p2,p3,p4,y,id,hu,h,V,U,r,vd,q1,q1a,q2,q2a,na,nb,nc,nr;
  
  nf=checknf(nf);  
  pol=(GEN)nf[1];N=lgef(pol)-3;id=idmat(N);
  na=gnorml2(a);nb=gnorml2(b);nc=gnorml2(c);
  U=gmul(idmat(3),gmodulcp(polun[varn(pol)],pol));
  if(gcmp(nc,nb)<0)
  {
    p1=c;c=b;b=p1;p1=nc;nc=nb;nb=p1;
    p1=(GEN)U[3];U[3]=U[2];U[2]=(long)p1;
  }
  if(gcmp(nc,na)<0)
  {
    p1=a;a=c;c=p1;p1=na;na=nc;nc=p1;
    p1=(GEN)U[1];U[1]=U[3];U[3]=(long)p1;
  }
  while(!gcmp0(gmin(na,nb)))
  {
    p1=cgetg(2*N+1,19);
    for(i=1;i<=N;i++)
    {
      p1[i]=(long)element_mul(nf,a,(GEN)id[i]);
      p1[i+N]=(long)element_mul(nf,b,(GEN)id[i]);
    }
    hu=hnfnew(p1);h=(GEN)hu[1];V=(GEN)hu[2];
    p2=(GEN)ker(concat(h,c))[1];p3=(GEN)p2[N+1];
    p4=cgetg(N+1,18);
    for(i=1;i<=N;i++) p4[i]=(long)ground(gdiv((GEN)p2[i],p3));
    r=gadd(c,gmul(h,p4));
    vd=cgetg(N+1,19);for(i=1;i<=N;i++) vd[i]=V[N+i];
    p2=gmul(vd,p4);
    q1=cgetg(N+1,18);q2=cgetg(N+1,18);
    for(i=1;i<=N;i++) {q1[i]=p2[i];q2[i]=p2[i+N];}
    q1a=basistoalg(nf,q1);q2a=basistoalg(nf,q2);
    U[3]=(long)gadd((GEN)U[3],gadd(gmul(q1a,(GEN)U[1]),gmul(q2a,(GEN)U[2])));
    nr=gnorml2(r);
    if(gcmp(nr,gmax(na,nb))>=0)
    {
      err(talker,"threetotwo does not work");
    }
    if(gcmp(na,nb)>=0)
    {
      c=a;nc=na;a=r;na=nr;p1=(GEN)U[1];U[1]=U[3];U[3]=(long)p1;
    }
    else
    {
      c=b;nc=nb;b=r;nb=nr;p1=(GEN)U[2];U[2]=U[3];U[3]=(long)p1;
    }
  }
  if(!gcmp0(na))
  {
    p1=a;a=b;b=p1;p1=(GEN)U[1];U[1]=U[2];U[2]=(long)p1;
  }
  tetpil=avma;y=cgetg(4,17);y[1]=lcopy(b);y[2]=lcopy(c);
  y[3]=(long)algtobasis(nf,U);return gerepile(av,tetpil,y);
}

GEN
twototwo(GEN nf, GEN a, GEN b)
/* Given 2 algebraic integers a,b in a number field nf given by their
   components on the integral basis, gives a three-components vector [d,e,U]
   whose first two component are algebraic integers d,e and the third a
   unimodular 2x2-matrix U such that [a,b]*U=[d,e], with d and e hopefully
   smaller than a and b. */
{
  long av=avma,tetpil,N,fl;
  GEN pol,p1,y,id,U,r,qr,qa,na,nb,nr;
  
  nf=checknf(nf);  
  pol=(GEN)nf[1];N=lgef(pol)-3;id=idmat(N);
  na=gnorml2(a);nb=gnorml2(b);
  U=gmul(idmat(2),gmodulcp(polun[varn(pol)],pol));
  if(gcmp(na,nb)>0)
  {
    p1=a;a=b;b=p1;p1=na;na=nb;nb=p1;
    p1=(GEN)U[2];U[2]=U[1];U[1]=(long)p1;
  }
  fl=1;
  while((!gcmp0(na))&&fl)
  {
    qr=nfdivres(nf,b,a);r=(GEN)qr[2];nr=gnorml2(r);
    if(gcmp(nr,na)<0)
    {
      b=a;a=r;nb=na;na=nr;qa=basistoalg(nf,(GEN)qr[1]);
      p1=gsub((GEN)U[2],gmul(qa,(GEN)U[1]));U[2]=U[1];U[1]=(long)p1;
    }
    else
    {
      fl=0;
      if(gcmp(nr,nb)<0)
      {
	qa=basistoalg(nf,(GEN)qr[1]);
	U[2]=lsub((GEN)U[2],gmul(qa,(GEN)U[1]));
      }
    }
  }
  tetpil=avma;y=cgetg(4,17);y[1]=lcopy(a);y[2]=lcopy(b);
  y[3]=(long)algtobasis(nf,U);return gerepile(av,tetpil,y);
}

GEN
element_mulvec(GEN nf, GEN x, GEN v)
{
  long lx=lg(v),i;
  GEN y;
  
  y=cgetg(lx,typ(v));
  for(i=1;i<lx;i++) y[i]=(long)element_mul(nf,x,(GEN)v[i]);
  return y;
}

GEN
element_mulvecrow(GEN nf, GEN x, GEN m, long i, long lim)
{
  long lx,j;
  GEN y;

  lx=min(lg(m),lim+1);
  y=cgetg(lx,17);
  for(j=1;j<lx;j++) y[j]=(long)element_mul(nf,x,gcoeff(m,i,j));
  return y;
}

GEN
element_reduce(GEN nf, GEN x, GEN ideal)
{
/* Given an element x and an ideal in HNF, gives an r such that
   x-r is in ideal and r is small. Pour l'instant a usage interne,
   donc pas de verifs. */

  long av=avma,tetpil,N,j;
  GEN p1,p2,p3;
  
  N=lg(x)-1;
  p1=cgetg(N+2,19);for(j=1;j<=N;j++) p1[j]=ideal[j];
  p1[N+1]=(long)x;
  p2=(GEN)ker(p1)[1];p3=(GEN)p2[N+1];
  p1=cgetg(N+1,18);
  for(j=1;j<=N;j++) p1[j]=lround(gdiv((GEN)p2[j],p3));
  p2=gmul(ideal,p1);
  tetpil=avma;return gerepile(av,tetpil,gadd(p2,x));
}

/* A torsion-free module M over Z_K will be given by a row vector [A,I]
   with two components. I=[\a_1,...,\a_k] is a row vector of k fractional
   ideals given in HNF. A is an nxk matrix (same k and n the rank of the
   module) such that if A_j is the j-th column of A then
   M=\a_1A_1+...\a_kA_k. We say that [A,I] is a pseudo-basis if k=n */

GEN
nfhermite(GEN nf, GEN x)
{
/* Given a torsion-free module x as above outputs a pseudo-basis for x in
 Hermite Normal Form */

  long av=avma,tetpil,i,j,def,k,m;
  GEN p1,p2,p3,p4,p5,p6,y,newid,A,I,J,u,v,den;

  nf=checknf(nf);  
  if((typ(x)!=17)||(lg(x)!=3)) err(talker,"not a module in nfhermite");
  A=(GEN)x[1];I=(GEN)x[2];
  if(typ(A)!=19) err(talker,"not a matrix in nfhermite");
  k=lg(A)-1;
  if((typ(I)!=17)||(lg(I)!=k+1))
    err(talker,"not a correct ideal list in nfhermite");
  if(!k) err(talker,"not a matrix of maximal rank in nfhermite");
  m=lg((GEN)A[1])-1;
  if(k<m) err(talker,"not a matrix of maximal rank in nfhermite");
  A=gcopy(A);I=gcopy(I);def=k+1;
  for(j=1;j<=k;j++)
    if(typ((GEN)I[j])!=19) I[j]=(long)idealhermite(nf,(GEN)I[j]);
  J=cgetg(k+1,17);for(j=1;j<=k;j++) J[j]=zero;
  for(i=m;i>=1;i--)
  {
    def--;
    for(j=def;(j>=1)&&gcmp0(gcoeff(A,i,j));j--);
    if(!j) err(talker,"not a matrix of maximal rank in nfhermite");
    if(j<def)
    {
      p1=(GEN)A[j];A[j]=A[def];A[def]=(long)p1;
      p1=(GEN)I[j];I[j]=I[def];I[def]=(long)p1;
    }
    p2=element_inv(nf,p3=gcoeff(A,i,def));
    A[def]=(long)element_mulvec(nf,p2,(GEN)A[def]);
    I[def]=(long)idealmul(nf,p3,(GEN)I[def]);
    for(j=def-1;j>=1;j--)
    {
      if(!gcmp0(p1=gcoeff(A,i,j)))
      {
	p2=idealmul(nf,p1,(GEN)I[j]);
	newid=idealadd(nf,p2,(GEN)I[def]);
	p3=idealinv(nf,newid);
	p4=idealmul(nf,p2,p3);p5=idealmul(nf,(GEN)I[def],p3);
	y=idealaddone(nf,p4,p5);
	u=element_div(nf,(GEN)y[1],p1);v=(GEN)y[2];
	p6=gsub((GEN)A[j],element_mulvec(nf,p1,(GEN)A[def]));
	A[def]=(long)gadd(element_mulvec(nf,u,(GEN)A[j]),element_mulvec(nf,v,(GEN)A[def]));
	A[j]=(long)p6;
	I[j]=(long)idealmul(nf,idealmul(nf,(GEN)I[j],(GEN)I[def]),p3);
	I[def]=(long)newid;
	den=denom((GEN)I[j]);
	if(!gcmp1(den))
	{
	  I[j]=lmul(den,(GEN)I[j]);A[j]=(long)gdiv((GEN)A[j],den);
	}
      }
    }
    p1=(GEN)I[def];
    J[def]=(long)idealinv(nf,p1);
    for(j=def+1;j<=k;j++)
    {
      p2=gsub(element_reduce(nf,gcoeff(A,i,j),idealmul(nf,p1,(GEN)J[j])),gcoeff(A,i,j));
      A[j]=(long)gadd((GEN)A[j],element_mulvec(nf,p2,(GEN)A[def]));
    }
  }
  tetpil=avma;y=cgetg(3,17);
  p1=cgetg(m+1,19);for(j=1;j<=m;j++) p1[j]=lcopy((GEN)A[j+k-m]);
  y[1]=(long)p1;
  p2=cgetg(m+1,17);for(j=1;j<=m;j++) p2[j]=lcopy((GEN)I[j+k-m]);
  y[2]=(long)p2;return gerepile(av,tetpil,y);
}

/* A torsion module M over Z_K will be given by a row vector [A,I,J]
   with three components. I=[\b_1,...,\b_n] is a row vector of k
   fractional ideals given in HNF, J=[\a_1,...,\a_n] is a row vector of
   n fractional ideals in HNF. A is an nxn matrix (same n)
   such that if A_j is the j-th column of A and e_n is the canonical
   basis of K^n, then M=(\b_1e_1+...+\b_ne_n)/(\a_1A_1+...\a_nA_n)
   */

GEN
nfsmith(GEN nf, GEN x)
{
/* We input a torsion module x=[A,I,J] as above, and output the
   smith normal form as K=[\c_1,...,\c_n] such that
   x=Z_K/\c_1+...+Z_K/\c_n.
   */
  
  long av,tetpil,i,j,k,l,lim,c,fl,n,m,N;
  GEN p1,p2,p3,p4,z,b,u,v,w,d,dinv,unnf,A,I,J;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if((typ(x)!=17)||(lg(x)!=4)) err(talker,"not a module in nfsmith");
  A=(GEN)x[1];I=(GEN)x[2];J=(GEN)x[3];
  if(typ(A)!=19) err(talker,"not a matrix in nfsmith");
  n=lg(A)-1;
  if((typ(I)!=17)||(lg(I)!=n+1)||(typ(J)!=17)||(lg(J)!=n+1))
    err(talker,"not a correct ideal list in nfsmith");
  if(!n) err(talker,"not a matrix of maximal rank in nfsmith");
  m=lg((GEN)A[1])-1;
  if(n<m) err(talker,"not a matrix of maximal rank in nfsmith");
  if(n>m) err(impl,"nfsmith for non square matrices");
  A=gcopy(A);I=gcopy(I);J=gcopy(J);
  for(j=1;j<=n;j++)
    if(typ((GEN)I[j])!=19) I[j]=(long)idealhermite(nf,(GEN)I[j]);
  for(j=1;j<=n;j++)
    if(typ((GEN)J[j])!=19) J[j]=(long)idealhermite(nf,(GEN)J[j]);
  lim=(avma+bot)>>1;
  av=avma;
  for(i=n;i>=2;i--)
  {
    do
    {
      c=0;
      for(j=i-1;j>=1;j--)
      {
	p1=gcoeff(A,i,j);
	if(!gcmp0(p1))
	{
	  p2=gcoeff(A,i,i);
	  d=nfbezout(nf,p2,p1,(GEN)J[i],(GEN)J[j],&u,&v,&w,&dinv);
	  if(!gcmp0(u))
	  {
	    if(!gcmp0(v))
	      b=gadd(element_mulvec(nf,u,(GEN)A[i]),element_mulvec(nf,v,(GEN)A[j]));
	    else b=element_mulvec(nf,u,(GEN)A[i]);
	  }
	  else b=element_mulvec(nf,v,(GEN)A[j]);
	  A[j]=lsub(element_mulvec(nf,p2,(GEN)A[j]),element_mulvec(nf,p1,(GEN)A[i]));
	  A[i]=(long)b;
	  J[j]=(long)w;J[i]=(long)d;
	}
      }
      for(j=i-1;j>=1;j--)
      {
	p1=gcoeff(A,j,i);
	if(!gcmp0(p1))
	{
	  p2=gcoeff(A,i,i);
	  d=nfbezout(nf,p2,p1,(GEN)I[i],(GEN)I[j],&u,&v,&w,&dinv);
	  if(!gcmp0(u))
	  {
	    if(!gcmp0(v))
	      b=gadd(element_mulvecrow(nf,u,A,i,i),element_mulvecrow(nf,v,A,j,i));
	    else b=element_mulvecrow(nf,u,A,i,i);
	  }
	  else b=element_mulvecrow(nf,v,A,j,i);
	  p3=gsub(element_mulvecrow(nf,p2,A,j,i),element_mulvecrow(nf,p1,A,i,i));
	  for(k=1;k<=i;k++) {coeff(A,j,k)=p3[k];coeff(A,i,k)=b[k];}
	  I[j]=(long)w;I[i]=(long)d;
	  c++;
	}
      }
      if(!c)
      {
	b=gcoeff(A,i,i);fl=1;
	if(!gcmp0(b))
	{
	  b=idealmul(nf,gtomat(b),idealmul(nf,(GEN)J[i],(GEN)I[i]));
	  for(k=1;(k<i)&&fl;k++)
	    for(l=1;(l<i)&&fl;l++)
	    {
	      p3=gcoeff(A,k,l);
	      if(!gcmp0(p3)) fl=gegal(idealadd(nf,b,idealmul(nf,gtomat(p3),idealmul(nf,(GEN)J[l],(GEN)I[k]))),b);
	    }
	  if(!fl)
	  {
	    k--;l--;
	    b=idealdiv(nf,(GEN)I[k],(GEN)I[i]);
	    p4=invmulmat(idealdiv(nf,(GEN)J[i],idealmul(nf,gtomat(p3),(GEN)J[l])),b);
	    for(l=1;(l<=N)&&gcmp1(denom((GEN)p4[l]));l++);
	    if(l>N) err(talker,"bug2 in nfsmith");
	    p3=element_mulvecrow(nf,(GEN)b[l],A,k,i);
	    for(l=1;l<=i;l++)
	      coeff(A,i,l)=(long)gadd(gcoeff(A,i,l),(GEN)p3[l]);
	  }
	}
      }
      if(avma<lim) {tetpil=avma;A=gerepile(av,tetpil,gcopy(A));}
    }
    while(c||(!fl));
  }
  unnf=cgetg(N+1,18);unnf[1]=un;for(i=2;i<=N;i++) unnf[i]=zero;
  p1=gtomat(gcoeff(A,1,1));coeff(A,1,1)=(long)unnf;
  J[1]=(long)idealmul(nf,p1,(GEN)J[1]);
  for(i=2;i<=n;i++)
    if(!gegal(gcoeff(A,i,i),unnf))
      err(talker,"bug in nfsmith");
  tetpil=avma;z=cgetg(n+1,17);
  for(i=1;i<=n;i++) z[i]=(long)idealmul(nf,(GEN)I[i],(GEN)J[i]);
  return gerepile(av,tetpil,z);
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*          ALGEBRE LINEAIRE DANS LES CORPS DE NOMBRES             */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

GEN
element_divmodpr(GEN nf, GEN x, GEN y, GEN prhall)
/* On ne peut PAS definir ca comme les autres par
#define element_divmodpr(nf,x,y,prhall) (nfreducemodpr(nf,element_div(nf,x,y),prhall))
car le element_div ne marche en general pas */
{
  GEN p1,z;
  long i,N=lgef((GEN)nf[1])-3,av=avma,tetpil;
  
  p1=lift(gdiv(gmodulcp(gmul((GEN)nf[7],lift(x)),(GEN)nf[1]),gmodulcp(gmul((GEN)nf[7],lift(y)),(GEN)nf[1])));
  z=cgetg(N+1,18);
  for(i=1;i<=N;i++) z[i]=(long)truecoeff(p1,i-1);
  p1=gmul((GEN)nf[8],z);tetpil=avma;
  return gerepile(av,tetpil,nfreducemodpr(nf,p1,prhall));
}

GEN
element_invmodpr(GEN nf, GEN y, GEN prhall)
{
  GEN p1,z;
  long i,N=lgef((GEN)nf[1])-3,av=avma,tetpil;

  p1=ginvmod(gmul((GEN)nf[7],lift(y)),(GEN)nf[1]);
  z=cgetg(N+1,18);for(i=1;i<=N;i++) z[i]=(long)truecoeff(p1,i-1);
  p1=gmul((GEN)nf[8],z);tetpil=avma;
  return gerepile(av,tetpil,nfreducemodpr(nf,p1,prhall));
}

GEN
element_powmodpr(GEN nf,GEN x,GEN k,GEN prhall)
{
  long i,f,av=avma,tetpil,N,s;
  GEN k1,y,z;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  s=signe(k);k1=(s>=0)?k:negi(k);z=x;f=1;y=cgetg(N+1,18);
  for(i=2;i<=N;i++) y[i]=zero;
  y[1]=un;
  while(f)
  {
    if(mpodd(k1)) y=element_mulmodpr(nf,z,y,prhall);
    k1=shifti(k1,-1);f=signe(k1);
    if(f) z=element_sqrmodpr(nf,z,prhall);
  }
  tetpil=avma;return gerepile(av,tetpil,(s>=0)?gcopy(y):element_invmodpr(nf,y,prhall));
}

GEN
nfker(GEN nf, GEN x, GEN prhall)
/* x est une matrice dont les coefficients sont des vecteurs dans
   la base d'entiers modulo un ideal premier prhall, sous forme reduite
   modulo prhall. */
{
  long i,j,k,r,t,n,n1,m,av,av1,av2,N;
  GEN c,d,y,unnf,munnf,unmodp,zeronf,p,pp,prh;
  
  av=avma;prh=(GEN)prhall[1];
  N=lgef(nf[1])-3;pp=gcoeff(prh,1,1);unmodp=gmodulcp(gun,pp);
  unnf=cgetg(N+1,18);for(i=2;i<=N;i++) unnf[i]=zero;unnf[1]=un;unnf=gmul(unnf,unmodp);
  munnf=cgetg(N+1,18);for(i=2;i<=N;i++) munnf[i]=zero;munnf[1]=lneg(gun);
  munnf=gmul(munnf,unmodp);
  zeronf=cgetg(N+1,18);for(i=1;i<=N;i++) zeronf[i]=zero;zeronf=gmul(unmodp,zeronf);
  if (typ(x)!=19) err(kerer1);
  n1=lg(x);n=n1-1;if(!n) return cgetg(1,19);
  m=lg((GEN)x[1])-1;x=gcopy(x);
  r=0;
  c=cgeti(m+1);for(k=1;k<=m;k++) c[k]=0;
  d=cgeti(n1);
  av1=avma;
  for(k=1;k<=n;k++)
  {
    j=1;
    while((j<=m)&&(c[j]||gcmp0(gcoeff(x,j,k)))) j++;
    if (j<=m)
    {
      p=element_divmodpr(nf,munnf,gcoeff(x,j,k),prhall);
      coeff(x,j,k)=(long)munnf;
      for(i=k+1;i<=n;i++)
	coeff(x,j,i)=(long)element_mulmodpr(nf,p,gcoeff(x,j,i),prhall);
      for(t=1;t<=m;t++)
	if(t!=j)
	{
	  p=gcoeff(x,t,k);
	  for(i=k+1;i<=n;i++)
           coeff(x,t,i)=(long)gadd(gcoeff(x,t,i),element_mulmodpr(nf,p,gcoeff(x,j,i),prhall));
	  coeff(x,t,k)=(long)zeronf;
	}
      c[j]=k;d[k]=j;
      av2=avma;
      x=gerepile(av1,av2,gcopy(x));
    }		  
    else {r++;d[k]=0;}
  }
  if(r)
  {
    av1=avma;
    y=cgetg(r+1,19);
    for(j=k=1;j<=r;j++,k++)
    {
      while(d[k]) k++;
      y[j]=(long)(p=cgetg(n1,18));
      for(i=1;i<k;i++) p[i]=d[i]? lcopy(gcoeff(x,d[i],k)):lcopy(zeronf);
      p[k]=lcopy(unnf);
      for(i=k+1;i<=n;i++) p[i]=lcopy(zeronf);
    }
    return gerepile(av,av1,y);
  }
  else {avma=av;y=cgetg(1,19);return y;}
}

GEN
nfsuppl(GEN nf, GEN x, long n, GEN prhall)
{
  long i,j,av=avma,tetpil,k,s,t,N;
  GEN y,p1,p2,pp,unmodp,unnf,prh;

  N=lgef(nf[1])-3;prh=(GEN)prhall[1];
  pp=gcoeff(prh,1,1);unmodp=gmodulcp(gun,pp);
  unnf=cgetg(N+1,18);for(i=2;i<=N;i++) unnf[i]=zero;unnf[1]=un;unnf=gmul(unnf,unmodp);
  if(typ(x)!=19) err(kerer1);
  k=lg(x)-1;if(k>n) err(suppler2);
  if(k&&((lg((GEN)x[1])-1)!=n)) err(talker,"incorrect dimension in nfsuppl");
  s=0;y=cgetg(n+1,19);
  for(j=1;j<=n;j++)
  {
    p1=cgetg(n+1,18);y[j]=(long)p1;
    for(i=1;i<=n;i++){if(i==j) p1[i]=(long)unnf;else p1[i]=(long)gmul(gzero,unnf);}
  }
  while(s<k)
  {
    s++;p1=nfgauss(nf,y,(GEN)x[s],prhall);t=s;
    while((t<=n)&&gcmp0((GEN)p1[t])) t++;
    if(t>n) err(suppler2);
    p2=(GEN)y[s];y[s]=x[s];if(s!=t) y[t]=(long)p2;
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(y));
}

GEN
nfgauss(GEN nf, GEN a, GEN b, GEN prhall)
/* a.x=b ou b est un vecteur (sans verif.) */

{
  long  nbli,nbco,i,j,k,av1,av2,av3,av4;
  GEN   aa,x,p,m,u;
  
  nbco=lg(a)-1;nbli=lg((GEN)a[1])-1;if (nbco!=nbli) err(gausser1);
  x=cgetg(nbli+1,18);av1=avma;
  for (j=1;j<=nbco;j++) x[j]=b[j];
  aa=cgetg(nbco+1,19);
  for (j=1;j<=nbco;j++)
  {
    aa[j]=lgetg(nbli+1,18);
    for (i=1;i<=nbli;i++) coeff(aa,i,j)=coeff(a,i,j);
  }
  for (i=1;i<nbli;i++)
  {
    p=gcoeff(aa,i,i);k=i;
    if (gcmp0(p))
    {
      for (k=i+1;(k<=nbli)&&gcmp0(gcoeff(aa,k,i));k++);
      if (k>nbco) err(matinv1);
      else
      {
	for (j=i;j<=nbco;j++)
	{
	  u=gcoeff(aa,i,j);coeff(aa,i,j)=coeff(aa,k,j);
	  coeff(aa,k,j)=(long)u;
	}
	u=(GEN)x[i];x[i]=x[k];x[k]=(long)u;
	p=gcoeff(aa,i,i);
      }
    }
    for (k=i+1;k<=nbli;k++)
    {
      m=gcoeff(aa,k,i);
      if (!gcmp0(m))
      {
	m=element_divmodpr(nf,m,p,prhall);
	for (j=i+1;j<=nbco;j++)
	  coeff(aa,k,j)=lsub(gcoeff(aa,k,j),element_mulmodpr(nf,m,gcoeff(aa,i,j),prhall));
	x[k]=lsub((GEN)x[k],element_mulmodpr(nf,m,(GEN)x[i],prhall));
      }
    }
  }
      /* Resolution systeme triangularise */
  av2=avma;
  p=gcoeff(aa,nbli,nbco);
  if (gcmp0(p)) err(matinv1);
  else
  {
    x[nbli]=(long)element_divmodpr(nf,(GEN)x[nbli],p,prhall);
    for (i=nbli-1;i>0;i--)
    {
      av3=avma;m=(GEN)x[i];
      for (j=i+1;j<=nbco;j++)
	m=gsub(m,element_mulmodpr(nf,gcoeff(aa,i,j),(GEN)x[j],prhall));
      av4=avma;
      x[i]=lpile(av3,av4,element_divmodpr(nf,m,gcoeff(aa,i,i),prhall));
    }
  }
  gerepile(av1,av2,(GEN)1);
  return x;
}

GEN
nfidealdet1(GEN nf, GEN a, GEN b)
/* Given two fractional ideals a and b, gives x in a, y in b, z in b^-1,
   t in a^-1 such that xt-yz=1. In the present version, z=1. */
{
  long av=avma,tetpil,i,N,fla,flb;
  GEN x,p1,p2,res,z,da,db;

  N=lgef((GEN)nf[1])-3;
  da=denom(a);db=denom(b);fla=gcmp1(da);flb=gcmp1(db);
  if(!fla) a=gmul(da,a);if(!flb) b=gmul(db,b);
  x=idealcoprimeinv(nf,a,b);p1=idealmul(nf,x,idealinv(nf,a));
  p2=idealaddone(nf,p1,b);tetpil=avma;res=cgetg(5,17);
  res[1]=fla?lcopy(x):(long)gdiv(x,da);
  res[2]=flb?lcopy((GEN)p2[2]):(long)gdiv((GEN)p2[2],db);
  z=cgetg(N+1,18);z[1]=flb?un:lcopy(db);
  for(i=2;i<=N;i++) z[i]=zero;res[3]=(long)z;
  res[4]=(long)element_div(nf,(GEN)p2[1],(GEN)res[1]);
  return gerepile(av,tetpil,res);
}

GEN
nfdetint(GEN nf,GEN pseudo)
/* Given a pseudo basis pseudo, outputs a multiple of the ideal determinant of pseudo */
{
  GEN pass,c,v,det1,piv,pivprec,vi,p1,x,I,unnf,zeronf,id,idprod;
  long i,j,k,rg,t,n,n1,m,m1,av=avma,av1,av3,tetpil,lim,dec,cm=0,N;
  
  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if((typ(pseudo)!=17)||(lg(pseudo)!=3)) err(talker,"not a module in nfdetint");
  x=(GEN)pseudo[1];I=(GEN)pseudo[2];
  if(typ(x)!=19) err(talker,"not a matrix in nfdetint");
  n=(n1=lg(x))-1;if(!n) return gun;
  m=(m1=lg((GEN)x[1]))-1;
  if((typ(I)!=17)||(lg(I)!=n1))
    err(talker,"not a correct ideal list in nfdetint");
  lim=(avma+bot)>>1;
  unnf=cgetg(N+1,18);unnf[1]=un;for(i=2;i<=N;i++) unnf[i]=zero;
  zeronf=cgetg(N+1,18);for(i=1;i<=N;i++) zeronf[i]=zero;
  c=cgeti(m1);for(k=1;k<=m;k++) c[k]=0;id=idmat(N);
  av1=avma;det1=gscalmat(gzero,N);piv=pivprec=unnf;pass=cgetg(m1,19);
  for(j=1;j<=m;j++) 
  {
    p1=cgetg(m1,18);pass[j]=(long)p1;
    for(i=1;i<=m;i++) p1[i]=(long)zeronf;
  }
  v=cgetg(m1,18);k=1;rg=0;
  while((k<=n)&&(rg<m))
  {
    for(t=0,i=1;i<=m;i++)
      if(!c[i])
      {
	vi=element_mul(nf,piv,gcoeff(x,i,k));
	for(j=1;j<=m;j++)
	  if(c[j]) vi=gadd(vi,element_mul(nf,gcoeff(pass,i,j),gcoeff(x,j,k)));
	v[i]=(long)vi;if(!t) if(!gcmp0(vi)) t=i;
      }
    if(t)
    {
      rg++;c[t]=k;pivprec=piv;piv=(GEN)v[t];
      if(rg<m)
      {
	for(i=1;i<=m;i++)
	  if(!c[i])
	    for(j=1;j<=m;j++)
	      if(c[j]&&(j!=t))
	      {
		p1=gsub(element_mul(nf,piv,gcoeff(pass,i,j)),element_mul(nf,(GEN)v[i],gcoeff(pass,t,j)));
		coeff(pass,i,j)=(rg>1)?(long)element_div(nf,p1,pivprec):(long)p1;
	      }
	for(i=1;i<=m;i++) if(!c[i]) coeff(pass,i,t)=lneg((GEN)v[i]);
      }
    }
    if(rg==m)
    {
      if(!cm)
      {
	idprod=id;for(i=1;i<=m;i++) if(i!=t) idprod=idealmul(nf,idprod,(GEN)I[c[i]]);
	cm=1;
      }
      p1=idealmul(nf,piv,(GEN)I[c[t]]);
      det1=idealadd(nf,p1,det1);piv=pivprec;c[t]=0;rg--;
    }
    k++;
    if(avma<lim) 
    {
      tetpil=avma;det1=gcopy(det1);piv=gcopy(piv);pivprec=gcopy(pivprec);
      pass=gcopy(pass);v=gcopy(v);idprod=gcopy(idprod);
      av3=avma;dec=lpile(av1,tetpil,0)>>TWOPOTBYTES_IN_LONG;
      if(adecaler(det1,tetpil,av3)) det1+=dec;
      if(adecaler(piv,tetpil,av3)) piv+=dec;
      if(adecaler(pivprec,tetpil,av3)) pivprec+=dec;
      if(adecaler(idprod,tetpil,av3)) idprod+=dec;
      if(adecaler(pass,tetpil,av3)) pass+=dec;
      if(adecaler(v,tetpil,av3)) v+=dec;
    }
  }
  if(rg+cm==m) {tetpil=avma;return gerepile(av,tetpil,idealmul(nf,idprod,det1));}
  else {avma=av;return gscalmat(gzero,N);}
}

GEN
nfcleanmod(GEN nf, GEN x, long lim, GEN detmat)
{
  long lx=lg(x),i;
  GEN y;

  if((lim<=0)||(lim>=lx)) lim=lx-1;
  y=cgetg(lx,18);
  for(i=1;i<=lim;i++) y[i]=(long)element_reduce(nf,(GEN)x[i],detmat);
  for(;i<lx;i++) y[i]=x[i];
  return y;
}

GEN
nfbezout(GEN nf, GEN a, GEN b, GEN ida, GEN idb, GEN *u, GEN *v, GEN *w, GEN *di)
{
/* a usage interne
   Given elements a,b, ideals ida, idb, outputs d=a.ida+b.idb and
   gives di=d^-1, w=ida.idb.di, u, v such that au+bv=1 and u in
   ida.di, v in idb.di. We assume ida, idb non-zero, but a and b can be
   zero. Error if a and b are both zero. */

  long av=avma,j,N,av1,tetpil,dec;
  GEN pa,pb,pab,pa1,pb1,pu,pv,pw,uv,d,dinv;

  if(gcmp0(a))
  {
    if(gcmp0(b)) err(talker,"both elements zero in nfbezout");
    pab=idealmul(nf,ida,idb);
    tetpil=avma;d=idealmulelt(nf,b,idb);dinv=idealinv(nf,d);
    N=lgef((GEN)nf[1])-3;pu=cgetg(N+1,18);
    for(j=1;j<=N;j++) pu[j]=zero;
    pw=idealmul(nf,pab,dinv);
    pv=element_inv(nf,b);
  }
  else
  {
    if(gcmp0(b))
    {
      pab=idealmul(nf,ida,idb);
      tetpil=avma;d=idealmulelt(nf,a,ida);dinv=idealinv(nf,d);
      N=lgef((GEN)nf[1])-3;pv=cgetg(N+1,18);
      for(j=1;j<=N;j++) pv[j]=zero;
      pw=idealmul(nf,pab,dinv);
      pu=element_inv(nf,a);
    }
    else
    {
      pa=idealmulelt(nf,a,ida);pb=idealmulelt(nf,b,idb);
      d=idealadd(nf,pa,pb);dinv=idealinv(nf,d);
      pa1=idealmullll(nf,pa,dinv);pb1=idealmullll(nf,pb,dinv);
      uv=idealaddone(nf,pa1,pb1);
      pab=idealmul(nf,ida,idb);
      tetpil=avma;
      pu=element_div(nf,(GEN)uv[1],a);
      pv=element_div(nf,(GEN)uv[2],b);
      pw=idealmul(nf,pab,dinv);
      d=gcopy(d);dinv=gcopy(dinv);
    }
  }
  av1=avma;dec=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;
  if(adecaler(pu,tetpil,av1)) pu+=dec;
  if(adecaler(pv,tetpil,av1)) pv+=dec;
  if(adecaler(pw,tetpil,av1)) pw+=dec;
  if(adecaler(d,tetpil,av1)) d+=dec;
  if(adecaler(dinv,tetpil,av1)) dinv+=dec;
  *u=pu;*v=pv;*w=pw;*di=dinv;return d;
}

GEN
idealmullll(GEN nf, GEN x, GEN y)
{
      /* A usage interne. Pas de verifs ni gestion de pile */
  
  GEN z,den;
  long fl;

  z=idealmul(nf,x,y);
  den=denom(z);fl=gcmp1(den);
  if(!fl) z=gmul(den,z);
  z=gmul(z,lllintpartial(z));
  if(!fl) z=gdiv(z,den);
  return z;
}

GEN
idealmulelt(GEN nf, GEN elt, GEN x)
{
      /* A usage interne. Pas de verifs ni gestion de pile */
  long lx=lg(x),j;
  GEN z;
  
  z=cgetg(lx,19);
  for(j=1;j<lx;j++) z[j]=(long)element_mul(nf,elt,(GEN)x[j]);
  return z;
}

GEN
idealdivlll(GEN nf, GEN x, GEN y)
{
      /* A usage interne. Pas de verifs ni gestion de pile */
  
  GEN z,den;
  long fl;

  z=idealdiv(nf,x,y);
  den=denom(z);fl=gcmp1(den);
  if(!fl) z=gmul(den,z);
  z=gmul(z,lllintpartial(z));
  if(!fl) z=gdiv(z,den);
  return z;
}

GEN
nfhermitemod(GEN nf, GEN pseudo, GEN detmat)
{
  long li,co,av,tetpil,i,j,jm1,def,ldef,lim,N;
  GEN b,q,w,p1,p2,y,d,u,v,den,x,I,J,dinv,unnf,wh;

  nf=checknf(nf);N=lgef((GEN)nf[1])-3;
  if((typ(pseudo)!=17)||(lg(pseudo)!=3)) err(talker,"not a module in nfhermitemod");
  x=(GEN)pseudo[1];I=(GEN)pseudo[2];
  if(typ(x)!=19) err(talker,"not a matrix in nfhermitemod");
  co=lg(x);
  if((typ(I)!=17)||(lg(I)!=co))
    err(talker,"not a correct ideal list in nfhermitemod");
  if(co==1) return cgetg(1,19);
  li=lg((GEN)x[1]);
  x=gcopy(x);I=gcopy(I);
  unnf=cgetg(N+1,18);unnf[1]=un;for(i=2;i<=N;i++) unnf[i]=zero;
  for(j=1;j<co;j++)
    if(typ((GEN)I[j])!=19) I[j]=(long)idealhermite(nf,(GEN)I[j]);
  if(typ(x)!=19) err(hnfer1);
  lim=(avma+bot)>>1;
  av=avma;den=denom(detmat);
  if(!gcmp1(den)) detmat=gmul(den,detmat);
  detmat=gmul(detmat,lllintpartial(detmat));
  y=gcopy(x);
  def=co;ldef=(li>co)?li-co+1:1;
  for(i=li-1;i>=ldef;i--)
  {
    def--;j=def-1;while(j&&gcmp0(gcoeff(y,i,j))) j--;
    while(j)
    {
      jm1=j-1;if(!jm1) jm1=def;
      d=nfbezout(nf,gcoeff(y,i,j),gcoeff(y,i,jm1),(GEN)I[j],(GEN)I[jm1],&u,&v,&w,&dinv);
      if(!gcmp0(u))
      {
	if(!gcmp0(v)) p1=gadd(element_mulvec(nf,u,(GEN)y[j]),element_mulvec(nf,v,(GEN)y[jm1]));
	else p1=element_mulvec(nf,u,(GEN)y[j]);
      }
      else p1=element_mulvec(nf,v,(GEN)y[jm1]);
      y[j]=lsub(element_mulvec(nf,gcoeff(y,i,j),(GEN)y[jm1]),element_mulvec(nf,gcoeff(y,i,jm1),(GEN)y[j]));
      y[j]=(long)nfcleanmod(nf,(GEN)y[j],i,idealdivlll(nf,detmat,w));
      y[jm1]=(long)nfcleanmod(nf,p1,i,idealmullll(nf,detmat,dinv));
      I[j]=(long)w;I[jm1]=(long)d;
      j--;while(j&&(gcmp0(gcoeff(y,i,j)))) j--;
      if(avma<lim) {tetpil=avma;y=gerepile(av,tetpil,gcopy(y));}
    }
  }
  b=detmat;
  wh=cgetg(li,19);def--;
  for(i=li-1;i>=1;i--)
  {
    d=nfbezout(nf,gcoeff(y,i,i+def),unnf,(GEN)I[i+def],b,&u,&v,&w,&dinv);
    wh[i]=(long)nfcleanmod(nf,element_mulvec(nf,u,(GEN)y[i+def]),i,idealmullll(nf,b,dinv));
    coeff(wh,i,i)=(long)unnf;
    I[i+def]=(long)d;
    if(i>1) b=idealmul(nf,b,dinv);
  }
  J=cgetg(li,17);J[1]=zero;
  for(j=2;j<li;j++) J[j]=(long)idealinv(nf,(GEN)I[j+def]);
  for(i=li-2;i>=1;i--)
  {
    for(j=i+1;j<li;j++)
    {
      q=idealmul(nf,(GEN)I[i+def],(GEN)J[j]);
      p1=gsub(element_reduce(nf,gcoeff(wh,i,j),q),gcoeff(wh,i,j));
      wh[j]=(long)gadd((GEN)wh[j],element_mulvec(nf,p1,(GEN)wh[i]));
    }
    if(avma<lim) {tetpil=avma;wh=gerepile(av,tetpil,gcopy(wh));}
  }
  tetpil=avma;p1=cgetg(3,17);p1[1]=lcopy(wh);
  p2=cgetg(li,17);p1[2]=(long)p2;
  for(j=1;j<li;j++) p2[j]=lcopy((GEN)I[j+def]);
  return gerepile(av,tetpil,p1);
}


