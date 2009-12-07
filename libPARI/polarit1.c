/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~                                                                   ~*/
/*~                     OPERATIONS ARITHMETIQUES                      ~*/
/*~                                                                   ~*/
/*~                         SUR LES POLYNOMES                         ~*/
/*~                                                                   ~*/
/*~                         (premiere partie)                         ~*/
/*~                                                                   ~*/
/*~                        copyright Babe Cool                        ~*/
/*~                                                                   ~*/
/*~                                                                   ~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

# include "genpari.h"
#define TRUE 1
#define FALSE 0

GEN gnorml1(GEN x, long PREC),laguer(GEN pol,long N,GEN y0,GEN EPS,long PREC);
GEN square_free_factorization(GEN pol);
GEN respm(GEN f1,GEN f2,GEN pm);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                           DIVISIBILITE                          */
/*                                                                 */
/*                 Renvoie 1 si y divise x, 0 sinon .              */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


int
gdivise(GEN x, GEN y)
{
  long    avmacourant,i;
  GEN   p1;

  avmacourant=avma;
  p1=gmod(x,y);i=gcmp0(p1);
  avma=avmacourant;
  return i;
}

int
poldivis(GEN x, GEN y, GEN *z)
{
  long av=avma;
  GEN p1,p2;

  p1=poldivres(x,y,&p2);
  if(signe(p2)) {avma=av;return 0;}
  cgiv(p2);*z=p1;return 1;
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                            REDUCTION                            */
/*                                                                 */
/*        Met sous forme de fraction une nfraction; sous           */
/*        forme de fr.rat une n.fr.rat; dans les autres cas        */
/*        creation d'une copie .                                   */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



GEN
gred(GEN x)
{
  long    tx,tetpil,l;
  GEN   y,p1,x1,x2,x3;

  tx=typ(x);l=avma;
  if ((tx==4)||(tx==5))
  {
    y=dvmdii((GEN)x[1],(GEN)x[2],&p1);
    if(!signe(p1)) cgiv(p1);
    else
    {
      p1=mppgcd((GEN)x[2],p1);tetpil=avma;
      y=cgetg(3,4);
      if(gcmp1(p1))
      {
        y[1]=lcopy((GEN)x[1]);
        y[2]=lcopy((GEN)x[2]);
      }
      else
      {
        y[1]=ldivii((GEN)x[1],p1);
        y[2]=ldivii((GEN)x[2],p1);
      }
      y=gerepile(l,tetpil,y);
    }
  }
  else if ((tx==13)||(tx==14))
  {
    x1=content((GEN)x[1]);x2=content((GEN)x[2]);x3=gdiv(x1,x2);
    x1=(gcmp0(x1))?(GEN)x[1]:gdiv((GEN)x[1],x1);
    x2=gdiv((GEN)x[2],x2);y=poldivres(x1,x2,&p1);
    if(gcmp0(p1)) {tetpil=avma;return gerepile(l,tetpil,gmul(x3,y));}
    else
    {
      p1=ggcd(x2,p1);y=cgetg(3,13);
      if(isscalar(p1)) {y[1]=(long)x1;y[2]=(long)x2;}
      else {y[1]=ldeuc(x1,p1);y[2]=ldeuc(x2,p1);}
      x1=numer(x3);x2=denom(x3);tetpil=avma;
      p1=cgetg(3,13);p1[1]=lmul(x1,(GEN)y[1]);p1[2]=lmul(x2,(GEN)y[2]);
      return gerepile(l,tetpil,p1);
    }
  }
  else y=gcopy(x);
  return y;
}

/*    REDUCTION SUR PLACE AVEC CHANGEMENT DE TYPE EVENTUEL  */

void
gredsp(GEN *px)
{
  long    tx,l,l1;
  GEN   x,y,p1;

  x= *px;tx=typ(x);
  if ((tx==4)||(tx==5))
  {
    l=avma;
    y=dvmdii((GEN)x[1],(GEN)x[2],&p1);
    if(!signe(p1))
    {
      cgiv(p1);l1=(long)(x+3);*px=gerepile(l1,l,y);
    }
    else
    {
      p1=mppgcd((GEN)x[2],p1);
      if(!gcmp1(p1))
      {
        mpdivz((GEN)x[1],p1,(GEN)x[1]);
        mpdivz((GEN)x[2],p1,(GEN)x[2]);
      }
      settyp(x,4);avma=l;
    }
  }
  else if ((tx==13)||(tx==14))
  {
    l=avma;
    y=poldivres((GEN)x[1],(GEN)x[2],&p1);
    if(!signe(p1))
    {
      cgiv(p1);l1=(long)(x+3);*px=gerepile(l1,l,y);
    }
    else
    {
      p1=primpart(ggcd((GEN)x[2],p1));
      if(isnonscalar(p1))
      {
        gdeucz((GEN)x[1],p1,(GEN)x[1]);
        gdeucz((GEN)x[2],p1,(GEN)x[2]);
      }
      settyp(y,13);avma=l;
    }
  }
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*               DIVISION EUCLIDIENNE DES POLYNOMES                */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

GEN
gdeuc(GEN x, GEN y)
{
  long  tx=typ(x),ty=typ(y),dx,dy,dz,vx,vy,i,j,l,tetpil,tetpil2,f;
  GEN z,p1,p2,p3;

  if(ty<10) return gdiv(x,y);
  if(tx<10) return gzero;
  if((tx!=10)||(ty!=10)) err(poler1);
  if(gcmp0(y)) err(poler4);
  dx=lgef(x)-3;dy=lgef(y)-3;vx=varn(x);vy=varn(y);
  if((vx>vy)||(dx<dy))
  {z=cgetg(3,10);z[1]=2;z[2]=zero;setvarn(z,vy);}
  else
  {
    if(vx<vy) z=gdiv(x,y);
    else
    {
      dz=dx-dy;
      z=cgetg(dz+3,10);z[1]=evalsigne(1)+evallgef(3+dz);setvarn(z,vx);
      f=gcmp1((GEN)y[dy+2]);
      if(f) z[dz+2]=lcopy((GEN)x[dx+2]);
      else
	z[dz+2]=ldiv((GEN)x[dx+2],(GEN)y[dy+2]);
      for(i=dx-1;i>=dy;--i)
      {
        l=avma;p1=((GEN)x[i+2]);tetpil2=0;
        for(j=i-dy+1;(j<=i)&&(j<=dz);j++)
        {
          p2=gmul((GEN)z[j+2],(GEN)y[i-j+2]);
          tetpil2=avma;p1=gsub(p1,p2);
        }
        tetpil=avma;
        if(f)
        {
          if(tetpil2) z[i-dy+2]=lpile(l,tetpil2,p1);
          else z[i-dy+2]=(long)p1;
        }
        else
        {
          p3=gdiv(p1,(GEN)y[dy+2]);
          z[i-dy+2]=lpile(l,tetpil,p3);
        }
      }
    }
  }
  return z;
}


GEN
gres(GEN x, GEN y)
{
  long  tx=typ(x),ty=typ(y),dx,dy,dz,i,j,k,f,l1,tetpil,vx,vy;
  GEN z,p1,p2,p4;

  if(ty<10) return gzero;
  if(ty!=10) err(poler2);
  if(gcmp0(y)) err(poler5);
  vy=varn(y);vx=gvar(x);if(vx>vy) return gcopy(x);
  if(tx!=10) err(poler2);
  dx=lgef(x)-3;dy=lgef(y)-3;
  if(dx<dy) return gcopy(x);
  if(vx<vy) {z=cgetg(3,10);z[1]=evallgef(2)+evalvarn(vx);z[2]=zero;return z;}
  dz=dx-dy;l1=avma;
  p4=cgetg(dz+3,10);p4[1]=evalsigne(1)+evallgef(3+dz);setvarn(p4,vx);
  p4[dz+2]=ldiv((GEN)x[dx+2],(GEN)y[dy+2]);
  for(i=dx-1;i>=dy;--i)
  {
    p1=((GEN)x[i+2]);
    for(j=i-dy+1;(j<=i)&&(j<=dz);j++)
    {
      p2=gmul((GEN)p4[j+2],(GEN)y[i-j+2]);
      p1=gsub(p1,p2);
    } 
    p4[i-dy+2]=(long)gdiv(p1,(GEN)y[dy+2]);
  }
  f=1;
  for(i=dy-1;(i>=0)&&f;--i)
  {
    p1=((GEN)x[i+2]);
    for(j=0;(j<=i)&&(j<=dz);j++)
    {
      p2=gmul((GEN)p4[j+2],(GEN)y[i-j+2]);
      p1=gsub(p1,p2);
    }
    f=gcmp0(p1);
  }
  if(f)
  {
    avma=l1;z=cgetg(3,10);z[1]=evallgef(2)+evalvarn(vx);z[2]=zero;
  }
  else
  {
    z=cgetg(i+4,10);z[1]=evalsigne(1)+evallgef(4+i)+evalvarn(vx);
    z[i+3]=(long)p1;
    for(k=i;k>=0;--k)
    {
      p1=((GEN)x[k+2]);
      for(j=0;(j<=k)&&(j<=dz);j++)
      {
	p2=gmul((GEN)p4[j+2],(GEN)y[k-j+2]);
	p1=gsub(p1,p2);
      }
      z[k+2]=(long)p1;
    }
    tetpil=avma;z=gerepile(l1,tetpil,gcopy(z));
  }
  return z;
}


GEN
poldivres(GEN x, GEN y, GEN *pr)
{
  long  tx=typ(x),ty=typ(y),dx,dy,dz,i,j,k,f,l,tetpil,vx,vy;
  GEN z,p1,p2,p3;

  if(ty<10) {*pr=gzero;return gdiv(x,y);}
  if(tx<10) {*pr=gcopy(x);return gzero;}
  if(tx!=10) err(poler3);
  vx=varn(x);vy=gvar9(y);
  if(vy>vx)
  {
    z=gdiv(x,y);p1=cgetg(3,10);*pr=p1;p1[1]=2;
    p1[2]=zero;setvarn(p1,vx);
  }
  else
  {
    if(typ(y)!=10) err(poler3);
    if(gcmp0(y)) err(poler6);
    dx=lgef(x)-3;dy=lgef(y)-3;
    if((vx>vy)||(dx<dy))
    {
      z=cgetg(3,10);z[1]=2;z[2]=zero;
      setvarn(z,vy);*pr=gcopy(x);
    }
    else
    {
      dz=dx-dy;
      z=cgetg(dz+3,10);z[1]=evalsigne(1)+evallgef(3+dz);setvarn(z,vx);
      z[dz+2]=ldiv((GEN)x[dx+2],(GEN)y[dy+2]);
      for(i=dx-1;i>=dy;--i)
      {
        l=avma;p1=((GEN)x[i+2]);
        for(j=i-dy+1;(j<=i)&&(j<=dz);j++)
        {
          p2=gmul((GEN)z[j+2],(GEN)y[i-j+2]);
          p1=gsub(p1,p2);
        }
        tetpil=avma;p3=gdiv(p1,(GEN)y[dy+2]);
        z[i-dy+2]=lpile(l,tetpil,p3);
      }
      l=avma;f=1;
      for(i=dy-1;(i>=0)&&f;--i)
      {
        l=avma;p1=((GEN)x[i+2]);
        for(j=0;(j<=i)&&(j<=dz);j++)
        {
          p2=gmul((GEN)z[j+2],(GEN)y[i-j+2]);
          tetpil=avma;p1=gsub(p1,p2);
        } 
	p1=gerepile(l,tetpil,p1);f=gcmp0(p1);
      }
      if(f)
      {
        avma=l;*pr=cgetg(3,10);(*pr)[1]=evallgef(2)+evalvarn(vx);(*pr)[2]=zero;
      }
      else
      {
        *pr=cgetg(i+4,10);(*pr)[1]=evalsigne(1)+evallgef(4+i)+evalvarn(vx);
        (*pr)[i+3]=(long)p1;
        for(k=i;k>=0;--k)
        {
          l=avma;p1=((GEN)x[k+2]);
          for(j=0;(j<=k)&&(j<=dz);j++)
          {
            p2=gmul((GEN)z[j+2],(GEN)y[k-j+2]);
            tetpil=avma;p1=gsub(p1,p2);
          }
          (*pr)[k+2]=lpile(l,tetpil,p1);
        }
      }
    }
  }
  return z;
}




/*********************************************************************/
/*********************************************************************/
/*                                                                   */
/*                         FONCTION BEZOUT                           */
/*                                                                   */
/*********************************************************************/
/*********************************************************************/

GEN
bezoutpolnun(GEN a, GEN b, GEN *u, GEN *v)
{
  GEN d,q,u1,u2,u3,w1,w2,w3,*bof;
  long  av,av1,av2,va,dec;

  if ((typ(a)!=10)||(typ(b)!=10)) err(bezoutpoler);
  if((va=varn(a))!=varn(b)) err(bezoutpoler);
  if (lgef(b)>lgef(a))
  {
    u1=b;b=a;a=u1;bof=u;u=v;v=bof;
  }
  if (signe(b))
  {
    w1=a;w2=b;u1=polun[va];u2=gzero;
    av=avma;
    do
    {
      q=poldivres(w1,w2,&w3);
      u3=gsub(u1,gmul(q,u2));
      u1=u2;u2=u3;w1=w2;w2=w3;
    }
    while(signe(w3));
    u3=gdiv(gsub(w1,gmul(u1,a)),b);av2=avma;
    d=gdiv(w1,w3=leadingterm(w1));u2=gdiv(u1,w3);
    u1=gdiv(u3,w3);av1=avma;dec=lpile(av,av2,0)>>TWOPOTBYTES_IN_LONG;
    *u=adecaler(u2,av2,av1)?u2+dec:u2;*v=adecaler(u1,av2,av1)?u1+dec:u1;
    if(adecaler(d,av2,av1)) d+=dec;
  } 
  else {d=gcopy(a);*u= *v=polun[va];}
  return d;
}

GEN
polinvinexact(GEN x, GEN y)
{
  long av=avma,tetpil,i,dx=lgef(x)-3,dy=lgef(y)-3,lz;
  GEN v,z;

  lz=dx+dy;
  v=cgetg(lz+1,18);for(i=1;i<lz;i++) v[i]=zero;
  v[lz]=un;
  v=gauss(sylvestermatrix(y,x),v);
  z=cgetg(dy+2,10);z[1]=evalsigne(1)+evalvarn(varn(y))+evallgef(dy+2);
  for(i=2;i<dy+2;i++) z[i]=v[lz-i+2];
  normalizepol(&z);
  tetpil=avma;return gerepile(av,tetpil,gcopy(z));
}

GEN   polinvmod(GEN x, GEN y)
{
  long  av,av1,fl,tx,vx=varn(x),vy=varn(y);
  GEN u,v,d,p1;

  while(vx!=vy)
  {
    if(vx>vy) {d=cgetg(3,13);d[1]=(long)polun[vx];d[2]=lcopy(x);return d;}
    if(vx<vy)
    {
      if(lgef(x)!=3) err(talker,"non-invertible polynomial in polinvmod");
      x=(GEN)x[2];vx=gvar(x);
    }
  }
  tx=typ(x);
  if((tx!=10)&&(tx!=13)&&(tx!=14)) err(talker,"incorrect type in polinvmod");
  if(tx==10)
  {
    if(isinexactfield(x)||isinexactfield(y))
      return polinvinexact(x,y);
    else
    {
      av=avma;
      d=subresext(x,y,&u,&v);
      if(gcmp0(d)) err(talker,"non-invertible polynomial in polinvmod");
      fl=((typ(d)==10)&&(varn(d)==vx));
      if(fl)
      {
	if(lgef(d)>3) err(talker,"non-invertible polynomial in polinvmod");
	else d=(GEN)d[2];
      }
      av1=avma;return gerepile(av,av1,gdiv(u,d));
    }
  }
  else
  {
    av=avma;
    p1=gmul((GEN)x[1],polinvmod((GEN)x[2],y));
    av1=avma;return gerepile(av,av1,gmod(p1,y));
  }
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                EVALUATION D'UN POLYNOME REEL                    */
/*                                                                 */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


GEN   poleval(GEN x, GEN y)
{
  long  l,av,tetpil,i,tx;
  GEN z,p1,p2,p3,r,s;
  
  if((tx=typ(x))<10) return gcopy(x);
  if(tx!=10) err(polevaler1);
  l=lgef(x);
  if (l==2) z=gzero;
  else
  {
    if (l==3) z=gcopy((GEN)x[2]);
    else
    {
      av=avma;p1=(GEN)x[l-1];
      if(typ(y)!=6)
      {
	for (i=l-1;i>3;--i)
	  p1=gadd(gmul(p1,y),(GEN)x[i-1]);
	p1=gmul(y,p1);tetpil=avma;
	z=gerepile(av,tetpil,gadd(p1,(GEN)x[2]));
      }
      else
      {
	p2=(GEN)x[l-2];r=gadd((GEN)y[1],(GEN)y[1]);
	s=gnorm(y);
	for(i=l-3;i>=2;--i)
	{
	  p3=gadd(p2,gmul(r,p1));
	  p2=gsub((GEN)x[i],gmul(s,p1));
	  p1=p3;
	}
	p1=gmul(y,p1);tetpil=avma;
	z=gerepile(av,tetpil,gadd(p1,p2));
      }
    }
  }
  return z;
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                         RACINES COMPLEXES                       */
/*                                                                 */
/*        l represente la longueur voulue pour les parties         */
/*            reelles et imaginaires des racines de x              */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


GEN   roots(GEN x, long l)
{
  long  av1=avma,av,i,j,f,g,fr,deg,l0,l1,l2,l3,l4,ln,dec;
  long  exc,expmin,m,deg0,k,ti,h,ii,e,e1,emax,v;
  GEN y,xc,xd0,xd,xdabs,p1,p2,p3,p4,p5,p6,p7,p8;
  GEN p9,p10,p11,p12,p14,p15,pa,pax,pb,pp,pq,ps;
  
  if(typ(x)!=10) err(poler7);
  else
  {
    v=varn(x);deg0=lgef(x)-3;expmin=((2-l)<<TWOPOTBITS_IN_LONG)+12;
    if(!signe(x)) err(poler8);
    y=cgetg(deg0+1,18);if(!deg0) return y;
    for(i=1;i<=deg0;i++)
    {
      p1=cgetg(3,6);p1[1]=lgetr(l);
      p1[2]=lgetr(l);y[i]=(long)p1;
      for(j=3;j<l;j++) ((GEN)p1[2])[j]=((GEN)p1[1])[j]=0;
    }
    g=1;f=1;
    for(i=2;i<=deg0+2;i++)
    {
      ti=typ((GEN)x[i]);
      if(ti==8)
      {
	p2=(GEN)((GEN)((GEN)x[i])[1])[2];
	if(gsigne(p2)>0) g=0;
      }
      else
	if(ti>5) g=0;
    }
    l1=avma;p2=cgetg(3,6);
    p2[1]=lmppi(DEFAULTPREC);
    p2[2]=ldivrs((GEN)p2[1],10);
    p11=cgetg(4,10);p11[1]=evalsigne(1)+evallgef(4);setvarn(p11,v);p11[3]=un;
    p12=cgetg(5,10);p12[1]=evalsigne(1)+evallgef(5);setvarn(p12,v);p12[4]=un;
    for(i=2;(i<=deg0+2)&&(gcmp0((GEN)x[i]));i++)
      gaffsg(0,(GEN)y[i-1]);k=i-2;
    if(k!=deg0)
    {
      if(k)
      {
	j=deg0+3-k;pax=cgetg(j,10);pax[1]=evalsigne(1)+evallgef(j);
	setvarn(pax,v);
	for(i=2;i<j;i++)
	  pax[i]=x[i+k];
      }
      else pax=x;
      xd0=deriv(pax,v);pp=ggcd(pax,xd0);m=1;pa=pax;
      h=isnonscalar(pp);if(h) pq=gdeuc(pax,pp);
      do
      {
	if(h)
	{
	  pa=pp;pb=pq;
	  pp=ggcd(pa,deriv(pa,v));h=isnonscalar(pp);
	  if(h) pq=gdeuc(pa,pp);else pq=pa;
	  ps=gdeuc(pb,pq);
	}
	else ps=pa;
	    /* calcul des racines d'ordre exactement m */
	deg=lgef(ps)-3;
	if(deg)
	{
	  l3=avma;e=gexpo((GEN)ps[deg+2]);emax=e;
	  for(i=2;i<deg+2;i++)
	  {
	    p3=(GEN)(ps[i]);
	    if(!gcmp0(p3))
	    {
	      e1=gexpo(p3);
	      if(e1>emax) emax=e1;
	    }
	  }
	  e=emax-e;if(e<0) e=0;avma=l3;
	  if(ps!=pax) xd0=deriv(ps,v);
	  xdabs=cgetg(deg+2,10);xdabs[1]=xd0[1];
	  for(i=2;i<deg+2;i++)
	  {
	    l3=avma;p3=(GEN)xd0[i];p4=gabs(greal(p3),l);
	    p5=gabs(gimag(p3),l);l4=avma;
	    xdabs[i]=lpile(l3,l4,gadd(p4,p5));
	  }
	  l0=avma;xc=gcopy(ps);xd=gcopy(xd0);l2=avma;
	  for(i=1;i<=deg;i++)
	  {
	    if(i==deg)
	    {
	      p1=(GEN)y[k+m*i];
	      gdivz(gneg((GEN)xc[2]),(GEN)xc[3],p1);
	      p14=(GEN)(p1[1]);p15=(GEN)(p1[2]);
	    }
	    else
	    {
	      p3=gshift(p2,e);p4=poleval(xc,p3);
	      p5=gnorm(p4);exc=0;
	      while(exc>= -20)
	      {
		p6=poleval(xd,p3);p7=gneg(gdiv(p4,p6));
		f=1;l3=avma;if(gcmp0(p5)) exc= -32;
		else exc=expo(gnorm(p7))-expo(gnorm(p3));
		avma=l3;
		for(j=1;(j<=10)&&f;j++)
		{
		  p8=gadd(p3,p7);p9=poleval(xc,p8);
		  p10=gnorm(p9);
		  f=(cmprr(p10,p5)>=0)&&(exc>= -20);
		  if(f) {gshiftz(p7,-2,p7);avma=l3;}
		}
		if(f)
		{
		  avma=av1;
		  if(DEBUGLEVEL)
		  {
		    fprintferr("too many iterations in roots(): ");
		    fprintferr("using roots2()\n");flusherr();
		  }
		  return roots2(x,l);
		}
		else
		{
		  av=avma;dec=lpile(l2,l3,0)>>TWOPOTBYTES_IN_LONG;
		  p3=adecaler(p8,l3,av)?p8+dec:p8;
		  p4=adecaler(p9,l3,av)?p9+dec:p9;
		  p5=adecaler(p10,l3,av)?p10+dec:p10;
		}
	      }
	      p1=(GEN)y[k+m*i];setlg((GEN)p1[1],3);
	      setlg((GEN)p1[2],3);gaffect(p3,p1);avma=l2;
	      p14=(GEN)(p1[1]);p15=(GEN)(p1[2]);
	      for(ln=4;ln<=l;ln=(ln<<1)-2)
	      {
		setlg(p14,ln);setlg(p15,ln);
		if(gcmp0(p14))
		{settyp(p14,1);p14[1]=2;}
		if(gcmp0(p15))
		{settyp(p15,1);p15[1]=2;}
		p4=poleval(xc,p1);p5=poleval(xd,p1);
		p6=gneg(gdiv(p4,p5));
		settyp(p14,2);settyp(p15,2);
		gaffect(gadd(p1,p6),p1);avma=l2;
	      }
	    }
	    setlg(p14,l);setlg(p15,l);
	    p7=gcopy(p1);
	    p14=(GEN)(p7[1]);p15=(GEN)(p7[2]);
	    setlg(p14,l+1);setlg(p15,l+1);
	    if(gcmp0(p14))
	    {settyp(p14,1);p14[1]=2;}
	    if(gcmp0(p15))
	    {settyp(p15,1);p15[1]=2;}
	    for(ii=1;ii<=5;ii++)
	    {
	      p4=poleval(ps,p7);p5=poleval(xd0,p7);
	      p6=gneg(gdiv(p4,p5));p7=gadd(p7,p6);
	      p14=(GEN)(p7[1]);p15=(GEN)(p7[2]);
	      if(gcmp0(p14))
	      {settyp(p14,1);p14[1]=2;}
	      if(gcmp0(p15))
	      {settyp(p15,1);p15[1]=2;}
	    }
	    gaffect(p7,p1);p4=poleval(ps,p7);
	    p6=gdiv(p4,poleval(xdabs,gabs(p7,l)));
	    if((!gcmp0(p6))&&(gexpo(p6)>=expmin))
	    {
	      avma=av1;
	      if(DEBUGLEVEL)
	      {
		fprintferr("internal error in roots(): using roots2()\n");flusherr();
	      }
	      return roots2(x,l);
	    }
	    avma=l2;
	    if((expo((GEN)p1[2])<expmin)&&g)
	    {
	      gaffect(gzero,(GEN)p1[2]);
	      for(j=1;j<m;j++)
		gaffect(p1,(GEN)y[k+(i-1)*m+j]);
	      p11[2]=lneg((GEN)p1[1]);l4=avma;
	      xc=gerepile(l0,l4,gdeuc(xc,p11));
	    }
	    else
	    {
	      for(j=1;j<m;j++)
		gaffect(p1,(GEN)y[k+(i-1)*m+j]);
	      if(g)
	      {
		p1=gconj(p1);
		for(j=1;j<=m;j++)
		  gaffect(p1,(GEN)y[k+i*m+j]);i++;
		p12[2]=lnorm(p1);
		p12[3]=lmulsg(-2,(GEN)p1[1]);
		l4=avma;
		xc=gerepile(l0,l4,gdeuc(xc,p12));
	      }
	      else
	      {
		p11[2]=lneg(p1);l4=avma;
		xc=gerepile(l0,l4,gdeuc(xc,p11));
	      }
	    }
	    xd=deriv(xc,v);l2=avma;
	  }
	  k=k+deg*m;
	}
	m++;
      }
      while (k!=deg0);
    }
    avma=l1;
    if(g&&(deg0>1))
    {
      for(j=2;j<=deg0;j++)
      {
	p1=(GEN)y[j];fr=gcmp0((GEN)p1[2]);
	i=j-1;
	if(fr)
	{
	  p2=(GEN)y[i];f=!gcmp0((GEN)p2[2]);
	  if(!f) f=(cmprr((GEN)p2[1],(GEN)p1[1])>0);
	  for(;(i>0)&&f;--i)
	  {
	    y[i+1]=y[i];
	    p2=(GEN)y[i];f=!gcmp0((GEN)p2[2]);
	    if(!f) f=(cmprr((GEN)p2[1],(GEN)p1[1])>0);
	  }
	}
	y[i+1]=(long)p1;
      }
    }
  }
  return y;
}

GEN   rootslong(GEN x, long l)
{
  long  av,av1=avma,i,j,f,g,fr,deg,l0,l1,l2,l3,l4,ln,dec;
  long  exc,expmin,m,deg0,k,ti,h,ii,e,e1,emax,v;
  GEN y,xc,xd0,xd,xdabs,p1,p2,p3,p4,p5,p6,p7,p8;
  GEN p9,p10,p11,p12,p14,p15,pa,pax,pb,pp,pq,ps;
  
  if(typ(x)!=10) err(poler7);
  else
  {
    v=varn(x);deg0=lgef(x)-3;expmin=((2-l)<<TWOPOTBITS_IN_LONG)+12;
    if(!signe(x)) err(poler8);
    y=cgetg(deg0+1,18);if(!deg0) return y;
    for(i=1;i<=deg0;i++)
    {
      p1=cgetg(3,6);p1[1]=lgetr(l);
      p1[2]=lgetr(l);y[i]=(long)p1;
      for(j=3;j<l;j++) ((GEN)p1[2])[j]=((GEN)p1[1])[j]=0;
    }
    g=1;f=1;
    for(i=2;i<=deg0+2;i++)
    {
      ti=typ((GEN)x[i]);
      if(ti==8)
      {
	p2=(GEN)((GEN)((GEN)x[i])[1])[2];
	if(gcmpgs(p2,0)>0) g=0;
      }
      else
	if(ti>5) g=0;
    }
    l1=avma;p2=cgetg(3,6);
    p2[1]=lmppi(l);
    p2[2]=ldivrs((GEN)p2[1],10);
    p11=cgetg(4,10);p11[1]=evalsigne(1)+evallgef(4);setvarn(p11,v);p11[3]=un;
    p12=cgetg(5,10);p12[1]=evalsigne(1)+evallgef(5);setvarn(p12,v);p12[4]=un;
    for(i=2;(i<=deg0+2)&&(gcmp0((GEN)x[i]));i++)
      gaffsg(0,(GEN)y[i-1]);k=i-2;
    if(k!=deg0)
    {
      if(k)
      {
	j=deg0+3-k;pax=cgetg(j,10);pax[1]=evalsigne(1)+evallgef(j);
	setvarn(pax,v);
	for(i=2;i<j;i++)
	  pax[i]=x[i+k];
      }
      else pax=x;
      xd0=deriv(pax,v);pp=ggcd(pax,xd0);m=1;pa=pax;
      h=isnonscalar(pp);if(h) pq=gdeuc(pax,pp);
      do
      {
	if(h)
	{
	  pa=pp;pb=pq;
	  pp=ggcd(pa,deriv(pa,v));h=isnonscalar(pp);
	  if(h) pq=gdeuc(pa,pp);else pq=pa;
	  ps=gdeuc(pb,pq);
	}
	else ps=pa;
	    /* calcul des racines d'ordre exactement m */
	deg=lgef(ps)-3;
	if(deg)
	{
	  l3=avma;e=gexpo((GEN)ps[deg+2]);emax=e;
	  for(i=2;i<deg+2;i++)
	  {
	    p3=(GEN)(ps[i]);
	    if(!gcmp0(p3))
	    {
	      e1=gexpo(p3);
	      if(e1>emax) emax=e1;
	    }
	  }
	  e=emax-e;if(e<0) e=0;avma=l3;
	  if(ps!=pax) xd0=deriv(ps,v);
	  xdabs=cgetg(deg+2,10);xdabs[1]=xd0[1];
	  for(i=2;i<deg+2;i++)
	  {
	    l3=avma;p3=(GEN)xd0[i];p4=gabs(greal(p3),l);
	    p5=gabs(gimag(p3),l);l4=avma;
	    xdabs[i]=lpile(l3,l4,gadd(p4,p5));
	  }
	  l0=avma;xc=gcopy(ps);xd=gcopy(xd0);l2=avma;
	  for(i=1;i<=deg;i++)
	  {
	    if(i==deg)
	    {
	      p1=(GEN)y[k+m*i];
	      gdivz(gneg((GEN)xc[2]),(GEN)xc[3],p1);
	      p14=(GEN)(p1[1]);p15=(GEN)(p1[2]);
	    }
	    else
	    {
	      p3=gshift(p2,e);p4=poleval(xc,p3);
	      p5=gnorm(p4);exc=0;
	      while(exc>= -20)
	      {
		p6=poleval(xd,p3);p7=gneg(gdiv(p4,p6));
		f=1;l3=avma;if(gcmp0(p5)) exc= -32;
		else exc=expo(gnorm(p7))-expo(gnorm(p3));
		avma=l3;
		for(j=1;(j<=50)&&f;j++)
		{
		  p8=gadd(p3,p7);p9=poleval(xc,p8);
		  p10=gnorm(p9);
		  f=(cmprr(p10,p5)>=0)&&(exc>= -20);
		  if(f) {gshiftz(p7,-2,p7);avma=l3;}
		}
		if(f) err(poler9);
		else
		{
		  av=avma;dec=lpile(l2,l3,0)>>TWOPOTBYTES_IN_LONG;
		  p3=adecaler(p8,l3,av)?p8+dec:p8;
		  p4=adecaler(p9,l3,av)?p9+dec:p9;
		  p5=adecaler(p10,l3,av)?p10+dec:p10;
		}
	      }
	      p1=(GEN)y[k+m*i];gaffect(p3,p1);avma=l2;
	      p14=(GEN)(p1[1]);p15=(GEN)(p1[2]);
	      for(ln=4;ln<=l;ln=(ln<<1)-2)
	      {
		if(gcmp0(p14))
		{settyp(p14,1);p14[1]=2;}
		if(gcmp0(p15))
		{settyp(p15,1);p15[1]=2;}
		p4=poleval(xc,p1);p5=poleval(xd,p1);
		p6=gneg(gdiv(p4,p5));
		settyp(p14,2);settyp(p15,2);
		gaffect(gadd(p1,p6),p1);avma=l2;
	      }
	    }
	    p7=gcopy(p1);
	    p14=(GEN)(p7[1]);p15=(GEN)(p7[2]);
	    setlg(p14,l+1);setlg(p15,l+1);
	    if(gcmp0(p14))
	    {settyp(p14,1);p14[1]=2;}
	    if(gcmp0(p15))
	    {settyp(p15,1);p15[1]=2;}
	    for(ii=1;ii<=max(32,((e<<TWOPOTBITS_IN_LONG)+2));ii<<=1)
	    {
	      p4=poleval(ps,p7);p5=poleval(xd0,p7);
	      p6=gneg(gdiv(p4,p5));p7=gadd(p7,p6);
	      p14=(GEN)(p7[1]);p15=(GEN)(p7[2]);
	      if(gcmp0(p14))
	      {settyp(p14,1);p14[1]=2;}
	      if(gcmp0(p15))
	      {settyp(p15,1);p15[1]=2;}
	    }
	    gaffect(p7,p1);p4=poleval(ps,p7);
	    p6=gdiv(p4,poleval(xdabs,gabs(p7,l)));
	    if((!gcmp0(p6))&&(gexpo(p6)>=expmin))
	    {
	      avma=av1;
	      if(DEBUGLEVEL)
	      {
		fprintferr("internal error in roots(): using roots2()\n");flusherr();
	      }
	      return roots2(x,l);
	    }
	    avma=l2;
	    if((expo((GEN)p1[2])<expmin)&&g)
	    {
	      gaffect(gzero,(GEN)p1[2]);
	      for(j=1;j<m;j++)
		gaffect(p1,(GEN)y[k+(i-1)*m+j]);
	      p11[2]=lneg((GEN)p1[1]);l4=avma;
	      xc=gerepile(l0,l4,gdeuc(xc,p11));
	    }
	    else
	    {
	      for(j=1;j<m;j++)
		gaffect(p1,(GEN)y[k+(i-1)*m+j]);
	      if(g)
	      {
		p1=gconj(p1);
		for(j=1;j<=m;j++)
		  gaffect(p1,(GEN)y[k+i*m+j]);i++;
		p12[2]=lnorm(p1);
		p12[3]=lmulsg(-2,(GEN)p1[1]);
		l4=avma;
		xc=gerepile(l0,l4,gdeuc(xc,p12));
	      }
	      else
	      {
		p11[2]=lneg(p1);l4=avma;
		xc=gerepile(l0,l4,gdeuc(xc,p11));
	      }
	    }
	    xd=deriv(xc,v);l2=avma;
	  }
	  k=k+deg*m;
	}
	m++;
      }
      while (k!=deg0);
    }
    avma=l1;
    if(g&&(deg0>1))
    {
      for(j=2;j<=deg0;j++)
      {
	p1=(GEN)y[j];fr=gcmp0((GEN)p1[2]);
	i=j-1;
	if(fr)
	{
	  p2=(GEN)y[i];f=!gcmp0((GEN)p2[2]);
	  if(!f) f=(cmprr((GEN)p2[1],(GEN)p1[1])>0);
	  for(;(i>0)&&f;--i)
	  {
	    y[i+1]=y[i];
	    p2=(GEN)y[i];f=!gcmp0((GEN)p2[2]);
	    if(!f) f=(cmprr((GEN)p2[1],(GEN)p1[1])>0);
	  }
	}
	y[i+1]=(long)p1;
      }
    }
  }
  return y;
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*      Recherche des racines modulo p ( par verif   f(x)=0 )      */
/*                                                                 */
/*      (retourne le vecteur horizontal dont les composantes       */
/*       sont les racines (eventuellement vecteur a 0 comp.)       */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

GEN     rootmod2(GEN f, GEN p)
{
  GEN g,y,z,ss;
  long vf,av,av1,av2,deg,s,nbrac,pasfini;
  
  if((typ(f)!=10)||gcmp0(f)) err(factmoder);
  y=(GEN)(f[2]);vf=varn(f);
  if((typ(y)!=3)||!gegal((GEN)y[1],p))
  {p=gcopy(p);av=avma;f=gmul(f,gmodulcp(gun,p));}
  else 
    av=avma;
  deg=lgef(f)-3;
  s=0;nbrac=0;
  y=cgetg(deg+1,17);
  av1=avma;
  do
  {
    if(av1==avma)
    {
      ss=stoi(s);
      pasfini=(gcmp(p,ss)>0);
    }
    else av1=avma;
    av2=avma;
    z=poleval(f,ss);
    if(gcmp0(z))
    {
      avma=av2;
      nbrac++;
      y[nbrac]=(long)gmodulcp(ss,p);
      f=gdiv(f,gsub(polx[vf],ss));
    }
    else
    {
      avma=av1;s++;
    }
  }
  while((nbrac<deg-1)&&pasfini);
  if (!nbrac) {avma=av;return cgetg(1,17);}
  if (nbrac==(deg-1))
  {
    nbrac++;y[nbrac]=lneg(gdiv((GEN)f[2],(GEN)f[3]));
  }
  g=cgetg(nbrac+1,17);
  for(s=1;s<=nbrac;s++) g[s]=y[s];
  av1=avma;return gerepile(av,av1,gcopy(g));
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*          Recherche intelligente des racines modulo p            */
/*                                                                 */
/*      (retourne le vecteur horizontal dont les composantes       */
/*       sont les racines (eventuellement vecteur a 0 comp.)       */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


GEN rootmod(GEN f, GEN p)
{
  GEN y,unmodp,pol,xun,g,a,b,d,e,u,h,q,p1;
  long av,tetpil,vf,n,i,j,la,flag,lb,rac[10];
  
  if((typ(f)!=10)||gcmp0(f)) err(factmoder);
  y=(GEN)(f[2]);vf=varn(f);
  if((typ(y)!=3)||!gegal((GEN)y[1],p))
  {p=gcopy(p);av=avma;f=gmul(f,gmodulcp(gun,p));}
  else av=avma;
  unmodp=gmodulcp(gun,p);
  if(gegal(p,gdeux))
  {
    j=0;if(gcmp0((GEN)f[2])) j++;
    if(gcmp0(gsubst(f,vf,unmodp))) j+=2;
    avma=av;
    switch(j)
    {
      case 0: y=cgetg(1,17);break;
      case 1: y=cgetg(2,17);y[1]=lmodulcp(gzero,p);break;
      case 2: y=cgetg(2,17);y[1]=lmodulcp(gun,p);break;
      case 3: y=cgetg(3,17);y[1]=lmodulcp(gzero,p);
	y[2]=lmodulcp(gun,p);
    }
    return y;
  }
  if(!gcmpgs(p,4))
  {
    j=0;if(gcmp0((GEN)f[2])) {j++;rac[j]=0;}
    p1=unmodp;
    for(i=1;i<=3;i++)
    {
      if(gcmp0(gsubst(f,vf,p1))) {j++;rac[j]=i;}
      if(i<3) p1=gadd(p1,unmodp);
    }
    avma=av;y=cgetg(j+1,17);
    for(i=1;i<=j;i++) y[i]=lmodulcp(stoi(rac[i]),p);
    return y;
  }
  pol=gmul(unmodp,polx[vf]);
  xun=gmodulcp(pol,f);g=ggcd((GEN)(gsub(gpui(xun,p,0),xun))[2],f);
  n=lgef(g)-3;if(!n) {avma=av;y=cgetg(1,17);return y;}
  y=cgetg(n+1,17);
  if(gcmp0((GEN)g[2]))
  {
    y[1]=zero;g=gdiv(g,polx[vf]);
    if(lgef(g)>3) y[2]=(long)g;
    j=2;
  }
  else {y[1]=(long)g;j=1;}
  while(j<=n)
  {
    a=(GEN)y[j];la=lgef(a)-3;
    if(la==1)
    {
      y[j]=(gneg(gdiv((GEN)a[2],(GEN)a[3])))[2];j++;
    }
    else if(la==2)
    {
      d=gsub(gmul((GEN)a[3],(GEN)a[3]),gmul2n(gmul((GEN)a[2],(GEN)a[4]),2));
      e=gsqrt(d,0);u=gdiv(gun,gmul2n((GEN)a[4],1));
      y[j]=(gmul(u,gsub(e,(GEN)a[3])))[2];y[j+1]=(gmul(u,gneg(gadd(e,(GEN)a[3]))))[2];
      j+=2;
    }
    else
    {
      flag=1;h=pol;q=shifti(subis(p,1),-1);
      while(flag)
      {
	p1=gmodulcp(h,a);b=ggcd((GEN)(gsub(gpui(p1,q,0),gun))[2],a);
	lb=lgef(b)-3;
	if(lb&&(lb<la))
	{
	  flag=0;y[j]=(long)b;y[j+lb]=ldiv(a,b);
	}
	else h=gadd(h,unmodp);
      }
    }
  }
  y=sort(y);tetpil=avma;return gerepile(av,tetpil,gmul(unmodp,y));
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                     FACTORISATION MODULO p                      */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

GEN     factmod2(GEN f, GEN p)
{
  long  i,j,k,d,e,vf,expos[100],psim,nbfact,nbf,av,pk,tetpil,calc;
  GEN   y,t[100],f1,f2,f3,df1,df2,g,g1,g2;
  GEN   xmod,u,v,pd,q;
  
  if((typ(f)!=10)||gcmp0(f)) err(factmoder);
  if((lgef(p)>3)||((lgef(p)==3)&&(cmpis(p,VERYBIGINT)>0))) 
    err(impl,"factmod2 for primes >2^31");
  av=avma;vf=varn(f);p=gcopy(p);f=gmul(f,gmodulcp(gun,p));
  if(lgef(f)==3)
  {avma=av;y=cgetg(3,19);y[1]=lgetg(1,18);y[2]=lgetg(1,18);return y;}
  
      /*  1ere etape : trouver les facteurs square-free produit
	  des facteurs premiers de meme multiplicite    */
  
  f1=f;e=0;nbfact=1;pk=1;df1=deriv(f1,vf);calc=0;
  
  do
  {
/*printf("debut 1re etape\n");*/
    e+=pk;
    while (gcmp0(df1))
    {
      pk *= (psim=itos(p));e=pk;
      j=(lgef(f1)-3)/psim+3;f2=cgetg(j,10);
      f2[1]=evalsigne(1)+evallgef(j);setvarn(f2,vf);
      for(i=2;i<j;i++)
	f2[i]=lcopy((GEN)f1[psim*(i-2)+2]);
      f1=f2;df1=deriv(f1,vf);calc=0;
    }
    if(calc) f2=f3;else f2=ggcd(f1,df1);calc=1;
    if (lgef(f2)<4) u=f1;
    else
    {
      g1=gdiv(f1,f2);
      if (gcmp0(df2=deriv(f2,vf)))
      {
	u=g1;f3=f2;
      }
      else
      {
	f3=ggcd(f2,df2);
	if (lgef(f3)<4) u=gdiv(g1,f2);
	else
	{
	  g2=gdiv(f2,f3);u=gdiv(g1,g2);
	}
      }
    }
    
	/*  2eme etape : 
	    Ici u est un polynome square-free
	    (produit des facteurs premiers de meme multiplicite  e :
	    trouver les facteurs (square-free) produit des facteurs de meme
	    degre d */
    
    d=0;pd=gun;
    xmod=gmodulcp(polx[vf],u);v=xmod;
    while(d<(lgef(u)-3)>>1)
    {
				/*printf("debut 2me etape\n");*/
      d++;
      pd=mulii(pd,p);
      q=shifti(subis(pd,1),-1);
      v=gpui(v,p,0);
      g=ggcd((GEN)(gsub(v,xmod))[2],u);
      
      if (lgef(g)>3)
      {
				/*printf("debut 3me etape\n");*/
	
/*  3eme etape :
    Ici g est produit de pol irreductibles ayant tous le meme degre d;*/
	
	t[nbfact]=g;j=nbfact+(lgef(g)-3)/d;
	split(p[2],t+nbfact,d,p[2],q);
/* le premier parametre est un entier variable m qui sera converti en un
   polynome w dont les coeff sont ses digits en base p (initialement m = p --> X)
   pour faire pgcd de g avec w^(p^d-1)/2 jusqu'a casser. */
	for(;nbfact<j;expos[nbfact++]=e);
	u=gdiv(u,g);v=gmodulcp((GEN)v[2],u);
      }
				/*printf("fin 3me etape\n");*/
      
    }                                 /*printf("fin 2me etape\n");*/
    if (lgef(u)>3) {t[nbfact]=u;expos[nbfact++]=e;}
    f1=f2;df1=df2;
  }
				/*printf("fin 1re etape\n");*/
  while(lgef(f1)>3);
  
				/*printf("fin de l'algorithme\n");*/
  nbf=nbfact;
  tetpil=avma;
  t[1]=gdiv((GEN)t[1],(GEN)((GEN)t[1])[lgef(t[1])-1]);
  for(j=2;j<nbfact;j++)
  {
    if (expos[j]) t[j]=gdiv((GEN)t[j],(GEN)((GEN)t[j])[lgef(t[j])-1]);
    for (k=1;k<j;k++)
      if (expos[k]&&polegal(t[j],t[k])) {expos[k]+=expos[j];expos[j]=0;nbf--;k=j;}
  }
  y=cgetg(3,19);
  u=cgetg(nbf,18);y[1]=(long)u;
  v=cgetg(nbf,18);y[2]=(long)v;
  for(j=1,k=0;j<nbfact;j++)
  {
    if (expos[j])
    {
      k++;
      u[k]=(long)t[j];
      v[k]=lstoi(expos[j]);
    }
  }
  return gerepile(av,tetpil,y);
}

GEN     factcantor(GEN f, GEN p)
{
  long  i,j,k,d,e,vf,expos[100],psim,nbfact,av,tetpil;
  GEN   y,t[100],f1,f2,f3,df1,g,g1;
  GEN   xmod,u,v,pd,q;
  
  if((typ(f)!=10)||gcmp0(f)||(typ(p)!=1)) err(factmoder);
  if((lgef(p)>3)||((lgef(p)==3)&&(cmpis(p,VERYBIGINT)>0))) 
    return factmod_gen(f,p);
  else
  {
    av=avma;vf=varn(f);p=gcopy(p);f=gmul(f,gmodulcp(gun,p));
    if(lgef(f)==3)
    {avma=av;y=cgetg(3,19);y[1]=lgetg(1,18);y[2]=lgetg(1,18);return y;}
    
	/*  1ere etape : trouver les facteurs square-free produit
	    des facteurs premiers de meme multiplicite    */
    
    f1=f;e=1;nbfact=1;
    psim=itos(p);
    while(lgef(f1)>3)
    {
	  /*printf("debut 1re etape\n");*/
      df1=deriv(f1,vf);f2=ggcd(f1,df1);g1=gdiv(f1,f2);k=0;
      while(lgef(g1)>3)
      {
	k++;if(!(k%psim)) {k++;f2=gdiv(f2,g1);}
	f3=ggcd(f2,g1);u=gdiv(g1,f3);g1=f3;f2=gdiv(f2,g1);
	if(lgef(u)>3)
	{
	      /*  2eme etape : 
		  Ici u est un polynome square-free
		  (produit des facteurs premiers de meme multiplicite  e*k :
		  trouver les facteurs (square-free) produit des facteurs de meme
		  degre d */
	  
	  d=0;pd=gun;
	  xmod=gmodulcp(polx[vf],u);v=xmod;
	  while(d<(lgef(u)-3)>>1)
	  {
		/*printf("debut 2me etape\n");*/
	    d++;
	    pd=mulii(pd,p);
	    q=shifti(subis(pd,1),-1);
	    v=gpui(v,p,0);
	    g=ggcd((GEN)(gsub(v,xmod))[2],u);
	    
	    if (lgef(g)>3)
	    {
		  /*printf("debut 3me etape\n");*/
	      
		  /*  3eme etape :
		      Ici g est produit de pol irreductibles ayant tous le meme degre d;*/
	      
	      t[nbfact]=g;j=nbfact+(lgef(g)-3)/d;
	      split(p[2],t+nbfact,d,p[2],q);
		  /* le premier parametre est un entier variable m qui sera converti en un
		     polynome w dont les coeff sont ses digits en base p (initialement m = p --> X)
		     pour faire pgcd de g avec w^(p^d-1)/2 jusqu'a casser. */
	      for(;nbfact<j;expos[nbfact++]=e*k);
	      u=gdiv(u,g);v=gmodulcp((GEN)v[2],u);
	    }
		/*printf("fin 3me etape\n");*/
	    
	  }                         /*printf("fin 2me etape\n");*/
	  if (lgef(u)>3) {t[nbfact]=u;expos[nbfact++]=e*k;}
	}
      }
      e*=psim;j=(lgef(f2)-3)/psim+3;f1=cgetg(j,10);
      f1[1]=evalsigne(1)+evallgef(j);setvarn(f1,vf);
      for(i=2;i<j;i++) f1[i]=lcopy((GEN)f2[psim*(i-2)+2]);
    }
    
	/*printf("fin 1re etape\n");*/
    
	/*printf("fin de l'algorithme\n");*/
    tetpil=avma;
    t[1]=gdiv((GEN)t[1],(GEN)((GEN)t[1])[lgef(t[1])-1]);
    for(j=2;j<nbfact;j++)
      if (expos[j]) t[j]=gdiv((GEN)t[j],(GEN)((GEN)t[j])[lgef(t[j])-1]);
    y=cgetg(3,19);
    u=cgetg(nbfact,18);y[1]=(long)u;
    v=cgetg(nbfact,18);y[2]=(long)v;
    for(j=1,k=0;j<nbfact;j++)
    {
      if (expos[j]) {k++;u[k]=(long)t[j];v[k]=lstoi(expos[j]);}
    }
    return gerepile(av,tetpil,y);
  }
}

GEN  factmod_gen(GEN f, GEN p)
{
  long  j,k,d,e,vf,expos[100],nbfact,av,tetpil;
  GEN   y,t[100],f1,f2,f3,df1,g,g1;
  GEN   xmod,u,v,pd,q;
  
  av=avma;vf=varn(f);p=gcopy(p);f=gmul(f,gmodulcp(gun,p));
  if(lgef(f)==3)
  {avma=av;y=cgetg(3,19);y[1]=lgetg(1,18);y[2]=lgetg(1,18);return y;}
  
      /*  1ere etape : trouver les facteurs square-free produit
	  des facteurs premiers de meme multiplicite    */
  
  f1=f;e=1;nbfact=1;
  
  while(lgef(f1)>3)
  { 
    df1=deriv(f1,vf);f2=ggcd(f1,df1);g1=gdiv(f1,f2);k=0;
    while(lgef(g1)>3)
    {
      k++;
      f3=ggcd(f2,g1);u=gdiv(g1,f3);g1=f3;f2=gdiv(f2,g1);
      if(lgef(u)>3)
      {   
	    /*  2eme etape : 
		Ici u est un polynome square-free
		(produit des facteurs premiers de meme multiplicite  e*k :
		trouver les facteurs (square-free) produit des facteurs de meme
		degre d */
	
	d=0;pd=gun;
	xmod=gmodulcp(polx[vf],u);v=xmod;
	while(d<(lgef(u)-3)>>1) 
	{
	  d++;
	  pd=mulii(pd,p);
	  q=shifti(subis(pd,1),-1);
	  v=gpui(v,p,0);
	  g=ggcd((GEN)(gsub(v,xmod))[2],u);
	  if (lgef(g)>3)
	  {
		/*  3eme etape :Ici g est produit de pol 
		    irreductibles ayant tous le meme degre d;*/
	    
	    t[nbfact]=g;
	    j=nbfact+(lgef(g)-3)/d;
	    splitgen(p,t+nbfact,d,p,q);
/* le premier parametre est un entier variable m qui sera converti en un
   polynome w dont les coeff sont ses digits en base p (initialement m = p --> X
   pour faire pgcd de g avec w^(p^d-1)/2 jusqu'a casser. */
	    
	    for(;nbfact<j;expos[nbfact++]=e*k);
	    u=gdiv(u,g);v=gmodulcp((GEN)v[2],u);
	  }
	}         
	if (lgef(u)>3) {t[nbfact]=u;expos[nbfact++]=e*k;}
      }
    }
    f1=gmul(polun[vf],(GEN)f2[2]);      
  }
  tetpil=avma;
  t[1]=gdiv((GEN)t[1],(GEN)((GEN)t[1])[lgef(t[1])-1]);
  for(j=2;j<nbfact;j++)
    if (expos[j]) t[j]=gdiv((GEN)t[j],(GEN)((GEN)t[j])[lgef(t[j])-1]);
  y=cgetg(3,19);
  u=cgetg(nbfact,18);y[1]=(long)u;
  v=cgetg(nbfact,18);y[2]=(long)v;
  for(j=1,k=0;j<nbfact;j++)
  {
    if (expos[j]) {k++;u[k]=(long)t[j];v[k]=lstoi(expos[j]);}
  }
  return gerepile(av,tetpil,y);
}

GEN     simplefactmod(GEN f, GEN p)
{
  long  i,j,k,d,e,vf,expos[200],t[200],psim,nbfact,av,tetpil;
  GEN   y,f1,f2,f3,df1,g,g1;
  GEN   xmod,u,v,pd,q;
  
  if((typ(f)!=10)||gcmp0(f)) err(factmoder);
  if((lgef(p)>3)||((lgef(p)==3)&&(cmpis(p,VERYBIGINT)>0))) 
    err(impl,"simplefactmod for primes >2^31");
  av=avma;vf=varn(f);p=gcopy(p);f=gmul(f,gmodulcp(gun,p));
  if(lgef(f)==3)
  {avma=av;y=cgetg(3,19);y[1]=lgetg(1,18);y[2]=lgetg(1,18);return y;}
  
      /*  1ere etape : trouver les facteurs square-free produit
	  des facteurs premiers de meme multiplicite    */
  
  f1=f;e=1;nbfact=1;
  psim=itos(p);
  while(lgef(f1)>3)
  {
/*printf("debut 1re etape\n");*/
    df1=deriv(f1,vf);f2=ggcd(f1,df1);g1=gdiv(f1,f2);k=0;
    while(lgef(g1)>3)
    {
      k++;if(!(k%psim)) {k++;f2=gdiv(f2,g1);}
      f3=ggcd(f2,g1);u=gdiv(g1,f3);g1=f3;f2=gdiv(f2,g1);
      if(lgef(u)>3)
      {
	    /*  2eme etape : 
		Ici u est un polynome square-free
		(produit des facteurs premiers de meme multiplicite  e*k :
		trouver les facteurs (square-free) produit des facteurs de meme
		degre d */
	
	d=0;pd=gun;
	xmod=gmodulcp(polx[vf],u);v=xmod;
	while(d<(lgef(u)-3)>>1)
	{
	      /*printf("debut 2me etape\n");*/
	  d++;
	  pd=mulii(pd,p);
	  q=shifti(subis(pd,1),-1);
	  v=gpui(v,p,0);
	  g=ggcd((GEN)(gsub(v,xmod))[2],u);
	  
	  if (lgef(g)>3)
	  {
		/*printf("debut 3me etape\n");*/
	    
		/*  3eme etape :
		    Ici g est produit de pol irreductibles ayant tous le meme degre d;*/
	    
	    j=nbfact+(lgef(g)-3)/d;
	    for(;nbfact<j;nbfact++) {t[nbfact]=d;expos[nbfact]=e*k;}
	    u=gdiv(u,g);v=gmodulcp((GEN)v[2],u);
	  }
				/*printf("fin 3me etape\n");*/
	  
	}                         /*printf("fin 2me etape\n");*/
	if (lgef(u)>3) {t[nbfact]=lgef(u)-3;expos[nbfact++]=e*k;}
      }
    }
    e*=psim;j=(lgef(f2)-3)/psim+3;f1=cgetg(j,10);
    f1[1]=evalsigne(1)+evallgef(j);setvarn(f1,vf);
    for(i=2;i<j;i++) f1[i]=lcopy((GEN)f2[psim*(i-2)+2]);
  }
				/*printf("fin 1re etape\n");*/
  
				/*printf("fin de l'algorithme\n");*/
  tetpil=avma;
  y=cgetg(3,19);
  u=cgetg(nbfact,18);y[1]=(long)u;
  v=cgetg(nbfact,18);y[2]=(long)v;
  for(j=1;j<nbfact;j++) {u[j]=lstoi(t[j]);v[j]=lstoi(expos[j]);}
  return gerepile(av,tetpil,y);
}


GEN stopoly(long m, long p, long v)            
/* renvoie un polynome de la variable v
   dont les coef sont les digits de m en base p */
{
  GEN y;long l=2,i,c[1000];
  do {c[l++]=m%p;m=m/p;} while(m);
  y=cgetg(l,10);for(i=2;i<l;i++) y[i]=lstoi(c[i]);
  y[1]=evalsigne(1)+evallgef(l);setvarn(y,v);
  return y;
}



GEN stopoly_gen(GEN m, GEN  p, long v)            
/* renvoie un polynome de la variable v
   dont les coef sont les digits de m en base p */
{
  GEN y,c;long l=2,i;
  
  c=cgetg(1001,17);
  do  {c[l++]=(long)gmod(m,p);m=gdivent(m,p);} while(!(gcmp0(m)));
  y=cgetg(l,10);for(i=2;i<l;i++) y[i]=c[i];
  y[1]=evalsigne(1)+evallgef(l);setvarn(y,v);
  return y;
}

void split(long m, GEN *t, long d, long p, GEN q)

/*----------------------------------------------------- 
  Programme recursif :
  Entree:
  m entier arbitraire ( converti en un polynome w )
   p nb premier; q=(p^d-1)/2
   t[0] polynome de degre k*d prod de k fact de deg d.
  Sortie:
   t[0],t[1]...t[k-1] contiennent les k facteurs de g
------------------------------------------------------*/

{  
  long j,l,v,av,av1,dv;
  GEN w,w0,wm,unmodp;
  
  if ((dv=lgef(*t)-3)==d) return;
  v=varn(*t);
  unmodp=gmodulcp(gun,stoi(p));
  do
  {
    av=avma;
    if(p==2)
    {
      w=gmul(gpuigs(polx[v],m-1),unmodp);m+=2;
      for(w0=w,j=1;j<d;j++) w=gmod(gadd(w0,gmul(w,w)),*t);
    }
    else
    {
      w=gmul(stopoly(m++,p,v),unmodp);
      wm=gpui(gmodulcp(w,*t),q,0);w=gsub((GEN)wm[2],unmodp);
    }
    av1=avma;
    w=gerepile(av,av1,ggcd(*t,w));
    l=lgef(w)-3;
  }
  while((!l)||(l==dv));
  t[j=l/d]=gdiv(*t,w);*t=w;
  split(m,t+j,d,p,q);split(m,t,d,p,q);
}

void splitgen(GEN m, GEN *t, long d, GEN  p, GEN q)

/*----------------------------------------------------- 
  Programme recursif :
  Entree:
  m entier arbitraire ( converti en un polynome w )
  p nb premier; q=(p^d-1)/2
  t[0] polynome de degre k*d prod de k fact de deg d.
  Sortie:
  t[0],t[1]...t[k-1] contiennent les k facteurs de g
  ------------------------------------------------------*/

{  
  long j,l,v,av,av1,av3,dv,dec;
  GEN w,wm,unmodp;
  
  if ((dv=lgef(*t)-3)==d) return;
  v=varn(*t);
  unmodp=gmodulcp(gun,p);
  do
  {
    av=avma;
    m=gadd(m,gun);
    w=gmul(stopoly_gen(m,p,v),unmodp);
    wm=gpui(gmodulcp(w,*t),q,0);
    w=gsub((GEN)wm[2],unmodp);
    av1=avma;
    m=gcopy(m);w=ggcd(*t,w);
    av3=avma;dec=lpile(av,av1,0)>>TWOPOTBYTES_IN_LONG;
    if(adecaler(m,av1,av3)) m+=dec;
    if(adecaler(w,av1,av3)) w+=dec;
    l=lgef(w)-3;
  }
  while((!l)||(l==dv));
  t[j=l/d]=gdiv(*t,w);*t=w;
  splitgen(m,t+j,d,p,q);splitgen(m,t,d,p,q);
}


GEN decpol(GEN x, long klim)

/* a modifier. Ecrit principalement pour le programme trace.c */
/* aucune verification */
/* klim=0 habituellement, sauf si l'on ne veut chercher que les facteurs de degre <= klim */

{
  short int pos[200];
  long av=avma,av1,tete,lx,k,kin,i,j,i1,i2,fl,d,nbfact;
  GEN res,p1,p2;
  
  kin=1;res=cgetg(lx=lgef(x)-2,17);nbfact=0;
  p1=roots(x,DEFAULTPREC);d=lg(p1)-1;if(!klim) klim=d;
  do
  {
    fl=1;
    for(k=kin;((k+k)<=d)&&fl&&(k<=klim);k++)
    {
      for(i=0;i<=k;i++) pos[i]=i;
      do
      {
	av1=avma;p2=gzero;j=0;
	for(i1=1;i1<=k;i1++) p2=gadd(p2,(GEN)p1[pos[i1]]);
	if((gexpo(gimag(p2))<-20)&&(gexpo(gsub(p2,ground(p2)))<-20))
	{
	  p2=gun;
	  for(i1=1;i1<=k;i1++) p2=gmul(p2,gsub(polx[0],(GEN)p1[pos[i1]]));p2=ground(p2);
	  if((gcmp0(gimag(p2)))&&(gcmp0(gmod(x,p2)))) 
	  {
	    res[++nbfact]=(long)p2;x=gdiv(x,p2);fl=0;kin=k;p2=cgetg(d-k+1,18);
	    i1=1;i2=1;for(i=1;i<=d;i++)
	    {
	      if((i1<=k)&&(i==pos[i1])) i1++;
	      else p2[i2++]=p1[i];
	    }
	    p1=p2;d-=k;
	  }
	}
	if(fl)
	{
	  avma=av1;pos[k]++;
	  while(pos[k-j]>(d-j))
	  {
	    j++;pos[k-j]++;
	  }
	  for(i=k-j+1;i<=k;i++) pos[i]=i+pos[k-j]-k+j;
	}
      }
      while((j<k)&&fl);
    }
  }
  while(((((k+k)<=d)&&(k<=klim))||(!fl))&&(lgef(x)>3));
  if(lgef(x)>3) res[++nbfact]=(long)x;
  setlg(res,nbfact+1);for(j=nbfact+1;j<lx;j++) res[j]=zero; /*necessaire pour gerepile */
  tete=avma;return gerepile(av,tete,greal(res));
}


GEN factpol2(GEN x, long klim)



/* klim=0 habituellement, sauf si l'on ne veut chercher que les facteurs de degre <= klim */

{
  long av=avma,av2,vv,k,i,j,i1,f,nbfac;
  GEN res,p1,p2,y,d,fa[30],a,ap,t,v,w;
  
  if((typ(x)!=10)||(!signe(x))) err(factpoler1);
  y=cgetg(3,19);if(lgef(x)==3) {y[1]=lgetg(1,18);y[2]=lgetg(1,18);return y;}
  if(lgef(x)==4)
  {
    p1=cgetg(2,18);y[1]=(long)p1;p1[1]=lcopy(x);
    p1=cgetg(2,18);y[2]=(long)p1;p1[1]=un;return y;
  }
  d=content(x);vv=varn(x);a=gdiv(x,d);ap=deriv(a,vv);t=ggcd(a,ap);v=gdiv(a,t);
  w=gdiv(ap,t);j=0;f=1;nbfac=0;
  while(f)
  {
    j++;w=gsub(w,deriv(v,vv));f=signe(w);
    if(f) {res=ggcd(v,w);v=gdiv(v,res);w=gdiv(w,res);}
    else res=v;
    fa[j]=(lgef(res)>3) ? decpol(res,klim) : cgetg(1,18);
    nbfac+=(lg(fa[j])-1);
  }
  av2=avma;y=cgetg(3,19);p1=cgetg(nbfac+1,18);y[1]=(long)p1;
  p2=cgetg(nbfac+1,18);y[2]=(long)p2;
  for(i=1,k=0;i<=j;i++)
    for(i1=1;i1<lg(fa[i]);i1++)
    {
      p1[++k]=lcopy((GEN)fa[i][i1]);p2[k]=lstoi(i);
    }
  return gerepile(av,av2,y);
}

GEN     factmod(GEN f, GEN p)
{
  long  i,j,k,e,vf,expos[100],psim,psim2,N,nbfact,av,tetpil,lb,ld,lf,r,kk;
  GEN   y,t[100],f1,f2,f3,df1,g1;
  GEN   xmod,u,v,w,polb,pold,polu,polt,puix,vker,vran,zmodp,unmodp;
  GEN   p1,p2,Q;
  
  if((typ(f)!=10)||gcmp0(f)||(typ(p)!=1)) err(factmoder);
  if(!cmpis(p,2)) return factcantor(f,p);
  if((lgef(p)>3)||((lgef(p)==3)&&(cmpis(p,VERYBIGINT)>0))) 
    return factmod_gen(f,p);
  else
  {
    av=avma;vf=varn(f);p=gcopy(p);unmodp=gmodulcp(gun,p);f=gmul(f,unmodp);
    zmodp=gmodulcp(gzero,p);
    lf=lgef(f)-3;
    puix=cgetg(lf+1,17);puix[1]=(long)polun[vf];
    for(i=1;i<lf;i++) puix[i+1]=lmul((GEN)puix[i],polx[vf]);
    if(!lf)
    {avma=av;y=cgetg(3,19);y[1]=lgetg(1,18);y[2]=lgetg(1,18);return y;}
    
	/*  1ere etape : trouver les facteurs square-free produit
	    des facteurs premiers de meme multiplicite    */
    
    f1=f;e=1;nbfact=1;
    psim=itos(p);psim2=psim>>1;
    while(lgef(f1)>3)
    {
      df1=deriv(f1,vf);f2=ggcd(f1,df1);g1=gdiv(f1,f2);k=0;
      while(lgef(g1)>3)
      {
	k++;if(!(k%psim)) {k++;f2=gdiv(f2,g1);}
	f3=ggcd(f2,g1);u=gdiv(g1,f3);g1=f3;f2=gdiv(f2,g1);
	N=lgef(u)-3;
	if(N)
	{
/*
  Ici u est un polynome square-free
  (produit des facteurs premiers de meme multiplicite  e*k)
  */
	  xmod=gmodulcp(gmul(unmodp,polx[vf]),u);Q=cgetg(N+1,19);
	  p1=cgetg(N+1,18);Q[1]=(long)p1;p1[1]=(long)unmodp;
	  for(i=2;i<=N;i++) p1[i]=(long)zmodp;
	  v=gpui(xmod,p,0);w=v;
	  for(j=2;j<=N;j++)
	  {
	    p1=cgetg(N+1,18);Q[j]=(long)p1;p2=(GEN)w[2];
	    for(i=1;i<=lgef(p2)-2;i++) p1[i]=p2[i+1];
	    for(;i<=N;i++) p1[i]=(long)zmodp;
	    if(j<N) w=gmul(w,v);
	  }
	  setlg(puix,N+1);vker=gmul(puix,ker(gsub(Q,idmat(N))));
	  r=lg(vker)-1;t[nbfact]=u;kk=1;
	  while(kk<r)
	  {
	    vran=cgetg(r+1,18);
	    for(i=1;i<=r;i++) vran[i]=lmodulcp(stoi(mymyrand()%psim),p);
	    polt=gmul(vker,vran);
	    for(i=1;(i<=kk)&&(kk<r);i++)
	    {
	      polb=t[nbfact+i-1];lb=lgef(polb)-3;
	      if(lb>1)
	      {
		polu=(GEN)gpuigs(gmodulcp(polt,polb),psim2)[2];
		pold=ggcd(polb,gaddgs(polu,-1));
		ld=lgef(pold)-3;
		if((ld>0)&&(ld<lb))
		{
		  t[nbfact+i-1]=pold;kk++;
		  t[nbfact+kk-1]=gdiv(polb,pold);
		}
	      }
	    }
	  }
	  for(i=nbfact;i<nbfact+r;i++) expos[i]=e*k;
	  nbfact+=r;
	}
      }
      e*=psim;j=(lgef(f2)-3)/psim+3;f1=cgetg(j,10);
      f1[1]=evalsigne(1)+evallgef(j);setvarn(f1,vf);
      for(i=2;i<j;i++) f1[i]=lcopy((GEN)f2[psim*(i-2)+2]);
    }
    
	/*printf("fin 1re etape\n");*/
    
	/*printf("fin de l'algorithme\n");*/
    setlg(puix,lf+1);
    tetpil=avma;
    t[1]=gdiv((GEN)t[1],(GEN)((GEN)t[1])[lgef(t[1])-1]);
    for(j=2;j<nbfact;j++)
      if (expos[j]) t[j]=gdiv((GEN)t[j],(GEN)((GEN)t[j])[lgef(t[j])-1]);
    y=cgetg(3,19);
    u=cgetg(nbfact,18);y[1]=(long)u;
    v=cgetg(nbfact,18);y[2]=(long)v;
    for(j=1,k=0;j<nbfact;j++)
    {
      if (expos[j]) {k++;u[k]=(long)t[j];v[k]=lstoi(expos[j]);}
    }
    return gerepile(av,tetpil,y);
  }
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                Recherche de racines  p-adiques                  */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/* a etant un p-adique, retourne le vecteur des racines p-adiques de f
   congrues a a modulo p dans le cas ou on suppose f(a) congru a 0 modulo p
(ou a 4 si p=2). */

#define gmaxval(x,y) (gcmp0(x)?BIGINT:ggval(x,y))

GEN
apprgen(GEN f, GEN a)
{
  GEN fp,fg,p1,p,pro,idiot,idiot2,u,ip,quatre;
  long av=avma,tetpil,v,vv,ps,i,j,k,lu,n,fl2;
  
  if((typ(f)!=10)||(typ(a)!=7)||gcmp0(f)) err(rootper1);
  v=varn(f);fp=deriv(f,v);fg=ggcd(f,fp);
  if(lgef(fg)>3) {f=gdiv(f,fg);fp=deriv(f,v);}
  p=(GEN)a[2];p1=poleval(f,a);if((vv=gmaxval(p1,p))<=0) err(rootper2);
  fl2=gegal(p,gdeux);if((vv==1)&&fl2) err(rootper2);
  if(fl2) quatre=stoi(4);vv=gmaxval(poleval(fp,a),p);
  if(!vv)
  {
    while(!gcmp0(p1)) {a=gsub(a,gdiv(p1,poleval(fp,a)));p1=poleval(f,a);}
    tetpil=avma;pro=cgetg(2,17);pro[1]=lcopy(a);
    return gerepile(av,tetpil,pro);
  }
  n=lgef(f)-3;pro=cgetg(n+1,17);j=0;
  p1=poleval(f,gadd(a,gmul(fl2?quatre:p,polx[v])));
  if(!gcmp0(p1)) p1=gdiv(p1,gpuigs(p,ggval(p1,p)));
  if(gcmpgs(p,VERYBIGINT)>0) err(impl,"apprgen for p>=2^31");
  ps=fl2?4:itos(p);idiot=gsub(a,a);idiot2=fl2 ? ggrandocp(p,2) : ggrandocp(p,1);
  for(i=0;i<ps;i++)
  {
    ip=stoi(i);
    if(gcmp0(poleval(p1,gadd(ip,idiot2))))
    {
      u=apprgen(p1,gadd(idiot,ip));
      lu=lg(u);
      for(k=1;k<lu;k++) {j++;pro[j]=ladd(a,gmul(fl2?quatre:p,(GEN)u[k]));}
    }
  }
  tetpil=avma;p1=cgetg(j+1,17);for(i=1;i<=j;i++) p1[i]=lcopy((GEN)pro[i]);
  return gerepile(av,tetpil,p1);
}

/* Retourne le vecteur des racines p-adiques de f en precision r */

GEN
rootpadic(GEN f, GEN p, long r)
{
  GEN y,fp,fg,pro,yi,p1,rac;
  long v,lx,i,j,k,n,av=avma,tetpil,fl2;

  if((typ(f)!=10)||gcmp0(f)) err(rootper1);
  if(r<=0) err(rootper4);
  v=varn(f);fp=deriv(f,v);fg=ggcd(f,fp);
  if(lgef(fg)>3) {f=gdiv(f,fg);fp=deriv(f,v);}
  fl2=gegal(p,gdeux);rac=(fl2&&(r>=2))? rootmod(f,stoi(4)) : rootmod(f,p);
  lx=lg(rac);
  if(r==1)
  {
    tetpil=avma;y=cgetg(lx,18);
    for(i=1;i<lx;i++)
    {
      p1=(GEN)rac[i];yi=cgetg(5,7);yi[2]=lclone(p);yi[3]=lcopy(p);
      yi[4]=lcopy((GEN)p1[2]);yi[1]=evalprecp(1)+HIGHVALPBIT;y[i]=(long)yi;
    }
    return gerepile(av,tetpil,y);
  }
  n=lgef(f)-3;pro=cgetg(n+1,18);j=0;
  for(i=1;i<lx;i++)
  {
    p1=(GEN)rac[i];yi=cgetg(5,7);yi[2]=lclone(p);
    if(signe((GEN)p1[2]))
    {
      if((!fl2)||mpodd((GEN)p1[2]))
      {
	yi[3]=lpuigs(p,r);yi[4]=p1[2];yi[1]=HIGHVALPBIT+evalprecp(r);
      }
      else
      {
	yi[3]=lpuigs(p,r);yi[4]=un;yi[1]=HIGHVALPBIT+1+evalprecp(r);
      }
    }
    else {yi[3]=un;yi[4]=zero;yi[1]=HIGHVALPBIT+r;}
    p1=apprgen(f,yi);for(k=1;k<lg(p1);k++) pro[++j]=p1[k];
  }
  tetpil=avma;p1=cgetg(j+1,17);for(i=1;i<=j;i++) p1[i]=lcopy((GEN)pro[i]);
  return gerepile(av,tetpil,p1);
}  

GEN
rootpadicfast(GEN f, GEN p, long r, long flall)
{
/* a usage interne. Pas de verifs ni de gestion de pile. On suppose que f est
   un polynome a coeffs dans Z de degre n ayant n racines distinctes mod p, et
   p>2, r>=2. On rend les n racines p-adiques en precision r si flall>0,
   1 seule si flall=0
   */

  long i,e,e1,n;
  GEN fa,fp,y,yi,p1;

  fa=rootmod(f,p);n=flall?lgef(f)-3:1;fp=deriv(f,varn(f));y=cgetg(n+1,17);
  for(i=1;i<=n;i++)
  {
    p1=(GEN)fa[i];yi=cgetg(5,7);yi[2]=lclone(p);
    yi[4]=p1[2];
    if(signe((GEN)p1[2])) {yi[3]=lmul(p,p);yi[1]=HIGHVALPBIT+evalprecp(2);}
    else {yi[3]=un;yi[1]=evalvalp(2);}
    e=1;e1=2;
    do
    {
      yi=gsub(yi,gdiv(poleval(f,yi),poleval(fp,yi)));
      e=e1;if(e<r) {e1<<=1;if(e1>r) e1=r;yi=gprec(yi,e1);}
    }
    while(e<r);
    y[i]=(long)yi;
  }
  return y;
}

/* a appartenant a une extension finie de Q_p, retourne le vecteur des racines
de f congrues a a modulo p dans le cas ou on suppose f(a) congru a 0 modulo p
(ou a 4 si p=2). */

GEN
apprgen9(GEN f, GEN a)
{
  GEN fp,fg,p1,p,pro,idiot,idiot2,u,ip,alpha,t,vecg,quatre;
  long av=avma,tetpil,v,vv,ps,i,j,k,lu,n,precpadique,d,va,flfin,fl2;

  if((typ(f)!=10)||gcmp0(f)) err(rootper1);
  if(typ(a)==7) return apprgen(f,a);
  if((typ(a)!=9)||(typ((GEN)a[2])!=10)) err(rootper1);
  v=varn(f);fp=deriv(f,v);fg=ggcd(f,fp);
  if(lgef(fg)>3) {f=gdiv(f,fg);fp=deriv(f,v);}
  alpha=(GEN)a[2];t=(GEN)a[1];precpadique=BIGINT;va=varn(t);
  for(i=2;i<lgef(alpha);i++)
  {
    pro=(GEN)alpha[i];
    if(typ(pro)==7) 
    {
      precpadique=min(precpadique,(signe((GEN)pro[4])?valp(pro)+precp(pro):valp(pro)));
      p=(GEN)pro[2];
    }
  }
  d=lgef(t)-3;
  for(i=2;i<=d+2;i++)
  {
    pro=(GEN)t[i];
    if(typ(pro)==7) 
    {
      precpadique=min(precpadique,(signe((GEN)pro[4])?valp(pro)+precp(pro):valp(pro)));
      p=(GEN)pro[2];
    }
  }
  if(precpadique==BIGINT) err(rootper1);
  p1=poleval(f,a);if((vv=gmaxval(lift(p1),p))<=0) err(rootper2);
  fl2=gegal(p,gdeux);if((vv==1)&&fl2) err(rootper2);
  if(fl2) quatre=stoi(4);vv=gmaxval(lift(poleval(fp,a)),p);
  if(!vv)
  {
    while(!gcmp0(p1)) {a=gsub(a,gdiv(p1,poleval(fp,a)));p1=poleval(f,a);}
    tetpil=avma;pro=cgetg(2,18);pro[1]=lcopy(a);
    return gerepile(av,tetpil,pro);
  }
  n=lgef(f)-3;pro=cgetg(n+1,18);j=0;
  p1=poleval(f,gadd(a,gmul(fl2?quatre:p,polx[v])));
  if(!gcmp0(p1)) p1=gdiv(p1,gpuigs(p,ggval(p1,p)));
  if(gcmpgs(p,VERYBIGINT)>0) err(impl,"apprgen9 for p>=2^31");
  ps=fl2?4:itos(p);idiot=gmodulcp(ggrandocp(p,precpadique),t);
  idiot2=gmodulcp(fl2 ? ggrandocp(p,2) : ggrandocp(p,1),t);
  vecg=cgetg(d+1,18);for(i=1;i<=d;i++) vecg[i]=zero;
  flfin=1;
  while(flfin)
  {
    ip=gmodulcp(gtopoly(vecg,va),t);
    if(gcmp0(poleval(p1,gadd(ip,idiot2))))
    {
      u=apprgen9(p1,gadd(ip,idiot));
      lu=lg(u);
      for(k=1;k<lu;k++) {j++;pro[j]=ladd(a,gmul(fl2?quatre:p,(GEN)u[k]));}
    }
    i=d;while(i&&(!cmpis((GEN)vecg[i],ps-1))) i--;
    if(i) vecg[i]=laddsi(1,(GEN)vecg[i]);
    else flfin=0;
  }
  tetpil=avma;p1=cgetg(j+1,18);for(i=1;i<=j;i++) p1[i]=lcopy((GEN)pro[i]);
  return gerepile(av,tetpil,p1);
}

/* GEN
   padicff(GEN f, GEN p, long r)
             
 interne donc aucune verification 
 En entree, res est un polynome a coeffs dans Z_p sans facteur carre 

{
  long av=avma,tetpil,lx,n,i,j,k,v,fl2;
  GEN rac,pro,p1,p2,y,yi,quatre;

  fl2=gegal(p,gdeux);
  if(fl2&&(r>=2)) err(impl,"padicfactor for p=2");
  rac=factmod(f,p);lx=lg((GEN)rac[1]);
  if(r==1)
    {
      tetpil=avma;y=cgetg(lx,18);p2=(GEN)rac[1];
      for(i=1;i<lx;i++) y[i]=(long)gcvtop((GEN)p2[i],p,1);
      return gerepile(av,tetpil,y);
    }
  if(fl2) quatre=stoi(4);
  n=lgef(f)-3;pro=cgetg(n+1,18);j=0;v=varn(f);p2=lift((GEN)rac[1]);
  for(i=1;i<lx;i++)
    {
      p1=(GEN)p2[i];
      if(lgef(p1)==4)
	{
	  p1=gcmp0((GEN)p1[2])?(GEN)p1[2]:gsub(fl2?quatre:p,(GEN)p1[2]);
	  yi=cgetg(5,7);yi[2]=lclone(p);
	  if(signe(p1))
	    {
	      if((!fl2)||mpodd(p1))
		{
		  yi[3]=lpuigs(p,r);yi[4]=(long)p1;yi[1]=HIGHVALPBIT+evalprecp(r);
		}
	      else
		{
		  yi[3]=lpuigs(p,r);yi[4]=un;yi[1]=HIGHVALPBIT+1+evalprecp(r);
		}
	    }
	  else {yi[3]=un;yi[4]=zero;yi[1]=HIGHVALPBIT+r;}
	  p1=apprgen(f,yi);
	  for(k=1;k<lg(p1);k++) pro[++j]=lsub(polx[v],(GEN)p1[k]);
	}
      else
	{
	  yi=gmodulcp(gcvtop(polx[v],p,r),gcvtop((GEN)p2[i],p,r));
	  p1=apprgen9(f,yi);
	  for(k=1;k<lg(p1);k++) pro[++j]=(long)caract((GEN)p1[k],v);
	}
    }
  tetpil=avma;p1=cgetg(j+1,18);for(i=1;i<=j;i++) p1[i]=lcopy((GEN)pro[i]);
  return gerepile(av,tetpil,p1);
}  
*/


/*****************************************/
/*  Factorisation p-adique d'un polynome */
/*****************************************/

GEN padicff2(GEN nf,GEN p,long pr)

/* Recoit toutes les donnees de
   nf (plus p et pr) et factorise le
 polynome T=nf[1] dans Zp avec la precision pr */

{ 
  long N=lgef((GEN)nf[1])-3,i,j,k,dim_pke,dim_a,av=avma,lgprdec,tetpil;
  GEN mat,mat_pe;
  GEN vecteur,pk,ident,dec_p,mat_pke,mat_a,theta,mat_theta;
  GEN mat_smith1,mat_smith2,facteur;
 
  pk=gpuigs(p,pr);ident=idmat(N);
  dec_p=primedec(nf,p);lgprdec=lg(dec_p)-1;
  facteur=cgetg(lgprdec+1,18);
  for(i=1;i<=lgprdec;i++)
  {
    vecteur=element_pow(nf,(GEN)((GEN)(dec_p[i]))[2],(GEN)((GEN)(dec_p[i]))[3]);
    mat_pe=cgetg(N+1,19);
    for(j=1;j<=N;j++) mat_pe[j]=(long)element_mul(nf,(GEN)ident[j],vecteur);
    mat_pe=hnfmodid(concat(mat_pe,gmul(p,ident)),p);
    mat_pke=mat_pe;
    for(j=1;j<=pr-1;j++) mat_pke=idealmul(nf,mat_pe,mat_pke);
    dim_pke=lg(mat_pke)-1;
    mat_smith1=smith(mat_pke);mat_smith2=(GEN)smith2(mat_pke)[1];
    for(j=1;(j<dim_pke)&&(gegal((GEN)mat_smith1[j],pk));j++);
    if(gegal((GEN)mat_smith1[j],pk)) dim_a=j;
    else dim_a=j-1;
    mat_a=ginv(mat_smith2);
    theta=(GEN)((GEN)nf[8])[2];
    mat=cgetg(dim_a+1,19);
    for(j=1;j<=dim_a;j++)
      mat[j]=(long)element_mul(nf,theta,(GEN)mat_a[j]);
    mat=inverseimage(mat_a,mat);mat_theta=cgetg(dim_a+1,19);
    for(j=1;j<=dim_a;j++) mat_theta[j]=lgetg(dim_a+1,18);
    for(j=1;j<=dim_a;j++)
      for(k=1;k<=dim_a;k++)
	coeff(mat_theta,j,k)=coeff(mat,j,k);
    facteur[i]=(long)caradj(mat_theta,0,0);
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(facteur));
}

GEN
padicff(GEN x,GEN p,long pr)
{
  GEN p1,p2,p3,p5,dx,nf,mat,un_p_adic;
  long N=lgef(x)-3,av=avma,tetpil,i,j;

  nf=cgetg(10,17);nf[1]=(long)x;dx=discsr(x);
  mat=cgetg(3,19);for(i=1;i<=2;i++) mat[i]=lgetg(3,18);
  coeff(mat,1,1)=(long)p;coeff(mat,1,2)=lstoi(ggval(dx,p));
  coeff(mat,2,1)=ldiv(dx,gpuigs(p,ggval(dx,p)));coeff(mat,2,2)=un;
  p3=allbase(x,(long)mat,(GEN*)(nf+3));
  if(!carrecomplet(divii(dx,(GEN)nf[3]),(GEN*)(nf+4))) err(initalgbug1);
  p1=cgetg(N+1,19);
  for(j=1;j<=N;j++)
  {
    p2=cgetg(N+1,18);p1[j]=(long)p2;
    for(i=1;i<=N;i++) p2[i]=(long)truecoeff((GEN)p3[j],i-1);
  }
  p5=cgetg(N*N+1,19);
  for(j=1;j<=N*N;j++) p5[j]=lgetg(N+1,18);
  for(j=1;j<=N*N;j++) 
    for(i=1;i<=N;i++) 
      coeff(p5,i,j)=(long)truecoeff(gmod(gmul((GEN)p3[(j-1)%N+1],(GEN)p3[((j-1)/N)+1]),x),i-1);
  nf[8]=(long)ginv(p1);
  nf[9]=lmul((GEN)nf[8],p5);nf[2]=nf[5]=nf[6]=zero;
  nf[7]=(long)p3;
  un_p_adic=cgetg(5,7);un_p_adic[2]=lclone(p);un_p_adic[3]=lpuigs(p,pr);
  un_p_adic[4]=un;un_p_adic[1]=HIGHVALPBIT+evalprecp(pr);
  p1=padicff2(nf,p,pr);tetpil=avma;
  return gerepile(av,tetpil,gmul(p1,un_p_adic));
}

GEN
factorpadic2(GEN x, GEN p, long r)
{
  long av=avma,av2,vv,k,i,j,i1,f,nbfac;
  GEN res,p1,p2,y,d,fa[100],a,ap,t,v,w;
  
  if((typ(x)!=10)||(!signe(x))) err(rootper1);
  if(r<=0) err(rootper4);
  y=cgetg(3,19);if(lgef(x)==3) {y[1]=lgetg(1,18);y[2]=lgetg(1,18);return y;}
  if(lgef(x)==4)
  {
    p1=cgetg(2,18);y[1]=(long)p1;p1[1]=lcopy(x);
    p1=cgetg(2,18);y[2]=(long)p1;p1[1]=un;return y;
  }
  d=content(x);vv=varn(x);a=gdiv(x,d);ap=deriv(a,vv);t=ggcd(a,ap);v=gdiv(a,t);
  w=gdiv(ap,t);j=0;f=1;nbfac=0;
  while(f)
  {
    j++;w=gsub(w,deriv(v,vv));f=signe(w);
    if(f)
    {
      res=ggcd(v,w);v=gdiv(v,res);w=gdiv(w,res);
    }
    else res=v;
    fa[j]=(lgef(res)>3) ? padicff(res,p,r) : cgetg(1,18);
    nbfac+=(lg(fa[j])-1);
  }
  av2=avma;y=cgetg(3,19);p1=cgetg(nbfac+1,18);y[1]=(long)p1;
  p2=cgetg(nbfac+1,18);y[2]=(long)p2;
  for(i=1,k=0;i<=j;i++)
    for(i1=1;i1<lg(fa[i]);i1++)
    {
      p1[++k]=lcopy((GEN)fa[i][i1]);p2[k]=lstoi(i);
    }
  return gerepile(av,av2,y);
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                  FACTORISATION P-adique avec ROUND 4            */
/*                                                                 */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

GEN
factorpadic4(GEN f,GEN p,long r)
{
  GEN w,g,matf,fx,resint,res,y,ynew,p1,unpadicr;
  long v=varn(f),n=lgef(f)-3,av=avma,tetpil,mfx,nbpoly,i,k,j,m;
  
  if((typ(f)!=10)||(!signe(f))) err(rootper1);
  if(r<=0) err(rootper4);
  
  y=cgetg(3,19);
  if(lgef(f)==3) {y[1]=lgetg(1,18);y[2]=lgetg(1,18);return y;}
  if(lgef(f)==4)
  {
    p1=cgetg(2,18);y[1]=(long)p1;p1[1]=lcopy(f);
    p1=cgetg(2,18);y[2]=(long)p1;p1[1]=un;
    return y;
  }
  j=1;
  matf=squarefree(f);
  nbpoly=lg((GEN)matf[1])-1;
  res=cgetg(3,19);res[1]=lgetg(n+1,18);res[2]=lgetg(n+1,18);
  for(i=1;i<=nbpoly;i++)
  {       
    fx=gcoeff(matf,i,1);
    mfx=ggval(discsr(fx),p);
    m=(r<=mfx)?mfx+1:r;
    w=factmod(fx,p);
    g=bestnu(w);
    if (lg((GEN)w[1])==2)
      resint=nilordpadic(p,m,fx,mfx,g);
    else
      resint=Decomppadic(p,m,fx,mfx,polx[v],fx,g);
    for(k=1;k<=lg((GEN)resint[1])-1;k++)
    {
      coeff(res,j,1)=coeff(resint,k,1);
      coeff(res,j,2)=lmulii(gcoeff(resint,k,2),gcoeff(matf,i,2));
      j++;
    }
  }
  y=cgetg(3,19);y[1]=lgetg(j,18);y[2]=lgetg(j,18);
  for(k=1;k<j;k++)
  {coeff(y,k,1)=coeff(res,k,1);coeff(y,k,2)=coeff(res,k,2);}
  if(r!=m)
  {
    unpadicr=cgetg(5,7);setprecp(unpadicr,r);setvalp(unpadicr,0);
    unpadicr[2]=(long)p;unpadicr[3]=(long)gpuigs(p,r);unpadicr[4]=un;
    tetpil=avma;ynew=cgetg(3,19);ynew[1]=lmul((GEN)y[1],unpadicr);
  }
  else {tetpil=avma;ynew=cgetg(3,19);ynew[1]=lcopy((GEN)y[1]);}
  ynew[2]=lcopy((GEN)y[2]);
  return gerepile(av,tetpil,ynew);
}

GEN
nilordpadic(GEN p,long r,GEN fx,long mf,GEN gx)
{
  
  long Da,Na,La,Ma,first,n,v=varn(fx),av=avma,tetpil; 
  GEN alpha,chi,nu,eta,w,phi,unpadic;
  GEN pmf,Dchi,unmodpmf,res;

  if(DEBUGLEVEL>=3)
  {
    fprintferr(" On entre dans Nilord_padic ");
    if(DEBUGLEVEL>=4)
    {
      fprintferr(" avec les parametres suivants \n ");
      fprintferr(" p=");bruterr(p,'g',-1);
      fprintferr(",  fx=");bruterr(fx,'g',-1);
      fprintferr(",  exposant=%ld,  gx= ",mf);bruterr(gx,'g',-1);
    }
    fprintferr("\n");
  } 
  
  pmf=gpuigs(p,mf+1);
  n=lgef(fx)-3;
  alpha=polx[v];
  first=TRUE;

  unmodpmf=gmodulcp(gun,pmf);
  while (1)
  {
    if (first) 
    {
      chi=fx;
      nu=gx;
      first=FALSE;
    }
    else
    {
      w=factcp(p,fx,alpha);
      chi=(GEN)w[1];
      nu=(GEN)w[2];
      if (cmpis((GEN)w[4],1)==1) 
      {
	tetpil=avma;
	return gerepile(av,tetpil,Decomppadic(p,r,fx,mf,alpha,chi,nu));
      }
    } 
    Da=lgef(nu)-3;
    Na=n/Da;

    if(mf+1<=padicprec(chi,p))
      Dchi=lift(gmul(discsr(lift(gmul(chi,unmodpmf))),unmodpmf));  
    else
      Dchi=discsr(chi);

    if (gcmp0(Dchi))
      Dchi=discsr(chi);
    if (gcmp0(Dchi))
      alpha=gadd(alpha,gmul(p,polx[v]));
    else
    {
      if (gcmp(vstar(p,chi),gzero)==1)
	alpha=gadd(alpha,gun);
      else
      { 	     
	w=setup(p,chi,polx[v],nu);
	eta=(GEN)w[2];
	La=itos((GEN)w[3]);
	Ma=itos((GEN)w[4]);
	if (La>1)
	  alpha=gadd(alpha,eleval(fx,eta,alpha));
	else
	{ 
	  if (Ma==Na)
	  { 
	    unpadic=cgetg(5,7);
	    setprecp(unpadic,r);setvalp(unpadic,0);unpadic[2]=(long)p;
	    unpadic[3]=(long)gpuigs(p,r);unpadic[4]=un;
	    tetpil=avma;
	    res=cgetg(3,19);res[1]=lgetg(2,18);res[2]=lgetg(2,18);
	    coeff(res,1,1)=lmul(fx,unpadic);coeff(res,1,2)=un;
	    if(DEBUGLEVEL>=3)
	    {
	      fprintferr(" On sort de Nilord_padic : Ce cas est fini ");
	      if(DEBUGLEVEL>=4)
	      {
		fprintferr(" avec les parametres suivants \n ");
		fprintferr(" p=");bruterr(p,'g',-1);
		fprintferr(",  fx=");bruterr(fx,'g',-1);
		fprintferr(",  alpha=");bruterr(alpha,'g',-1);
		fprintferr(",  chi=");bruterr(chi,'g',-1);
		
	      }
	      
	      fprintferr("\n");
	    }     
	    return gerepile(av,tetpil,res);
	  }
	  else
	  {
	    w=bsrch(p,chi,ggval(Dchi,p),eta,Ma);
	    phi=eleval(fx,(GEN)w[2],alpha);
	    if (gcmp1((GEN)w[1]))
	    {
	      tetpil=avma;
	      return gerepile(av,tetpil,Decomppadic(p,r,fx,mf,phi,(GEN)w[3],(GEN)w[4]));
	    }
	    else alpha=phi;
	  }
	}
      }
    }
  }
}

GEN
Decomppadic(GEN p,long r,GEN f,long mf,GEN theta,GEN chi,GEN nu)
{
  long n1,n2,j,i,av=avma,tetpil, v=varn(f);
  GEN unmodp,unmodpdrp,unmodpkdr,unmodpr,unmodpdr,unpadic;
  GEN pdr,pk,ph,pr;
  GEN b1,b2,b3,a2,a1,e,f1,f2;
  GEN res;
  long valk;
  
  if(DEBUGLEVEL>=3)
  {
    fprintferr(" On entre dans Decomp_padic ");
    if(DEBUGLEVEL>=4)
    {
      fprintferr(" avec les parametres suivants \n ");
      fprintferr(" p=");bruterr(p,'g',-1);
      fprintferr(" precision=%ld",r);
      fprintferr(",  f=");bruterr(f,'g',-1);
      fprintferr(",  exposant=%ld ",mf);
    }
    fprintferr("\n");
  }

  unmodp=gmodulcp(gun,p);
  pr=gpuigs(p,r);
  unmodpr=gmodulcp(gun,pr);
  
  pdr = (GEN)respm(f, deriv(f,v), gpuigs(p, mf));
  unmodpdr = gmodulcp(gun,pdr);
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
  e=eleval(f,lift(gmul(a1,b2)),theta);

  if(1<=padicprec(e,p))
    e=gdiv(lift(gmul(gmul(pdr,e),unmodpdrp)),pdr);

  pk=p;
  ph=mulii(pdr, pr);
  valk = 1;  
  
 /*    E(t)- e(t) belongs to p^k Op, which is contained in p^(k-df)*Zp[xi]  */

  while (cmpii(pk,ph)==-1)
  {
    e=gmod(gmul(e,gmul(e,gsubsg(3,gmulsg(2,e)))),f);
    pk=gmul(pk,pk);
    valk = 2*valk;  
    unmodpkdr=gmodulcp(gun,mulii(pk,pdr));
    if(valk<=padicprec(e,p))
      e=gdiv(lift(gmul(gmul(pdr,e),unmodpkdr)),pdr);
  }  
  f1=gcdpm(f,gmul(pdr,gsubsg(1,e)),mulii(pr,pdr));
  f1=lift(gmul(gmod(f1,f),unmodpr));
  f2=gdivent(f,f1);                         f2=lift(gmul(gmod(f2,f),unmodpr));

  if(DEBUGLEVEL>=4)
  {
    fprintferr(" Decomp : On considere deux nouveaux polynomes : ");
    fprintferr(" f1=");bruterr(f1,'g',-1);
    fprintferr(",  f2=");bruterr(f2,'g',-1);
    fprintferr("\n");
  }     
  n1=lgef(f1)-3;  b1=factorpadic4(f1,p,r);
  n2=lgef(f2)-3;  b2=factorpadic4(f2,p,r);
  
  unpadic=cgetg(5,7);
  setprecp(unpadic,r);setvalp(unpadic,0);unpadic[2]=(long)p;
  unpadic[3]=(long)gpuigs(p,r);unpadic[4]=un;
  tetpil=avma;
  res=cgetg(3,19);
  for(j=1;j<=2;j++) {res[j]=lgetg(lg((GEN)b1[1])+lg((GEN)b2[1])-1,18);}
  for(i=1;i<=lg((GEN)b1[1])-1;i++)
  {
    coeff(res,i,1)=lmul(gcoeff(b1,i,1),unpadic);
    coeff(res,i,2)=lcopy(gcoeff(b1,i,2));
  }
  for(i=1;i<=lg((GEN)b2[1])-1;i++)
  {
    coeff(res,i+lg((GEN)b1[1])-1,1)=lmul(gcoeff(b2,i,1),unpadic);
    coeff(res,i+lg((GEN)b1[1])-1,2)=lcopy(gcoeff(b2,i,2));
  }
  return gerepile(av,tetpil,res);
}

GEN
squarefree(GEN f)
{
  GEN T,V,W,A,B;
  long n=lgef(f)-3,v=varn(f),i,j,k,av=avma,tetpil;
  
  T=ggcd(deriv(f,v),f);V=gdiv(f,T);
  A=cgetg(3,19);A[1]=lgetg(n+1,18);A[2]=lgetg(n+1,18);
  k=1;i=1;
  do
  {
    W=ggcd(T,V);T=gdiv(T,W); 
    if(lgef(V)!=lgef(W)) {coeff(A,i,1)=ldiv(V,W);coeff(A,i,2)=lstoi(k);i++;}
    k++;V=W;
  }
  while(lgef(V)>3);
  tetpil=avma;
  B=cgetg(3,19);B[1]=lgetg(i,18);B[2]=lgetg(i,18);
  for(j=1;j<=i-1;j++)
  {coeff(B,j,1)=lcopy(gcoeff(A,j,1));coeff(B,j,2)=lcopy(gcoeff(A,j,2));}
  return gerepile(av,tetpil,B);
}

  
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                     FACTORISATION DANS F_q                      */
/*                                                                 */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


GEN
factmod9(GEN f, GEN p, GEN a)
{
  long  i,j,k,d,e,vf,va,expos[100],psim,nbfact,nbf,av,pk,tetpil,calc;
  GEN   y,t[100],f1,f2,f3,df1,df2,g,g1,g2;
  GEN   xmod,u,v,pd,q,qq,unfp,unfq;
  
  if((typ(a)!=10)||(typ(f)!=10)||gcmp0(f)||gcmp0(a)) err(factmoder);
  if(lgef(f)==3)
  {y=cgetg(3,19);y[1]=lgetg(1,18);y[2]=lgetg(1,18);return y;}
  vf=varn(f);va=varn(a);
  if(va<=vf) err(talker,"polynomial variable must be of higher priority than finite field variable\nin factfq");
  av=avma;unfp=gmodulcp(gun,p);
  a=gmul(unfp,a);unfq=gmodulcp(gmul(unfp,polun[va]),a);
  f=gmul(unfq,f);qq=gpuigs(p,lgef(a)-3);
  
      /*  1ere etape : trouver les facteurs square-free produit
	  des facteurs premiers de meme multiplicite    */

  f1=f;e=0;nbfact=1;pk=1;df1=deriv(f1,vf);calc=0;

  do
  {
/*printf("debut 1re etape\n");*/
    e+=pk;
    while (gcmp0(df1))
    {
      pk *= (psim=itos(p));e=pk;
      j=(lgef(f1)-3)/psim+3;f2=cgetg(j,10);
      f2[1]=evalsigne(1)+evallgef(j);setvarn(f2,vf);
      for(i=2;i<j;i++)
	f2[i]=lcopy((GEN)f1[psim*(i-2)+2]);
      f1=f2;df1=deriv(f1,vf);calc=0;
    }
    if(calc) f2=f3;else f2=ggcd(f1,df1);calc=1;
    if (lgef(f2)<4) u=f1;
    else
    {
      g1=gdiv(f1,f2);
      if (gcmp0(df2=deriv(f2,vf)))
      {
	u=g1;f3=f2;
      }
      else
      {
	f3=ggcd(f2,df2);
	if (lgef(f3)<4) u=gdiv(g1,f2);
	else
	{
	  g2=gdiv(f2,f3);u=gdiv(g1,g2);
	}
      }
    }
      
	/*  2eme etape : 
	    Ici u est un polynome square-free
	    (produit des facteurs premiers de meme multiplicite  e :
	    trouver les facteurs (square-free) produit des facteurs de meme
	    degre d */
      
    d=0;pd=gun;
    xmod=gmodulcp(polx[vf],u);v=xmod;
    while(d<(lgef(u)-3)>>1)
    {
				/*printf("debut 2me etape\n");*/
      d++;
      pd=mulii(pd,qq);
      q=shifti(subis(pd,1),-1);
      v=gpui(v,qq,0);
      g=ggcd((GEN)(gsub(v,xmod))[2],u);
	  
      if (lgef(g)>3)
      {
				/*printf("debut 3me etape\n");*/
	      
/*  3eme etape :
    Ici g est produit de pol irreductibles ayant tous le meme degre d;*/

	t[nbfact]=g;j=nbfact+(lgef(g)-3)/d;
	split9(qq,t+nbfact,d,p[2],q,unfq,qq,a);
/* le premier parametre est un entier variable m qui sera converti en un
   polynome w dont les coeff sont ses digits en base p (initialement m = p --> X)
   pour faire pgcd de g avec w^(qq^d-1)/2 jusqu'a casser. */
	for(;nbfact<j;expos[nbfact++]=e);
	u=gdiv(u,g);v=gmodulcp((GEN)v[2],u);
      }
				/*printf("fin 3me etape\n");*/

    }                                 /*printf("fin 2me etape\n");*/
    if (lgef(u)>3) {t[nbfact]=u;expos[nbfact++]=e;}
    f1=f2;df1=df2;
  }
				/*printf("fin 1re etape\n");*/
  while(lgef(f1)>3);
  
				/*printf("fin de l'algorithme\n");*/
  nbf=nbfact;
  tetpil=avma;
  t[1]=gdiv((GEN)t[1],(GEN)((GEN)t[1])[lgef(t[1])-1]);
  for(j=2;j<nbfact;j++)
  {
    if (expos[j]) t[j]=gdiv((GEN)t[j],(GEN)((GEN)t[j])[lgef(t[j])-1]);
    for (k=1;k<j;k++)
      if (expos[k]&&polegal(t[j],t[k])) {expos[k]+=expos[j];expos[j]=0;nbf--;k=j;}
  }
  y=cgetg(3,19);
  u=cgetg(nbf,18);y[1]=(long)u;
  v=cgetg(nbf,18);y[2]=(long)v;
  for(j=1,k=0;j<nbfact;j++)
  {
    if (expos[j])
    {
      k++;
      u[k]=(long)t[j];
      v[k]=lstoi(expos[j]);
    }
  }
  return gerepile(av,tetpil,y);
}

GEN
stopoly9(GEN m, long p, GEN qq, long v, GEN a)

/* renvoie un polynome de la variable v
dont les coef sont les digits de m en base qq */

{
  GEN y,p1,p2,r,c[1000],d[1000];
  long l=2,l1,i,j,va=varn(a);

  do
  {
    m=dvmdii(m,qq,&r);c[l++]=r;
  }
  while(signe(m));
  y=cgetg(l,10);y[1]=evalsigne(1)+evallgef(l);setvarn(y,v);
  for(i=2;i<l;i++) 
  {
    p1=c[i];l1=2;
    do
    {
      p1=dvmdis(p1,p,&r);d[l1++]=r;
    }
    while(signe(p1));
    p2=cgetg(l1,10);p2[1]=evalsigne(1)+evallgef(l1);setvarn(p2,va);
    for(j=2;j<l1;j++) p2[j]=(long)d[j];
    y[i]=lmodulcp(p2,a);
  }
  return y;
}

GEN
stopoly92(long d1, long v, GEN a, GEN *ptres)

/* renvoie un polynome aleatoire de la variable v
de degre inferieur ou egal a 2*d1-1 */

{
  GEN y,p2;
  long p1,l1,i,j,d2,l,va=varn(a),c[1000],d[1000],k=lgef(a)-3,nsh;

/* qqs=2^k */

  c[1]=1;d2=d1+d1+1;
  nsh=BITS_IN_RANDOM-1-k;if(nsh<=0) nsh=1;
  do
  {
    for(l=2;l<=d2;l++) c[l]=(mymyrand())>>nsh;
    l=d2;while(!c[l]) l--;
  }
  while(l<=2);
  l++;
  y=cgetg(l,10);y[1]=evalsigne(1)+evallgef(l);setvarn(y,v);
  for(i=2;i<l;i++) 
  {
    p1=c[i];l1=2;
    do {d[l1++]=p1&1;p1>>=1;} while(p1);
    p2=cgetg(l1,10);p2[1]=evalsigne(1)+evallgef(l1);setvarn(p2,va);
    for(j=2;j<l1;j++) p2[j]=lstoi(d[j]);
    y[i]=lmodulcp(p2,a);
  }
  p1=(mymyrand())>>nsh;l1=2;
  do {d[l1++]=p1&1;p1>>=1;} while(p1);
  p2=cgetg(l1,10);p2[1]=evalsigne(1)+evallgef(l1);setvarn(p2,va);
  for(j=2;j<l1;j++) p2[j]=lstoi(d[j]);
  *ptres=gmodulcp(p2,a);
  return y;
}

void
split9(GEN m, GEN *t, long d, long p, GEN q, GEN unfq, GEN qq, GEN a)

/*----------------------------------------------------- 
   Programme recursif :
  Entree:
   m entier arbitraire ( converti en un polynome w )
   p nb premier; q=(qq^d-1)/2
   t[0] polynome de degre k*d prod de k fact de deg d.
  Sortie:
   t[0],t[1]...t[k-1] contiennent les k facteurs de g
------------------------------------------------------*/

{  
  long j,l,v,av,av1,av2,lim,dec,dv;
  GEN w,w0,wm,res;

  if ((dv=lgef(*t)-3)==d) return;
  v=varn(*t);lim=(bot+avma)>>1;
  do
  {
    av=avma;
    if(p==2)
    {
      w=gmul(stopoly92(d,v,a,&res),unfq);
      for(w0=w,j=1;j<d;j++) w=gmod(gadd(w0,gpui(w,qq,0)),*t);
      w=gsub(w,res);
    }
    else
    {
      w=gmul(stopoly9(m,p,qq,v,a),unfq);m=addsi(1,m);
      wm=gpui(gmodulcp(w,*t),q,0);w=gsub((GEN)wm[2],unfq);
    }
    av2=avma;w=ggcd(*t,w);l=lgef(w)-3;
    if(avma<lim)
    {
      m=gcopy(m);av1=avma;dec=lpile(av,av2,0)>>TWOPOTBYTES_IN_LONG;
      if(adecaler(w,av2,av1)) w+=dec;
      if(adecaler(m,av2,av1)) m+=dec;
    }
  }
  while((!l)||(l==dv));
  t[j=l/d]=gdiv(*t,w);*t=w;
  split9(m,t+j,d,p,q,unfq,qq,a);split9(m,t,d,p,q,unfq,qq,a);
}

GEN
roots2(GEN pol,long PREC)
{
  long av,tetpil,N,flagexactpol,flagrealpol,flagrealrac,ti,i,j,jj;
  long nbpol,k,av1,multiqol,deg,nbroot;
  GEN unp,p1,p2,rr,EPS,qol,qolbis,x,b,c,*ad,v,tabqol;

  av=avma;
  if(typ(pol)!=10) err(poler7);
  if(!signe(pol)) err(poler8);
  N=lgef(pol)-3;
  if(!N) return cgetg(1,18);
  if(N==1)
  {
    unp=cgetr(PREC);affsr(1,unp);p1=gmul(unp,(GEN)pol[3]);p2=gneg(gdiv((GEN)pol[2],p1));
    tetpil=avma;return gerepile(av,tetpil,gcopy(p2));
  }
  EPS=cgetr(3);affsr(1,EPS);flagrealpol=1;flagexactpol=1;
  EPS=gmul2n(EPS,(2-PREC)*BITS_IN_LONG+12);
  for(i=0;i<=N;i++)
  {
    ti=typ((GEN)pol[i+2]);if((ti==2)||(ti>5)) flagexactpol=0;
    if(ti==8){p1=(GEN)((GEN)((GEN)pol[i+2])[1])[2];if(gsigne(p1)>0) flagrealpol=0;}
    else if(ti>5) flagrealpol=0;
  }
  rr=cgetg(N+1,18);
  for(i=1;i<=N;i++)
  {rr[i]=lgetg(3,6);((GEN)(rr[i]))[1]=lgetr(PREC);((GEN)(rr[i]))[2]=lgetr(PREC);}
  tabqol=square_free_factorization(pol);
  nbpol=lg((GEN)tabqol[1])-1;nbroot=0;
  for(k=1;k<=nbpol;k++)
  {
    av1=avma;
    qol=(GEN)((GEN)tabqol[2])[k];qolbis=gcopy(qol);
    multiqol=itos((GEN)((GEN)tabqol[1])[k]);
    deg=lgef(qol)-3;
    for(j=deg;j>=1;j--)
    {
      x=gzero;flagrealrac=0;
      if(j==1) x=gneg(gdiv((GEN)qolbis[2],(GEN)qolbis[3]));
      else x=laguer(qolbis,j,x,EPS,PREC);
      if(flagexactpol)
      {
	x=gprec(x,(long)((PREC-1)*K));
	x=laguer(qol,deg,x,gmul2n(EPS,-32),PREC+1);
      }
      else x=laguer(qol,deg,x,EPS,PREC);
      if((typ(x)==6)&&(gcmp(gabs(gimag(x),PREC),gmul(gdeux,gmul(EPS,gabs(greal(x),PREC))))<=0))
      {x[2]=zero;flagrealrac=1;}
      else if(typ(x)<6) flagrealrac=1;
      for(i=1;i<=multiqol;i++) gaffect(x,(GEN)rr[nbroot+i]);
      nbroot+=multiqol;
      if((!flagrealpol)||flagrealrac)
      {
	ad=(GEN*)newbloc(j+1);
	for(i=0;i<=j;i++) ad[i]=(GEN)qolbis[i+2];
	b=(GEN)ad[j];
	for(jj=j-1;jj>=0;jj--)
	{
	  c=(GEN)ad[jj];
	  ad[jj]=b;
	  b=gadd(gmul((GEN)rr[nbroot],b),c);
	}
	v=cgetg(j+1,17);
	for(i=1;i<=j;i++) v[i]=(long)ad[j-i];
	qolbis=gtopoly(v,varn(qolbis));
	if(flagrealpol)
	  for(i=0;i<=j-1;i++) if(typ((GEN)qolbis[i+2])==6) ((GEN)qolbis[i+2])[2]=zero;
	killbloc((GEN)ad);
      }
      else
      {
	ad=(GEN*)newbloc(j-2);
	ad[j-2]=(GEN)qolbis[j+2];
	p1=gmulsg(2,greal((GEN)rr[nbroot]));p2=gnorm((GEN)rr[nbroot]);
	ad[j-3]=gadd((GEN)qolbis[j+1],gmul(p1,ad[j-2]));
	for(i=j-2;i>=2;i--)
	{ad[i-2]=gadd((GEN)qolbis[i+2],gsub(gmul(p1,ad[i-1]),gmul(p2,ad[i])));}
	v=cgetg(j,17);
	for(i=1;i<=j-1;i++) v[i]=(long)ad[j-1-i];
	qolbis=gtopoly(v,varn(qolbis));
	for(i=0;i<=j-2;i++) if(typ((GEN)qolbis[i+2])==6) ((GEN)qolbis[i+2])[2]=zero;
	for(i=1;i<=multiqol;i++)
	{p1=gconj((GEN)rr[nbroot]);gaffect(p1,(GEN)rr[nbroot+i]);}
	nbroot+=multiqol;
	j--;
	killbloc((GEN)ad);
      }
    }
    avma=av1;
  }
  for(j=2;j<=N;j++)
  {
    x=(GEN)rr[j];
    for(i=j-1;i>=1;i--)
    {
      if(gcmp(greal((GEN)rr[i]),greal(x))<=0) break;
      rr[i+1]=rr[i];
    }
    rr[i+1]=(long)x;
  }
  tetpil=avma;
  return gerepile(av,tetpil,gcopy(rr));
}

#define MR 8
#define MT 10


GEN
laguer(GEN pol,long N,GEN y0,GEN EPS,long PREC)
{
  long av,tetpil,MAXIT,iter,i,j;
  GEN rac,erre,I,x,abx,abp,abm,dx,x1,b,d,f,g,h,sq,gp,gm,g2,*ffrac;
  
  av=avma;
  MAXIT=MR*MT;
  rac=cgetg(3,6);rac[1]=lgetr(PREC);rac[2]=lgetr(PREC);
  I=cgetg(3,6);I[1]=un;I[2]=un;ffrac=(GEN*)newbloc(MR+1);
  for(i=0;i<=MR;i++){ffrac[i]=cgetr(PREC);}
  affrr(dbltor(0.0),ffrac[0]);affrr(dbltor(0.5),ffrac[1]);
  affrr(dbltor(0.25),ffrac[2]);affrr(dbltor(0.75),ffrac[3]);
  affrr(dbltor(0.13),ffrac[4]);affrr(dbltor(0.38),ffrac[5]);
  affrr(dbltor(0.62),ffrac[6]);affrr(dbltor(0.88),ffrac[7]);
  affrr(dbltor(1.0),ffrac[8]);
  x=y0;
  for(iter=1;iter<=MAXIT;iter++)
  {
    b=(GEN)pol[N+2];
    erre=gnorml1(b,PREC);
    d=gzero;
    f=gzero;
    abx=gnorml1(x,PREC);
    for(j=N-1;j>=0;j--)
    {
      f=gadd(gmul(x,f),d);
      d=gadd(gmul(x,d),b);
      b=gadd(gmul(x,b),(GEN)pol[j+2]);
      erre=gadd(gnorml1(b,PREC),gmul(abx,erre));
    }
    erre=gmul(erre,EPS);
    if(gcmp(gnorml1(b,PREC),erre)<=0)
    {
      tetpil=avma;
      killbloc((GEN)ffrac);
      gaffect(x,rac);
      return gerepile(av,tetpil,gcopy(rac));
    }
    g=gdiv(d,b);
    g2=gpuigs(g,2);
    h=gsub(g2,gmul(gdeux,gdiv(f,b)));
    sq=gsqrt(gmul(stoi(N-1),gsub(gmul(stoi(N),h),g2)),PREC);
    gp=gadd(g,sq);
    gm=gsub(g,sq);
    abp=gnorm(gp);
    abm=gnorm(gm);
    if(gcmp(abp,abm)<0) gp=gcopy(gm);
    if(gsigne(max(abp,abm))==1) dx=gdivsg(N,gp);
    else dx=gmul(gadd(gun,abx),gexp(gmulgs(I,iter),PREC));
    x1=gsub(x,dx);
    if(gcmp(gnorml1(gsub(x,x1),PREC),EPS)<0)
    {
      tetpil=avma;
      killbloc((GEN)ffrac);
      gaffect(x,rac);
      return gerepile(av,tetpil,gcopy(rac));
    }
    if(iter%MT) x=gcopy(x1);
    else x=gsub(x,gmul(ffrac[iter/MT],dx));
  }
  avma=av;killbloc((GEN)ffrac);
  err(poler9);
  return gnil;
}

GEN
gnorml1(GEN x,long PREC)
{
  long av,tetpil,lx,i;
  GEN p1,p2,s;
  av=avma;
  switch(typ(x))
  {
    case 1:case 2:case 4: case 5:return gabs(x,PREC);
    case 3:case 7:case 9:case 10:case 11:case 13:case 14:case 15:case 16:return gcopy(x);
    case 6:p1=gabs((GEN)x[1],PREC);p2=gabs((GEN)x[2],PREC);tetpil=avma;
      return gerepile(av,tetpil,gadd(p1,p2));
    case 8:p1=gabs((GEN)x[2],PREC);p2=gabs((GEN)x[3],PREC);tetpil=avma;
      return gerepile(av,tetpil,gadd(p1,p2));
    case 17:case 18:case 19:lx=lg(x);s=gzero;
      for(i=1;i<lx;i++) s=gadd(s,gnorml1((GEN)x[i],PREC));tetpil=avma;
      return gerepile(av,tetpil,gcopy(s));
    default: err(talker,"not a PARI object in gnorml1");
      return gnil;
  }
}

#undef MR
#undef MT

GEN
square_free_factorization(GEN pol)
/* retourne une matrice a deux colonnes: la 1ere contient les i tels que A_i non
   constant, la deuxieme les A_i, telle que pol=A_i1^i1.A_i2^i2...A_in^in.
   Si pol est constant, retourne la matrice vide. */
{
  long av,tetpil,deg,i,j,va,m;
  GEN p1,p2,x,t1,v1,t,v,*A;

  if(typ(pol)!=10) err(poler7);
  deg=lgef(pol)-3;
  if(deg<1) return cgetg(1,19);
  if(deg==1)
  {
    x=cgetg(3,19);x[1]=lgetg(2,18);x[2]=lgetg(2,18);
    p1=(GEN)x[1];p1[1]=un;p2=(GEN)x[2];p2[1]=lcopy(pol);return x;
  }
  av=avma;va=varn(pol);t1=ggcd(pol,deriv(pol,va));
  if(isscalar(t1))
  {
    avma=av;
    x=cgetg(3,19);x[1]=lgetg(2,18);x[2]=lgetg(2,18);
    p1=(GEN)x[1];p1[1]=un;p2=(GEN)x[2];p2[1]=lcopy(pol);return x;
  }
  A=(GEN*)newbloc(deg+1);
  v1=gdeuc(pol,t1);v=v1;i=0;
  while(lgef(v)>3)
  {
    v=ggcd(t1,v1);i++;
    A[i]=gdeuc(v1,v);
    t=gdeuc(t1,v);
    v1=v;t1=t;
  }
  m=0;for(j=1;j<=i;j++) if(isnonscalar(A[j])) m++;
  x=cgetg(3,19);x[1]=lgetg(m+1,18);x[2]=lgetg(m+1,18);m=0;
  for(j=1;j<=i;j++)
  {
    if(isnonscalar(A[j]))
    {
      m++;
      p1=(GEN)x[1];p1[m]=lstoi(j);
      p2=(GEN)x[2];p2[m]=lcopy(A[j]);
    }
  }
  killbloc((GEN)A);
  tetpil=avma;return gerepile(av,tetpil,gcopy(x));
}
