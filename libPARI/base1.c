/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                            */
/*                   CORPS DE NOMBRES                         */
/*                                                            */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

# include "genpari.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~                                                          ~*/
/*~		   HERMITE NORMAL FORM REDUCTION	     ~*/
/*~							     ~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

GEN
hnf(GEN x)
{
  long li,co,av,av0,tetpil,i,j,k,def,ldef,lim;
  GEN p1,p2,y,z,u,v,d,denx;

  if(typ(x)!=19) err(hnfer1);
  lim=(avma+bot)>>1;
  co=lg(x);if(co==1) return gcopy(x);
  li=lg((GEN)x[1]);av0=avma;denx=denom(x);
  av=avma;y=gcmp1(denx)?gcopy(x):gmul(denx,x);
  def=co;ldef=(li>co)?li-co+1:1;
  for(i=li-1;i>=ldef;i--)
  {
    def--;j=def-1;while(j&&(!signe(gcoeff(y,i,j)))) j--;
    while(j>1)
    {
      d=bezout(gcoeff(y,i,j),gcoeff(y,i,j-1),&u,&v);
      if(DEBUGLEVEL>5) {outerr(u);outerr(v);}
      p1=gadd(gmul(u,(GEN)y[j]),gmul(v,(GEN)y[j-1]));
      y[j]=lsub(gmul(divii(gcoeff(y,i,j),d),(GEN)y[j-1]),gmul(divii(gcoeff(y,i,j-1),d),(GEN)y[j]));
      y[j-1]=(long)p1;
      j--;while(j&&(!signe(gcoeff(y,i,j)))) j--;
    }
    if(j==1)
    {
      d=bezout(gcoeff(y,i,1),gcoeff(y,i,def),&u,&v);
      if(DEBUGLEVEL>5) {outerr(u);outerr(v);}
      p1=gadd(gmul(u,(GEN)y[1]),gmul(v,(GEN)y[def]));
      y[1]=lsub(gmul(divii(gcoeff(y,i,1),d),(GEN)y[def]),gmul(divii(gcoeff(y,i,def),d),(GEN)y[1]));
      y[def]=(long)p1;
    }
    p1=gcoeff(y,i,def);
    if(signe(p1)<0) {y[def]=lneg((GEN)y[def]);p1=gcoeff(y,i,def);}
    if(signe(p1))
    {
      for(j=def+1;j<co;j++)
      {
	p2=negi(gdivent(gcoeff(y,i,j),p1));
	y[j]=ladd((GEN)y[j],gmul(p2,(GEN)y[def]));
      }
    }
    else def++;
    if(avma<lim) {tetpil=avma;y=gerepile(av,tetpil,gcopy(y));}
  }
  for(i=0,j=1;j<co;j++) if(!gcmp0((GEN)y[j])) i++;
  tetpil=avma;z=cgetg(i+1,19);
  for(k=0,j=1;j<co;j++) if(!gcmp0((GEN)y[j])) z[++k]=lcopy((GEN)y[j]);
  if(gcmp1(denx)) return gerepile(av0,tetpil,z);
  else {tetpil=avma;return gerepile(av0,tetpil,gdiv(z,denx));}
}

GEN
fasthnf(GEN x,GEN detmat)
{
  long li,co,av,tetpil,i,j,ii,jj,def,lim,jm1;
  GEN p1,p2,y,w,u,v,d,dms2,b;

/* usage interne pas de verification. */
  lim=(avma+bot)>>1;
  av=avma;co=lg(x);li=lg((GEN)x[1]);dms2=shifti(detmat,-1);y=x;
  def=co;
  for(i=li-1;i>=1;i--)
  {
    def--;j=co-li;while(j&&(!signe(gcoeff(y,i,j)))) j--;
    if(j)
    {
      ii=i-1;while(ii&&(!signe(gcoeff(y,ii,def)))) ii--;
      if(!ii)
      {
	p1=gcoeff(y,i,def);
	if(gcmp1(p1)) 
	{
	  for(jj=j;jj;jj--) coeff(y,i,jj)=zero;
	  j=0;
	}
	else
	{
	  for(jj=j;jj;jj--) coeff(y,i,jj)=lmodii(gcoeff(y,i,jj),p1);
	  while(j&&(!signe(gcoeff(y,i,j)))) j--;
	}
      }
    }
    while(j)
    {
      jm1=(j>1)?j-1:def;
      d=bezout(gcoeff(y,i,j),gcoeff(y,i,jm1),&u,&v);
      if(signe(u))
      {
	if(signe(v)) p1=gadd(gmul(u,(GEN)y[j]),gmul(v,(GEN)y[jm1]));
	else p1=gmul(u,(GEN)y[j]);
      }
      else p1=gmul(v,(GEN)y[jm1]);
      y[j]=lsub(gmul(divii(gcoeff(y,i,j),d),(GEN)y[jm1]),gmul(divii(gcoeff(y,i,jm1),d),(GEN)y[j]));
      y[j]=(long)cleanmod((GEN)y[j],i,detmat,dms2);
      y[jm1]=(long)cleanmod(p1,i,detmat,dms2);
      j--;while(j&&(!signe(gcoeff(y,i,j)))) j--;
      if(avma<lim) {tetpil=avma;y=gerepile(av,tetpil,gcopy(y));}
    }
    if(avma<lim) {tetpil=avma;y=gerepile(av,tetpil,gcopy(y));}
  }
  b=detmat;w=cgetg(li,19);def--;
  for(i=li-1;i>=1;i--)
  {
    d=bezout(gcoeff(y,i,i+def),b,&u,&v);w[i]=lmod(gmul(u,(GEN)y[i+def]),b);
    if(!signe(gcoeff(w,i,i))) coeff(w,i,i)=(long)d;
    if(i>1) b=divii(b,d);
  }
  for(i=li-2;i>=1;i--)
  {
    for(j=i+1;j<li;j++)
    {
      p2=gdivent(gcoeff(w,i,j),gcoeff(w,i,i));w[j]=lsub((GEN)w[j],gmul(p2,(GEN)w[i]));
    }
    if(avma<lim) {tetpil=avma;w=gerepile(av,tetpil,gcopy(w));}
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(w));
}

GEN
cleanmod(GEN x,long lim,GEN detmat,GEN detmatsur2)
{
  long lx=lg(x),i;
  GEN y,p1;

  if((lim<=0)||(lim>=lx)) lim=lx-1;
  y=cgetg(lx,18);
  for(i=1;i<=lim;i++)
  {
    p1=modii((GEN)x[i],detmat);
    y[i]=(cmpii(p1,detmatsur2)>0) ? lsubii(p1,detmat) : (long)p1;
  }
  for(;i<lx;i++) y[i]=x[i];
  return y;
}

GEN
allhnfmod(GEN x,GEN detmat,long all)
{
  long li,co,av,tetpil,i,j,jm1,def,ldef,lim;
  GEN b,q,w,p1,y,d,u,v,dms2;

  if(typ(x)!=19) err(hnfer1);
  if(all&&gcmp0(detmat)) return hnf(x);
  lim=(avma+bot)>>1;
  av=avma;co=lg(x);if(co==1) return cgetg(1,19);
  li=lg((GEN)x[1]);dms2=shifti(detmat,-1);y=gcopy(x);
  def=co;ldef=(li>co)?li-co+1:1;
  if(DEBUGLEVEL>6) {fprintferr("entering hnfmod");flusherr();}
  for(i=li-1;i>=ldef;i--)
  {
    if(DEBUGLEVEL>6) {fprintferr("\ni = %ld:",i);flusherr();}
    def--;j=def-1;while(j&&(!signe(gcoeff(y,i,j)))) j--;
    while(j)
    {
      if(DEBUGLEVEL>8) {fprintferr(" %ld",j);flusherr();}
      jm1=(j>1)?j-1:def;
      d=bezout(gcoeff(y,i,j),gcoeff(y,i,jm1),&u,&v);
      if(signe(u))
      {
	if(signe(v)) p1=gadd(gmul(u,(GEN)y[j]),gmul(v,(GEN)y[jm1]));
	else p1=gmul(u,(GEN)y[j]);
      }
      else p1=gmul(v,(GEN)y[jm1]);
      y[j]=lsub(gmul(divii(gcoeff(y,i,j),d),(GEN)y[jm1]),gmul(divii(gcoeff(y,i,jm1),d),(GEN)y[j]));
      y[j]=(long)cleanmod((GEN)y[j],i,detmat,dms2);
      y[jm1]=(long)cleanmod(p1,i,detmat,dms2);
      j--;while(j&&(!signe(gcoeff(y,i,j)))) j--;
      if(avma<lim) {tetpil=avma;y=gerepile(av,tetpil,gcopy(y));}
    }
  }
  if(DEBUGLEVEL>6) {fprintferr("\nCalcul de la diagonale:");flusherr();}
  b=detmat;
  w=cgetg(li,19);def--;
  for(i=li-1;i>=1;i--)
  {
    if(DEBUGLEVEL>6) {fprintferr(" %ld",i);flusherr();}
    d=bezout(gcoeff(y,i,i+def),b,&u,&v);w[i]=lmod(gmul(u,(GEN)y[i+def]),b);
    if(!signe(gcoeff(w,i,i))) coeff(w,i,i)=(long)d;
    if((i>1)&&all) b=divii(b,d);
  }
  if(DEBUGLEVEL>6) {fprintferr("\nNettoyage triangle");flusherr();}  
  for(i=li-2;i>=1;i--)
  {
    for(j=i+1;j<li;j++)
    {
      q=gdivent(gcoeff(w,i,j),gcoeff(w,i,i));w[j]=lsub((GEN)w[j],gmul(q,(GEN)w[i]));
    }
    if(avma<lim) {tetpil=avma;w=gerepile(av,tetpil,gcopy(w));}
  }
  if(DEBUGLEVEL>6) {fprintferr("\nFin hnfmod\n");flusherr();}    
  tetpil=avma;return gerepile(av,tetpil,gcopy(w));
}

GEN
hnfmod(GEN x, GEN detmat)
{
  return allhnfmod(x,detmat,1);
}

GEN
hnfmodid(GEN x, GEN p)
{
  return allhnfmod(x,p,0);
}

GEN
hnfhavas(GEN x)
  /* rend [y,U,V] tel que y=V.x.U , V vecteur de permutation, U matrice 
    unimodulaire, y forme HNF de x */
{
  long li,co,av,av0,av1,tetpil,i,j,k,def,ldef,lim,fl,fl1,imin,jmin,vpk;
  long jpro,com,dec,vi;
  GEN p1,p2,p3,p12,p13,y,z,u,denx,vperm,mat1,col2,lil2,s,pmin,apro,bpro,cpro;

  if(DEBUGLEVEL>6)
  {fprintferr("A l'entree de hnfhavas : AVMA = %ld\n",avma);flusherr();}
  if(typ(x)!=19) err(hnfer1);
  av0=avma;co=lg(x);if(co==1) return gcopy(x);
  li=lg((GEN)x[1]);denx=denom(x);
  vperm=cgeti(li);for(i=1;i<li;i++) vperm[i]=i;
  lim=(avma+bot)>>1;av=avma;u=idmat(co-1);
  y=gcmp1(denx)?gcopy(x):gmul(denx,x);def=co;ldef=(li>co)?li-co+1:1;
  for(i=li-1;i>=ldef;i--)
  {
    def--;av1=avma;mat1=cgetg(def+1,19);col2=cgetg(def+1,18);
    for(j=1;j<=def;j++)
    {
      p1=cgetg(i+1,18);mat1[j]=(long)p1;s=gzero;
      for(k=1;k<=i;k++)
      {p2=gsqr(gcoeff(y,vperm[k],j));p1[k]=(long)p2;s=gadd(s,p2);}
      col2[j]=(long)s;
    }
    lil2=cgetg(i+1,18);
    for(k=1;k<=i;k++)
    {
      s=gzero;
      for(j=1;j<=def;j++) s=gadd(s,gcoeff(mat1,k,j));
      lil2[k]=(long)s;
    }
    fl=1;fl1=1;
    for(k=i;k>=1;k--)
    {
      vpk=vperm[k];
      if(signe((GEN)lil2[k]))
      {
	fl=0;
	if(fl1||(cmpii((GEN)lil2[k],pmin)<=0))
	{
	  j=1;while(!signe(gcoeff(y,vpk,j))) j++;
	  if(fl1)
	  {
	    imin=k;jmin=j;pmin=mulii((GEN)lil2[k],(GEN)col2[j]);
	    cpro=absi(gcoeff(y,vpk,j));fl1=0;
	  }
	  jpro=j;apro=absi(gcoeff(y,vpk,j));j++;
	  for(;j<=def;j++)
	  {
	    if(signe(gcoeff(y,vpk,j))&&(com=cmpii((GEN)col2[j],(GEN)col2[jpro])<=0))
	    {
	      if(com<0) {jpro=j;apro=absi(gcoeff(y,vpk,j));}
	      else
	      {
		bpro=absi(gcoeff(y,vpk,j));
		if(cmpii(bpro,apro)<0) {jpro=j;apro=bpro;}
	      }
	    }
	  }
	  com=cmpii(p1=mulii((GEN)lil2[k],(GEN)col2[jpro]),pmin);
	  if((com<0)||((com==0)&&(cmpii(apro,cpro)<0)))
	  {pmin=p1;imin=k;jmin=jpro;cpro=apro;}
	}
      }
    }
    avma=av1;
    if(fl) goto comterm;
    else
    {
      if(jmin!=def) 
      {
	p1=(GEN)y[def];y[def]=y[jmin];y[jmin]=(long)p1;
	p1=(GEN)u[def];u[def]=u[jmin];u[jmin]=(long)p1;
      }
      if(imin!=i) {fl=vperm[i];vperm[i]=vperm[imin];vperm[imin]=fl;}
      vi=vperm[i];fl=1;
      while(fl)
      {
	if(signe(gcoeff(y,vi,def))<0) 
	{
	  y[def]=lneg((GEN)y[def]);u[def]=lneg((GEN)u[def]);
	}
	p1=gcoeff(y,vi,def);p12=shifti(p1,-1);p13=negi(p12);
	for(j=1;j<def;j++)
	{
	  p2=dvmdii(gcoeff(y,vi,j),p1,&p3);
	  if(cmpii(p3,p13)<0) p2=addis(p2,-1);
	  else {if(cmpii(p3,p12)>0) p2=addis(p2,1);}
	  if(DEBUGLEVEL>5) outerr(p2);
	  y[j]=ladd((GEN)y[j],gmul(p2=negi(p2),(GEN)y[def]));
	  u[j]=ladd((GEN)u[j],gmul(p2,(GEN)u[def]));
	}
	j=1;while(!signe(gcoeff(y,vi,j))) j++;
	if(j<def) 
	{
	  pmin=gnorml2((GEN)y[j]);jmin=j;apro=absi(gcoeff(y,vi,j));
	  j++;
	  for(;j<def;j++)
	  {
	    if(signe(gcoeff(y,vi,j)))
	    {
	      p1=gnorml2((GEN)y[j]);com=cmpii(p1,pmin);
	      if((com<0)||((com==0)&&(cmpii(bpro=absi(gcoeff(y,vi,j)),apro)<0)))
	      {
		pmin=p1;jmin=j;apro=bpro;
	      }
	    }
	  }
	  p1=(GEN)y[def];y[def]=y[jmin];y[jmin]=(long)p1;
	  p1=(GEN)u[def];u[def]=u[jmin];u[jmin]=(long)p1;
	}
	else fl=0;
      }
    }
    vi=vperm[i];p1=gcoeff(y,vi,def);
    for(j=def+1;j<co;j++)
    {
      p2=negi(gdivent(gcoeff(y,vi,j),p1));
      if(DEBUGLEVEL>5) outerr(p2);
      y[j]=ladd((GEN)y[j],gmul(p2,(GEN)y[def]));
      u[j]=ladd((GEN)u[j],gmul(p2,(GEN)u[def]));
    }
    if(avma<lim)
    {
      tetpil=avma;y=gcopy(y);u=gcopy(u);
      dec=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;
      y+=dec;u+=dec;
    }
  }
  comterm:
  tetpil=avma;z=cgetg(4,17);
  fl=gcmp1(denx);p1=cgetg(co,19);
  for(j=1;j<co;j++)
  {
    p2=cgetg(li,18);p1[j]=(long)p2;
    for(i=1;i<li;i++) p2[i]=fl?lcopy(gcoeff(y,vperm[i],j)):ldiv(gcoeff(y,vperm[i],j),denx);
  }
  z[1]=(long)p1;z[2]=lcopy(u);
  p1=cgetg(li,17);for(i=1;i<li;i++) p1[i]=lstoi(vperm[i]);
  z[3]=(long)p1;return gerepile(av0,tetpil,z);
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~							      ~*/
/*~      	    SMITH NORMAL FORM REDUCTION	              ~*/
/*~							      ~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

GEN
smith(GEN x)
{
  long li,av,tetpil,i,j,k,l,lim,c,fl,n;
  GEN p1,p2,p3,p4,y,z,b,u,v,d;

  if(typ(x)!=19) err(hnfer1);
  lim=(avma+bot)>>1;
  av=avma;n=lg(x)-1;if(!n) return cgetg(1,17);
  li=lg((GEN)x[1])-1;y=gcopy(x);
  if(li!=n) err(hnfer2);
  for(i=n;i>=2;i--)
  {
    do
    {
      c=0;
      for(j=i-1;j>=1;j--)
      {
	p1=gcoeff(y,i,j);
	if(signe(p1))
	{
	  p2=gcoeff(y,i,i);
	  if(gegal(p1,p2)) {d=p1;u=gun;v=gzero;p3=gun;p4=gun;}
	  else if(!signe(addii(p1,p2))) 
	  {
	    d=absi(p1);u=(signe(p2)>0)?gun:gneg(gun);
	    v=gzero;p3=u;p4=gneg(u);
	  }
	  else {d=bezout(p2,p1,&u,&v);p3=divii(p2,d);p4=divii(p1,d);}
	  for(k=1;k<=i;k++)
	  {
	    b=addii(mulii(u,gcoeff(y,k,i)),mulii(v,gcoeff(y,k,j)));
	    coeff(y,k,j)=lsubii(mulii(p3,gcoeff(y,k,j)),mulii(p4,gcoeff(y,k,i)));
	    coeff(y,k,i)=(long)b;
	  }
	}
      }
      for(j=i-1;j>=1;j--)
      {
	p1=gcoeff(y,j,i);
	if(signe(p1))
	{
	  p2=gcoeff(y,i,i);
	  if(gegal(p1,p2)) {d=p1;u=gun;v=gzero;p3=gun;p4=gun;}
	  else if(!signe(addii(p1,p2))) 
	  {
	    d=absi(p1);u=(signe(p2)>0)?gun:gneg(gun);
	    v=gzero;p3=u;p4=gneg(u);
	  }
	  else {d=bezout(p2,p1,&u,&v);p3=divii(p2,d);p4=divii(p1,d);}
	  for(k=1;k<=i;k++)
	  {
	    b=addii(mulii(u,gcoeff(y,i,k)),mulii(v,gcoeff(y,j,k)));
	    coeff(y,j,k)=lsubii(mulii(p3,gcoeff(y,j,k)),mulii(p4,gcoeff(y,i,k)));
	    coeff(y,i,k)=(long)b;
	  }
	  c++;
	}
      }
      if(!c)
      {
	b=gcoeff(y,i,i);fl=1;
	if(signe(b))
	{
	  for(k=1;(k<i)&&fl;k++)
	    for(l=1;(l<i)&&fl;l++)
	      fl= !signe(modii(gcoeff(y,k,l),b));
	  if(!fl)
	  {
	    k--;
	    for(l=1;l<=i;l++)
	      coeff(y,i,l)=laddii(gcoeff(y,i,l),gcoeff(y,k,l));
	  }
	}
      }
      if(avma<lim) {tetpil=avma;y=gerepile(av,tetpil,gcopy(y));}
    }
    while(c||(!fl));
  }
  tetpil=avma;z=cgetg(n+1,17);
  for(j=0,k=1;k<=n;k++) if(!signe(gcoeff(y,k,k))) z[++j]=zero;
  for(k=1;k<=n;k++) if(signe(gcoeff(y,k,k))) z[++j]=(long)gabs(gcoeff(y,k,k),0);
  return gerepile(av,tetpil,z);
}

GEN
smith2(GEN x)
{
  long li,av,tetpil,i,j,k,l,lim,c,fl,n,dec;
  GEN p1,p2,p3,p4,y,z,b,u,v,d,ml,mr;

  if(typ(x)!=19) err(hnfer1);
  lim=(avma+bot)>>1;
  av=avma;n=lg(x)-1;
  if(!n) {z=cgetg(3,17);z[1]=lgetg(1,19);z[2]=lgetg(1,19);return z;}
  li=lg((GEN)x[1])-1;y=gcopy(x);
  if(li!=n) err(hnfer2);
  ml=idmat(n);mr=idmat(n);
  for(i=n;i>=2;i--)
  {
    do
    {
      c=0;
      for(j=i-1;j>=1;j--)
      {
	p1=gcoeff(y,i,j);
	if(signe(p1))
	{
	  p2=gcoeff(y,i,i);
	  if(gegal(p1,p2)) {d=p1;u=gun;v=gzero;p3=gun;p4=gun;}
	  else if(!signe(addii(p1,p2))) 
	  {
	    d=absi(p1);u=(signe(p2)>0)?gun:gneg(gun);
	    v=gzero;p3=u;p4=gneg(u);
	  }
	  else {d=bezout(p2,p1,&u,&v);p3=divii(p2,d);p4=divii(p1,d);}
	  for(k=1;k<=i;k++)
	  {
	    b=addii(mulii(u,gcoeff(y,k,i)),mulii(v,gcoeff(y,k,j)));
	    coeff(y,k,j)=lsubii(mulii(p3,gcoeff(y,k,j)),mulii(p4,gcoeff(y,k,i)));
	    coeff(y,k,i)=(long)b;
	  }
	  b=gadd(gmul(u,(GEN)mr[i]),gmul(v,(GEN)mr[j]));
	  mr[j]=lsub(gmul(p3,(GEN)mr[j]),gmul(p4,(GEN)mr[i]));
	  mr[i]=(long)b;
	}
      }
      for(j=i-1;j>=1;j--)
      {
	p1=gcoeff(y,j,i);
	if(signe(p1))
	{
	  p2=gcoeff(y,i,i);
	  if(gegal(p1,p2)) {d=p1;u=gun;v=gzero;p3=gun;p4=gun;}
	  else if(!signe(addii(p1,p2))) 
	  {
	    d=absi(p1);u=(signe(p2)>0)?gun:gneg(gun);
	    v=gzero;p3=u;p4=gneg(u);
	  }
	  else {d=bezout(p2,p1,&u,&v);p3=divii(p2,d);p4=divii(p1,d);}
	  for(k=1;k<=i;k++)
	  {
	    b=addii(mulii(u,gcoeff(y,i,k)),mulii(v,gcoeff(y,j,k)));
	    coeff(y,j,k)=lsubii(mulii(p3,gcoeff(y,j,k)),mulii(p4,gcoeff(y,i,k)));
	    coeff(y,i,k)=(long)b;
	  }
	  b=gadd(gmul(u,(GEN)ml[i]),gmul(v,(GEN)ml[j]));
	  ml[j]=lsub(gmul(p3,(GEN)ml[j]),gmul(p4,(GEN)ml[i]));
	  ml[i]=(long)b;
	  c++;
	}
      }
      if(!c)
      {
	b=gcoeff(y,i,i);fl=1;
	if(signe(b))
	{
	  for(k=1;(k<i)&&fl;k++)
	    for(l=1;(l<i)&&fl;l++)
	      fl= !signe(modii(gcoeff(y,k,l),b));
	  if(!fl)
	  {
	    k--;
	    for(l=1;l<=i;l++)
	      coeff(y,i,l)=laddii(gcoeff(y,i,l),gcoeff(y,k,l));
	    ml[i]=ladd((GEN)ml[i],(GEN)ml[k]);
	  }
	}
      }
      if(avma<lim) 
      {
	tetpil=avma;y=gcopy(y);ml=gcopy(ml);mr=gcopy(mr);
	dec=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;y+=dec;ml+=dec;mr+=dec;
      }
    }
    while(c||(!fl));
  }
  for(k=1;k<=n;k++) if(signe(gcoeff(y,k,k))<0) mr[k]=lneg((GEN)mr[k]);
  ml=gtrans(ml);tetpil=avma;z=cgetg(3,17);
  z[1]=lcopy(ml);z[2]=lcopy(mr);
  return gerepile(av,tetpil,z);
}

GEN
gsmith(GEN x)
{
  long li,av,tetpil,i,j,k,l,lim,c,fl,n;
  GEN p1,p2,p3,p4,y,z,b,u,v,d;

  if(typ(x)!=19) err(hnfer1);
  lim=(avma+bot)>>1;
  av=avma;n=lg(x)-1;if(!n) return cgetg(1,17);
  li=lg((GEN)x[1])-1;y=gcopy(x);
  if(li!=n) err(hnfer2);
  for(i=n;i>=2;i--)
  {
    do
    {
      c=0;
      for(j=i-1;j>=1;j--)
      {
	p1=gcoeff(y,i,j);
	if(signe(p1))
	{
	  p2=gcoeff(y,i,i);v=gdiventres(p1,p2);
	  if(gcmp0((GEN)v[2])) {d=p2;p4=(GEN)v[1];v=gzero;p3=gun;u=gun;}
	  else {d=gbezout(p2,p1,&u,&v);p3=gdiv(p2,d);p4=gdiv(p1,d);}
	  for(k=1;k<=i;k++)
	  {
	    b=gadd(gmul(u,gcoeff(y,k,i)),gmul(v,gcoeff(y,k,j)));
	    coeff(y,k,j)=lsub(gmul(p3,gcoeff(y,k,j)),gmul(p4,gcoeff(y,k,i)));
	    coeff(y,k,i)=(long)b;
	  }
	}
      }
      for(j=i-1;j>=1;j--)
      {
	p1=gcoeff(y,j,i);
	if(signe(p1))
	{
	  p2=gcoeff(y,i,i);v=gdiventres(p1,p2);
	  if(gcmp0((GEN)v[2])) {d=p2;p4=(GEN)v[1];v=gzero;p3=gun;u=gun;}
	  else {d=gbezout(p2,p1,&u,&v);p3=gdiv(p2,d);p4=gdiv(p1,d);}
	  for(k=1;k<=i;k++)
	  {
	    b=gadd(gmul(u,gcoeff(y,i,k)),gmul(v,gcoeff(y,j,k)));
	    coeff(y,j,k)=lsub(gmul(p3,gcoeff(y,j,k)),gmul(p4,gcoeff(y,i,k)));
	    coeff(y,i,k)=(long)b;
	  }
	  c++;
	}
      }
      if(!c)
      {
	b=gcoeff(y,i,i);fl=1;
	if(signe(b))
	{
	  for(k=1;(k<i)&&fl;k++)
	    for(l=1;(l<i)&&fl;l++)
	      fl= !signe(gmod(gcoeff(y,k,l),b));
	  if(!fl)
	  {
	    k--;
	    for(l=1;l<=i;l++)
	      coeff(y,i,l)=ladd(gcoeff(y,i,l),gcoeff(y,k,l));
	  }
	}
      }
      if(avma<lim) {tetpil=avma;y=gerepile(av,tetpil,gcopy(y));}
    }
    while(c||(!fl));
  }
  tetpil=avma;z=cgetg(n+1,17);
  for(j=0,k=1;k<=n;k++) if(!signe(gcoeff(y,k,k))) z[++j]=zero;
  for(k=1;k<=n;k++)
    if(signe(p1=gcoeff(y,k,k)))
    {
      if(typ(p1)==1) z[++j]=(long)gabs(p1,0);
      else
      {
	if(typ(p1)==10)
	{
	  p2=(GEN)p1[lgef(p1)-1];
	  if((typ(p2)==1)&&(signe(p2)<0)) z[++j]=lneg(p1);
	  else z[++j]=lcopy(p1);
	}
	else z[++j]=lcopy(p1);
      }
    }
  return gerepile(av,tetpil,z);
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~									~*/
/*~			    GALOIS OPERATIONS				~*/
/*~									~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

GEN
transroot(GEN x, int i, int j)
{
  long n=lg(x),k;
  GEN y;

  y=cgetg(n,18);
  for(k=1;k<n;k++) y[k]=((k==i)||(k==j))?x[i+j-k]:x[k];
  return y;
}

GEN
tschirnhaus(GEN x)
{
  long av=avma,tetpil,v,n,a,b,c;
  GEN u;

  if(typ(x)!=10) err(galer1);
  n=lgef(x)-3;if(n<=0) err(galer1);v=varn(x);
  if(v) {u=gcopy(x);setvarn(u,0);x=u;}
  do
  {
    a=mymyrand()&3;if(!a) a=1;b=mymyrand()&7;if(b>=4) b-=8;
    c=mymyrand()&7;if(c>=4) c-=8;
    u=caract(gmodulcp(gaddsg(c,gmul(polx[0],gaddsg(b,gmulsg(a,polx[0])))),x),v);
  }
  while(lgef(polgcd(u,deriv(u,v)))>=4);
  tetpil=avma;return gerepile(av,tetpil,gcopy(u));
}

int
gpolcomp(GEN p1, GEN p2)
{
  int d,j;

  d=lgef(p1)-3;if((lgef(p2)-3)!=d) err(gpolcompbug1);
  j=d+1;while((j>=2)&&gegal(absi((GEN)p1[j]),absi((GEN)p2[j]))) j--;
  if(j==1) return 0;
  return gcmp(absi((GEN)p1[j]),absi((GEN)p2[j]));
}

GEN
galois(GEN x, long prec)
{

  long av=avma,av1,i,j,k,n,f,l,l2,e,e1;
  GEN x1,p1,p2,p3,p4,p5,p6,y,z,el;
  static int ind5[20]={2,5,3,4,1,3,4,5,1,5,2,4,1,2,3,5,1,4,2,3};
  static int ind6[60]={3,5,4,6,2,6,4,5,2,3,5,6,2,4,3,6,2,5,3,4,1,4,5,6,1,5,3,6,1,6,3,4,1,3,4,5,1,6,2,5,1,2,4,6,1,5,2,4,1,3,2,6,1,2,3,5,1,4,2,3};

  if(typ(x)!=10) err(galer1);
  n=lgef(x)-3;if(n<=0) err(galer1);
  if(n>7) err(impl,"galois of degree higher than 7");
  x=gdiv(x,content(x));
  for(i=2;i<=n+2;i++) if(typ((GEN)x[i])!=1) err(galer2);
  p1=(GEN)x[n+2];
  if(!gcmp1(p1))
  {
    x1=cgetg(n+3,10);x1[1]=x[1];x1[n+2]=un;p2=gun;
    for(i=n+1;i>=2;i--) {x1[i]=lmul((GEN)x[i],p2);if(i>2) p2=gmul(p1,p2);}
    x=x1;
  }
  p1=factor(x);
  if((lg((GEN)p1[1])>2)||(!gcmp1(gcoeff(p1,1,2)))) err(impl,"galois of reducible polynomial");
  x1=gcopy(x);av1=avma;
  for(;;)
  {
    switch(n)
    {
      case 1: avma=av;y=cgetg(4,17);y[1]=y[3]=un;y[2]=lneg(gun);return y;
      case 2: avma=av;y=cgetg(4,17);y[1]=deux;y[3]=un;y[2]=lneg(gun);return y;
      case 3: f=carreparfait(discsr(x));avma=av;y=cgetg(4,17);y[3]=un;
	if(f) {y[2]=un;y[1]=lstoi(3);return y;}
	else {y[2]=lneg(gun);y[1]=lstoi(6);return y;}
      case 4: 
	do
	{
	  p1=rootslong(x,prec);p2=p1;
	  p3=gzero;for(i=1;i<=4;i++) p3=gadd(p3,gmul((GEN)p2[i],gsqr((GEN)p2[(i&3)+1])));
	  p4=gsub(polx[0],p3);p2=transroot(p1,1,2);
	  p3=gzero;for(i=1;i<=4;i++) p3=gadd(p3,gmul((GEN)p2[i],gsqr((GEN)p2[(i&3)+1])));
	  p4=gmul(p4,gsub(polx[0],p3));p2=transroot(p1,1,3);
	  p3=gzero;for(i=1;i<=4;i++) p3=gadd(p3,gmul((GEN)p2[i],gsqr((GEN)p2[(i&3)+1])));
	  p4=gmul(p4,gsub(polx[0],p3));p2=transroot(p1,1,4);
	  p3=gzero;for(i=1;i<=4;i++) p3=gadd(p3,gmul((GEN)p2[i],gsqr((GEN)p2[(i&3)+1])));
	  p4=gmul(p4,gsub(polx[0],p3));p2=transroot(p1,2,3);
	  p3=gzero;for(i=1;i<=4;i++) p3=gadd(p3,gmul((GEN)p2[i],gsqr((GEN)p2[(i&3)+1])));
	  p4=gmul(p4,gsub(polx[0],p3));p2=transroot(p1,3,4);
	  p3=gzero;for(i=1;i<=4;i++) p3=gadd(p3,gmul((GEN)p2[i],gsqr((GEN)p2[(i&3)+1])));
	  p4=gmul(p4,gsub(polx[0],p3));p5=grndtoi(greal(p4),&e);
	  e=max(e,gexpo(gimag(p4)));
	  if(e> -10) prec=(prec<<2)-2;
	}
	while(e> -10);
	p6=ggcd(p5,deriv(p5,0));
	f=(typ(p6)==10)&&(lgef(p6)>3);
	if(f) goto tchi;
	p1=factor(p5);p2=(GEN)p1[1];l=lg(p2)-1;
	switch(l)
	{
	  case 1: f=carreparfait(discsr(x));avma=av;y=cgetg(4,17);y[3]=un;
	    if(f) {y[2]=un;y[1]=lstoi(12);return y;}
	    else {y[2]=lneg(gun);y[1]=lstoi(24);return y;}
	  case 2: avma=av;y=cgetg(4,17);y[3]=un;y[2]=lneg(gun);y[1]=lstoi(8);return y;
	  case 3: l2=lgef((GEN)p2[1])-3;
	    if(l2==2) {avma=av;y=cgetg(4,17);y[3]=y[2]=un;y[1]=lstoi(4);return y;}
	    else {avma=av;y=cgetg(4,17);y[3]=un;y[2]=lneg(gun);y[1]=lstoi(4);return y;}
	  default: err(galbug1);
	}
      case 5:
	do
	{
	  do
	  {
	    p1=rootslong(x,prec);z=cgetg(7,17);el=cgeti(7);
	    for(l=1;l<=5;l++)
	    {
	      p2=(l==1)?p1:transroot(p1,1,l);
	      p3=gzero;k=0;for(i=1;i<=5;i++)
	      {
		p5=gadd(gmul((GEN)p2[ind5[k]],(GEN)p2[ind5[k+1]]),gmul((GEN)p2[ind5[k+2]],(GEN)p2[ind5[k+3]]));
		p3=gadd(p3,gmul(gsqr((GEN)p2[i]),p5));k+=4;
	      }
	      z[l]=lrndtoi(greal(p3),&e);
	      el[l]=max(e,gexpo(gimag(p3)));
	      p4=(l==1)?gsub(polx[0],p3):gmul(p4,gsub(polx[0],p3));
	    }
	    p2=transroot(p1,2,5);
	    p3=gzero;k=0;for(i=1;i<=5;i++)
	    {
	      p5=gadd(gmul((GEN)p2[ind5[k]],(GEN)p2[ind5[k+1]]),gmul((GEN)p2[ind5[k+2]],(GEN)p2[ind5[k+3]]));
	      p3=gadd(p3,gmul(gsqr((GEN)p2[i]),p5));k+=4;
	    }
	    z[6]=lrndtoi(greal(p3),&e);
	    el[6]=max(e,gexpo(gimag(p3)));
	    p4=gmul(p4,gsub(polx[0],p3));
	    p5=grndtoi(greal(p4),&e);
	    e=max(e,gexpo(gimag(p4)));
	    if(e> -10) prec=(prec<<2)-2;
	  }
	  while(e> -10);
	  p6=ggcd(p5,deriv(p5,0));
	  f=(typ(p6)==10)&&(lgef(p6)>3);
	  if(f) goto tchi;
	  p3=factor(p5);l=lg((GEN)p3[1])-1;
	  f=carreparfait(discsr(x));
	  if(l==1)
	  {
	    avma=av;y=cgetg(4,17);y[3]=un;
	    if(f) {y[2]=un;y[1]=lstoi(60);return y;}
	    else {y[2]=lneg(gun);y[1]=lstoi(120);return y;}
	  }
	  else
	  {
	    if(f) 
	    {
	      l=1;while((l<=6)&&((el[l]>-((prec-2)<<(TWOPOTBITS_IN_LONG-1)))||(!gcmp0(poleval(p5,(GEN)z[l]))))) l++;
	      if(l>6) err(galbug4);
	      p2=(l==6)?transroot(p1,2,5):transroot(p1,1,l);
	      p3=gzero;
	      for(i=1;i<=5;i++)
	      {
		j=(i%5)+1;
		p3=gadd(p3,gmul(gmul((GEN)p2[i],(GEN)p2[j]),gsub((GEN)p2[j],(GEN)p2[i])));
	      }
	      p5=gmul(p3,p3);p4=grndtoi(greal(p5),&e1);
	      e1=max(e1,gexpo(gimag(p5)));
	      if(e1<= -10)
	      {
		if(gcmp0(p4)) goto tchi;
		f=carreparfait(p4);
		avma=av;y=cgetg(4,17);y[3]=y[2]=un;y[1]=lstoi(f?5:10);return y;
	      }
	      else prec=(prec<<2)-2;
	    }
	    else
	    {
	      avma=av;y=cgetg(4,17);y[3]=un;y[2]=lneg(gun);y[1]=lstoi(20);return y;
	    }
	  }
	}
	while(e1> -10);
      case 6: 
	do
	{
	  do
	  {
	    p1=rootslong(x,prec);
	    for(l=1;l<=6;l++)
	    {
	      p2=(l==1)?p1:transroot(p1,1,l);
	      p3=gzero;k=0;for(i=1;i<=5;i++) for(j=i+1;j<=6;j++)
	      {
		p5=gadd(gmul((GEN)p2[ind6[k]],(GEN)p2[ind6[k+1]]),gmul((GEN)p2[ind6[k+2]],(GEN)p2[ind6[k+3]]));
		p3=gadd(p3,gmul(gsqr(gmul((GEN)p2[i],(GEN)p2[j])),p5));k+=4;
	      }
	      p4=(l==1)?gsub(polx[0],p3):gmul(p4,gsub(polx[0],p3));
	    }
	    p5=grndtoi(greal(p4),&e);
	    e=max(e,gexpo(gimag(p4)));
	    if(e> -10) prec=(prec<<2)-2;
	  }
	  while(e> -10);
	  p6=ggcd(p5,deriv(p5,0));
	  f=(typ(p6)==10)&&(lgef(p6)>3);
	  if(f) goto tchi;
	  p3=factor(p5);p2=(GEN)p3[1];l=lg(p2)-1;
	  switch(l)
	  {
	    case 1:	p3=gadd(gmul(gmul((GEN)p1[1],(GEN)p1[2]),(GEN)p1[3]),gmul(gmul((GEN)p1[4],(GEN)p1[5]),(GEN)p1[6]));
	      p4=gsub(polx[0],p3);
	      for(i=1;i<=3;i++)
		for(j=4;j<=6;j++)
		{
		  p2=transroot(p1,i,j);
		  p3=gadd(gmul(gmul((GEN)p2[1],(GEN)p2[2]),(GEN)p2[3]),gmul(gmul((GEN)p2[4],(GEN)p2[5]),(GEN)p2[6]));
		  p4=gmul(p4,gsub(polx[0],p3));
		}
	      p5=grndtoi(greal(p4),&e1);
	      e1=max(e1,gexpo(gimag(p4)));
	      if(e1<= 10)
	      {
		p6=ggcd(p5,deriv(p5,0));
		f=(typ(p6)==10)&&(lgef(p6)>3);
		if(f) goto tchi;
		p3=factor(p5);p2=(GEN)p3[1];l=lg(p2)-1;f=carreparfait(discsr(x));
		avma=av;y=cgetg(4,17);y[3]=un;
		if(l==1)
		{
		  if(f) {y[2]=un;y[1]=lstoi(360);return y;}
		  else {y[2]=lneg(gun);y[1]=lstoi(720);return y;}
		}
		else
		{
		  if(f) {y[2]=un;y[1]=lstoi(36);return y;}
		  else {y[2]=lneg(gun);y[1]=lstoi(72);return y;}
		}
	      }
	      else prec=(prec<<2)-2;
	      break;
		  
	    case 2: l2=lgef((GEN)p2[1])-3;if(l2>3) l2=6-l2;
	      switch(l2)
	      {
		case 1: f=carreparfait(discsr(x));avma=av;y=cgetg(4,17);y[3]=un;
		  if(f) {y[2]=un;y[1]=lstoi(60);return y;}
		  else {y[2]=lneg(gun);y[1]=lstoi(120);return y;}
		case 2: f=carreparfait(discsr(x));
		  if(f) {avma=av;y=cgetg(4,17);y[3]=y[2]=un;y[1]=lstoi(24);return y;}
		  else
		  {
		    p3=(lgef((GEN)p2[1])==5)?(GEN)p2[2]:(GEN)p2[1];
		    f=carreparfait(discsr(p3));avma=av;y=cgetg(4,17);y[2]=lneg(gun);
		    if(f) {y[1]=lstoi(24);y[3]=deux;return y;}
		    else {y[1]=lstoi(48);y[3]=un;return y;}
		  }
		case 3: f=carreparfait(discsr((GEN)p2[1]))||carreparfait(discsr((GEN)p2[2]));
		  avma=av;y=cgetg(4,17);y[3]=un;y[2]=lneg(gun);y[1]=lstoi(f?18:36);
		  return y;
	      }
	    case 3: for(l2=1;l2<=3;l2++) if(lgef((GEN)p2[l2])>=6) p3=(GEN)p2[l2];
	      if(lgef(p3)==6)
	      {
		f=carreparfait(discsr(p3));avma=av;y=cgetg(4,17);y[2]=lneg(gun);
		y[3]=un;y[1]=f?lstoi(6):lstoi(12);return y;
	      }
	      else
	      {
		f=carreparfait(discsr(x));avma=av;y=cgetg(4,17);y[3]=un;
		if(f) {y[2]=un;y[1]=lstoi(12);return y;}
		else {y[2]=lneg(gun);y[1]=lstoi(24);return y;}
	      }
	    case 4: avma=av;y=cgetg(4,17);y[1]=lstoi(6);y[2]=lneg(gun);
	      y[3]=deux;return y;
	    default: err(galbug3);
	  }
	}
	while(e1> -10);
	  
      case 7: 
	do
	{
	  p1=rootslong(x,prec);p4=gun;
	  for(i=1;i<=5;i++)
	    for(j=i+1;j<=6;j++)
	      for(k=j+1;k<=7;k++)
		p4=gmul(p4,gsub(polx[0],gadd(gadd((GEN)p1[i],(GEN)p1[j]),(GEN)p1[k])));
	  p5=grndtoi(greal(p4),&e);e=max(e,gexpo(gimag(p4)));
	  if(e> -10) prec=(prec<<2)-2;
	}
	while(e> -10);
	p6=ggcd(p5,deriv(p5,0));
	f=(typ(p6)==10)&&(lgef(p6)>3);
	if(f) goto tchi;
	p1=factpol(p5,0,7);p2=(GEN)p1[1];l=lg(p2)-1;
	switch(l)
	{
	  case 1: f=carreparfait(discsr(x));avma=av;y=cgetg(4,17);y[3]=un;
	    if(f) {y[2]=un;y[1]=lstoi(2520);return y;}
	    else {y[2]=lneg(gun);y[1]=lstoi(5040);return y;}
	  case 2: f=lgef((GEN)p2[1])-3;avma=av;y=cgetg(4,17);y[3]=un;
	    if((f==7)||(f==28)) {y[2]=un;y[1]=lstoi(168);return y;}
	    else {y[2]=lneg(gun);y[1]=lstoi(42);return y;}
	  case 3: avma=av;y=cgetg(4,17);y[3]=y[2]=un;y[1]=lstoi(21);return y;
	  case 4: avma=av;y=cgetg(4,17);y[3]=un;y[2]=lneg(gun);y[1]=lstoi(14);return y;
	  case 5: avma=av;y=cgetg(4,17);y[3]=y[2]=un;y[1]=lstoi(7);return y;
	  default: err(galbug2);
	}
    }
    tchi:
    avma=av1;x=tschirnhaus(x1);
    if(DEBUGLEVEL>=2)
    {
      fprintferr("transformation de Tschirnhaus: nouveau polynome ");
      outerr(x);flusherr();
    }
  }
}

long
computehenselbound(GEN nf, GEN p, long prec)
/* a usage interne de galoisconj */
{
  long n,r1,r2,ru,i,j,e1,ep,e,av=avma;
  GEN p1,m,mm,pt5,pt6,mi,pmax;
  
  n=lgef((GEN)nf[1])-3;p1=(GEN)nf[2];r1=itos((GEN)p1[1]);r2=itos((GEN)p1[2]);
  ru=r1+r2;pt5=(GEN)nf[5];mm=(GEN)pt5[1];
  m=cgetg(n+1,19);
  for(j=1;j<=n;j++)
  {
    p1=cgetg(n+1,18);m[j]=(long)p1;
    for(i=1;i<=n;i++)
      p1[i]=(i<=ru)?coeff(mm,i,j):(long)gconj(gcoeff(mm,i-r2,j));
  }
  mi=gmul(m,(GEN)pt5[6]);
  pmax=gzero;for(j=1;j<=n;j++) pmax=gmax(pmax,gnorml2((GEN)mi[j]));
  pt6=(GEN)nf[6];
  p1=gzero;for(i=1;i<=r1;i++) p1=gadd(p1,gnorm((GEN)pt6[i]));
  for(;i<=ru;i++) p1=gadd(p1,gmul2n(gnorm((GEN)pt6[i]),1));
  p1=gaddsg(1,gmulsg(n,gsqr(gmul(p1,pmax))));
  e1=gexpo(p1);ep=expi(p);
  avma=av;e=(e1+n+1)/ep+3; /* a modifier */
  return e;
}
  
GEN
galoisconj(GEN nf, long prec)
{
  long av=avma,tetpil,i,j,k,n,v,c,e,fl;
  GEN x,y,p1,p2,p,a,a1,m,fa,pe,dp,bas,mred;
  byteptr pt=diffptr;

  if(DEBUGLEVEL>1) {fprintferr("Entree dans galoisconj()\n");flusherr();}
  if(typ(nf)==10)
    err(talker,"galoisconj wants a number field, not a polynomial;\n        please apply initalg first");
  nf=checknf(nf);x=(GEN)nf[1];
  if(typ(x)!=10) err(galer1);
  n=lgef(x)-3;if(n<=0) return cgetg(1,17);
  v=varn(x);p=cgeti(3);p[1]=evallgef(3)+evalsigne(1);p[2]=0;fl=1;
  dp=mulii((GEN)nf[3],(GEN)nf[4]);
  while((c= *pt++)&&fl)
  {
    p[2]+=c;
    if(cmpis(p,2)&&(!divise(dp,p)))
    {fa=simplefactmod(x,p);fl=(lg((GEN)fa[1])<(n+1));}
  }
  if(!c) err(primer1);
  e=computehenselbound(nf,p,prec);
  fa=gtrunc(rootpadicfast(x,p,e,1));
  pe=gpuigs(p,e);
  m=idmat(n+1);
  coeff(m,1,1)=(long)pe;
  a=gmodulcp((GEN)fa[1],pe);a1=centerlift(gneg(gsubst((GEN)nf[7],v,a)));
  for(j=2;j<=n;j++) coeff(m,1,j)=a1[j];
  mred=cgetg(n+1,19);for(j=1;j<=n;j++) mred[j]=m[j];
  mred=gmul(mred,lllint(mred));for(j=1;j<=n;j++) m[j]=mred[j];
  y=cgetg(n+1,17);y[1]=(long)polx[v];
  for(i=2;i<=n;i++)
  {
    coeff(m,1,n+1)=(long)centerlift(gneg(gmodulcp((GEN)fa[i],pe)));
    mred=(GEN)gmul(m,lllint(m));fl=0;
    for(k=1;(k<=n)&&(!fl);k++)
    {
      p1=(GEN)mred[k];
      p2=(GEN)p1[n+1];if(signe(p2)>0) {p1=gneg(p1);p2=(GEN)p1[n+1];}
      if(gcmp_1(p2))
      {
	bas=(GEN)nf[7];
	p2=gzero;for(j=1;j<=n;j++) p2=gadd(p2,gmul((GEN)p1[j],(GEN)bas[j]));
	fl=gcmp0(gsubst(x,v,gmodulcp(p2,x)));
      }
    }
    y[i]=fl?(long)p2:zero;
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(y));
}

GEN
isisomfastall(GEN nf1, GEN nf2, long fliso, long prec)
{
  long av=avma,tetpil,i,j,k,n,n1,v,c,e,fl,count;
  GEN x1,x2,y,p1,p2,p,a,a1,m,fa1,fa2,pe,dp,bas,mred;
  byteptr pt=diffptr;

  if((typ(nf1)==10)||(typ(nf2)==10))
    err(talker,"isisomfast and isinclfast want two number fields, not polynomials;\n        please apply initalg first");
  nf1=checknf(nf1);nf2=checknf(nf2);x1=(GEN)nf1[1];x2=(GEN)nf2[1];
  if((typ(x1)!=10)||(typ(x2)!=10)) err(nfisoer1);
  if(fliso)
  {
    if(!gegal((GEN)nf1[2],(GEN)nf2[2])) {avma=av;return gzero;}
    if(!gegal((GEN)nf1[3],(GEN)nf2[3])) {avma=av;return gzero;}
    n1=n=lgef(x1)-3;if(n<=0) err(nfisoer2);
  }
  else
  {
    n=lgef(x2)-3;n1=lgef(x1)-3;if((n<=0)||(n1<=0)) err(nfisoer2);
    if((n%n1)||(!divise((GEN)nf2[3],gpuigs((GEN)nf1[3],n/n1))))
    {avma=av;return gzero;}
  }
  v=varn(x1);p=cgeti(3);p[1]=evallgef(3)+evalsigne(1);p[2]=0;fl=1;
  dp=mulii(mulii((GEN)nf2[3],(GEN)nf1[4]),(GEN)nf2[4]);
  while((c= *pt++)&&fl)
  {
    p[2]+=c;
    if(cmpis(p,2)&&(!divise(dp,p)))
    {
      fa2=simplefactmod(x2,p);
      fl=(lg((GEN)fa2[1])<(n+1));
      if((!fl)&&(lg((GEN)simplefactmod(x1,p)[1])<(n1+1)))
      {avma=av;return gzero;}
    }
  }
  if(!c) err(primer1);
  e=max(computehenselbound(nf2,p,prec),computehenselbound(nf1,p,prec));
      /* a modifier nf1 et nf2 doivent etre pris en compte a la fois */
  fa2=gtrunc(rootpadicfast(x2,p,e,0));
  fa1=gtrunc(rootpadicfast(x1,p,e,1));  
  pe=gpuigs(p,e);
  m=idmat(n+1);
  coeff(m,1,1)=(long)pe;
  a=gmodulcp((GEN)fa2[1],pe);a1=centerlift(gneg(gsubst((GEN)nf2[7],v,a)));
  for(j=2;j<=n;j++) coeff(m,1,j)=a1[j];
  mred=cgetg(n+1,19);for(j=1;j<=n;j++) mred[j]=m[j];
  mred=gmul(mred,lllint(mred));for(j=1;j<=n;j++) m[j]=mred[j];
  y=cgetg(n1+1,17);count=0;
  for(i=1;i<=n1;i++)
  {
    coeff(m,1,n+1)=(long)centerlift(gneg(gmodulcp((GEN)fa1[i],pe)));
    mred=(GEN)gmul(m,lllint(m));fl=0;
    for(k=1;(k<=n)&&(!fl);k++)
    {
      p1=(GEN)mred[k];
      p2=(GEN)p1[n+1];if(signe(p2)>0) {p1=gneg(p1);p2=(GEN)p1[n+1];}
      if(gcmp_1(p2))
      {
	bas=(GEN)nf2[7];
	p2=gzero;for(j=1;j<=n;j++) p2=gadd(p2,gmul((GEN)p1[j],(GEN)bas[j]));
	fl=gcmp0(gsubst(x1,v,gmodulcp(p2,x2)));
      }
    }
    if(fl) {y[i]=(long)p2;count++;} else y[i]=zero;
  }
  if(!count) {avma=av;return gzero;}
  else
  {
    tetpil=avma;p1=cgetg(count+1,17);
    for(count=0,i=1;i<=n1;i++)
      if(signe((GEN)y[i])) p1[++count]=lcopy((GEN)y[i]);
    return gerepile(av,tetpil,p1);
  }
}

GEN
isisomfast(GEN nf1, GEN nf2, long prec)
{
  return isisomfastall(nf1,nf2,1,prec);
}

GEN
isinclfast(GEN nf1, GEN nf2, long prec)
{
  return isisomfastall(nf1,nf2,0,prec);
}

GEN
henselstep(GEN x, GEN polr, long e)
{
  long i,n;
  GEN y,yi,xp;
  
  n=lgef(x)-3;xp=deriv(x,varn(x));
  y=cgetg(n+1,17);
  for(i=1;i<=n;i++)
  {
    yi=gprec((GEN)polr[i],e<<1);
    y[i]=lsub(yi,gdiv(poleval(x,yi),poleval(xp,yi)));
  }
  return y;
}

GEN
galoisconjforce(GEN nf, long prec)
{
  long av=avma,tetpil,i,j,n,v,c,e,fl;
  GEN x,y,p1,p2,p,a,a1,m,fa,pe,dp,bas,mred,polr;
  byteptr pt=diffptr;

  if(DEBUGLEVEL>1) {fprintferr("Entree dans galoisconjforce()\n");flusherr();}
  if(typ(nf)==10)
    err(talker,"galoisconjforce wants a number field, not a polynomial;\n        please apply initalg first");
  nf=checknf(nf);x=(GEN)nf[1];
  if(typ(x)!=10) err(galer1);
  n=lgef(x)-3;if(n<=0) return cgetg(1,17);
  v=varn(x);p=cgeti(3);p[1]=evallgef(3)+evalsigne(1);p[2]=0;fl=1;
  dp=mulii((GEN)nf[3],(GEN)nf[4]);
  while((c= *pt++)&&fl)
  {
    p[2]+=c;
    if(cmpis(p,2)&&(!divise(dp,p)))
    {fa=simplefactmod(x,p);fl=(lg((GEN)fa[1])<(n+1));}
  }
  if(!c) err(primer1);
  e=computehenselbound(nf,p,prec);
  polr=rootpadicfast(x,p,e,1);pe=gpuigs(p,e);
  y=cgetg(n+1,17);y[1]=(long)polx[v];for(i=2;i<=n;i++) y[i]=zero;
  do
  {
    fa=gtrunc(polr);m=idmat(n+1);coeff(m,1,1)=(long)pe;
    a=gmodulcp((GEN)fa[1],pe);a1=centerlift(gneg(gsubst((GEN)nf[7],v,a)));
    for(j=2;j<=n;j++) coeff(m,1,j)=a1[j];
    mred=cgetg(n+1,19);for(j=1;j<=n;j++) mred[j]=m[j];
    mred=gmul(mred,lllint(mred));for(j=1;j<=n;j++) m[j]=mred[j];
    for(fl=1,i=2;(i<=n)&&fl;i++)
    {
      if(!signe((GEN)y[i]))
      {
	coeff(m,1,n+1)=(long)centerlift(gneg(gmodulcp((GEN)fa[i],pe)));
	mred=(GEN)gmul(m,lllint(m));
	p1=(GEN)mred[1];
	p2=(GEN)p1[n+1];if(signe(p2)>0) {p1=gneg(p1);p2=(GEN)p1[n+1];}
	if(gcmp_1(p2))
	{
	  bas=(GEN)nf[7];
	  p2=gzero;for(j=1;j<=n;j++) p2=gadd(p2,gmul((GEN)p1[j],(GEN)bas[j]));
	  fl=gcmp0(gsubst(x,v,gmodulcp(p2,x)));
	  if(fl) y[i]=(long)p2;
	}
	else fl=0;
      }
    }
    if(!fl) {polr=henselstep(x,polr,e);e<<=1;pe=gsqr(pe);}
  }
  while(!fl);
  tetpil=avma;return gerepile(av,tetpil,gcopy(y));
}

GEN
galoisconj1(GEN nf, long prec)
{
  long av=avma,tetpil,i,j,n,v;
  GEN x,y,w,polr,p1,p2,b;

  if(DEBUGLEVEL>1) {fprintferr("Entree dans galoisconj1()\n");flusherr();}
  if(typ(nf)==10)
    err(talker,"galoisconj1 wants a number field, not a polynomial;\n        please apply initalg first");
  nf=checknf(nf);x=(GEN)nf[1];
  if(typ(x)!=10) err(galer1);
  n=lgef(x)-3;if(n<=0) return cgetg(1,17);
  v=varn(x);polr=(GEN)nf[6];p1=(GEN)polr[1];b=(GEN)nf[7];
  p2=(GEN)((GEN)nf[5])[1];
  w=cgetg(n+1,17);for(i=1;i<=n;i++) w[i]=coeff(p2,1,i);
  y=cgetg(n+1,17);y[1]=(long)polx[v];
  for(i=2;i<=n;i++)
  {
    p1=lindep2(concat(w,(GEN)polr[i]),(prec-2)*2*BYTES_IN_LONG);
    if(gcmp0((GEN)p1[n+1])) y[i]=zero;
    else
    {
      p2=gzero;
      for(j=1;j<=n;j++) p2=gadd(p2,gmul((GEN)p1[j],(GEN)b[j]));
      p2=gdiv(p2,gneg((GEN)p1[n+1]));
      if(gcmp0(gmod(gsubst(x,v,p2),x))) y[i]=(long)p2;else y[i]=zero;
    }
    if(DEBUGLEVEL>1) outerr((GEN)y[i]);
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(y));
}

GEN
galoisconj2(GEN x, long prec)
{
  long av=avma,tetpil,i,j,n,v;
  GEN y,w,polr,p1,p2;

  if(DEBUGLEVEL>1){fprintferr("Entree dans galoisconj2()\n");flusherr();}
  if(typ(x)!=10) err(galer1);
  n=lgef(x)-3;if(n<=0) return cgetg(1,17);
  v=varn(x);p1=factor(x);
  if((lg((GEN)p1[1])>2)||(!gcmp1(gcoeff(p1,1,2)))) err(galcer1);
  polr=rootslong(x,prec);p1=(GEN)polr[1];
  w=cgetg(n+1,17);w[1]=un;for(i=2;i<=n;i++) w[i]=lmul(p1,(GEN)w[i-1]);
  y=cgetg(n+1,17);y[1]=(long)polx[v];
  if(DEBUGLEVEL>1){fprintferr("Calcul des conjugues...\n");flusherr();}
  for(i=2;i<=n;i++)
  {
    p1=lindep2(concat(w,(GEN)polr[i]),(prec-2)*2*BYTES_IN_LONG);
/*  p1=lindep(concat(w,(GEN)polr[i]),prec);*/
    if(gcmp0((GEN)p1[n+1])) y[i]=zero;
    else
    {
      p2=gzero;
      for(j=n;j;j--) p2=gadd((GEN)p1[j],gmul(p2,polx[v]));
      p2=gdiv(p2,gneg((GEN)p1[n+1]));
      if(gcmp0(gmod(gsubst(x,v,p2),x))) y[i]=(long)p2;else y[i]=zero;
    }
    if(DEBUGLEVEL>1) outerr((GEN)y[i]);
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(y));
}

GEN
galoisapply(GEN nf, GEN aut, GEN x)
{
  long av=avma,tetpil,tx=typ(x),lx,j,N;
  GEN p1,y,pol,unmod;
  
  nf=checknf(nf);pol=(GEN)nf[1];
  if(typ(aut)==10) aut=gmodulcp(aut,pol);
  else
  {
    if((typ(aut)!=9)||(!gegal((GEN)aut[1],pol))||gcmp0((GEN)aut[1]))
      err(talker,"incorrect galois automorphism in galoisapply");
  }
  switch(tx)
  {
    case 1: case 3: case 4: case 5: case 7:
      tetpil=avma;return gerepile(av,tetpil,gcopy(x));
    case 9:
      tetpil=avma;p1=gsubst((GEN)x[2],varn(pol),aut);
      if((typ(p1)!=9)||(!gegal((GEN)p1[1],pol))) {tetpil=avma;p1=gmodulcp(p1,pol);}
      return gerepile(av,tetpil,p1);
    case 10:
      tetpil=avma;p1=gsubst(x,varn(pol),aut);
      if((typ(p1)!=9)||(!gegal((GEN)p1[1],pol))) {tetpil=avma;p1=gmodulcp(p1,pol);}
      return gerepile(av,tetpil,p1);
    case 17:
      if(lg(x)==3)
      {
	tetpil=avma;y=cgetg(3,17);y[1]=(long)galoisapply(nf,aut,(GEN)x[1]);
	y[2]=lcopy((GEN)x[2]);return gerepile(av,tetpil,y);
      }
      if(lg(x)!=6) err(talker,"incorrect type in galoisapply");
      y=cgetg(6,17);y[1]=x[1];unmod=gmodulcp(gun,(GEN)x[1]);
      p1=centerlift(gmul(unmod,algtobasis(nf,galoisapply(nf,aut,(GEN)x[2]))));
      if(gcmp1((GEN)x[3]))
      {
	if(ggval(subres(gmul((GEN)nf[7],p1),pol),(GEN)x[1])>itos((GEN)x[4]))
	{
	  if(signe((GEN)p1[1])>0) p1[1]=lsub((GEN)p1[1],(GEN)x[1]);
	  else p1[1]=ladd((GEN)p1[1],(GEN)x[1]);
	}
      }
      y[2]=(long)p1;y[3]=x[3];y[4]=x[4];
      y[5]=(long)centerlift(gmul(unmod,algtobasis(nf,galoisapply(nf,aut,(GEN)x[5]))));
      tetpil=avma;return gerepile(av,tetpil,gcopy(y));
    case 18:
      N=lgef(pol)-3;
      if(lg(x)!=N+1) err(talker,"incorrect type in galoisapply");
      p1=gmul((GEN)nf[7],x);tetpil=avma;
      return gerepile(av,tetpil,galoisapply(nf,aut,p1));
    case 19:
      lx=lg(x);if(lx==1) return gcopy(x);
      N=lgef(pol)-3;
      if(lg((GEN)x[1])!=N+1) err(talker,"incorrect type in galoisapply");
      tetpil=avma;p1=cgetg(lx,19);
      for(j=1;j<lx;j++) p1[j]=(long)galoisapply(nf,aut,(GEN)x[j]);
      if(lg(x)!=N+1) return gerepile(av,tetpil,p1);
      else {tetpil=avma;return gerepile(av,tetpil,idealhermite(nf,p1));}
    default:
      err(talker,"incorrect type in galoisapply");return gnil;
  }
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~									~*/
/*~			       INITALG					~*/
/*~									~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

GEN
initalgredall(GEN x, long flun, long prec)
{
  GEN y,yy,p1,p2,p3,p4,p5,p6,p7,p10,p11,p12,p20,fieldd,polr,ptrace;
  GEN dx,adx,polmax,s,a,dxn,adxn,sn,phimax,ind,fa;
  long tx=typ(x),n=lgef(x)-3,i,j,av=avma,av1,av2,av3,v,k,j1,j2,lgp,r1,r2,ru,imax,numb,flc,tetpil;

  if(DEBUGLEVEL) timer2();
  if((tx!=10)||(n<=0)) err(poltyper);
  p1=content(x);p1=gcmp1(p1) ? x: gdiv(x,content(x));
  for(k=2;k<=n+2;k++) if(typ((GEN)p1[k])!=1) err(impl,"general algebraic extension");
  if(gcmp1(p2=(GEN)p1[n+2])) x=p1;
  else
  {
    x=cgetg(n+3,10);x[1]=p1[1];x[n+2]=un;x[n+1]=p1[n+1];p3=p2;
    for(k=n;k>=2;k--)
    {
      x[k]=lmulii(p3,(GEN)p1[k]);
      if(k>2) p3=mulii(p2,p3);
    }
  }
  p1=factor(x);
  if(DEBUGLEVEL) {fprintferr("temps factpol: ");fprintferr("%ld\n",timer2());flusherr();}
  if(lgef(gcoeff(p1,1,1))!=n+3) err(redper1);
  r1=sturm(x);p4=allbase4(x,0,&fieldd,&fa);
  if(DEBUGLEVEL) {fprintferr("temps round4: ");fprintferr("%ld\n",timer2());flusherr();}
  if(r1<n)
  {
    polr=rootslong(x,prec);p3=cgetg(n+1,19);
    for(i=1;i<=n;i++)
    {
      p1=cgetg(n+1,18);p3[i]=(long)p1;
      for(j=1;j<=n;j++)
	p1[j]=lsubst((GEN)p4[i],varn((GEN)p4[i]),(GEN)polr[j]);
    }
    p2=greal(gmul(gconj(gtrans(p3)),p3));
    p1=lllgram(p2,prec);
    for(s=gzero,i=1;i<=n;i++) s=gadd(s,gnorm((GEN)polr[i]));
  }
  else
  {
    ptrace=cgetg(n+1,17);ptrace[1]=lstoi(n);
    for(k=1;k<n;k++) 
    {
      p3=gmulsg(k,(GEN)x[n-k+2]);
      for(i=1;i<k;i++) p3=gadd(p3,gmul((GEN)x[n-i+2],(GEN)ptrace[k-i+1]));
      ptrace[k+1]=lneg(p3);
    }
    p2=cgetg(n+1,19);
    for(i=1;i<=n;i++)
    {
      p1=cgetg(n+1,18);p2[i]=(long)p1;
      for(j=1;j<i;j++) p1[j]=lcopy(gcoeff(p2,i,j));
      for(j=i;j<=n;j++)
      {
	p5=gres(gmul((GEN)p4[i],(GEN)p4[j]),x);p6=gzero;
	for(k=0;k<=lgef(p5)-3;k++) p6=gadd(p6,gmul((GEN)p5[k+2],(GEN)ptrace[k+1]));
	p1[j]=(long)p6;
      }
    }
    p1=lllgramint(p2);
    s=(n>1)?gsub(gsqr((GEN)x[n+1]),gmul2n((GEN)x[n],1)):gsqr((GEN)x[2]);
  }
  if(DEBUGLEVEL) {fprintferr("temps matrice T2: ");fprintferr("%ld\n",timer2());flusherr();}
  v=varn(x);dx=discsr(x);adx=absi(dx);polmax=x;imax=0;
  a=cgetg(n+1,18);for(i=1;i<=n;i++) a[i]=lmul(p4,(GEN)p1[i]);
  for(numb=0,i=1;i<=n;i++)
  {
    av1=avma;p3=gmodulcp((GEN)a[i],x);p7=content((GEN)p3[2]);
    if(gcmp1(p7)) p3=caract(p3,v);
    else
    {
      p3=caract(gdiv(p3,p7),v);
      p3=gmul(gpuigs(p7,lgef(p3)-3),gsubst(p3,v,gdiv(polx[v],p7)));
    }
    p5=ggcd(deriv(p3,v),p3);
    if(lgef(p5)==3)
    {
      dxn=discsr(p3);adxn=absi(dxn);flc=gcmp(adxn,adx);numb++;
      if(flc<=0)
      {
	if(r1<n) for(sn=gzero,j=1;j<=n;j++) sn=gadd(sn,gnorm(poleval((GEN)a[i],(GEN)polr[j])));
	else sn=(n>1)?gsub(gsqr((GEN)p3[n+1]),gmul2n((GEN)p3[n],1)):gsqr((GEN)p3[2]);
	if(flc<0) {dx=dxn;adx=adxn;s=sn;polmax=p3;imax=i;}
	else 
	{
	  flc=gcmp(sn,s);
	  if((flc<0)||((!flc)&&(gpolcomp(p3,polmax)<0)))
	  {dx=dxn;adx=adxn;s=sn;polmax=p3;imax=i;}
	}
      }
    }
  }
  if(!numb) err(polreder1);
  phimax=imax?(GEN)a[imax]:polx[v];
  j=n+1;
  while((j>=2)&&(!signe((GEN)polmax[j]))) j-=2;
  if((j>=2)&&(signe((GEN)polmax[j])>0))
  {
    if(polmax==x) polmax=gcopy(x);
    for(;j>=2;j-=2) setsigne((GEN)polmax[j],-signe((GEN)polmax[j]));
    phimax=gneg(phimax);
  }
  if(DEBUGLEVEL) {fprintferr("temps polmax: ");fprintferr("%ld\n",timer2());flusherr();}
  p2=gmodulcp(phimax,x);p20=polymodrecip(p2);
  p2=gcmp0(p2)?p4:lift(gsubst(p4,v,p20));
  p3=cgetg(n+1,19);
  for(j=1;j<=n;j++)
  {
    p4=cgetg(n+1,18);p3[j]=(long)p4;
    for(i=1;i<=n;i++) p4[i]=(long)truecoeff((GEN)p2[j],i-1);
  }
  p4=denom(p3);p1=gmul(p4,p3);
  p2=gdiv(hnfmod(p1,detint(p1)),p4);p3=cgetg(n+1,17);
  for(j=1;j<=n;j++) 
  {
    p1=gzero;for(i=n;i;i--) p1=gadd(gcoeff(p2,i,j),gmul(p1,polx[v]));
    p3[j]=(long)p1;
  }
  if(!carrecomplet(divii(dx,fieldd),&ind)) err(initalgbug1);
  p1=cgetg(n+1,19);
  for(j=1;j<=n;j++)
  {
    p2=cgetg(n+1,18);p1[j]=(long)p2;
    for(i=1;i<=n;i++) p2[i]=(long)truecoeff((GEN)p3[j],i-1);
  }
  p5=cgetg(n*n+1,19);
  for(j=1;j<=n*n;j++) p5[j]=lgetg(n+1,18);
  p11=cgetg(n+1,19);for(j=1;j<=n;j++) p11[j]=lgetg(n+1,18);
  p12=polsym(polmax,lgef(polmax)-4);
  for(i=1;i<=n;i++)
    for(j=i;j<=n;j++)
    {
      p10=gmod(gmul((GEN)p3[j],(GEN)p3[i]),polmax);
      j1=j+(i-1)*n;j2=i+(j-1)*n;lgp=lgef(p10)-2;
      for(k=1;k<=lgp;k++) coeff(p5,k,j1)=coeff(p5,k,j2)=p10[k+1];
      for(;k<=n;k++) coeff(p5,k,j1)=coeff(p5,k,j2)=zero;
      coeff(p11,i,j)=coeff(p11,j,i)=(long)trace9(p10,p12);
    }
  if(DEBUGLEVEL) {fprintferr("temps table de mult: ");fprintferr("%ld\n",timer2());flusherr();}
  p2=rootslong(polmax,prec);
  if(DEBUGLEVEL) {fprintferr("temps racines: ");fprintferr("%ld\n",timer2());flusherr();}
  for(i=1;i<=r1;i++) p2[i]=lreal((GEN)p2[i]);
  tetpil=avma;y=cgetg(10,17);y[1]=lcopy(polmax);y[8]=(long)ginv(p1);
  y[9]=lmul((GEN)y[8],p5);p1=cgetg(3,17);
  p1[1]=lstoi(r1);p1[2]=lstoi(r2=(n-r1)>>1);y[2]=(long)p1;ru=r1+r2;
  y[3]=lcopy(fieldd);y[4]=lcopy(ind);
  p4=cgetg(n+1,19);
  for(j=1;j<=n;j++)
  {
    p6=cgetg(ru+1,18);p4[j]=(long)p6;
    for(i=1;i<=ru;i++) p6[i]=(long)gsubst((GEN)p3[j],v,(GEN)p2[max(i,(i<<1)-r1)]);
  }
  av2=avma;p7=cgetg(8,17);p7[1]=(long)p4;
  p5=cgetg(ru+1,19);
  for(j=1;j<=ru;j++)
  {
    p6=cgetg(n+1,18);p5[j]=(long)p6;
    for(i=1;i<=n;i++) 
      p6[i]=(j<=r1)?lconj(gcoeff(p4,j,i)):lmul2n(gconj(gcoeff(p4,j,i)),1);
  }
  p7[2]=(long)p5;p7[3]=(long)greal(gmul(p5,p4));
  p7[4]=lcopy(p11);p7[5]=zero;
  p7[6]=lmul(ginv(p11), fieldd);
         /* discriminant * (inverse of y[5][4] = p7[4]) */
  p7[7]=(long)hnfmod((GEN)p7[6],detint((GEN)p7[6]));
         /* Ideal basis for discriminant * (inverse of different) */
  av3=avma;
  y[5]=lpile(av2,av3,gcopy(p7));p4=cgetg(ru+1,17);y[6]=(long)p4;
  for(i=1;i<=ru;i++) p4[i]=lcopy((GEN)p2[max(i,(i<<1)-r1)]);
  if(DEBUGLEVEL) {fprintferr("temps matrices: ");fprintferr("%ld\n",timer2());flusherr();}
  y[7]=lcopy(p3);
  ((GEN)y[5])[5]=(long)differente(y,fa);
/*  av2=avma;p1=idealinv(y,(GEN)(((GEN)y[5])[7]));av3=avma;
  ((GEN)y[5])[5]=lpile(av2,av3,gmul(fieldd,p1)); */
  if(DEBUGLEVEL) {fprintferr("temps fin differente: ");fprintferr("%ld\n",timer2());flusherr();}
  if(flun) {yy=cgetg(3,17);yy[1]=(long)y;yy[2]=lcopy(p20);}
  else yy=y;
  return gerepile(av,tetpil,yy);
}

GEN
initalgred(GEN x, long prec)
{
  return initalgredall(x,0,prec);
}

GEN
initalgred2(GEN x, long prec)
{
  return initalgredall(x,1,prec);
}

GEN
initalgall0(GEN x, long fldif, long prec)
{
  GEN y,p1,p2,p3,p4,p5,p6,p7,p10,p11,p12,fieldd,dx,ind,fa;
  long tx=typ(x),n=lgef(x)-3,i,j,av=avma,av2,av3,v,k,r1,r2,ru,tetpil,j1,j2,lgp;
  long PRECREG,PRECREGINT,parikeep;

  if(DEBUGLEVEL) timer2();
  if((tx!=10)||(n<=0)) err(poltyper);
  for(k=2;k<=n+2;k++) if(typ((GEN)x[k])!=1) err(initer1);
  if(!gcmp1((GEN)x[n+2])) err(initer2);
  p1=factor(x);
  if(DEBUGLEVEL) {fprintferr("temps factpol: ");fprintferr("%ld\n",timer2());flusherr();}
  if(lgef(gcoeff(p1,1,1))!=n+3) err(redper1);
  r1=sturm(x);p3=allbase4(x,0,&fieldd,&fa);
  if(DEBUGLEVEL) {fprintferr("temps round4: ");fprintferr("%ld\n",timer2());flusherr();}
  PRECREGINT=(gexpo(fieldd)>>TWOPOTBITS_IN_LONG)+n+3;
  PRECREG=prec+PRECREGINT;
  v=varn(x);dx=discsr(x);
  if(!carrecomplet(divii(dx,fieldd),&ind)) err(initalgbug1);
  p1=cgetg(n+1,19);
  for(j=1;j<=n;j++)
  {
    p2=cgetg(n+1,18);p1[j]=(long)p2;
    for(i=1;i<=n;i++) p2[i]=(long)truecoeff((GEN)p3[j],i-1);
  }
  p5=cgetg(n*n+1,19);
  for(j=1;j<=n*n;j++) p5[j]=lgetg(n+1,18);
  p11=cgetg(n+1,19);for(j=1;j<=n;j++) p11[j]=lgetg(n+1,18);
  p12=polsym(x,lgef(x)-4);
  for(i=1;i<=n;i++)
    for(j=i;j<=n;j++)
    {
      p10=gmod(gmul((GEN)p3[j],(GEN)p3[i]),x);
      j1=j+(i-1)*n;j2=i+(j-1)*n;lgp=lgef(p10)-2;
      for(k=1;k<=lgp;k++) coeff(p5,k,j1)=coeff(p5,k,j2)=p10[k+1];
      for(;k<=n;k++) coeff(p5,k,j1)=coeff(p5,k,j2)=zero;
      coeff(p11,i,j)=coeff(p11,j,i)=(long)trace9(p10,p12);
    }
  if(DEBUGLEVEL) {fprintferr("temps table de mult: ");fprintferr("%ld\n",timer2());flusherr();}
  p2=rootslong(x,PRECREG);
  if(DEBUGLEVEL) {fprintferr("temps racines: ");fprintferr("%ld\n",timer2());flusherr();}
  for(i=1;i<=r1;i++) p2[i]=lreal((GEN)p2[i]);
  tetpil=avma;y=cgetg(10,17);y[1]=lcopy(x);y[8]=(long)ginv(p1);
  y[9]=lmul((GEN)y[8],p5);p1=cgetg(3,17);
  p1[1]=lstoi(r1);p1[2]=lstoi(r2=(n-r1)>>1);y[2]=(long)p1;ru=r1+r2;
  y[3]=lcopy(fieldd);y[4]=lcopy(ind);
  p4=cgetg(n+1,19);
  for(j=1;j<=n;j++)
  {
    p6=cgetg(ru+1,18);p4[j]=(long)p6;
    for(i=1;i<=ru;i++) p6[i]=(long)gsubst((GEN)p3[j],v,(GEN)p2[max(i,(i<<1)-r1)]);
  }
  av2=avma;p7=cgetg(8,17);p7[1]=(long)p4;
  p5=cgetg(ru+1,19);
  for(j=1;j<=ru;j++)
  {
    p6=cgetg(n+1,18);p5[j]=(long)p6;
    for(i=1;i<=n;i++) 
      p6[i]=(j<=r1)?lconj(gcoeff(p4,j,i)):lmul2n(gconj(gcoeff(p4,j,i)),1);
  }
  p7[2]=(long)p5;p7[3]=(long)greal(gmul(p5,p4));
  p7[4]=lcopy(p11);p7[5]=zero;
  p7[6]=lmul(ginv(p11), fieldd);
         /* discriminant * (inverse of y[5][4] = p7[4]) */
  p7[7]=(long)hnfmod((GEN)p7[6],detint((GEN)p7[6]));
         /* Ideal basis for discriminant * (inverse of different) */
  av3=avma;
  y[5]=lpile(av2,av3,gcopy(p7));p4=cgetg(ru+1,17);y[6]=(long)p4;
  for(i=1;i<=ru;i++) p4[i]=lcopy((GEN)p2[max(i,(i<<1)-r1)]);
  if(DEBUGLEVEL) {fprintferr("temps matrices: ");fprintferr("%ld\n",timer2());flusherr();}
  y[7]=lcopy(p3);
  if(fldif)
  {
    parikeep=pari_randseed;((GEN)y[5])[5]=(long)differente(y,fa);
    pari_randseed=parikeep;
/*    
    av2=avma;p1=idealinv(y,(GEN)(((GEN)y[5])[7]));av3=avma;
    ((GEN)y[5])[5]=lpile(av2,av3,gmul(fieldd,p1)); */
    if(DEBUGLEVEL) {fprintferr("temps fin differente: ");fprintferr("%ld\n",timer2());flusherr();}
  }
  else ((GEN)y[5])[5]=zero;
  return gerepile(av,tetpil,y);
}

GEN
initalg(GEN x, long prec)
{
  return initalgall0(x,1,prec);
}

GEN
oldinitalgred2(GEN x, long prec)
{
  GEN p1,p2,p3,p4,p5,p6,p7,fieldd,polr,ptrace,dx,adx,polmax,s,a,dxn,adxn,sn,phimax,fa;
  long tx=typ(x),n=lgef(x)-3,i,j,av=avma,av1,v,k,r1,imax,numb,flc,tetpil;

  if((tx!=10)||(n<=0)) err(poltyper);
  p1=content(x);p1=gcmp1(p1) ? x: gdiv(x,content(x));
  for(k=2;k<=n+2;k++) if(typ((GEN)p1[k])!=1) err(impl,"general algebraic extension");
  if(gcmp1(p2=(GEN)p1[n+2])) x=p1;
  else
  {
    x=cgetg(n+3,10);x[1]=p1[1];x[n+2]=un;x[n+1]=p1[n+1];p3=p2;
    for(k=n;k>=2;k--)
    {
      x[k]=lmulii(p3,(GEN)p1[k]);
      if(k>2) p3=mulii(p2,p3);
    }
  }
  p1=factor(x);
  if(lgef(gcoeff(p1,1,1))!=n+3) err(redper1);
  r1=sturm(x);p4=allbase4(x,0,&fieldd,&fa);
  if(r1<n)
  {
    polr=rootslong(x,prec);p3=cgetg(n+1,19);
    for(i=1;i<=n;i++)
    {
      p1=cgetg(n+1,18);p3[i]=(long)p1;
      for(j=1;j<=n;j++)
	p1[j]=lsubst((GEN)p4[i],varn((GEN)p4[i]),(GEN)polr[j]);
    }
    p2=greal(gmul(gconj(gtrans(p3)),p3));
    p1=lllgram(p2,prec);
    for(s=gzero,i=1;i<=n;i++) s=gadd(s,gnorm((GEN)polr[i]));
  }
  else
  {
    ptrace=cgetg(n+1,17);ptrace[1]=lstoi(n);
    for(k=1;k<n;k++) 
    {
      p3=gmulsg(k,(GEN)x[n-k+2]);
      for(i=1;i<k;i++) p3=gadd(p3,gmul((GEN)x[n-i+2],(GEN)ptrace[k-i+1]));
      ptrace[k+1]=lneg(p3);
    }
    p2=cgetg(n+1,19);
    for(i=1;i<=n;i++)
    {
      p1=cgetg(n+1,18);p2[i]=(long)p1;
      for(j=1;j<i;j++) p1[j]=lcopy(gcoeff(p2,i,j));
      for(j=i;j<=n;j++)
      {
	p5=gres(gmul((GEN)p4[i],(GEN)p4[j]),x);p6=gzero;
	for(k=0;k<=lgef(p5)-3;k++) p6=gadd(p6,gmul((GEN)p5[k+2],(GEN)ptrace[k+1]));
	p1[j]=(long)p6;
      }
    }
    p1=lllgramint(p2);
    s=(n>1)?gsub(gsqr((GEN)x[n+1]),gmul2n((GEN)x[n],1)):gsqr((GEN)x[2]);
  }
  v=varn(x);dx=discsr(x);adx=absi(dx);polmax=x;imax=0;
  a=cgetg(n+1,18);for(i=1;i<=n;i++) a[i]=lmul(p4,(GEN)p1[i]);
  for(numb=0,i=1;i<=n;i++)
  {
    av1=avma;p3=gmodulcp((GEN)a[i],x);p7=content((GEN)p3[2]);
    if(gcmp1(p7)) p3=caract(p3,v);
    else
    {
      p3=caract(gdiv(p3,p7),v);
      p3=gmul(gpuigs(p7,lgef(p3)-3),gsubst(p3,v,gdiv(polx[v],p7)));
    }
    p5=ggcd(deriv(p3,v),p3);
    if(lgef(p5)==3)
    {
      dxn=discsr(p3);adxn=absi(dxn);flc=gcmp(adxn,adx);numb++;
      if(flc<=0)
      {
	if(r1<n) for(sn=gzero,j=1;j<=n;j++) sn=gadd(sn,gnorm(poleval((GEN)a[i],(GEN)polr[j])));
	else sn=(n>1)?gsub(gsqr((GEN)p3[n+1]),gmul2n((GEN)p3[n],1)):gsqr((GEN)p3[2]);
	if(flc<0) {dx=dxn;adx=adxn;s=sn;polmax=p3;imax=i;}
	else 
	{
	  flc=gcmp(sn,s);
	  if((flc<0)||((!flc)&&(gpolcomp(p3,polmax)<0)))
	  {dx=dxn;adx=adxn;s=sn;polmax=p3;imax=i;}
	}
      }
    }
  }
  if(!numb) err(polreder1);
  phimax=imax?(GEN)a[imax]:polx[v];
  j=n+1;
  while((j>=2)&&(!signe((GEN)polmax[j]))) j-=2;
  if((j>=2)&&(signe((GEN)polmax[j])>0))
  {
    if(polmax==x) polmax=gcopy(x);
    for(;j>=2;j-=2) setsigne((GEN)polmax[j],-signe((GEN)polmax[j]));
    phimax=gneg(phimax);
  }
  p2=gmodulcp(phimax,x);
  tetpil=avma;return gerepile(av,tetpil,polymodrecip(p2));
}

GEN
rootsof1(GEN nf)
{
  long av=avma,tetpil,N,ld,fl,k,i,j;
  GEN algun,p1,y,R1,d,list,w;

  nf=checknf(nf);  
  N=lgef((GEN)nf[1])-3;R1=(GEN)((GEN)nf[2])[1];
  algun=(GEN)((GEN)nf[8])[1];
  if(signe(R1)) 
  {
    avma=av;y=cgetg(3,17);y[1]=deux;y[2]=lneg(algun);
    return y;
  }
  y=minim((GEN)((GEN)nf[5])[3],N,1000);
  if(itos((GEN)y[2])!=N) err(rootofer2);
  w=(GEN)y[1];
  if(!cmpii(w,gdeux)) 
  {
    avma=av;y=cgetg(3,17);y[1]=deux;y[2]=lneg(algun);
    return y;
  }
  d=divisors(w);ld=lg(d)-1;
  list=(GEN)y[3];k=lg(list);fl=1;
  for(i=1;(i<k)&&fl;i++)
  {
    p1=(GEN)list[i];j=1;
    while((j<ld-1)&&(!gegal(element_pow(nf,p1,(GEN)d[j]),algun))) j++;
    if(j<ld-1) p1=gneg(p1);j=1;
    while((j<ld-1)&&(!gegal(element_pow(nf,p1,(GEN)d[j]),algun))) j++;
    fl=(j<ld-1);
  }
  if(fl) err(rootsof1er);
  tetpil=avma;y=cgetg(3,17);y[1]=lcopy(w);y[2]=lcopy(p1);
  return gerepile(av,tetpil,y);
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~								     ~*/
/*~		      COMPOSITUM OF TWO NUMBER FIELDS                ~*/
/*~								     ~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#define sign2(x) (x>=0?0:1)

GEN
compositum(GEN pol1, GEN pol2)
{
  long av=avma,tetpil,i,j,a,b,l,l1,l2,fl,c,degmax,a1;
  GEN m,sol,pro0,pro1,pro2,pro3,k1,p;

  if((typ(pol1)!=10)||(typ(pol2)!=10)) err(talker,"incorrect type in compositum");
  if((varn(pol1)!=0)||(varn(pol2)!=0)) err(talker,"not the x variable in compositum");
  l1=lgef(pol1)-3;l2=lgef(pol2)-3;
  if((l1<=0)||(l2<=0)) err(talker,"constant polynomial in compositum");
  l=l1*l2;m=cgetg(l+2,19);
  for(j=1;j<=l+1;j++) m[j]=lgetg(l+1,18);
  fl=1;sol=polun[0];degmax=0;
  for(c=1,a=1;(c<=l)&&fl;a=sign2(a)-a)
  {
    pro0=cgetg(3,10);pro0[1]=evallgef(3)+evalsigne(1);pro0[2]=(long)polun[MAXVARN];
    pro1=gadd(polx[0],gmulsg(a,polx[MAXVARN]));
    for(b=0;b<=l;b++)
    {
      for(i=0;i<l1;i++)
      {
	pro3=truecoeff(pro0,i);setvarn(pro3,0);
	pro2=gmod(pro3,pol2);
	for(j=0;j<l2;j++)	coeff(m,i*l2+j+1,l+1-b)=(long)truecoeff(pro2,j);
      }
      if(b<l) pro0=gmod(gmul(pro0,pro1),pol1);
    }
    k1=kerint(m);c++;
    for(a1=1;(a1<lg(k1))&&fl;a1++)
    {
      p=gtopoly((GEN)k1[a1],0);
      pro0=(GEN)factor(p)[1];
      for(i=1;(i<lg(pro0))&&fl;i++) 
	if(lgef((GEN)pro0[i])>lgef(sol)) 
	{
	  sol=(GEN)pro0[i];degmax=lgef(sol)-3;fl=(degmax<l);
	}
    }
  }
  tetpil=avma;return gerepile(av,tetpil,gdiv(sol,(GEN)sol[lgef(sol)-1]));
}
/*====================================================================
            Forme Normale d'Hermite avec permutations
====================================================================*/

GEN
hnfperm(GEN a)
{
  GEN b,U,c,l,perm,s,d,u,v,p,q,x,y,x1,y1;
  long r,t,i,j,j1,k,m,n,av=avma,tetpil;


  n=lg(a)-1;if(!n) return cgetg(1,19);
  m=lg((GEN)a[1])-1;
  c=cgeti(m+1);l=cgeti(n+1);perm=cgeti(m+1);
  for(i=1;i<=m;i++) c[i]=0;
  for(j=1;j<=n;j++) l[j]=0;
  U=idmat(n);b=gcopy(a);
/* U Matrice de Passage :  a*U=b tout au long algo */

  for(i=1;i<=m;i++) c[i]=0;
  for(r=0,k=1;k<=n;k++)
  {
    for(j=1;j<k;j++) if((t=l[j]))
    {
      y=gcoeff(b,t,k);

      if(signe(y))
      {
	p=gcoeff(b,t,j);
	d=bezout(p,y,&u,&v);
	x1=divii(p,d);y1=divii(y,d);
	for(i=1;i<=m;i++)
	{
	  x=gcoeff(b,i,j);y=gcoeff(b,i,k);
	  coeff(b,i,j)=(long)addii(mulii(u,x),mulii(v,y));
	  coeff(b,i,k)=(long)subii(mulii(x1,y),mulii(y1,x));
	}
	for(i=1;i<=n;i++)
	{
	  x=gcoeff(U,i,j);y=gcoeff(U,i,k);
	  coeff(U,i,j)=(long)addii(mulii(u,x),mulii(v,y));
	  coeff(U,i,k)=(long)subii(mulii(x1,y),mulii(y1,x));
	}
	for(j1=1;j1<j;j1++) if(l[j1])
	{
	  d=gcoeff(b,t,j);
	  q=dvmdii(gcoeff(b,t,j1),d,&s);
	  if(signe(s)<0) q=addsi(-1,q);
	  if(signe(q))
	  {
	    for(i=1;i<=m;i++)
	      coeff(b,i,j1)=(long)subii(gcoeff(b,i,j1),mulii(q,gcoeff(b,i,j)));
	    for(i=1;i<=n;i++)
	      coeff(U,i,j1)=(long)subii(gcoeff(U,i,j1),mulii(q,gcoeff(U,i,j)));
	  }
	}
      }
    }
    for(t=m;t&&((c[t])||(gcmp0(gcoeff(b,t,k))));t--);
    if(t)
    {
      p=gabs(gcoeff(b,t,k),0);
      for(i=t-1;i;i--)
	if(signe(gcoeff(b,i,k))&&(gcmp(p,q=gabs(gcoeff(b,i,k),0))>0)) {p=q;t=i;}
      perm[++r]=l[k]=t;c[t]=k;
      if(signe(gcoeff(b,t,k))<0)
      {
	for(i=1;i<=m;i++) coeff(b,i,k)= (long)negi(gcoeff(b,i,k));
	for(i=1;i<=n;i++) coeff(U,i,k)= (long)negi(gcoeff(U,i,k));
	p=gcoeff(b,t,k);
      }
      for(j=1;j<k;j++) if(l[j])
      {
	q=dvmdii(gcoeff(b,t,j),p,&s);
	if(signe(s)<0) q=addsi(-1,q);
	if(signe(q))
	{
	  for(i=1;i<=m;i++)
	    coeff(b,i,j)=(long)subii(gcoeff(b,i,j),mulii(q,gcoeff(b,i,k)));
	  for(i=1;i<=n;i++)
	    coeff(U,i,j)=(long)subii(gcoeff(U,i,j),mulii(q,gcoeff(U,i,k)));
	}
      }
    }
    else l[k]=0;
  }

/* On a :    a*U=b   ( matrice (m,n) )
   U  matrice (n,n) tq  |det(U)|=1    ( U dans GL(n) )
   Les colonnes de b  telles que l[j]<>0 : base de Im(a)  ( il y en a r )
   Les colonnes de U telles que l[j]=0  : base de Ker(a) ( il y en a n-r )   */
   

  tetpil=avma;
  y=cgetg(4,17);
  p=cgetg(r+1,19);u=cgetg(n+1,19);
  for(t=1,k=r,j=1;j<=n;j++) if(l[j]) 
  {
	/*    p[++k]=lcopy((GEN)b[j]);
	      k=0 au depart: pour Matrice hnf lignes non permutees */
    q=(GEN)(p[k]=lgetg(m+1,18));
    for(i=1;i<=m;i++) q[i]=lcopy(gcoeff(b,perm[m-i+1],j));
    u[k+n-r]=lcopy((GEN)U[j]);
    k--;
  }
  else u[t++]=lcopy((GEN)U[j]);
  y[1]=(long)p;y[2]=(long)u;
  q=(GEN)(y[3]=lgetg(m+1,17));
  for(i=1;i<=m;i++) q[m-i+1]=lstoi(perm[i]);
  return gerepile(av,tetpil,y);
}

/*====================================================================
	    Forme Normale d'Hermite (Version par colonnes 31/01/94)
====================================================================*/


GEN
hnfnew(GEN a)	
{
  GEN b,c,h,x,y,u,v,x1,y1,p,q,s,d,U;
  long m,n,r,i,j,j1,k,li,ii,rg,z,av=avma,av1,tetpil,lim,dec;

  n=lg(a)-1;if(!n) return cgetg(1,19);
  m=lg((GEN)a[1])-1;
  lim=(bot+avma)>>1;
  c=cgeti(m+1);h=cgeti(n+1);
  for(i=1;i<=m;i++) c[i]=0;
  for(j=1;j<=n;j++) h[j]=m;
  av1=avma;
  b=gcopy(a);U=idmat(n);
  for(r=n+1,li=m;li;li--)
  {
    j=1;
    do
    {
      for(i=h[j];i>li;i--)
	if(signe(y=gcoeff(b,i,j)))
	{
	  k=c[i];x=gcoeff(b,i,k);   /* annuler bij a l'aide de p=bik */
	  d=bezout(x,y,&u,&v);
	  if(DEBUGLEVEL>5) {outerr(u);outerr(v);}
	  x1=divii(x,d);y1=divii(y,d);
	  for(ii=1;ii<=i;ii++)
	  {
	    x=gcoeff(b,ii,k);y=gcoeff(b,ii,j);
	    coeff(b,ii,k)=(long)addii(mulii(u,x),mulii(v,y));
	    coeff(b,ii,j)=(long)subii(mulii(x1,y),mulii(y1,x));
	  }
	  for(ii=1;ii<=n;ii++)
	  {
	    x=gcoeff(U,ii,k);y=gcoeff(U,ii,j);
	    coeff(U,ii,k)=(long)addii(mulii(u,x),mulii(v,y));
	    coeff(U,ii,j)=(long)subii(mulii(x1,y),mulii(y1,x));
	  }
	  for(j1=k+1;j1<=n;j1++)
	  {
	    q=dvmdii(gcoeff(b,i,j1),d,&s);
	    if(signe(s)<0) q=addsi(-1,q);
	    if(signe(q))
	    {
	      for(ii=1;ii<=i;ii++)
		coeff(b,ii,j1)=(long)subii(gcoeff(b,ii,j1),mulii(q,gcoeff(b,ii,k)));
	      for(ii=1;ii<=n;ii++)
		coeff(U,ii,j1)=(long)subii(gcoeff(U,ii,j1),mulii(q,gcoeff(U,ii,k)));
	    }
	  }
	}
      if(avma<lim)
      {
	tetpil=avma;b=gcopy(b);U=gcopy(U);
	dec=lpile(av1,tetpil,0)>>TWOPOTBYTES_IN_LONG;b+=dec;U+=dec;
      }	      
      x=gcoeff(b,li,j);if(!signe(x)) {h[j]=li-1;j++;}
    }
    while((j<r)&&(!signe(x)));
    if(j<r)
    {
      r--;
      if(j<r)
      {
	z=b[j];b[j]=b[r];b[r]=z;
	z=U[j];U[j]=U[r];U[r]=z;
	h[j]=h[r];h[r]=li;c[li]=r;
      }
      if(signe(gcoeff(b,li,r))<0)
      {
	for(i=1;i<=li;i++) coeff(b,i,r)=(long)negi(gcoeff(b,i,r));
	for(i=1;i<=n;i++) coeff(U,i,r)=(long)negi(gcoeff(U,i,r));
      }
      p=gcoeff(b,li,r);
      for(j=r+1;j<=n;j++)
      {
	q=dvmdii(gcoeff(b,li,j),p,&s);
	if(signe(s)<0) q=addsi(-1,q);
	if(signe(q)) 
	{
	  for(i=1;i<=li;i++)
	    coeff(b,i,j)=(long)subii(gcoeff(b,i,j),mulii(q,gcoeff(b,i,r)));
	  for(i=1;i<=n;i++)
	    coeff(U,i,j)=(long)subii(gcoeff(U,i,j),mulii(q,gcoeff(U,i,r)));
	}
      }
    }
  }
  r--;
  for(j=1;j<=r;j++)   /* les r 1ers vec sont dans l'espace engendre
			 par les n-r derniers (base sur R de im(a)) */
    for(i=h[j];i;i--)
      if(signe(y=gcoeff(b,i,j)))
      {
	k=c[i];x=gcoeff(b,i,k);
	d=bezout(x,y,&u,&v);
	if(DEBUGLEVEL>5) {outerr(u);outerr(v);}
	x1=divii(x,d);y1=divii(y,d);
	for(ii=1;ii<=i;ii++)
	{
	  x=gcoeff(b,ii,k);y=gcoeff(b,ii,j);
	  coeff(b,ii,k)=(long)addii(mulii(u,x),mulii(v,y));
	  coeff(b,ii,j)=(long)subii(mulii(x1,y),mulii(y1,x));
	}
	for(ii=1;ii<=n;ii++)
	{
	  x=gcoeff(U,ii,k);y=gcoeff(U,ii,j);
	  coeff(U,ii,k)=(long)addii(mulii(u,x),mulii(v,y));
	  coeff(U,ii,j)=(long)subii(mulii(x1,y),mulii(y1,x));
	}
	for(j1=k+1;j1<=n;j1++)
	{
	  q=dvmdii(gcoeff(b,i,j1),d,&s);
	  if(signe(s)<0) q=addsi(-1,q);
	  if(signe(q))
	  {
	    for(ii=1;ii<=i;ii++)
	      coeff(b,ii,j1)=(long)subii(gcoeff(b,ii,j1),mulii(q,gcoeff(b,ii,k)));
	    for(ii=1;ii<=n;ii++)
	      coeff(U,ii,j1)=(long)subii(gcoeff(U,ii,j1),mulii(q,gcoeff(U,ii,k)));
	  }
	}
      }
  rg=n-r;
  tetpil=avma;
  y=cgetg(3,17);
  h=cgetg(rg+1,19);
  for(j=1;j<=rg;j++) h[j]=lcopy((GEN)b[r+j]);
  y[1]=(long)h;y[2]=lcopy(U);
  return gerepile(av,tetpil,y);
}
	  
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                    FONCTION ZETA DE DEDEKIND                    */
/*                                                                 */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#define T ((GEN)(nf[1]))
#define D ((GEN)(nf[3]))

GEN
initzeta(GEN pol, long prec)
{
  GEN nf,p_i,tab_prem,alpha,beta,mu,gr1,gr2;
  GEN cst,A0,c0,c1,c2,eps,lneps,limx,bg,ballg,resi,zet;
  GEN c_pair,ck_pair,c_impair,ck_impair;
  GEN serie_pair,serie_impair,serie_exp,cik,vectf;
  GEN coeflog,cstlog,logpro,racpi,gN0,pg,borneps,ftest;
  GEN aij,tabj,tabcstn,tabcstni,nfz,pnfz,pnfz2,tmp;
  long imin,itest,imax,r1,r2,N,i,j,k,n,av,tetpil,N0,lv;
  long dec,expolim,limk,f_j,p,q,limpile,prec2;
  long* coef;
  long* coef2;
  byteptr ptpi;

  expolim=-((prec-2)<<TWOPOTBITS_IN_LONG)-6;
  eps=gmul2n(gun,expolim);
  prec2=prec+1;
  prec=(prec<<1)-1;
  constpi(prec);consteuler(prec);
  racpi=gsqrt(gpi,prec);

      /*************** Calcul du residu et des constantes ***************/

      /* Nb de classes et regulateur */

  ballg=buchinit(pol,dbltor(0.5),dbltor(0.5),prec2);
  bg=(GEN)ballg[8];nf=(GEN)ballg[7];N=lgef(T)-3;
  av=avma;
  gr1=(GEN)((GEN)nf[2])[1];gr2=(GEN)((GEN)nf[2])[2];
  r1=itos(gr1);r2=itos(gr2);
  resi=gdiv(gmul(gmul(gpuigs(gdeux,r1),(GEN)((GEN)bg[1])[1]),
		 (GEN)bg[2]),(GEN)((GEN)bg[4])[1]);

      /* Calcul de N0 */

  mu=gadd(gmul2n(gr1,-1),gr2);
  alpha=gmul2n(stoi(r1+r2+3),-1);
  beta=gpui(gdeux,gmul2n(gr1,-1),DEFAULTPREC);
  A0=gmul2n(gun,r1);
  A0=gmul(A0,gpuigs(mu,r1+r2+2));
  A0=gmul(A0,gpuigs(gmul2n(gpi,1),1-r1-r2));
  A0=gsqrt(A0,DEFAULTPREC);
  c1=gmul(mu,gpui(beta,ginv(mu),DEFAULTPREC));
  c0=gdiv(gmul(A0,gpuigs(gmul(gdeux,gpi),r1+r2-1)),mu);
  c0=gmul(c0,gpui(c1,gsub(gun,alpha),DEFAULTPREC));
  c2=gdiv(gsub(alpha,gun),mu);
  lneps=glog(gdiv(c0,eps),DEFAULTPREC);
  limx=gdiv(gsub(glog(lneps,DEFAULTPREC),glog(c1,DEFAULTPREC)),gadd(c2,gdiv(lneps,mu)));
  limx=gmul(gpui(gdiv(c1,lneps),mu,DEFAULTPREC),gadd(gun,gmul(c2,gmul(mu,limx))));
  cst=gdiv(gsqrt(gabs(D,0),prec),gmul(gpuigs(racpi,N),gpuigs(gdeux,r2)));
  gN0=gfloor(gdiv(cst,limx));
  if(cmpis(gN0,MAXHALFULONG)>=0)
    err(talker,"discriminant too large for initzeta, sorry");
  N0=itos(gN0);
  if(DEBUGLEVEL>=2) 
  {
    fprintferr("N0 = %ld\n",N0);
    flusherr();
  }

      /* Calcul de imax */

  imin=1;imax=1400;
  borneps=gmul(gpuigs(gmul2n(racpi,1),r2),gpuigs(stoi(5),r1));
  borneps=gdiv(borneps,gmul(gmul(gsqrt(limx,DEFAULTPREC),gmul2n(eps,4)),
			    gpuigs(racpi,3)));
  while((imax-imin)>=4)
  {
    itest=(imax+imin)>>1;
    ftest=gmul(gpuigs(mpfactr(itest,DEFAULTPREC),r2),gpuigs(limx,itest));
    ftest=gmul(ftest,gpuigs(mpfactr(itest/2,DEFAULTPREC),r1));
    if(gcmp(ftest,borneps)>=0) imax=itest;
    else imin=itest;
  }
  imax=(imax/2)*2;
  if(DEBUGLEVEL>=2) 
  {
    fprintferr("imax = %ld\n",imax);
    flusherr();
  }

      /* Tableau des i/cst (i=1 a N0) */

  cstlog=glog(cst,prec);
  tabcstn=cgetg(N0+1,17);
  tabcstni=cgetg(N0+1,17);
  for(i=1;i<=N0;i++) 
  {
    tabcstn[i]=ldivsg(i,cst);
    tabcstni[i]=un;
  }
  tetpil=avma;
  resi=gcopy(resi);cst=gcopy(cst);cstlog=gcopy(cstlog);
  tabcstn=gcopy(tabcstn);tabcstni=gcopy(tabcstni);
  dec=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;
  resi+=dec;cst+=dec;cstlog+=dec;tabcstn+=dec;tabcstni+=dec;

      /********** Calcul des coefficients a(i,j) independants de s **********/

  if(DEBUGLEVEL>=2) {fprintferr("Calcul des a(i,j)\n");flusherr();}
  av=avma;
  zet=cgetg(r1+r2+2,17);
  zet[1]=(long)geuler;
  for(i=2;i<=r1+r2+1;i++)
    zet[i]=(long)gzeta(stoi(i),prec);
  aij=cgetg(imax+1,17);
  for(i=1;i<=imax;i++)
    aij[i]=lgetg(r1+r2+2,17);
  affsr(1,c_pair=cgetr(prec));c_pair=gmul2n(c_pair,r1);
  c_impair=gmul(c_pair,gpuigs(racpi,r1));
  if((r1+r2)&1) c_impair=gneg(c_impair);
  ck_pair=cgetg(r1+r2+2,17);ck_impair=cgetg(r2+2,17);
  for(k=1;k<=r1+r2+1;k++)
  {
    ck_pair[k]=lmul((GEN)zet[k],gadd(gr2,gdiv(gr1,gpuigs(gdeux,k))));
    if(k&1) ck_pair[k]=lneg((GEN)ck_pair[k]);
  }
  for(k=1;k<=r2+1;k++)
  {
    ck_impair[k]=lmul((GEN)zet[k],gadd(gr2,gmul(gr1,gsub(gun,gmul2n(gun,-k)))));
    if(k&1) ck_impair[k]=lneg((GEN)ck_impair[k]);
    ck_impair[k]=laddsg(r1+r2,(GEN)ck_impair[k]);
  }
  ck_impair[1]=lsub((GEN)ck_impair[1],gmul(gr1,glog(gdeux,prec)));
  serie_pair=cgetg(r1+r2+3,11);
  serie_pair[1]=evalsigne(1)+evalvalp(1);
  serie_impair=cgetg(r2+3,11);
  serie_impair[1]=evalsigne(1)+evalvalp(1);
  i=0;
  while(i<imax/2)
  {
    for(k=1;k<=r1+r2+1;k++)
      serie_pair[k+1]=ldivgs((GEN)ck_pair[k],k);
    serie_exp=gmul(c_pair,gexp(serie_pair,prec));
    for(j=1;j<=r1+r2+1;j++)
      ((GEN)aij[2*i+1])[j]=serie_exp[r1+r2+3-j];
    for(k=1;k<=r2+1;k++)
      serie_impair[k+1]=ldivgs((GEN)ck_impair[k],k);
    serie_exp=gmul(c_impair,gexp(serie_impair,prec));
    for(j=1;j<=r2+1;j++)
      ((GEN)aij[2*i+2])[j]=serie_exp[r2+3-j];
    for(j=r2+2;j<=r1+r2+1;j++)
      ((GEN)aij[2*i+2])[j]=zero;
    i++;
    c_pair=gdiv(c_pair,gmul(gpuigs(stoi(i),r1+r2),gpuigs(stoi(4*i-2),r2)));
    if(r1&1) c_pair=gneg(c_pair);
    c_impair=gdiv(c_impair,gmul(gpuigs(stoi(i),r2),gpuigs(stoi(2*i+1),r1+r2)));
    c_impair=gmul2n(c_impair,r1-r2);if(r1&1) c_impair=gneg(c_impair);
    for(k=1;k<=r1+r2+1;k++)
      ck_pair[k]=ladd((GEN)ck_pair[k],gadd(gdiv(gr2,gpuigs(stoi(2*i-1),k)),
					   gdiv(stoi(r1+r2),gpuigs(stoi(2*i),k))));
    for(k=1;k<=r2+1;k++)
      ck_impair[k]=ladd((GEN)ck_impair[k],gsub(gmul(stoi(r1+r2),
						    gadd(ginv(gpuigs(stoi(2*i),k)),ginv(gpuigs(stoi(2*i+1),
											       k)))),gdiv(gr1,gpuigs(stoi(2*i),k))));
  }
  tetpil=avma;aij=gerepile(av,tetpil,gcopy(aij));
  nfz=cgetg(10,17);pnfz=cgetg(5,17);nfz[1]=(long)pnfz;nfz[2]=(long)resi;
  nfz[5]=(long)cst;nfz[6]=(long)cstlog;nfz[7]=(long)aij;pnfz[1]=lstoi(r1);
  pnfz[2]=lstoi(r2);pnfz[3]=lstoi(imax);pnfz[4]=(long)ballg;

      /************* Calcul du nombre d'ideaux de norme donnee *************/

  if(DEBUGLEVEL>=2) {fprintferr("Calcul des a(n)\n");flusherr();}
  av=avma;
  ptpi=diffptr;gN0=stoi(N0);
  coef=newbloc(N0+1);coef[0]=evaltyp(2)+evalpere(MAXPERE)+evallg(N0+1);
  coef2=newbloc(N0+1);coef2[0]=evaltyp(2)+evalpere(MAXPERE)+evallg(N0+1);
      /* coef[0] et coef2[0] necessaires pour la fonction taille */
  for(i=2;i<=N0;i++) coef[i]=0;
  coef[1]=1;
  p_i=cgeti(3);p_i[1]=evalsigne(1)+evallgef(3);
  p_i[2]=*ptpi++;
  while((*ptpi)&&(p_i[2]<=N0))
  {
    if(signe(modii((GEN)nf[4],p_i))) 
    {vectf=(GEN)simplefactmod(T,p_i)[1];lv=lg(vectf);}
    else
    {
      tab_prem=primedec(nf,p_i);lv=lg(tab_prem);
      vectf=cgetg(lv,18);
      for(j=1;j<lv;j++) vectf[j]=((GEN)tab_prem[j])[4];
    }
    for(j=1;j<lv;j++)
    {
      f_j=itos((GEN)vectf[j]);
      pg=gpuigs(p_i,f_j);
      if(gcmp(pg,gN0)<=0)
      {
	p=itos(pg);
	for(k=1;k<=N0;k++) coef2[k]=coef[k];
	q=p;limk=N0/q;
	while(q<=N0)
	{
	  for(k=1;k<=limk;k++) coef2[k*q]+=coef[k];
	  q*=p;limk/=p;
	}
	tmp=coef;coef=coef2;coef2=tmp;
      }
    }
    p_i[2]+=*ptpi++;
  }
  if(!(*ptpi)) err(primer1);
  coeflog=cgetg(N0+1,17);for(i=1;i<=N0;i++) coeflog[i]=zero;
  for(i=2;i<=N0;i++)
    if(coef[i]!=0)
      coeflog[i]=(long)glog(stoi(i),prec);
  tabj=cgetg(N0+1,19);
  for(i=1;i<=N0;i++) 
    tabj[i]=lgetg(r1+r2+2,18);    /* tabj[n,j]=coef(n)*ln(c/n)^(j-1)/(j-1)! */
  for(i=1;i<=N0;i++)
  {
    if(coef[i]!=0)
    {
      coeff(tabj,1,i)=un;logpro=gneg(glog((GEN)tabcstn[i],prec));
      for(j=2;j<=r1+r2+1;j++)
	coeff(tabj,j,i)=ldivgs(gmul(gcoeff(tabj,j-1,i),logpro),j-1);
      tabj[i]=lmulgs((GEN)tabj[i],coef[i]);
    }
    else
      for(j=1;j<=r1+r2+1;j++) coeff(tabj,j,i)=zero;
  }
  tetpil=avma;coeflog=gneg(coeflog);tabj=gcopy(tabj);
  dec=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;
  coeflog+=dec;tabj+=dec;nfz[3]=(long)tabj;
  pnfz2=cgetg(N0+1,17);for(i=1;i<=N0;i++) pnfz2[i]=lstoi(coef[i]);
  nfz[8]=(long)pnfz2;nfz[9]=(long)coeflog;
  killbloc(coef);killbloc(coef2);

      /******************** Calcul des coefficients Cik ********************/

  if(DEBUGLEVEL>=2) {fprintferr("Calcul des Cik\n");flusherr();}
  av=avma;limpile=(avma+bot)>>1;
  cik=cgetg(r1+r2+1,19);
  for(k=1;k<=r1+r2;k++) 
  {
    cik[k]=lgetg(imax+1,18);
    for(i=1;i<=imax;i++) coeff(cik,i,k)=zero;
  }      
  for(i=1;i<=imax;i++)
  {
    for(k=0;k<=r1+r2-1;k++)
    {
      for(n=N0;n>=1;n--)
	if(signe(gcoeff(tabj,1,n)))
	  for(j=1;j<=r1+r2-k;j++)
	    coeff(cik,i,k+1)=ladd(gcoeff(cik,i,k+1),gmul(gmul((GEN)((GEN)aij[i])[j+k+1],gcoeff(tabj,j,n)),(GEN)tabcstni[n]));
    }
    for(n=1;n<=N0;n++) tabcstni[n]=lmul((GEN)tabcstni[n],(GEN)tabcstn[n]);
    if(avma<limpile)
    {
      tetpil=avma;cik=gcopy(cik);tabcstni=gcopy(tabcstni);
      dec=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;cik+=dec;tabcstni+=dec;
    }
  }
  tetpil=avma;cik=gerepile(av,tetpil,gcopy(cik));
  nfz[4]=(long)cik;return nfz;
}

GEN
gzetakall(GEN nfz, GEN s, long flag, long prec)

{
  GEN pnfz,resi,tabj,cik,cst,cstlog,aij,coeflog,varint,cs,gcoef;
  GEN lambd,zetak,gammas,gammaunmoins,gammas2,gammaunmoins2,var1,var2;
  GEN unmoins,gexpro,gar,nf,val,valm,valk,valkm,val2,valm2,valk2,valkm2;
  long i,j,k,N,N0,av,tetpil,sl;
  long r1,r2,imax,precbig;

  if((typ(nfz)!=17)||(lg(nfz)!=10))
    err(talker,"not a zeta number field in zetakall");
  precbig=(prec<<1)-1;prec++;
  s=gprec(s,(precbig-2)*10);
  if ((typ(s)==2)&&(!signe(gfrac(s)))) s=gtrunc(s);
  resi=(GEN)nfz[2];tabj=(GEN)nfz[3];cik=(GEN)nfz[4];cst=(GEN)nfz[5];
  cstlog=(GEN)nfz[6];aij=(GEN)nfz[7];gcoef=(GEN)nfz[8];N0=lg(gcoef)-1;
  coeflog=(GEN)nfz[9];pnfz=(GEN)nfz[1];
  r1=itos((GEN)pnfz[1]);r2=itos((GEN)pnfz[2]);N=r1+(r2<<1);
  imax=itos((GEN)pnfz[3]);nf=(GEN)((GEN)pnfz[4])[7];
  av=avma;
  unmoins=gsub(gun,s);
  if(typ(s)==1)
  {
    sl=itos(s);
    if((sl<0)&&((!signe(modii(s,gdeux)))||(r2!=0))) return(gzero);
    if(sl==0)
    {
      if((r1==1)&&(r2==0)) return(gneg(ghalf));
      if((r1==0)&&(r2==1)) return(gneg(resi));
      if(!((r1==1)&&(r2==0))&&!((r1==0)&&(r2==1))) return(gzero);
    }
    if((sl<0)&&(r2==0))
    {
      lambd=gdiv(resi,gmul(s,gsub(s,gun)));
      gammas2=ggamma(gdiv(s,gdeux),prec);gar=gpuigs(gammas2,r1);
      cs=gexp(gmul(cstlog,s),prec);
      gammaunmoins2=ggamma(gdiv(unmoins,gdeux),prec);
      var1=var2=gun;
      for(i=2;i<=N0;i++)
	if(signe((GEN)gcoef[i]))
	{
	  var1=gadd(var1,gmulsg(((GEN)gcoef[i])[2],
				gexpro=gexp(gmul((GEN)coeflog[i],s),prec)));
	  var2=gadd(var2,gdivsg(((GEN)gcoef[i])[2],gmulsg(i,gexpro)));
	}
      lambd=gadd(lambd,gmul(gmul(var1,cs),gar));
      lambd=gadd(lambd,gmul(gmul(var2,gdiv(cst,cs)),
			    gpuigs(gammaunmoins2,r1)));
      val=gsub(s,gdeux);valm=gsub(unmoins,gdeux);
      for(i=1;i<=imax;i+=2)
      {
	val=gadd(val,gdeux);valm=gadd(valm,gdeux);
	valk=val;valkm=valm;
	for(k=0;k<=r1+r2-1;k++)	
	{
	  lambd=gsub(lambd,gdiv(gcoeff(cik,i,k+1),valk));
	  lambd=gsub(lambd,gdiv(gcoeff(cik,i,k+1),valkm));
	  valk=gmul(val,valk);valkm=gmul(valm,valkm);
	}
      }
      zetak=flag?lambd:gdiv(lambd,gmul(gar,cs));
    }
    else
    {
      lambd=gdiv(resi,gmul(s,gsub(s,gun)));
      gammas=ggamma(s,prec);gammas2=ggamma(gdiv(s,gdeux),prec);
      gar=gmul(gpuigs(gammas,r2),gpuigs(gammas2,r1));
      cs=gexp(gmul(cstlog,s),prec);	  
      var1=var2=gzero;
      for(i=1;i<=N0;i++)
	if(signe((GEN)gcoef[i]))
	{
	  varint=gzero;
	  var1=gadd(var1,gmulsg(((GEN)gcoef[i])[2],
				gexpro=gexp(gmul((GEN)coeflog[i],s),prec)));
	  for(j=1;j<=r1+r2+1;j++)
	    varint=gadd(varint,gmul((GEN)((GEN)aij[sl])[j],
				    gcoeff(tabj,j,i)));
	  var2=gadd(var2,gdiv(varint,gmulsg(i,gexpro)));
	}
      lambd=gadd(lambd,gmul(gmul(var1,cs),gar));
      lambd=gadd(lambd,gmul(var2,gdiv(cst,cs)));
      val=s;valm=unmoins;
      for(i=1;i<=imax;i++)
      {
	valk=val;valkm=valm;
	for(k=0;k<=r1+r2-1;k++)
	{	      
	  lambd=gsub(lambd,gdiv(gcoeff(cik,i,k+1),valk));
	  if(i!=sl)
	    lambd=gsub(lambd,gdiv(gcoeff(cik,i,k+1),valkm));
	  valk=gmul(val,valk);valkm=gmul(valm,valkm);
	}
	val=gadd(val,gun);valm=gadd(valm,gun);
      }
      zetak=flag?lambd:gdiv(lambd,gmul(gar,cs));
    }
  }
  else
  {
    lambd=gdiv(resi,gmul(s,gsub(s,gun)));
    gammas=ggamma(s,prec);gammas2=ggamma(gmul2n(s,-1),prec);
    gar=gmul(gpuigs(gammas,r2),gpuigs(gammas2,r1));
    cs=gexp(gmul(cstlog,s),prec);
    var1=gmul(gpi,s);
    gammaunmoins=gdiv(gpi,gmul(gsin(var1,prec),gammas));
    gammaunmoins2=gdiv(gmul(gmul(gsqrt(gpi,prec),gpui(gdeux,gsub(s,gun),prec)),gammas2),gmul(gcos(gmul2n(var1,-1),prec),gammas));
    var1=var2=gun;
    for(i=2;i<=N0;i++)
      if(signe((GEN)gcoef[i]))
      {
	var1=gadd(var1,gmulsg(((GEN)gcoef[i])[2],
			      gexpro=gexp(gmul((GEN)coeflog[i],s),prec)));
	var2=gadd(var2,gdivsg(((GEN)gcoef[i])[2],gmulsg(i,gexpro)));
      }
    lambd=gadd(lambd,gmul(gmul(var1,cs),gar));
    lambd=gadd(lambd,gmul(gmul(gmul(var2,gdiv(cst,cs)),
			       gpuigs(gammaunmoins,r2)),gpuigs(gammaunmoins2,r1)));
    val=gsub(s,gdeux);valm=gsub(unmoins,gdeux);
    for(i=1;i<=imax;i+=2)
    {
      val=gadd(val,gdeux);valm=gadd(valm,gdeux);
      valk=val;valkm=valm;
      if(r2!=0)
      {
	val2=gadd(val,gun);valm2=gadd(valm,gun);
	valk2=val2;valkm2=valm2;
      } 
      for(k=0;k<=r1+r2-1;k++)
      {
	lambd=gsub(lambd,gdiv(gcoeff(cik,i,k+1),valk));
	lambd=gsub(lambd,gdiv(gcoeff(cik,i,k+1),valkm));
	valk=gmul(val,valk);valkm=gmul(valm,valkm);
      }
      if(r2!=0)
	for(k=0;k<=r1+r2-1;k++)
	{
	  lambd=gsub(lambd,gdiv(gcoeff(cik,i+1,k+1),valk2));
	  lambd=gsub(lambd,gdiv(gcoeff(cik,i+1,k+1),valkm2));
	  valk2=gmul(val2,valk2);valkm2=gmul(valm2,valkm2);
	}
    }
    zetak=flag?lambd:gdiv(lambd,gmul(gar,cs));
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(zetak));
}

GEN
gzetak(GEN nfz, GEN s, long prec)
{
  return gzetakall(nfz,s,0,prec);
}

GEN
glambdak(GEN nfz, GEN s, long prec)
{
  return gzetakall(nfz,s,1,prec);
}


