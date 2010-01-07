/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@                                                               @*/
/*@             PROGRAMMES D'ENTREES-SORTIES DES GEN              @*/
/*@                                                               @*/
/*@                      copyright Babe Cool                      @*/
/*@                                                               @*/
/*@                                                               @*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

/*******************************************************************/
/*******************************************************************/
/*                                                                 */
/*                 LISTE DES TYPES GENERIQUES                      */
/*                 ~~~~~~~~~~~~~~~~~~~~~~~~~~                      */
/*                                                                 */
/*  1  :entier long     [ cod1 ] [ cod2 ] [ man1 ] ... [ manl ]    */
/*  2  :reel            [ cod1 ] [ cod2 ] [ man1 ] ... [ manl ]    */
/*  3  :entier modulo   [ code ] [ mod  ] [ entier modulo ]        */
/*  4  :fraction        [ code ] [ num. ] [ den. ]                 */
/*  5  :nfraction       [ code ] [ num. ] [ den. ]                 */
/*  6  :complexe        [ code ] [ reel ] [ imag ]                 */
/*  7  :p-adique        [ cod1 ] [ cod2 ] [ p ] [ p^r ] [ entier]  */
/*  8  :quadrat         [ cod1 ] [ mod  ] [ reel ] [ imag ]        */
/*  9  :poly mod        [ code ] [ mod  ] [ polynome  mod ]        */
/* --------------------------------------------------------------- */
/*  10 :polynome        [ cod1 ] [ cod2 ] [ man1 ] ... [ manl ]    */
/*  11 :serie           [ cod1 ] [ cod2 ] [ man1 ] ... [ manl ]    */
/*  13 :fr.rat          [ code ] [ num. ] [ den. ]                 */
/*  14 :n.fr.rat        [ code ] [ num. ] [ den. ]                 */
/*  16 :forme quadrat   [ code ] [  a   ] [  b   ] [  c   ]        */
/*  17 :vecteur ligne   [ code ] [  x1  ] ... [  xl  ]             */
/*  18 :vecteur colonne [ code ] [  x1  ] ... [  xl  ]             */
/*  19 :matrice         [ code ] [ col1 ] ... [ coll ]             */
/*                                                                 */
/*******************************************************************/
/*******************************************************************/

# include       "genpari.h"
static void monome(long v, long deg),texnome(long v, long deg);

/********************************************************************/
/********************************************************************/
/**                                                                **/
/**                     FILTRAGE D'ENTREE                          **/
/**                                                                **/
/********************************************************************/
/********************************************************************/

void
filtre(char *s)
{
  char c, *s1 = s;
  int outer = 1;
  
  while((c = *s++))
  {
    if (outer)
      if (isspace(c)) continue; else *s1++ = isupper(c) ? tolower(c) : c;
    else
      *s1++ = c;
    if (c == '"') outer = !outer;
  }
  *s1 = 0;
}

/********************************************************************/
/********************************************************************/
/**                                                                **/
/**             UTILITAIRES GENERAUX D'IMPRESSION                  **/
/**                                                                **/
/********************************************************************/
/********************************************************************/

PariOUT *pariOut;
PariERR *pariErr;

typedef struct outString {
  char *s;
  uLong len;
  uLong size;
} outString;

outString pariOutStr;

#define INISIZE 1024
#define STEPSIZE 1024

void
outstrAppC(char c)
{
  if (pariOutStr.len+1 >= pariOutStr.size) {
#if __MWERKS__
    pariOutStr.s = (char*) macrealloc(pariOutStr.s,pariOutStr.size,pariOutStr.size+STEPSIZE);
    pariOutStr.size+=STEPSIZE;
#else
    pariOutStr.s = (char*) realloc(pariOutStr.s,pariOutStr.size+=STEPSIZE);
#endif
    if (!pariOutStr.s) {
      err(memer);
    }
  }
  pariOutStr.s[pariOutStr.len]=c;
  pariOutStr.s[++pariOutStr.len]=0;
}


void
outstrAppS(char* s)
{
  long len=strlen(s);
  if (pariOutStr.len+len >= pariOutStr.size) {
#if __MWERKS__
    pariOutStr.s = (char*) macrealloc(pariOutStr.s,pariOutStr.size,pariOutStr.size+STEPSIZE+len);
    pariOutStr.size+=STEPSIZE+len;
#else
    pariOutStr.s = (char*)realloc(pariOutStr.s, pariOutStr.size+=STEPSIZE+len);
#endif
    if (!pariOutStr.s) {
      err(memer);
    }
  }
  strcpy(pariOutStr.s+pariOutStr.len,s);
  pariOutStr.len+=len;
}

PariOUT pariOut2Str = {outstrAppC, outstrAppS};

/* This function is non-reenterable,
 * but as far as I can understand,
 * this should not be a problem */

char *
gen2str(GEN x)
{
  pariOut=&pariOut2Str;
  pariOutStr.s = (char*)malloc(INISIZE);
  if (!pariOutStr.s) {
    err(memer);
  }
  pariOutStr.s[0]=0;
  pariOutStr.len=0;
  pariOutStr.size=INISIZE;
  brute(x,'g',-1);
  pariOut=NULL;
  return pariOutStr.s;
}

void
pariputc(char c)
{
  if ((errfile==outfile? pariErr && pariErr->putc: pariOut && pariOut->putc)) {
    (errfile==outfile? pariErr->putc: pariOut->putc)(c);
    return;
  }
  putc(c, outfile);
  if (logfile) putc(c, logfile);
}

void
pariputs(char *s)
{
  if ((errfile==outfile? pariErr && pariErr->puts: pariOut && pariOut->puts)) {
    (errfile==outfile? pariErr->puts: pariOut->puts)(s);
    return;
  }
  fputs(s, outfile);
  if (logfile) fputs(s, logfile);
}

static void
blancs(long nb)
{
  while(nb-->0) pariputc(' ');
}

static void
zeros(long nb)
{
  while(nb-->0) pariputc('0');
}

static long
coinit(long grandmot)
{
  char cha[10], *p = cha + 9;
  *p = 0;
  do {*--p = grandmot%10 + '0'; grandmot /= 10;} while (grandmot);
  pariputs(p);
  return (cha - p) + 9;
}

static void
comilieu(long grandmot)
{
  char cha[10], *p = cha + 9;
  
  for(*p = 0; p > cha; grandmot /= 10) *--p = grandmot%10 + '0';
  pariputs(cha);
}

static void
cofin(long grandmot, long decim)
{
  char cha[10], *p = cha + 9;
  for(; p > cha; grandmot /= 10) *--p = grandmot%10 + '0';
  cha[decim] = 0;
  pariputs(cha);
}

static long
nbdch(long l)
{
  if (l<10) return 1;
  if (l<100) return 2;
  if (l<1000) return 3;
  if (l<10000) return 4;
  if (l<100000) return 5;
  if (l<1000000) return 6;
  if (l<10000000) return 7;
  if (l<100000000) return 8;
  if (l<1000000000) return 9;
  return 10; /* ne doit pas se produire */
}

/********************************************************************/
/********************************************************************/
/**                                                                **/
/**                     ECRIRE UN NOMBRE                           **/
/**                                                                **/
/********************************************************************/
/********************************************************************/


void
ecrireall(GEN x, char format, long dec, long chmp, long flbl)
{
  int     typy,sgn,i;
  GEN     enti,frac,modifie,p1,dix;
  long    avmacourant,d ,longueur,e,f,ex;
  long    nbch,*res,*re,decmax,deceff,arrondi[3];
  char    thestring[20];
  
  typy=typ(x);
  sgn=signe(x);
  
  if (typy==1)
	/* ecriture d'un entier */
  {
    if (! sgn) pariputc('0');
    else
    {
      re=res=(long *)convi(x);
      nbch=nbdch(*--re);
      while (*--re!= -1) nbch+=9;
      if (sgn!=1) nbch++;
      blancs(chmp-nbch);
      if (sgn!=1) pariputc('-');
      coinit(*--res);
      while (*--res!= -1) comilieu(*res);
    }
  }
  else
	/* ecriture d'un reel */
    switch (format)
    {
      case 'f':
        if (! sgn)
	      /*  reel 0 */
	{
	  pariputs("0.");
	  longueur=1+((-expo(x))>>TWOPOTBITS_IN_LONG);
	  if (longueur<0) longueur=0;
	  if (dec<0) dec=(long)(K*longueur);
	  zeros(dec);
	}
        else
	      /* reel non nul */
	{
	  avmacourant=avma;
	  if (sgn!=1) pariputc('-');
	      /*  on arrondit si il y a lieu */
	  for (i=0;i<=2;i++) arrondi[i]=x[i];
	  setlg(arrondi,3);
	  if (dec>0)
	  {
	    arrondi[1]=(long)(arrondi[1]-((double)BITS_IN_LONG/K)*dec-2);
	    modifie=mpadd(x,arrondi);
	  }
	  else  modifie=x;
	  if(expo(modifie)>=((lg(modifie)-2)<<TWOPOTBITS_IN_LONG))
	    ecrireall(x,'e',dec,chmp,flbl);
	  else
	  {
		/* partie entiere */
	    enti=gcvtoi(modifie,&e);
	    if(e>0) err(talker,"bug in ecrireall");
	    res=(long *)convi(enti);
	    d=coinit(*(--res));
	    while (*(--res)!= -1)
	    {
	      d=d+9;comilieu(*res);
	    }
	    pariputc('.');
		/* partie fractionnaire */
	    frac=subri(modifie,enti);
	    if(!signe (frac))
	    {
	      if (dec<0) dec=(long)(-expo(frac)*L2SL10+1);
	      dec=dec-d;
	      if (dec>0) zeros(dec);
	    }
	    else
	    {
	      if(!signe(enti))
	      {
		d=0;
		do
		{
		  p1=mulsr(1000000000,frac);
		  if((f=(expo(p1)<0)))
		  {
		    zeros(9);frac=p1;
		  }
		}
		while(f);
		do
		{
		  p1=mulsr(10,frac);
		  if((f=(expo(p1)<0)))
		  {
		    zeros(1);frac=p1;
		  }
		}
		while(f);
	      }
	      res=(long *)confrac(frac);
	      decmax= *res+++d;
	      if (dec<0) dec=decmax;
	      deceff=dec-decmax;
	      dec=dec-d;
	      while (dec>8)
	      {
		if (dec>deceff) comilieu(*res++);
		else zeros(9);
		dec=dec-9;
	      }
	      if (dec>0)
	      {
		if (dec>deceff) cofin(*res,dec);
		else zeros(dec);
	      }
	    }
	    avma=avmacourant;
	  }
	}
        break;
      case 'e':
	    /* impression d'un reel en format exponentiel */
        ex=expo(x);ex=(ex>=0)?(long)(ex*L2SL10) : (long)(-(-ex*L2SL10)-1);
        if (! sgn) 
	{
	  sprintf(thestring,flbl?" 0.E%ld":"0.E%ld",ex+1);
	  pariputs(thestring);
	}
        else
	{
	  avmacourant=avma;
	  affsr(10,dix=cgetr(lg(x)+1));
	  p1=(ex>0)?gdiv(x,gpuigs(dix,ex)):gmul(x,gpuigs(dix,-ex));
	  if(gcmpgs(gabs(p1,0),10)>=0)
	  {p1=gdivgs(p1,10);ex++;}
	  ecrireall(p1,'f',dec,chmp,flbl);
	  sprintf(thestring,flbl?" E%ld":"E%ld",ex);
	  pariputs(thestring);avma=avmacourant;
	}     
        break;
      case 'g':
	    /* impression d'un reel en format 'f',sauf s'il est trop petit */
        if(expo(x)>= -32) ecrireall(x,'f',dec,chmp,flbl);
        else ecrireall(x,'e',dec,chmp,flbl);
        break;
      default: err(formater);
    }
}

void
ecrire(GEN x, char format, long dec, long chmp)
{
  ecrireall(x,format,dec,chmp,1);
}

char *
gitoascii(GEN g, char* buf)
/*
  Convert g (which must be a PARI integer) to a (decimal) string in buf.
  Return buf.  The array buf must be large enough to hold g.
*/
{
  const long ltop = avma;

  GEN gleft = gabs(g,0);
  GEN gbillion = stoi(1000000000L);	/* Do 9 digits at a time */
  long i, irem, nchar = 0;

  if (typ(g) != 1) err(talker, "non-integer argument to gitoascii");

  while (1)
  {
    GEN gqr = gdiventres(gleft, gbillion);

    irem = itos((GEN)gqr[2]);	/* Remainder */
    gleft = (GEN)gqr[1];	/* Quotient */

    if (!signe(gleft)) break;

    for (i = 0; i <= 8; i++)
    {
      buf[nchar++] = (char)(irem % 10) + '0';
      irem /= 10;
    }
  } /* while */

  do
  {
    buf[nchar++] = (char)(irem % 10) + '0';
    irem /= 10;
  }
  while (irem != 0);

  if (signe(g) < 0)
  {
    buf[nchar++] = '-';			/* Check for negative number */
    nchar++;
  }
  buf[nchar] = '\0';			/* String terminator */

  for (i = 0; 2*i < nchar; i++)
  {
    char ch = buf[i];
    buf[i] = buf[nchar-1-i];
    buf[nchar-1-i] = ch;		/* Reverse order of digits */
  }

  avma = ltop;
  return buf;
} 

/********************************************************************/
/********************************************************************/
/**                                                                **/
/**                      SORTIE HEXADECIMALE                       **/
/**                                                                **/
/********************************************************************/
/********************************************************************/

void
sorstring(char* b, long x)
{
  char thestring[20];
#ifdef LONG_IS_64BIT
  sprintf(thestring,b,(uLong)x>>32,x&MAXHALFULONG);pariputs(thestring);
#else
  sprintf(thestring,b,x);pariputs(thestring);
#endif
}

static void
voir2(GEN x, long nb, long bl)
{
  long    tx=typ(x),i,j,e,dx,nb2,lx=lg(x);
  char    thestring[20];
#ifndef LONG_IS_64BIT
#define VOIR_STRING1 "[&=%08x] "
#define VOIR_STRING2 "%08x  "
#define VOIR_STRING3 "  %08x  :  "
#else
#define VOIR_STRING1 "[&=%08x%08x] "
#define VOIR_STRING2 "%08x%08x  "
#define VOIR_STRING3 "  %08x%08x  :  "
#endif
  
  bl+=2;
  sorstring(VOIR_STRING1,(uLong)x);
  if (nb<0)  nb2=lg(x);else nb2=nb;
  switch(tx)
  {
    case 1 : nb2=lgef(x);
    case 2 : for(i=0;i<nb2;i++) sorstring(VOIR_STRING2,x[i]);
      pariputc('\n');
      break;
    case 3 :
    case 9 : for(i=0;i<3;i++) sorstring(VOIR_STRING2,x[i]);
      pariputc('\n');
      blancs(bl);pariputs("mod = ");voir2((GEN)x[1],lgef((GEN)x[1]),bl);
      blancs(bl);
      if(tx==3) pariputs("int = ");
      else pariputs("pol = ");voir2((GEN)x[2],lgef((GEN)x[2]),bl);
      break;
    case 4 :
    case 5 :
    case 13:
    case 14: for(i=0;i<3;i++) sorstring(VOIR_STRING2,x[i]);
      pariputc('\n');
      blancs(bl);pariputs("num = ");voir2((GEN)x[1],lgef((GEN)x[1]),bl);
      blancs(bl);pariputs("den = ");voir2((GEN)x[2],lgef((GEN)x[2]),bl);
      break;
    case 6 : for(i=0;i<3;i++) sorstring(VOIR_STRING2,x[i]);
      pariputc('\n');
      blancs(bl);pariputs("real = ");voir2((GEN)x[1],nb,bl);
      blancs(bl);pariputs("imag = ");voir2((GEN)x[2],nb,bl);
      break;
    case 7 : for(i=0;i<5;i++) sorstring(VOIR_STRING2,x[i]);
      pariputc('\n');
      blancs(bl);pariputs("  p : ");voir2((GEN)x[2] ,lgef((GEN)x[2]),bl);
      blancs(bl);pariputs("p^l : ");voir2((GEN)x[3] ,lgef((GEN)x[3]),bl);
      blancs(bl);pariputs("  I : ");voir2((GEN)x[4] ,lgef((GEN)x[3]),bl);
      break;
    case 8 : for(i=0;i<4;i++) sorstring(VOIR_STRING2,x[i]);
      pariputc('\n');
      blancs(bl);pariputs("polynomial=");voir2((GEN)x[1],nb,bl);
      blancs(bl);pariputs("real = ");voir2((GEN)x[2],nb,bl);
      blancs(bl);pariputs("imag = ");voir2((GEN)x[3],nb,bl);
      break;
    case 10: for(i=0;i<lgef(x);i++) sorstring(VOIR_STRING2,x[i]);
      pariputc('\n');
      for(i=2;i<lgef(x);i++)
      {
	blancs(bl);
	sprintf(thestring, "coef of degree %ld = ",i-2);
	pariputs(thestring);voir2((GEN)x[i],nb, bl);
      }
      break;
    case 11: for(i=0;i<lx;i++) sorstring(VOIR_STRING2,x[i]);
      pariputc('\n');
      e=valp(x);
      if(signe(x))
        for(i=2;i<lx;i++)
	{
	  blancs(bl);
	  sprintf(thestring, "coef of degree %ld = ",e+i-2);
	  pariputs(thestring);voir2((GEN)x[i],nb, bl);
	}
      break;
    case 15:
    case 16:
    case 17:
    case 18: for(i=0;i<lx;i++) sorstring(VOIR_STRING2,x[i]);
      pariputc('\n');
      for(i=1;i<lx;i++)
      {
	blancs(bl);
	sprintf(thestring, "%ld-th component = ",i);
	pariputs(thestring);
	voir2((GEN)x[i],nb,bl);
      }
      break;
    case 19: for(i=0;i<lx;i++) sorstring(VOIR_STRING2,x[i]);
      pariputc('\n');
      if(lx>1)
      {
	dx=lg((GEN)x[1]);
	for (i=1;i<dx;i++)
	  for (j=1;j<lx;j++)
	  {
	    blancs(bl);
	    sprintf(thestring, "mat(%ld,%ld) = ",i,j);
	    pariputs(thestring) ;
	    voir2(gcoeff(x,i,j) ,nb, bl);
	  }
      }
  }
}

void
voir(GEN x, long nb)
{
  voir2(x,nb,0);
  pariputc('\n');
}

/********************************************************************/
/********************************************************************/
/**                                                                **/
/**                         SORTIE FORMATEE                        **/
/**                                                                **/
/********************************************************************/
/********************************************************************/

void
printvargp(long v)
{
  pariputs(varentries[v]->name);
}

static void
printvar(long v)
{
  (*printvariable)(v);
}

static void
sori(GEN g, char fo, long dd, long chmp)
{
  long  typy,sig,v,i,j,i0,e,l,l1,l2,n;
  long  dx,lx,av;
  char  thestring[50];
  GEN   p,a1,b1,a,b;
  
  typy=typ(g);if((typy==4)||(typy==5)) sig=gsigne(g);
  if ((typy>3)&&(typy<18))  chmp=0;
  if (gcmp0(g)&&(typy<17))
  {
    switch(typy)
    {
      case 2 : ecrire(g,fo,dd,chmp);break;
      case 3 :
      case 9 : pariputs("(0 mod ");sori((GEN)g[1],fo,-1,chmp);
	pariputc(')');break;
      case 7 :
	pariputs(" 0+O(");ecrire((GEN)g[2],fo,dd,chmp);
	sprintf(thestring, "^%ld",valp(g));
	pariputs(thestring);
	pariputc(')');break;
      case 11: pariputs(" 0+O(");printvar(ordvar[varn(g)]);
	sprintf(thestring, "^%ld)\n",valp(g)); pariputs(thestring);break;
      default: blancs(chmp-1);pariputc('0');
    }
	/*    if (typy>9) pariputc('\n');*/
  }
  else  if (gcmp1(g))
  {
    switch(typy)
    {
      case 2 : ecrire(g,fo,dd,chmp);break;
      case 3 :
      case 9 : pariputs("(1 mod ");sori((GEN)g[1],fo,-1,chmp);
	pariputc(')');break;
      case 7 : pariputs("1+O(");ecrire((GEN)g[2],fo,dd,chmp);
	sprintf(thestring, "^%ld",precp(g)); pariputs(thestring);pariputc(')');break;
      case 11: pariputs("1+O(");printvar(ordvar[varn(g)]);
	sprintf(thestring, "^%ld)\n",lg(g)-2); pariputs(thestring);break;
      default: blancs(chmp-1);pariputc('1');
    }
	/*    if (typy>9) pariputc('\n'); */
  }
  else
    if (((typy==4)||(typy==5))&&gcmp1((GEN)g[2])) ecrire((GEN)g[1],fo,dd,chmp);
    else
    {
      if ((typy>2)&&(typy<15))
      {
	if (((typy==4) || (typy==5))&&(sig<0)) pariputc('-');
	if ((typy!=13)&&(typy!=14)) pariputc('(');
      }
      switch(typy)
      {
	case 1 :
	case 2 : ecrire(g,fo,dd,chmp);break;
	case 3 :
	  if (signe((GEN)g[2])<0)
	  {
	    l=avma;sori(addii((GEN)g[2],(GEN)g[1]),fo,-1,chmp);
	    avma=l;
	  }
	  else sori((GEN)g[2],fo,dd,chmp);
	  pariputs(" mod ");
	  sori((GEN)g[1],fo,dd,chmp);
	  break;
	    
	case 9 :
	  sori((GEN)g[2],fo,dd,chmp);
	  pariputs(" mod ");
	  sori((GEN)g[1],fo,dd,chmp);
	  break;
	    
	case 4 :
	case 5 :
	  a=(GEN)g[1];
	  if (sig<0)
	  {setsigne(a,1);ecrire(a,fo,dd,chmp);setsigne(a,-1);}
	  else ecrire(a,fo ,dd,chmp);
	  if (!gcmp1((GEN)g[2]))
	  {pariputs(" /");ecrire((GEN)g[2],fo,dd,chmp);}
	  break;
	    
	case 6 :
	  a=(GEN)g[1];b=(GEN)g[2];
	  if (!gcmp0(a)) sori(a,fo,dd,chmp);
	  if((signe(b)>0)&&!gcmp0(a)) pariputc('+');
	  else pariputc(' ');
	  if (!gcmp0(b))
	  {
	    if (gcmp1(b)) pariputs(" i");
	    else
	    {
	      if(gcmp_1(b)) pariputs("-i");
	      else {sori(b,fo,dd,chmp);pariputs(" i");}
	    }
	  }
	  break;
	    
	case 7 :                          /*  ecrire un p-adique  */
	  e=valp(g);l=precp(g);
	  av=avma;
	  a1=gcopy((GEN)g[4]);p=(GEN)g[2];
	  for (i=0;i<l;i++)
	  {
	    a1=dvmdii(a1,p,&b1);
	    if (signe(b1))
	    {
	      if (!(e+i) || (!gcmp1(b1)))
	      {
		ecrire(b1,fo,dd,chmp);
		if((e+i)) pariputc('*');else pariputc(' ');
	      }
	      if (e+i==1) {ecrire(p,fo,dd,chmp);pariputc(' ');}
	      else if (e+i) {ecrire(p,fo,dd,chmp);sprintf(thestring, "^%ld ",e+i); pariputs(thestring);}
	      pariputc('+');
	    }
	  }
	  pariputs(" O(");
	  if (!(e+l)) pariputs(" 1");
	  else {ecrire(p,fo,dd,chmp);if((e+l)!=1) sprintf(thestring, "^%ld",e+l); pariputs(thestring);}
	  pariputc(')');
	  avma=av;
	  break;
	    
	case  8 :
	  a=(GEN)g[2];b=(GEN)g[3];
	  if (!gcmp0(a)) sori(a,fo,dd,chmp);
	  if((signe(b)>0)&&!gcmp0(a)) pariputs(" +");
	  else pariputc(' ');
	  if (!gcmp0(b))
	  {
	    if (gcmp1(b)) pariputs(" w");
	    else
	    {
	      if(gcmp_1(b)) pariputs("-w");
	      else {sori(b,fo,dd,chmp);pariputs(" w");}
	    }
	  }
	  break;
	    
	case 10 :                                 /* sortir un polynome */
	  i0=gval(g,varn(g))+2;l=lgef(g)-1;v=ordvar[varn(g)];
	  for (i=l;i>=i0;i--)
	  {
	    a=(GEN)g[i];
	    if (!gcmp0(a))
	    {
	      if ((i==l)&&gcmp_1(a) &&
		  (l>2)&&(typ(a)!=3)&&(typ(a)!=9)) pariputc('-');
	      if ((!gcmp1(a)&&!gcmp_1(a)) || (i==2) || (typ(a)==3)
		  || (typ(a)==9)) sori(a,fo,dd,chmp);
	      if (i==3) {pariputc(' ');printvar(v);pariputc(' ');}
	      if (i>3) {pariputc(' ');printvar(v);sprintf(thestring, "^%ld ",i-2); pariputs(thestring);}
	    }
	    if (i>i0)
	    {
	      b=(GEN)g[i-1];if(!gcmp0(b))
	      {
		if ((i>3)&&gcmp_1(b)&&(typ(b)!=3)&&(typ(b)!=9)) 
		  pariputc('-');
		else if (((signe(b)>0)||(typ(b)==3)||(typ(b)>5)))
		  pariputc('+');
	      }
	    }
	  }
	  break;
	    
	case 11 :                                 /* serie    */
	  e=valp(g)-2;l=lg(g);v=ordvar[varn(g)];
	  for (i=2;i<l;i++)
	  {
	    a=(GEN)g[i];
	    if (!gcmp0(a))
	    {
	      if (!(e+i) || (!gcmp1(a)&&!gcmp_1(a)) || (typ(a)==3)
		  || (typ(a)==9))
	      {
		sori(a,fo,dd,chmp);
		if(!(e+i)) pariputc(' ');
	      }
	      else if (gcmp_1(a)) pariputc('-');
	      if (e+i==1) {pariputc(' ');printvar(v);pariputc(' ');}
	      if (e+i>1) {pariputc(' ');printvar(v);sprintf(thestring, "^%ld ",e+i); pariputs(thestring);}
	      if (e+i<0) {pariputc(' ');printvar(v);sprintf(thestring, "^(%ld); pariputs(thestring) ",e+i);}
	    }
	    b=(GEN)g[i+1];
	    if ((i<l-1)&&((typ(b)==3) || (typ(b)>5) ||(signe(b)>0)))
	      pariputc('+');
	  }
	  if (!(e+l)) pariputs("+ O(1)");
	  else if (e+l==1){pariputs("+ O(");printvar(v);pariputc(')');}
	  else {pariputs("+ O(");printvar(v);sprintf(thestring, "^%ld)",e+l); pariputs(thestring);}
	  break;
	    
	case 13 :
	    
	case 14 :  pariputs("\n\n");
	  l1=lg((GEN)g[1]);l2=lg((GEN)g[2]);
	  l=(l1>l2) ? l1-2 : l2-2;
	  sori((GEN)g[1],fo,dd,chmp);pariputc('\n');
	  for (n=1;n<l;n++)
	    pariputs("----------");pariputc('\n');
	  sori((GEN)g[2],fo,dd,chmp);
	  break;
	    
	case 15: pariputc('{');sori((GEN)g[1],fo,dd,chmp);pariputc(','); 
	  sori((GEN)g[2],fo,dd,chmp);pariputc(',');sori((GEN)g[3],fo,dd,chmp); 
	  pariputc(',');sori((GEN)g[4],fo,dd,chmp);pariputs("}\n"); 
	  break; 
	    
	case 16: pariputc('{');sori((GEN)g[1],fo,dd,chmp);pariputc(',');
	  sori((GEN)g[2],fo,dd,chmp);pariputc(',');sori((GEN)g[3],fo,dd,chmp);
	  pariputs("}\n");
	  break;
	    
	case 17 :                                    /* vecteur ligne  */
	  pariputc('[');
	  for (i=1;i<lg(g);i++)
	  {
	    sori((GEN)g[i],fo,dd,chmp);
	    if (i<lg(g)-1) pariputc(',');
	  }
	  pariputs("]\n");
	  break;
	    
	case 18 :                                   /* vecteur colonne */
	  if(lg(g)==1) pariputs("[]\n");
	  else
	    for (i=1;i<lg(g);i++)
	    {
	      pariputc('[');
	      sori((GEN)g[i],fo,dd,chmp);
	      pariputs("]\n");
	    }
	  break;
	    
	case 19 :
	  pariputc('\n');lx=lg(g);dx=(lx>1)?lg((GEN)g[1]):2;
	  for (i=1;i<dx;i++)
	  {
	    pariputc('[');
	    for (j=1;j<lx;j++)
	    {
	      sori(gcoeff(g,i,j),fo,dd,chmp);
	      if(j<lx-1) pariputc(' ');
	    }
	    if(i<dx-1) pariputs("]\n\n");else pariputs("]\n");
	  }
	  break;
	default: sorstring(VOIR_STRING2,*g);
      }
      if ((typy>2)&&(typy<13)) pariputc(')');
    }                             /* fin du else */
}

void
sor(GEN g, char fo, long dd, long chmp)
{
  long av=avma;
  if(varchanged) sori(changevar(g,polvar), fo, dd, chmp);
  else sori(g, fo, dd, chmp);
  avma = av;
}

void
etatpile(unsigned int n)
{
  long  nu,i,l,m;
  GEN adr,adr1;
  double r;
  char thestring[80];
  
  nu=(top-avma)/BYTES_IN_LONG;
  l=(top-bot)/BYTES_IN_LONG;
  r=100.0*nu/l;
  sprintf(thestring, "\n Top : %lx   Bottom : %lx   Current stack : %lx\n",top,bot,avma); pariputs(thestring);
  sprintf(thestring, " Used :                         %ld  long words  (%ld K)\n",nu,nu/1024*BYTES_IN_LONG); pariputs(thestring);
  sprintf(thestring, " Available :                    %ld  long words  (%ld K)\n",(l-nu),(l-nu)/1024*BYTES_IN_LONG); pariputs(thestring);
  sprintf(thestring, " Occupation of the PARI stack : %6.2f percent\n",r); pariputs(thestring);
  getheapaux(&m, &l);
  sprintf(thestring, " %ld objects on heap occupy %ld long words\n\n", m, l); pariputs(thestring);
  
  if (n)
  {
    if (n>nu) n=nu;
    adr=(GEN)avma;adr1=adr+n;
    while (adr<adr1)
    {
      sorstring(VOIR_STRING3,(uLong)adr);
      l=lg(adr);m=(adr==polvar) ? MAXVAR+1 : 0;
      for (i=0;(i<l)&&(adr<adr1);i++,adr++) sorstring(VOIR_STRING2,*adr);
      pariputc('\n');if(m) adr=polvar+m;
    }
    pariputc('\n');
  }
}

/********************************************************************/
/********************************************************************/
/**                                                                **/
/**                          SORTIE BRUTE                          **/
/**                                                                **/
/********************************************************************/
/********************************************************************/

static long
isnull(GEN g)
{
  long i;
  switch (typ(g))
  {
    case 1: return !signe(g);
    case 6: return isnull((GEN)g[1])&&isnull((GEN)g[2]);
    case 8: return isnull((GEN)g[2])&&isnull((GEN)g[3]);
    case 4:
    case 5:
    case 13:
    case 14: return isnull((GEN)g[1]);
    case 10: for (i=lgef(g)-1;i>1;i--) if (!isnull((GEN)g[i])) return 0;
      return 1;
    case 9: return isnull((GEN)g[2]);
    default: return 0;
  }
}

static long
isone(GEN g) /* renvoie 1 ou-1 si g est 1 ou-1,0 sinon */
{
  long i,sig;
  switch (typ(g))
  {
    case 1: if(!signe(g)) return 0;
    else return (g[2]==1)&&(lgef(g)==3) ? signe(g) : 0;
    case 6: return isnull((GEN)g[2]) * isone((GEN)g[1]);
    case 8: return isnull((GEN)g[3]) * isone((GEN)g[2]);
    case 4:
    case 5:
    case 13:
    case 14: return isone((GEN)g[1])*isone((GEN)g[2]);
    case 10: if(!signe(g)) return 0;
      if (!(sig=isone((GEN)g[2]))) return 0;
    for (i=lgef(g)-1;i>2;i--) if (!isnull((GEN)g[i])) return 0;
    return sig;
    default: return 0;
  }
}

static long
isfactor(GEN g) /* si g est un monome,renvoie son signe,0 sinon */
{
  long i,deja=0,sig=1;
  switch(typ(g))
  {
    case 1:
    case 2: return signe(g)<0 ?-1 : 1;
    case 4:
    case 5:
    case 13:
    case 14: return isfactor((GEN)g[1]);
    case 6: if (isnull((GEN)g[1])) return isfactor((GEN)g[2]);
      return isnull((GEN)g[2]) ? isfactor((GEN)g[1]) : 0;
    case 7: return !signe((GEN)g[4]);
    case 8: if (isnull((GEN)g[2])) return isfactor((GEN)g[3]);
      return isnull((GEN)g[3]) ? isfactor((GEN)g[2]) : 0;
    case 10:
      for (i=lgef(g)-1;i>1;i--)
        if (!isnull((GEN)g[i]))
	{
	  if (deja) return 0;
	  sig=isfactor((GEN)g[i]);
	  deja=1;
	}
      return sig ? sig : 1;
    case 11: if(!signe(g)) return 1;
      for (i=lg(g)-1;i>1;i--) if (!isnull((GEN)g[i])) return 0;
    default: return 1;
  }
}

static long
isdenom(GEN g) /* renvoie 1 si g est un truc... */
{
  long i,deja=0;
  switch(typ(g))
  {
    case 4:
    case 5:
    case 13:
    case 14: return 0;
    case 6: return isnull((GEN)g[2]);
    case 7: return !signe((GEN)g[4]);
    case 8: return isnull((GEN)g[3]);
    case 10:
      for (i=lgef(g)-1;i>1;i--)
        if (!isnull((GEN)g[i]))
	{
	  if (deja) return 0;
	  if (i==2) return isdenom((GEN)g[2]);
	  if (!isone((GEN)g[i])) return 0;
	  deja=1;
	}
      return 1;
    case 11: if(!signe(g)) return 1;
      for (i=lg(g)-1;i>1;i--) if (!isnull((GEN)g[i])) return 0;
    default: return 1;
  }
}

#define putsigne(x) pariputs(x>0?" + ":" - ")
#define putsignb(x,f) pariputs(x>0?(f?" + ":"+"):(f?" - ":"-"))

static void
monome(long v, long deg)
{
  char thestring[20];
  if (deg)
  {
    printvar(v);
    if (deg!=1) {sprintf(thestring, "^%ld",deg); pariputs(thestring);}
  }
  else pariputc('1');
}

static void
brutiall(GEN g, char format, long dec, long sanssigne, long flbl)
{
  long e,l,sig,i,j,r,v,av=avma;
  GEN a1,b1,p;
  char thestring[20];
  
  if (isnull(g)) pariputc('0');
  else if((sig=isone(g))) {if (!sanssigne&&(sig<0)) pariputc('-');pariputc('1');}
  else switch(typ(g))
  {
    case 1:
    case 2: if (sanssigne&&(signe(g)<0)) g=mpabs(g);
      ecrireall(g,format,dec,0,flbl);break;
    case 3:
    case 9: pariputs("mod(");brutiall((GEN)g[2],format,dec,0,flbl);
      pariputs(flbl?", ":",");
      brutiall((GEN)g[1],format,dec,0,flbl);pariputc(')');break;
    case 4:
    case 5:
    case 13:
    case 14:
      if (!(sig=isfactor((GEN)g[1]))) pariputc('(');
      brutiall((GEN)g[1],format,dec,sanssigne,flbl);
      if (!sig) pariputc(')');
      pariputc('/');
      if (!(sig=isdenom((GEN)g[2]))) pariputc('(');
      brutiall((GEN)g[2],format,dec,0,flbl);
      if (!sig) pariputc(')');
      break;
    case 6:
    case 8:
      r=(typ(g)==8);
      if (isnull((GEN)g[r+1])) 
        if ((sig=isone((GEN)g[r+2])))
	{
	  if (!sanssigne&&(sig<0)) pariputc('-');
	  pariputc(r ? 'w' : 'i');
	}
        else
	{
	  if (!(sig=isfactor((GEN)g[r+2]))) pariputc('(');
	  brutiall((GEN)g[r+2],format,dec,sanssigne,flbl);
	  if (!sig) pariputc(')');
	  pariputc('*');
	  pariputc(r ? 'w' : 'i');
	}
      else
      {
	brutiall((GEN)g[r+1],format,dec,sanssigne,flbl);
	if (!isnull((GEN)g[r+2]))
	  if ((sig=isone((GEN)g[r+2])))
	  {
	    putsignb(sig,flbl);pariputc(r ? 'w' : 'i');
	  }
	  else
	  {
	    if ((sig=isfactor((GEN)g[r+2]))) putsignb(sig,flbl);
	    else pariputs(flbl?" + (":"+(");
	    brutiall((GEN)g[r+2],format,dec,1,flbl);
	    if (!sig) pariputc(')');
	    pariputc('*');
	    pariputc(r ? 'w' : 'i');
	  }
      }
      break;
    case 10:
      v=ordvar[varn(g)];for (i=lgef(g)-1;isnull((GEN)g[i]);i--);
      if ((sig=isone((GEN)g[i]))) 
      {
	if (!sanssigne&&(sig<0)) pariputc('-');
	monome(v,i-2);
      }
      else
      {
	if (isfactor((GEN)g[i])) 
	  brutiall((GEN)g[i],format,dec,sanssigne,flbl);
	else
	{
	  pariputc('(');
	  brutiall((GEN)g[i],format,dec,0,flbl);
	  pariputc(')');
	}
	if (i>2) {pariputc('*');monome(v,i-2);}
      }
      for(;--i>1;) if (!isnull((GEN)g[i]))
        if ((sig=isone((GEN)g[i]))) {putsignb(sig,flbl);monome(v,i-2);}
        else
	{
	  if ((sig=isfactor((GEN)g[i]))) putsignb(sig,flbl);
	  else pariputs(flbl?" + (":"+(");
	  brutiall((GEN)g[i],format,dec,sig,flbl);
	  if (!sig) pariputc(')');
	  if (i>2) {pariputc('*');monome(v,i-2);}
	}
      break;
    case 7:
      e=valp(g);l=precp(g);
      a1=(GEN)g[4];p=(GEN)g[2];
      for (i=0;i<l;i++)
      {
	a1=dvmdii(a1,p,&b1);
	if (signe(b1))
	{
	  if (!(e+i)||!gcmp1(b1))
	  {
	    ecrireall(b1,format,0,0,flbl);
	    if (e+i) pariputc('*');
	  }
	  if(e+i)
	  {
	    ecrireall(p,format,0,0,flbl);
	    if ((e+i)!=1) 
	    {
	      sprintf(thestring,flbl?"^%d ":"^%d",e+i); 
	      pariputs(thestring);
	    }
	  }
	  pariputs(flbl?" + ":"+");
	}
      }
      pariputs("O(");
      ecrireall(p,format,0,0,flbl);
      if ((e+l)!=1) {sprintf(thestring,"^%ld",e+l); pariputs(thestring);}
      pariputc(')');
      break;
    case 11:
      e=valp(g)-2;v=ordvar[varn(g)];
      if (signe(g))
      {
	l=lg(g);
	if ((sig=isone((GEN)g[2]))) {if (sig<0) pariputc('-');monome(v,2+e);}
	else
	{
	  if (!(sig=isfactor((GEN)g[2]))) pariputc('(');
	  brutiall((GEN)g[2],format,dec,sanssigne,flbl);
	  if (!sig) pariputc(')');
	  if (valp(g)) {pariputc('*');monome(v,valp(g));}
	}
	for(i=3;i<l;i++) if (!isnull((GEN)g[i]))
	  if ((sig=isone((GEN)g[i]))) {putsignb(sig,flbl);monome(v,i+e);}
	  else
	  {
	    if ((sig=isfactor((GEN)g[i]))) putsignb(sig,flbl);
	    else pariputs(flbl?" + (":"+(");
	    brutiall((GEN)g[i],format,dec,sig,flbl);
	    if (!sig) pariputc(')');
	    if ((i+e)!=0) {pariputc('*');monome(v,i+e);}
	  }
	pariputs(flbl?" + ":"+");
      }
      else  l=2;
      pariputs("O(");
      printvar(v);
      if ((e+l)!=1) {sprintf(thestring,"^%ld",e+l);pariputs(thestring);}
      pariputc(')');
      break;
    case 15: pariputs("qfr(");brutiall((GEN)g[1],format,dec,0,flbl);
      pariputs(flbl?", ":",");brutiall((GEN)g[2],format,dec,0,flbl);
      pariputs(flbl?", ":",");brutiall((GEN)g[3],format,dec,0,flbl);
      pariputs(flbl?", ":",");brutiall((GEN)g[4],format,dec,0,flbl);
      pariputc(')');
      break;
    case 16: pariputs("qfi(");brutiall((GEN)g[1],format,dec,0,flbl);
      pariputs(flbl?", ":",");brutiall((GEN)g[2],format,dec,0,flbl);
      pariputs(flbl?", ":",");brutiall((GEN)g[3],format,dec,0,flbl); 
      pariputc(')'); 
      break;
    case 17:
    case 18:
      pariputc('[');
      for(i=1;i<lg(g);i++)
      {
	brutiall((GEN)g[i],format,dec,0,flbl);
	if (i<lg(g)-1) pariputs(flbl?", ":",");
      }
      pariputc(']');
      if (typ(g)==18) pariputc('~');
      break;
    case 19:
      pariputc('[');
      if (lg(g)>1) for(i=1;i<lg((GEN)g[1]);i++)
      {
	for(j=1;j<lg(g);j++)
	{
	  brutiall((GEN)((GEN)g[j])[i],format,dec,0,flbl);
	  if (j<lg(g)-1) pariputs(flbl?", ":",");
	}
	if (i<lg((GEN)g[1])-1) pariputs(flbl?"; ":";");
      }
      pariputc(']');break;
    default: sorstring(VOIR_STRING2,*g);
  }
  avma=av;
}

void
bruteall(GEN g, char format, long dec, long flbl)
{
  long av=avma;
  if(varchanged) brutiall(changevar(g,polvar),format,dec,0,flbl);
  else brutiall(g,format,dec,0,flbl);
  avma=av;
}

void
brute(GEN g, char format, long dec)
{
  bruteall(g,format,dec,1);
}

static void
matbruti(GEN g, char format, long dec, long sanssigne)

/* for this function the only difference with bruti is in type 19 */

{
  long i,j,av=avma,lx,dx;
  
  if(typ(g)<19) brutiall(g,format,dec,sanssigne,1);
  else 
  {
    pariputc('\n');lx=lg(g);dx=(lx>1)?lg((GEN)g[1]):2;
    for (i=1;i<dx;i++)
    {
      pariputc('[');
      for (j=1;j<lx;j++)
      {
	brutiall((GEN)((GEN)g[j])[i],format,dec,0,1);
	if(j<lx-1) pariputc(' ');
      }
      if(i<dx-1) pariputs("]\n\n");else pariputs("]\n");
    }
  }
  avma=av;
}

void
matbrute(GEN g, char format, long dec)
{
  long av=avma;
  if(varchanged) matbruti(changevar(g,polvar),format,dec,0);
  else matbruti(g,format,dec,0);
  avma=av;
}


/********************************************************************/
/********************************************************************/
/**                                                                **/
/**                          FORMATTAGE TeX                        **/
/**                                                                **/
/********************************************************************/
/********************************************************************/

static void
texnome(long v, long deg)
{
  char    thestring[20];
  if (deg)
  {
    printvar(v);
    if (deg!=1) {sprintf(thestring, "^{%ld}",deg); pariputs(thestring);}
  }
  else pariputc('1');
}

static void
texi(GEN g, char format, long dec, long sanssigne)
{
  long e,l,sig,i,j,r,v,av=avma;
  GEN a1,b1,p;
  char thestring[20];
  
  pariputc('{');
  if (isnull(g)) pariputc('0');
  else if ((sig=isone(g))) {if (!sanssigne&&(sig<0)) pariputc('-');pariputc('1');}
  else switch(typ(g))
  {
    case 1:
    case 2:
      if (sanssigne&&(signe(g)<0)) g=mpabs(g);
      ecrire(g,format,dec,0);break;
    case 3:
    case 9:
      texi((GEN)g[2],format,dec,0);pariputs("mod");
      texi((GEN)g[1],format,dec,0);break;
    case 4:
    case 5:
    case 13:
    case 14:
      texi((GEN)g[1],format,dec,sanssigne);
      pariputs("\\over");
      texi((GEN)g[2],format,dec,0);
      break;
    case 6:
    case 8:
      r=(typ(g)==8);
      if (isnull((GEN)g[r+1])) 
        if ((sig=isone((GEN)g[r+2]))) {if (!sanssigne&&(sig<0)) pariputc('-');pariputc(r ? 'w' : 'i');}
        else
	{ 
	  if (!(sig=isfactor((GEN)g[r+2]))) pariputc('('); 
	  texi((GEN)g[r+2],format,dec,sanssigne); 
	  if (!sig) pariputc(')'); 
	  pariputc(r ? 'w' : 'i'); 
	}
      else
      { 
	texi((GEN)g[r+1],format,dec,sanssigne); 
	if (!isnull((GEN)g[r+2])) 
	  if ((sig=isone((GEN)g[r+2]))) {putsigne(sig);pariputc(r ? 'w' : 'i');} 
	  else 
	  { 
	    if ((sig=isfactor((GEN)g[r+2]))) putsigne(sig); 
	    else pariputs("+("); 
	    texi((GEN)g[r+2],format,dec,1); 
	    if (!sig) pariputc(')'); 
	    pariputc(r ? 'w' : 'i'); 
	  } 
      }
      break;
    case 10:
      v=ordvar[varn(g)];for (i=lgef(g)-1;isnull((GEN)g[i]);i--);
      if ((sig=isone((GEN)g[i]))) {if (!sanssigne&&(sig<0)) pariputc('-');texnome(v,i-2);}
      else
      { 
/*	  if (!(sig=isfactor((GEN)g[i]))) pariputc('('); 
	  texi((GEN)g[i],format,dec,sanssigne); 
	  if (!sig) pariputc(')'); 
	  if (i>2) texnome(v,i-2);  */
	if (isfactor((GEN)g[i])) texi((GEN)g[i],format,dec,sanssigne); 
	else
	{
	  pariputc('('); 
	  texi((GEN)g[i],format,dec,0); 
	  pariputc(')');
	}
	if (i>2) texnome(v,i-2); 
      }
      for(;--i>1;) if (!isnull((GEN)g[i]))
        if ((sig=isone((GEN)g[i]))) {putsigne(sig);texnome(v,i-2);}
        else
	{ 
	  if ((sig=isfactor((GEN)g[i]))) putsigne(sig);else pariputs("+("); 
	  texi((GEN)g[i],format,dec,sig); 
	  if (!sig) pariputc(')'); 
	  if (i>2) texnome(v,i-2); 
	}
      break;
    case 7:
      e=valp(g);l=precp(g);
      a1=(GEN)g[4];p=(GEN)g[2];
      for (i=0;i<l;i++)
      { 
	a1=dvmdii(a1,p,&b1); 
	if (signe(b1)) 
	{ 
	  if (!(e+i) || !gcmp1(b1)) 
	  { 
	    ecrire(b1,format,0,0); 
	    if (e+i) pariputs("\\cdot"); 
	  } 
	  if (e+i) 
	  { 
	    ecrire(p,format,0,0); 
	    if ((e+i)!=1) {sprintf(thestring, "^{%ld}",e+i); pariputs(thestring);} 
	  } 
	  pariputc('+'); 
	} 
      }
      pariputs("O(");
      ecrire(p,format,0,0);if ((e+l)!=1) {sprintf(thestring, "^{%ld}",e+l); pariputs(thestring);}
      pariputc(')');
      break;
    case 11:
      e=valp(g)-2;v=ordvar[varn(g)];
      if (signe(g))
      { 
	l=lg(g); 
	if ((sig=isone((GEN)g[2]))) {if (sig<0) pariputc('-');texnome(v,2+e);} 
	else 
	{ 
	  if (!(sig=isfactor((GEN)g[2]))) pariputc('('); 
	  texi((GEN)g[2],format,dec,sanssigne); 
	  if (!sig) pariputc(')'); 
	  if (valp(g)) texnome(v,valp(g)); 
	} 
	for(i=3;i<l;i++) if (!isnull((GEN)g[i])) 
	  if ((sig=isone((GEN)g[i]))) {putsigne(sig);texnome(v,i+e);} 
	  else 
	  { 
	    if ((sig=isfactor((GEN)g[i]))) putsigne(sig);else pariputs(" + ("); 
	    texi((GEN)g[i],format,dec,sig); 
	    if (!sig) pariputc(')'); 
	    if (i+e) texnome(v,i+e); 
	  } 
	pariputc('+'); 
      }
      else  l=2;
      pariputs("O(");
      printvar(v);if ((e+l)!=1) {sprintf(thestring, "^{%ld}",e+l); pariputs(thestring);}
      pariputc(')');
      break;
    case 15: pariputs("qfr(");texi((GEN)g[1],format,dec,0);pariputs(", ");
      texi((GEN)g[2],format,dec,0);pariputs(", ");texi((GEN)g[3],format,dec,0);
      pariputs(", ");texi((GEN)g[4],format,dec,0);
      pariputc(')');
      break;
    case 16: pariputs("qfi(");texi((GEN)g[1],format,dec,0);pariputs(", "); 
      texi((GEN)g[2],format,dec,0);pariputs(", ");texi((GEN)g[3],format,dec,0); 
      pariputc(')'); 
      break; 
      
    case 17:
      pariputs("\\pmatrix{");
      for(i=1;i<lg(g);i++)
      {
	texi((GEN)g[i],format,dec,0);
	if (i<lg(g)-1) pariputc('&');
      }
      pariputs("\\cr}");
      break;
    case 18:
      pariputs("\\pmatrix{");
      for(i=1;i<lg(g);i++)
      {
	texi((GEN)g[i],format,dec,0);
	pariputs("\\cr");
      }
      pariputc('}');
      break;
    case 19:
      pariputs("\\pmatrix{");
      if (lg(g)>1) for(i=1;i<lg((GEN)g[1]);i++)
      {
	for(j=1;j<lg(g);j++)
	{
	  texi((GEN)((GEN)g[j])[i],format,dec,0);
	  if (j<lg(g)-1) pariputc('&');
	}
	pariputs("\\cr");
      }
      pariputc('}');
  }
  avma=av;
  pariputc('}');
}

void
texe(GEN g, char format, long dec)
{
  long av=avma;
  if(varchanged) texi(changevar(g,polvar),format,dec,0);
  else texi(g,format,dec,0);
  avma=av;
}

/********************************************************************/
/********************************************************************/
/**                                                                **/
/**                GESTION DES FICHIERS IN, OUT ET LOG             **/
/**                                                                **/
/********************************************************************/
/********************************************************************/

void
switchin(char *name)
{
#ifdef __cplusplus
  const int MAXFILES = 5;
#else
# define MAXFILES 5
#endif
  static FILE *stack[MAXFILES];
  static long depth = 0;
  if (name)
  {
    char namecopy[100];
    if (depth >= MAXFILES) err(includer1);
    stack[depth++] = infile;
    if((infile = fopen(name, "r"))) return;
    name=strcpy(namecopy, name);
    name=strcat(name,".gp");
    if((infile = fopen(name, "r"))) return;
    infile = stack[--depth];
    err(inputer1);
  }
  if(!depth) quitting_pari=1;
  else {fclose(infile); infile = stack[--depth];}
}

void
switchout(char *name)
{
  if (name)
  {
    FILE *glou = fopen(name, "a");
    if (!glou) err(outputer1);
    outfile = glou;
  }
  else
  {
    fclose(outfile);
    outfile = stdout;
  }
}

void
fliplog(void)
{
  if (logfile)
  {
    fclose(logfile);
    logfile = NULL;
    pariputs("    logging off\n");
  }
  else
  {
    logfile = fopen("pari.log", "w");
    if (!logfile) err(outloger1);
    pariputs("    logging on\n");
  }
}

void
imprimer(GEN x)
{
  bruteall(x,'g',-1,1);printf("\n");fflush(stdout);
}

void
outerr(GEN x)
{
  FILE *tmpfile;
  
  tmpfile=outfile;outfile=errfile;output(x);outfile=tmpfile;
}

void
fprintferr(char* pat, ...)
{
    va_list args;

    va_start(args, pat);

    if (pariErr && pariErr->puts) { 
      char buff[1024];
      vsprintf(buff,pat,args); 
      (pariErr->puts)(buff);
    } else vfprintf(errfile,pat,args);

    va_end(args);
}

void
flusherr()
{
  if (pariErr && pariErr->flush) {
    (pariErr->flush)();
  } else {
    fflush(errfile);
  }
}

void
bruterr(GEN x,char format,long dec)
{
  FILE *tmpfile;

  tmpfile=outfile;outfile=errfile;bruteall(x,format,dec,1);outfile=tmpfile;
}

void
outbeauterr(GEN x)
{
  FILE *tmpfile;

  tmpfile=outfile;outfile=errfile;outbeaut(x);outfile=tmpfile;
}

void
printversion(void)
{
  pariputs("            GP/PARI CALCULATOR Version 1.39.03\n");
  pariputs(pariversion);
}

void
addhelp(entree *ep, char *s)
{
    if (ep->help && !(ep->valence & EpSTATIC)) {
	free(ep->help);
	ep->valence &= ~EpSTATIC;
    }
    ep->help = (char*)malloc(strlen(s) + 1);
    strcpy(ep->help, s);
}

void
freeep(entree *ep)
{
    if (foreignFuncFree && ep->code && (*ep->code == 'x')) {
	/* The function is created by foreign interpreter, ask it
	 * to free the function. */
	(*foreignFuncFree)(ep);
    }
    if (ep->help && !(ep->valence & EpSTATIC)) free(ep->help);
    free(ep); 
}

GEN
killep(entree *ep)
{
    GEN res = gnil;
    int i;
    entree *ep1;

    if (EpVALENCE(ep)<100) err(killer1);
    killvalue(ep);
    if (EpVALENCE(ep) == 200) res = (GEN)ep->value;
    else {
	for(i = 0; i<TBLSZ; i++) {
	    if (hashtable[i] == ep) {
		hashtable[i] = ep->next;
		freeep(ep); 
		break;
	    } else {
		for(ep1 = hashtable[i]; ep1; ep1 = ep1->next) {
		    if (ep1->next == ep) {
			ep1->next = ep->next; 
			freeep(ep); 
			break;
		    }
		}
	    }
	}
    }
    return res;
}
