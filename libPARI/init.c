/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                                                                 */
/*              PROGRAMME D'INITIALISATION DU SYSTEME              */
/*                                                                 */
/*                    ET TRAITEMENT DES ERREURS                    */
/*                                                                 */
/*                       copyright Babe Cool                       */
/*                                                                 */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include        "genpari.h"
#include "rect.h"

const int STACKSIZE = 5000;  /* nombre de gn possibles */
const int TBLSZ = 135;  /* taille de la table de hashcodes */
const int NUMPRTBELT = 20; /* taille table de premiers prives */

const double K = 9.632959862*(BYTES_IN_LONG/4);  /* 32*log(2)/log(10)  */
const double K1 = 0.103810253/(BYTES_IN_LONG/4); /* log(10)/(32*log(2))*/
const double K2 = 1.1239968;               /* 1/(1-(log(2)/(2*pi)))    */
const double K4 = 17.079468445347/BITS_IN_LONG;  /* 2*e*pi/32          */
const double LOG2 = 0.69314718055994531;     /* log(2)                   */
const double L2SL10 = 0.301029995663981;   /* log(2)/log(10)           */
#ifndef  PI
const double PI = 3.141592653589;          /* pi                       */
#endif
const double rac5 = 2.23606797749;         /* racine de 5              */
const double C1 = 0.9189385332;            /* log(2*pi)/2              */
const double C2 = 22.18070978*(BYTES_IN_LONG/4);  /* 32*log(2)         */
const double C3 = 0.0216950598/(BYTES_IN_LONG/4); /* log((1+sqrt(5))/2)/(32*log(2)) */
#ifdef LONG_IS_64BIT
const double C31 = 9223372036854775808.0;  /* 2^63 */
#else
const double C31 = 2147483648.0;           /* 2^31                     */
#endif

#ifndef LONG_IS_64BIT
const long BIGINT = 32767;                 /* 2^15-1                   */
const long EXP220 = 1048576;               /* 2^20                     */
const long VERYBIGINT = 2147483647;        /* 2^31-1                   */

#else

const long BIGINT = 2147483647;                 /* 2^31-1              */
const long EXP220 = 1099511627776;              /* 2^40                */
const long VERYBIGINT = 9223372036854775807;    /* 2^63-1              */
#endif

/*      Variables statiques communes :          */

unsigned long top,bot,avma;
long    avloc;
#ifndef LONG_IS_64BIT
long    prec=5;
#else
long    prec=4;
#endif
long    precdl=16, defaultpadicprecision=16;
long    tglobal,paribuffsize=30000,pariecho=0;
long    compact_arrays;
long    quitting_pari=0;
jmp_buf environnement;
unsigned long init_opts = INIT_JMPm | INIT_SIGm;
FILE    *outfile;
FILE    *errfile;
FILE    *logfile;
FILE    *infile;
long    nvar = 0;
GEN     gnil,gzero,gun,gdeux,ghalf,polvar,gi,RAVYZARC;
GEN     gpi=(GEN)0;
GEN     geuler=(GEN)0;
GEN     bernzone=(GEN)0;
GEN     premierbloc=(GEN)0;
entree  **varentries, **hashtable;
GEN     *polun, *polx, *g;
long    *ordvar,varchanged=0;
long    nextbloc = 0;
#ifdef LONG_IS_64BIT
long    glbfmt[]={'g',0,38};
#else
long    glbfmt[]={'g',0,28};
#endif
Rect    **rectgraph;
long    pari_randseed;
long    DEBUGLEVEL = 0;

byteptr diffptr;
GEN     primetab; /* nombres premiers prives */

#ifdef LONG_IS_64BIT
long    lontyp[30]={0,0x100000000,0x100000000,1,1,1,1,2,1,1,2,2,0,1,1,1,1,1,1,1};
long    lontyp2[30]={0,0x100000000,0x100000000,2,1,1,1,3,2,2,2,2,0,1,1,1,1,1,1,1};     
#else
long    lontyp[30]={0,0x10000,0x10000,1,1,1,1,2,1,1,2,2,0,1,1,1,1,1,1,1};
long    lontyp2[30]={0,0x10000,0x10000,2,1,1,1,3,2,2,2,2,0,1,1,1,1,1,1,1};     
#endif
void    (*printvariable)(long);
void	*foreignHandler;	/* Handler for foreign commands. */
GEN	(*foreignExprHandler)(char*); /* Handler for foreign expressions. */
char	foreignExprSwitch = 3;	/* Just some unprobable char. */
long	(*foreignAutoload)(char*, long); /* What to call if unknown
					    function is found. */
void	(*foreignFuncFree)(entree *); /* How to free external enree. */

     
     /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
     /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
     /*                                                                 */
     /*                      INITIALISATION DU SYSTEME                  */
     /*                                                                 */
     /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
     /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

const long dummy=0; /* Ne pas enlever !!! */

void
catchinterrupt(int truc)
{
  truc = truc; /* Only way to keep quiet all ANSI compilers... */
  if (infile != stdin) switchin(NULL);
  signal(SIGINT,catchinterrupt);
  err(interrupter);
}

#ifdef SIGBUS
void
catchbus(int truc)
{
  truc = truc; 
  signal(SIGBUS,catchbus);
  err(talker,"bus error: bug in GP or in calling program");
}
#endif

void
catchsegv(int truc)
{
  truc = truc; 
  signal(SIGSEGV,catchsegv);
  err(talker,"segmentation fault: bug in GP or in calling program");
}

#ifdef __cplusplus
void
init(long parisize, long maxprime, void (*printvar)(long))
#else
void
init(long parisize, long maxprime)
#endif
{
  long v, n;
  char *p;
  GEN p1;
  Rect *e;
  
#ifdef __cplusplus
  printvariable=printvar;
#else
  printvariable=printvargp;
#endif

  outfile = stdout;errfile = stderr;logfile = NULL;infile = stdin;
  if (INIT_JMP && setjmp(environnement))
  {
    fprintferr( "\n  ***   Error in the PARI system. End of the program.\n");
    exit(1);
  }
  if (INIT_SIG) {
    signal(SIGINT,catchinterrupt);
#ifdef SIGBUS
    signal(SIGBUS,catchbus);
#endif
    signal(SIGSEGV,catchsegv);
  }
  compact_arrays=1;
  pari_randseed=1;v=parisize&(BYTES_IN_LONG-1);
  if(v) parisize+=(BYTES_IN_LONG-v);
  if (!(diffptr=initprimes(maxprime))) err(memer);
#if __MWERKS__
  {
    OSErr resultCode; Handle newHand = MFTempNewHandle(parisize,&resultCode);
    if (!newHand) err(memer);
    HLock(newHand);
    bot=(long)*newHand;
  }
#else
  if (!(bot=(long)malloc(parisize))) err(memer);
#endif
  top=avma=bot+parisize;
  if (!(varentries=(entree **)malloc(sizeof(entree*)*MAXVAR))) err(memer);
  if (!(hashtable=(entree **)malloc(sizeof(entree*)*TBLSZ))) err(memer);
  if (!(ordvar=(long *)malloc(sizeof(long)*MAXVAR))) err(memer);
  if (!(polun=(GEN *)malloc(sizeof(GEN)<<MAXSHIFTVAR))) err(memer);
  if (!(polx=(GEN *)malloc(sizeof(GEN)<<MAXSHIFTVAR))) err(memer);
  if (!(g=(GEN *)malloc(sizeof(GEN)*STACKSIZE))) err(memer);
  if (!(rectgraph=(Rect**)malloc(sizeof(Rect*)*NUMRECT))) err(memer);
  for(n=0;n<NUMRECT;n++) 
  {
    if(!(e=rectgraph[n]=(Rect*)malloc(sizeof(Rect)))) err(memer);
    e->head=e->tail=NULL;
    e->sizex=e->sizey=0;
    e->cursorx=cgetr(3);
    e->cursory=cgetr(3);
    e->xscale=cgetr(3);
    e->xshift=cgetr(3);
    e->yscale=cgetr(3);
    e->yshift=cgetr(3);
  }
  for(n = 0; n < TBLSZ; n++) hashtable[n] = NULL;
  for(v = 0; v < NUMFUNC; v++)
  {
    fonctions[v].help = helpmessage[v];
    fonctions[v].valence |= EpSTATIC;
    for(n = 0, p = fonctions[v].name; *p; p++) n = n << 1 ^ *p;
    if (n < 0) n = -n; n %= TBLSZ;
    fonctions[v].next = hashtable[n];
    hashtable[n] = fonctions + v;
  }
  gnil = cgeti(2);gnil[1]=2; setpere(gnil,MAXPERE);
  gzero = cgeti(2);gzero[1]=2; setpere(gzero, MAXPERE);
  gun = stoi(1); setpere(gun, MAXPERE);
  gdeux = stoi(2); setpere(gdeux, MAXPERE);
  ghalf = cgetg(3,4);ghalf[1]=un;ghalf[2]=deux; setpere(ghalf, MAXPERE);
  gi = cgetg(3,6); gi[1] = zero; gi[2] = un; setpere(gi, MAXPERE);
  p1=cgetg(4,10);p1[1]=evalsigne(1)+evalvarn(MAXVARN)+evallgef(4);
  p1[2]=zero;p1[3]=un;polx[MAXVARN]=p1;
  p1=cgetg(3,10);p1[1]=evalsigne(1)+evalvarn(MAXVARN)+evallgef(3);
  p1[2]=un;polun[MAXVARN]=p1;
  for(v=0; v < MAXVAR; v++) ordvar[v] = v;
  polvar = cgetg(MAXVAR + 1,17); setlg(polvar,1); setpere(polvar, MAXPERE);
  for(v=1;v<=MAXVAR;v++) polvar[v]=evaltyp(17)+evalpere(MAXPERE)+evallg(1);
  primetab = cgetg(NUMPRTBELT+2,17);
  for(v = 1; v <= NUMPRTBELT+1; v++) primetab[v]=un;
  for(v = 0; v < STACKSIZE; v++) g[v] = gzero;
  lisseq("x");avloc=avma;
}

void
killall()
{
  long i,n;
  char *p;
  Rect *e;
  
  for(n=0;n<NUMRECT;n++) {e=(rectgraph[n]);if(RHead(e)) killrect(n);}
  for(i=1;i<=NUMPRTBELT;i++)
    if(!gcmp1((GEN)primetab[i])) {killbloc((GEN)primetab[i]);primetab[i]=un;}
  premierbloc=gpi=geuler=bernzone=(GEN)0;
  for(i=0;i<STACKSIZE;i++) g[i]=gzero;
  for(i=0;i<MAXVAR;i++) ordvar[i]=i;
  setlg(polvar,1);
  for(i=1;i<=MAXVAR;i++) polvar[i]=evaltyp(17)+evalpere(MAXPERE)+evallg(1);  
  for(n=0;n<TBLSZ;n++) hashtable[n]=NULL;
  for(i=0;i<NUMFUNC;i++)
  {
    for(n=0,p=fonctions[i].name;*p;p++) n=n<<1^*p;
    if (n<0) n=-n;n%=TBLSZ;
    fonctions[i].next=hashtable[n];
    hashtable[n]=fonctions+i;
  }
  precdl=16;defaultpadicprecision=16;
  pariecho=nvar=varchanged=nextbloc=0;
  lisseq("x");
}  

void
freeall()
{
  int n;
  long i;
  Rect *e;
  entree *ep, *ep1;
  
  for(n=0;n<NUMRECT;n++)
  {
    e=(rectgraph[n]);if(RHead(e)) killrect(n);
    free((void *)e);
  }
  free((void *)rectgraph);
  free((void *)g);free((void *)polx);
  free((void *)polun);free((void *)ordvar);
  while(premierbloc) killbloc(premierbloc);

  for(i = 0; i < TBLSZ; i++)
    for(ep = hashtable[i]; ep; ep = ep1)
    {
      int m = ep - fonctions;
      ep1 = ep->next;
      if (ep->help && !(ep->valence & EpSTATIC)) {
	  free(ep->help);
	  ep->help = NULL;	/* To avoid doublicate freeing in freeep. */
      }
      if (m < 0 || m >= NUMFUNC) freeep(ep);
    }
   
  free((void *)hashtable);free((void *)varentries);free((void *)bot);
  free((void *)diffptr);
}
  
GEN
geni(void)
{
  return gi;
}

long
marklist(void)
{
  return nextbloc;
}

GEN
newbloc(long n)
{
  long *x;
  x = (long *)malloc((n << TWOPOTBYTES_IN_LONG) + 4*BYTES_IN_LONG);
  if (!x) err(memer);
  *x++ = 0;
  *x++ = (long)premierbloc;
  *x++ = nextbloc++;
  *x++ = 0;
  if (premierbloc) premierbloc[-4] = (long)x;
  return premierbloc = x;
}

void
killbloc(GEN x)
{
  if (!x || isonstack(x)) return;
  if (x[-4]) ((GEN)x[-4])[-3] = x[-3]; else premierbloc = (GEN)x[-3];
  if (x[-3]) ((GEN)x[-3])[-4] = x[-4];
  free((void *)(x-4));
}

void
newvalue(entree *ep, GEN val)
{
  GEN y = gclone(val);
  y[-1] = (long) ep->value;
  ep->value = (void *)y;
}

void
changevalue(entree *ep, GEN val)
{
  GEN y = gclone(val);
  GEN x = (GEN)ep->value;
  ep->value = (void *)y;
  if ((long)x - (long)ep == sizeof(entree)) 
  {
    y[-1] = (long)x;
    return;
  }
  y[-1] = x[-1];
  killbloc(x);
}

void
killvalue(entree *ep)
{
  GEN x = (GEN)ep->value;
  if ((long)x - (long)ep == sizeof(entree)) return;
  ep->value = (void *)x[-1];
  killbloc(x);
}


entree *
install(void *f, char *name, int valence)
{
  if ((valence < 0) || (valence > 3)) err(valencer1);
  return installep(f, name, valence, NULL, NULL);
}

entree *
installep(void *f, char *name, int valence, char *code, char *help)
{
  int n;
  entree *ep;
  char *p;
  
  for(n = 0, p = name; *p; p++) n = n << 1 ^ *p;
  if (n < 0) n = -n; n %= TBLSZ;
  for(ep = hashtable[n]; ep; ep = ep->next)
    if (!strcmp(name, ep->name)) err(nomer1);
  ep = (entree *)malloc(sizeof(entree) + strlen(name) + 1);
  ep->name = (char *)ep + sizeof(entree); strcpy(ep->name, name);
  ep->value = (void *)f;
  ep->valence = valence;
  ep->menu = 0;
  ep->next = hashtable[n];
  ep->help = help;
  ep->code = code;
  hashtable[n] = ep;
  return ep;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*              TRAITEMENT DES ERREURS                             */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef __GNUC__
__volatile__
#endif
void
err(long numerr, ...)
{
  char c;
  va_list poer;
  char *ch;
  GEN noninv;

  va_start(poer,numerr);
  if (numerr != memer) {
    flusherr();fflush(outfile);
  } else {
    pariErr=NULL;		/* Try to avoid loop */
  }
  fprintferr( "  ***   %s",errmessage[numerr]);
  switch (numerr)
  {
    case matcher1:
      ch=va_arg(poer, char*);c = *ch++;
      fprintferr( "'%c'\n  ***   instead of: '%s'", c, ch); break;
    case impl: ch=va_arg(poer, char*);
      fprintferr( " %s is not yet implemented.",ch); break;
    case talker: ch=va_arg(poer, char*);
      fprintferr( "%s.",ch); break;
    case varer1:
    case unknowner1:
    case caracer1: ch=va_arg(poer, char*);fprintferr( "'%s'",ch);break;
    case invmoder: ch=va_arg(poer, char*);noninv=va_arg(poer, GEN);
      fprintferr(": ");bruterr(noninv,'g',-1);break;
    case errpile: 
      if (!pariErr || !pariErr->putc) {
	fprintferr("\n");
      }
      allocatemoremem(0);break;
  }
  if (pariErr && pariErr->die) /* empty */;
  else flusherr();
  outfile=stdout;errfile=stderr;
  if (!pariErr || !pariErr->putc) fprintferr("\n");
  va_end(poer);
  if (pariErr && pariErr->die) pariErr->die();
  if (environnement) {
    longjmp(environnement, numerr);
  } else {
    exit(1);
  }
}

void
recover(long listloc)
{
  long  m, n;
  GEN x, y;
  entree *ep, *ep2;

  for (n = 0; n < TBLSZ; n++)
    for (ep = hashtable[n]; ep;)
      if (EpVALENCE(ep) >= 100)
      {
        x = (GEN)ep->value;
        if ((long)x - (long)ep == sizeof(entree))
        {
          if (EpVALENCE(ep) == 200) ep = ep->next;
          else
            if (ep == hashtable[n])
            {
              hashtable[n] = ep->next;
              freeep(ep);
              ep = hashtable[n];
            }
            else
            {
              for(ep2 = hashtable[n]; ep2->next != ep; ep2 = ep2->next);
              ep2->next = ep->next;
              freeep(ep); ep = ep2->next;
            }
          continue;
        }
	m = x[-2];
        if ((m < listloc) || (m >= nextbloc)) ep=ep->next;
        else killvalue(ep);
      }
      else ep = ep->next;
  for(x = premierbloc; x && x[-2] >= listloc; x = y)
  {
    y = (GEN)x[-3];
    if (x != gpi && x != geuler) killbloc(x);
  }
}

void
allocatemoremem(unsigned long newsize)
{
  long av,declg,declg2,tl,parisize,v;
  GEN ll,pp,l1,l2,l3;
  unsigned long topold,avmaold,botold;

  avmaold=avloc;topold=top;botold=bot;
  if(newsize<3) parisize=(topold-botold)<<1;
  else 
  {
    if(newsize<(topold-avmaold)) 
      err(talker,"required stack memory too small");
    else parisize=newsize+16-(((newsize-1)&15)+1);
  }
  if (!(bot=(long)malloc(parisize))) err(nomer2);
  if(!newsize)
  {
    fprintferr( "  *** Warning: doubling the stack size; new stack = %ld\n",parisize);
    fprintferr( "  *** Please reissue the same command if you are under GP\n");
  }
  top=avma=bot+parisize;
  declg=(long)top-(long)topold;declg2=declg>>TWOPOTBYTES_IN_LONG;
  for(ll=(GEN)top,pp=(GEN)topold;pp>(GEN)avmaold;) *--ll= *--pp;
  av=(long)ll;
  while(ll<(GEN)top)
  {
    l2=ll+lontyp[tl=typ(ll)];
    if(tl==10) {l3=ll+lgef(ll);ll+=lg(ll);if(l3>ll) l3=l2;}
    else {ll+=lg(ll);l3=ll;} 
    for(;l2<l3;l2++) 
    {
      l1=(GEN)(*l2);
      if((l1<(GEN)topold)&&(l1>=(GEN)avmaold)) *l2+=declg;
    }
  }
  gnil+=declg2;gzero+=declg2;gun+=declg2;gdeux+=declg2;ghalf+=declg2;
  gi+=declg2;polx[MAXUBYTE]+=declg2;polun[MAXUBYTE]+=declg2;polvar+=declg2;
  for(v=0;v<=tglobal;v++)
    if((g[v]<(GEN)topold)&&(g[v]>=(GEN)avmaold)) g[v]+=declg2;
  free((void *)botold);avloc=avma=av;
}

GEN
allocatemem(unsigned long newsize)
{
#if __MWERKS__
  newsize = newsize;
  err(talker, "Not implemented in this version, sorry");
#else
  allocatemoremem(newsize);
  longjmp(environnement,errpile);
#endif
  return gnil; /*inutile mais ca fait plaisir a des compilos */
}

#if __MWERKS__
void *macrealloc(void *p, size_t oldsize, size_t newsize)
{
  char *q = malloc(newsize);
  char *qq = q, *pp = p;
  int l = newsize > oldsize ? oldsize : newsize;
  while (l--) *qq++ = *pp++;
  free(p);
  return q;
}
#endif

/* In case if dynamic linking may lead to several mallocs defined. */

#ifdef MALLOC_PROCS
#undef malloc
#undef realloc
#undef free
Malloc_procs malloc_procs = {malloc, realloc, free};
#endif
