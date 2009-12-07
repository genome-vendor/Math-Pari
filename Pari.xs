#ifdef __cplusplus
extern "C" {
#endif 

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#ifdef __cplusplus
}
#endif 

/* This should not be defined at this moment, but in 5.001n is. */
#ifdef coeff
#  undef coeff
#endif

#include <pari.h>
#include <language/anal.h>
#include <gp/gp.h>			/* init_opts */

/* 	$Id: Pari.xs,v 1.7 1995/01/23 18:50:58 ilya Exp ilya $	 */
/* dFUNCTION should be the last declaration! */

#ifdef __cplusplus
  #define VARARG ...
#else
  #define VARARG
#endif

#define dFUNCTION(retv)  retv (*FUNCTION)(VARARG) = \
            (retv (*)(VARARG)) XSANY.any_dptr

#define DO_INTERFACE(inter) subaddr = CAT2(XS_Math__Pari_interface, inter)
#define CASE_INTERFACE(inter) case inter: \
                   DO_INTERFACE(inter); break

/* Here is the rationals for managing SVs which keep GENs: when newly
   created SVs from GENs on stack, the same moved to heap, and
   originally from heap. We assume that we do not need to free stuff
   that was originally on heap. However, we need to free the stuff we
   moved from the stack ourself.
   
   Here is how we do it: The variables that were initially off stack
   have SvPVX == GENheap. 
   
   The variables that were moved from the stack have SvPVX ==
   GENmovedOffStack.

   If the variable is on stack, and it is the oldest one which is on
   stack, then SvPVX == GENfirstOnStack.

   Otherwise SvPVX is the next older SV that refers to a GEN on stack.

   In the last two cases SvCUR is the offset on stack of the stack
   frame on the entry into the function for which SV is the argument.
*/

#define GENmovedOffStack ((char*) 1) /* Just an atom. */
#define GENfirstOnStack ((char*) 2) /* Just an atom. */
#define GENheap NULL
#define ifact mpfact

typedef entree * PariVar;		/* For loop variables. */
typedef entree * PariName;		/* For changevalue.  */
typedef char * PariExpr;
typedef GEN * GEN_Ptr;

/* We make a "fake" PVAV, not enough entries.  */

#define setSVpari(sv, in, oldavma) do {				\
    sv_setref_pv(sv, "Math::Pari", (void*)in);			\
    if (typ(in) >= 17 && SvTYPE(SvRV(sv)) != SVt_PVAV) {	\
	make_PariAV(sv);					\
    }								\
    if (isonstack(in)) {					\
	SvCUR(SvRV(sv)) = oldavma - bot;			\
	SvPVX(SvRV(sv)) = (char*)PariStack;			\
	PariStack = SvRV(sv);					\
	perlavma = avma;					\
	onStack_inc;						\
    }								\
    SVnum_inc;							\
} while (0)

SV* PariStack;			/* PariStack keeps the latest SV that
				 * keeps a GEN on stack. */
long perlavma;				/* How much stack is needed
					   for GENs in Perl variables. */
long sentinel;				/* How much stack was used
					   when Pari called Perl. */

#ifdef DEBUG_PARI

long SVnum;
long SVnumtotal;
long onStack;
long offStack;

#  define SVnum_inc (SVnum++, SVnumtotal++)
#  define SVnum_dec (SVnum--)
#  define onStack_inc (onStack++)
#  define onStack_dec (onStack--)
#  define offStack_inc (offStack++)
#else  /* !defined DEBUG_PARI */
#  define SVnum_inc 
#  define SVnum_dec
#  define onStack_inc
#  define onStack_dec
#  define offStack_inc
#endif /* !defined DEBUG_PARI */

void
make_PariAV(SV *sv)
{
    AV *av = (AV*)SvRV(sv);
    char *s = SvPVX(av);
    IV i = SvIVX(av);
    SvUPGRADE((SV*)av, SVt_PVAV);    
    SvPVX(av)	= s;
    SvIVX(av)	= i;
    sv_magic((SV*)av, sv, 'P', Nullch, 0);
}

SV*
wrongT(SV *sv, char *file, int line)
{
  if (SvTYPE(sv) != SVt_PVCV || SvTYPE(sv) != SVt_PVGV) {
    croak("Got the type 0x%x instead of CV=0x%x or GV=0x%x in %s, %i", 
	  SvTYPE(sv), SVt_PVCV, SVt_PVGV, file, line);      
  } else {
    croak("Something very wrong  in %s, %i", file, line);
  }
  return NULL;				/* To pacify compiler. */
}

HV *pariStash;				/* For quick id. */
HV *pariEpStash;

/* Copied from anal.c. */
static entree *
installep(void *f, char *name, int len, int valence, int add, entree **table)
{
  entree *ep = (entree *) gpmalloc(sizeof(entree) + add + len+1);
  const entree *ep1 = initial_value(ep);
  char *u = (char *) ep1 + add;

  ep->name    = u; strncpy(u, name,len); u[len]=0;
  ep->args    = NULL; ep->help = NULL; ep->code = NULL;
  ep->value   = f? f: (void *) ep1;
  ep->next    = *table;
  ep->valence = valence;
  ep->menu    = 0;
  return *table = ep;
}
static void
changevalue(entree *ep, GEN val)
{
  GEN y = gclone(val), x = (GEN)ep->value;

  ep->value = (void *)y;
  if (x == (GEN) initial_value(ep) || !isclone(x))
  {
    y[-1] = (long)x; /* push new value */
    return;
  }
  y[-1] = x[-1]; /* save initial value */
  killbloc(x);   /* destroy intermediate one */
}
static long
numvar(GEN x)
{
  if (typ(x) != t_POL || lgef(x) != 4 || 
    !gcmp0((GEN)x[2]) || !gcmp1((GEN)x[3])) 
      croak("Corrupted data: should be variable");
  return varn(x);
}


static SV *
PARIvar(char *s)
{
  char *olds = s, *u, *v;
  long hash;
  SV *sv;
  GEN p1;
  entree *ep = is_entry_intern(s, functions_hash, &hash);

  if (ep) {
      if (EpVALENCE(ep) != EpVAR)
	  croak("Got a function name instead of a variable");
  } else {
      ep = installep(NULL, s, strlen(s), EpVAR, 7*sizeof(long),
		     functions_hash + hash);
      manage_var(0,ep);
#if 0
      ep = (entree *)malloc(sizeof(entree) + 7*BYTES_IN_LONG 
			    + s - olds + 1);
      ep->name = (char *)ep + sizeof(entree) + 7*BYTES_IN_LONG;
      for (u = ep->name, v = olds; v < s;) *u++ = *v++; *u = 0;
      ep->value = (void *)((char *)ep + sizeof(entree));
      ep->code = ep->help = NULL;
      ep->next = hashtable[n];
      hashtable[n] = ep;
      p1 = (GEN)ep->value;
      if (nvar == MAXVAR) err(trucer1);
      ep->valence = 200;
      p1[0] = evaltyp(10)+evalpere(1)+evallg(4);
      p1[1] = evalsigne(1)+evallgef(4)+evalvarn(nvar);
      p1[2] = zero; p1[3] = un;
      polx[nvar] = p1;
      polvar[nvar+1] = (long)p1;
      p1 += 4;
      p1[0] = evaltyp(10)+evalpere(1)+evallg(3);
      p1[1] = evalsigne(1)+evallgef(3)+evalvarn(nvar); p1[2] = un;
      polun[nvar] = p1;
      varentries[nvar++] = ep;
      setlg(polvar, nvar+1);    
#endif
  }
  
  found:
  sv = NEWSV(909,0);
  sv_setref_pv(sv, "Math::Pari::Ep", (void*)ep);
  make_PariAV(sv);
  return sv;
}

static entree *
findVariable(SV *sv, int generate)
{
    /* There may be 4 important cases:
       a) we got a 'word' string, which we interpret as the name of
          the variable to use;
       b1) It is a pari value containing a polynomial 0+1*v, we use it;
       b2) It is other pari value, we ignore it;
       c) it is a string containing junk, same as 'b';
       d) It is an ep value => typo (same iterator in two loops).
       In any case we localize the value.
     */
  char *s;
  char *s1, *u, *v;
  long hash;
  GEN p1;
  entree *ep;
  char name[50];

  if (SvROK(sv)) {
      SV* tsv = SvRV(sv);
      if (SvOBJECT(tsv)) {
	  if (SvSTASH(tsv) == pariStash) {
	      GEN x = (GEN)SvIV(tsv);
	      if (typ(x)==10 		/* Polynomial. */
		  && lgef(x)==4		/* 2 terms */
		  && (gcmp0((GEN)x[2]))	/* Free */
		  && (gcmp1((GEN)x[3]))) { /* Leading */
		  s = varentries[ordvar[varn(x)]]->name;
		  goto repeat;
	      }
	      goto ignore;
	  } else if (SvSTASH(tsv) == pariEpStash) {
	      /* It is not good to croak: $v=PARIvar 'v'; vector(3,$v,'v'); */
	      if (generate)
		  /*croak("Same iterator in embedded PARI loop construct")*/;
	      return (entree*) SvIV(tsv);
	  }
      }
  }
  if (!SvOK(sv))
      goto ignore;
  s = SvPV(sv,na);
  repeat:
  s1 = s;
  while (isalnum(*s1)) 
      s1++;
  if (*s1 || s1 == s || !isalpha(*s)) {
      static int depth;

    ignore:
      if (!generate)
	  croak("Bad PARI variable name \"%s\" specified",s);
      SAVEINT(depth);
      sprintf(name, "intiter%i",depth++);
      s = name;
      goto repeat;
  }
  
  ep = is_entry_intern(s, functions_hash, &hash);

  if (ep) {
      if (EpVALENCE(ep) != EpVAR)
	  croak("Got a function name instead of a variable");
  } else {
      ep = installep(NULL, s, s1 - s, EpVAR, 7*sizeof(long),
		     functions_hash + hash);
      manage_var(0,ep);
  }

#if 0
  olds = s;
  for (n = 0; isalnum(*s); s++) n = n << 1 ^ *s;
  if (n < 0) n = -n; n %= TBLSZ;
  for(ep = hashtable[n]; ep; ep = ep->next)
  {
    for(u = ep->name, v = olds; (*u) && *u == *v; u++, v++);
    if (!*u && !*v) {
      if (EpVALENCE(ep) != 200)
	croak("Got a function name instead of a variable");
      return ep;
    }
  }
  ep = (entree *)malloc(sizeof(entree) + 7*BYTES_IN_LONG 
			+ s - olds + 1);
  ep->name = (char *)ep + sizeof(entree) + 7*BYTES_IN_LONG;
  for (u = ep->name, v = olds; v < s;) *u++ = *v++; *u = 0;
  ep->value = (void *)((char *)ep + sizeof(entree));
  ep->code = ep->help = NULL;
  ep->next = hashtable[n];
  hashtable[n] = ep;
  p1 = (GEN)ep->value;
  if (nvar == MAXVAR) err(trucer1);
  ep->valence = 200;
  p1[0] = evaltyp(10)+evalpere(1)+evallg(4);
  p1[1] = evalsigne(1)+evallgef(4)+evalvarn(nvar);
  p1[2] = zero; p1[3] = un;
  polx[nvar] = p1;
  polvar[nvar+1] = (long)p1;
  p1 += 4;
  p1[0] = evaltyp(10)+evalpere(1)+evallg(3);
  p1[1] = evalsigne(1)+evallgef(3)+evalvarn(nvar); p1[2] = un;
  polun[nvar] = p1;
  varentries[nvar++] = ep;
  setlg(polvar, nvar+1);
#endif
  return ep;
}

static PariVar
bindVariable(SV *sv)
{
    /* There may be 4 important cases:
       a) we got a 'word' string, which we interpret as the name of
          the variable to use;
       b1) It is a pari value containing a polynomial 0+1*v, we use it;
       b2) It is other pari value, we ignore it;
       c) it is a string containing junk, same as 'b';
       d) It is an ep value => typo (same iterator in two loops).
       In any case we localize the value.
     */
  char *s;
  char *olds, *u, *v;
  long n, override = 0;
  GEN p1;
  entree *ep;
  char name[50];

  if (!SvREADONLY(sv)) {
      save_item(sv);			/* Localize it. */
      override = 1;
  }
  ep = findVariable(sv, 1);
  if (override) {
      sv_setref_pv(sv, "Math::Pari::Ep", (void*)ep);
      make_PariAV(sv);
  }
  return ep;
}

static int
not_here(char *s)
{
    croak("%s not implemented on this architecture", s);
    return -1;
}

SV* worksv;
SV* workErrsv;

void
svputc(char c)
{
  sv_catpvn(worksv,&c,1);
}


void
svputs(char* p)
{
  sv_catpv(worksv,p);
}

void
svErrputc(char c)
{
  sv_catpvn(workErrsv,&c,1);
}


void
svErrputs(char* p)
{
  sv_catpv(workErrsv,p);
}

void
svOutflush()
{
    /* EMPTY */
}

void
svErrflush()
{
    STRLEN len;

    if (SvPV(workErrsv, len) && len) {
	warn("PARI: %s",SvPV(workErrsv,na));
	sv_setpv(workErrsv,"");
    }
}

void
svErrdie()
{
  SV *errSv = newSVsv(workErrsv);

  sv_setpv(workErrsv,"");
  croak("PARI: %s", SvPV(errSv,na));
}


PariOUT perlOut={svputc, svputs, svOutflush, NULL};
PariOUT perlErr={svErrputc, svErrputs, svErrflush, svErrdie};

GEN
sv2pari(SV* sv)
{
  if (SvROK(sv)) {
      SV* tsv = SvRV(sv);
      if (SvOBJECT(tsv)) {
	  if (SvSTASH(tsv) == pariStash) {
	      IV tmp = SvIV(tsv);
	      return (GEN) tmp;
	  } else if (SvSTASH(tsv) == pariEpStash) {
	      IV tmp = SvIV(tsv);
	      return (GEN)(((entree*) tmp)->value);
	  }
      }
      {
	  int type = SvTYPE(tsv);
	  if (type==SVt_PVAV) {
	      AV* av=(AV*) tsv;
	      I32 len=av_len(av);	/* Length-1 */
	      GEN ret=cgetg(len+2,17);
	      int i;
	      for (i=0;i<=len;i++) {
		  SV** svp=av_fetch(av,i,0);
		  if (!svp) croak("Internal error in perl2pari!");
		  ret[i+1]=(long)sv2pari(*svp);
	      }
	      return ret;
	  } else {
	      return lisexpr(SvPV(sv,na)); /* For overloading */
	  }
      }
  }
  else if (SvIOKp(sv)) return stoi(SvIV(sv));
  else if (SvNOKp(sv)) return dbltor(SvNV(sv));
  else if (SvPOK(sv)) return lisexpr(SvPV(sv,na));
  else if (!SvOK(sv)) return stoi(0);
  else
    croak("Variable in perl2pari is not of known type");  
}

GEN
sv2parimat(SV* sv)
{
  GEN in=sv2pari(sv);
  if (typ(in)==17) {
    long len=lg(in)-1;
    long t;
    long l=lg(in[1]);
    for (;len;len--) {
      if ((t=typ(in[len])) == 17) {
	settyp(in[len],18);
      } else if (t!=18) {
	croak("Not a vector where column of a matrix expected");
      }
      if (lg(in[len])!=l) {
	croak("Columns of input matrix are of different height");
      }
    }
    settyp(in,19);
    return in;
  } else if (typ(in)==19) {
    return in;
  } else {
    croak("Not a matrix where matrix expected");
  }
}


SV*
pari2iv(GEN in)
{
  return newSViv((IV)gtolong(in));
}

SV*
pari2nv(GEN in)
{
  return newSVnv(gtodouble(in));
}

SV*
pari2pv(GEN in)
{
  PariOUT *oldOut = pariOut;
  pariOut = &perlOut;
  worksv = newSVpv("",0);
  bruteall(in,'g',-1,0);		/* 0: compact pari-readable form */
  pariOut = oldOut;
  return worksv;
}

#ifdef LONG_IS_64BIT
#define fmt_nb 38
#else
#define fmt_nb 28
#endif

SV*
pari_print(GEN in)
{
  PariOUT *oldOut = pariOut;
  pariOut = &perlOut;
  worksv = newSVpv("",0);
  brute(in, 'g', fmt_nb);
  pariOut = oldOut;
  return worksv;
}

SV*
pari_pprint(GEN in)
{
  PariOUT *oldOut = pariOut;
  pariOut = &perlOut;
  worksv = newSVpv("",0);
  sor(in, 'g'/*fmt.format*/, fmt_nb, 0/*fmt.field*/);
  pariOut = oldOut;
  return worksv;
}

SV*
pari_texprint(GEN in)
{
  PariOUT *oldOut = pariOut;
  pariOut = &perlOut;
  worksv = newSVpv("",0);
  texe(in, 'g', fmt_nb);
  pariOut = oldOut;
  return worksv;
}

/* Should be used for conversion of procedure arguments only. */
SV*
pari2mortalsv(GEN in, long oldavma)
{				/* Oldavma should keep the value of
				 * avma when entering a function call. */
    SV *sv = sv_newmortal();

    setSVpari(sv, in, oldavma);
    return sv;
}

static const 
unsigned char defcode[] = "\06xD0,G,D0,G,D0,G,D0,G,D0,G,D0,G,";

static int doing_PARI_autoload = 0;

entree *
installPerlFunctionCV(SV* cv, char *name, I32 numargs, char *help)
{
    char *code, *s;
    I32 req = numargs, opt = 0;
    entree *ep;

    if(SvROK(cv))
	cv = SvRV(cv);

    if (numargs < 0 && SvPOK(cv) && (s = SvPV(cv,na))) {
	/* Get number of arguments. */
	req = opt = 0;
	while (*s == '$')
	    req++, s++;
	if (*s == ';') 
	    s++;
	while (*s == '$')
	    opt++, s++;
	if (*s == '@') {
	    opt += 6;			/* Max 6 optional arguments. */
	    s++;
	}
	if (*s == 0) {			/* Got it! */
	    numargs = req + opt;
	}
    }
    
    if (numargs < 0) {		/* Variable number of arguments. */
	/* Install something hairy with <= 6 args */
	code = (char*)defcode + 1;		/* Remove constness. */
	numargs = code[-1];
    } else if (numargs >= 256) {
	croak("Import of Perl function with too many arguments");
    } else {
	code = (char *)malloc(numargs*6 - req*5 + 2);
	code[0] = 'x';
	memset(code + 1, 'G', req);
	s = code + 1 + req;
	while (opt--) {
	    strcpy(s, "D0,G,");
	    s += 6;
	}
	*s = '\0';
    }
    ((CV*)cv)->sv_any->xof_off = numargs;	/* XXXX Nasty of us... */
    SAVEINT(doing_PARI_autoload);
    doing_PARI_autoload = 1;
    ep = install((void*)SvREFCNT_inc(cv), name, code);
    doing_PARI_autoload = 0;
    if (code != (char*)defcode + 1)
	free(code);
    ep->help = help;
    return ep;
}

void
freePerlFunction(entree *ep)
{
    if (!ep->code || (*ep->code != 'x')) {
	croak("Attempt to ask Perl to free PARI function not installed from Perl");
    }
    if (ep->code != (char *)defcode + 1)
	free(ep->code - 1);
    if (ep->help)
	free(ep->help);
    SvREFCNT_dec((SV*)ep->value);
}

long
moveoffstack_newer_than(SV* sv)
{
  SV* sv1;
  SV* nextsv;
  long ret=0;
  
  for (sv1 = PariStack; sv1 != sv; sv1 = nextsv) {
    ret++;
    nextsv = (SV *) SvPVX(sv1);
    SvPVX(sv1) = GENmovedOffStack; /* Mark as moved off stack. */
    SvIVX(sv1) = (IV) gclone((GEN)SvIV(sv1));
    onStack_dec;
    offStack_inc;
  }
  PariStack = sv;
  return ret;
}

GEN
callPerlFunction(entree *ep, ...)
{
    va_list args;
    char *s = ep->code;
    SV *cv = (SV*) ep->value;
    int numargs = ((CV*)cv)->sv_any->xof_off;	/* XXXX Nasty of us... */
    GEN res;
    int i;
    dSP;
    int count ;
    long oldavma = avma;
    SV *oPariStack = PariStack;
    SV *sv;

    va_start(args, ep);
    ENTER ;
    SAVETMPS;
    SAVEINT(sentinel);
    sentinel = avma;
    PUSHMARK(sp);
    EXTEND(sp, numargs + 1);
    for (i = 0; i < numargs; i++) {
	PUSHs(pari2mortalsv(va_arg(args, GEN), oldavma));
    }
    va_end(args);
    PUTBACK;
    count = perl_call_sv(cv, G_SCALAR);

    SPAGAIN;
    if (count != 1)
	croak("Perl function exported into PARI did not return a value");

    sv = SvREFCNT_inc(POPs);		/* Preserve the guy. */

    PUTBACK ;
    FREETMPS ;
    LEAVE ;
    /* Now PARI data created inside this subroutine sits above
       oldavma, but the caller is going to unwind the stack: */
    if (PariStack != oPariStack)
	moveoffstack_newer_than(oPariStack);
    /* Now, when everything is moved off stack, and avma is reset, we
       can get the answer: */
    res = sv2pari(sv);			/* XXXX When to decrement the count? */
    /* We need to copy it back to stack, otherwise we cannot decrement
     the count.  XXXX not necessary! */
    avma -= taille(res)<<TWOPOTBYTES_IN_LONG;
    brutcopy(res, (GEN)avma);
    SvREFCNT_dec(sv);
    
    return (GEN)avma;
}

/* Currently with <=6 arguments only! */

entree *
autoloadPerlFunction(char *s, long len)
{
    CV *cv;
    SV* name;
    HV* converted;

    if (doing_PARI_autoload)
	return 0;
    converted = perl_get_hv("Math::Pari::converted",TRUE);
    if (hv_fetch(converted, s, len, FALSE)) 
	return 0;

    name = sv_2mortal(newSVpv(s, len));

    cv = perl_get_cv(SvPVX(name), FALSE);
    if (cv == Nullcv) {
	return 0;
    }
    /* Got it! */
    return installPerlFunctionCV((SV*)cv, SvPVX(name), -1, NULL); /* -1 gives variable. */
}

GEN
exprHandler_Perl(char *s)
{
    SV* dummy;
    SV* cv = (SV*)(s - LSB_in_U32 - 
		   ((char*)&(dummy->sv_flags) - ((char*)dummy)));
    GEN res;
    long count;
    dSP;
    SV *sv;
    SV *oPariStack = PariStack;

    ENTER ;
    SAVETMPS;
    PUSHMARK(sp);
    SAVEINT(sentinel);
    sentinel = avma;
    count = perl_call_sv(cv, G_SCALAR);

    SPAGAIN;
    sv = SvREFCNT_inc(POPs);		/* Preserve it through FREETMPS */

    PUTBACK ;
    FREETMPS ;
    LEAVE ;

    /* Now PARI data created inside this subroutine sits above
       oldavma, but the caller is going to unwind the stack: */
    if (PariStack != oPariStack)
	moveoffstack_newer_than(oPariStack);
    /* Now, when everything is moved off stack, and avma is reset, we
       can get the answer: */
    res = sv2pari(sv);
    /* We need to copy it back to stack, otherwise we cannot decrement
     the count. */
    avma -= taille(res)<<TWOPOTBYTES_IN_LONG;
    brutcopy(res, (GEN)avma);
    SvREFCNT_dec(sv);
    
    return (GEN)avma;
}


static GEN
Arr_FETCH(GEN g, I32 n) 
{
    I32 l = lg(g) - 1;

    if (typ(g) < 17)
	croak("Access to elements of not-a-vector");
    if (n >= l || n < 0)
	croak("Array index %i out of range", n);
    return (GEN)g[n + 1];
}

typedef int (*FUNC_PTR)();
#define set_gnuterm(a,b) set_term_funcp((FUNC_PTR)(a),(struct termentry *)(b))

MODULE = Math::Pari PACKAGE = Math::Pari PREFIX = Arr_

PROTOTYPES: ENABLE

GEN
Arr_FETCH(g,n)
    long	oldavma=avma;
    GEN g
    I32 n

MODULE = Math::Pari PACKAGE = Math::Pari

PROTOTYPES: ENABLE

GEN
sv2pari(sv)
long	oldavma=avma;
     SV *	sv

GEN
sv2parimat(sv)
long	oldavma=avma;
     SV *	sv


SV *
pari2iv(in)
long	oldavma=avma;
     GEN	in


SV *
pari2nv(in)
long	oldavma=avma;
     GEN	in


SV *
pari2num(in)
long	oldavma=avma;
     GEN	in
   CODE:
     if (typ(in)==1) {
       RETVAL=pari2iv(in);
     } else {
       RETVAL=pari2nv(in);
     }
   OUTPUT:
     RETVAL

SV *
pari2pv(in,...)
long	oldavma=avma;
     GEN	in
   CODE:
     RETVAL=pari2pv(in);
   OUTPUT:
     RETVAL

GEN
PARI(...)
long	oldavma=avma;
   CODE:
     if (items==1) {
       RETVAL=sv2pari(ST(0));
     } else {
	int i;

	RETVAL=cgetg(items+1,17);
	for (i=0;i<items;i++) {
	  RETVAL[i+1]=(long)sv2pari(ST(i));
	}
     }
   OUTPUT:
     RETVAL

GEN
PARIcol(...)
long	oldavma=avma;
   CODE:
     if (items==1) {
       RETVAL=sv2pari(ST(0));
     } else {
	int i;

	RETVAL=cgetg(items+1,17);
	for (i=0;i<items;i++) {
	  RETVAL[i+1]=(long)sv2pari(ST(i));
	}
     }
     settyp(RETVAL,18);
   OUTPUT:
     RETVAL

GEN
PARImat(...)
long	oldavma=avma;
   CODE:
     if (items==1) {
       RETVAL=sv2parimat(ST(0));
     } else {
	int i;

	RETVAL=cgetg(items+1,17);
	for (i=0;i<items;i++) {
	  RETVAL[i+1]=(long)sv2pari(ST(i));
	  settyp(RETVAL[i+1],18);
	}
     }
     settyp(RETVAL,19);
   OUTPUT:
     RETVAL

void
installPerlFunctionCV(cv, name, numargs = 1, help = NULL)
SV*	cv
char   *name
I32	numargs
char   *help

# In what follows if a function returns long, we do not need anything
# on the stack, thus we add a cleanup section.

GEN
interface0()
long	oldavma=avma;
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(prec);
  }
 OUTPUT:
   RETVAL

GEN
interface00()
long	oldavma=avma;
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION();
  }
 OUTPUT:
   RETVAL

GEN
interface1(arg1)
long	oldavma=avma;
GEN	arg1
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,prec);
  }
 OUTPUT:
   RETVAL

# with fake arguments for overloading

GEN
interface199(arg1,arg2,inv)
long	oldavma=avma;
GEN	arg1
GEN	arg2 = NO_INIT
long	inv = NO_INIT
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,prec);
  }
 OUTPUT:
   RETVAL


long
interface10(arg1)
long	oldavma=avma;
GEN	arg1
 CODE:
  {
    dFUNCTION(long);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1);
  }
 OUTPUT:
   RETVAL
 CLEANUP:
   avma=oldavma;

# With fake arguments for overloading

long
interface109(arg1,arg2,inv)
long	oldavma=avma;
GEN	arg1
GEN	arg2 = NO_INIT
long	inv = NO_INIT
 CODE:
  {
    dFUNCTION(long);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1);
  }
 OUTPUT:
   RETVAL
 CLEANUP:
   avma=oldavma;

GEN
interface11(arg1)
long	oldavma=avma;
long	arg1
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1);
  }
 OUTPUT:
   RETVAL

long
interface15(arg1)
long	oldavma=avma;
long	arg1
 CODE:
  {
    dFUNCTION(long);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1);
  }
 OUTPUT:
   RETVAL
 CLEANUP:
   avma=oldavma;

GEN
interface18(arg1)
long	oldavma=avma;
GEN	arg1
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1);
  }
 OUTPUT:
   RETVAL

GEN
interface2(arg1,arg2)
long	oldavma=avma;
GEN	arg1
GEN	arg2
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2);
  }
 OUTPUT:
   RETVAL

# With fake arguments for overloading

GEN
interface299(arg1,arg2,inv)
long	oldavma=avma;
GEN	arg1
GEN	arg2
bool	inv
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL = inv? FUNCTION(arg2,arg1): FUNCTION(arg1,arg2);
  }
 OUTPUT:
   RETVAL

long
interface20(arg1,arg2)
long	oldavma=avma;
GEN	arg1
GEN	arg2
 CODE:
  {
    dFUNCTION(long);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2);
  }
 OUTPUT:
   RETVAL
 CLEANUP:
   avma=oldavma;

# With fake arguments for overloading and comparison to gun for speed

long
interface2099(arg1,arg2,inv)
long	oldavma=avma;
GEN	arg1
GEN	arg2
bool	inv
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL = (inv? FUNCTION(arg2,arg1): FUNCTION(arg1,arg2)) == gun;
  }
 OUTPUT:
   RETVAL
 CLEANUP:
   avma=oldavma;

# With fake arguments for overloading

long
interface209(arg1,arg2,inv)
long	oldavma=avma;
GEN	arg1
GEN	arg2
bool	inv
 CODE:
  {
    dFUNCTION(long);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL = inv? FUNCTION(arg2,arg1): FUNCTION(arg1,arg2);
  }
 OUTPUT:
   RETVAL
 CLEANUP:
   avma=oldavma;

GEN
interface29(arg1,arg2)
long	oldavma=avma;
GEN	arg1
GEN	arg2
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2,prec);
  }
 OUTPUT:
   RETVAL

GEN
interface3(arg1,arg2,arg3)
long	oldavma=avma;
GEN	arg1
GEN	arg2
GEN	arg3
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2,arg3);
  }
 OUTPUT:
   RETVAL

long
interface30(arg1,arg2,arg3)
long	oldavma=avma;
GEN	arg1
GEN	arg2
GEN	arg3
 CODE:
  {
    dFUNCTION(long);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2,arg3);
  }
 OUTPUT:
   RETVAL
 CLEANUP:
   avma=oldavma;

GEN
interface4(arg1,arg2,arg3,arg4)
long	oldavma=avma;
GEN	arg1
GEN	arg2
GEN	arg3
GEN	arg4
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2,arg3,arg4);
  }
 OUTPUT:
   RETVAL

GEN
interface5(arg1,arg2,arg3,arg4,arg5)
long	oldavma=avma;
GEN	arg1
GEN	arg2
GEN	arg3
GEN	arg4
GEN	arg5
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2,arg3,arg4,prec);
  }
 OUTPUT:
   RETVAL

GEN
interface12(arg1,arg2)
long	oldavma=avma;
GEN	arg1
GEN	arg2
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,numvar(arg2), precdl);
  }
 OUTPUT:
   RETVAL


GEN
interface13(arg1, arg2=0, arg3=gzero)
long	oldavma=avma;
GEN	arg1
long	arg2
GEN	arg3
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, arg3);
  }
 OUTPUT:
   RETVAL


GEN
interface14(arg1,arg2=0)
long	oldavma=avma;
GEN	arg1
GEN	arg2
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2 ? numvar(arg2) : -1);
  }
 OUTPUT:
   RETVAL


GEN
interface21(arg1,arg2)
long	oldavma=avma;
GEN	arg1
long	arg2
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2);
  }
 OUTPUT:
   RETVAL


GEN
interface22(arg1,arg2,arg3)
long	oldavma=avma;
GEN	arg1
PariVar	arg2
PariExpr	arg3
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, arg3);
  }
 OUTPUT:
   RETVAL

GEN
interface23(arg1,arg2)
long	oldavma=avma;
GEN	arg1
long	arg2
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2);
  }
 OUTPUT:
   RETVAL

GEN
interface24(arg1,arg2)
long	oldavma=avma;
long	arg1
GEN	arg2
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2);
  }
 OUTPUT:
   RETVAL

GEN
interface25(arg1,arg2,arg3=0)
long	oldavma=avma;
GEN	arg1
GEN	arg2
long	arg3
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2,arg3);
  }
 OUTPUT:
   RETVAL

GEN
interface26(arg1,arg2,arg3)
long	oldavma=avma;
GEN	arg1
GEN	arg2
GEN	arg3
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, numvar(arg2), arg3);
  }
 OUTPUT:
   RETVAL

GEN
interface27(arg1,arg2,arg3)
long	oldavma=avma;
PariVar	arg1
GEN	arg2
PariExpr	arg3
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, arg3, prec);
  }
 OUTPUT:
   RETVAL

GEN
interface28(arg1,arg2)
long	oldavma=avma;
GEN	arg1
GEN	arg2
 CODE:
  {
    long junk;
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, &junk);
  }
 OUTPUT:
   RETVAL

long
interface29_old(arg1,arg2)
long	oldavma=avma;
GEN	arg1
long	arg2
 CODE:
  {
    dFUNCTION(long);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2);
  }
 OUTPUT:
   RETVAL
 CLEANUP:
   avma=oldavma;

GEN
interface31(arg1,arg2,arg3=0,arg4=0)
long	oldavma=avma;
GEN	arg1
GEN	arg2
GEN	arg3
GEN	arg4
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, arg3, arg4 ? &arg4 : NULL);
  }
 OUTPUT:
   RETVAL

GEN
interface32(arg1,arg2,arg3)
long	oldavma=avma;
GEN	arg1
GEN	arg2
long	arg3
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2,arg3);
  }
 OUTPUT:
   RETVAL

GEN
interface33(arg1,arg2,arg3,arg4=0)
long	oldavma=avma;
GEN	arg1
GEN	arg2
GEN	arg3
long	arg4
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2,arg3,arg4,prec);
  }
 OUTPUT:
   RETVAL

GEN
interface34(arg1,arg2,arg3)
long	oldavma=avma;
long	arg1
long	arg2
long	arg3
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, arg3);
  }
 OUTPUT:
   RETVAL

GEN
interface35(arg1,arg2,arg3)
long	oldavma=avma;
long	arg1
GEN	arg2
GEN	arg3
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1,arg2,arg3);
  }
 OUTPUT:
   RETVAL

GEN
interface37(arg1,arg2,arg3,arg4)
long	oldavma=avma;
PariVar	arg1
GEN	arg2
GEN	arg3
PariExpr	arg4
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, arg3, arg4, prec);
  }
 OUTPUT:
   RETVAL

GEN
interface47(arg1,arg2,arg3,arg4,arg0=gun)
long	oldavma=avma;
GEN	arg0
PariVar	arg1
GEN	arg2
GEN	arg3
PariExpr	arg4
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, arg3, arg4, arg0);
  }
 OUTPUT:
   RETVAL

GEN
interface48(arg1,arg2,arg3,arg4,arg0=gzero)
long	oldavma=avma;
GEN	arg0
PariVar	arg1
GEN	arg2
GEN	arg3
PariExpr	arg4
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, arg3, arg4, arg0);
  }
 OUTPUT:
   RETVAL

GEN
interface49(arg0,arg00,arg1,arg2,arg3)
long	oldavma=avma;
GEN	arg0
GEN	arg00
PariVar	arg1
PariVar	arg2
PariExpr	arg3
 CODE:
  {
    dFUNCTION(GEN);
# arg1 and arg2 may finish to be the same entree*, like after $x=$y=PARIvar 'x'
    if (arg1 == arg2) {
	if (ST(2) == ST(3)) 
	    croak("Same iterator for a double loop");
# ST(3) is localized now
	sv_unref(ST(3));
	arg2 = findVariable(ST(3),1);
	sv_setref_pv(ST(3), "Math::Pari::Ep", (void*)arg2);
    }
    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg0, arg00, arg1, arg2, arg3);
  }
 OUTPUT:
   RETVAL

GEN
interface83(arg1,arg2,arg3,arg4)
long	oldavma=avma;
PariVar	arg1
GEN	arg2
GEN	arg3
PariExpr	arg4
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, arg3, arg4);
  }
 OUTPUT:
   RETVAL

GEN
interface84(arg1,arg2,arg3)
long	oldavma=avma;
GEN	arg1
PariVar	arg2
PariExpr	arg3
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, arg3);
  }
 OUTPUT:
   RETVAL

# These interfaces were automatically generated:

long
interface16(arg1)
long	oldavma=avma;
    char * arg1
 CODE:
  {
    dFUNCTION(long);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1);
  }
 OUTPUT:
   RETVAL
 CLEANUP:
   avma=oldavma;


void
interface19(arg1, arg2)
long	oldavma=avma;
    long arg1
    long arg2
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    FUNCTION(arg1, arg2);
  }


GEN
interface44(arg1, arg2, arg3, arg4)
long	oldavma=avma;
    long arg1
    long arg2
    long arg3
    long arg4
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, arg3, arg4);
  }
 OUTPUT:
   RETVAL


GEN
interface45(arg1, arg2, arg3, arg4)
long	oldavma=avma;
    long arg1
    GEN arg2
    GEN arg3
    long arg4
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, arg3, arg4, prec);
  }
 OUTPUT:
   RETVAL


GEN
interface59(arg1, arg2, arg3, arg4, arg5)
long	oldavma=avma;
    long arg1
    GEN arg2
    GEN arg3
    GEN arg4
    GEN arg5
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, arg3, arg4, arg5);
  }
 OUTPUT:
   RETVAL


GEN
interface73(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
long	oldavma=avma;
    long arg1
    PariVar arg2
    GEN arg3
    GEN arg4
    PariExpr arg5
    long arg6
    long arg7
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION(arg1, arg2, arg3, arg4, arg5, prec, arg6, arg7);
  }
 OUTPUT:
   RETVAL


void
interface86(arg1, arg2, arg3, arg4, arg5)
long	oldavma=avma;
    PariVar arg1
    GEN arg2
    GEN arg3
    GEN arg4
    PariExpr arg5
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    FUNCTION(arg1, arg2, arg3, arg4, arg5);
  }


void
interface87(arg1, arg2, arg3, arg4=0)
long	oldavma=avma;
    PariVar arg1
    GEN arg2
    PariExpr arg3
    long arg4
 CODE:
  {
    dFUNCTION(GEN);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    FUNCTION(arg1, arg2, arg3, arg4);
  }


bool
_2bool(arg1,arg2,inv)
long	oldavma=avma;
     GEN	arg1
     GEN	arg2 = NO_INIT
     long	inv = NO_INIT
   CODE:
     RETVAL=!gcmp0(arg1);
   OUTPUT:
     RETVAL

bool
pari2bool(arg1)
long	oldavma=avma;
     GEN	arg1
   CODE:
     RETVAL=!gcmp0(arg1);
   OUTPUT:
     RETVAL

CV *
loadPari(name)
     char *	name
   CODE:
     {
       char *olds = name;
       entree *ep=NULL;
       long hash, valence;
       void (*func)(void*)=NULL;
       void (*unsupported)(void*) = (void (*)(void*)) not_here;

       if (*name=='g') {
	   switch (name[1]) {
	   case 'a':
	       if (strEQ(name,"gadd")) {
		   valence=2;
		   func=(void (*)(void*)) gadd;
	       } else if (strEQ(name,"gand")) {
		   valence=2;
		   func=(void (*)(void*)) gand;
	       }
	       break;
	   case 'c':
	       if (strEQ(name,"gcmp0")) {
		   valence=10;
		   func=(void (*)(void*)) gcmp0;
	       } else if (strEQ(name,"gcmp1")) {
		   valence=10;
		   func=(void (*)(void*)) gcmp1;
	       } else if (strEQ(name,"gcmp_1")) {
		   valence=10;
		   func=(void (*)(void*)) gcmp_1;
	       } else if (strEQ(name,"gcmp")) {
		   valence=20;
		   func=(void (*)(void*)) gcmp;
	       }
	       break;
	   case 'd':
	       if (strEQ(name,"gdiv")) {
		   valence=2;
		   func=(void (*)(void*)) gdiv;
	       } else if (strEQ(name,"gdivent")) {
		   valence=2;
		   func=(void (*)(void*)) gdivent;
	       } else if (strEQ(name,"gdivround")) {
		   valence=2;
		   func=(void (*)(void*)) gdivround;
	       }
	       break;
	   case 'e':
	       if (strEQ(name,"geq")) {
		   valence=2;
		   func=(void (*)(void*)) geq;
	       } else if (strEQ(name,"gegal")) {
		   valence=20;
		   func=(void (*)(void*)) gegal;
	       }
	       break;
	   case 'g':
	       if (strEQ(name,"gge")) {
		   valence=2;
		   func=(void (*)(void*)) gge;
	       } else if (strEQ(name,"ggt")) {
		   valence=2;
		   func=(void (*)(void*)) ggt;
	       } 
	       break;
	   case 'l':
	       if (strEQ(name,"gle")) {
		   valence=2;
		   func=(void (*)(void*)) gle;
	       } else if (strEQ(name,"glt")) {
		   valence=2;
		   func=(void (*)(void*)) glt;
	       } 
	       break;
	   case 'm':
	       if (strEQ(name,"gmul")) {
		   valence=2;
		   func=(void (*)(void*)) gmul;
	       } else if (strEQ(name,"gmod")) {
		   valence=2;
		   func=(void (*)(void*)) gmod;
	       } 
	       break;
	   case 'n':
	       if (strEQ(name,"gneg")) {
		   valence=1;
		   func=(void (*)(void*)) gneg;
	       } else if (strEQ(name,"gne")) {
		   valence=2;
		   func=(void (*)(void*)) gne;
	       } 
	       break;
	   case 'o':
	       if (strEQ(name,"gor")) {
		   valence=2;
		   func=(void (*)(void*)) gor;
	       }
	       break;
	   case 'p':
	       if (strEQ(name,"gpui")) {
		   valence=2;
		   func=(void (*)(void*)) gpui;
	       }
	       break;
	   case 's':
	       if (strEQ(name,"gsub")) {
		   valence=2;
		   func=(void (*)(void*)) gsub;
	       }
	       break;
	   }
       } else if (*name=='_') {
	   if (name[1] == 'g') {
	       switch (name[2]) {
	       case 'a':
		   if (strEQ(name,"_gadd")) {
		       valence=299;
		       func=(void (*)(void*)) gadd;
		   } else if (strEQ(name,"_gand")) {
		       valence=2099;
		       func=(void (*)(void*)) gand;
		   } 
		   break;
	       case 'c':
		   if (strEQ(name,"_gcmp")) {
		       valence=209;
		       func=(void (*)(void*)) gcmp;
		   } else if (strEQ(name,"_gcmp0")) {
		       valence=109;
		       func=(void (*)(void*)) gcmp0;
		   }
		   break;
	       case 'd':
		   if (strEQ(name,"_gdiv")) {
		       valence=299;
		       func=(void (*)(void*)) gdiv;
		   }
		   break;
	       case 'e':
		   if (strEQ(name,"_geq")) {
		       valence=2099;
		       func=(void (*)(void*)) geq;
		   }
		   break;
	       case 'g':
		   if (strEQ(name,"_gge")) {
		       valence=2099;
		       func=(void (*)(void*)) gge;
		   } else if (strEQ(name,"_ggt")) {
		       valence=2099;
		       func=(void (*)(void*)) ggt;
		   }
		   break;
	       case 'l':
		   if (strEQ(name,"_gle")) {
		       valence=2099;
		       func=(void (*)(void*)) gle;
		   } else if (strEQ(name,"_glt")) {
		       valence=2099;
		       func=(void (*)(void*)) glt;
		   }
		   break;
	       case 'm':
		   if (strEQ(name,"_gmul")) {
		       valence=299;
		       func=(void (*)(void*)) gmul;
		   } else if (strEQ(name,"_gmod")) {
		       valence=299;
		       func=(void (*)(void*)) gmod;
		   }
		   break;
	       case 'n':
		   if (strEQ(name,"_gneg")) {
		       valence=199;
		       func=(void (*)(void*)) gneg;
		   } else if (strEQ(name,"_gne")) {
		       valence=2099;
		       func=(void (*)(void*)) gne;
		   }
		   break;
	       case 'o':
		   if (strEQ(name,"_gor")) {
		       valence=2099;
		       func=(void (*)(void*)) gor;
		   }
		   break;
	       case 'p':
		   if (strEQ(name,"_gpui")) {
		       valence=299;
		       func=(void (*)(void*)) gpui;
		   }
		   break;
	       case 's':
		   if (strEQ(name,"_gsub")) {
		       valence=299;
		       func=(void (*)(void*)) gsub;
		   } 
		   break;
	       }
	   } else {
	       switch (name[1]) {
	       case 'a':
		   if (strEQ(name,"_abs")) {
		       valence=199;
		       func=(void (*)(void*)) gabs;
		   } 
		   break;
	       case 'c':
		   if (strEQ(name,"_cos")) {
		       valence=199;
		       func=(void (*)(void*)) gcos;
		   } 
		   break;
	       case 'e':
		   if (strEQ(name,"_exp")) {
		       valence=199;
		       func=(void (*)(void*)) gexp;
		   } 
		   break;
	       case 'l':
		   if (strEQ(name,"_lex")) {
		       valence=209;
		       func=(void (*)(void*)) lexcmp;
		   } else if (strEQ(name,"_log")) {
		       valence=199;
		       func=(void (*)(void*)) glog;
		   } 
		   break;
	       case 's':
		   if (strEQ(name,"_sin")) {
		       valence=199;
		       func=(void (*)(void*)) gsin;
		   } else if (strEQ(name,"_sqrt")) {
		       valence=199;
		       func=(void (*)(void*)) gsqrt;
		   }
		   break;
	       }
	   }
       }
       if (!func) {
	   SAVEINT(doing_PARI_autoload);
	   doing_PARI_autoload = 1;
	   ep = is_entry_intern(name, functions_hash, &hash);
	   doing_PARI_autoload = 0;
#if 0
	 for (n = 0; *name; name++) n = n << 1 ^ *name;
	 if (n < 0) n = -n; n %= TBLSZ;
	 for(ep = hashtable[n]; ep; ep = ep->next) {
	   if (strEQ(olds,ep->name)) { /* Name in the symbol table */
	     break;
	   }
	 }
#endif
	 if (!ep) {
#if 0					/* findentry() is static. */
	     ep = findentry(name,strlen(name),funct_old_hash[hash]);
#endif
	     if (!ep)
		 croak("`%s' is not a Pari function name",name);
	     else
		 warn("`%s' is an obsolete Pari function name", name);
	 }
	 if (ep && (EpVALENCE(ep) < EpUSER 
		    /* && ep>=fonctions && ep < fonctions+NUMFUNC) */)) {
	     /* Builtin */
	   valence=EpVALENCE(ep);
	   func=(void (*)(void*)) ep->value;
	   if (!func) {
	     func = unsupported;
	   }
	 }
       }
       if (func == unsupported) {
	 croak("Do not know how to work with Pari control structure `%s'",
	       olds);
       } else if (func) {
	 void (*subaddr)(CV*);
	 char* file = __FILE__, *proto = NULL;
	 char subname[276]="Math::Pari::";
	 char buf[10];
	 CV *protocv;
	 
	 strcpy(subname+12,"interface");
	 sprintf(buf, "%d", valence);
	 strcpy(subname+12+9,buf);
	 protocv = perl_get_cv(subname, FALSE);
	 if (protocv) {
	     proto = SvPV((SV*)protocv,na);
	 }
	 
	 strcpy(subname+12,olds);
	 switch (valence) {
	 case 0:
	     if (ep->code[0] == 'p' && ep->code[1] == 0) {
		 DO_INTERFACE(0);
	     } else if (ep->code[0] == 0) {
		 DO_INTERFACE(00);
	     } else {
		 croak("Unsupported interface %d for a Pari function %s with code \"%s\"",
		       valence, olds, ep->code);
	     }
	     break;
	   CASE_INTERFACE(1);
	   CASE_INTERFACE(10);
	   CASE_INTERFACE(199);
	   CASE_INTERFACE(109);
	   CASE_INTERFACE(11);
	   CASE_INTERFACE(15);
	   CASE_INTERFACE(2);
	   CASE_INTERFACE(20);
	   CASE_INTERFACE(299);
	   CASE_INTERFACE(209);
	   CASE_INTERFACE(2099);
	   CASE_INTERFACE(3);
	   CASE_INTERFACE(30);
	   CASE_INTERFACE(4);
	   CASE_INTERFACE(5);
	   CASE_INTERFACE(21);
	   CASE_INTERFACE(23);
	   CASE_INTERFACE(24);
	   CASE_INTERFACE(25);
	   CASE_INTERFACE(29);
	   CASE_INTERFACE(32);
	   CASE_INTERFACE(33);
	   CASE_INTERFACE(35);
	   CASE_INTERFACE(12);
	   CASE_INTERFACE(13);
	   CASE_INTERFACE(14);
	   CASE_INTERFACE(26);
	   CASE_INTERFACE(28);
	   CASE_INTERFACE(31);
	   CASE_INTERFACE(34);
	   CASE_INTERFACE(22);
	   CASE_INTERFACE(27);
	   CASE_INTERFACE(37);
	   CASE_INTERFACE(47);
	   CASE_INTERFACE(48);
	   CASE_INTERFACE(49);
	   CASE_INTERFACE(83);
	   CASE_INTERFACE(84);
	   CASE_INTERFACE(18);
	   /* These interfaces were automatically generated: */
	   CASE_INTERFACE(16);
	   CASE_INTERFACE(19);
	   CASE_INTERFACE(44);
	   CASE_INTERFACE(45);
	   CASE_INTERFACE(59);
	   CASE_INTERFACE(73);
	   CASE_INTERFACE(86);
	   CASE_INTERFACE(87);

	 default: croak("Unsupported interface %d for a Pari function %s",
			valence, olds);
	 }
	 RETVAL = newXS(subname,subaddr,file);
	 if (proto)
	     sv_setpv((SV*)RETVAL, proto);
	 CvXSUBANY(RETVAL).any_dptr = func;
       } else {
	 croak("Cannot load a Pari macro `%s'", olds);
       }
     }
   OUTPUT:
     RETVAL


# Tag is menu entry, or -1 for all.

SV *
listPari(tag)
     int tag
   PPCODE:
     {
       long v, valence;
       entree *ep;

       for(ep = functions_basic; ep->name; ep++)  {
	   valence = EpVALENCE(ep);
	   if (tag == -1 || ep->menu == tag) {
	       switch (valence) {
		   case 0:
		       if ((ep->code[0] != 0) 
			   && ((ep->code[0] != 'p' || ep->code[1] != 0)))
			   break;
		       /* FALL THROUGH */
	           case 1:
		   case 10:
		   case 199:
		   case 109:
		   case 11:
		   case 15:
		   case 2:
		   case 20:
		   case 299:
		   case 209:
		   case 2099:
		   case 3:
		   case 30:
		   case 4:
		   case 5:
		   case 21:
		   case 23:
		   case 24:
		   case 25:
		   case 29:
		   case 32:
		   case 33:
		   case 35:
		   case 12:
		   case 13:
		   case 14:
		   case 26:
		   case 28:
		   case 31:
		   case 34:
		   case 22:
		   case 27:
		   case 37:
		   case 47:
		   case 48:
		   case 49:
		   case 83:
		   case 84:
		   case 18:
		       /* These interfaces were automatically generated: */
	           case 16:
		   case 19:
		   case 44:
		   case 45:
		   case 59:
		   case 73:
		   case 86:
		   case 87:
		   XPUSHs(sv_2mortal(newSVpv(ep->name, 0)));
	       }
	   }
       }
     }

BOOT:
{
   SV *mem = perl_get_sv("Math::Pari::initmem", FALSE);
   SV *pri = perl_get_sv("Math::Pari::initprimes", FALSE);
   if (!mem || !SvOK(mem)) {
       croak("$Math::Pari::initmem not defined!");
   }
   if (!pri || !SvOK(pri)) {
       croak("$Math::Pari::initprimes not defined!");
   }
   INIT_JMP_off;
   INIT_SIG_off;
   /* These guys are new in 2.0. */
   init_defaults(1);
   pari_addfunctions(&pari_modules, functions_highlevel,helpmessages_highlevel);
   init_graph();

   init(SvIV(mem),SvIV(pri)); /* Default: take four million bytes of
			       * memory for the stack, calculate
			       * primes up to 500000. */
   PariStack = (SV *) GENfirstOnStack;
   workErrsv = newSVpv("",0);
   pariErr = &perlErr;
   foreignHandler = (void*)&callPerlFunction;
   foreignAutoload = &autoloadPerlFunction;
   foreignExprSwitch = (char)SVt_PVCV;
   foreignExprHandler = &exprHandler_Perl;
   foreignFuncFree = &freePerlFunction;
   pariStash = gv_stashpv("Math::Pari", TRUE);
   pariEpStash = gv_stashpv("Math::Pari::Ep", TRUE);
   perlavma = sentinel = avma;
}

void
memUsage()
PPCODE:
#ifdef DEBUG_PARI
  EXTEND(sp, 3);		/* Got cv + 0, return 4. */
  PUSHs(sv_2mortal(newSViv(SVnumtotal)));
  PUSHs(sv_2mortal(newSViv(SVnum)));
  PUSHs(sv_2mortal(newSViv(onStack)));
  PUSHs(sv_2mortal(newSViv(offStack)));
#endif  
  

MODULE = Math::Pari PACKAGE = Math::Pari

void
DESTROY(rv)
     SV *	rv
   CODE:
     {
	 /* PariStack keeps the latest SV that keeps a GEN on stack. */
	 SV* sv = SvRV(rv);
	 char* type = SvPVX(sv);	/* The value of PariStack when the
					 * variable was created, thus the
					 * previous SV that keeps a GEN from
					 * stack, or some atoms. */
	 long oldavma = SvCUR(sv) + bot; /* The value of avma on the entry
					  * to function having the SV as
					  * argument. */
	 long howmany;
       
	 SvPVX(sv) = GENheap;		/* To avoid extra free() in moveoff.... */
	 switch ((IV)type) {
	 case (IV)GENheap:	/* Leave it alone? XXXX */
	     break;
	 case (IV)GENmovedOffStack:	/* Know that it _was temporary. */
	     killbloc((GEN)SvIV(sv));
	     break;
	 default:		/* Still on stack */
	     if (type != (char*)PariStack) { /* But not the newest one. */
		 howmany=moveoffstack_newer_than(sv);
		 DEBUG_u( deb("%li items moved off stack\n", howmany) );
	     }
	     /* Now fall through: */
/* case (IV)GENfirstOnStack: */
	     /* Now sv is the newest one on stack. */
	     onStack_dec;
	     perlavma = oldavma;
	     if (oldavma > sentinel) {
		 avma = sentinel;	/* Mark the space on stack as free. */
	     } else {
		 avma = oldavma;	/* Mark the space on stack as free. */
	     }
	     PariStack = (SV*)type; /* The same on the Perl/PARI side. */
	 }
	 SVnum_dec;
     }


SV *
pari_print(in)
    GEN in

SV *
pari_pprint(in)
    GEN in

SV *
pari_texprint(in)
    GEN in

I32
typ(in)
    GEN in

SV *
PARIvar(in)
    char *in

GEN
ifact(arg1)
long	oldavma=avma;
long	arg1

void
changevalue(name, val)
    PariName name
    GEN val

void
set_gnuterm(a,b)
    IV a
    IV b
