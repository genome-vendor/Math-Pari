#  include <pari.h>
#  include <language/anal.h>
#  include <gp/gp.h>			/* init_opts */

/* On some systems /usr/include/sys/dl.h attempts to declare
   ladd which pari.h already defined with a different meaning.

   It is not clear whether this is a correct fix...
 */
#undef ladd

#define PERL_POLLUTE			/* We use older varnames */

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

#ifdef warner
#  undef warner
#endif

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

#ifndef XSINTERFACE_FUNC_SET		/* Not in 5.004_04 */
#  define XSINTERFACE_FUNC_SET(cv,f)	\
		CvXSUBANY(cv).any_dptr = (void (*) (void*))(f)
#endif

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

#if defined(MYMALLOC) && defined(EMBEDMYMALLOC) && defined(UNEMBEDMYMALLOC)

Malloc_t
malloc(register size_t nbytes)
{
    return Perl_malloc(nbytes);
}

Free_t
free(void *mp)
{
    Perl_mfree(mp);			/* What to return? */
}

Malloc_t
realloc(void *mp, size_t nbytes)
{
    return Perl_realloc(mp, nbytes);
}

#endif

/* We make a "fake" PVAV, not enough entries.  */

#define setSVpari(sv, in, oldavma) do {				\
    sv_setref_pv(sv, "Math::Pari", (void*)in);			\
    if (is_matvec_t(typ(in)) && SvTYPE(SvRV(sv)) != SVt_PVAV) {	\
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

#define pari_version_exp() PARI_VERSION_EXP

#if PARI_VERSION_EXP >= 2000018

GEN
_gbitneg(GEN g)
{
    return gbitneg(g,-1);
}

#endif	/* PARI_VERSION_EXP >= 2000018 */ 

void
make_PariAV(SV *sv)
{
    AV *av = (AV*)SvRV(sv);
    char *s = SvPVX(av);
    IV i = SvIVX(av);
    MAGIC *mg;
    SV *new = newRV_noinc((SV*)av);	/* cannot use sv, it may be 
					   sv_restore()d */

    SvUPGRADE((SV*)av, SVt_PVAV);    
    SvPVX(av)	= s;
    SvIVX(av)	= i;
    sv_magic((SV*)av, new, 'P', Nullch, 0);
    SvREFCNT_dec(new);			/* now RC(new)==1 */
	/* We avoid an reference loop, so should be careful on DESTROY */
#if 0
    if ((mg = SvMAGIC((SV*)av)) && mg->mg_type == 'P' /* be extra paranoid */
	&& (mg->mg_flags & MGf_REFCOUNTED)) {
/*      mg->mg_flags &= ~MGf_REFCOUNTED;	*/
/*	SvREFCNT_dec(sv);			*/
	sv_2mortal((SV*)av);		/* We restore refcount on DESTROY */
    }
#endif
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

#if PARI_VERSION_EXP <= 2002000		/* Global after 2.2.0 */
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
#endif

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
      ep = (entree *)gpmalloc(sizeof(entree) + 7*BYTES_IN_LONG 
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
	      if (typ(x) == t_POL	/* Polynomial. */
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
  ep = (entree *)gpmalloc(sizeof(entree) + 7*BYTES_IN_LONG 
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

unsigned long
longword(GEN x, long n)
{
  if (n < 0 || n >= lg(x))
      croak("The longword %ld ordinal out of bound", n);
  return x[n];
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
  if (SvGMAGICAL(sv)) mg_get(sv); /* MAYCHANGE in perlguts.pod - bug in perl */
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
	      GEN ret=cgetg(len+2, t_VEC);
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
  else if (SvNOKp(sv)) {
      /* Need more voodoo, since sv_true sv_false are NOK, but not IOK */
      if (/*SvIOKp(sv) && */ SvPOKp(sv) && SvCUR(sv)) {
	  char *s = SvPVX(sv), *e = SvEND(sv);

	  if (*s == '+' || *s == '-')
	      s++;
	  while (1) {
	      if (s == e)
		  return stoi(SvIV(sv)); /* In fact, is integer */
	      if (*s < '0' || *s > '9')	/* Not an integer */
		  break;
	      s++;
	  }
      }
      return dbltor(SvNV(sv));
  } else if (SvPOK(sv)) return lisexpr(SvPV(sv,na));
  else if (!SvOK(sv)) return stoi(0);
  else
    croak("Variable in perl2pari is not of known type");  
}

GEN
sv2parimat(SV* sv)
{
  GEN in=sv2pari(sv);
  if (typ(in)==t_VEC) {
    long len=lg(in)-1;
    long t;
    long l=lg(in[1]);
    for (;len;len--) {
      if ((t=typ(in[len])) == t_VEC) {
	settyp(in[len], t_COL);
      } else if (t != t_COL) {
	croak("Not a vector where column of a matrix expected");
      }
      if (lg(in[len])!=l) {
	croak("Columns of input matrix are of different height");
      }
    }
    settyp(in, t_MAT);
    return in;
  } else if (typ(in) == t_MAT) {
    return in;
  } else {
    croak("Not a matrix where matrix expected");
  }
}

SV*
pari2iv(GEN in)
{
#ifdef SvIsUV
#  define HAVE_UVs 1
    UV uv;
#else
#  define HAVE_UVs 0
    IV uv;
#endif
    int overflow = 0;

    if (typ(in) != t_INT) 
	return newSViv((IV)gtolong(in));
    switch (lgef(in)) {
    case 2:
	uv = 0;
	break;
    case 3:
	uv = in[2];
	if (sizeof(long) >= sizeof(IV) && in[2] < 0)
	    overflow = 1;
	break;
    case 4:
	if ( 2 * sizeof(long) > sizeof(IV)
	     || (2 * sizeof(long) == sizeof(IV)) && !HAVE_UVs && in[2] < 0 )
	    goto do_nv;
	uv = in[2];
	uv = (uv << TWOPOTBYTES_IN_LONG) + in[3];
	break;
    default:
	goto do_nv;
    }
    if (overflow) {
#ifdef SvIsUV
	if (signe(in) > 0) {
	    SV *sv = newSViv((IV)uv);

	    SvIsUV_on(sv);
	    return sv;
	} else
#endif
	    goto do_nv;
    }
    return newSViv(signe(in) > 0 ? (IV)uv : -(IV)uv);
do_nv:
    return newSVnv(gtodouble(in));	/* XXXX to NV, not to double? */
}

SV*
pari2nv(GEN in)
{
  return newSVnv(gtodouble(in));
}

SV*
pari2pv(GEN in)
{
    if (typ(in) == t_STR)		/* Puts "" around without special-casing */
	return newSVpv(GSTR(in),0);
    {    
	PariOUT *oldOut = pariOut;
	pariOut = &perlOut;
	worksv = newSVpv("",0);
	bruteall(in,'g',-1,0);		/* 0: compact pari-readable form */
	pariOut = oldOut;
	return worksv;
    }
}

int fmt_nb;

#ifdef LONG_IS_64BIT
#  define def_fmt_nb 38
#else
#  define def_fmt_nb 28
#endif

long
setprecision(long digits)
{
  long m = fmt_nb;

  if(digits>0) {fmt_nb = digits; prec = (long)(digits*pariK1 + 3);}
  return m;
}

long
setseriesprecision(long digits)
{
  long m = precdl;

  if(digits>0) {precdl = digits;}
  return m;
}

static IV primelimit;
static UV parisize;

IV
setprimelimit(IV n)
{
    byteptr ptr;
    IV o = primelimit;

    if (n != 0) {
	ptr = initprimes(n);
	free(diffptr);
	diffptr = ptr;
	primelimit = n;
    }
    return primelimit;
}

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
	code = (char *)gpmalloc(numargs*6 - req*5 + 2);
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

void
detach_stack(void)
{
    moveoffstack_newer_than((SV *) GENfirstOnStack);
}

UV
allocatemem(UV newsize)
{
    if (newsize) {
	detach_stack();
	parisize = allocatemoremem(newsize);
    }
    return parisize;
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

    if (!is_matvec_t(typ(g)))
	croak("Access to elements of not-a-vector");
    if (n >= l || n < 0)
	croak("Array index %i out of range", n);
#if 0
    warn("fetching %d-th element of type %d", n, typ((GEN)g[n + 1]));
#endif
    return (GEN)g[n + 1];
}

static void
Arr_STORE(GEN g, I32 n, GEN elt)
{
    I32 l = lg(g) - 1, docol = 0;
    GEN old;

    if (!is_matvec_t(typ(g)))
	croak("Access to elements of not-a-vector");
    if (n >= l || n < 0)
	croak("Array index %i out of range", n);
#if 0
    warn("storing %d-th element of type %d", n, typ((GEN)g[n + 1]));
#endif	/* 0 */

    if (typ(g) == t_MAT) {
	long len = lg(g);
	long l   = lg(g[1]);
	if (typ(elt) != t_COL) {
	    if (typ(elt) != t_VEC)
		croak("Not a vector where column of a matrix expected");
	    docol = 1;
	}
	if (lg(elt)!=l && len != 2)
	    croak("Assignment of a columns into a matrix of incompatible height");
    }
    
    old = (GEN)g[n + 1];
    /* It is not clear whether we need to clone if the elt is offstack */
    elt = gclone(elt);
    if (docol)
	settyp(elt, t_COL);

    /* anal.c is optimizing inspection away around here... */
    if (isclone(old)) killbloc0(old, 1);
    g[n + 1] = (long)elt;
}

#define Arr_FETCHSIZE(g)  (lg(g) - 1)

#define DFT_VAR (GEN)-1
#define DFT_GEN (GEN)NULL

static void
check_pointer(unsigned int ptrs, GEN argvec[])
{
  unsigned int i;
  for (i=0; ptrs; i++,ptrs>>=1)
    if (ptrs & 1) *((GEN*)argvec[i]) = gclone(*((GEN*)argvec[i]));
}

#define RETTYPE_VOID 0
#define RETTYPE_LONG 1
#define RETTYPE_GEN 2

static void
fill_argvect(entree *ep, char *s, long *has_pointer, GEN *argvec,
	     long *rettype, SV **args, int items)
{
    entree *ep1;
    int i = 0, j = 0;
    long fake;
    PariExpr expr;

    if (!ep)
	croak("XSUB call through interface did not provide *function");
    if (!s)
	croak("XSUB call through interface with a NULL code");

    while (*s) {
	if (i >= 8)
	    croak("Too many args for a flexible-interface function");
	switch (*s++)
	    {
	    case 'G': /* GEN */
		argvec[i++] = sv2pari(args[j++]);
		break;

	    case 'L': /* long */
		argvec[i++] = (GEN) (long)SvIV(args[j]);
		j++;
		break;

	    case 'n': /* var number */
		argvec[i++] = (GEN) numvar(sv2pari(args[j++]));
		break;

	    case 'V': /* variable */
		ep1 = bindVariable(args[j++]);
		argvec[i++] = (GEN)ep1;
		if (EpVALENCE(ep1) != EpVAR && *(s-1) == 'V')
		    croak("Did not get a variable");
		break;
	    case 'S': /* symbol */
		ep1 = bindVariable(args[j++]);
		argvec[i++] = (GEN)ep1;
		break;
	    case '&': /* *GEN */
		ep1 = bindVariable(args[j++]);
		if (ep1->value == (void*)initial_value(ep1))
		    changevalue(ep1, gzero); /* don't overwrite initial value */
		*has_pointer |= (1 << i);
		argvec[i++] = (GEN) &(ep1->value); 
		break;
	    case  'E': /* Input position - subroutine */
	    case  'I': /* Input position - subroutine */
		if (SvROK(args[j]) && SvTYPE(SvRV(args[j])) == SVt_PVCV) {
		    expr = ((char*)&(SvRV(args[j])->sv_flags)) + LSB_in_U32;
		} else expr = (char *)SvPV(args[j],na);
		argvec[i++] = (GEN) expr;
		j++;
		break;

	    case 'r': /* raw */
	    case 's': /* expanded string; empty arg yields "" */
		argvec[i++] = (GEN) SvPV(args[j],na);
		j++;
		break;

	    case 'p': /* precision */
		argvec[i++] = (GEN) prec; 
		break;

	    case '=':
	    case ',':
		break;

	    case 'D': /* Has a default value */
		if (j >= items || !SvOK(args[j]))
		    {
			char *pre = s;

			if (j < items)
			    j++;

			if ( *s == 'G' || *s == '&' 
			     || *s == 'E' || *s == 'I' || *s == 'V') { 
			    argvec[i++]=DFT_GEN; s++; 
			    break; 
			}
			if (*s == 'n')              { 
			    argvec[i++]=DFT_VAR; s++; 
			    break; 
			}
			while (*s++ != ',');
			switch (*s) {
			case 'r': case 's':
			    if (pre[0] == '\"' && pre[1] == '\"' 
				&& s - pre == 3) {
				argvec[i++] = (GEN) "";
				break;
			    }
			    goto unknown;
			case 'L': /* long */
			    argvec[i++] = (GEN) atol(pre);
			    break;
			case 'G':
			    if ((*pre == '1' || *pre == '0') && pre[1]==',') {
				argvec[i++] = (*pre == '1'
					       ? gun : gzero);
				break;
			    }
			default:
			  unknown:
			    croak("Cannot process default argument %.*s of type %.1s",
				  s - pre - 1, pre, s);
			}
			s++;			/* Skip ',' */
		    }
		else
		    if (*s == 'G' || *s == '&' || *s == 'n'
			|| *s == 'E' || *s == 'I' || *s == 'V') 
			break;
		while (*s++ != ',');
		break;

	    case 'P': /* series precision */
		argvec[i++] = (GEN) precdl; 
		break;

	    case 'f': /* Fake *long argument */
		argvec[i++] = (GEN) &fake; 
		break;

	    case 'x': /* Foreign function */
		croak("Calling Perl via PARI with an unknown interface: avoiding loop");
		break;

	    case 'l': /* Return long */
		*rettype = RETTYPE_LONG; break;

	    case 'v': /* Return void */
		*rettype = RETTYPE_VOID; break;

	    default: 
		croak("Unsupported code '%.1s' in signature of a PARI function", s-1);
	    }
	if (j > items)
	    croak("Too few args %d for PARI function %s", items, ep->name);
    }
    if (j < items)
	croak("%d unused args for PARI function %s", items - j, ep->name);
#if PURIFY
    for ( ; i<9; i++) argvec[i]=NULL; 
#endif
}

#define _to_int(in,dummy1,dummy2) to_int(in)

static GEN
to_int(GEN in)
{
    long sign = gcmp(in,gzero);
    long dummy;
    

    if (!sign)
	return gzero;
    switch (typ(in)) {
    case t_INT:
    case t_SMALL:
	return in;
    case t_INTMOD:
	return lift0(in, -1);		/* -1: not as polmod */
    default:
	return gtrunc(in);
    }
}

typedef int (*FUNC_PTR)();
typedef void (*TSET_FP)(char *s);
#if PARI_VERSION_EXP < 2000013
#  define set_gnuterm(a,b,c) \
	set_term_funcp((FUNC_PTR)(a),(struct termentry *)(b))
#else	/* !( PARI_VERSION_EXP < 2000013 ) */ 
#  define set_gnuterm(a,b,c) \
	set_term_funcp3((FUNC_PTR)(a),(struct termentry *)(b), (TSET_FP)(c))
#endif	/* PARI_VERSION_EXP < 2000013 */

#define int_set_term_ftable(a) (v_set_term_ftable((void*)a))
extern  void v_set_term_ftable(void *a);

#define s_type_name(x) type_name(typ(x));


MODULE = Math::Pari PACKAGE = Math::Pari PREFIX = Arr_

PROTOTYPES: ENABLE

GEN
Arr_FETCH(g,n)
    long	oldavma=avma;
    GEN g
    I32 n

void
Arr_STORE(g,n,elt)
    long	oldavma=avma;
    GEN g
    I32 n
    GEN elt

I32
Arr_FETCHSIZE(g)
    GEN g

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
     if (typ(in) == t_INT) {
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
_to_int(in, dummy1, dummy2)
long	oldavma=avma;
    GEN in
    SV  *dummy1
    SV  *dummy2

GEN
PARI(...)
long	oldavma=avma;
   CODE:
     if (items==1) {
       RETVAL=sv2pari(ST(0));
     } else {
	int i;

	RETVAL=cgetg(items+1, t_VEC);
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

	RETVAL=cgetg(items+1, t_VEC);
	for (i=0;i<items;i++) {
	  RETVAL[i+1]=(long)sv2pari(ST(i));
	}
     }
     settyp(RETVAL, t_COL);
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

	RETVAL=cgetg(items+1, t_VEC);
	for (i=0;i<items;i++) {
	  RETVAL[i+1]=(long)sv2pari(ST(i));
	  settyp(RETVAL[i+1], t_COL);
	}
     }
     settyp(RETVAL, t_MAT);
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

void
interface_flexible_void(...)
long	oldavma=avma;
 CODE:
  {
    entree *ep = (entree *) XSANY.any_dptr, *ep1;
    void (*FUNCTION_real)(VARARG)
	= (void (*)(VARARG))ep->value;
    GEN argvec[9];
    long rettype = RETTYPE_GEN;
    long has_pointer = 0;		/* XXXX ?? */

    fill_argvect(ep, ep->code, &has_pointer, argvec, &rettype, &ST(0), items);

    if (rettype != RETTYPE_VOID)
	croak("Expected VOID return type, got code '%s'", ep->code);
    
    (FUNCTION_real)(argvec[0], argvec[1], argvec[2], argvec[3],
	          argvec[4], argvec[5], argvec[6], argvec[7], argvec[8]);
    if (has_pointer) 
	check_pointer(has_pointer,argvec);
  }

GEN
interface_flexible_gen(...)
long	oldavma=avma;
 CODE:
  {
    entree *ep = (entree *) XSANY.any_dptr, *ep1;
    GEN (*FUNCTION_real)(VARARG)
	= (GEN (*)(VARARG))ep->value;
    GEN argvec[9];
    long rettype = RETTYPE_GEN;
    long has_pointer = 0;		/* XXXX ?? */

    fill_argvect(ep, ep->code, &has_pointer, argvec, &rettype, &ST(0), items);

    if (rettype != RETTYPE_GEN)
	croak("Expected GEN return type, got code '%s'", ep->code);
    
    RETVAL=(FUNCTION_real)(argvec[0], argvec[1], argvec[2], argvec[3],
	          argvec[4], argvec[5], argvec[6], argvec[7], argvec[8]);
    if (has_pointer) 
	check_pointer(has_pointer,argvec);
  }
 OUTPUT:
   RETVAL

long
interface_flexible_long(...)
long	oldavma=avma;
 CODE:
  {
    entree *ep = (entree *) XSANY.any_dptr, *ep1;
    long (*FUNCTION_real)(VARARG)
	= (long (*)(VARARG))ep->value;
    GEN argvec[9];
    long rettype = RETTYPE_GEN;
    long has_pointer = 0;		/* XXXX ?? */

    fill_argvect(ep, ep->code, &has_pointer, argvec, &rettype, &ST(0), items);

    if (rettype != RETTYPE_LONG)
	croak("Expected long return type, got code '%s'", ep->code);
    
    RETVAL=FUNCTION_real(argvec[0], argvec[1], argvec[2], argvec[3],
	          argvec[4], argvec[5], argvec[6], argvec[7], argvec[8]);
    if (has_pointer) 
	check_pointer(has_pointer,argvec);
  }
 OUTPUT:
   RETVAL

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

    RETVAL = FUNCTION(arg1, arg2, arg3);
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
interface28(arg1,arg2=0,arg3=0)
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

    RETVAL = FUNCTION(arg1, arg2, arg3);
  }
 OUTPUT:
   RETVAL

GEN
interface28_old(arg1,arg2)
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

void
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

    FUNCTION(arg1, arg2, arg3);
  }

void
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

    FUNCTION(arg1,arg2,arg3);
  }

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
interface49(arg0,arg00,arg1=0,arg2=0,arg3=0)
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
    if (arg1 == arg2 && arg1) {
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

void
interface83(arg1,arg2,arg3,arg4)
long	oldavma=avma;
PariVar	arg1
GEN	arg2
GEN	arg3
PariExpr	arg4
 CODE:
  {
    dFUNCTION(void);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    FUNCTION(arg1, arg2, arg3, arg4);
  }

void
interface84(arg1,arg2,arg3)
long	oldavma=avma;
GEN	arg1
PariVar	arg2
PariExpr	arg3
 CODE:
  {
    dFUNCTION(void);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    FUNCTION(arg1, arg2, arg3);
  }


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
interface45(arg1, arg2, arg3=0)
long	oldavma=avma;
    long arg1
    GEN arg2
    long arg3
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


void
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

    FUNCTION(arg1, arg2, arg3, arg4, arg5);
  }


GEN
interface73(arg1, arg2, arg3, arg4, arg5, arg6=0, arg7=0)
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
#if PARI_VERSION_EXP >= 2000018
	       case 'b':
		   if (strEQ(name,"_gbitand")) {
		       valence=299;
		       func=(void (*)(void*)) gbitand;
		   } else if (strEQ(name,"_gbitor")) {
		       valence=299;
		       func=(void (*)(void*)) gbitor;
		   } else if (strEQ(name,"_gbitxor")) {
		       valence=299;
		       func=(void (*)(void*)) gbitxor;
		   } else if (strEQ(name,"_gbitneg")) {
		       valence=199;
		       func=(void (*)(void*)) _gbitneg;
		   } 
		   break;
#endif
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
	 XS((*subaddr));
	 char* file = __FILE__, *proto = NULL;
	 char subname[276]="Math::Pari::";
	 char buf[64], *s, *s1;
	 CV *protocv;
	 int flexible = 0;
	 
	 sprintf(buf, "%d", valence);
	 switch (valence) {
	 case 0:
	     if (!ep->code) {
		 croak("Unsupported Pari function %s, interface 0 code NULL");
	     } else if (ep->code[0] == 'p' && ep->code[1] == 0) {
		 DO_INTERFACE(0);
	     } else if (ep->code[0] == 0) {
		 DO_INTERFACE(00);
	     } else {
		 goto flexible;
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

	 default: 
	     if (!ep->code)
		 croak("Unsupported interface %d and no code for a Pari function %s",
		       valence, olds);
	   flexible:
	     s1 = s = ep->code;
	     if (*s1 == 'x')
		 s1++;
	     if (*s1 == 'v') {
		 strcpy(buf, "_flexible_void");
		 DO_INTERFACE(_flexible_void);
	     }
	     else if (*s1 == 'l') {
		 strcpy(buf, "_flexible_long");
		 DO_INTERFACE(_flexible_long);
	     }
	     else {
		 strcpy(buf, "_flexible_gen");
		 DO_INTERFACE(_flexible_gen);
	     }
	     
	     flexible = 1;
	 }
	 strcpy(subname+12,"interface");
	 strcpy(subname+12+9,buf);
	 protocv = perl_get_cv(subname, FALSE);
	 if (protocv) {
	     proto = SvPV((SV*)protocv,na);
	 }
	 
	 strcpy(subname+12,olds);
	 RETVAL = newXS(subname,subaddr,file);
	 if (proto)
	     sv_setpv((SV*)RETVAL, proto);
	 XSINTERFACE_FUNC_SET(RETVAL, flexible ? (void*)ep : (void*)func);
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
       entree *ep, *table = functions_basic;
       int i=-1;

       while (++i <= 1) {
	   if (i==1)
	       table = functions_highlevel;
	   
	   for(ep = table; ep->name; ep++)  {
	       valence = EpVALENCE(ep);
	       if (tag == -1 || ep->menu == tag) {
		   switch (valence) {
		   default:
		   case 0:
		       if (ep->code == 0)
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

   primelimit = SvIV(pri);
   parisize = SvIV(mem);
   init(parisize, primelimit); /* Default: take four million bytes of
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
   fmt_nb = def_fmt_nb;
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
#if 1
	 if (SvMAGICAL(sv) && SvTYPE(sv) == SVt_PVAV) {
	     MAGIC *mg = mg_find(sv, 'P');
	     SV *obj;

		/* Be extra paranoid: is refcount is artificially low? */
	     if (mg && (obj = mg->mg_obj) && SvROK(obj) && SvRV(obj) == sv) {
		 mg->mg_flags &= ~MGf_REFCOUNTED;
		 SvREFCNT_inc(sv);
		 SvREFCNT_dec(obj);
	     }
	 }
#endif
	 SvPVX(sv) = GENheap;		/* To avoid extra free() in moveoff.... */
	 if (type == GENheap)	/* Leave it alone? XXXX */
	 /* break */ ;
	 else if (type == GENmovedOffStack) {/* Know that it _was temporary. */
	     killbloc((GEN)SvIV(sv));	     
	 } else {
	 		/* Still on stack */
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
set_gnuterm(a,b,c=0)
    IV a
    IV b
    IV c

long
setprecision(digits=0)
    long digits

long
setseriesprecision(digits=0)
    long digits

IV
setprimelimit(n = 0)
    IV n

void
int_set_term_ftable(a)
    IV a

long
pari_version_exp()

# Cannot do this: it is xsubpp which needs the typemap entry for UV,
# and it needs to convert *all* the branches.
#/* #if defined(PERL_VERSION) && (PERL_VERSION >= 6)*//* 5.6.0 has UV in the typemap */

#if 0
#UV
#allocatemem(newsize = 0)
#UV newsize

#else	/* !( HAVE_UVs ) */

unsigned long
allocatemem(newsize = 0)
    unsigned long newsize

#endif	/* !( HAVE_UVs ) */

long
lgef(x)
    GEN x

long
lgefint(x)
    GEN x

long
lg(x)
    GEN x

unsigned long
longword(x,n)
    GEN x
    long n

MODULE = Math::Pari PACKAGE = Math::Pari	PREFIX = s_

char *
s_type_name(x)
    GEN x
