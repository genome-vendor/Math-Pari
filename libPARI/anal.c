/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*            Analyseur syntactique pour la calculette             */
/*                                                                 */
/*                       copyright Babe Cool                       */
/*                                                                 */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "genpari.h"

GEN seq(void), expr(void), facteur(void), truc(void), identifier(void), constante();
void skipseq(void), skipexpr(void), skipfacteur(void), skiptruc(void), skipconstante(void), skipidentifier(void);
entree *findentry(void), *skipentry(void);

static char *analyseurs,*labellist[100];
static long analyseurtetpil;

GEN
lisexpr(char *t)
{
  GEN res;
  long av, oldtetpil = analyseurtetpil;
  char *olds = analyseurs;

  if (foreignExprHandler && *t == foreignExprSwitch) {
      return (*foreignExprHandler)(t);
  }
  analyseurs = t; analyseurtetpil = av = avma;
  res = expr();
  res = gerepile(av, analyseurtetpil, res);      
  analyseurtetpil = oldtetpil; analyseurs = olds;
  return res;
}

GEN
readexpr(char **c)
{
  char *olds = analyseurs, *oldc = *c;
  analyseurs = oldc; skipexpr();
  if ((*analyseurs) && !separe(*analyseurs)) err(caracer1, analyseurs);
  *c = analyseurs; analyseurs = olds;
  return lisexpr(oldc);
}

GEN
lisseq(char *t)
{
  GEN res;
  long av, oldtetpil = analyseurtetpil;
  char *olds = analyseurs;

  if (foreignExprHandler && *t == foreignExprSwitch) {
      return (*foreignExprHandler)(t);
  }
  analyseurs = t; analyseurtetpil = av = avma;
  res = seq();
  res = gerepile(av, analyseurtetpil, res);      
  analyseurtetpil = oldtetpil; analyseurs = olds;
  return res;
}

GEN
readseq(char **c)
{
  long i;
  char *olds = analyseurs, *oldc = *c;
  for(i=0;i<100;i++) labellist[i]=(char*)0;
  analyseurs = oldc; skipseq();
  *c = analyseurs; analyseurs = olds;
  return lisseq(oldc);
}

/*
 * Values of different codes:
 *
 * p Argument in C is precision
 * l Return long
 * v void return
 * | Interchange last two args
 * ~ Interchange last two pairs of args
 * f Fake *long argument
 * F Fake *GEN argument
 * P Argument in C is series precision
 * = Separator is =
 * G Argument is GEN
 * L Argument is long, XXX 64bit?
 * Z Argument must be converted to long
 * n Argument is ordinal of var
 * V Argument is variable, old value is irrelevant.
 * S Argument is a symbol
 * I Input position (should be processed with lisexpr or lisseq only).
 * s Argument is string
 * D Has a default value (after it follows ",", default value, followed by
 *   ',', followed by the type designator and ",", like "D,567,Z,")
 * d Has a default value which is the value of the previous argument,
 *   like dG. (Make a GEN, if not present, copy the previous argument).
 * x Installed foreign function. Put the ep of the function as an argument,
 *   call installedHandler instead.
 * 
 * If the code is NULL, the previous mechanism is used.
 *
 * When new codes are added, changes should be made to identifier and
 * skipidentifier.
 *
 * Currently the following functions have no code word:
 * O 50
 * goto 61
 * if 80
 * label 60
 * o 50
 * pprint 54
 * pprint1 52
 * print 53
 * print1 51
 * read 56
 * string 57
 * texprint 55
 * until 82
 * while 81
 *
 * Unsupported types: 50-57, 60-61, 80-82.
 *
 * Type 99 reserved for functions without the optimized indentifier code.
 */

entree fonctions[]={
{"O",50,0,7,0,NULL,NULL},
{"abs",1,(void *)gabs,3,0,"Gp",NULL},
{"acos",1,(void *)gacos,3,0,"Gp",NULL},
{"acosh",1,(void *)gach,3,0,"Gp",NULL},
{"addell",3,(void *)addell,5,0,"GGGp",NULL},
{"addhelp",99,(void *)addhelp,11,0,"vSs",NULL},
{"addprimes",1,(void *)addprimestotable,4,0,"Gp",NULL},
{"adj",1,(void *)adj,8,0,"Gp",NULL},
{"agm",2,(void *)agm,3,0,"GGp",NULL},
{"akell",2,(void *)akell,5,0,"GGp",NULL},
{"algdep",23,(void *)algdep,8,0,"GLp",NULL},
{"algdep2",33,(void *)algdep2,8,0,"GLLp",NULL},
{"algtobasis",2,(void*)algtobasis,6,0,"GGp",NULL},
{"allocatemem",11,(void *)allocatemem,2,0,"Lp",NULL},
{"anell",23,(void *)anell,5,0,"GLp",NULL},
{"apell",2,(void *)apell,5,0,"GGp",NULL},
{"apell2",2,(void *)apell2,5,0,"GGp",NULL},
{"apprpadic",2,(void *)apprgen9,7,0,"GGp",NULL},
{"arg",1,(void *)garg,3,0,"Gp",NULL},
{"asin",1,(void *)gasin,3,0,"Gp",NULL},
{"asinh",1,(void *)gash,3,0,"Gp",NULL},
{"assmat",1,(void *)assmat,8,0,"Gp",NULL},
{"atan",1,(void *)gatan,3,0,"Gp",NULL},
{"atanh",1,(void *)gath,3,0,"Gp",NULL},
{"basis",13,(void *)base,6,0,"Gf",NULL},
{"basis2",13,(void *)base2,6,0,"Gf",NULL},
{"basistoalg",2,(void*)basistoalg,6,0,"GGp",NULL},  
{"bernreal",11,(void *)bernreal,3,0,"Lp",NULL},
{"bernvec",11,(void *)bernvec,3,0,"Lp",NULL},
{"bestappr",2,(void *)bestappr,4,0,"GGp",NULL},
{"bezout",2,(void *)vecbezout,4,0,"GGp",NULL},
{"bezoutres",2,(void *)vecbezoutres,4,0,"GGp",NULL},
{"bigomega",1,(void *)gbigomega,4,0,"Gp",NULL},
{"bilhell",3,(void *)bilhell,5,0,"GGGp",NULL},
{"bin",21,(void *)binome,4,0,"GL",NULL},
{"binary",1,(void *)binaire,2,0,"Gp",NULL},
{"bittest",2,(void *)gbittest,2,0,"GGp",NULL},
{"boundcf",21,(void *)gboundcf,4,0,"GL",NULL},
{"boundfact",21,(void *)boundfact,4,0,"GL",NULL},
{"box",35,(void *)rectbox,10,0,"LGG",NULL},
{"buchcertify",10,(void *)certifybuchall,6,0,"Gl",NULL},
{"buchfu",1,(void *)buchfu,6,0,"Gp",NULL},
{"buchgen",92,(void *)buchall,6,0,"GD,0.3,G,dGD,5,G,D,1,G,D,4,L,D,3,L,D,0,L,p",NULL},
{"buchgenforcefu",95,(void *)buchall,6,0,"GD,0.3,G,dGD,5,G,D,1,G,D,4,L,D,3,L,D,3,L,p",NULL},
{"buchgenfu",94,(void *)buchall,6,0,"GD,0.3,G,dGD,5,G,D,1,G,D,4,L,D,3,L,D,2,L,p",NULL},
{"buchimag",96,(void *)buchimag,4,0,"GD,0.1,G,dGD,5,G,",NULL},
{"buchinit",91,(void *)buchall,6,0,"GD,0.3,G,dGD,5,G,D,1,G,D,4,L,D,3,L,D,-1,L,p",NULL},
{"buchinitforcefu",89,(void *)buchall,6,0,"GD,0.3,G,dGD,5,G,D,1,G,D,4,L,D,3,L,D,-3,L,p",NULL},
{"buchinitfu",90,(void *)buchall,6,0,"GD,0.3,G,dGD,5,G,D,1,G,D,4,L,D,3,L,D,-2,L,p",NULL},
{"buchnarrow",1,(void *)buchnarrow,6,0,"Gp",NULL},
{"buchray",2,(void *)buchray,6,0,"GGp",NULL},
{"buchrayinit",2,(void *)buchrayinit,6,0,"GGp",NULL},
{"buchreal",97,(void *)buchreal,4,0,"GD,0,G,D,0.1,G,dGD,5,G,p",NULL},
{"bytesize",10,(void *)taille2,2,0,"Gl",NULL},
{"ceil",1,(void *)gceil,2,0,"Gp",NULL},
{"centerlift",1,(void *)centerlift,2,0,"Gp",NULL},
{"cf",1,(void *)gcf,4,0,"Gp",NULL},
{"cf2",2,(void *)gcf2,4,0,"GGp",NULL},
{"changevar",2,(void *)changevar,2,0,"GGp",NULL},
{"char",14,(void *)caradj0,8,0,"Gn",NULL},
{"char1",14,(void *)caract,8,0,"Gn",NULL},
{"char2",14,(void *)carhess,8,0,"Gn",NULL},
{"chell",2,(void *)coordch,5,0,"GGp",NULL},
{"chinese",2,(void *)chinois,4,0,"GGp",NULL},
{"chptell",2,(void *)pointch,5,0,"GGp",NULL},
{"classno",1,(void *)classno,4,0,"Gp",NULL},
{"classno2",1,(void *)classno2,4,0,"Gp",NULL},
{"coeff",21,(void *)truecoeff,2,0,"GL",NULL},
{"compimag",2,(void *)compimag,4,0,"GGp",NULL},
{"compo",21,(void *)compo,2,0,"GL",NULL},
{"compositum",2,(void *)compositum,6,0,"GGp",NULL},
{"comprealraw",2,(void *)comprealraw,4,0,"GGp",NULL},
{"concat",2,(void *)concat,8,0,"GGp",NULL},
{"conj",1,(void *)gconj,2,0,"Gp",NULL},
{"conjvec",1,(void *)conjvec,2,0,"Gp",NULL},
{"content",1,(void *)content,4,0,"Gp",NULL},
{"convol",2,(void *)convol,7,0,"GGp",NULL},
{"cos",1,(void *)gcos,3,0,"Gp",NULL},
{"cosh",1,(void *)gch,3,0,"Gp",NULL},
{"cursor",11,(void*)rectcursor,10,0,"Lp",NULL},
{"cvtoi",13,(void *)gcvtoi,2,0,"Gf",NULL},
{"cyclo",11,(void *)cyclo,7,0,"Lp",NULL},
{"denom",1,(void *)denom,2,0,"Gp",NULL},
{"deplin",1,(void *)deplin,8,0,"Gp",NULL},
{"deriv",14,(void *)deriv,7,0,"Gn",NULL},
{"det",1,(void *)det,8,0,"Gp",NULL},
{"det2",1,(void *)det2,8,0,"Gp",NULL},
{"detint",1,(void *)detint,8,0,"Gp",NULL},
{"detr",1,(void *)detreel,8,0,"Gp",NULL},
{"dilog",1,(void *)dilog,3,0,"Gp",NULL},
{"dirdiv",2,(void *)dirdiv,7,0,"GGp",NULL},
{"dirmul",2,(void *)dirmul,7,0,"GGp",NULL},
{"dirzetak",2,(void *)dirzetak,6,0,"GGp",NULL},
{"disc",1,(void *)discsr,7,0,"Gp",NULL},
{"discf",1,(void *)discf,6,0,"Gp",NULL},
{"discf2",1,(void *)discf2,6,0,"Gp",NULL},
{"divisors",1,(void *)divisors,4,0,"Gp",NULL},
{"divres",2,(void *)gdiventres,1,0,"GGp",NULL},
{"divsum",22,(void *)divsomme,9,0,"GV|I",NULL},
{"draw",1,(void*)rectdraw,10,0,"Gp",NULL},
{"eigen",1,(void *)eigen,8,0,"Gp",NULL},
{"eint1",1,(void *)eint1,3,0,"Gp",NULL},
{"erfc",1,(void *)gerfc,3,0,"Gp",NULL},
{"eta",1,(void *)eta,3,0,"Gp",NULL},
{"euler",0,(void *)mpeuler,3,0,"p",NULL},
{"eval",1,(void *)geval,7,0,"Gp",NULL},
{"exp",1,(void *)gexp,3,0,"Gp",NULL},
{"extract",2,(void *)extract,8,0,"GGp",NULL},
{"fact",11,(void *)mpfactr,4,0,"Lp",NULL},
{"factcantor",2,(void *)factcantor,4,0,"GGp",NULL},
{"factfq",3,(void *)factmod9,4,0,"GGGp",NULL},
{"factmod",2,(void *)factmod,4,0,"GGp",NULL},
{"factor",1,(void *)factor,4,0,"Gp",NULL},
{"factoredbasis",28,(void *)factoredbase,6,0,"GGf",NULL},
{"factoreddiscf",2,(void *)factoreddiscf,6,0,"GGp",NULL},
{"factoredpolred",2,(void *)factoredpolred,6,0,"GGp",NULL},
{"factoredpolred2",2,(void *)factoredpolred2,6,0,"GGp",NULL},
{"factornf",2,(void *)polfnf,6,0,"GGp",NULL},
{"factorpadic",32,(void *)factorpadic4,7,0,"GGL",NULL},
{"factorpadic2",32,(void *)factorpadic2,7,0,"GGL",NULL},
{"factpol",33,(void *)factpol,7,0,"GLLp",NULL},
{"factpol2",21,(void *)factpol2,7,0,"GL",NULL},
{"fibo",11,(void *)fibo,4,0,"Lp",NULL},
{"floor",1,(void *)gfloor,2,0,"Gp",NULL},
{"for",83,(void *)forpari,10,0,"V=GGI",NULL},
{"fordiv",84,(void *)fordiv,10,0,"GV|I",NULL},
{"forprime",83,(void *)forprime,10,0,"V=GGI",NULL},
{"forstep",86,(void *)forstep,10,0,"V=GGGI",NULL},
{"forvec",87,(void *)forvec,10,0,"V=GIp",NULL},
{"frac",1,(void *)gfrac,2,0,"Gp",NULL},
{"galois",1,(void *)galois,6,0,"Gp",NULL},
{"galoisapply",3,(void *)galoisapply,6,0,"GGGp",NULL},
{"galoisconj",1,(void *)galoisconj,6,0,"Gp",NULL},
{"galoisconj1",1,(void *)galoisconj1,6,0,"Gp",NULL},
{"galoisconjforce",1,(void *)galoisconjforce,6,0,"Gp",NULL},
{"gamh",1,(void *)ggamd,3,0,"Gp",NULL},
{"gamma",1,(void *)ggamma,3,0,"Gp",NULL},
{"gauss",2,(void *)gauss,8,0,"GGp",NULL},
{"gcd",2,(void *)ggcd,4,0,"GGp",NULL},
{"getheap",0,(void *)getheap,2,0,"p",NULL},
{"getrand",0,(void *)getrand,2,0,"p",NULL},
{"getstack",0,(void *)getstack,2,0,"p",NULL},
{"gettime",0,(void *)gettime,2,0,"p",NULL},
{"globalred",1,(void *)globalreduction,5,0,"Gp",NULL},
{"goto",61,0,11,0,NULL,NULL},
{"hclassno",1,(void *)classno3,4,0,"Gp",NULL},
{"hell",2,(void *)ghell,5,0,"GGp",NULL},
{"hell2",2,(void *)ghell2,5,0,"GGp",NULL},
{"hermite",1,(void *)hnf,8,0,"Gp",NULL},
{"hermitebatut",1,(void *)hnfnew,8,0,"Gp",NULL},
{"hermitehavas",1,(void *)hnfhavas,8,0,"Gp",NULL},
{"hermitemod",2,(void *)hnfmod,8,0,"GGp",NULL},
{"hermitemodid",2,(void *)hnfmodid,8,0,"GGp",NULL},
{"hermiteperm",1,(void *)hnfperm,8,0,"Gp",NULL},
{"hess",1,(void *)hess,8,0,"Gp",NULL},
{"hilb",30,(void *) hil,4,0,"GGGl",NULL},
{"hilbert",11,(void *)hilb,8,0,"Lp",NULL},
{"hilbp",20,(void *) hil,4,0,"GGl",NULL},
{"hvector",22,(void *)vecteur,9,0,"GV|I",NULL},
{"hyperu",3,(void *)hyperu,3,0,"GGGp",NULL},
{"i",0,(void *)geni,2,0,"p",NULL},
{"idealadd",3,(void *)idealadd,6,0,"GGGp",NULL},
{"idealaddone",3,(void *)idealaddone,6,0,"GGGp",NULL},
{"idealaddmultone",2,(void *)idealaddmultone,6,0,"GGp",NULL},  
{"idealappr",2,(void *)idealappr,6,0,"GGp",NULL},  
{"idealapprfact",2,(void *)idealapprfact,6,0,"GGp",NULL},
{"idealchinese",3,(void *)idealchinese,6,0,"GGGp",NULL},  
{"idealcoprime",3,(void *)idealcoprime,6,0,"GGGp",NULL},  
{"idealdiv",3,(void *)idealdiv,6,0,"GGGp",NULL},
{"idealdivexact",3,(void *)idealdivexact,6,0,"GGGp",NULL},
{"idealfactor",2,(void *)idealfactor,6,0,"GGp",NULL},
{"idealhermite",2,(void *)idealhermite,6,0,"GGp",NULL},  
{"idealhermite2",3,(void *)idealhermite2,6,0,"GGGp",NULL},  
{"idealintersect",3,(void *)idealintersect,6,0,"GGGp",NULL},
{"idealinv",2,(void *)idealinv,6,0,"GGp",NULL},
{"idealinv2",2,(void *)oldidealinv,6,0,"GGp",NULL},
{"ideallllred",3,(void *)ideallllred,6,0,"GGGp",NULL},
{"idealmul",3,(void *)idealmul,6,0,"GGGp",NULL},
{"idealmulred",3,(void *)idealmulred,6,0,"GGGp",NULL},
{"idealnorm",2,(void *)idealnorm,6,0,"GGp",NULL},
{"idealpow",3,(void *)idealpow,6,0,"GGGp",NULL},
{"idealpowred",3,(void *)idealpowred,6,0,"GGGp",NULL},
{"idealtwoelt",2,(void *)ideal_two_elt,6,0,"GGp",NULL},
{"idealtwoelt2",3,(void *)ideal_two_elt2,6,0,"GGGp",NULL},
{"idealval",30,(void *)idealval,6,0,"GGGl",NULL},
{"idmat",11,(void *)idmat,8,0,"Lp",NULL},
{"if",80,0,11,0,NULL,NULL},
{"imag",1,(void *)gimag,2,0,"Gp",NULL},
{"image",1,(void *)image,8,0,"Gp",NULL},
{"image2",1,(void *)image2,8,0,"Gp",NULL},
{"imagecompl",1,(void *)imagecompl,8,0,"Gp",NULL},
{"imager",1,(void *)imagereel,8,0,"Gp",NULL},
{"incgam",2,(void *)incgam,3,0,"GGp",NULL},
{"incgam1",2,(void *)incgam1,3,0,"GGp",NULL},
{"incgam2",2,(void *)incgam2,3,0,"GGp",NULL},
{"incgam3",2,(void *)incgam3,3,0,"GGp",NULL},
{"incgam4",3,(void *)incgam4,3,0,"GGGp",NULL},
{"indexrank",1,(void *)indexrank,8,0,"Gp",NULL},
{"indsort",1,(void *)indexsort,8,0,"Gp",NULL},
{"initalg",1,(void *)initalg,6,0,"Gp",NULL},
{"initalgred",1,(void *)initalgred,6,0,"Gp",NULL},
{"initalgred2",1,(void *)initalgred2,6,0,"Gp",NULL},
{"initell",1,(void *)initell,5,0,"Gp",NULL},
{"initell2",1,(void *)initell2,5,0,"Gp",NULL},
{"initrect",34,(void*)initrect,10,0,"ZZZ",NULL},
{"initzeta",1,(void *)initzeta,6,0,"Gp",NULL},
{"integ",14,(void *)integ,7,0,"Gn",NULL},
{"intersect",2,(void *)intersect,8,0,"GGp",NULL},
{"intgen",37,(void *)rombint,9,0,"V=GGIp",NULL},
{"intinf",37,(void *)qromi,9,0,"V=GGIp",NULL},
{"intnum",37,(void *)qromb,9,0,"V=GGIp",NULL},
{"intopen",37,(void *)qromo,9,0,"V=GGIp",NULL},
{"inverseimage",2,(void *)inverseimage,8,0,"GGp",NULL},
{"isfund",1,(void *)gisfundamental,4,0,"Gp",NULL},
{"isideal",20,(void *)isideal,6,0,"GGl",NULL},
{"isincl",2,(void *)nfincl,6,0,"GGp",NULL},
{"isinclfast",2,(void *)isinclfast,6,0,"GGp",NULL},  
{"isirreducible",1,(void *)gisirreducible,7,0,"Gp",NULL},
{"isisom",2,(void *)nfiso,6,0,"GGp",NULL},
{"isisomfast",2,(void *)isisomfast,6,0,"GGp",NULL},  
{"isoncurve",20,(void *)oncurve,5,0,"GGl",NULL},
{"isprime",1,(void *)gisprime,4,0,"Gp",NULL},
{"isprincipal",2,(void *)isprincipal,6,0,"GGp",NULL},
{"isprincipalgen",2,(void *)isprincipalgen,6,0,"GGp",NULL},
{"isprincipalray",2,(void *)isprincipalray,6,0,"GGp",NULL},
{"isprincipalraygen",2,(void *)isprincipalraygen,6,0,"GGp",NULL},
{"ispsp",1,(void *)gispsp,4,0,"Gp",NULL},
{"isqrt",1,(void *)racine,4,0,"Gp",NULL},
{"isset",10,(void *)isvecset,8,0,"Gl",NULL},
{"issqfree",1,(void *)gissquarefree,4,0,"Gp",NULL},
{"issquare",1,(void *)gcarreparfait,4,0,"Gp",NULL},
{"isunit",2,(void *)isunit,6,0,"GGp",NULL},
{"jacobi",1,(void *)jacobi,8,0,"Gp",NULL},
{"jbesselh",2,(void *)jbesselh,3,0,"GGp",NULL},
{"jell",1,(void *)jell,3,0,"Gp",NULL},
{"karamul",32,(void *)karamul,7,0,"GGL",NULL},
{"kbessel",2,(void *)kbessel,3,0,"GGp",NULL},
{"kbessel2",2,(void *)kbessel2,3,0,"GGp",NULL},
{"ker",1,(void *)ker,8,0,"Gp",NULL},
{"keri",1,(void *)keri,8,0,"Gp",NULL},
{"kerint",1,(void *)kerint,8,0,"Gp",NULL},
{"kerint1",1,(void *)kerint1,8,0,"Gp",NULL},
{"kerint2",1,(void *)kerint2,8,0,"Gp",NULL},
{"kerr",1,(void *)kerreel,8,0,"Gp",NULL},
{"kill",85,(void*)killep,11,0,"S",NULL},
{"killrect",11,(void *)killrect,10,0,"Lp",NULL},
{"kro",2,(void *)gkronecker,4,0,"GGp",NULL},
{"label",60,0,11,0,NULL,NULL},
{"lambdak",2,(void *)glambdak,6,0,"GGp",NULL},
{"laplace",1,(void *)laplace,7,0,"Gp",NULL},
{"lcm",2,(void *)glcm,4,0,"GGp",NULL},
{"legendre",11,(void *)legendre,7,0,"Lp",NULL},
{"length",1,(void *)glength,2,0,"Gp",NULL},
{"lex",20,(void *)lexcmp,2,0,"GGl",NULL},
{"lexsort",1,(void *)lexsort,8,0,"Gp",NULL},
{"lift",1,(void *)lift,2,0,"Gp",NULL},
{"lindep",1,(void *)lindep,8,0,"Gp",NULL},
{"lindep2",23,(void *)lindep2,8,0,"GLp",NULL},
{"line",35,(void *)rectline,10,0,"LGG",NULL},
{"lines",35,(void *)rectlines,10,0,"LGG",NULL},
{"lll",1,(void *)lll,8,0,"Gp",NULL},
{"lll1",1,(void *)lll1,8,0,"Gp",NULL},
{"lllgen",1,(void *)lllgen,8,0,"Gp",NULL},
{"lllgram",1,(void *)lllgram,8,0,"Gp",NULL},
{"lllgram1",1,(void *)lllgram1,8,0,"Gp",NULL},
{"lllgramgen",1,(void *)lllgramgen,8,0,"Gp",NULL},
{"lllgramint",1,(void *)lllgramint,8,0,"Gp",NULL},
{"lllgramkerim",1,(void *)lllgramkerim,8,0,"Gp",NULL},
{"lllgramkerimgen",1,(void *)lllgramkerimgen,8,0,"Gp",NULL},
{"lllint",1,(void *)lllint,8,0,"Gp",NULL},
{"lllintpartial",1,(void *)lllintpartial,8,0,"Gp",NULL},
{"lllkerim",1,(void *)lllkerim,8,0,"Gp",NULL},
{"lllkerimgen",1,(void *)lllkerimgen,8,0,"Gp",NULL},
{"lllrat",1,(void *)lllrat,8,0,"Gp",NULL},
{"ln",1,(void *)glog,3,0,"Gp",NULL},
{"lngamma",1,(void *)glngamma,3,0,"Gp",NULL},
{"localred",2,(void *)localreduction,5,0,"GGp",NULL},
{"log",1,(void *)glog,3,0,"Gp",NULL},
{"logagm",1,(void *)glogagm,3,0,"Gp",NULL},
{"lseriesell",4,(void *)lseriesell,5,0,"GGGGp",NULL},
{"mat",1,(void *)gtomat,8,0,"Gp",NULL},
{"matextract",3,(void *)matextract,8,0,"GGGp",NULL},
{"mathell",2,(void *)mathell,5,0,"GGp",NULL},
{"matinvr",1,(void *)invmatreel,8,0,"Gp",NULL},
{"matrix",49,(void *)matrice,9,0,"GGVV~I",NULL},
{"matrixqz",2,(void *)matrixqz,8,0,"GGp",NULL},
{"matrixqz2",1,(void *)matrixqz2,8,0,"Gp",NULL},
{"matrixqz3",1,(void *)matrixqz3,8,0,"Gp",NULL},
{"matsize",1,(void *)matsize,2,0,"Gp",NULL},
{"max",2,(void *)gmax,1,0,"GGp",NULL},
{"min",2,(void *)gmin,1,0,"GGp",NULL},
{"minideal",3,(void *)minideal,6,0,"GGGp",NULL},
{"minim",33,(void *)minim,8,0,"GLLp",NULL},
{"minim2",23,(void *)minim2,8,0,"GLp",NULL},
{"mod",25,(void *)gmodulcp,2,0,"GG",NULL},
{"modp",25,(void *)gmodulo,2,0,"GG",NULL},
{"modreverse",1,(void *)polymodrecip,6,0,"Gp",NULL},
{"move",35,(void*)rectmove,10,0,"LGG",NULL},
{"mu",1,(void *)gmu,4,0,"Gp",NULL},
{"newtonpoly",2,(void *)newtonpoly,6,0,"GGp",NULL},
{"nextprime",1,(void *)bigprem,4,0,"Gp",NULL},
{"nfdetint",2,(void *)nfdetint,6,0,"GGp",NULL},
{"nfdiv",3,(void *)element_div,6,0,"GGGp",NULL},
{"nfdiveuc",3,(void *)nfdiveuc,6,0,"GGGp",NULL},
{"nfdivres",3,(void *)nfdivres,6,0,"GGGp",NULL},
{"nfhermite",2,(void *)nfhermite,6,0,"GGp",NULL},
{"nfhermitemod",3,(void *)nfhermitemod,6,0,"GGGp",NULL},
{"nfmod",3,(void *)nfmod,6,0,"GGGp",NULL},
{"nfmul",3,(void *)element_mul,6,0,"GGGp",NULL},
{"nfpow",3,(void *)element_pow,6,0,"GGGp",NULL},
{"nfreduce",3,(void *)element_reduce,6,0,"GGGp",NULL},
{"nfsmith",2,(void *)nfsmith,6,0,"GGp",NULL},  
{"nfval",30,(void *)element_val,6,0,"GGGl",NULL},
{"norm",1,(void *)gnorm,2,0,"Gp",NULL},
{"norml2",1,(void *)gnorml2,2,0,"Gp",NULL},
{"nucomp",3,(void *)nucomp,4,0,"GGGp",NULL},
{"numdiv",1,(void *)numbdiv,4,0,"Gp",NULL},
{"numer",1,(void *)numer,2,0,"Gp",NULL},
{"nupow",2,(void *)nupow,4,0,"GGp",NULL},
{"o",50,0,7,0,NULL,NULL},
{"omega",1,(void *)gomega,4,0,"Gp",NULL},
{"ordell",2,(void *)ordell,5,0,"GGp",NULL},
{"order",1,(void *)order,4,0,"Gp",NULL},
{"orderell",2,(void *)orderell,5,0,"GGp",NULL},
{"ordred",1,(void *)ordred,6,0,"Gp",NULL},
{"padicprec",20,(void *)padicprec,2,0,"GGl",NULL},
{"pascal",11,(void *)pasc,8,0,"Lp",NULL},
{"perf",10,(void *)perf,8,0,"Gl",NULL},
{"permutation",24,(void *)permute,2,0,"LGp",NULL},
{"permutation2num",1,(void *)permuteInv,2,0,"Gp",NULL},
{"pf",2,(void *)primeform,4,0,"GGp",NULL},
{"phi",1,(void *)phi,4,0,"Gp",NULL},
{"pi",0,(void *)mppi,3,0,"p",NULL},
{"plot",37,(void *)plot,10,0,"V=GGIp",NULL},
{"ploth",37,(void *)ploth,10,0,"V=GGIp",NULL},
{"ploth2",37,(void *)ploth2,10,0,"V=GGIp",NULL},
{"plothmult",37,(void *)plothmult,10,0,"V=GGIp",NULL},
{"plothraw",2,(void *)plothraw,10,0,"GGp",NULL},
{"plothrawlines",2,(void *)plothrawlines,10,0,"GGp",NULL},
{"plothsizes",0,(void *)plothsizes,10,0,"p",NULL},
{"plotterm",16,(void *)term_set,10,0,"sl",NULL},
{"pnqn",1,(void *)pnqn,4,0,"Gp",NULL},
{"point",35,(void *)rectpoint,10,0,"LGG",NULL},
{"pointell",2,(void *)pointell,5,0,"GGp",NULL},
{"points",35,(void *)rectpoints,10,0,"LGG",NULL},
{"polint",31,(void *)polint,7,0,"GGGF",NULL},
{"polred",1,(void *)polred,6,0,"Gp",NULL},
{"polred2",1,(void *)polred2,6,0,"Gp",NULL},
{"polredabs",1,(void *)polredabs,6,0,"Gp",NULL},
{"polredabs2",1,(void *)polredabs2,6,0,"Gp",NULL},
{"polredabsall",1,(void *)polredabsall,6,0,"Gp",NULL},
{"polredabsfast",1,(void *)polredabsfast,6,0,"Gp",NULL},
{"polsym",21,(void *)polsym,7,0,"GL",NULL},
{"polvar",1,(void *)gpolvar,2,0,"Gp",NULL},
{"poly",14,(void *)gtopoly,2,0,"Gn",NULL},
{"polylog",24,(void *)gpolylog,3,0,"LGp",NULL},
{"polylogd",24,(void *)polylogd,3,0,"LGp",NULL},
{"polylogdold",24,(void *)polylogdold,3,0,"LGp",NULL},
{"polylogp",24,(void *)polylogp,3,0,"LGp",NULL},
{"polyrev",14,(void *)gtopolyrev,2,0,"Gn",NULL},
{"postdraw",1,(void *)postdraw,10,0,"Gp",NULL},
{"postploth",37,(void *)postploth,10,0,"V=GGIp",NULL},
{"postploth2",37,(void *)postploth2,10,0,"V=GGIp",NULL},
{"postplothraw",2,(void *)postplothraw,10,0,"GGp",NULL},
{"powell",3,(void *)powell,5,0,"GGGp",NULL},
{"powrealraw",23,(void *)powrealraw,4,0,"GLp",NULL},
{"pprint",54,0,11,0,NULL,NULL},
{"pprint1",52,0,11,0,NULL,NULL},
{"prec",21,(void *)gprec,2,0,"GL",NULL},
{"prime",11,(void *)prime,4,0,"Lp",NULL},
{"primedec",2,(void *)primedec,6,0,"GGp",NULL},
{"primes",11,(void *)primes,4,0,"Lp",NULL},
{"primroot",1,(void *)gener,4,0,"Gp",NULL},
{"principalideal",2,(void *)principalideal,6,0,"GGp",NULL},
{"principalidele",2,(void *)principalidele,6,0,"GGp",NULL},
{"print",53,0,11,0,NULL,NULL},
{"print1",51,0,11,0,NULL,NULL},
{"prod",48,(void *)produit,9,0,"GV|=GGIp",NULL},
{"prodeuler",37,(void *)prodeuler,9,0,"V=GGIp",NULL},
{"prodinf",27,(void *)prodinf,9,0,"V=GIp",NULL},
{"prodinf1",27,(void *)prodinf1,9,0,"V=GIp",NULL},
{"psi",1,(void *)gpsi,3,0,"Gp",NULL},
{"qfi",3,(void *)qfi,4,0,"GGGp",NULL},
{"qfr",4,(void *)qfr,4,0,"GGGGp",NULL},
{"quaddisc",1,(void *)quaddisc,4,0,"Gp",NULL},   
{"quadgen",1,(void *)quadgen,2,0,"Gp",NULL},   
{"quadpoly",1,(void *)quadpoly,2,0,"Gp",NULL},   
{"random",0,(void *)genrand,2,0,"p",NULL},
{"rank",10,(void *)rank,8,0,"Gl",NULL},
{"rbox",35,(void *)rectrbox,10,0,"LGG",NULL},
{"rcopy",44,(void *)rectcopy,10,0,"LLLL",NULL},
{"read",56,0,11,0,NULL,NULL},
{"real",1,(void *)greal,2,0,"Gp",NULL},
{"recip",1,(void *)polrecip,7,0,"Gp",NULL},
{"redimag",1,(void *)redimag,4,0,"Gp",NULL},      
{"redreal",1,(void *)redreal,4,0,"Gp",NULL},      
{"redrealnod",2,(void *)redrealnod,4,0,"GGp",NULL},      
{"regula",1,(void *)regula,4,0,"Gp",NULL}, 
{"reorder",1,(void *)reorder,11,0,"Gp",NULL}, 
{"resultant",2,(void *)subres,7,0,"GGp",NULL},    
{"resultant2",2,(void *)resultant2,7,0,"GGp",NULL},    
{"reverse",1,(void *)recip,7,0,"Gp",NULL}, 
{"rhoreal",1,(void *)rhoreal,4,0,"Gp",NULL},      
{"rhorealnod",2,(void *)rhorealnod,4,0,"GGp",NULL},      
{"rline",35,(void *)rectrline,10,0,"LGG",NULL},
{"rlines",35,(void *)rectlines,10,0,"LGG",NULL},
{"rlinetype",19,(void *)rectlinetype,10,0,"LL",NULL},      
{"rmove",35,(void *)rectrmove,10,0,"LGG",NULL},      
{"rndtoi",13,(void *)grndtoi,2,0,"Gf",NULL},
{"rnfbasis",2,(void *)rnfbasis,6,0,"GGp",NULL},  
{"rnfdiscf",2,(void *)rnfdiscf,6,0,"GGp",NULL},
{"rnfhermitebasis",2,(void *)rnfhermitebasis,6,0,"GGp",NULL},  
{"rnfisfree",20,(void *)rnfisfree,6,0,"GGl",NULL},
{"rnfpseudobasis",2,(void *)rnfpseudobasis,6,0,"GGp",NULL},
{"rnfsteinitz",2,(void *)rnfsteinitz,6,0,"GGp",NULL},
{"rootmod",2,(void *)rootmod,7,0,"GGp",NULL},
{"rootmod2",2,(void *)rootmod2,7,0,"GGp",NULL},
{"rootpadic",32,(void *)rootpadic,7,0,"GGL",NULL},
{"roots",1,(void *)roots,7,0,"Gp",NULL},
{"roots2",1,(void *)roots2,7,0,"Gp",NULL},
{"rootslong",1,(void *)rootslong,7,0,"Gp",NULL},
{"rootsof1",1,(void *)rootsof1,6,0,"Gp",NULL},
{"round",1,(void *)ground,2,0,"Gp",NULL},
{"rounderror",10,(void *)rounderror,2,0,"Gl",NULL},
{"rploth",73,(void *)recplothmultin,10,0,"LV=GGIpLL",NULL},      
{"rplothraw",45,(void *)recplothraw,10,0,"LGGLp",NULL},      
{"rpoint",35,(void *)rectrpoint,10,0,"LGG",NULL},      
{"rpoints",35,(void *)rectpoints,10,0,"LGG",NULL},      
{"rpointtype",19,(void *)rectpointtype,10,0,"LL",NULL},      
{"scale",59,(void *)rectscale,10,0,"LGGGG",NULL},
{"series",14,(void *)gtoser,2,0,"Gn",NULL},
{"set",1,(void *)gtoset,8,0,"Gp",NULL},
{"setintersect",2,(void *)setintersect,8,0,"GGp",NULL},
{"setminus",2,(void *)setminus,8,0,"GGp",NULL},
{"setprecision",15,(void *)setprecr,2,0,"Ll",NULL},
{"setrand",11,(void *)setrand,2,0,"Lp",NULL},
{"setsearch",20,(void *)setsearch,8,0,"GGl",NULL},
{"setserieslength",15,(void *)setserieslength,2,0,"Ll",NULL},
{"settype",21,(void *)gsettype,2,0,"GL",NULL},
{"setunion",2,(void *)setunion,8,0,"GGp",NULL},
{"shift",21,(void *)gshift,1,0,"GL",NULL},
{"shiftmul",21,(void *)gmul2n,1,0,"GL",NULL},
{"sigma",1,(void *)sumdiv,4,0,"Gp",NULL},
{"sigmak",24,(void *)sumdivk,4,0,"LGp",NULL},
{"sign",10,(void *)gsigne,1,0,"Gl",NULL},
{"signat",1,(void *)signat,8,0,"Gp",NULL},
{"signunit",1,(void *)signunit,6,0,"Gp",NULL},
{"simplefactmod",2,(void *)simplefactmod,4,0,"GGp",NULL},
{"simplify",1,(void *)simplify,2,0,"Gp",NULL},
{"sin",1,(void *)gsin,3,0,"Gp",NULL},
{"sinh",1,(void *)gsh,3,0,"Gp",NULL},
{"size",10,(void *)gsize,2,0,"Gl",NULL},
{"smallbasis",13,(void *)smallbase,6,0,"Gf",NULL},
{"smalldiscf",1,(void *)smalldiscf,6,0,"Gp",NULL},
{"smallfact",1,(void *)smallfact,4,0,"Gp",NULL},
{"smallinitell",1,(void *)smallinitell,5,0,"Gp",NULL},
{"smallpolred",1,(void *)smallpolred,6,0,"Gp",NULL},
{"smallpolred2",1,(void *)smallpolred2,6,0,"Gp",NULL},
{"smith",1,(void *)smith,8,0,"Gp",NULL},
{"smith2",1,(void *)smith2,8,0,"Gp",NULL},
{"smithpol",1,(void *)gsmith,8,0,"Gp",NULL},
{"solve",37,(void *)zbrent,9,0,"V=GGIp",NULL},
{"sort",1,(void *)sort,8,0,"Gp",NULL},
{"sqr",1,(void *)gsqr,3,0,"Gp",NULL},
{"sqred",1,(void *)sqred,8,0,"Gp",NULL},
{"sqrt",1,(void *)gsqrt,3,0,"Gp",NULL},
{"srgcd",2,(void *)srgcd,7,0,"GGp",NULL},
{"string",57,(void*)rectstring,10,0,NULL,NULL},
{"sturm",10,(void *)sturm,7,0,"Gl",NULL},
{"sturmpart",30,(void *)sturmpart,7,0,"GGGl",NULL},
{"subcyclo",2,(void *)subcyclo,6,0,"GGp",NULL},
{"subell",3,(void *)subell,5,0,"GGGp",NULL},
{"subst",26,(void *)gsubst,7,0,"GnG",NULL},
{"sum",48,(void *)somme,9,0,"GV|=GGIp",NULL},
{"sumalt",27,(void *)sumalt,9,0,"V=GIp",NULL},
{"sumalt2",27,(void *)sumalt2,9,0,"V=GIp",NULL},
{"suminf",27,(void *)suminf,9,0,"V=GIp",NULL},
{"sumpos",27,(void *)sumpos,9,0,"V=GIp",NULL},
{"supplement",1,(void *)suppl,8,0,"Gp",NULL},
{"sylvestermatrix",2,(void *)sylvestermatrix,7,0,"GGp",NULL},
{"tan",1,(void *)gtan,3,0,"Gp",NULL},
{"tanh",1,(void *)gth,3,0,"Gp",NULL},
{"taniyama",1,(void *)taniyama,5,0,"Gp",NULL},
{"taylor",12,(void *)tayl,7,0,"GnP",NULL},
{"tchebi",11,(void *)tchebi,7,0,"Lp",NULL},
{"teich",1,(void *)teich,3,0,"Gp",NULL},
{"texprint",55,0,11,0,NULL,NULL},
{"theta",2,(void *)theta,3,0,"GGp",NULL},
{"thetanullk",21,(void *)thetanullk,3,0,"GL",NULL},
{"threetotwo",4,(void *)threetotwo,6,0,"GGGGp",NULL},
{"threetotwo2",4,(void *)threetotwo2,6,0,"GGGGp",NULL},
{"torsell",1,(void *)torsell,5,0,"Gp",NULL},
{"trace",1,(void *)trace,8,0,"Gp",NULL},
{"trans",1,(void *)gtrans,8,0,"Gp",NULL},
{"trunc",1,(void *)gtrunc,2,0,"Gp",NULL},
{"tschirnhaus",1,(void *)tschirnhaus,6,0,"Gp",NULL},
{"twototwo",3,(void *)twototwo,6,0,"GGGp",NULL},  
{"type",1,(void *)gtype,2,0,"Gp",NULL},
{"unit",1,(void *)fundunit,4,0,"Gp",NULL},
{"until",82,0,11,0,NULL,NULL},
{"valuation",20,(void *)ggval,2,0,"GGl",NULL},
{"vec",1,(void *)gtovec,2,0,"Gp",NULL},
{"vecmax",1,(void *)vecmax,2,0,"Gp",NULL},
{"vecmin",1,(void *)vecmin,2,0,"Gp",NULL},
{"vecsort",2,(void *)vecsort,8,0,"GGp",NULL},
{"vector",22,(void *)vecteur,9,0,"GV|I",NULL},
{"vvector",22,(void *)vvecteur,9,0,"GV|I",NULL},
{"wf",1,(void *)wf,3,0,"Gp",NULL},
{"wf2",1,(void *)wf2,3,0,"Gp",NULL},
{"while",81,0,11,0,NULL,NULL},
{"zell",2,(void *)zell,5,0,"GGp",NULL},
{"zeta",1,(void *)gzeta,3,0,"Gp",NULL},
{"zetak",2,(void *)gzetak,6,0,"GGp",NULL},
{"zideallog",3,(void *)zideallog,6,0,"GGGp",NULL},
{"zidealstar",2,(void *)zidealstar,6,0,"GGp",NULL},
{"zidealstarinit",2,(void *)zidealstarinit,6,0,"GGp",NULL},
{"znstar",1,(void *)znstar,4,0,"Gp",NULL},
{"zzzz",2,(void *)zidealstarinitold,6,0,"GGp",NULL},
{"zzzzz",10,(void *)certifybuchall,8,0,"Gl",NULL}
};

long    NUMFUNC=sizeof(fonctions)/sizeof(entree);

static void
matcherr(char c)
{
  static char reste[100];
  char *p;
  long i;
  
  for(analyseurs--, p=reste, *p++=c, i=0; i<97; i++) *p++ = *analyseurs++;
  *p = 0;err(matcher1,reste);
}

#define match(c)  if(*analyseurs++ != c) matcherr(c)

GEN
seq(void)
{
  GEN res=gnil;
  for(;;)
  {
    while(separe(*analyseurs)) analyseurs++;
    if ((!*analyseurs) || (*analyseurs == ')') || (*analyseurs == ',')) return res;
    res = expr();
    if(!separe(*analyseurs)) return res;
  }
}

GEN
expr(void)
{
#ifdef __cplusplus
  typedef GEN (*PFGEN)(...);
#else
  typedef GEN (*PFGEN)();
#endif
  PFGEN func[4];
  GEN aux,e,e1,e2,e3;
  long niveau;
  
  for(niveau=0;niveau<4;niveau++) func[niveau]=NULL;
  e1=e2=e3=(GEN)0;
  niveau=3;
  for(;;)
    switch(niveau)
    {
      case 3: aux=facteur();
	if(func[3]) {analyseurtetpil=avma;e3=((GEN (*)(GEN,GEN))func[3])(e3,aux);}
	else e3=aux;
	switch(*analyseurs)
	{
	  case '*': analyseurs++;func[3]=(PFGEN)&gmul;break;
	  case '/': analyseurs++;func[3]=(PFGEN)&gdiv;break;
	  case '\\': analyseurs++;
	    if((*analyseurs)=='/') {analyseurs++;func[3]=(PFGEN)&gdivround;}
	    else func[3]=(PFGEN)&gdivent;
	    break;
	  case '%': analyseurs++;func[3]=(PFGEN)&gmod;break;
	  default: niveau--;func[3]=NULL; 
	}
	break;
      case 2: 
	if(!e3) {niveau++;break;}
	if(func[2]) {analyseurtetpil=avma;e2=((GEN (*)(GEN,GEN))func[2])(e2,e3);}
	else e2=e3;
	e3=(GEN)0;
	switch(*analyseurs)
	{
	  case '+': analyseurs++;func[2]=(PFGEN)&gadd;niveau++;break;
	  case '-': analyseurs++;func[2]=(PFGEN)&gsub;niveau++;break;
	  default: niveau--;func[2]=NULL;
	}
	break;
      case 1: 
	if(!e2) {niveau++;break;}
	if(func[1]) {analyseurtetpil=avma;e1=((GEN (*)(GEN,GEN))func[1])(e1,e2);}
	else e1=e2;
	e2=(GEN)0;
	switch(*analyseurs)
	{
	  case '<': analyseurs++;
	    switch(*analyseurs)
	    {
	      case '=': analyseurs++;func[1]=(PFGEN)&gle;break;
	      case '>': analyseurs++;func[1]=(PFGEN)&gne;break;
	      default : func[1]=(PFGEN)&glt;
	    }
	    niveau++;break;
	  case '>': analyseurs++;
	    if((*analyseurs)=='=') {analyseurs++;func[1]=(PFGEN)&gge;}
	    else func[1]=(PFGEN)&ggt;
	    niveau++;break;
	  case '=': 
	    if((analyseurs[1])=='=') {analyseurs+=2;func[1]=(PFGEN)&geq;niveau++;}
	    break;
	  case '!': 
	    if((analyseurs[1])=='=') {analyseurs+=2;func[1]=(PFGEN)&gne;niveau++;}
	    break;
	  default: niveau--;func[1]=NULL;
	}
	break;
      case 0: 
	if(!e1) {niveau++;break;}
	if(func[0]) {analyseurtetpil=avma;e=((GEN (*)(GEN,GEN))func[0])(e,e1);}
	else e=e1;
	e1=(GEN)0;
	switch(*analyseurs)
	{
	  case '&': analyseurs++;if(*analyseurs=='&') analyseurs++;
	    func[0]=(PFGEN)&gand;niveau++;break;
	  case '|': analyseurs++;if(*analyseurs=='|') analyseurs++;
	    func[0]=(PFGEN)&gor;niveau++;break;
	  default: return e;
	}
    }
}

GEN
facteur(void)
{
  GEN tru,p1,arg,arg1;
  long tx,c,e,av2,flcol,flrow;
  long plus = (*analyseurs =='+')||(*analyseurs =='-')?(*analyseurs++=='+'):1;
  tru=truc();
  for (;;) switch(*analyseurs)
  {
    case '^': analyseurs++;p1=facteur();analyseurtetpil=avma; tru=gpui(tru,p1,prec);break;
    case '~': analyseurs++;analyseurtetpil=avma;tru=gtrans(tru);break;
    case '_': analyseurs++;analyseurtetpil=avma;tru=gconj(tru);break;
    case '\'': analyseurs++;analyseurtetpil=avma;tru=deriv(tru,gvar9(tru));break;
    case '[': 
      tx=typ(p1=tru);
      if((tx<17)||(tx>19)) err(caracer1,analyseurs);
      analyseurs++;av2=avma;flcol=flrow=0;
      if(tx<19)
      {
	arg=expr();if(typ(arg)!=1) err(caseer);
	c=itos(arg);if((c<1)||(c>=lg(p1))) err(arrayer1);
      }
      else
      {
	if(lg(p1)==1) err(arrayer1);
	if(*analyseurs==',')
	{
	  analyseurs++;arg=expr();if(typ(arg)!=1) err(caseer);
	  c=itos(arg);if((c<1)||(c>=lg(p1))) err(arrayer1);
	  flcol=1;
	}
	else
	{
	  arg=expr();if(typ(arg)!=1) err(caseer);
	  e=itos(arg);if((e<1)||(e>=lg((GEN)p1[1]))) err(arrayer1);
	  match(',');
	  if(*analyseurs==']') flrow=1;
	  else
	  {
	    arg1=expr();if(typ(arg1)!=1) err(caseer);
	    c=itos(arg1);
	    if((c<1)||(c>=lg(p1))) err(arrayer1);
	  }
	}
      }
      match(']'); analyseurtetpil=avma=av2;
      if((tx<19)||flcol) tru=gcopy((GEN)p1[c]);
      else
      {
	if(flrow)
	{
	  tru=cgetg(lg(p1),17);
	  for(c=1;c<lg(p1);c++) tru[c]=lcopy((GEN)((GEN)p1[c])[e]);
	}
	else tru = gcopy((GEN)((GEN)p1[c])[e]);
      }
      break;
    case '!': analyseurs++;if((*analyseurs)!='=') {analyseurtetpil=avma;tru=mpfact(itos(tru));break;} else analyseurs--;
    default: if(plus) return tru; else {analyseurtetpil=avma;return gneg(tru);}
  }
}

GEN
truc(void)
{
  long i,n=0,j,p=0,m=1;
  GEN *table,p1;
  
  if (isalpha(*analyseurs)) return identifier();
  if (isdigit(*analyseurs) || (*analyseurs=='.')) return constante();
  switch(*analyseurs++)
  {
    case '(': p1=expr();match(')');return p1;
    case '[':
      table = (GEN *)newbloc(paribuffsize>>TWOPOTBYTES_IN_LONG);
      if (*analyseurs!=']') 
      {do table[++n]=expr();while (*analyseurs++==',');analyseurs--;}
      switch (*analyseurs++)
      {
	case ']': analyseurtetpil=avma;p1=cgetg(n+1,17);
	  for (i=1;i<=n;i++) p1[i]=lcopy(table[i]);
	  break;
	case ';': m=n;do table[++n]=expr();while (*analyseurs++!=']');
	  if (n % m) err(recter1);
	  p=n/m;analyseurtetpil=avma;p1=cgetg(m+1,19);
	  for (i=1;i<=m;i++) p1[i]=(long)cgetg(p+1,18);
	  for (j=1;j<=m;j++)
	    for(i=1;i<=p;i++)
	      ((GEN)p1[j])[i]=lcopy(table[(i-1)*m+j]);
	  break;
	default: err(vectmater1);
      }
      killbloc((GEN)table);
      return p1;
    case '%':
      p=0;while((*analyseurs)=='`') {analyseurs++;p++;}
      if(p>tglobal) err(referer1);
      if(p) return g[tglobal-p];
      while (isdigit(*analyseurs)) p = 10*p + *analyseurs++ - '0';
      if(p>tglobal) err(referer2);
      return g[p];
  }
  err(caracer1,analyseurs-1);return gnil;
}

int
numvar(GEN x)
{
  if(typ(x)!=10) err(numvarer);
  if(lgef(x)!=4) err(numvarer);
  if((!gcmp0((GEN)x[2])) || (!gcmp1((GEN)x[3]))) err(numvarer);
  return varn(x);
}

GEN brutcopy(GEN x, GEN y);

#define TOTAL_STRING_ARG 4096		/* Total lenght of string args */
#define RET_INT 1
#define RET_VOID 2
#define RET_GEN 0

#define match_comma() if (matchcomma) {match(',');} else matchcomma = 1

GEN
identifier(void)
{
  long c,e,va,m,nparam,i,av,av2,tx,flrow,flcol;
  static long yatileugoto;
  GEN arg,arg1,arg2,arg3,arg4,argvec[10],res=gnil,p1;
#ifdef __cplusplus
  GEN (*f)(...);
#else
  GEN (*f)();
#endif
  char *ch1, *ch2, *readbuffer;
  entree *ep, *ep1, **p;
  
  ep = findentry();
  if (ep->code) {
      int ret = RET_GEN, matchcomma = 0;
      static char buffer[TOTAL_STRING_ARG]; /* Do not recurse into string! */
      char *bp=buffer, *sinit, delim, *s = ep->code;
      char *oldanalyseurs = NULL;
      long fake;
      GEN *fakepGEN = NULL;
      GEN fakeGEN;
      void *call = ep->value;
      entree *ep1;
        
      i=0;
      if (!*(s+1) && (*s == 'p') && (*analyseurs != '('))
	  return ((GEN (*)(long))ep->value)(prec);
      match('(');
      /* Optimized for G and p. */
      while (*s == 'G') {
	  match_comma();
	  argvec[i++] = expr();
	  s++;
      }
      if (*s == 'p') {
	  argvec[i++] = (GEN) prec;
	  s++;
      }
      while (*s) {
	  switch (*s++) {
	  default:
	      err(talker,"Unknown code for function arguments");
	  case 'G':		/* Argument is GEN */
	      match_comma();
	      argvec[i++] = expr();
	      break;
	  case 'L':		/* Argument is long, XXX 64bit? */
	      match_comma();
	      arg = expr();
	      if(typ(arg)!=1) err(caseer);
	      argvec[i++] = (GEN) itos(arg);
	      break;
	  case 'Z':		/* Argument must be converted to long */
	      match_comma();
	      argvec[i++] = (GEN) gtolong(expr());
	      break;
	  case 'n':		/* Argument is ordinal of var */
	      match_comma();
	      argvec[i++] = (GEN) numvar(expr());
	      break;
	  case 'V':		/* Argument is a variable */
	  case 'S':		/* Argument is a symbol */
	      match_comma();
	      if(!isalpha(*analyseurs)) err(varer1,analyseurs);
	      ep1 = findentry();
	      if (EpVALENCE(ep1) != 200 && *(s-1) == 'V') 
		  err(varer1,analyseurs);
	      argvec[i++] = (GEN) ep1;
	      break;
	  case  'I':		/* Input position */
	      match_comma();
	      argvec[i++] = (GEN) analyseurs;
	      skipexpr();
	      break;
	  case 's':		/* Argument is string */
	      match_comma();
	      delim = *analyseurs++;
	      if (delim != '"' && delim != '\'') {
		  err(talker,"String should be delimited by \" or \'");
	      }
	      sinit = bp;
	      while ((*analyseurs) && (*analyseurs != delim) 
		     && ((bp - buffer) < TOTAL_STRING_ARG - 1)) 
		  *bp++ = *analyseurs++;
	      match(delim);
	      *bp = 0;
	      bp++;
	      if ((bp - buffer) >= TOTAL_STRING_ARG - 1) 
		  err(talker,"Too long string");
	      argvec[i++] = (GEN) sinit;
	      break;
	  case 'p':		/* Argument in C is precision */
	      argvec[i++] = (GEN) prec;
	      break;
	  case 'l':		/* Return long */
	      ret = RET_INT;
	      break;
	  case '=':
	      match('=');
	      matchcomma = 0;
	      break;
	  case '|':		/* Interchange last two args */
	  {
	      GEN tmp = argvec[i-1];

	      argvec[i-1] = argvec[i-2];
	      argvec[i-2] = tmp;
	      break;
	  }
	  case '~':		/* Interchange last two pairs
				 * of args */
	  {
	      GEN tmp = argvec[i-1];

	      argvec[i-1] = argvec[i-3];
	      argvec[i-3] = tmp;
	      tmp = argvec[i-2];
	      argvec[i-2] = argvec[i-4];
	      argvec[i-4] = tmp;
	      break;
	  }
	  case 'D':		/* Has a default value */
	      if (*analyseurs == ')') {
		  oldanalyseurs = analyseurs;

		  analyseurs = matchcomma ? s : s + 1 ;
		  s++;
		  while (*s != ',') s++;
		  s++;
	      } else {
		  s++;
		  while (*s != ',') s++;
		  s++;
	      }
	      break;
	  case 'd':		/* Has a default value: prev arg */
	      if (*analyseurs == ')') {
		  argvec[i] = argvec[i-1];
		  i++; s++;
	      }
	      break;
	  case 'P':		/* Argument in C is series precision */
	      argvec[i++] = (GEN) precdl;
	      break;
	  case 'f':		/* Fake *long argument */
	      argvec[i++] = (GEN) &fake;
	      break;
	  case 'F':		/* Fake *GEN argument */
	      argvec[i++] = (GEN) &fakeGEN;
	      fakepGEN = &fakeGEN;
	      break;
	  case 'x':		/* Foreign function */
	      argvec[i++] = (GEN) ep;
	      call = foreignHandler;
	      break;
	  case 'v':		/* Void return */
	      ret = RET_VOID;
	      break;
	  case ',':		/* Clean up default */
	      if (oldanalyseurs) {
		  analyseurs = oldanalyseurs;
		  oldanalyseurs = NULL;
	      }
	      break;
	  }
      }
      switch (ret) {
      case RET_INT: {
	  long resint;
      
#ifdef __cplusplus
	  long (*f1)(...) = (long (*)(...))call;
#else
	  long (*f1)() = (long (*)())call;
#endif
	  resint = (*f1)(argvec[0], argvec[1], argvec[2], argvec[3], argvec[4], 
			 argvec[5], argvec[6], argvec[7], argvec[8], argvec[9]);
	  analyseurtetpil=avma;
	  res = stoi(resint);
	  break;
      }
      case RET_VOID: {
#ifdef __cplusplus
	  void (*f1)(...) = (void (*)(...))call;
#else
	  void (*f1)() = (void (*)())call;
#endif
	  (*f1)(argvec[0], argvec[1], argvec[2], argvec[3], argvec[4], 
		argvec[5], argvec[6], argvec[7], argvec[8], argvec[9]);
	  analyseurtetpil=avma;
	  res = gnil;
	  break;
      }
      case RET_GEN:
	  analyseurtetpil=avma;
#ifdef __cplusplus
	  f = (GEN (*)(...))call;
#else
	  f = (GEN (*)())call;
#endif
	  res = (*f)(argvec[0], argvec[1], argvec[2], argvec[3], argvec[4], 
		     argvec[5], argvec[6], argvec[7], argvec[8], argvec[9]);
	  break;
      }
      if (fakepGEN) cgiv(fakeGEN);
      match(')');
      return res;
  }
  if (EpVALENCE(ep) < 100) /* fonctions predefinies */
  {
#ifdef __cplusplus
    f = (GEN (*)(...))ep->value;
#else
    f = (GEN (*)())ep->value;
#endif
    if (!EpVALENCE(ep) && (*analyseurs != '(')) return (*f)(prec);
    match('(');
    switch(EpVALENCE(ep))
    {
/*       case 0: res=(*f)(prec);break; */
/*       case 1: arg=expr();analyseurtetpil=avma; */
/* 	res=(*f)(arg,prec);break; */
/*       case 2: arg=expr();match(',');arg1=expr(); */
/* 	analyseurtetpil=avma;res=(*f)(arg,arg1,prec);break; */
/*       case 3: arg=expr();match(',');arg1=expr(); */
/* 	match(',');arg2=expr();analyseurtetpil=avma; */
/* 	res=(*f)(arg,arg1,arg2,prec);break; */
/*       case 4: arg=expr();match(',');arg1=expr(); */
/* 	match(',');arg2=expr();match(',');arg3=expr();analyseurtetpil=avma; */
/* 	res=(*f)(arg,arg1,arg2,arg3,prec);break; */
/*       case 5: arg=expr();match(',');arg1=expr();match(',');arg2=expr(); */
/* 	match(',');arg3=expr();match(',');arg4=expr();analyseurtetpil=avma; */
/* 	res=(*f)(arg,arg1,arg2,arg3,arg4,prec);break; */
/*       case 10: */
/* #ifdef __cplusplus	 */
/* 	p1=(GEN)(*(long(*)(...))f)(expr()); */
/* #else */
/* 	p1=(GEN)(*(long(*)())f)(expr()); */
/* #endif */
/* 	analyseurtetpil=avma;res=stoi((long)p1);break; */
/*       case 11: arg=expr();if(typ(arg)!=1) err(caseer); */
/* 	analyseurtetpil=avma;res=(*f)(itos(arg),prec);break; */
/*       case 12: arg=expr();match(',');arg1=expr();va=numvar(arg1); */
/* 	analyseurtetpil=avma;res=(*f)(arg,va,precdl);break; */
/*       case 13: arg=expr();analyseurtetpil=avma; */
/* 	res=(*f)(arg,&e);break; */
/*       case 14: arg=expr();match(',');arg1=expr();va=numvar(arg1); */
/* 	analyseurtetpil=avma;res=(*f)(arg,va);break; */
/*       case 15: arg=expr();if(typ(arg)!=1) err(caseer); */
/* 	analyseurtetpil=avma; */
/* #ifdef __cplusplus */
/* 	res=stoi((*(long(*)(...))f)(itos(arg))); */
/* #else	 */
/* 	res=stoi((*(long(*)())f)(itos(arg))); */
/* #endif	 */
/* 	break; */
/*       case 16: match('"'); */
/* 	{   char buffer[80], *s=buffer; */
/* 	    long r; */

/* 	    while ((*analyseurs)&&(*analyseurs!='"') && ((s - buffer)<79))  */
/* 		*s++ = *analyseurs++; */
/* 	    match('"'); */
/* 	    *s = 0; */
/* #ifdef __cplusplus	 */
/* 	    r=(*(long(*)(...))f)(buffer); */
/* #else */
/* 	    r=(*(long(*)())f)(buffer); */
/* #endif */
/* 	    res=stoi(r); */
/* 	    break; */
/* 	} */
/*       case 19:  */
/* 	arg=expr(); */
/* 	if(typ(arg)!=1) err(caseer); */
/* 	match(','); */
/* 	arg1=expr(); */
/* 	if(typ(arg1)!=1) err(caseer); */
/* 	analyseurtetpil=avma; */
/* 	res=f(itos(arg),itos(arg1)); */
/* 	break; */
/*       case 20: arg=expr();match(','); */
/* #ifdef __cplusplus */
/* 	p1=(GEN)(*(long(*)(...))f)(arg,expr()); */
/* #else */
/* 	p1=(GEN)(*(long(*)())f)(arg,expr()); */
/* #endif	 */
/* 	analyseurtetpil=avma;res=stoi((long)p1);break; */
/*       case 21: arg=expr();match(',');arg1=expr();if(typ(arg1)!=1) err(caseer); */
/* 	analyseurtetpil=avma;res=(*f)(arg,itos(arg1));break; */
/*       case 22: arg=expr();match(','); */
/* 	if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = findentry();if (EpVALENCE(ep)!=200) err(varer1,analyseurs);match(','); */
/* 	analyseurtetpil=avma; */
/* 	res=(*f)(ep,arg,analyseurs); skipexpr(); break; */
/*       case 23: arg=expr();match(',');arg1=expr();if(typ(arg1)!=1) err(caseer); */
/* 	analyseurtetpil=avma;res=(*f)(arg,itos(arg1),prec);break; */
/*       case 24: arg=expr();if(typ(arg)!=1) err(caseer); */
/* 	match(',');arg1=expr();analyseurtetpil=avma; */
/* 	res=(*f)(itos(arg),arg1,prec);break; */
/*       case 25: arg=expr();match(',');arg1=expr();analyseurtetpil=avma; */
/* 	res=(*f)(arg,arg1);break; */
/*       case 26: arg=expr();match(',');arg1=expr(); */
/* 	va=numvar(arg1);match(',');arg2=expr(); */
/* 	analyseurtetpil=avma;res=(*f)(arg,va,arg2);break; */
/*       case 27: if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = findentry();match('='); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	arg=expr(); match(','); analyseurtetpil=avma; */
/* 	res=(*f)(ep,arg,analyseurs,prec); skipexpr(); break; */
/*       case 28: arg=expr();match(',');arg1=expr(); */
/* 	analyseurtetpil=avma;res=(*f)(arg,arg1,&e);break; */
/*       case 29: arg=expr();match(',');arg1=expr();if(typ(arg1)!=1) err(caseer); */
/* #ifdef __cplusplus */
/* 	p1=(GEN)(*(long(*)(...))f)(arg,itos(arg1)); */
/* #else */
/* 	p1=(GEN)(*(long(*)())f)(arg,itos(arg1)); */
/* #endif	 */
/* 	analyseurtetpil=avma;res=stoi((long)p1);break; */
/*       case 30: arg=expr();match(',');arg1=expr();match(','); */
/* #ifdef __cplusplus */
/* 	p1=(GEN)(*(long(*)(...))f)(arg,arg1,expr()); */
/* #else */
/* 	p1=(GEN)(*(long(*)())f)(arg,arg1,expr()); */
/* #endif	 */
/* 	analyseurtetpil=avma;res=stoi((long)p1);break; */
/*       case 31: arg=expr();match(',');arg1=expr();match(','); */
/* 	analyseurtetpil=avma;res=(*f)(arg,arg1,expr(),&arg2);cgiv(arg2); */
/* 	break; */
/*       case 32: arg=expr();match(',');arg1=expr();match(',');arg2=expr(); */
/* 	if(typ(arg2)!=1) err(caseer); */
/* 	analyseurtetpil=avma;res=(*f)(arg,arg1,itos(arg2));break; */
/*       case 33: arg=expr();match(',');arg1=expr();match(',');arg2=expr(); */
/* 	if((typ(arg2)!=1)||(typ(arg1)!=1)) err(caseer); */
/* 	analyseurtetpil=avma;res=(*f)(arg,itos(arg1),itos(arg2),prec);break; */
/*       case 34: arg1=expr();match(',');arg2=expr();match(',');arg3=expr(); */
/* 	analyseurtetpil=avma; */
/* 	res=(*f)(gtolong(arg1),gtolong(arg2),gtolong(arg3)); */
/* 	break; */
/*       case 35: arg=expr();match(',');arg1=expr();match(',');arg2=expr(); */
/* 	if(typ(arg)!=1) err(caseer); */
/* 	analyseurtetpil=avma;res=(*f)(itos(arg),arg1,arg2);break; */
/*       case 37: if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = findentry();match('='); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	arg=expr(); match(','); arg1=expr(); match(',');  */
/* 	analyseurtetpil=avma; */
/* 	res=(*f)(ep,arg,arg1,analyseurs,prec); skipexpr(); break; */
/*       case 44: arg=expr();match(','); */
/* 	if(typ(arg)!=1) err(caseer); */
/* 	arg1=expr();match(','); */
/* 	if(typ(arg1)!=1) err(caseer); */
/* 	arg2=expr();match(','); */
/* 	if(typ(arg2)!=1) err(caseer); */
/* 	arg3=expr(); */
/* 	if(typ(arg3)!=1) err(caseer); */
/* 	analyseurtetpil=avma; */
/* 	res=(*f)(itos(arg),itos(arg1),itos(arg2),itos(arg3)); */
/* 	break; */
/*       case 45: arg=expr();match(','); */
/* 	if(typ(arg)!=1) err(caseer); */
/* 	arg1=expr();match(',');arg2=expr();match(',');arg3=expr(); */
/* 	if(typ(arg3)!=1) err(caseer); */
/* 	analyseurtetpil=avma; */
/* 	res=(*f)(itos(arg),arg1,arg2,itos(arg3),prec); */
/* 	break; */
/*       case 48: arg=expr(); match(','); */
/* 	if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = findentry();match('='); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	arg1=expr(); match(','); arg2=expr(); match(','); */
/* 	analyseurtetpil=avma; */
/* 	res=(*f)(ep,arg,arg1,arg2,analyseurs,prec); skipexpr(); break; */
/*       case 49: arg=expr();match(','); arg1=expr(); match(','); */
/* 	if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = findentry();if (EpVALENCE(ep)!=200) err(varer1,analyseurs);match(','); */
/* 	if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep1 = findentry();if (EpVALENCE(ep1)!=200) err(varer1,analyseurs);match(','); */
/* 	analyseurtetpil=avma; */
/* 	res=(*f)(ep,ep1,arg,arg1,analyseurs); skipexpr(); break; */
      case 50: p1=truc();
	if (*analyseurs++=='^') 
	{
	  arg=facteur();if(typ(arg)!=1) err(caseer);
	  e=itos(arg);
	}
	else {e = 1; analyseurs--;}
	analyseurtetpil=avma; res = ggrando(p1,e); break;
      case 51: case 52: case 53: case 54: case 55:
	if (*analyseurs != ')') for(;;)
	{
	  if (*analyseurs == '"')
	  {
	    analyseurs++;
	    while ((*analyseurs)&&(*analyseurs!='"')) pariputc(*analyseurs++);
	    match('"');
	  }
	  else
	  {
	    analyseurtetpil=avma;res=expr();
	    switch(EpVALENCE(ep))
	    {
	      case 51: case 53:
		brute(res,(char)glbfmt[0],glbfmt[2]);break;
	      case 52: case 54:
		sor(res,(char)glbfmt[0],glbfmt[2],glbfmt[1]);break;
	      case 55: texe(res,(char)glbfmt[0],glbfmt[2]);break;
	    }
	  }
	  if (*analyseurs == ')') break;
	  match(',');
	}
	if (EpVALENCE(ep)>52) pariputc('\n'); 
	fflush(outfile); if (logfile) fflush(logfile); break;
      case 56: 
	readbuffer = (char *)newbloc(paribuffsize>>TWOPOTBYTES_IN_LONG);
	while(!fgets(readbuffer, paribuffsize, infile)) switchin(NULL);
	if (pariecho) pariputs(readbuffer);
	else if (logfile) fputs(readbuffer, logfile);
	res=lisseq(readbuffer);killbloc((GEN)readbuffer);break;
      case 57:
	arg=expr();if(typ(arg)!=1) err(caseer);
	match(',');
	if(*analyseurs!='"') arg1=expr();
	else
	{
	  match('"');ch1=(char*)malloc(256);m=0;
	  while ((m<256)&&(*analyseurs)&&(*analyseurs!='"'))
	    ch1[m++]=*analyseurs++;
	  match('"');arg1=cgetg(m+1,17);
	  for(i=1;i<=m;i++) arg1[i]=lstoi((long)ch1[i-1]);
	  free(ch1);
	}
	analyseurtetpil=avma;res=(*f)(itos(arg),arg1);
	break;
/*       case 59: arg=expr();match(',');arg1=expr();match(','); */
/* 	arg2=expr();match(',');arg3=expr();match(',');arg4=expr(); */
/* 	if(typ(arg)!=1) err(caseer); */
/* 	analyseurtetpil=avma;res=(*f)(itos(arg),arg1,arg2,arg3,arg4); */
/* 	break;	   */
      case 60: arg=expr();if(typ(arg)!=1) err(caseer);
	m=itos(arg);if((m>=100)||(m<0)) err(labeler);
	labellist[m]=analyseurs;break;
      case 61: arg=expr();if(typ(arg)!=1) err(caseer);
	m=itos(arg);if((m>=100)||(m<0)||(!labellist[m])) err(labeler);
	analyseurs=labellist[m];yatileugoto=1;break;
/*       case 73:  */
/* 	arg=expr();if(typ(arg)!=1) err(caseer); */
/* 	match(',');  */
/* 	if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = findentry(); */
/* 	if (EpVALENCE(ep)!=200) err(varer1,analyseurs);match('='); */
/* 	arg1=expr(); match(','); */
/* 	arg2=expr(); match(','); */
/* 	arg3=expr(); if(typ(arg3)!=1) err(caseer); match(',');  */
/* 	arg4=expr(); if(typ(arg4)!=1) err(caseer); match(','); */
/* 	analyseurtetpil=avma; */
/* 	res = (*f)(itos(arg),ep,arg1,arg2,analyseurs, */
/* 		   prec,itos(arg3),itos(arg4)); */
/* 	skipexpr(); */
/* 	break; */
      case 80: av = avma; c=gcmp0(expr()); analyseurtetpil = avma = av; match(',');
	if (c) {skipseq();match(',');res = seq();}
	else 
	{
	  yatileugoto=0;res = seq();
	  if(!yatileugoto) {match(',');skipseq();}
	}
	break;
      case 81: analyseurtetpil = av = avma; ch1 = analyseurs;
	while (!gcmp0(expr()))
	{
	  analyseurtetpil = avma = av; match(',');
	  yatileugoto=0;seq();
	  if(!yatileugoto) analyseurs = ch1;else break;
	}
	if(!yatileugoto) {match(','); skipseq();}
	break;
      case 82: av = avma; ch1 = analyseurs;
	skipexpr();
	do 
	{
	  analyseurtetpil = avma = av; match(','); 
	  yatileugoto=0;seq();
	  if(!yatileugoto) analyseurs = ch1;else break;
	}
	while (gcmp0(expr()));
	if(!yatileugoto) {match(','); skipseq();}
	break;
/*       case 83: if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = findentry();match('='); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	arg=expr(); match(','); arg1=expr(); match(',');  */
/* 	analyseurtetpil=avma; */
/* 	res=(*f)(ep,arg,arg1,analyseurs); skipseq(); break; */
/*       case 84: arg=expr();match(','); */
/* 	if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = findentry();if (EpVALENCE(ep)!=200) err(varer1,analyseurs);match(','); */
/* 	analyseurtetpil=avma; */
/* 	res=(*f)(ep,arg,analyseurs); skipseq(); break; */
/*       case 85: if(!isalpha(*analyseurs)) err(killer1); */
/* 	ep = findentry(); if (EpVALENCE(ep)<100) err(killer1); */
/* 	killvalue(ep); */
/* 	if (EpVALENCE(ep) == 200) res = (GEN)ep->value; */
/* 	else */
/* 	{ */
/* 	  for(i = 0; i<TBLSZ; i++) */
/* 	    if (hashtable[i] == ep) { */
/* 		hashtable[i] = ep->next; */
/* 		if (ep->help && !(ep->valence & EpSTATIC)) free(ep->help); */
/* 		free(ep); break; */
/* 	    } */
/* 	    else */
/* 	      for(ep1 = hashtable[i]; ep1; ep1 = ep1->next) */
/* 		if (ep1->next == ep) { */
/* 		    ep1->next = ep->next;  */
/* 		    if (ep->help && !(ep->valence & EpSTATIC)) free(ep->help); */
/* 		    free(ep); break; */
/* 		} */
/* 	} */
/* 	break; */
/*       case 86: if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = findentry();match('='); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	arg=expr();match(',');arg1=expr();match(',');arg2=expr();match(','); */
/* 	analyseurtetpil=avma; */
/* 	res=(*f)(ep,arg,arg1,arg2,analyseurs); skipseq(); break; */
/*       case 87: if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = findentry();match('='); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	arg=expr(); match(','); analyseurtetpil=avma; */
/* 	res=(*f)(ep,arg,analyseurs,prec); skipseq(); break; */
/*       case 89: case 90: case 91: case 92: case 94: case 95: */
/* 	i=0; */
/* 	do */
/* 	{ */
/* 	  if(i) match(','); */
/* 	  argvec[i++]=expr(); */
/* 	} */
/* 	while((i<=6)&&(*analyseurs!=')')); */
/* 	switch(i) */
/* 	{ */
/* 	  case 1: argvec[1]=dbltor(0.3); /* cbach */ 
/* 	  case 2: argvec[2]=argvec[1]; /* cbach2 */ 
/* 	  case 3: argvec[3]=stoi(5); /* nrelsup */ 
/* 	  case 4: argvec[4]=gun; /* gborne pour petite norme */ 
/* 	  case 5: argvec[5]=stoi(4); /* nombre de relations par ideal */ 
/* 	  case 6: argvec[6]=stoi(3); /* cardinal minimal de la sous-fb */ 
/* 	  default: break; */
/* 	} */
/* 	analyseurtetpil=avma; */
/* 	res=(*f)(argvec[0],argvec[1],argvec[2],argvec[3],argvec[4],itos(argvec[5]),itos(argvec[6]),(EpVALENCE(ep))-92,prec); */
/* 	break; */
/*       case 96: /* buchimag */ 
/* 	i=0; */
/* 	do */
/* 	{ */
/* 	  if(i) match(','); */
/* 	  argvec[i++]=expr(); */
/* 	} */
/* 	while((i<=3)&&(*analyseurs!=')')); */
/* 	switch(i) */
/* 	{ */
/* 	  case 1: argvec[1]=dbltor(0.1); */
/* 	  case 2: argvec[2]=argvec[1]; */
/* 	  case 3: argvec[3]=stoi(5); */
/* 	  default: break; */
/* 	} */
/* 	analyseurtetpil=avma; */
/* 	res=(*f)(argvec[0],argvec[1],argvec[2],argvec[3]); */
/* 	break; */
/*       case 97: /* buchreal */ 
/* 	i=0; */
/* 	do */
/* 	{ */
/* 	  if(i) match(','); */
/* 	  argvec[i++]=expr(); */
/* 	} */
/* 	while((i<=4)&&(*analyseurs!=')')); */
/* 	switch(i) */
/* 	{ */
/* 	  case 1: argvec[1]=gzero; */
/* 	  case 2: argvec[2]=dbltor(0.1); */
/* 	  case 3: argvec[3]=argvec[2]; */
/* 	  case 4: argvec[4]=stoi(5); */
/* 	  default: break; */
/* 	} */
/* 	analyseurtetpil=avma; */
/* 	res=(*f)(argvec[0],argvec[1],argvec[2],argvec[3],argvec[4],prec); */
/* 	break; */

      default: err(valencer1);
    }
    match(')');return res;
  }
  switch (EpVALENCE(ep))
  {
    case 200: /* variables */
      if((*analyseurs)=='[')
      {
	tx=typ(p1=(GEN)ep->value);
	if((tx<17)||(tx>19)) err(caracer1,analyseurs);
	analyseurs++;av2=avma;flcol=flrow=0;
	if(tx<19)
	{
	  arg=expr();if(typ(arg)!=1) err(caseer);
	  c=itos(arg);if((c<1)||(c>=lg(p1))) err(arrayer1);
	}
	else
	{
	  if(lg(p1)==1) err(arrayer1);
	  if(*analyseurs==',')
	  {
	    analyseurs++;arg=expr();if(typ(arg)!=1) err(caseer);
	    c=itos(arg);if((c<1)||(c>=lg(p1))) err(arrayer1);
	    flcol=1;
	  }
	  else
	  {
	    arg=expr();if(typ(arg)!=1) err(caseer);
	    e=itos(arg);if((e<1)||(e>=lg((GEN)p1[1]))) err(arrayer1);
	    match(',');
	    if(*analyseurs==']') flrow=1;
	    else
	    {
	      arg1=expr();if(typ(arg1)!=1) err(caseer);
	      c=itos(arg1);
	      if((c<1)||(c>=lg(p1))) err(arrayer1);
	    }
	  }
	}
	match(']'); avma=av2;
	if(((*analyseurs)=='=')&&(*(analyseurs+1)!='=')) 
	{
	  long good=0, diff;
	  GEN oldx;
	  analyseurs++;res=expr();
	  if((tx==19)&&(!flcol))
	  {
	    if(flrow)
	    {
	      int diff;
	      if((typ(res)!=17)||(lg(res)!=lg(p1))) err(caseer2);
	      oldx=(GEN)p1[c]; good=2;
	      for(c=1;c<lg(p1);c++)
	      {
		if((diff=mateltsize(p1,c,e)-taille((GEN)res[c]))<0)
		{good=0; break;}		
		if (diff>0) good=1;
	      }
	      if (good==2 && !compact_arrays) good=0;
	      if (!good)
		for(c=1;c<lg(p1);c++) ((GEN)p1[c])[e]=(long)res[c];
		  /* maybe lcopy(res[c]) instead ? */
	      else for(c=1;c<lg(p1);c++)
		brutcopy((GEN)res[c], (GEN)(((GEN)p1[c])[e]));
	    }
	    else
	    {
	      diff=mateltsize(p1,c,e)-taille(res);
	      if (diff==0 || (!compact_arrays && diff>0)) good=1;
	      if (good) brutcopy(res, (GEN)(((GEN)p1[c])[e]));
	      else ((GEN)p1[c])[e]=(long)res;
	    }
	  }
	  else 
	  {
	    if(flcol)
	    {
	      if((typ(res)!=18)||(lg(res)!=lg((GEN)p1[1]))) err(caseer2);
	      diff=matcolsize(p1,c)-(taille(res)-lg(res));
	      if (diff==0 || (!compact_arrays && diff>0)) good=1;
	      if (good)
	      {
		GEN first=(GEN)((GEN)p1[c])[1];
		brutcopy((GEN)res[1], first);
		for(c=1;c<lg(p1)-1;)
		{
		  first += taille((GEN)((res)[c++]));
		  brutcopy((GEN)res[c], first);
		  p1[c]=(long)first;
		}
	      }
	      else p1[c]=(long)res;
	    }
	    else
	    {
	      diff=aryeltsize(p1,c)-taille(res);
	      if (diff==0 || (!compact_arrays && diff>0)) good=1;
	      if (good) brutcopy(res, (GEN) p1[c]);
	      else p1[c]=(long)res;
	    }
	  }
	  if (!good) {changevalue(ep, p1);p1=(GEN)ep->value;}
	}
	analyseurtetpil=avma;
	if((tx<19)||flcol) return (GEN)p1[c];
	else
	{
	  if(flrow) 
	  {
	    res=cgetg(lg(p1),17);
	    for(c=1;c<lg(p1);c++) res[c]=((GEN)p1[c])[e];
/* maybe lcopy() instead */
	    return res;
	  }
	  else return (GEN)((GEN)p1[c])[e];
	}
      }
      if(((*analyseurs)=='=')&&(*(analyseurs+1)!='=')) 
      {
	analyseurs++;changevalue(ep, expr()); 
      }
      analyseurtetpil=avma;return (GEN)ep->value;
      
    case 100: /* fonctions utilisateur */
      ch1 = analyseurs;
      match('(');
      p = (entree **)ep->value;
      nparam = (long)*p++;
      arg1 = arg = cgetg(nparam+1, 17);
      for(i = 0; (i < nparam) && (*analyseurs != ')'); i++)
      {
	if (i) match(',');
	*++arg = (long)expr();
      }
      if ((*analyseurs==')') && ((analyseurs[1] != '=') || (analyseurs[2] == '=')))
      {
	analyseurs++;
	while(i++ < nparam) *++arg = zero;
	analyseurtetpil = avma;
	for(i=0; i<nparam; i++) newvalue(*p++,(GEN)(*++arg1));
	res = lisseq((char *)p);
	res = forcecopy(res);
	for(i = 0; i < nparam; i++) 
	  killvalue(*--p);
	return res;
      }
      while (*analyseurs == ',') {analyseurs++; skipexpr();}
      match(')');
      if ((*analyseurs != '=') || (analyseurs[1] == '=')) err(nparamer1);
      analyseurs = ch1;
      killbloc((GEN)ep->value);
      
    case 101: /* nouvelle fonction */
      
      match('(');
      ch1 = analyseurs;
      for(nparam = 0; *analyseurs != ')'; nparam++)
      {
	if (nparam) match(',');
	if (!isalpha(*analyseurs)) err(paramer1);
	if (EpVALENCE(skipentry()) != 200) err(paramer1);
      }
      match(')'); match('='); ch2 = analyseurs; skipseq(); 
      p = (entree **)newbloc(nparam + ((analyseurs - ch2)>>TWOPOTBYTES_IN_LONG) + 2);
      p[-1] = (entree *)ep->value;
      ep->value = (void *)p;
      *p++ = (entree *)nparam;
      ch2 = analyseurs; analyseurs = ch1;
      for(i = 0; i < nparam; i++)
      {
	if (i) match(',');
	*p++ = ep1 = findentry();
	if (EpVALENCE(ep1) != 200) err(paramer1);
      }      
      match(')'); match('=');
      strncpy((char *)p, analyseurs, ch2 - analyseurs);
      *((char *)p + (ch2 - analyseurs)) = 0;
      ep->valence = 100;
      ep->menu = 0;
      analyseurs = ch2;
      return gnil;
      
    default: err(valencer1);return gnil;
  }
}

static long
word(long *nb)
{
  int m=0;
  for(*nb = 0; (*nb < 9) && isdigit(*analyseurs); (*nb)++)
    m = 10 * m + *analyseurs++-'0';
  return m;
}

GEN
constante()
{
  static long pw10[] = {1, 10, 100, 1000, 10000, 100000,
			1000000, 10000000, 100000000, 1000000000};
  long l,m=0,n=0,plus=1,nb, av = avma, limite=(avma + bot)/2;
  GEN z,y;
  
  analyseurtetpil=avma;
  y = stoi(word(&nb));
  while (isdigit(*analyseurs))
  {
    m = word(&nb); y = mulsi(pw10[nb], y);
    analyseurtetpil = avma;
    y = addsi(m, y);
    if (avma < limite)
    {
      y = gerepile(av, analyseurtetpil, y);
      analyseurtetpil = av;
    }
  }
  if ((*analyseurs!='.')&&(*analyseurs!='e')&&(*analyseurs!='E')) return y;
  if (*analyseurs=='.') 
  {
    analyseurs++;
    while (isdigit(*analyseurs))
    {
      m = word(&nb); y = mulsi(pw10[nb], y);
      analyseurtetpil = avma;
      y = addsi(m, y);
      if (avma < limite)
      {
	y = gerepile(av, analyseurtetpil, y);
	analyseurtetpil = av;
      }
      n -= nb;
    }
  }
  l=lgef(y);if(l<prec) l=prec;
  analyseurtetpil=avma;
  z=cgetr(l);affir(y,z);
  if ((*analyseurs=='e') || (*analyseurs=='E'))
  {
    analyseurs++;
    if (((*analyseurs)=='+') || ((*analyseurs)=='-')) plus=(*analyseurs++=='+');
    m = word(&nb);
    if(isdigit(*analyseurs)) err(expter1);
    if (plus) n += m;else n -= m;
  }
  if (n)
  {
    affsr(10, y = cgetr(l));
    y = gpuigs(y, labs(n));
    analyseurtetpil=avma;
    z = n > 0 ?  mulrr(z, y) : divrr(z, y);
  }
  return z;
}

entree *
findentry(void)
{
  char *olds = analyseurs, *u, *v;
  long sv, n;
  GEN p1;
  entree *ep;
  
  for (n = 0; isalnum(*analyseurs); analyseurs++) n = n << 1 ^ *analyseurs;
  if (n < 0) n = -n; n %= TBLSZ;
  for(ep = hashtable[n]; ep; ep = ep->next)
  {
      for(u = ep->name, v = olds; (*u) && *u == *v; u++, v++);
      if (!*u && (v == analyseurs)) return ep;
  }
  if (foreignAutoload && *analyseurs == '(') { /* Try to autoload. */
      if ((*foreignAutoload)(olds, analyseurs - olds)) {
	  for(ep = hashtable[n]; ep; ep = ep->next)  {
	      for(u = ep->name, v = olds; (*u) && *u == *v; u++, v++);
	      if (!*u && (v == analyseurs)) return ep;
	  }
	  /* Here we got an error. */
	  err(talker, "foreignAutoload reported success, but did not install a function");
      }
  }
  sv = (*analyseurs == '(') ? 0 : 7*BYTES_IN_LONG;
  ep = (entree *)malloc(sizeof(entree) + sv + analyseurs - olds + 1);
  ep->name = (char *)ep + sizeof(entree) + sv;
  for (u = ep->name, v = olds; v < analyseurs;) *u++ = *v++; *u = 0;
  ep->value = (void *)((char *)ep + sizeof(entree));
  ep->next = hashtable[n];
  hashtable[n] = ep;
  p1 = (GEN)ep->value;
  if (*analyseurs == '(') ep->valence = 101;
  else {
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
  }
  ep->help = ep->code = NULL;
  return ep;
}



void
skipseq(void)
{
  for(;;)
  {
    while(separe(*analyseurs)) analyseurs++;
    if ((!*analyseurs) || (*analyseurs == ')') || (*analyseurs == ',')) return;
    skipexpr();
    if(!separe(*analyseurs)) return;
  }
}

void
skipexpr(void)
{
  long niveau=3,e1,e2,e3;
  
  e1=e2=e3=0;
  for(;;)
    switch(niveau)
    {
      case 3: e3=1;skipfacteur();
	switch(*analyseurs)
	{
	  case '*': 
	  case '/':
	  case '%': analyseurs++;break;
	  case '\\': analyseurs++;if((*analyseurs)=='/') analyseurs++;
	    break;
	  default: niveau--;
	}
	break;
      case 2:
	if(!e3) {niveau++;break;}
	e3=0;e2=1;
	switch(*analyseurs)
	{
	  case '+':
	  case '-': analyseurs++;niveau++;break;
	  default: niveau--;
	}
	break;
      case 1: 
	if(!e2) {niveau++;break;}
	e2=0;e1=1;
	switch(*analyseurs)
	{
	  case '<': analyseurs++;
	    switch(*analyseurs)
	    {
	      case '=':
	      case '>': analyseurs++;niveau++;break;
	      default : niveau++;break;
	    }
	    break;
	  case '>': analyseurs++;
	    if((*analyseurs)=='=') analyseurs++;
	    niveau++; break;
	  case '=': 
	  case '!': 
	    if((analyseurs[1])=='=') {analyseurs+=2;niveau++;}
	    break;
	  default: niveau--;
	}
	break;
      case 0: 
	if(!e1) {niveau++;break;}
	e1=0;
	switch(*analyseurs)
	{
	  case '&': analyseurs++;if(*analyseurs=='&') analyseurs++;niveau++;break;
	  case '|': analyseurs++;if(*analyseurs=='|') analyseurs++;niveau++;break;
	  default: return;
	}
    }
}

void
skipfacteur(void)
{
  if (((*analyseurs)=='+') || ((*analyseurs)=='-')) analyseurs++;
  skiptruc();
  for (;;) switch(*analyseurs)
  {
    case '^': analyseurs++;skipfacteur(); break;
    case '~': 
    case '_':
    case '\'': analyseurs++;break;
    case '[': 
      analyseurs++;
      if(*analyseurs == ',') {analyseurs++;skipexpr();}
      else
      {
	skipexpr();
	if(*analyseurs==',')
	{
	  analyseurs++;if(*analyseurs != ']') skipexpr();
	}
      }
      match(']');break;
    case '!': analyseurs++;if((*analyseurs)!='=') break; else analyseurs--;
    default: return;
  }
}

void
skiptruc(void)
{
  long n=0,p=0,m=1;
  
  if (isalpha(*analyseurs)) {skipidentifier(); return;}
  if (isdigit(*analyseurs) || (*analyseurs=='.')) {skipconstante(); return;}
  switch(*analyseurs++)
  {
    case '(': skipexpr();match(')');return;
    case '[': if (*analyseurs!=']') 
    {do {n++; skipexpr();} while (*analyseurs++==',');analyseurs--;}
    switch (*analyseurs++)
    {
      case ']': return;
      case ';': m=n;do {n++; skipexpr();} while (*analyseurs++!=']');
	if (n % m) err(recter1);
	return;
      default: err(vectmater1);
    }
    case '%':
      p=0;while((*analyseurs)=='`') {analyseurs++;p++;}
    if(p>tglobal) err(referer1);
    if(p) return;
    while (isdigit(*analyseurs)) p = 10*p + *analyseurs++ - '0';
    if(p>tglobal) err(referer1);
    return;
  }
  err(caracer1,analyseurs-1);
}

void
skipidentifier(void)
{
  long nparam, i, m;
  entree *ep, **p;
  char *ch1;
  GEN arg;

  ep = skipentry();
  if (ep->code) {
    int matchcomma = 0;
    char delim, *s = ep->code;
    entree *ep1;
    
    if (!*(s+1) && (*s == 'p') && (*analyseurs != '(')) return;
    match('(');
    /* Optimized for G and p. */
    while (*s == 'G') {
	if (matchcomma) {
	    match(',');
	} else matchcomma = 1;
	skipexpr();
	s++;
    }
    if (*s == 'p') s++;
    while (*s) {
	switch (*s++) {
	default:
	    err(talker,"Unknown type designator during skip");
	case 'G':				/* Argument is GEN */
	case 'n':				/* Argument is ordinal of var */
	case 'L':				/* Argument is long, XXX 64bit? */
	case 'Z':				/* Argument must be converted
						 * to long */
	case 'I':				/* Input argument */
	    if (matchcomma) {
		match(',');
	    } else matchcomma = 1;
	    skipexpr();
	    break;
	case 'V':				/* Argument is variable */
	case 'S':				/* Argument is a symbol */
	    if (matchcomma) {
		match(',');
	    } else matchcomma = 1;
	    if(!isalpha(*analyseurs)) err(varer1,analyseurs);
	    ep1 = skipentry();
	    if (EpVALENCE(ep1) != 200 && *(s-1) == 'V') 
		err(varer1,analyseurs);
	    break;
	case 's':				/* Argument is string */
	    if (matchcomma) {
		match(',');
	    } else matchcomma = 1;
	    delim = *analyseurs++;
	    if (delim != '"' && delim != '\'') {
		err(talker,"String should be delimited by \" or \'");
	    }
	    while ((*analyseurs) && (*analyseurs != delim))
		analyseurs++;
	    match(delim);
	    break;
	case 'p':
	case 'P':
	case 'l':
	case 'v':
	case 'f':
	case 'F':
	case '|':
	case '~':
	case 'x':
	    break;
	case 'D':
	    if (*analyseurs == ')') {
		match(')');
		return;
	    }
	    s++;
	    while (*s != ',') s++;
	    s++;		/* Comma. */
	    break;
	case 'd':
	    if (*analyseurs == ')') {
		match(')');
		return;
	    }
	    break;			/* Just visual separation. */
	case '=':
	    match('=');
	    matchcomma = 0;
	    break;
	case ',':
	    break;
	}
    }
    match(')');
    return;
  }
  if (EpVALENCE(ep) < 100) /* fonctions predefinies */
  {
    if (!EpVALENCE(ep) && (*analyseurs != '(')) return;
    match('(');
    switch(EpVALENCE(ep))
    {
/*       case 0: */
      case 56: 
	  break;
/*       case 1: */
/*       case 10: */
/*       case 11: */
/*       case 13: */
/*       case 15: */
      case 61: skipexpr(); break;
      case 60: arg=expr();if(typ(arg)!=1) err(caseer);
	m=itos(arg);if((m>=100)||(m<0)) err(labeler);
	labellist[m]=analyseurs;break;
      case 51: case 52: case 53: case 54: case 55: case 16:
	if (*analyseurs != ')') for(;;)
	{
	  if (*analyseurs == '"')
	  {
	    analyseurs++;
	    while ((*analyseurs)&&(*analyseurs!='"')) analyseurs++;
	    match('"');
	  }
	  else skipexpr();
	  if (*analyseurs == ')') break;
	  match(',');
	}
	break;
      case 57:
	skipexpr();match(',');if(*analyseurs!='"') skipexpr();
	else
	{
	  match('"');m=0;
	  while ((m<256)&&(*analyseurs)&&(*analyseurs!='"')) 
	  {m++;analyseurs++;}
	  match('"');
	}
	break;
/*       case 2: */
/*       case 12: */
/*       case 14: */
/*       case 20: */
/*       case 21: */
/*       case 23: */
/*       case 24: */
/*       case 25: */
/*       case 28: */
/*       case 19: */
/*       case 29: skipexpr(); match(','); skipexpr(); break; */
/*       case 22: skipexpr(); match(','); */
/* 	if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = skipentry(); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	match(','); skipexpr(); break; */
/*       case 27: if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = skipentry(); match('='); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	skipexpr(); match(','); skipexpr(); break; */
/*       case 3: */
/*       case 26:  */
/*       case 30: */
/*       case 31: */
/*       case 32: */
/*       case 33: */
/*       case 34: */
/*       case 35: skipexpr();match(',');skipexpr();match(',');skipexpr(); */
/* 	break; */
/*       case 37: if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = skipentry(); match('='); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	skipexpr(); match(','); skipexpr(); match(','); skipexpr();break; */
/*       case 4: skipexpr();match(',');skipexpr();match(',');skipexpr(); */
/* 	match(',');skipexpr();break; */
/*       case 48: skipexpr(); match(','); */
/* 	if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = skipentry(); match('='); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	skipexpr(); match(','); skipexpr(); match(','); skipexpr();break; */
/*       case 49: skipexpr(); match(','); skipexpr(); match(','); */
/* 	if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = skipentry(); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	match(',');if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = skipentry(); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	match(','); skipexpr(); break; */
/*       case 44:  */
/*       case 45:  */
/* 	skipexpr();match(','); */
/* 	skipexpr();match(','); */
/* 	skipexpr();match(','); */
/* 	skipexpr(); */
/* 	break; */
/*       case 73: skipexpr(); match(','); */
/* 	if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = skipentry(); match('='); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	skipexpr(); match(','); skipexpr(); match(',');  */
/* 	skipexpr(); match(','); skipexpr(); match(',');  */
/* 	skipexpr(); */
/* 	break; */
      case 50: skiptruc();
	if (*analyseurs++=='^') skipfacteur();else analyseurs--;
	break;
/*       case 5: */
/*       case 59: skipexpr();match(',');skipexpr();match(',');skipexpr(); */
/* 	match(',');skipexpr();match(',');skipexpr(); */
/* 	break; */
      case 80: skipexpr(); match(','); skipseq(); match(','); skipseq(); 
	  break;
      case 81:
      case 82: skipexpr(); match(','); skipseq(); 
	  break;
/*       case 83: if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = skipentry(); match('='); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	skipexpr(); match(','); skipexpr(); match(','); skipseq();  */
/* 	break; */
/*       case 84: skipexpr(); match(','); */
/* 	if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = skipentry(); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	match(','); skipseq();  */
/* 	break; */
/*       case 85: if(!isalpha(*analyseurs)) err(killer1); */
/* 	ep = skipentry(); if (EpVALENCE(ep)<100) err(killer1); */
/* 	break; */
/*       case 86: if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = skipentry(); match('='); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	skipexpr();match(',');skipexpr();match(',');skipexpr();match(',');skipseq(); */
/* 	break; */
/*       case 87: if(!isalpha(*analyseurs)) err(varer1,analyseurs); */
/* 	ep = skipentry(); match('='); if (EpVALENCE(ep)!=200) err(varer1,analyseurs); */
/* 	skipexpr(); match(','); skipseq();  */
/* 	break; */
/*       case 89: case 90: case 91: case 92: case 94: case 95: case 97: */
/* 	i=0; */
/* 	do */
/* 	{ */
/* 	  if(i) match(','); */
/* 	  i++;skipexpr(); */
/* 	} */
/* 	while((i<=6)&&(*analyseurs!=')')); */
/* 	break; */
/*       case 96: */
/* 	i=0; */
/* 	do */
/* 	{ */
/* 	  if(i) match(','); */
/* 	  i++;skipexpr(); */
/* 	} */
/* 	while((i<=3)&&(*analyseurs!=')')); */
/* 	break; */
      default: err(valencer1);
    }
    match(')');
    return;
  }
  switch (EpVALENCE(ep))
  {
    case 200: /* variables */
      if((*analyseurs)=='[')
      {
	analyseurs++;
	if(*analyseurs == ',') {analyseurs++;skipexpr();}
	else
	{
	  skipexpr();
	  if(*analyseurs == ',')
	  {
	    analyseurs++;if(*analyseurs != ']') skipexpr();
	  }
	}
	match(']');
      }
      if(((*analyseurs)=='=')&&(*(analyseurs+1)!='=')) 
      {
	analyseurs++;skipexpr(); 
      }
      return;
      
    case 100: /* fonctions utilisateur */
      ch1 = analyseurs;
      match('(');
      p = (entree **)ep->value;
      nparam = (long)*p++;
      i = 0;
      for(i = 0; (i < nparam) && (*analyseurs != ')'); i++)
      {
	if (i) match(',');
	skipexpr();
      }
      if ((*analyseurs==')') && ((analyseurs[1] != '=') || (analyseurs[2] == '='))) {analyseurs++; return;}
      while (*analyseurs == ',') {analyseurs++; skipexpr();}
      match(')');
      if ((*analyseurs != '=') || (analyseurs[1] == '=')) err(nparamer1);
      analyseurs = ch1;
      
    case 101: /* nouvelle fonction */
      
      match('(');
      for(nparam = 0; *analyseurs != ')'; nparam++)
      {
	if (nparam) match(',');
	skipexpr();
      };
      match(')');
      if (*analyseurs == '=') {analyseurs++; skipseq();}
      return;
      
    default: err(valencer1);
  }
}

void
skipconstante(void)
{  
  while (isdigit(*analyseurs)) analyseurs++;
  if ((*analyseurs!='.')&&(*analyseurs!='e')&&(*analyseurs!='E')) return;
  if (*analyseurs=='.') analyseurs++;
  while (isdigit(*analyseurs)) analyseurs++;
  if ((*analyseurs=='e') || (*analyseurs=='E'))
  {
    analyseurs++;
    if (((*analyseurs)=='+') || ((*analyseurs)=='-')) analyseurs++;
    while (isdigit(*analyseurs)) analyseurs++;
  }
}

entree fake101 = {"",101,0,0,0,NULL,NULL};
entree fake200 = {"",200,0,0,0,NULL,NULL};

entree *
skipentry(void)
{
  char *u, *v, *olds = analyseurs;
  long n;
  entree *ep;
  
  for(n = 0; isalnum(*analyseurs); analyseurs++) n = n << 1 ^ *analyseurs;
  if (n < 0) n = -n; n %= TBLSZ;
  for(ep = hashtable[n]; ep; ep = ep->next)
  {
    for(u = ep->name, v = olds; (*u) && *u == *v; u++, v++);
    if (!*u && (v == analyseurs)) return ep;
  }
  return (*analyseurs == '(') ? &fake101 : &fake200;
}
