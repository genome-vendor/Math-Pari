/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*            Algorithmes SOUS-EXPONENTIELS DE CALCUL              */
/*            DU GROUPE DE CLASSES ET DU REGULATEUR                */
/*                   CORPS DE NOMBRES GENERAUX                     */
/*                      (McCURLEY, BUCHMANN)                       */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "genpari.h"

long CBUCHG; /* de la forme 2^k-1 */

GEN buchall_for_degree_one_pol(long v, long flun);
GEN initalgall0(GEN x,long fldif,long prec);
GEN subfactorbasegen(long N,long m,GEN vectbase,long minsfb,long *vperm,long *ptss);
GEN **powsubfactgen(GEN nf,GEN w,long a,long PRECREG,long PRECREGINT);
GEN idealmulprimered(GEN nf,GEN x,GEN vp,long PRECREG);
GEN idealmulprimeredall(GEN nf,GEN x,GEN vp,long PRECREG,long PRECREGINT);
GEN getfu(GEN nf,GEN *ptxarch,GEN reg,long flun,long *pte,long PRECREG);
GEN cleancol(GEN x,long N,long RU,long PRECREG);
GEN quad_form(GEN *cbase,GEN ma,GEN ideal,long PRECREG);
GEN compute_regulator(GEN mit,long col,GEN matarch,long RU,long R1,long N,long PRECREG,GEN *xarch,GEN *xreal,GEN *image_mdet,long *ptsreg);
GEN compute_class_number(GEN mit,GEN *met,GEN *u1);
GEN compute_check(GEN mit,long col,GEN image_mdet,GEN xreal,long RU,GEN *reg,GEN z,GEN zu,GEN clh,long PRECREG,GEN *parch);
GEN class_group_generators(long N,long RU,long PRECREG,long PRECREGINT,GEN nf,GEN met,GEN clh,GEN u1,GEN vectbase,long *vperm);
long factorisegen(GEN nf,GEN ideal,long kcz,long limp,long *primfact,long *expoprimfact,long *primfactorbase,GEN *idealbase,long *numideal,long *numprimfactorbase);
GEN compute_matt2(long RU,GEN nf,long *nbmatt2);
long factorisealpha(GEN nf,GEN alpha,long kcz,long limp,long *primfact,long *expoprimfact,long *primfactorbase,GEN *idealbase,long *numideal,long *numprimfactorbase);
long factorbasegen(GEN nf,long n2,long n,long **ptnumprim,long **ptprim,long **ptnum,GEN **ptideal,long *ptkc,long *ptkcz,long *ptkcz2,GEN *ptlfun);
long small_norm_for_buchall(long s,long **mat,GEN matarch,long KC,long KCCO,long KCZ,long LIMC,long N,GEN D,long RU,long R1,long PRECREG,GEN P,GEN nf,GEN vectbase,long *prim,GEN *idealbase,long *numideal,long *numprim,GEN IDN,GEN gborne,long nbrelpid);
long random_relation(long t,long q,long LIMC,long KC,long KCZ,long RU,long PRECREG,long PRECREGINT,GEN nf,long *firstime,long *secondtime,long *iz,long *jideal,long *bouext,long *bou,long *slim,long **ptvinvperm,long lgsub,GEN *idealbase,long *numideal,long *prim,long *numprim,long *vperm,long *subvperm,long **mat,GEN matarch,long **matcopy,GEN lmatt2,long nbmatt2,GEN mit,long *ex,GEN IDN,GEN **vp,GEN vectbase,long nrelsup,long nlze,long *phase,long *extrarel,long *cmptglob,long *nbtest,GEN *extramat,GEN *extramatc);
long already_found_relation_of_long (long **mat,long s,long KC);
void desallocate(GEN primsubfactorbase,long cmptglob,GEN **vp,long *numprim,long *prim,long *numideal,GEN *idealbase,long *ex,long **matcopy,long *vperm,long *subvperm,long *vinvperm,long flvp,long phase);
long be_honest(long KCZ,long KCZ2,long RU,long PRECREGINT,GEN IDN,long s,GEN nf,long lgsub,GEN *idealbase,long *numideal,long *prim,long *numprim,long q,GEN **vp,long *ex,GEN exu);
GEN compute_prhall(GEN nf,GEN pr);

#define MAXRELSUP (min(50,4*KC))

GEN
buchall(GEN P,GEN gcbach,GEN gcbach2,GEN gRELSUP,GEN gborne,long nbrelpid,long minsfb,long flun,long prec)
{
  long av,av0,av1,av2,tetpil,dec,i,j,k,ss,lgsub,n1,*p1,*ex,q,s,nbtest,resuhonest,c;
  long N,R1,R2,RU,CO,PRECREG,PRECREGINT,LIMC,LIMC2,lim,KC,KC2,KCZ,KCZ2,KCCO,KCCOPRO;
  long ip,extrarel,iz,ncz,col,k0,nlze,sreg,colnew,bou,bouext,cmptglob,limpile;
  long **mat,**matcopy,*vperm,*vinvperm,*subvperm;
  long *numprim,*prim,*numideal,nrelsup,nreldep,phase,sfact,slim,jideal;
  long firstime,secondtime,nbmatt2,flT2,flvp,flHONEST;
  double cbach,cbach2,drc,LOGD,rootdisc;
  GEN pgen1,pgen2,pgen3,matgenpro,lmatt2,**vp,*idealbase,matarch,exu,fu,zu;
  GEN nf,D,F,BINV,primsubfactorbase,vectbase,xarch,xreal,met,mit,reg,lfun,z,clh;
  GEN matalpha,u1,RES,extramat,extramatc,c_1,image_mdet,pdep,parch,IDN;

  if(DEBUGLEVEL) timer2();
  if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  N=lgef(P)-3;if(N<=1) return buchall_for_degree_one_pol(varn(P),flun);
  av=avma;pgen1=factor(P);pgen2=(GEN)pgen1[1];pgen3=(GEN)pgen1[2];
  if((lg(pgen2)>2)||(cmpis((GEN)pgen3[1],1)>0)) err(bucher10);avma=av;
  av=avma;
  if(DEBUGLEVEL)
  {
    fprintferr("temps trivialites: ");fprintferr("%ld\n",timer2());
    flusherr();
  }
FORCEFU:

/* Initialisations */  

  if(typ(gRELSUP)!=1) gRELSUP=gtrunc(gRELSUP);if(gsigne(gRELSUP)<=0) err(bucher5);
  IDN=idmat(N);nf=initalgall0(P,flun>=0?0:1,max(8,prec));
  if(DEBUGLEVEL)
  {
    fprintferr("temps initalg: ");fprintferr("%ld\n",timer2());
    flusherr();
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  }
  if(abs(flun)>1) RES=cgetg(11,17);else RES=flun?cgetg(9,17):cgetg(8,17);
  RES[1]=nf[1];RES[2]=nf[2];pgen1=cgetg(3,17);RES[3]=(long)pgen1;pgen1[1]=nf[3];
  pgen1[2]=nf[4];RES[4]=nf[7];
  R1=itos((GEN)((GEN)nf[2])[1]);R2=(N-R1)/2;RU=R1+R2;F=(GEN)nf[4];D=(GEN)nf[3];
  BINV=(GEN)nf[8];zu=rootsof1(nf);CO=itos(gRELSUP)+RU-1;
  drc=gtodouble(D);LOGD=log(fabs(drc));
  lim=max((long)(exp(-(N+0.))*sqrt(2*PI*N*fabs(drc))*pow(4/PI,R2+0.)),3);
  cbach=gtodouble(gcbach);cbach2=gtodouble(gcbach2);rootdisc=exp(LOGD/N);
  if(DEBUGLEVEL)
  {
    fprintferr("temps rootsof1 et/ou trivialites: ");
    fprintferr("%ld\n",timer2());flusherr();
  }
  av0=avma;
INCREASEGEN:
  if(DEBUGLEVEL) fprintferr("cbach = %f\n",cbach);
  nreldep=nrelsup=0;LIMC=(long)(cbach*LOGD*LOGD);LIMC2=max(3*N,(long)(cbach2*LOGD*LOGD));
  if(LIMC>LIMC2) LIMC2=LIMC;

/* Calcul de la factor-base, de la sub-factor-base et initialisation de vperm */  

  if(DEBUGLEVEL)
  {
    fprintferr("\n N = %ld , R1 = %ld , R2 = %ld , RU = %ld",N,R1,R2,RU);
    fprintferr("\n D^(1/N) = %f",rootdisc);
    fprintferr(" , LIMC = %ld , LIMC2 = %ld\n",LIMC,LIMC2);
    fprintferr(" D = ");outerr(D);flusherr();
  }
  KC2=factorbasegen(nf,LIMC2,LIMC,&numprim,&prim,&numideal,&idealbase,&KC,&KCZ,&KCZ2,&lfun);
  if(!KC)
  {
    free(numprim);free(prim);free(numideal);free(idealbase);
    avma=av0;cbach=2*cbach;cbach2=max(cbach2,cbach);goto INCREASEGEN;
  }

  if(DEBUGLEVEL)
  {
    fprintferr("temps creation de la factor base: ");
    fprintferr("%ld\n",timer2());flusherr();
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  }
  if(DEBUGLEVEL>6)
  {for(i=1;i<=KCZ;i++){fprintferr("++ idealbase[%ld] = ",i);outerr(idealbase[i]);}}
  if(DEBUGLEVEL>7)
  {
    fprintferr("########## FACTORBASE ##########\n\n");
    fprintferr("KC2 = %ld , KC = %ld , KCZ = %ld , KCZ2 = %ld\n",KC2,KC,KCZ,KCZ2);
    for(i=1;i<=KCZ;i++){fprintferr("++ idealbase[%ld] = ",i);outerr(idealbase[i]);}
  }
  vectbase=cgetg(KC+1,18);
  for(i=1;i<=KCZ;i++) 
  {
    ip=numideal[prim[i]];pgen1=idealbase[i];n1=lg(pgen1);
    for(j=1;j<n1;j++) vectbase[ip+j]=pgen1[j];
  }
  z=gmul(gmul(gmul(gpuigs(gdeux,-R1),gpuigs(gmul2n(mppi(DEFAULTPREC),1),-R2)),gsqrt(gabs(D,0),DEFAULTPREC)),lfun);
  vperm=(long*)malloc(sizeof(long)*(KC+1));for(i=1;i<=KC;i++) vperm[i]=i;
  primsubfactorbase=subfactorbasegen(N,min(lim,LIMC2),vectbase,minsfb,vperm,&ss);
  lgsub=lg(primsubfactorbase);
  if(DEBUGLEVEL)
  {
    fprintferr("temps creation de la sous-factor base: ");
    fprintferr("%ld\n",timer2());flusherr();
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  }
  if(DEBUGLEVEL>3)
  {
    fprintferr("\n***** LISTE DES IDEAUX DE LA FACTEUR BASE *****\n\n");
    for(i=1;i<=KC;i++){fprintferr("no %ld = ",i);outerr((GEN)vectbase[i]);}
    fprintferr("\n***** LISTE DES IDEAUX DE LA SOUS-FACTEUR BASE *****\n\n");
    outerr(primsubfactorbase);
  }
  if(DEBUGLEVEL>3)
  {
    fprintferr("\n***** PERMUTATION D'ORIGINE *****\n\n");
    for(i=1;i<=KC;i++) fprintferr("vperm[%ld] = %ld\n",i,vperm[i]);
    fprintferr("\n");flusherr();
  }
  if(ss==-1) 
  {
    free(vperm);free(numprim);free(prim);free(numideal);free(idealbase);
    avma=av0;cbach=2*cbach;cbach2=max(cbach2,cbach);
    goto INCREASEGEN;
  }

/* Calcul de la precision, creation de matarch, affectation des relations triviales */  

  PRECREGINT=((gexpo(D)*(lgsub-2)+((N*N)>>2))>>TWOPOTBITS_IN_LONG)+4;
  PRECREG=max(prec+1,PRECREGINT);CBUCHG=15;KCCO=KC+RU-1+max(ss,itos(gRELSUP));
  if(DEBUGLEVEL) fprintferr("sous-factor base a %ld elements\n",lgsub-1);
  if(DEBUGLEVEL)
  {
    fprintferr("\nCO = %ld , nbrelsup = %ld , ss = %ld , ",CO,itos(gRELSUP),ss);
    fprintferr("KCZ = %ld , KC = %ld , KCCO = %ld \n",KCZ,KC,KCCO);
    flusherr();
  }
  subvperm=(long*)malloc(sizeof(long)*lgsub);for(i=1;i<lgsub;i++) subvperm[i]=vperm[i];
  ex=(long*)malloc(lgsub<<TWOPOTBYTES_IN_LONG);exu=cgetg(RU+1,17);
  q=BITS_IN_RANDOM-1-(long)ceil(log((double)CBUCHG)/log(2.0));s=0;
  mat=(long**)malloc(sizeof(long*)*(KCCO+1));matarch=cgetg(KCCO+1,19);
  for(i=1;i<=KCCO;i++)
  {
    p1=(long *)malloc(sizeof(long)*(KC+1));mat[i]=p1;
    for(j=1;j<=KC;j++) p1[j]=0;pgen1=cgetg(RU+1,18);matarch[i]=(long)pgen1;
    for(j=1;j<=RU;j++) 
    {
      pgen2=cgetg(3,6);pgen1[j]=(long)pgen2;
      pgen2[1]=lgetr(PRECREG);pgen2[2]=lgetr(PRECREG);
    }
  }
  for(i=1;i<=KCZ;i++)
  {
    ip=numideal[prim[i]];pgen1=idealbase[i];n1=lg(pgen1);k=0;
    for(j=1;j<n1;j++) k+=(itos((GEN)((GEN)pgen1[j])[3])*itos((GEN)((GEN)pgen1[j])[4]));
    if(k==N)
    {
      s++;for(j=1;j<n1;j++) mat[s][ip+j]=itos((GEN)((GEN)pgen1[j])[3]);
      for(j=1;j<=RU;j++) gaffsg(0,(GEN)((GEN)matarch[s])[j]);
    }
  }
  nbtest=ss=s;
  if(DEBUGLEVEL)
  {
    fprintferr("temps creation mat et matarch: ");
    fprintferr("%ld\n",timer2());flusherr();
    fprintferr("Apres relations triviales, s = nbtest = %ld\n",s);flusherr();
  }

/* Recherche de relations par elements de petite norme */

  if(DEBUGLEVEL)
  {
    fprintferr("\n#### Recherche des relations avec elements de petite norme\n");
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  }
  if(gsigne(gborne)>0)
    s=small_norm_for_buchall(s,mat,matarch,KC,KCCO,KCZ,LIMC,N,D,RU,R1,PRECREG,P,nf,vectbase,prim,idealbase,numideal,numprim,IDN,gborne,nbrelpid);

/* Recherche de relations aleatoires */
 
  flT2=flvp=flHONEST=0;
  nlze=slim=KCCO;firstime=secondtime=cmptglob=phase=0;
  iz=jideal=bouext=bou=1;
  limpile=(avma+bot)>>1;
  if(s==KCCO) goto HNFSPEC;
  if(DEBUGLEVEL)
  {
    fprintferr("\n#### Recherche de relations aleatoires\n");flusherr();
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  }
LABELINT:
  if(!flT2){flT2=1;lmatt2=compute_matt2(RU,nf,&nbmatt2);}
LABELVP:
  if(!flvp)
  {
    if(DEBUGLEVEL)
    {
      fprintferr("Calcul des puissances de la sub-factor base :\n");
      if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
      flusherr();
    }
    vp=powsubfactgen(nf,primsubfactorbase,CBUCHG+1,PRECREG,PRECREGINT);flvp=1;
    if(DEBUGLEVEL)
    {
      fprintferr("temps powsubfactgen: ");fprintferr("%ld\n",timer2());
      flusherr();
      if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
    }
    if(DEBUGLEVEL>7)
    {
      fprintferr("\n***** LISTE DES PUISSANCES DE LA SOUS-FACTEUR BASE *****\n\n");
      for(i=1;i<lgsub;i++)
      {
	for(j=0;j<=CBUCHG+1;j++){fprintferr("vp[%ld][%ld] = ",i,j);outerr(vp[i][j]);}
	fprintferr("\n");
      }
    }
  }
  if(flHONEST) goto LABELHONEST;
  av1=avma;
  s=random_relation(s,q,LIMC,KC,KCZ,RU,PRECREG,PRECREGINT,nf,&firstime,&secondtime,&iz,&jideal,&bouext,&bou,&slim,&vinvperm,lgsub,idealbase,numideal,prim,numprim,vperm,subvperm,mat,matarch,matcopy,lmatt2,nbmatt2,mit,ex,IDN,vp,vectbase,nrelsup,nlze,&phase,&extrarel,&cmptglob,&nbtest,&extramat,&extramatc);
  if(s==-1)
  {
    for(j=1;j<=KCCO;j++) free(mat[j]);free(mat);
    desallocate(primsubfactorbase,cmptglob,vp,numprim,prim,numideal,idealbase,ex,matcopy,vperm,subvperm,vinvperm,flvp,phase);
    if(cbach>11.99) 
      err(talker,"sorry, buchxxx is not able to compute this field PLEASE REPORT!!!");
    else
    {
      avma=av0;cbach=min(2*cbach,12);cbach2=max(cbach2,cbach);
      goto INCREASEGEN;
    }
  }
HNFSPEC:

/* Reductions des matrices des relations */  

  if(!phase)
  {
    if(DEBUGLEVEL>2)
    {
      fprintferr("A la sortie : iz=%ld,",iz);
      fprintferr("jideal=%ld,bouext=%ld,bou=%ld,",jideal,bouext,bou);
      fprintferr("phase=%ld\n",phase);flusherr();
      if(DEBUGLEVEL>3){fprintferr("***** AVMA = %ld\n",avma);flusherr();}
      av2=avma;matgenpro=cgetg(s+1,19);
      for(j=1;j<=s;j++)
      {
	pgen1=cgetg(KC+1,18);matgenpro[j]=(long)pgen1;
	for(i=1;i<=KC;i++) pgen1[i]=lstoi(mat[j][i]);
      }
      fprintferr("\nRang  = %ld, temps = %ld\n",rank(matgenpro),timer2());
      flusherr();avma=av2;
      if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
    }
    k0=lgsub-1;
    if(DEBUGLEVEL>3)
    {
      fprintferr("mat = \n");
      for(i=1;i<=KC;i++)
      {
	fprintferr("[");
	for(j=1;j<=KCCO;j++) fprintferr("%4ld",mat[j][i]);fprintferr("]\n");
      }
      fprintferr("\nmatarch = ");outerr(matarch);
    }
    matcopy=(long**)malloc(sizeof(long*)*(10*KCCO+1+MAXRELSUP));
    for(j=1;j<=KCCO;j++)
    {
      p1=(long *)malloc(sizeof(long)*(KC+1));matcopy[j]=p1;
      for(i=1;i<=KC;i++) p1[i]=mat[j][i];
    }
    mit=hnfspec(mat,&pdep,&matarch,vperm,&matalpha,KCCO,KC,k0,&nlze,&col);
    if(DEBUGLEVEL){fprintferr("temps hnfspec = %ld\n",timer2());flusherr();}
    for(j=1;j<=KCCO;j++) free(mat[j]);free(mat);
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
    cmptglob=KCCO;KCCOPRO=KCCO;sfact=s;secondtime=0;
    phase=1;firstime=1;
    if(nlze)
    {
      if(avma<limpile)
      {
	tetpil=avma;mit=gcopy(mit);
	pdep=gcopy(pdep);matarch=gcopy(matarch);matalpha=gcopy(matalpha);
	dec=lpile(av1,tetpil,0)>>TWOPOTBYTES_IN_LONG;
	mit+=dec;matalpha+=dec;matarch+=dec;pdep+=dec;
      }
      goto LABELINT;
    }
  }
  else
  {
    if(nrelsup) nlze=0;
    if(DEBUGLEVEL>6) 
    {
      fprintferr("avant hnfadd :\n");fprintferr("vectbase[vperm[]] = \n");
      fprintferr("[");
      for(i=1;i<=KC;i++)
      {bruterr((GEN)vectbase[vperm[i]],'g',-1);if(i<KC) fprintferr(",");}
      fprintferr("]~\n");flusherr();
    }
    mit=hnfadd(mit,&pdep,&matarch,vperm,&matalpha,KCCOPRO,KC,col,&nlze,extramat,extramatc);
    if(DEBUGLEVEL){fprintferr("temps hnfadd = %ld\n",timer2());flusherr();}
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
    KCCOPRO+=extrarel;col=KCCOPRO-(lg(matalpha)-1);free(vinvperm);
    if(nlze)
    {
      nreldep++;
      if(nreldep>MAXRELSUP) 
      {
	desallocate(primsubfactorbase,cmptglob,vp,numprim,prim,numideal,idealbase,ex,matcopy,vperm,subvperm,vinvperm,flvp,phase);
	if(cbach>11.99) 
	 err(talker,"sorry, buchxxx is not able to compute this field PLEASE REPORT!!!");
	else
	{avma=av0;cbach=min(2*cbach,12);cbach2=max(cbach2,cbach);goto INCREASEGEN;}
      }
      else
      {
	phase=1;firstime=1;
	if(avma<limpile)
	{
	  tetpil=avma;mit=gcopy(mit);
	  pdep=gcopy(pdep);matarch=gcopy(matarch);matalpha=gcopy(matalpha);
	  dec=lpile(av1,tetpil,0)>>TWOPOTBYTES_IN_LONG;
	  mit+=dec;matalpha+=dec;matarch+=dec;pdep+=dec;
	}
	goto LABELINT;
      }
    }
  }
  phase=1;firstime=1;

/* Calcul du regulateur provisoire pour le check */

  if(DEBUGLEVEL)
  {
    fprintferr("\n#### Calcul du regulateur\n");flusherr();
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  }
  reg=compute_regulator(mit,col,matarch,RU,R1,N,PRECREG,&xarch,&xreal,&image_mdet,&sreg);
  if(gexpo(reg)<=-4)
  {
    nrelsup++;
    if(nrelsup>MAXRELSUP) 
    {
      desallocate(primsubfactorbase,cmptglob,vp,numprim,prim,numideal,idealbase,ex,matcopy,vperm,subvperm,vinvperm,flvp,phase);
      if(cbach>11.99) 
	err(talker,"sorry, buchxxx is not able to compute this field PLEASE REPORT!!!");
      else{avma=av0;cbach=min(2*cbach,12);cbach2=max(cbach2,cbach);goto INCREASEGEN;}
    }
    else 
    {
      if(DEBUGLEVEL) fprintferr("regulateur nul : \n");
      if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
      nlze=1;phase=1;firstime=1;
      if(avma<limpile)
      {
	tetpil=avma;mit=gcopy(mit);
	pdep=gcopy(pdep);matarch=gcopy(matarch);matalpha=gcopy(matalpha);
	dec=lpile(av1,tetpil,0)>>TWOPOTBYTES_IN_LONG;
	mit+=dec;matalpha+=dec;matarch+=dec;pdep+=dec;
      }
      goto LABELINT;
    }
  }

/* Calcul du nombre de classes */  

  if(DEBUGLEVEL)
  {
    fprintferr("\n#### Calcul du nombre de classes\n");flusherr();
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  }
  clh=compute_class_number(mit,&met,&u1);
  u1=reducemodmatrix(u1,mit);
  if(DEBUGLEVEL)
  {
    fprintferr("temps smith/groupe de classes: ");fprintferr("%ld\n",timer2());
    flusherr();
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  }

/* Calcul du check */
  
  if(DEBUGLEVEL)
  {
    fprintferr("\n#### Calcul du check\n");flusherr();
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  }
  c_1=compute_check(mit,col,image_mdet,xreal,RU,&reg,z,zu,clh,PRECREG,&parch);
  if(gcmp0(c_1))
  {
    fprintferr("\nYou are probably being dishonest, doubling the Bach constant\n");
    flusherr();
    avma=av0;cbach=min(2*cbach,12);cbach2=max(cbach2,cbach);goto INCREASEGEN;
  }
  if(gcmpgs(gmul2n(c_1,1),3)<0) {c_1=stoi(20);nrelsup=10*MAXRELSUP;}
  if(gcmpgs(c_1,3)>0) 
  {
    nrelsup++;
    if(nrelsup>MAXRELSUP) 
    {
      desallocate(primsubfactorbase,cmptglob,vp,numprim,prim,numideal,idealbase,ex,matcopy,vperm,subvperm,vinvperm,flvp,phase);
      if(cbach>11.99) 
	fprintferr("\n  ***   Warning: check is greater than 1.5, suggest increasing extra relations\n");
      else{avma=av0;cbach=min(2*cbach,12);cbach2=max(cbach2,cbach);goto INCREASEGEN;}
    }
    else 
    {
      if(DEBUGLEVEL) fprintferr("\n ***** check = %f\n\n",gtodouble(c_1)/2);
      if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
      nlze=1;phase=1;firstime=1;
      if(avma<limpile)
      {
	tetpil=avma;mit=gcopy(mit);
	pdep=gcopy(pdep);matarch=gcopy(matarch);matalpha=gcopy(matalpha);
	dec=lpile(av1,tetpil,0)>>TWOPOTBYTES_IN_LONG;
	mit+=dec;matalpha+=dec;matarch+=dec;pdep+=dec;
      }
      goto LABELINT;
    }
  }
  s=sfact;
  if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}

/* Phase "be honest" */
 
  if(DEBUGLEVEL)
  {
    if(KCZ2>KCZ)
    {
      fprintferr("be honest for primes from %ld to %ld\n",prim[KCZ+1],prim[KCZ2]);
      if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
    }
  }
  if((KCZ2>KCZ)&&(!flvp)){flHONEST=1;goto LABELVP;}
LABELHONEST:
  resuhonest=be_honest(KCZ,KCZ2,RU,PRECREGINT,IDN,s,nf,lgsub,idealbase,numideal,prim,numprim,q,vp,ex,exu);
  if(!resuhonest)
  {
    desallocate(primsubfactorbase,cmptglob,vp,numprim,prim,numideal,idealbase,ex,matcopy,vperm,subvperm,vinvperm,flvp,phase);
    if(cbach>11.99) 
      err(talker,"sorry, buchxxx is not able to compute this field PLEASE REPORT!!!");
    else{avma=av0;cbach=min(2*cbach,12);cbach2=max(cbach2,cbach);goto INCREASEGEN;}
  }
  if(DEBUGLEVEL)
  {
    fprintferr("temps be honest: ");fprintferr("%ld\n",timer2());
    flusherr();
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  }
  
/* Calcul des generateurs du groupe de classes */

  if(DEBUGLEVEL)
  {
    fprintferr("\n#### Calcul des generateurs du groupe des classes\n");
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  }
  RES[5]=(long)class_group_generators(N,RU,PRECREG,PRECREGINT,nf,met,clh,u1,vectbase,vperm);
  if(DEBUGLEVEL)
  {
    fprintferr("temps generateurs du groupe de classes: ");
    fprintferr("%ld\n",timer2());flusherr();
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  }
  
/* Affectation du regulateur, des racines de l'unite et des unites fondamentales */

  xarch=cleancol(gmul(xarch,parch),N,RU,PRECREG);
  RES[6]=(long)reg;z=gdiv(gmul(reg,clh),gmul((GEN)zu[1],z));RES[7]=(long)z;
  if(flun)
  {
    pgen1=cgetg(3,17);RES[8]=(long)pgen1;
    pgen1[1]=zu[1];pgen1[2]=lmul((GEN)nf[7],(GEN)zu[2]);
  }
  if(abs(flun)>1)
  {
    if(DEBUGLEVEL)
    {
      fprintferr("\n#### Calcul des unites fondamentales\n");
      if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
    }
    fu=getfu(nf,&xarch,reg,flun,&c,PRECREG);
    if(DEBUGLEVEL)
    {
      fprintferr("temps getfu: ");fprintferr("%ld\n",timer2());
      flusherr();
      if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
    }
    if((abs(flun)>2)&&(!c))
    {
      desallocate(primsubfactorbase,cmptglob,vp,numprim,prim,numideal,idealbase,ex,matcopy,vperm,subvperm,vinvperm,flvp,phase);
      prec=(prec<<1)-2;avma=av;
      if(DEBUGLEVEL)
      {
	fprintferr("Precision insuffisante pour fu, on double la precision.\n");
	fprintferr("Nouvelle precision prec = %ld\n",prec);
      }
      goto FORCEFU;
    }
    RES[9]=c?lmul((GEN)nf[7],fu):(long)fu;RES[10]=lstoi(c);
  }
 
/* Nettoyage et fin */
 
  if(DEBUGLEVEL)
  {
    fprintferr("\n#### Phase finale\n");
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  }
  tetpil=avma;
  if(flun>=0) {z=cgetg(2,19);z[1]=lcopy(RES);}
  else
  {
    z=cgetg(9,17);z[1]=lcopy(mit);z[2]=lcopy(matalpha);z[3]=lcopy(xarch);
    colnew=lg(matarch)-1-sreg;pgen1=cgetg(colnew+1,19);z[4]=(long)pgen1;
    for(j=1;j<=colnew;j++) pgen1[j]=(long)cleancol((GEN)matarch[j+sreg],N,RU,PRECREG);
    z[5]=lcopy(vectbase);pgen1=cgetg(KC+1,18);z[6]=(long)pgen1;
    for(i=1;i<=KC;i++) pgen1[i]=lstoi(vperm[i]);z[7]=lcopy(nf);
    ncz=(flun==-1)?4:6;pgen1=cgetg(ncz+1,17);z[8]=(long)pgen1;
    for(i=1;i<=ncz;i++) pgen1[i]=lcopy((GEN)RES[i+4]);
  }
  desallocate(primsubfactorbase,cmptglob,vp,numprim,prim,numideal,idealbase,ex,matcopy,vperm,subvperm,vinvperm,flvp,phase);
  if(DEBUGLEVEL)
  {
    fprintferr("temps free finaux: ");fprintferr("%ld\n",timer2());
    flusherr();
  }
  if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  return gerepile(av,tetpil,z);
}

/* Calcul de la base des sous-facteurs pour generer les relations entre ideaux. */

GEN
subfactorbasegen(long N, long m, GEN vectbase, long minsfb, long *vperm, long* ptss)
{
  long av=avma,tetpil,i,j,k,s,s1,s2,e,nbidp,ss=0;
  GEN y,y1,y2,y3,p1,perm,perm1;
  double prod;
 
  k=lg(vectbase)-1;
  y=cgetg(k+1,18);y1=cgetg(k+1,18);y3=cgetg(k+1,18);
  s=0;s1=0;
  for(i=1;i<=k;i++) 
  {
    p1=(GEN)vectbase[i];y3[i]=(long)gpui((GEN)p1[1],(GEN)p1[4],0);
    if((i>1)&&cmpii((GEN)p1[1],(GEN)((GEN)vectbase[i-1])[1]))
    {
      if(s==N) y1[i-1]=zero;
      s=0;
      if(s1==N) ss++;
      s1=0;
    }
    s2=(e=itos((GEN)p1[3]))*itos((GEN)p1[4]);
    s+=s2;s1+=s2;
    if(e>1) {y1[i]=zero;s=0;}
    else y1[i]=y3[i];
  }
  if(s==N) y1[k]=zero;
  if(s1==N) ss++;
  perm=indexsort(y1);
  i=1;while((i<=k)&&(!signe((GEN)y1[itos((GEN)(perm[i]))]))) i++;
  if(i>k+1-minsfb) {*ptss= -1;avma=av;return gzero;}
  nbidp=0;prod=1.0;
  while(((nbidp<=(k-i))&&(prod<m+0.5))||(nbidp<minsfb))
  {nbidp++;prod*=gtodouble((GEN)y1[itos((GEN)perm[nbidp+i-1])]);}
  if(prod<m) {*ptss= -1;avma=av;return gzero;}
  for(j=1;j<=nbidp;j++) y3[itos((GEN)(perm[j+i-1]))]=zero;
  perm1=indexsort(y3);
  for(j=1;j<=nbidp;j++) vperm[j]=itos((GEN)(perm[j+i-1]));
  for(j=nbidp+1;j<=k;j++) vperm[j]=itos((GEN)(perm1[j]));
  tetpil=avma;y2=cgetg(nbidp+1,18);
  for(j=1;j<=nbidp;j++) y2[j]=lcopy((GEN)vectbase[itos((GEN)perm[j+i-1])]);
  *ptss=ss;return gerepile(av,tetpil,y2);
}

GEN
new_subfactorbasegen(long N, long m, GEN vectbase, long minsfb, long *vperm, long* ptss)
{
  long av=avma,tetpil,i,j,k,s,s1,e,nbidp,ss=0;
  GEN y,y1,y2,y3,p1,perm,perm1;
  double prod;
 
  k=lg(vectbase)-1;
  y=cgetg(k+1,18);y1=cgetg(k+1,18);y3=cgetg(k+1,18);
  s=0;s1=0;
  for(i=1;i<=k;i++) 
  {
    p1=(GEN)vectbase[i];y3[i]=(long)gpui((GEN)p1[1],(GEN)p1[4],0);
    if(!gcmp1((GEN)p1[4])) y1[i]=zero;
    else
    {
      e=itos((GEN)p1[3]);
      if((i>1)&&cmpii((GEN)p1[1],(GEN)((GEN)vectbase[i-1])[1]))	y1[i]=(e>1)?zero:p1[1];
      else
      {
	if(i>1) y1[i]=zero;
	else y1[i]=(e>1)?zero:p1[1];
      }
    }
  }
  perm=indexsort(y1);
  i=1;while((i<=k)&&(!signe((GEN)y1[itos((GEN)(perm[i]))]))) i++;
  if(i>k+1-minsfb) {*ptss= -1;avma=av;return gzero;}
  nbidp=0;prod=1.0;
  while(((nbidp<=(k-i))&&(prod<m+0.5))||(nbidp<minsfb))
  {nbidp++;prod*=gtodouble((GEN)y1[itos((GEN)perm[nbidp+i-1])]);}
  if(prod<m) {*ptss= -1;avma=av;return gzero;}
  for(j=1;j<=nbidp;j++) y3[itos((GEN)(perm[j+i-1]))]=zero;
  perm1=indexsort(y3);
  for(j=1;j<=nbidp;j++) vperm[j]=itos((GEN)(perm[j+i-1]));
  for(j=nbidp+1;j<=k;j++) vperm[j]=itos((GEN)(perm1[j]));
  tetpil=avma;y2=cgetg(nbidp+1,18);
  for(j=1;j<=nbidp;j++) y2[j]=lcopy((GEN)vectbase[itos((GEN)perm[j+i-1])]);
  *ptss=ss;return gerepile(av,tetpil,y2);
}

/*Calcul la table des puissances des ideaux premiers de la base des sous-facteurs ;
  w est la base des sous-facteurs, a est l'exposant maximum calcule : 
  renvoie une matrice x de GEN , x[j][i] contenant
  l'ideal (P_i)^j, sous forme HNF */

GEN **
powsubfactgen(GEN nf, GEN w, long a, long PRECREG,long PRECREGINT)
{
  long i,j,n=lg(w)-1,N=lgef((GEN)nf[1])-3,R1,R2,RU;
  GEN **x,pgen1,id;
 
  R1=itos((GEN)((GEN)nf[2])[1]);R2=(N-R1)/2;RU=R1+R2;
  x=(GEN**)malloc(sizeof(GEN*)*(n+1));
  for(i=1;i<=n;i++) x[i]=(GEN*)malloc(sizeof(GEN)*(a+1));
  id=cgetg(3,17);id[1]=(long)idmat(N);pgen1=cgetg(RU+1,17);
  id[2]=(long)pgen1;for(i=1;i<=RU;i++) pgen1[i]=zero;
  for(i=1;i<=n;i++) 
  {
    x[i][0]=id;
    if(DEBUGLEVEL>=2) {fprintferr(" 0");flusherr();}
    for(j=1;j<=a;j++) 
    {
      x[i][j]=(j==1)?idealmulprime(nf,id,(GEN)w[i]):idealmulprimeredall(nf,x[i][j-1],(GEN)w[i],PRECREG,PRECREGINT);
      if(DEBUGLEVEL>=2) {fprintferr(" %ld",j);flusherr();}
    }
    if(DEBUGLEVEL>=2) {fprintferr("\n");flusherr();}
  }
  return x;
}

/*Calcul de la base de facteurs : n2 est la borne pour les nombres premiers
  qui vont etre testes pour obtenir toutes les relations et la norme des ideaux
  choisis, n est la borne des nombres premiers qui vont etre testes pour les 
  relations primaires (constante de Bach).
  Cette fonction cree et affecte un certain nombre de variables :
  numprimfactorbase[i] est l'indice k tel que primfactorbase[k]=i (si i n'est
  pas premier, numprimfactorbase[i]=0), ; primfactorbase[i] contient le i-eme
  nombre premier utilise pour construire la base de facteurs ; numideal[i] est
  l'indice k tel que idealbase[k]=i ; idealbase[i] contient les ideaux premiers
  de norme convenable et au dessus du nombre premier numero i ;
  KCZ contient le nombre de nombres premiers utilises pour construire la base 
  de facteurs jusqu'a la constante de Bach, KC contient le nombre d'ideaux 
  premiers jusqu'a la constante de Bach ; KCZ2 contient le nombre de nombres 
  premiers de la base de facteurs au total ; enfin, la fonction renvoie le 
  nombre d'ideaux premiers utilises au total. On n'utilise que des nombres
  premiers ne divisant pas l'index F, et non inertes.*/

long
factorbasegen(GEN nf, long n2, long n, long **ptnumprim, long **ptprim, long **ptnum, GEN **ptideal, long *ptkc, long *ptkcz, long *ptkcz2, GEN *ptlfun)
{
  byteptr delta=diffptr;
  long av1,tetpil,i,j,k,pp,fl,lon,ip,*numprimfactorbase,*primfactorbase;
  long *numideal,ip0;
  GEN prim,pgen1,pgen2,pgen3,pgen4,*idealbase,lfun;
 
  numprimfactorbase=(long*)malloc(sizeof(long)*(n2+1));
  primfactorbase=(long*)malloc(sizeof(long)*(n2+1));
  numideal=(long*)malloc(sizeof(long)*(n2+1));
  idealbase=(GEN*)malloc(sizeof(GEN)*(n2+1));
  lfun=cgetr(DEFAULTPREC);affsr(1,lfun);*ptlfun=lfun;
  i=0;pp=*delta++;fl=1;ip=0;*ptkc=0;
  while(pp<=n2)
  {
    av1=avma;
    prim=stoi(pp);
    if(DEBUGLEVEL>=2) {fprintferr(" %ld",pp);flusherr();}
    pgen1=primedec(nf,prim);
    lon=lg(pgen1);divriz(mulir(subis(prim,1),lfun),prim,lfun);
    if((lon>2)||cmpis((GEN)((GEN)pgen1[1])[3],1)>0)
    {
      ip0=ip;
      pgen2=cgetg(lon,18);
      for(j=0,k=1;k<lon;k++)
      {
	pgen3=(GEN)pgen1[k];pgen4=gpui(prim,(GEN)pgen3[4],0);
	if(cmpis(pgen4,n2)<=0)
	{
	  ip++;j++;pgen2[j]=(long)pgen3;
	  divriz(mulir(pgen4,lfun),subis(pgen4,1),lfun);
	}
      }
      i++;numprimfactorbase[pp]=i;primfactorbase[i]=pp;
      numideal[pp]=ip0;tetpil=avma;pgen3=cgetg(j+1,17);
      for(k=1;k<=j;k++) pgen3[k]=lcopy((GEN)pgen2[k]);
      idealbase[i]=gerepile(av1,tetpil,pgen3);
    }
    else
    {
      pgen4=gpui(prim,(GEN)((GEN)pgen1[1])[4],0);
      if(cmpis(pgen4,n2)<=0) 
	divriz(mulir(pgen4,lfun),subis(pgen4,1),lfun);
    }
    pp+=*delta++;
    if((pp>n)&&fl) {*ptkc=ip;fl=0;*ptkcz=i;}
  }
  *ptkcz2=i;
  *ptnumprim=numprimfactorbase;*ptprim=primfactorbase;*ptnum=numideal;*ptideal=idealbase;
  if(DEBUGLEVEL>=2) fprintferr("\n");
  return ip;
}

long
factorisegen(GEN nf, GEN idealvec, long kcz, long limp, long *primfact, long *expoprimfact, long *primfactorbase, GEN *idealbase, long *numideal, long *numprimfactorbase)
{
  long sr,i,j,n1,ip,v,p,k,fl=1,av1,lo,ifinal;
  GEN x,q,r,pg1,p1,listexpo,ideal=(GEN)idealvec[1];
 
  av1=avma;lo=0;
  for(x=gun,i=1;i<lg(ideal);i++)
  {p1=gcoeff(ideal,i,i);if(!gcmp1(p1)) x=mulii(x,p1);}
  x=gdiv((GEN)idealvec[3],x);
  if(gcmp1(x)) {avma=av1;primfact[0]=0;return 1;}
  listexpo=cgeti(kcz+1);
  for(i=1;(i<=kcz)&&fl;i++)
  {
    p=primfactorbase[i];q=dvmdis(x,p,&r);k=0;
    if((sr=(!signe(r))))
    {
      x=q;while(sr) {k++;q=dvmdis(x,p,&r);if((sr=(!signe(r)))) x=q;}
    }
    listexpo[i]=k;
    fl=(cmpis(q,p)>0);
  }
  if(fl) {avma=av1;return 0;}
  else
  {
    ifinal=i-1;
    if(cmpis(x,limp)>0) {avma=av1;return 0;}
    for(i=1;i<=ifinal;i++)
    {
      p=primfactorbase[i];
      if((k=listexpo[i]))
      {
	pg1=idealbase[numprimfactorbase[p]];
	n1=lg(pg1);ip=numideal[p];
	for(j=1;(j<n1)&&k;j++)
	{
	  v=idealval(nf,ideal,(GEN)pg1[j])-element_val2(nf,(GEN)idealvec[2],(GEN)idealvec[3],(GEN)pg1[j]);
	  if(v) 
	  {
	    primfact[++lo]=ip+j;expoprimfact[lo]=v;
	    k+=v*itos((GEN)((GEN)pg1[j])[4]);
	  }
	}
	if(k) {avma=av1;return 0;}
      }
    }
    if(gcmp1(x)) {avma=av1;primfact[0]=lo;return 1;}
    else
    {
      p=itos(x);
      pg1=idealbase[numprimfactorbase[p]];
      n1=lg(pg1);ip=numideal[p];k=1;
      for(j=1;(j<n1)&&k;j++)
      {
	v=idealval(nf,ideal,(GEN)pg1[j])-element_val2(nf,(GEN)idealvec[2],(GEN)idealvec[3],(GEN)pg1[j]);
	if(v) 
	{
	  primfact[++lo]=ip+j;expoprimfact[lo]=v;
	  k+=v*itos((GEN)((GEN)pg1[j])[4]);
	}
      }
      if(k) {avma=av1;return 0;}
      avma=av1;primfact[0]=lo;return 1;
    }
  }
}

long
factorisealpha(GEN nf, GEN alpha, long kcz, long limp, long *primfact, long *expoprimfact, long *primfactorbase, GEN *idealbase, long *numideal, long *numprimfactorbase)
{
  long sr,i,j,n1,ip,v,p,k,fl=1,av1,lo,ifinal;
  GEN d,x,q,r,pg1,listexpo;
 
  av1=avma;lo=0;
  d=x=absi(subres(gmul((GEN)nf[7],alpha),(GEN)nf[1]));
  if(gcmp1(x)) {avma=av1;primfact[0]=0;return 1;}
  listexpo=cgeti(kcz+1);
  for(i=1;(i<=kcz)&&fl;i++)
  {
    p=primfactorbase[i];q=dvmdis(x,p,&r);k=0;
    if((sr=(!signe(r))))
    {
      x=q;while(sr) {k++;q=dvmdis(x,p,&r);if((sr=(!signe(r)))) x=q;}
    }
    listexpo[i]=k;
    fl=(cmpis(q,p)>0);
  }
  if(fl) {avma=av1;return 0;}
  else
  {
    ifinal=i-1;
    if(cmpis(x,limp)>0) {avma=av1;return 0;}
    for(i=1;i<=ifinal;i++)
    {
      p=primfactorbase[i];
      if((k=listexpo[i]))
      {
	pg1=idealbase[numprimfactorbase[p]];
	n1=lg(pg1);ip=numideal[p];
	for(j=1;(j<n1)&&k;j++)
	{
	  v=element_val2(nf,alpha,d,(GEN)pg1[j]);
	  if(v) 
	  {
	    primfact[++lo]=ip+j;expoprimfact[lo]=v;
	    k-=v*itos((GEN)((GEN)pg1[j])[4]);
	  }
	}
	if(k) {avma=av1;return 0;}
      }
    }
    if(gcmp1(x)) {avma=av1;primfact[0]=lo;return 1;}
    else
    {
      p=itos(x);
      pg1=idealbase[numprimfactorbase[p]];
      n1=lg(pg1);ip=numideal[p];k=1;
      for(j=1;(j<n1)&&k;j++)
      {
	v=element_val2(nf,alpha,d,(GEN)pg1[j]);
	if(v) 
	{
	  primfact[++lo]=ip+j;expoprimfact[lo]=v;
	  k-=v*itos((GEN)((GEN)pg1[j])[4]);
	}
      }
      if(k) {avma=av1;return 0;}
      avma=av1;primfact[0]=lo;return 1;
    }
  }
}

GEN
idealmulprimeredall(GEN nf, GEN x, GEN vp, long PRECREG,long PRECREGINT)
{
  long av=avma,tetpil;
  GEN y;
 
  y=idealmulprime(nf,x,vp);tetpil=avma;
  return gerepile(av,tetpil,ideallllredall(nf,y,gzero,PRECREG,PRECREGINT));
}

GEN
idealmulprimered(GEN nf, GEN x, GEN vp, long PRECREG)
{
  long av=avma,tetpil;
  GEN y;
 
  y=idealmulprime(nf,x,vp);tetpil=avma;
  return gerepile(av,tetpil,ideallllredall(nf,y,gzero,PRECREG,PRECREG));
}

GEN
cleancol(GEN x,long N,long RU,long PRECREG)
{
  long i,j,av=avma,tetpil,tx=typ(x),R1=RU+RU-N;
  GEN s,s2,p1,p2,p3,p4,y;
 
  if(tx<17) err(talker,"not a vector/matrix in cleancol");
  if(tx==19) 
  {
    y=cgetg(lg(x),tx);
    for(j=1;j<lg(x);j++) y[j]=(long)cleancol((GEN)x[j],N,RU,PRECREG);
    return y;
  }
  p1=greal(x);s=gzero;for(i=1;i<=RU;i++) s=gadd(s,(GEN)p1[i]);
  s=gdivgs(s,-N);if(N>R1) s2=gmul2n(s,1);
  p2=gmul2n(mppi(PRECREG),2);p3=gimag(x);
  tetpil=avma;y=cgetg(RU+1,tx);
  for(i=1;i<=RU;i++)
  {
    p4=cgetg(3,6);y[i]=(long)p4;
    p4[1]=(i<=R1)?ladd((GEN)p1[i],s):ladd((GEN)p1[i],s2);
    p4[2]=(long)gmod((GEN)p3[i],p2);
  }
  return gerepile(av,tetpil,y);
}

GEN
getfu(GEN nf, GEN *ptxarch, GEN reg, long flun, long *pte, long PRECREG)
{
  long av=avma,tetpil,i,j,RU,N=lgef((GEN)nf[1])-3,e,R1,R2,dec;
  GEN p1,p2,p3,y,matep,s,xarch;
 
  R1=itos((GEN)((GEN)nf[2])[1]);R2=(N-R1)>>1;RU=R1+R2;xarch=*ptxarch;
  if(RU==1) return cgetg(1,19);
  if(gexpo(reg)<-8) 
  {
    if(abs(flun)==2) fprintferr("\n  ***   Warning: not enough relations for fundamental units, not given\n");
    return cgetg(1,19);
  }
  matep=cgetg(RU,19);
  for(j=1;j<RU;j++)
  {
    s=gzero;for(i=1;i<=RU;i++) s=gadd(s,greal(gcoeff(xarch,i,j)));
    s=gdivgs(s,N);
    p1=cgetg(N+1,18);matep[j]=(long)p1;
    for(i=1;i<=R1;i++)
      p1[i]=lsub(gcoeff(xarch,i,j),s);
    for(i=R1+1;i<=RU;i++)
    {
      p1[i]=lsub(gmul2n(gcoeff(xarch,i,j),-1),s);
      p1[i+R2]=lconj((GEN)p1[i]);
    }
  }
  matep=gexp(gmul(matep,p1=lll(greal(matep),PRECREG)),PRECREG);
  xarch=gmul(xarch,p1);
  p1=(GEN)((GEN)nf[5])[1];
  p2=cgetg(N+1,19);
  for(j=1;j<=N;j++)
  {
    p3=cgetg(N+1,18);p2[j]=(long)p3;
    for(i=1;i<=R1;i++) p3[i]=(long)coeff(p1,i,j);
    for(i=R1+1;i<=RU;i++)
    {
      p3[i]=(long)coeff(p1,i,j);p3[i+R2]=lconj((GEN)p3[i]);
    }
  }
  y=greal(grndtoi(invmulmatreel(p2,matep),&e));
  if(e>=0)
  {
    if(abs(flun)==2) fprintferr("\n  ***   Warning: insufficient precision for fundamental units, not given\n");
    avma=av;*pte=0;return cgetg(1,19);
  }
  *pte= -e;
  for(j=1;j<lg(y);j++)
  {
    p1=(GEN)y[j];p2=(GEN)(ginv(gmodulcp(gmul((GEN)nf[7],p1),(GEN)nf[1]))[2]);
    p3=cgetg(N+1,18);for(i=1;i<lgef(p2)-1;i++) p3[i]=p2[i+1];
    for(;i<=N;i++) p3[i]=zero;
    p2=gmul((GEN)nf[8],p3);
    if(gcmp(gnorml2(p2),gnorml2(p1))<0) 
    {
      p1=p2;xarch[j]=lneg((GEN)xarch[j]);
    }
    for(i=N;(i>=1)&&gcmp0((GEN)p1[i]);i--);
    if(gsigne((GEN)p1[i])<0)
    {
      y[j]=lneg(p1);
      p1=cgetg(RU+1,18);p2=gmul(gi,mppi(PRECREG));
      for(i=1;i<=R1;i++) p1[i]=(long)p2;
      for(;i<=RU;i++) p1[i]=lmul2n(p2,1);
      xarch[j]=(long)gadd((GEN)xarch[j],p1);
    }
    else y[j]=(long)p1;
  }
  p1=gmul((GEN)nf[7],y);
  for(j=1;j<lg(y);j++) 
    if(!gcmp1(gabs(gnorm(gmodulcp((GEN)p1[j],(GEN)nf[1])),0)))
    {
      if(abs(flun)==2) fprintferr("\n  ***   Warning: fundamental units too large, not given\n");
      avma=av;*pte=0;return cgetg(1,17);
    }
  tetpil=avma;y=gcopy(y);xarch=gcopy(xarch);
  dec=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;*ptxarch=xarch+dec;return y+dec;
}

GEN
buchfu(GEN bnf)
{
  GEN nf,xarch,reg,res,fu,y,sig;
  long av=avma,tetpil,c,RU;

  nf=(GEN)bnf[7];sig=(GEN)nf[2];RU=itos((GEN)sig[1])+itos((GEN)sig[2]);
  res=(GEN)bnf[8];
  if((lg(res)==7)&&(lg((GEN)res[5])==RU))
  {y=cgetg(3,17);y[1]=lcopy((GEN)res[5]);y[2]=lcopy((GEN)res[6]);return y;}
  else
  {
    xarch=(GEN)bnf[3];reg=(GEN)res[2];
    fu=getfu(nf,&xarch,reg,2,&c,gprecision(xarch));
    tetpil=avma;y=cgetg(3,17);
    y[1]=c?lmul((GEN)nf[7],fu):lcopy(fu);y[2]=lstoi(c);
    return gerepile(av,tetpil,y);
  }
}

GEN
gcdrealnoer(GEN a, GEN b, long *pte)
{
  long av,tetpil,e;
  GEN k1,r;
 
  if(typ(a)==1)
  {
    if(typ(b)==1) return mppgcd(a,b);
    k1=cgetr(lg(b));affir(a,k1);a=k1;
  }
  else if(typ(b)==1) {k1=cgetr(lg(a));affir(b,k1);b=k1;}
  if(expo(a)<-5) return absr(b);
  if(expo(b)<-5) return absr(a);
  av=avma;a=absr(a);b=absr(b);
  while((expo(b)>=(-5))&&(signe(b)))
  {
    k1=gcvtoi(divrr(a,b),&e);if(e>0) {avma=av;*pte=e;return gzero;}
    r=subrr(a,mulir(k1,b));a=b;b=r;
  }
  tetpil=avma;*pte=expo(b);return gerepile(av,tetpil,gabs(a,0));
}

long
factorisegensimple(GEN nf, GEN ideal, long *primfact, long *expoprimfact, GEN vectbase)
{
  long i,p,v,fl=1,av1,av2,lo;
  GEN x,p1;
 
  av1=avma;lo=0;
  for(x=gun,i=1;i<lg(ideal);i++) x=mulii(x,gcoeff(ideal,i,i));
  if(gcmp1(x)) {avma=av1;primfact[0]=0;return 1;}
  av2=avma;fl=1;
  for(i=1;(i<lg(vectbase))&&fl;i++)
  {
    p1=(GEN)vectbase[i];p=itos((GEN)p1[1]);
    if(!signe(modis(x,p)))
    {
      v=idealval(nf,ideal,p1);
      if(v) 
      {
	primfact[++lo]=i;expoprimfact[lo]=v;
	x=divii(x,gpui((GEN)p1[1],mulsi(v,(GEN)p1[4]),0));
	if(gcmp1(x)) fl=0;
      }
    }
  }
  primfact[0]=lo;return 1-fl;
}

#define MAXITER 100

GEN
isprincipalall(GEN bignf, GEN x, long flall)
{
  long expoprimfact[500],primfact[500],av=avma,tetpil,*vinvperm;
  long i,j,fpc,colmit,colnew,k,N,R1,R2,RU,KC,e,pr,c,ss,nbtest,bou,q,tx=typ(x);
  GEN xinit,id,xalpha,yalpha,mit,matalpha,matunit,matalphac,vectbase,vperm,nf,a,RES;
  GEN u1u2,u1,u2,y,p1,p2,p3,p4,p5,s,s1,s2,om1,om2,om3,xar,pgen1,pgen2,ideal,vdir,exu;
  GEN detp5,p5min,orderp5,detpgen1,cycliclist;
  long *ex,expomax;
 
  if((typ(bignf)!=17)||(lg(bignf)!=9))
    err(talker,"not a big number field vector in isprincipal");
  mit=(GEN)bignf[1];matalpha=(GEN)bignf[2];matunit=(GEN)bignf[3];
  matalphac=gcopy((GEN)bignf[4]);vectbase=(GEN)bignf[5];vperm=(GEN)bignf[6];
  nf=(GEN)bignf[7];if((typ(nf)!=17)||(lg(nf)<10)) err(idealer1);
  N=lgef((GEN)nf[1])-3;id=idmat(N);
  RES=(GEN)bignf[8];KC=lg(vperm)-1;
  if((tx==9)||(tx==10))
  {
    if(tx==9)
    {
      if(!gegal((GEN)x[1],(GEN)nf[1]))
	err(talker,"not the same number field in isprincipal");
      x=(GEN)x[2];
    }
    avma=av;y=cgetg(4,17);c=lg((GEN)((GEN)RES[1])[2])-1;
    p1=cgetg(c+1,18);y[1]=(long)p1;for(i=1;i<=c;i++) p1[i]=zero;
    y[2]=lmul((GEN)nf[8],x);y[3]=lstoi(BIGINT);
    return y;
  }
  if((tx==17)&&(lg(x)==3)) x=(GEN)x[1];
  if((typ(x)==17)&&(lg(x)==6)) x=idealmulprime(nf,id,x);
  if((typ(x)!=19)||(lg(x)!=(N+1))||(lg((GEN)x[1])!=(N+1)))
    err(talker,"not an ideal in isprincipal");
  if(N==1)
  {y=cgetg(4,17);y[1]=lgetg(1,18);y[2]=x[1];y[3]=lstoi(BIGINT);return y;}
  R1=itos((GEN)((GEN)nf[2])[1]);R2=(N-R1)>>1;RU=R1+R2;
  pr=gprecision((GEN)matalphac[1]);if(!pr) pr=8;
  colmit=lg(mit)-1;colnew=lg(matalpha)-1;
  a=content(x);x=gdiv(x,a);xinit=x=hnfmod(x,detint(x));
  vinvperm=(long*)malloc(lg(vectbase)<<TWOPOTBYTES_IN_LONG);
  for(i=1;i<lg(vectbase);i++) vinvperm[itos((GEN)vperm[i])]=i;
  s=gun;for(i=1;i<=N;i++) s=mulii(s,gcoeff(x,i,i));
  xar=cgetg(RU+1,17);
  for(i=1;i<=RU;i++)
  {
    p2=cgetg(3,6);xar[i]=(long)p2;
    p3=cgetr(pr);p2[1]=(long)p3;
    p3=cgetr(pr);p2[2]=(long)p3;
    affsr(0,(GEN)p2[1]);affsr(0,(GEN)p2[2]);
  }
  fpc=factorisegensimple(nf,x,primfact,expoprimfact,vectbase);
  nbtest=0;ex=(long*)malloc(3*sizeof(long));exu=cgetg(RU+1,17);
  q=BITS_IN_RANDOM-1-(long)ceil(log((double)CBUCHG)/log(2.0));
  if(!fpc) 
  {
    p3=cgetg(3,17);p3[1]=(long)x;p3[2]=(long)xar;
    p1=ideallllred(nf,p3,gzero,pr);
    x=(GEN)p1[1];xar=cleancol((GEN)p1[2],N,RU,pr);
    fpc=factorisegensimple(nf,x,primfact,expoprimfact,vectbase);
    if(!fpc) 
    {
      vdir=cgetg(RU+1,17);for(i=1;i<=RU;i++) vdir[i]=zero;
      for(i=1;(i<=RU)&&(!fpc);i++)
      {
	vdir[i]=lstoi(10);if(i>1) vdir[i-1]=zero;
	p1=ideallllred(nf,p3,vdir,pr);
	x=(GEN)p1[1];xar=cleancol((GEN)p1[2],N,RU,pr);
	fpc=factorisegensimple(nf,x,primfact,expoprimfact,vectbase);
      }
      while((!fpc)&&(nbtest<MAXITER))
      {
	for(i=1;i<=2;i++) ex[i]=mymyrand()>>q;
	ideal=idealpowprime(nf,(GEN)vectbase[itos((GEN)vperm[1])],stoi(ex[1]),pr);
	pgen1=idealpowprime(nf,(GEN)vectbase[itos((GEN)vperm[2])],stoi(ex[2]),pr);
	ideal=idealmulh(nf,idealmul(nf,p3,ideal),pgen1);
	for(bou=1;(bou<=RU)&&(!fpc);bou++)
	{
	  if(bou==1) 
	  {for(i=1;i<=RU;i++) exu[i]=lstoi(mymyrand()>>q);}
	  else
	  {
	    for(i=1;i<=RU;i++) exu[i]=zero;
	    exu[bou]=lstoi(10);
	  }
	  p1=ideallllred(nf,ideal,exu,pr);nbtest++;
	  if(DEBUGLEVEL>2)
	  {
	    fprintferr("nbtest = %ld, ideal = ",nbtest);
	    outerr((GEN)p1[1]);
	  }
	  x=(GEN)p1[1];xar=cleancol((GEN)p1[2],N,RU,pr);
	  fpc=factorisegensimple(nf,x,primfact,expoprimfact,vectbase);
	}
      }
      if(nbtest>=MAXITER) 
      {
	free(ex);free(vinvperm);
	err(talker,"too many iterations in isprincipal");
      }
    }
  }
  xalpha=cgetg(colmit+1,18);yalpha=cgetg(colnew+1,18);
  for(i=1;i<=colmit;i++) xalpha[i]=zero;
  for(i=1;i<=colnew;i++) yalpha[i]=zero;
  if(nbtest)
  {
    for(i=1;i<=2;i++) 
    {
      if(i<=colmit) xalpha[i]=lstoi(-ex[i]);
      else yalpha[i-colmit]=lstoi(-ex[i]);
    }
  }
  free(ex);
  for(i=1;i<=primfact[0];i++)
  {
    k=vinvperm[primfact[i]];
    if(k<=colmit) xalpha[k]=laddsg(expoprimfact[i],(GEN)xalpha[k]);
    else yalpha[k-colmit]=laddsg(expoprimfact[i],(GEN)yalpha[k-colmit]);
  }
  u1u2=smith2(mit);u1=(GEN)u1u2[1];u2=(GEN)u1u2[2];
  cycliclist=(GEN)((GEN)RES[1])[2];
  p1=gmul(u1,gsub(xalpha,gmul(matalpha,yalpha)));
  c=lg(cycliclist)-1;u1=reducemodmatrix(ginv(u1),mit);
  p4=cgetg(colmit+colnew+1,18);p2=cgetg(c+1,18);
  for(i=1;i<=c;i++) 
  {
    p4[i]=(long)dvmdii((GEN)p1[i],(GEN)cycliclist[i],(GEN*)(p2+i));
    if(signe((GEN)p2[i])<0) 
    {
      p2[i]=ladd((GEN)p2[i],(GEN)cycliclist[i]);
      p4[i]=(long)gaddgs((GEN)p4[i],-1);
    }
  }
  if(flall)
  {
    for(;i<=colmit;i++) p4[i]=p1[i];
    for(;i<=colmit+colnew;i++) p4[i]=yalpha[i-colmit];
    p5=cgetg(colmit+1,19);for(i=1;i<=colmit;i++) p5[i]=matalphac[i];
    p3=gmul(p5,u2);for(i=1;i<=colmit;i++) matalphac[i]=p3[i];
    settyp(xar,18);p1=gsub(gmul(matalphac,p4),xar);
    p4=cgetg(c+1,19);
    for(j=1;j<=c;j++)
    {
      p5=(GEN)idealpowprime(nf,(GEN)vectbase[itos((GEN)(vperm[1]))],pgen1=gcoeff(u1,1,j),pr);
      if(signe(pgen1)<0) p5[1]=(long)numer((GEN)p5[1]);
      for(i=2;i<=colmit;i++)
      {
	pgen1=gcoeff(u1,i,j);ss=signe(pgen1);
	if(ss)
	{
	  if(ss>0)
	    p5=idealmulh(nf,p5,(GEN)idealpowprime(nf,(GEN)vectbase[itos((GEN)(vperm[i]))],pgen1,pr));
	  else
	    p5=idealmulh(nf,p5,numer((GEN)idealpowprime(nf,(GEN)vectbase[itos((GEN)(vperm[i]))],pgen1,pr)));
	  p5=ideallllred(nf,p5,gzero,pr);
	}
      }
      detp5=det((GEN)p5[1]);p5min=p5;pgen1=p5;
      pgen2=idealinv(nf,pgen1);pgen2[1]=(long)numer((GEN)pgen2[1]);
      detpgen1=det((GEN)pgen2[1]);
      if(gcmp(detpgen1,detp5)<0){detp5=detpgen1;p5min=pgen2;}
      orderp5=(GEN)cycliclist[j];
      expomax=itos(gmin(orderp5,stoi(50)));
      for(i=2;i<expomax;i++)
      {
	pgen1=idealmulred(nf,pgen1,p5,pr);
	if(gcmp1(mppgcd(stoi(i),orderp5)))
	{
	  detpgen1=det((GEN)pgen1[1]);
	  if(gcmp(detpgen1,detp5)<0){detp5=detpgen1;p5min=pgen1;}
	  pgen2=idealinv(nf,pgen1);pgen2[1]=(long)numer((GEN)pgen2[1]);
	  detpgen1=det((GEN)pgen2[1]);
	  if(gcmp(detpgen1,detp5)<0){detp5=detpgen1;p5min=pgen2;}
	}
      }
      p5=p5min;
      if(!gegal((GEN)p5[1],(GEN)((GEN)((GEN)RES[1])[3])[j]))
      {free(vinvperm);err(talker,"bug1 in isprincipal");}
      p4[j]=lneg((GEN)p5[2]);settyp((GEN)p4[j],18);
    }
    p1=cleancol(c?gadd(p1,gmul(p4,p2)):p1,N,RU,pr);
    if(RU>1)
    {
      s2=gzero;
      p4=cgetg(RU+1,19);
      for(j=1;j<RU;j++)
      {
	p5=cgetg(RU+1,18);p4[j]=(long)p5;
	s1=gzero;
	for(i=1;i<RU;i++) 
	{
	  p5[i]=(long)greal(gcoeff(matunit,i,j));s1=gadd(s1,gmul((GEN)p5[i],(GEN)p5[i]));
	}
	p5[RU]=zero;if(gcmp(s1,s2)>0) s2=s1;
      }
      p5=cgetg(RU+1,18);p4[RU]=(long)p5;
      for(i=1;i<RU;i++) p5[i]=(long)greal((GEN)p1[i]);
      s2=gsqrt(gmul2n(s2,RU+1),pr);if(gcmpgs(s2,100000000)<0) s2=stoi(100000000);
      p5[RU]=(long)s2;
      p4=(GEN)lll(p4,pr)[RU];
      if(signe((GEN)p4[RU])<0) p4=gneg(p4);
      if(!gcmp1((GEN)p4[RU]))
      {free(vinvperm);err(talker,"bug2 in isprincipal");}
      setlg(p4,RU);
      p1=gadd(p1,gmul(matunit,p4));
      setlg(p4,RU+1);
    }
    s2=gun;
    for(j=1;j<=c;j++) 
    {
      p5=(GEN)((GEN)((GEN)RES[1])[3])[j];
      s1=gun;for(i=1;i<=N;i++) s1=mulii(s1,gcoeff(p5,i,i));
      if(signe((GEN)p2[j])) s2=mulii(s2,gpui(s1,(GEN)p2[j],0));
    }
    s=gdivgs(glog(gdiv(s,s2),pr),N);
    p4=cgetg(N+1,18);
    for(i=1;i<=R1;i++) p4[i]=(long)gexp(gadd(s,(GEN)p1[i]),pr);
    for(i=R1+1;i<=RU;i++)
    {
      p4[i]=(long)gexp(gadd(s,gmul2n((GEN)p1[i],-1)),pr);;
      p4[i+R2]=lconj((GEN)p4[i]);
    }
    om1=(GEN)((GEN)nf[5])[1];
    om2=cgetg(N+1,19);
    for(j=1;j<=N;j++)
    {
      om3=cgetg(N+1,18);om2[j]=(long)om3;
      for(i=1;i<=R1;i++) om3[i]=(long)coeff(om1,i,j);
      for(i=R1+1;i<=RU;i++)
      {
	om3[i]=(long)coeff(om1,i,j);om3[i+R2]=lconj((GEN)om3[i]);
      }
    }
    p1=gdiv(grndtoi(gmul(s2,greal(gauss(om2,p4))),&e),s2);
    if(e<0)
    {
      p3=cgetg(2,19);p3[1]=(long)p1;
      p5=(GEN)((GEN)RES[1])[3];
      if(!c) p3=idealmul(nf,p3,id);
      else for(j=1;j<=c;j++) p3=idealmul(nf,p3,idealpow(nf,(GEN)p5[j],(GEN)p2[j]));
      if(!gegal(xinit,p3)) e=0;
    }
    if(e>=0)
    {
      fprintferr("\n ***   Warning: insufficient precision for generators, not given\n");
      tetpil=avma;y=cgetg(4,17);y[1]=lcopy(p2);y[2]=lgetg(1,18);y[3]=zero;
    }
    else{tetpil=avma;y=cgetg(4,17);y[1]=lcopy(p2);y[2]=lmul(a,p1);y[3]=lstoi(-e);}
  }
  else {tetpil=avma;y=gcopy(p2);}
  free(vinvperm);return gerepile(av,tetpil,y);
}

GEN
isprincipal(GEN bignf, GEN x)
{
  return isprincipalall(bignf,x,0);
}

GEN
isprincipalgen(GEN bignf, GEN x)
{
  return isprincipalall(bignf,x,1);
}

#undef MAXRELSUP
#undef MAXITER

GEN
isunit(GEN bignf, GEN x)
{
  long av=avma,tetpil,tx,i,R1,R2,RU,nru;
  GEN RES,matunit,y,p1,p2,p3,nf,ro1;
 
  if((typ(bignf)!=17)||(lg(bignf)!=9))
    err(talker,"not a big number field vector in isunit");
  tx=typ(x);
  if((tx==4)||(tx==5)) return cgetg(1,17);
  matunit=(GEN)bignf[3];RU=lg(matunit);
  RES=(GEN)bignf[8];nf=(GEN)bignf[7];
  ro1=(GEN)RES[4];nru=itos((GEN)ro1[1]);
  if(tx==1)
  {
    if(!gcmp1(absi(x))) return cgetg(1,17);
    y=cgetg(RU+1,17);
    for(i=1;i<RU;i++) y[i]=zero;
    y[RU]=(signe(x)>0)?lmodulcp(gzero,(GEN)ro1[1]):lmodulcp(stoi(nru>>1),(GEN)ro1[1]);
    return y;
  }
  if(tx!=10) 
  {
    if(tx!=9) err(talker,"not an algebraic number in isunit");
    if(!gegal((GEN)nf[1],(GEN)x[1])) err(talker,"not the same number field in isunit");
  }
  p1=(GEN)nf[2];R1=itos((GEN)p1[1]);R2=itos((GEN)p1[2]);
  p1=cgetg(RU+1,18);for(i=1;i<=R1;i++) p1[i]=un;for(;i<=RU;i++) p1[i]=deux;
  p3=(GEN)principalidele(nf,x)[2];
  p1=concat(matunit,p1);p2=ground(gauss(greal(p1),greal(p3)));
  if(!gcmp0((GEN)p2[RU])) err(talker,"insufficient precision (1) in isunit");
  p1=gsub(p3,gmul(p1,p2));
  p3=(GEN)principalidele(nf,(GEN)ro1[2])[2];
  p1=gdiv((GEN)p1[1],(GEN)p3[1]);
  p1=gmod(ground(greal(p1)),(GEN)ro1[1]);
  tetpil=avma;y=cgetg(RU+1,17);for(i=1;i<RU;i++) y[i]=lcopy((GEN)p2[i]);
  y[RU]=lmodulcp(p1,(GEN)ro1[1]);
  return gerepile(av,tetpil,y);
}

GEN
signunit(GEN bignf)
{
  long av,i,j,r,R1,RU;
  GEN matunit,y,p1,p2,nf;
 
  if((typ(bignf)!=17)||(lg(bignf)!=9))
    err(talker,"not a big number field vector in isunit");
  matunit=(GEN)bignf[3];RU=lg(matunit);nf=(GEN)bignf[7];
  p1=(GEN)nf[2];R1=itos((GEN)p1[1]);constpi(5);
  y=cgetg(RU,19);
  for(j=1;j<RU;j++)
  {
    p1=cgetg(R1+1,18);y[j]=(long)p1;
    for(i=1;i<=R1;i++)
    {
      av=avma;p2=ground(gdiv(gimag(gcoeff(matunit,i,j)),gpi));
      r=signe(p2)?p2[lgef(p2)-1]&1:0;avma=av;
      p1[i]=r?lneg(gun):un;
    }
  }
  return y;
}

GEN
buchall_for_degree_one_pol(long v, long flun)
{
  GEN RES,z,pgen1;
  long av,tetpil,i,ncz;

  av=avma;
  if(abs(flun)>1) RES=cgetg(11,17);else RES=flun?cgetg(9,17):cgetg(8,17);
  RES[1]=(long)polx[v];pgen1=cgetg(3,17);RES[2]=(long)pgen1;pgen1[1]=un;
  pgen1[2]=zero;pgen1=cgetg(3,17);RES[3]=(long)pgen1;pgen1[1]=un;pgen1[2]=un;
  pgen1=cgetg(2,17);RES[4]=(long)pgen1;pgen1[1]=un;
  pgen1=cgetg(4,17);RES[5]=(long)pgen1;pgen1[1]=un;pgen1[2]=lgetg(1,17);
  pgen1[3]=lgetg(1,17);RES[6]=un;RES[7]=un;
  if(flun){pgen1=cgetg(3,17);RES[8]=(long)pgen1;pgen1[1]=deux;pgen1[2]=lneg(gun);}
  if(abs(flun)>1) {RES[9]=lgetg(1,19);RES[10]=lstoi(EXP220);}
  if(flun>=0)
  {z=cgetg(2,19);z[1]=(long)RES;tetpil=avma;return gerepile(av,tetpil,gcopy(z));}
  else
  {
    z=cgetg(9,17);z[1]=lgetg(1,19);z[2]=lgetg(1,19);z[3]=lgetg(1,19);
    z[4]=lgetg(1,19);z[5]=lgetg(1,18);z[6]=lgetg(1,17);
    z[7]=(long)initalg(polx[v],MEDDEFAULTPREC);ncz=(flun==-1)?4:6;
    pgen1=cgetg(ncz+1,17);z[8]=(long)pgen1;
    for(i=1;i<=ncz;i++) pgen1[i]=lcopy((GEN)RES[i+4]);return z;
  }
}

GEN
quad_form(GEN *cbase, GEN ma, GEN ideal, long PRECREG)
{
  long av,tetpil,dec;
  GEN pg3,pg5,rrr;
 
  av=avma;pg3=gmul(gtrans(ideal),gmul(ma,ideal));*cbase=lllgram(pg3,PRECREG);
  if(DEBUGLEVEL>6)
  {
    fprintferr(" matrice d'entree dans lllgram :\n");
    outerr(pg3);
    fprintferr(" matrice donnee par lllgram : PRECREG = %ld\n",PRECREG);
    outerr(*cbase);
  }
  pg5=gmul(gtrans(*cbase),gmul(pg3,*cbase));
  rrr=sqred1(pg5);
  tetpil=avma;rrr=gcopy(rrr);*cbase=gcopy(*cbase);
  dec=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;
  rrr+=dec;*cbase+=dec;return rrr;
}

long
small_norm_for_buchall(long s,long **mat,GEN matarch,long KC,long KCCO,long KCZ,long LIMC,long N,GEN D,long RU,long R1,long PRECREG,GEN P,GEN nf,GEN vectbase,long *prim,GEN *idealbase,long *numideal,long *numprim,GEN IDN,GEN gborne,long nbrelpid)
{
  long av1,av2,av3,n2,nbsmallnorm,nbsmallfact,nbrel,nbrelideal;
  long *x,i,j,k,l,lll,noideal,fpc,expoprimfact[500],primfact[500];
  long t,trouve,nonnul,ra,R2;
  double eps,MINKOVSKI_BOUND,*y,*zz,**qq;
  double IDEAL_BOUND,*vv,normideal,p,quotfact;
  GEN pgen1,pgen2,pgen3,pgen4,pgen5,pg3,pg4;
  GEN matpro,vecindep,cq,bq,ma,ideal,rrr,cbase,gnormideal;

  av1=avma;nbrel=KCCO;n2=N+1;eps=0.000001;t=s;R2=(N-R1)>>1;pgen2=mppi(PRECREG);
  pgen1=gdiv(gsqrt(gmulsg(N,gmul2n(pgen2,1)),PRECREG),gexp(stoi(N),PRECREG));
  pgen1=gmulsg(N,gpui(pgen1,dbltor(2./(double)N),PRECREG));
  pgen3=gpui(gdivsg(4,pgen2),gdiv(stoi(2*R2),stoi(N)),PRECREG);pgen1=gmul(pgen1,pgen3);
  MINKOVSKI_BOUND=gtodouble(gmul(pgen1,gpui(gabs(D,0),dbltor(1./(double)N),PRECREG)));
  if(DEBUGLEVEL)
  {fprintferr("Borne pour les normes = %.0f\n",MINKOVSKI_BOUND);flusherr();}
  nbsmallnorm=nbsmallfact=0;
  x=(long*)malloc(n2*sizeof(long));y=(double*)malloc(n2*sizeof(double));
  zz=(double*)malloc(n2*sizeof(double));vv=(double*)malloc(n2*sizeof(double));
  qq=(double**)malloc(n2*sizeof(double*));
  for(k=1;k<=N;k++) qq[k]=(double*)malloc(n2*sizeof(double));
  noideal=0;ma=(GEN)((GEN)nf[5])[3];
  while((t<nbrel)&&(noideal<KC))
  {
    noideal++;nbrelideal=0;av2=avma;ideal=(GEN)(vectbase[KC+1-noideal]);
    normideal=gtodouble(gnormideal=gpui((GEN)ideal[1],(GEN)ideal[4],0));
    if(DEBUGLEVEL>1)
    {
      fprintferr("\n*** Ideal no %ld : S = %ld , ",noideal,t);
      fprintferr("nombre premier = %ld , ",itos((GEN)ideal[1]));
      fprintferr("ideal = ");outerr(ideal);flusherr();
    }
    ideal=idealmul(nf,ideal,IDN);rrr=quad_form(&cbase,ma,ideal,PRECREG);
    for(k=1;k<=N;k++)
    {
      vv[k]=gtodouble(gcoeff(rrr,k,k));
      for(l=1;l<k;l++) qq[l][k]=gtodouble(gcoeff(rrr,l,k));
    }
    if(DEBUGLEVEL>3)
    {
      for(k=1;k<=N;k++) fprintferr("vv[%ld]=%.0f ",k,vv[k]);
      fprintferr("\n");flusherr();
    }
    IDEAL_BOUND=gtodouble(gborne)*pow(normideal,2./(double)N)*MINKOVSKI_BOUND;
    if(DEBUGLEVEL>1){fprintferr("IDEAL_BOUND  = %.0f\n",IDEAL_BOUND);flusherr();}
    k=N;y[N]=zz[N]=0;x[N]=(long)(sqrt(IDEAL_BOUND/vv[N]+eps));
    do
    {
      do
      {
	if(k>1)
	{
	  k--;zz[k]=0;for(l=k+1;l<=N;l++) zz[k]=zz[k]+qq[k][l]*x[l];
	  p=x[k+1]+zz[k+1];y[k]=y[k+1]+p*p*vv[k+1];
	  x[k]=(long)floor(sqrt((IDEAL_BOUND-y[k]+eps)/vv[k])-zz[k]);
	}
	while(vv[k]*(x[k]+zz[k])*(x[k]+zz[k])>IDEAL_BOUND-y[k]+eps){k++;x[k]--;}
      }
      while(k>1);
      if((nonnul=(x[1]||(y[1]>eps))))
      {
	if(ccontent(x,N)==1)
	{
	  av3=avma;pg3=cgetg(N+1,18);for(l=1;l<=N;l++) pg3[l]=lstoi(x[l]);pg4=pg3;
	  pg3=gmul(ideal,gmul(cbase,pg3));j=N;while((j>=2)&&(!signe((GEN)pg3[j]))) --j;
	  if(j>1)
	  {
	    if(DEBUGLEVEL>6)
	    {
	      outerr(gdiv(absi(subres(gmul((GEN)nf[7],pg3),(GEN)nf[1])),gnormideal));
	      bq=gzero;
	      for(l=1;l<=N;l++)
	      {
		cq=gzero;
		for(lll=l+1;lll<=N;lll++)
		  cq=gadd(cq,gmul(gcoeff(rrr,l,lll),(GEN)pg4[lll]));
		cq=gadd(cq,(GEN)pg4[l]);bq=gadd(bq,gmul(gmul(cq,cq),gcoeff(rrr,l,l)));
	      }
	      outerr(bq);
	    }
	    if((DEBUGLEVEL>1)&&(DEBUGLEVEL<7)){fprintferr(".");flusherr();}
	    fpc=factorisealpha(nf,pg3,KCZ,LIMC,primfact,expoprimfact,prim,idealbase,numideal,numprim);
	    nbsmallnorm++;
	    if(fpc)
	    {
	      t++;
	      for(l=1;l<=primfact[0];l++) mat[t][primfact[l]]=expoprimfact[l];
	      trouve=already_found_relation_of_long(mat,t,KC);
	      if(trouve)
	      {
		if(DEBUGLEVEL>1)
		{
		  fprintferr("rel. elim: %ld",trouve);
		  if(DEBUGLEVEL>2)
		  {
		    fprintferr(" : [");
		    for(l=1;l<=N;l++){fprintferr(" ");bruterr((GEN)(pg3[l]),'g',-1);}
		    fprintferr(" ]");
		  }
		  fprintferr("\n");flusherr();
		  if(DEBUGLEVEL>1)
		  {
		    fprintferr("[");
		    for(l=1;l<=KC;l++) fprintferr("%4ld",mat[t][l]);
		    fprintferr("]\n");flusherr();
		  }
		}
		for(l=1;l<=KC;l++) mat[t][l]=0;t--;
	      }
	      else
	      {
		nbrelideal++;nbsmallfact++;
		if(DEBUGLEVEL==1){fprintferr("%4ld",t);flusherr();}
		if(DEBUGLEVEL>1)
		{
		  fprintferr("t = %ld",t);
		  if(DEBUGLEVEL>2)
		  {
		    fprintferr("[");
		    for(l=1;l<=N;l++){fprintferr(" ");bruterr((GEN)(pg3[l]),'g',-1);}
		    fprintferr(" ]");
		  }
		  fprintferr("\n");flusherr();
		  if(DEBUGLEVEL>1)
		  {
		    fprintferr("[");for(l=1;l<=KC;l++)fprintferr("%4ld",mat[t][l]);
		    fprintferr("]\n");flusherr();
		  }
		}
		pgen5=gmul((GEN)((GEN)nf[5])[1],pg3);pgen4=cgetg(RU+1,17);
		for(l=1;l<=R1;l++) pgen4[l]=(long)glog((GEN)pgen5[l],PRECREG);
		for(l=R1+1;l<=RU;l++) pgen4[l]=lmul2n(glog((GEN)pgen5[l],PRECREG),1);
		for(l=1;l<=RU;l++) gaffect((GEN)(pgen4)[l],(GEN)((GEN)matarch[t])[l]);
		if(t>=nbrel) goto LABEL2;if(nbrelideal==nbrelpid) goto LABEL3;
	      }
	    }
	  }
	  avma=av3;
	}
	x[k]--;
      }
    }
    while(nonnul);
  LABEL3:
    avma=av2;
    if(DEBUGLEVEL>1)
    {
      fprintferr("temps pour cet ideal : ");fprintferr("%ld\n",timer2());
      flusherr();
    }
  }
LABEL2:
  avma=av1;
  if(DEBUGLEVEL>1)
  {
    fprintferr("\nLes elements de petite norme ont donne %ld relations.\n",t);
    fprintferr("\nCalcul du rang :");flusherr();avma=av1;matpro=cgetg(t+1,19);
    for(j=1;j<=t;j++)
    {
      pgen3=cgetg(KC+1,18);matpro[j]=(long)pgen3;
      for(i=1;i<=KC;i++) pgen3[i]=lstoi(mat[j][i]);
    }
    vecindep=(GEN)(indexrank(matpro)[2]);ra=lg(vecindep)-1;
    fprintferr("rang = %ld ; colonnes independantes :\n",ra);flusherr();
    for(i=1;i<=ra;i++) fprintferr("%4ld",itos((GEN)vecindep[i]));avma=av1;
    fprintferr("\ntemps du calcul du rang = %ld\n",timer2());flusherr();
  }
  if(DEBUGLEVEL)
  {
    fprintferr("\n");quotfact=((double)nbsmallfact)/((double)nbsmallnorm);
    fprintferr("nb. fact./nb. small norm = %ld/%ld = %f\n",nbsmallfact,nbsmallnorm,quotfact);
    if(DEBUGLEVEL<4)
    {fprintferr("Temps des relations de petite norme = %ld\n",timer2());}
    flusherr();
  }
  free(x);free(y);free(zz);free(vv);for(l=1;l<=N;l++) free(qq[l]);free(qq);
  return t;
}

#define MAXRELSUP (min(50,4*KC))

long
random_relation(long t,long q,long LIMC,long KC,long KCZ,long RU,long PRECREG,long PRECREGINT,GEN nf,long *firstime,long *secondtime,long *iz,long *jideal,long *bouext,long *bou,long *slim,long **ptvinvperm,long lgsub,GEN *idealbase,long *numideal,long *prim,long *numprim,long *vperm,long *subvperm,long **mat,GEN matarch,long **matcopy,GEN lmatt2,long nbmatt2,GEN mit,long *ex,GEN IDN,GEN **vp,GEN vectbase,long nrelsup,long nlze,long *phase,long *extrarel,long *cmptglob,long *nbtest,GEN *extramat,GEN *extramatc)
{
  long av,av1,i,j,k,pz,ss,fl,sinit,s,exsup;
  long mitcol,fpc,primfact[500],expoprimfact[500],ep,trouve,cptzer;
  long *p1,*vinvperm;
  GEN ideal,idealpro,ideal2,idealcourant,pgen1,pgen2;
  
  av1=avma;cptzer=0;s=t;
  if(((*phase)==1)&&(*firstime)&&(!(*secondtime)))
  {
    *secondtime=1;
    s=0;*slim=*extrarel=nlze;
  }
  else if((*phase)==1){s=0;*slim=(*extrarel)=nlze;}if((*phase)!=1) *iz=1;
  for(;s<(*slim);(*iz)++)
  {
    idealcourant=idealbase[numprim[pz=prim[(*iz-1)%KCZ+1]]];
    if((*phase)!=1) *jideal=1;
    for(;((*jideal)<lg(idealcourant))&&(s<(*slim));(*jideal)++) 
    {
      if(*phase)
      {
	if((*phase)==1)
	{
	  vinvperm=(long*)malloc(sizeof(long)*(KC+1));
	  *ptvinvperm=vinvperm;
	  for(i=1;i<=KC;i++) vinvperm[vperm[i]]=i;
	}
	ideal=idealmulprime(nf,vp[1][0],(GEN)idealcourant[*jideal]);exsup=1;
      } 
      else
      {ideal=idealmulprime(nf,vp[1][0],(GEN)idealcourant[*jideal]);exsup=1;}
      if(((*phase)!=1)||(*secondtime)) for(i=1;i<lgsub;i++) ex[i]=mymyrand()>>q;
      for (i=1;i<lgsub;i++) ideal=idealmulh(nf,ideal,(GEN)vp[i][ex[i]]);
      fl=gegal((GEN)ideal[1],IDN);av=avma;sinit=s;if((*phase)!=1) *bouext=1;
      for(;((*bouext)<=nbmatt2)&&(s==sinit)&&(s<(*slim));(*bouext)++)
      {
	if((*phase)!=1) *bou=*bouext;
	for(;((*bou)<=nbmatt2)&&(s==sinit)&&(s<(*slim));(*bou)++)
	{
	  if((*phase)==1)
	  {
	    (*bou)++;(*phase)++;mitcol=lg(mit)-1;
	    *nbtest=0;*extramat=cgetg((*extrarel)+1,19);
	    if(DEBUGLEVEL)
            fprintferr("\n(recherche de %ld relations supplementaires)\n",*extrarel);
	    *extramatc=cgetg((*extrarel)+1,19);
	    for(j=1;j<=*extrarel;j++) 
	    {
	      (*extramat)[j]=lgetg(KC+1,18);pgen1=cgetg(RU+1,18);
	      (*extramatc)[j]=(long)pgen1;
	      for(i=1;i<=RU;i++) 
	      {
		pgen2=cgetg(3,6);pgen1[i]=(long)pgen2;
		pgen2[1]=lgetr(PRECREG);pgen2[2]=lgetr(PRECREG);
	      }
	    }
	    av=avma;sinit=s;
	  }
	  if((*bou)<=nbmatt2)
	  {
	    if(DEBUGLEVEL>2)
	    {
	      fprintferr("iz=%ld,jideal=%ld,",*iz,*jideal);
	      fprintferr("bouext=%ld,bou=%ld,phase=%ld,",*bouext,*bou,*phase);
	      fprintferr("random=%ld\n",pari_randseed);
	      if(DEBUGLEVEL>3)
	      {fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
	      flusherr();
	    }
	    idealpro=ideallllredpart1spec(nf,(GEN)ideal[1],gcoeff(lmatt2,*bouext,*bou),1,PRECREGINT);
	    if(DEBUGLEVEL>5)
	    {fprintferr("idealpro = ");output(idealpro);flusherr();}
	    fpc=factorisegen(nf,idealpro,KCZ,LIMC,primfact,expoprimfact,prim,idealbase,numideal,numprim);
	    if(DEBUGLEVEL>1){fprintferr(".");flusherr();}
	    (*nbtest)++;
	    if (fpc&&(s<(*slim)))
	    {
	      ideal2=ideallllredpart2(nf,(GEN)ideal[2],(GEN)idealpro[2],PRECREG);
	      s++;
	      if(*phase)
	      {
		pgen1=(GEN)(*extramat)[s];for(i=1;i<=KC;i++) pgen1[i]=zero;
		for(i=1;i<lgsub;i++) pgen1[vinvperm[subvperm[i]]]=lstoi(-ex[i]);
		ss=vinvperm[numideal[pz]+(*jideal)];
		pgen1[ss]=laddsi(-exsup,(GEN)pgen1[ss]);
		
		for(i=1;i<=primfact[0];i++)
		{
		  k=vinvperm[primfact[i]];ep=expoprimfact[i];
		  pgen1[k]=laddsg(ep,(GEN)pgen1[k]);
		}
		(*cmptglob)++;p1=(long *)malloc(sizeof(long)*(KC+1));
		matcopy[*cmptglob]=p1;
		for(i=1;i<=KC;i++) p1[vperm[i]]=itos((GEN)pgen1[i]);
		if(gcmp0(pgen1)){s--;free(p1);(*cmptglob)--;}
		else 
		{
		  for(i=1;i<=RU;i++)
		    gaffect((GEN)ideal2[i],(GEN)((GEN)(*extramatc)[s])[i]);
		  trouve=already_found_relation_of_long(matcopy,*cmptglob,KC);
		  if(trouve) 
		  {
		    if(DEBUGLEVEL>1)
		    {
		      fprintferr("rel. elim. phase 2: %ld ",trouve);
		      if(DEBUGLEVEL>3)
		      {
			fprintferr("(iz=%ld,",*iz);
			fprintferr("jideal=%ld,bouext=%ld,",*jideal,*bouext);
			fprintferr("bou=%ld,phase=%ld)",*bou,*phase);
			if(DEBUGLEVEL>3)
			{fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
		      }
		      fprintferr("\n[");
		      for(i=1;i<=KC;i++) fprintferr("%4ld",matcopy[(*cmptglob)][i]);
		      fprintferr("]\n");flusherr();
		    }
		    pgen1=(GEN)(*extramat)[s];for(i=1;i<=KC;i++) pgen1[i]=zero;
		    s--;free(p1);(*cmptglob)--;cptzer++;
		    if(cptzer>MAXRELSUP){avma=av1;return -1;}
		  }
		  else
		  {
		    cptzer=0;
		    if(DEBUGLEVEL)
		    {
		      fprintferr("\n++++ cmptglob = %ld , s = %ld ",*cmptglob,s);
		      fprintferr(":nouvelle relation:");
		      if(DEBUGLEVEL>3)
		      {
			fprintferr("(iz=%ld,",*iz);
			fprintferr("jideal=%ld,bouext=%ld,",*jideal,*bouext);
			fprintferr("bou=%ld,phase=%ld)",*bou,*phase);
		      }
		      fprintferr("\n[");
		      for(i=1;i<=KC;i++)fprintferr("%4ld",matcopy[(*cmptglob)][i]);
		      fprintferr("]\n");
		      if(DEBUGLEVEL>3)
		      fprintferr("temps pour cette relation : %ld\n",timer2());
		      if(DEBUGLEVEL>6)
		      {
			fprintferr("partie archimedienne =\n");
			for(i=1;i<=RU;i++) outerr((GEN)((GEN)(*extramatc)[s])[i]);
			for(i=1;i<=KC;i++) fprintferr("vperm[%ld]=%ld,",i,vperm[i]);
			fprintferr("\n");flusherr();
			fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();
		      }
		    }
		  }
		}
	      }
	      else
	      {
		for(i=1;i<lgsub;i++) mat[s][subvperm[i]]= -ex[i];
		for(i=1;i<=primfact[0];i++) mat[s][primfact[i]]+=expoprimfact[i];
		mat[s][numideal[pz]+(*jideal)]--;
		for(i=1;i<=RU;i++) gaffect((GEN)ideal2[i],(GEN)((GEN)matarch[s])[i]);
		trouve=already_found_relation_of_long(mat,s,KC);
		if(trouve) 
		{
		  if(DEBUGLEVEL>1)
		  {
		    fprintferr("rel. elim. phase 0: %ld ",trouve);
		    if(DEBUGLEVEL>3)
		    {
		      fprintferr("(iz=%ld,",*iz);
		      fprintferr("jideal=%ld,bouext=%ld,",*jideal,*bouext);
		      fprintferr("bou=%ld,phase=%ld)",*bou,*phase);
		      if(DEBUGLEVEL>6)
		      {fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
		    }
		    fprintferr("\n[");
		    for(i=1;i<=KC;i++) fprintferr("%4ld",mat[s][i]);
		    fprintferr("]\n");flusherr();
		  }
		  for(i=1;i<=KC;i++) mat[s][i]=0;s--;cptzer++;
		  if(cptzer>MAXRELSUP){avma=av1;return -1;}
		}
		else
                {
		  if(DEBUGLEVEL)
		  {
		    fprintferr("\n++++ s = %ld :nouvelle relation:",s);
		    if(DEBUGLEVEL>3)
		    {
		      fprintferr("(iz=%ld,",*iz);
		      fprintferr("jideal=%ld,bouext=%ld,",*jideal,*bouext);
		      fprintferr("bou=%ld,phase=%ld)",*bou,*phase);
		    }
		    fprintferr("\n[");
		    for(i=1;i<=KC;i++)fprintferr("%4ld",mat[s][i]);
		    fprintferr("]\n");
		    if(DEBUGLEVEL>3)
		    fprintferr("temps pour cette relation : %ld\n",timer2());
		    if(DEBUGLEVEL>6)
		    {
		      fprintferr("partie archimedienne =\n");
		      for(i=1;i<=RU;i++) outerr((GEN)((GEN)(matarch)[s])[i]);
		      for(i=1;i<=KC;i++) fprintferr("vperm[%ld]=%ld,",i,vperm[i]);
		      fprintferr("\n");flusherr();
		      fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();
		    }
		  }
		  cptzer=0;
		}
	      }
	    }
	  }
	  if(fl) goto LABBUCH;
	}
      }
    LABBUCH:    
      if(((*phase)!=1)&&(sinit==s)&&(s<(*slim))) (*jideal)--;
      if((sinit==s)||(!(*phase))) avma=av;
    }
  }
  if(s>=(*slim))
  {
    if((*phase)!=2)
    {
      if((*iz)>1) (*iz)--;
      if((*jideal)>1) (*jideal)--;if((*bouext)>1) (*bouext)--;if((*bou)>1) (*bou)--;
    }
    else
    {
      *bouext=1;*bou=1;
      if(*jideal>=lg(idealcourant)) *jideal=1;
      else (*iz)--;
    }
  }
  return s;
}

#undef MAXRELSUP

long
already_found_relation_of_long(long **mat,long s,long KC)
{
  long non_trouve,l,ll,lll,cca,ccb,rr,*pts,*ptl,fff;
  non_trouve=1;l=s-1;
  while(non_trouve&&(l>=1))
  {
    lll=1;while((!mat[s][lll])&&(!mat[l][lll])) lll++;
    cca=mat[s][lll];ccb=mat[l][lll];
    if((!cca)||(!ccb)) l--;
    else
    {
      if(!(abs(cca)%abs(ccb))){rr=cca/ccb;pts=mat[s];ptl=mat[l];fff=1;}
      else if(!(abs(ccb)%abs(cca))){rr=ccb/cca;pts=mat[l];ptl=mat[s];fff=1;}
      else{fff=0;l--;}
      if(fff)
      {ll=lll+1;while((ll<=KC)&&(pts[ll]==(rr*ptl[ll]))) ll++;non_trouve=(ll<=KC);l--;}
    }
  }
  if(non_trouve) return 0;else return l+1;
}

void
desallocate(GEN primsubfactorbase,long cmptglob,GEN **vp,long *numprim,long *prim,long *numideal,GEN *idealbase,long *ex,long **matcopy,long *vperm,long *subvperm,long *vinvperm,long flvp,long phase)
{
  long i;
  free(numprim);free(prim);free(numideal);free(idealbase);free(ex);
  free(vperm);free(subvperm);
  if(flvp){for(i=1;i<lg(primsubfactorbase);i++) free(vp[i]);free(vp);}
  if(phase)
  {
    for(i=1;i<=cmptglob;i++) free(matcopy[i]);free(matcopy);
  }
}

GEN
compute_regulator(GEN mit,long col,GEN matarch,long RU,long R1,long N,long PRECREG,GEN *xarch,GEN *xreal,GEN *image_mdet,long *ptsreg)
{
  long sizeofmit,sreg,i,j,av,tetpil,dec;
  GEN vei,mdet,reg;
 
  av=avma;sizeofmit=lg(mit)-1;sreg=col-sizeofmit;
  *xarch=cgetg(sreg+1,19);for(j=1;j<=sreg;j++) (*xarch)[j]=matarch[j];
  (*xreal)=greal(*xarch);
  vei=cgetg(RU+1,18);for(i=1;i<=RU;i++) vei[i]=(i<=R1)?un:deux;
  mdet=cgetg(sreg+2,19);for(j=2;j<=sreg+1;j++) mdet[j]=(*xreal)[j-1];
  mdet[1]=(long)vei;
  if(DEBUGLEVEL)
  {
    fprintferr("temps trivialites xarch, xreal, mdet: ");
    fprintferr("%ld\n",timer2());flusherr();
  }
  *image_mdet=imagereel(mdet,PRECREG);
  if(DEBUGLEVEL)
  {
    fprintferr("temps imagereel: ");fprintferr("%ld\n",timer2());
    flusherr();
  }
  if(lg(*image_mdet)!=(RU+1)) affsr(0,reg=cgetr(PRECREG));
  else reg=gdivgs(gabs(detreel(*image_mdet),PRECREG),N);
  if(DEBUGLEVEL)
  {
    fprintferr("temps detreel: ");fprintferr("%ld\n",timer2());
    flusherr();
  }
  tetpil=avma;*xarch=gcopy(*xarch);*xreal=gcopy(*xreal);*image_mdet=gcopy(*image_mdet);
  reg=gcopy(reg);dec=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;
  (*xarch)+=dec;(*xreal)+=dec;(*image_mdet)+=dec;reg+=dec;
  *ptsreg=sreg;return reg;
}

GEN
compute_class_number(GEN mit,GEN *met,GEN *u1)
{
  long sizeofmit,i,av,tetpil,dec;
  GEN pgen1,clh,u1u2,u2;

  sizeofmit=lg(mit)-1;av=avma;
  if(sizeofmit)
  {
    u1u2=smith2(mit);*u1=(GEN)u1u2[1];u2=(GEN)u1u2[2];
    pgen1=gmul(*u1,gmul(mit,u2));*u1=ginv(*u1);*met=pgen1;
    clh=gun;for(i=1;i<=sizeofmit;i++) clh=mulii(clh,gcoeff(pgen1,i,i));
  }
  else{clh=gun;*met=cgetg(1,19);*u1=cgetg(1,19);}
  tetpil=avma;*met=gcopy(*met);*u1=gcopy(*u1);clh=gcopy(clh);
  dec=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;
  (*met)+=dec;(*u1)+=dec;clh+=dec;return clh;
}

GEN
compute_check(GEN mit,long col,GEN image_mdet,GEN xreal,long RU,GEN *reg,GEN z,GEN zu,GEN clh,long PRECREG,GEN *parch)
{
  long sizeofmit,sreg,i,j,av,tetpil,dec;
  GEN lambda,sublambda_1,pgen2,pgen3,c_1,den;

  sizeofmit=lg(mit)-1;sreg=col-sizeofmit;av=avma;
  lambda=invmulmatreel(image_mdet,xreal);sublambda_1=cgetg(sreg+1,19);
  for(i=1;i<=sreg;i++)
  {
    pgen2=cgetg(RU,18);sublambda_1[i]=(long)pgen2;pgen3=(GEN)lambda[i];
    for(j=1;j<RU;j++) pgen2[j]=pgen3[j+1];
  }
  if(DEBUGLEVEL)
  {
    fprintferr("temps trivialites lambda: ");fprintferr("%ld\n",timer2());
    flusherr();
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
  }
  c_1=gmul2n(gdiv(*reg,gdiv(gmul(z,(GEN)zu[1]),clh)),1);
  sublambda_1=bestappr(sublambda_1,c_1);den=denom(sublambda_1);
  if(gcmp(den,c_1)>0){avma=av;return gzero;}
  *parch=(GEN)lllkerim(gmul(sublambda_1,den))[2];
  *reg=gabs(gmul(*reg,det(gmul(sublambda_1,*parch))),PRECREG);
  c_1=gmul2n(gdiv(*reg,gdiv(gmul(z,(GEN)zu[1]),clh)),1);
  if(DEBUGLEVEL)
  {
    fprintferr("temps bestappr/regulateur: ");fprintferr("%ld\n",timer2());
    flusherr();
  }
  tetpil=avma;*reg=gcopy(*reg);*parch=gcopy(*parch);c_1=gcopy(c_1);
  dec=lpile(av,tetpil,0)>>TWOPOTBYTES_IN_LONG;
  (*reg)+=dec;(*parch)+=dec;c_1+=dec;return c_1;
}

long
be_honest(long KCZ,long KCZ2,long RU,long PRECREGINT,GEN IDN,long s,GEN nf,long lgsub,GEN *idealbase,long *numideal,long *prim,long *numprim,long q,GEN **vp,long *ex,GEN exu)
{
  long av,av1,i,j,iz1,pz,sinit,nbtest,bouh,fpc,primfact[500],expoprimfact[500];
  GEN pgen1,ideal,idealpro;

  av=avma;
  for(iz1=KCZ+1;iz1<=KCZ2;iz1++)
  {
    pgen1=idealbase[numprim[pz=prim[iz1]]];
    if(DEBUGLEVEL>1) fprintferr(" %ld,",pz);
    for(j=1;j<lg(pgen1);j++) 
    {
      sinit=s;nbtest=0;
      do
      {
	for(i=1;i<lgsub;i++) ex[i]=mymyrand()>>q;av1=avma;
	ideal=idealmulprime(nf,IDN,(GEN)pgen1[j]);
	for (i=1;i<lgsub;i++) ideal=idealmulh(nf,ideal,(GEN)vp[i][ex[i]]);
	for(bouh=1;(bouh<=RU)&&(s==sinit);bouh++)
	{
	  if(bouh==1) {for(i=1;i<=RU;i++) exu[i]=lstoi(mymyrand()>>q);}
	  else{for(i=1;i<=RU;i++) exu[i]=zero;exu[bouh]=lstoi(10);}
	  idealpro=ideallllredpart1(nf,ideal,exu,1,PRECREGINT);nbtest++;
	  fpc=factorisegen(nf,idealpro,iz1-1,prim[iz1-1],primfact,expoprimfact,prim,idealbase,numideal,numprim);
	  if(fpc) s++;
	}
	avma=av1;
      }
      while((sinit==s)&&(nbtest<20));
      if(sinit==s){avma=av;return 0;}
    }
  }
  avma=av;return 1;
}

GEN
class_group_generators(long N,long RU,long PRECREG,long PRECREGINT,GEN nf,GEN met,GEN clh,GEN u1,GEN vectbase,long *vperm)
{
  long av,tetpil,i,j,c,lo,s,expomax;
  GEN res,basecl,p3,pgen1,pgen2,mot,detp3,p3min,orderp3,detpgen1;
  
  av=avma;lo=lg(met)-1;
  c=0;for(i=1;i<=lo;i++) if(!gcmp1(gcoeff(met,i,i))) c++;basecl=cgetg(c+1,17);
  if(DEBUGLEVEL)
  {
    fprintferr("temps trivialites cleancol: ");fprintferr("%ld\n",timer2());
    flusherr();
  }
  for(j=1;j<=c;j++)
  {
    p3=(GEN)idealpowprime(nf,(GEN)vectbase[vperm[1]],pgen1=gcoeff(u1,1,j),PRECREG)[1];
    if(signe(pgen1)<0) p3=numer(p3);
    for(i=2;i<=lo;i++)
    {
      pgen1=gcoeff(u1,i,j);s=signe(pgen1);
      if(s)
      {
	if(s>0)
	  p3=idealmulh(nf,p3,(GEN)idealpowprime(nf,(GEN)vectbase[vperm[i]],pgen1,PRECREGINT)[1]);
	else
	  p3=idealmulh(nf,p3,numer((GEN)idealpowprime(nf,(GEN)vectbase[vperm[i]],pgen1,PRECREGINT)[1]));
	p3=ideallllredall(nf,p3,gzero,PRECREGINT,PRECREGINT);
      }
    }
    detp3=det(p3);p3min=p3;pgen1=p3;
    pgen2=numer(idealinv(nf,pgen1));
    detpgen1=det(pgen2);
    if(gcmp(detpgen1,detp3)<0){detp3=detpgen1;p3min=pgen2;}
    orderp3=gcoeff(met,j,j);
    expomax=itos(gmin(orderp3,stoi(50)));
    for(i=2;i<expomax;i++)
    {
      pgen1=idealmulred(nf,pgen1,p3,PRECREGINT);
      if(gcmp1(mppgcd(stoi(i),orderp3)))
      {
	detpgen1=det(pgen1);
	if(gcmp(detpgen1,detp3)<0){detp3=detpgen1;p3min=pgen1;}
	pgen2=numer(idealinv(nf,pgen1));
	detpgen1=det(pgen2);
	if(gcmp(detpgen1,detp3)<0){detp3=detpgen1;p3min=pgen2;}
      }
    }
    basecl[j]=(long)p3min;
  }
  mot=cgetg(c+1,17);for(i=1;i<=c;i++) mot[i]=coeff(met,i,i);
  pgen1=cgetg(4,17);res=pgen1;pgen1[1]=(long)clh;
  pgen1[2]=(long)mot;pgen1[3]=(long)basecl;tetpil=avma;
  return gerepile(av,tetpil,gcopy(res));
}

GEN
compute_matt2(long RU,GEN nf,long *nbmatt2)
{
  long av,tetpil,i,j,k;
  GEN lmatt2,pgen1,pgen2,pgen3,pgen4;

  if(DEBUGLEVEL)
  {
    fprintferr("Calcul des matrices T2 avec poids : ");
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
    flusherr();
  }
  if(RU<10) *nbmatt2=RU;else *nbmatt2=9;av=avma;
  lmatt2=cgetg((*nbmatt2)+1,19);pgen1=(GEN)nf[5];pgen3=(GEN)pgen1[2];
  for(j=1;j<=(*nbmatt2);j++)
  {
    pgen2=cgetg(*nbmatt2+1,18);lmatt2[j]=(long)pgen2;
    for(i=1;i<=(*nbmatt2);i++)
    {
      if(j<i) pgen2[i]=zero;
      else
      {
	pgen4=cgetg(RU+1,19);for(k=1;k<=RU;k++) pgen4[k]=pgen3[k];
	pgen4[i]=lmul2n((GEN)pgen3[i],20);if(i!=j) pgen4[j]=lmul2n((GEN)pgen3[j],20);
	pgen2[i]=lreal(gmul(pgen4,(GEN)pgen1[1]));
      }
    }
  }
  if(DEBUGLEVEL)
  {
    fprintferr("temps matrices T2 avec poids : ");fprintferr("%ld\n",timer2());
    if(DEBUGLEVEL>7){fprintferr("$$$$$ AVMA = %ld\n",avma);flusherr();}
    flusherr();
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(lmatt2));
}

GEN
buchnarrow(GEN bignf)
{
  GEN nf,dataclassgroup,zeromod2,unmod2,rac,cyclic,generator,v,matsign,genplus;
  GEN dataunit,p1,p2,p3,p4,h,v2,vecsign,clh,basecl,met,u1,mot,id,v3,alpha,bas,lambda;
  long R1,R,i,j,va,ngen,sizeh,t,lo,c,s,av=avma,tetpil,N,fl,r,rr,rankinit,ranknew,k,limr;
  long av1,kk;

  if((typ(bignf)!=17)||(lg(bignf)!=9))
    err(talker,"not a big number field vector in buchnarrow");
  nf=(GEN)bignf[7];
  if((typ(nf)!=17)||(lg(nf)<10)) err(idealer1);
  R1=itos((GEN)((GEN)nf[2])[1]);if(!R1) return gcopy((GEN)((GEN)bignf[8])[1]);
  N=lgef((GEN)nf[1])-3;id=idmat(N);
  unmod2=gmodulcp(gun,gdeux);zeromod2=gmodulcp(gzero,gdeux);
  va=varn((GEN)nf[1]);rac=(GEN)nf[6];
  dataclassgroup=(GEN)((GEN)bignf[8])[1];
  cyclic=(GEN)dataclassgroup[2];generator=(GEN)dataclassgroup[3];ngen=lg(cyclic)-1;
  v=cgetg(R1+1,18);for(i=1;i<=R1;i++) v[i]=(long)unmod2;
  matsign=signunit(bignf);R=lg(matsign);
  dataunit=cgetg(R+1,19);
  for(j=1;j<R;j++)
  {
    p1=cgetg(R1+1,18);dataunit[j]=(long)p1;
    for(i=1;i<=R1;i++)
    {p1[i]=(signe(gcoeff(matsign,i,j))>0)? (long)zeromod2:(long)unmod2;}
  }
  dataunit[R]=(long)v;
  v2=image(dataunit);t=lg(v2)-1;
  vecsign=cgetg(R1+1,18);sizeh=ngen+R1-t;
  genplus=cgetg(sizeh+1,18);for(i=1;i<=ngen;i++) genplus[i]=generator[i];
  r=1;rr=3;bas=(GEN)nf[7];rankinit=t;j=ngen;fl=(t<R1);
  lambda=cgeti(N+1);
  while(fl)
  {
    limr=(itos(gmin(stoi(BIGINT),gpuigs(stoi(rr),N)))-1)>>1;
    for(k=rr;(k<=limr)&&fl;k++)
    {
      kk=k;for(i=1;i<=N;i++){lambda[i]=(kk+r)%rr-r;kk/=rr;}
      av1=avma;
      alpha=gzero;for(i=1;i<=N;i++) alpha=gadd(alpha,gmulsg(lambda[i],(GEN)bas[i]));
      for(i=1;i<=R1;i++)
	vecsign[i]=(signe(gsubst(alpha,va,(GEN)rac[i]))>0)?(long)zeromod2:(long)unmod2;
      v3=concat(v2,gtomat(vecsign));ranknew=rank(v3);
      if(ranknew>rankinit)
      {
	v2=v3;rankinit++;fl=(rankinit<R1);j++;genplus[j]=(long)alpha;
      }
      else avma=av1;
    }
    if(fl){r++;rr+=2;}
  }
  generator=genplus;
  v=ginv(v2);
  h=cgetg(sizeh+1,19);
  for(j=1;j<=ngen;j++)
  {
    p1=cgetg(sizeh+1,18);h[j]=(long)p1;
    p4=(GEN)cyclic[j];
    p2=(GEN)generator[j];
    p3=gmul((GEN)nf[7],(GEN)isprincipalgen(bignf,idealpow(nf,p2,p4))[2]);
    for(i=1;i<=R1;i++)
      vecsign[i]=(signe(gsubst(p3,va,(GEN)rac[i]))>0)?(long)zeromod2:(long)unmod2;
    p3=gmul(v,vecsign);
    for(i=1;i<=ngen;i++) p1[i]=(i==j)?cyclic[j]:zero;
    for(i=ngen+1;i<=sizeh;i++) p1[i]=(long)lift((GEN)p3[i+t-ngen]);
  }
  for(j=ngen+1;j<=sizeh;j++)
  {
    p1=cgetg(sizeh+1,18);h[j]=(long)p1;
    for(i=1;i<=sizeh;i++) p1[i]=(i==j)?deux:zero;
  }
  clh=compute_class_number(h,&met,&u1);
  u1=reducemodmatrix(u1,h);
  lo=lg(met)-1;
  c=0;for(i=1;i<=lo;i++) if(!gcmp1(gcoeff(met,i,i))) c++;basecl=cgetg(c+1,17);
  for(j=1;j<=c;j++)
  {
    p3=(GEN)idealpow(nf,(GEN)generator[1],p1=gcoeff(u1,1,j));
    if(signe(p1)<0) p3=numer(p3);
    for(i=2;i<=lo;i++)
    {
      p1=gcoeff(u1,i,j);s=signe(p1);
      if(s)
      {
	p3=idealmul(nf,p3,(GEN)idealpow(nf,(GEN)generator[i],p1));
	p3=gdiv(p3,content(p3));
      }
    }
    basecl[j]=(long)p3;
  }
  tetpil=avma;p1=cgetg(4,17);p1[1]=lcopy(clh);
  mot=cgetg(c+1,17);for(i=1;i<=c;i++) mot[i]=lcopy(gcoeff(met,i,i));
  p1[2]=(long)mot;p1[3]=lcopy(basecl);
  return gerepile(av,tetpil,p1);
}

GEN
buchrayall(GEN bignf,GEN ideal,long flun)
{
  GEN nf,dataclassgroup,cyclic,generator,genplus;
  GEN dataunit,p1,p2,p3,p4,p5,h,clh,basecl,met,u1,mot,id;
  GEN racunit,bigres,bid,resbid,resbid2,resbid3,x,y,funits,hmat,vecel;
  long R1,R,R3,i,j,ngen,sizeh,lo,c,s,av=avma,tetpil,N;
  long flid;

  if((typ(bignf)!=17)||(lg(bignf)!=9))
    err(talker,"not a big number field vector in buchray");
  nf=(GEN)bignf[7];bigres=(GEN)bignf[8];
  if(lg(bigres)!=7)
    err(talker,"buchray needs the units, please apply buchinitfu instead of buchinit");
  if((typ(nf)!=17)||(lg(nf)<10)) err(idealer1);
  dataclassgroup=(GEN)bigres[1];
  cyclic=(GEN)dataclassgroup[2];generator=(GEN)dataclassgroup[3];
  ngen=lg(cyclic)-1;
  R1=itos((GEN)((GEN)nf[2])[1]);N=lgef((GEN)nf[1])-3;id=idmat(N);
  bid=zidealstarinitold(nf,ideal);
  if((typ(ideal)==17)&&(lg(ideal)==3)) {flid=1;x=(GEN)ideal[1];}
  else {flid=0;x=ideal;}
  x=idealhermite(nf,x);
  resbid=(GEN)bid[2];resbid2=(GEN)resbid[2];resbid3=(GEN)resbid[3];
  R3=lg(resbid2)-1;
  if(!R3)
  {
    if(flun>0){avma=av;return gcopy((GEN)bigres[1]);}
    tetpil=avma;
    y=cgetg(6,17);y[1]=lcopy(bignf);y[2]=lcopy(bid);
    p1=cgetg(ngen+1,17);y[3]=(long)p1;
    for(i=1;i<=ngen;i++) p1[i]=un;y[4]=(long)idmat(ngen);y[5]=lcopy((GEN)bigres[1]);
    return gerepile(av,tetpil,y);
  }
  funits=(GEN)bigres[5];R=lg(funits);
  if(R!=((N+R1)>>1))
    err(talker,"incorrect unit length in buchray");
  racunit=(GEN)((GEN)bigres[4])[2];
  dataunit=cgetg(R+R3+1,19);
  dataunit[1]=(long)zideallog(nf,racunit,bid);
  for(j=2;j<=R;j++) dataunit[j]=(long)zideallog(nf,(GEN)funits[j-1],bid);
  for(j=R+1;j<=R+R3;j++)
  {
    p1=cgetg(R3+1,18);dataunit[j]=(long)p1;
    for(i=1;i<=R3;i++) p1[i]=(i==(j-R))?resbid2[i]:zero;
  }
  hmat=hnfmod(dataunit,(GEN)resbid[1]);
  sizeh=ngen+R3;h=cgetg(sizeh+1,19);
  genplus=cgetg(sizeh+1,17);
  vecel=cgetg(ngen+1,17);
  for(j=1;j<=ngen;j++)
  {
    p1=cgetg(sizeh+1,18);h[j]=(long)p1;
    p4=(GEN)cyclic[j];p2=(GEN)generator[j];
    vecel[j]=(long)idealcoprime(nf,p2,x);
    p3=idealmul(nf,(GEN)vecel[j],p2);genplus[j]=(long)p3;
    p5=(GEN)isprincipalgen(bignf,idealpow(nf,p3,p4))[2];
    p2=zideallog(nf,p5,bid);
    for(i=1;i<=ngen;i++) p1[i]=(i==j)?cyclic[j]:zero;
    for(i=ngen+1;i<=sizeh;i++) p1[i]=lnegi((GEN)p2[i-ngen]);
  }
  for(j=ngen+1;j<=sizeh;j++)
  {
    p1=cgetg(sizeh+1,18);h[j]=(long)p1;
    for(i=1;i<=ngen;i++) p1[i]=zero;
    for(i=ngen+1;i<=sizeh;i++) p1[i]=coeff(hmat,i-ngen,j-ngen);
    genplus[j]=resbid3[j-ngen];
  }
  generator=genplus;
  clh=compute_class_number(h,&met,&u1);
  u1=reducemodmatrix(u1,h);
  lo=lg(met)-1;
  c=0;for(i=1;i<=lo;i++) if(!gcmp1(gcoeff(met,i,i))) c++;basecl=cgetg(c+1,17);
  for(j=1;j<=c;j++)
  {
    p3=idmat(N);p4=p3;
    for(i=1;i<=lo;i++)
    {
      p1=gcoeff(u1,i,j);s=signe(p1);
      if(s)
      {
	if(s>0)	p3=idealmul(nf,p3,(GEN)idealpow(nf,(GEN)generator[i],p1));
	else p4=idealmul(nf,p4,(GEN)idealpow(nf,(GEN)generator[i],negi(p1)));
      }
    }
    p5=(GEN)ideleaddone(nf,p4,ideal)[1];
    basecl[j]=(long)idealdivexact(nf,idealmul(nf,p5,p3),p4);
  }
  tetpil=avma;p1=cgetg(4,17);p1[1]=lcopy(clh);
  mot=cgetg(c+1,17);for(i=1;i<=c;i++) mot[i]=lcopy(gcoeff(met,i,i));
  p1[2]=(long)mot;p1[3]=lcopy(basecl);
  if(flun>0) return gerepile(av,tetpil,p1);
  tetpil=avma;
  y=cgetg(6,17);y[1]=lcopy(bignf);y[2]=lcopy(bid);y[3]=lcopy(vecel);
  y[4]=(long)ginv(u1);y[5]=lcopy(p1);
  return gerepile(av,tetpil,y);
}

GEN
buchrayinit(GEN bignf, GEN ideal)
{
  return buchrayall(bignf,ideal,-1);
}

GEN
buchray(GEN bignf, GEN ideal)
{
  return buchrayall(bignf,ideal,1);
}

GEN
isprincipalrayall(GEN bignfray, GEN x, long flall)
{
  long av=avma,tetpil,i,c,N,ngen,ngzk;
  GEN bnf,nf,bid,vecel,vecep,matu,ep,generator,p1,p2,p3,beta,idep,y,rayclass;
  GEN divray,genray;

  if((typ(bignfray)!=17)||(lg(bignfray)!=6))
    err(talker,"incorrect bignfray in isprincipalray");
  bnf=(GEN)bignfray[1];bid=(GEN)bignfray[2];vecel=(GEN)bignfray[3];matu=(GEN)bignfray[4];
  rayclass=(GEN)bignfray[5];nf=(GEN)bnf[7];ngen=lg(vecel)-1;
  idep=isprincipalall(bnf,x,1);
  if(lg((GEN)idep[1])!=(ngen+1))
    err(talker,"incorrect generator length in isprincipalray");
  N=lgef((GEN)nf[1])-3;p2=cgetg(N+1,18);p2[1]=un;for(i=2;i<=N;i++) p2[i]=zero;
  ep=(GEN)idep[1];
  for(i=1;i<=ngen;i++)
    if(typ((GEN)vecel[i])>1)
      p2=element_mul(nf,p2,element_pow(nf,(GEN)vecel[i],(GEN)ep[i]));
  beta=element_div(nf,(GEN)idep[2],p2);
  p3=zideallog(nf,beta,bid);ngzk=lg(p3)-1;
  vecep=cgetg(ngen+ngzk+1,18);
  for(i=1;i<=ngen;i++) vecep[i]=ep[i];for(i=ngen+1;i<=ngen+ngzk;i++) vecep[i]=p3[i-ngen];
  p1=gmul(matu,vecep);
  p3=(GEN)bid[2];generator=(GEN)p3[3];
  divray=(GEN)rayclass[2];genray=(GEN)rayclass[3];c=lg(divray)-1;
  tetpil=avma;y=cgetg(c+1,18);
  for(i=1;i<=c;i++) y[i]=(long)modii((GEN)p1[i],(GEN)divray[i]);
  if(!flall) return gerepile(av,tetpil,y);
  else
  {
    p1=idmat(N);
    for(i=1;i<=c;i++) p1=idealmul(nf,idealpow(nf,(GEN)genray[i],(GEN)y[i]),p1);
    p2=isprincipalall(bnf,idealdiv(nf,x,p1),1);
    if(!gcmp0((GEN)p2[1])) err(talker,"bug in isprincipalray");
    tetpil=avma;
    p1=cgetg(4,17);p1[1]=lcopy(y);p1[2]=lcopy((GEN)p2[2]);
    p1[3]=(long)gmin((GEN)idep[3],(GEN)p2[3]);
    return gerepile(av,tetpil,p1);
  }
}
  
GEN
isprincipalray(GEN bignfray, GEN x)
{
  return isprincipalrayall(bignfray,x,0);
}

GEN
isprincipalraygen(GEN bignfray, GEN x)
{
  return isprincipalrayall(bignfray,x,1);
}

GEN
computezimmertbound(long N,long R1,GEN DK) /* DK = |dK| */
{
  long av,tetpil,i,R2;
  double **c;
  GEN w,p1,minkowskibound;
  
  if(N<2) return gun;
  else if(N<21)
  {
    av=avma;
    c=(double**)malloc(sizeof(double*)*21);
    for(i=1;i<=20;i++) c[i]=(double*)malloc(sizeof(double)*21);
    c[2][2]=-0.6931;c[2][0]=-0.45158;
    c[3][3]=-1.71733859;c[3][1]=-1.37420604;
    c[4][4]=-2.91799837;c[4][2]=-2.50091538;c[4][0]=-2.11943331;
    c[5][5]=-4.22701425;c[5][3]=-3.75471588;c[5][1]=-3.31196660;
    c[6][6]=-5.61209925;c[6][4]=-5.09730381;c[6][2]=-4.60693851;c[6][0]=-4.14303665;
    c[7][7]=-7.05406203;c[7][5]=-6.50550021;c[7][3]=-5.97735406;c[7][1]=-5.47145968;
    c[8][8]=-8.54052636;c[8][6]=-7.96438858;c[8][4]=-7.40555445;c[8][2]=-6.86558259;
    c[8][0]=-6.34608077;
    c[9][9]=-10.0630022;c[9][7]=-9.46382812;c[9][5]=-8.87952524;c[9][3]=-8.31139202;
    c[9][1]=-7.76081149;
    c[10][10]=-11.6153797;c[10][8]=-10.9966020;c[10][6]=-10.3907654;c[10][4]=-9.79895170;
    c[10][2]=-9.22232770;c[10][0]=-8.66213267;
    c[11][11]=-13.1930961;c[11][9]=-12.5573772;c[11][7]=-11.9330458;c[11][5]=-11.3210061;
    c[11][3]=-10.7222412;c[11][1]=-10.1378082;
    c[12][12]=-14.7926394;c[12][10]=-14.1420915;c[12][8]=-13.5016616;
    c[12][6]=-12.8721114;c[12][4]=-12.2542699;c[12][2]=-11.6490374;c[12][0]=-11.0573775;
    c[13][13]=-16.4112395;c[13][11]=-15.7475710;c[13][9]=-15.0929680;
    c[13][7]=-14.4480777;c[13][5]=-13.8136054;c[13][3]=-13.1903162;c[13][1]=-12.5790381;
    c[14][14]=-18.0466672;c[14][12]=-17.3712806;c[14][10]=-16.7040780;
    c[14][8]=-16.0456127;c[14][6]=-15.3964878;c[14][4]=-14.7573587;
    c[14][2]=-14.1289364;c[14][0]=-13.5119848;
    c[15][15]=-19.6970961;c[15][13]=-19.0111606;c[15][11]=-18.3326615;
    c[15][9]=-17.6620757;
    c[15][7]=-16.9999233;c[15][5]=-16.3467686;c[15][3]=-15.7032228;c[15][1]=-15.0699480;
    c[16][16]=-21.3610081;c[16][14]=-20.6655103;c[16][12]=-19.9768082;
    c[16][10]=-19.2953176;
    c[16][8]=-18.6214885;c[16][6]=-17.9558093;c[16][4]=-17.2988108;c[16][2]=-16.6510652;
    c[16][0]=-16.0131906;
    c[17][17]=-23.0371259;c[17][15]=-22.3329066;c[17][13]=-21.6349299;
    c[17][11]=-20.9435607;
    c[17][9]=-20.2591899;c[17][7]=-19.5822454;c[17][5]=-18.9131878;c[17][3]=-18.2525157;
    c[17][1]=-17.6007672;
    c[18][18]=-24.7243611;c[18][16]=-24.0121449;c[18][14]=-23.3056902;
    c[18][12]=-22.6053167;
    c[18][10]=-21.9113705;c[18][8]=-21.2242247;c[18][6]=-20.5442836;c[18][4]=-19.8719830;
    c[18][2]=-19.2077941;c[18][0]=-18.5522234;
    c[19][19]=-26.4217792;c[19][17]=-25.7021950;c[19][15]=-24.9879497;
    c[19][13]=-24.2793271;
    c[19][11]=-23.5766321;c[19][9]=-22.8801952;c[19][7]=-22.1903709;c[19][5]=-21.5075437;
    c[19][3]=-20.8321263;c[19][1]=-20.1645647;
    c[20][20]=-28.1285704;c[20][18]=-27.4021674;c[20][16]=-26.6807314;
    c[20][14]=-25.9645140;
    c[20][12]=-25.2537867;c[20][10]=-24.5488420;c[20][8]=-23.8499943;
    c[20][6]=-23.1575823;
    c[20][4]=-22.4719720;c[20][2]=-21.7935548;c[20][0]=-21.1227537;
    w=gexp(dbltor(c[N][R1]),6);
    for(i=1;i<=20;i++) free(c[i]);free(c);
    p1=gmul(gsqrt(DK,6),w);
    tetpil=avma;return gerepile(av,tetpil,ground(p1));
  }
  else
  {
    av=avma;R2=(N-R1)>>1;
    p1=gdiv(mpfact(N),gpuigs(stoi(N),N));
    minkowskibound=ground(gmul(gmul(p1,gpuigs(gdivsg(4,mppi(6)),R2)),gsqrt(DK,6)));
    if(cmpis(minkowskibound,500000)>0)
    {
      err(talker,"***** The field has degree more than 20 and the Minkowski bound is larger than 500 000 : it is unrealistic to certify this field *****");
      return gnil;
    }
    else
    {
      tetpil=avma;return gerepile(av,tetpil,gcopy(minkowskibound));
    }
  }
}

long
isinfactorbase(GEN factorbase,GEN ideal) /* l'ideal est premier sous forme 5 composantes
					    et la factor base est ordonnee */
{
  long av,i,j,lgfb;
  GEN p,alpha;

  av=avma;lgfb=lg(factorbase)-1;
  p=(GEN)ideal[1];alpha=(GEN)ideal[2];
  i=1;
  while((i<=lgfb)&&(gcmp(p,(GEN)((GEN)factorbase[i])[1])>0))
  {
    i++;
  }
  if(i>lgfb){avma=av;return 0;}
  else
  {
    if(!gegal(p,(GEN)((GEN)factorbase[i])[1])){avma=av;return 0;}
    else
    {
      if(gegal(alpha,(GEN)((GEN)factorbase[i])[2])){avma=av;return 1;}
      else
      {
	j=i+1;
	while((j<=lgfb)&&(gegal(p,(GEN)((GEN)factorbase[j])[1])))
	{
	  if(gegal(alpha,(GEN)((GEN)factorbase[j])[2])){avma=av;return 1;}
	  else j++;
	}
	avma=av;return 0;
      }
    }
  }
}

long
testprimeuptominkowskibound(GEN bnf,GEN minkowskibound)
{
  byteptr delta;
  long av,pp,i,nbideal,N,lgfb,fl;
  GEN nf,f,p1,vectpp,different,factorbase,pmax,dK;
  
  if(DEBUGLEVEL)
  {
   fprintferr("\nPHASE 1 : prime up to Zimmert bound are in computed class group ?\n\n");
   flusherr();
  }
  av=avma;
  nf=(GEN)bnf[7];N=lgef((GEN)nf[1])-3;f=(GEN)nf[4];dK=(GEN)nf[3];
  if(!gcmp1(f))
  {
    different=(GEN)((GEN)nf[5])[5];
    if(DEBUGLEVEL)
    {
      fprintferr("***** Testing the Different = ");outerr(different);flusherr();
    }
    p1=isprincipalall(bnf,different,0);
    if(DEBUGLEVEL>1)
    {
      fprintferr("      The answer is ");outerr(p1);flusherr();
    }
  }
  avma=av;
  factorbase=(GEN)bnf[5];lgfb=lg(factorbase)-1;pmax=(GEN)((GEN)factorbase[lgfb])[1];
  delta=diffptr;pp=*delta++;
  while(cmpsi(pp,minkowskibound)<1)
  {
    if(DEBUGLEVEL)
    {
      fprintferr("*** Testing prime p = %ld\n",pp);flusherr();
    }
    vectpp=primedec(bnf,stoi(pp));nbideal=lg(vectpp)-1;
    if(gdivise(dK,stoi(pp))) fl=1;else fl=0;
    for(i=1;i<(nbideal+fl);i++)
    {
      if(DEBUGLEVEL)
      {
	fprintferr("..... Testing ideal P = ");outerr((GEN)vectpp[i]);flusherr();
      }
      if((cmpis((GEN)((GEN)vectpp[i])[4],N))&&(cmpii(idealnorm(bnf,(GEN)vectpp[i]),minkowskibound)<1))
      {
	if((cmpsi(pp,pmax)<=0)&&(isinfactorbase(factorbase,(GEN)vectpp[i])))
	{
	  if(DEBUGLEVEL>1)
	  {
	    fprintferr("      The prime ideal is in the factor base\n");flusherr();
	  }
	}
	else
	{
	  p1=isprincipalall(bnf,(GEN)vectpp[i],0);
	  if(DEBUGLEVEL>1)
	  {
	    fprintferr("      The answer is ");outerr(p1);flusherr();
	  }
	}
/* Il vaudrait mieux attraper le message d'erreur en cas d'echec !! */
      }
      else
      {
	if(DEBUGLEVEL>1)
	{
	  fprintferr("      Norm ideal > Zimmert bound\n");flusherr();
	}
      }
    }
    pp+=*delta++;
  }
  avma=av;
  if(DEBUGLEVEL)
  {
    fprintferr("\nEnd of PHASE 1.\n\n\n");
    flusherr();
  }
  return 1;
}

GEN
mynfshanks(GEN nf,GEN x,GEN g0,GEN pr,GEN prhall)
/* rend le plus petit entier positif n tel que g0^n=x=1 modulo pr */
{
  long av,tetpil,i,N,lbaby,fl;
  GEN pf1,pfqr,p1,smalltable,p2,giant,smalltable2,permtable,pp,unmodp,p3,kgen;

  av=avma;
  pp=(GEN)pr[1];pf1=addsi(-1,gpui((GEN)pr[1],(GEN)pr[4],0));pfqr=racine(pf1);
  unmodp=gmodulcp(gun,pp);
  if(cmpis(pfqr,65535)>=0) err(talker,"module is too large in mynfshanks");
  lbaby=itos(pfqr);N=lgef((GEN)nf[1])-3;
  p1=g0;p2=cgetg(N+1,18);p2[1]=un;for(i=2;i<=N;i++) p2[i]=zero;
  smalltable=cgetg(lbaby+1,17);fl=0;
  for(i=1;(i<=lbaby)&&(!fl);i++)
  {
    smalltable[i]=(long)lift(element_mulmodpr(nf,p1,p2,prhall));p2=(GEN)smalltable[i];
    p3=lift(gmul(unmodp,gsub(p2,x)));fl=gcmp0(p3);
  }
  if(fl){avma=av;return stoi(i-1);}
  permtable=indexlexsort(smalltable);smalltable2=cgetg(lbaby+1,17);
  for(i=1;i<=lbaby;i++) smalltable2[i]=smalltable[itos((GEN)permtable[i])];
  giant=(GEN)smalltable[lbaby];
  smalltable=smalltable2;
  kgen=stoi(lbaby);p2=giant;
  while(1)
  {
    kgen=addis(kgen,lbaby);
    p2=lift(element_mulmodpr(nf,giant,p2,prhall));
    i=nfsearch(smalltable,p2);
    if(i)
    {tetpil=avma;return gerepile(av,tetpil,subii(kgen,(GEN)permtable[i]));}
  }
}

GEN
hermiteconstant(long n) /* rend constante d'Hermite^n si connue, sinon une borne sup */
{
  long av=avma,tetpil;
  GEN h;
  
  switch(n)
  {
    case 1: h=gun;break;
    case 2: h=gdivgs(stoi(4),3);break;
    case 3: h=gdeux;break;
    case 4: h=stoi(4);break;
    case 5: h=stoi(8);break;
    case 6: h=gdivgs(stoi(64),3);break;
    case 7: h=stoi(64);break;
    case 8: h=stoi(256);break;
    default: h=gmul(gpuigs(gdiv(gdeux,mppi(5)),n),gsqr(ggamma(gdivgs(stoi(n+4),2),5)));
  }
  tetpil=avma;return gerepile(av,tetpil,gcopy(h));
}

long
isprimitive(GEN nf) /* 1 primitif, 0 s'il est peut etre imprimitif... */
{
  long N,firstfactor,i,av,lfdK,maxexpo;
  GEN dK,fN,fN1,fdK,fdK2;

  av=avma;N=lgef((GEN)nf[1])-3;
  dK=(GEN)nf[3];
  fN=factor(stoi(N));fN1=(GEN)fN[1];
  firstfactor=itos((GEN)fN1[1]);
  if(firstfactor==N){avma=av;return 1;}
  fdK=factor(absi(dK));fdK2=(GEN)fdK[2];lfdK=lg(fdK2)-1;
  if(mpodd(dK)){maxexpo=1;i=1;}else{maxexpo=itos((GEN)fdK2[1])>>1;i=2;}
  for(;i<=lfdK;i++) maxexpo=max(maxexpo,itos((GEN)fdK2[i]));
  avma=av;return (long)(maxexpo<firstfactor);
}

GEN
badboundforregulator(GEN bnf,GEN units)
{
  long av,tetpil,N,R1,R2,R;
  GEN nf,dKa,defaultbound,p1,signature,c1;

  av=avma;
  nf=(GEN)bnf[7];
  N=lgef((GEN)nf[1])-3;
  defaultbound=dbltor(0.2);
  if(!isprimitive(nf))
  {
    if(DEBUGLEVEL){fprintferr("Default bound for regulator: 0.2\n");flusherr();}
    av=avma;return defaultbound;
  }
  else
  {
    dKa=absi((GEN)nf[3]);
    signature=(GEN)nf[2];R1=itos((GEN)signature[1]);R2=itos((GEN)signature[2]);
    R=R1+R2-1;
    if((!R2)&&(N<12)) c1=gpuigs(stoi(4),N>>1);else c1=gpuigs(stoi(N),N);
    if(cmpii(dKa,c1)<=0)
    {
      if(DEBUGLEVEL){fprintferr("Default bound for regulator: 0.2\n");flusherr();}
      av=avma;return defaultbound;
    }
    else
    {
      p1=gmax(defaultbound,gsqrt(gdiv(gdivgs(gmul2n(gpuigs(gdivgs(gmulgs(gsqr(glog(gdiv(dKa,c1),DEFAULTPREC)),3),N*(N*N-1)-6*R2),R),R2),N),hermiteconstant(R)),DEFAULTPREC));
      if(DEBUGLEVEL)
      {fprintferr("Mahler bound for regulator: ");outerr(p1);flusherr();}
      tetpil=avma;return gerepile(av,tetpil,gcopy(p1));
    }
  }
}

#define NBMAX 5000

GEN
minimforunits(GEN nf, long borne, long stockmax)
{
  GEN u,r,unr,S,S1,a,base,p1;
  long n1,n,av,av1,nonnul,i,j,k,s,norme,normax,flg,*x,fl1,cmpt;
  double p,b,c;
  double **q,*v,*y,*z;
  double eps=0.000001;

  av=avma;
  if(DEBUGLEVEL)
  {fprintferr("Searching minimum of T2-form on units: ");flusherr();}
  if(DEBUGLEVEL>2){fprintferr("\n   borne = %ld\n",borne);flusherr();}
  nf=checknf(nf);
  a=(GEN)((GEN)nf[5])[3];n1=lg(a);
  n=n1-1;
  x=(long*)malloc(n1*sizeof(long));
  y=(double*)malloc(n1*sizeof(double));
  z=(double*)malloc(n1*sizeof(double));
  v=(double*)malloc(n1*sizeof(double));
  q=(double**)malloc(n1*sizeof(double*));
  for(j=1;j<=n;j++) q[j]=(double*)malloc(n1*sizeof(double));
  affsr(1,unr=cgetr(10));
  u=lllgram(a,BIGDEFAULTPREC);base=gmul((GEN)nf[7],u);
  a=gmul(gtrans(u),gmul(a,u));a=gmul(a,unr);
  r=sqred1(a);
  for(j=1;j<=n;j++)
  {
    v[j]=rtodbl(gcoeff(r,j,j));
    for(i=1;i<j;i++)
      q[i][j]=rtodbl(gcoeff(r,i,j));
  }
  if(!(flg=borne))   /* flg=0 <==> chercher les vec min */
  {
    b=rtodbl(gcoeff(a,1,1));
    for(i=2;i<=n;i++)
      if((c=rtodbl(gcoeff(a,i,i)))<b) b=c;
    borne=(long)(b+eps);
  }
  normax=0;
  if(stockmax) S=cgetg(stockmax+1,19);
  s=0;k=n;cmpt=0;
  y[n]=z[n]=0;
  x[n]=(long)(sqrt(borne/v[n]+eps));
  do
  {
    do
    {
      if(k>1)
      {
	k--;
	z[k]=0;
	for(j=k+1;j<=n;j++) z[k]=z[k]+q[k][j]*x[j];
	p=x[k+1]+z[k+1];
	y[k]=y[k+1]+p*p*v[k+1];
	x[k]=(long)floor(sqrt((borne-y[k]+eps)/v[k])-z[k]);
      }
      while(v[k]*(x[k]+z[k])*(x[k]+z[k])>borne-y[k]+eps)
      {k++;x[k]--;}
    }
    while(k>1);
    if((nonnul=(x[1]||(y[1]>eps))))
    {
      cmpt++;
      if(cmpt>NBMAX)
      {
	free(x);free(y);free(z);free(v);for(j=1;j<=n;j++) free(q[j]);free(q);
	av=avma;return gzero;
      }
      if(DEBUGLEVEL>8){fprintferr(".");flusherr();}
      norme=(long)(y[1]+v[1]*(x[1]+z[1])*(x[1]+z[1])+eps);
      if(!flg&&(norme<borne)) {s=0; borne=norme;}
      if(norme>normax) normax=norme;
      av1=avma;p1=gzero;for(i=1;i<=n;i++) p1=gadd(p1,gmulsg(x[i],(GEN)base[i]));
      fl1=gcmp1(gabs(subres(p1,(GEN)nf[1]),0));avma=av1;
      if(fl1)
      {
	if(DEBUGLEVEL){fprintferr("*");flusherr();}
	s++;cmpt=0;
	if(s<=stockmax)
	{
	  p1=cgetg(n+1,18);
	  for(i=1;i<=n;i++) p1[i]=lstoi(x[i]);
	  S[s]=lmul(u,p1);
	}
      }
      x[k]--;
    }
  }
  while(nonnul);
  free(x);free(y);free(z);free(v);
  for(j=1;j<=n;j++) free(q[j]);free(q);
  if(DEBUGLEVEL){fprintferr("\n");flusherr();}
  if(stockmax) 
  {
    av1=avma;
    k=(s<stockmax)? s:stockmax;
    S1=cgetg(k+1,19);for(j=1;j<=k;j++) S1[j]=lcopy((GEN)S[j]);
    S=gerepile(av,av1,S1);
  }
  else {avma=av;S=cgetg(1,19);}
  u=cgetg(4,17);
  u[1]=lstoi(s<<1);
  u[2]=(flg)?lstoi(normax):lstoi(borne);
  u[3]=(long)S;
  return u;
}

#undef NBMAX

GEN
lowerboundforregulator(GEN bnf,GEN units)
{
  long av,N,R1,R2,R,tetpil,i,nbrootsofone,nbmin;
  GEN rootsofone,nf,signature,M0,M,m,col,T2,bound,minunit,newminunit,unitsonbasis;
  GEN vecminim,colalg,p1,pol,y;
    
  av=avma;
  rootsofone=(GEN)((GEN)bnf[8])[4];nbrootsofone=itos((GEN)rootsofone[1]);
  nf=(GEN)bnf[7];
  T2=(GEN)((GEN)nf[5])[3];
  N=lgef(nf[1])-3;
  signature=(GEN)nf[2];
  R1=itos((GEN)signature[1]);R2=itos((GEN)signature[2]);R=R1+R2-1;
  unitsonbasis=algtobasis(bnf,units);
  minunit=gmul(gmul(gtrans((GEN)unitsonbasis[1]),T2),(GEN)unitsonbasis[1]);
  for(i=2;i<=R;i++)
  {
    newminunit=gmul(gmul(gtrans((GEN)unitsonbasis[i]),T2),(GEN)unitsonbasis[i]);
    if(gcmp(newminunit,minunit)<0) minunit=newminunit;
  }
  if(gcmpgs(minunit,1000000000)>0){avma=av;return badboundforregulator(bnf,units);}
  vecminim=minimforunits(nf,itos(gceil(minunit)),10000);
  if(gcmp0(vecminim)){av=avma;return badboundforregulator(bnf,units);}
  m=(GEN)vecminim[3];nbmin=lg(m)-1;
  if(nbmin==10000){avma=av;return badboundforregulator(bnf,units);}
  bound=gaddgs(minunit,1);
  for(i=1;i<=nbmin;i++)
  {
    col=(GEN)m[i];colalg=basistoalg(nf,col);
    if(!gcmp1(lift(gpuigs(colalg,nbrootsofone))))
    {
      newminunit=gmul(gmul(gtrans(col),T2),col);
      if(gcmp(newminunit,bound)<0) bound=newminunit;
    }
  }
  if(gcmp(bound,minunit)>0) err(talker,"bug in lowerboundforregulator");
  if(DEBUGLEVEL){fprintferr("M* = ");outerr(bound);flusherr();}
  p1=polx[0];pol=gaddgs(gsub(gpuigs(p1,N),gmul(bound,p1)),N-1);
  if(DEBUGLEVEL>1){fprintferr("pol = ");outerr(pol);flusherr();}
  p1=roots(pol,6);
  if(N&1) y=greal((GEN)p1[3]);else y=greal((GEN)p1[2]);
  if(DEBUGLEVEL>1){fprintferr("y = ");outerr(y);flusherr();}
  M0=gmul2n(gmulsg(N*(N-1),gsqr(glog(y,6))),-2);
  if(DEBUGLEVEL){fprintferr("M0 = ");outerr(M0);flusherr();}
  M=gmul2n(gdivgs(gdiv(gpuigs(M0,R),hermiteconstant(R)),N),R2);
  if(DEBUGLEVEL){fprintferr("M = ");outerr(gsqrt(M,6));flusherr();}
  tetpil=avma;return gerepile(av,tetpil,gsqrt(M,6));
}

GEN
computeprimelement(GEN bnf,GEN pr,GEN prhall)
/* calcule un generateur du groupe cyclique (Zk/pr)* */
{
  long N,f,i,j,*vec,*vecindex,av,av1,tetpil,fl,r;
  GEN nf,pp,pfm1,ideal,bun,lambda,blambda,order;

  av=avma;
  nf=(GEN)bnf[7];N=lgef((GEN)nf[1])-3;pp=(GEN)pr[1];f=itos((GEN)pr[4]);
  pfm1=addis(gpuigs(pp,f),-1);
  ideal=idealmul(nf,pr,idmat(N));
  vec=(long*)malloc(sizeof(long)*(N+1));for(i=1;i<=N;i++) vec[i]=0;
  vecindex=(long*)malloc(sizeof(long)*(f+1));
  j=0;for(i=1;i<=N;i++) if(gegal(gcoeff(ideal,i,i),pp)) vecindex[++j]=i;
  bun=cgetg(N+1,18);bun[1]=un;for(i=2;i<=N;i++) bun[i]=zero;
  fl=r=1;
  for(i=1;i<=f;i++) vec[vecindex[i]]=r;
  av1=avma;
  do
  {
    avma=av1;
    lambda=cgetg(N+1,17);for(i=1;i<=N;i++) lambda[i]=(long)stoi(vec[i]);
    blambda=gtrans(lambda);
    order=mynfshanks(nf,bun,blambda,pr,prhall);
    if(DEBUGLEVEL>1)
    {
      fprintferr("               The element ");outerr(blambda);
      fprintferr("               has order ");outerr(order);
      flusherr();
    }
    if(gegal(order,pfm1)) fl=0;
    if(fl)
    {
      for(j=f;(vec[vecindex[j]]+r)==0;j--);
      vec[vecindex[j]]--;
      for(i=j+1;i<=f;i++) vec[vecindex[i]]=r;
      for(j=1;(j<f)&&(!vec[vecindex[j]]);j++);
      if(!vec[vecindex[j]])
      {
	r++;for(i=1;i<=f;i++) vec[vecindex[i]]=r;
	if(cmpis(pp,(r<<1))<0) err(talker,"bug in computeprimelement");
      }
    }
  }
  while(fl);
  free(vec);free(vecindex);tetpil=avma;
  return gerepile(av,tetpil,gcopy(blambda));
}
  
long
primecertify(GEN bnf,GEN beta,long pp,long nbalphapro,GEN big)
/* calcule une matrice carree de rang lg(beta) associee a une famille d'ideaux premiers
   P_i tels que : 1<=i<=nombre de beta ; N(P_i) congru a 1 mod pp
   v_(P_i)(beta[j])=0 pour tout 1<=j<=nbalphapro et 1<=i<=lg(beta) */

{
  long av,i,j,qq,nbcol,sizeofmat,flag,nbqq,ra,N;
  GEN nf,pgen,mat,mat1,qgen,decqq,newcol,eltgen,qrhall,ideal,eltgenpro;
  
  av=avma;
  pgen=stoi(pp);nbcol=0;nf=(GEN)bnf[7];N=lgef((GEN)nf[1])-3;
  sizeofmat=lg(beta)-1;mat=cgetg(1,19);
  qq=2*pp+1;
  qgen=stoi(qq);flag=1;
  while(flag)
  {
    if(isprime(qgen))
    {
      if(!gdivise(big,qgen))
      {
	decqq=primedec(bnf,qgen);nbqq=lg(decqq)-1;
	for(i=1;i<=nbqq;i++)
	{
	  ideal=(GEN)decqq[i];
	  if(itos((GEN)ideal[4])==1)
	  {
	    qrhall=compute_prhall(nf,ideal);
	    nbcol++;
	    newcol=cgetg(sizeofmat+1,18);
	    if(DEBUGLEVEL>1)
	    {
	      fprintferr("          The prime number q : %ld\n",qq);
	      fprintferr("          The prime ideal Q : ");outerr(ideal);
	      flusherr();
	    }
	    eltgenpro=lift(gener(qgen));
	    eltgen=cgetg(N+1,18);eltgen[1]=(long)eltgenpro;
	    for(j=2;j<=N;j++) eltgen[j]=zero;
	    for(j=1;j<=sizeofmat;j++)
	    {
	      newcol[j]=(long)nfshanks(nf,(GEN)beta[j],eltgen,ideal,qrhall);
	    }
	    if(DEBUGLEVEL>1)
	    {
	      fprintferr("          The generator of (Zk/Q)^* : ");outerr(eltgen);
	      fprintferr("          The column #%ld of the matrix log(b_j/Q) : ",nbcol);
	      outerr(newcol);
	      flusherr();
	    }
	    mat1=concat(mat,newcol);
	    ra=rank(mat1);
	    if(DEBUGLEVEL>1)
	    {
	      fprintferr("          The new rank of the matrix : %ld\n\n",ra);
	      flusherr();
	    }
	    if(ra!=nbcol) nbcol--;
	    else
	    {
	      mat=mat1;
	      if(nbcol==sizeofmat){i=nbqq+1;flag=0;}
	    }
	  }
	}
      }
    }
    qq+=2*pp;qgen=stoi(qq);
  }
  avma=av;return 1;
}

long
certifybuchall(GEN bnf)
{
  byteptr delta;
  long nbgenerator,i,j,pp,N,R1,R2,R,av,av1,nfa,nbf1,flag,nbalpha,nbalphapro,fl,testprime;
  GEN nf,cl,reg,rootsofone,fundunits,generator,p1,boundforprime;
  GEN zimmertbound,DK,alpha,factfd1,f1,classnumber;
  GEN cycliccomponent,nbrootsofone,pgen,beta,big;

  av=avma;
  nf=(GEN)bnf[7];N=lgef((GEN)nf[1])-3;
  R1=itos((GEN)((GEN)nf[2])[1]);R2=itos((GEN)((GEN)nf[2])[2]);R=R1+R2-1;
  DK=absi((GEN)nf[3]);
  zimmertbound=computezimmertbound(N,R1,DK);
  if(DEBUGLEVEL){fprintferr("Zimmert's bound = ");outerr(zimmertbound);flusherr();}
  testprimeuptominkowskibound(bnf,zimmertbound);
  cl=(GEN)((GEN)bnf[8])[1];classnumber=(GEN)cl[1];
  cycliccomponent=(GEN)cl[2];nbgenerator=lg(cycliccomponent)-1;generator=(GEN)cl[3];
  reg=(GEN)((GEN)bnf[8])[2];
  rootsofone=(GEN)((GEN)bnf[8])[4];nbrootsofone=(GEN)rootsofone[1];
  fundunits=(GEN)buchfu(bnf)[1];
  if(DEBUGLEVEL)
  {
    fprintferr("Class number = ");outerr(classnumber);
    fprintferr("Cyclic components = ");outerr(cycliccomponent);
    fprintferr("Generators = ");outerr(generator);
    fprintferr("Regulator = ");outerr(reg);
    fprintferr("Roots of one = ");outerr(rootsofone);
    fprintferr("Fundamental units = ");outerr(fundunits);
    flusherr();
  }
  alpha=cgetg(nbgenerator+1,17);
  for(i=1;i<=nbgenerator;i++)
  {
    p1=idealpow(nf,(GEN)generator[i],(GEN)cycliccomponent[i]);
    alpha[i]=isprincipalgen(bnf,p1)[2];
  }
  boundforprime=ground(gdiv(reg,lowerboundforregulator(bnf,fundunits)));
  if(DEBUGLEVEL)
  {
    fprintferr("\nPHASE 2 : are all primes good ?\n\n");
    flusherr();
  }
  if(DEBUGLEVEL)
  {
    fprintferr("  Bound for testing primes : B = ");outerr(boundforprime);flusherr();
  }
  nfa=0;
  if(nbgenerator)
  {
    factfd1=factor((GEN)cycliccomponent[1]);
    nbf1=lg((GEN)factfd1[1])-1;f1=(GEN)factfd1[1];
    nfa=0;for(i=1;i<=nbf1;i++){if(cmpii((GEN)f1[i],boundforprime)>0) nfa++;}
  }
  if(DEBUGLEVEL)
  {
    fprintferr("  Number of primes > boundforprime = %ld\n",nfa);flusherr();
  }
  delta=diffptr;pp=*delta++;flag=1;pgen=stoi(pp);
  if(DEBUGLEVEL)
  {
    fprintferr("\n  Testing primes <= B\n\n");flusherr();
  }
  while(cmpsi(pp,boundforprime)<1)
  {
    if(DEBUGLEVEL)
    {
      fprintferr("\n***** Testing prime p = %ld\n",pp);flusherr();
    }
    av1=avma;
    nbalphapro=0;
    if(gdivise(classnumber,pgen))
    {
      if(DEBUGLEVEL>1)
      {
	fprintferr("      p divides cl(k)\n");flusherr();
      }
      nbalphapro=nbgenerator;fl=1;
      while(fl&&(nbalphapro>0))
      {if(!gdivise((GEN)cycliccomponent[nbalphapro],pgen)) nbalphapro--;else fl=0;}
    }
    if(gdivise(nbrootsofone,pgen))
    {
      if(DEBUGLEVEL>1)
      {
	fprintferr("      p divides w(k)\n");flusherr();
      }
      fl=1;nbalpha=nbalphapro+1+R;
    }
    else{fl=0;nbalpha=nbalphapro+R;}
    if(DEBUGLEVEL>1)
    {
      fprintferr("      t+r+e = %ld\n",nbalpha);flusherr();
    }
    beta=cgetg(nbalpha+1,17);
    for(i=1;i<=nbalphapro;i++) beta[i]=alpha[i];
    if(fl) beta[nbalphapro+1]=(long)algtobasis(nf,(GEN)rootsofone[2]);
    for(i=1;i<=R;i++) beta[i+nbalphapro+fl]=(long)algtobasis(nf,(GEN)fundunits[i]);
    if(DEBUGLEVEL>1)
    {
      fprintferr("      Beta list = ");outerr(beta);fprintferr("\n");flusherr();
    }
    big=gun;for(i=1;i<=nbgenerator;i++) big=gmul(big,idealnorm(nf,(GEN)generator[i]));
    testprime=primecertify(bnf,beta,pp,nbalphapro,big);
    avma=av1;pp+=*delta++;pgen=stoi(pp);
  }
  for(j=1;j<=nfa;j++)
  {
    if(DEBUGLEVEL)
    {
      fprintferr("  Testing primes > B\n\n");flusherr();
    }
    av1=avma;
    pgen=(GEN)f1[nbf1-j+1];pp=itos(pgen);
    if(DEBUGLEVEL)
    {
      fprintferr("***** Testing prime p = %ld\n",itos(pgen));flusherr();
    }
    nbalphapro=nbgenerator;fl=1;
    while(fl&&(nbalphapro>0))
    {if(!gdivise((GEN)cycliccomponent[nbalphapro],pgen)) nbalphapro--;else fl=0;}
    if(gdivise(nbrootsofone,pgen)){fl=1;nbalpha=nbalphapro+1+R;}
    else{fl=0;nbalpha=nbalphapro+R;}
    if(DEBUGLEVEL>1)
    {
      fprintferr("      t+r+e = %ld\n",nbalpha);flusherr();
    }
    beta=cgetg(nbalpha+1,17);
    for(i=1;i<=nbalphapro;i++) beta[i]=alpha[i];
    if(fl) beta[nbalphapro+1]=(long)algtobasis(nf,(GEN)rootsofone[2]);
    for(i=1;i<=R;i++) beta[i+nbalphapro+fl]=(long)algtobasis(nf,(GEN)fundunits[i]);
    if(DEBUGLEVEL>1)
    {
      fprintferr("      Liste des beta = ");outerr(beta);flusherr();
    }
    big=gun;for(i=1;i<=nbgenerator;i++) big=gmul(big,idealnorm(nf,(GEN)generator[i]));
    testprime=primecertify(bnf,beta,pp,nbalphapro,big);
    avma=av1;
  }
  avma=av;return 1;
}

