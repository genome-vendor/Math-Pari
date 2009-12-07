#include <stdio.h>
#include <stdlib.h>
#include <setjmp.h>

#ifdef __cplusplus
  extern "C" {
#endif

/* CAT2:
 *      This macro catenates 2 tokens together.
 */
/* STRINGIFY:
 *      This macro surrounds its token with double quotes.
 */
#ifndef CAT2
# if 42 == 1
#  define CAT2(a,b)a/**/b
#  define CAT3(a,b,c)a/**/b/**/c
#  define CAT4(a,b,c,d)a/**/b/**/c/**/d
#  define CAT5(a,b,c,d,e)a/**/b/**/c/**/d/**/e
#  define STRINGIFY(a)"a"
                /* If you can get stringification with catify, tell me how! */
# endif
# if 42 == 42
#  define CAT2(a,b)a ## b
#  define CAT3(a,b,c)a ## b ## c
#  define CAT4(a,b,c,d)a ## b ## c ## d
#  define CAT5(a,b,c,d,e)a ## b ## c ## d ## e
#  define StGiFy(a)# a
#  define STRINGIFY(a)StGiFy(a)
#  define SCAT2(a,b)StGiFy(a) StGiFy(b)
#  define SCAT3(a,b,c)StGiFy(a) StGiFy(b) StGiFy(c)
#  define SCAT4(a,b,c,d)StGiFy(a) StGiFy(b) StGiFy(c) StGiFy(d)
#  define SCAT5(a,b,c,d,e)StGiFy(a) StGiFy(b) StGiFy(c) StGiFy(d) StGiFy(e)
# endif
# ifndef CAT2
#   include "Bletch: How does this C preprocessor catenate tokens?"
# endif
#endif /* CAT2 */



#ifndef USE_JUNK

char *
alloc(unsigned long size, char *name)
{
  return (char *)malloc((size_t)size);
}

void *
const_express() {return NULL;}

extern int term;
int term;
float xsize=1.0, ysize=1.0, pointsize=1.0;		/* During test! */

#ifndef BITS_IN_HALFULONG /* In pari it is already defined. */
  FILE *outfile = stdout;
#endif

char term_options[4] = "";
#define MAX_ID_LEN 50
extern char     default_font[]; 
char            default_font[MAX_ID_LEN+1] = "\0"; /* Entry added by DJL */
extern char outstr[];
char        outstr[MAX_ID_LEN+1] = "STDOUT";
extern double ticscale; /* scale factor for tic marks (was (0..1])*/
double        ticscale = 1.0; /* scale factor for tic mark */

jmp_buf env;

char *input_line;
int inline_num;          /* from command.c */
int interactive;    /* from plot.c */
char *infile_name;       /* from plot.c */

/* Not used: */

char *token;
long c_token, num_tokens;

#endif /* !defined(USE_JUNK) */


/* Cannot pull the whole plot.h, too many contradictions. */

#ifdef __ZTC__
typedef int (*FUNC_PTR)(...);
#else
typedef int (*FUNC_PTR)();
#endif

struct TERMENTRY {
        char *name;
#if defined(_Windows) && !defined(WIN32)
        char GPFAR description[80];     /* to make text go in FAR segment */
#else
        char *description;
#endif
        unsigned int xmax,ymax,v_char,h_char,v_tic,h_tic;
        FUNC_PTR options,init,reset,text,scale,graphics,move,vector,linetype,
                put_text,text_angle,justify_text,point,arrow;
};

#ifdef _Windows
#  define termentry TERMENTRY far
#else
#  define termentry TERMENTRY
#endif

extern struct termentry term_tbl[];
 
#define RETVOID 
#define RETINT , 1

#define F_0 void(*)()
#define F_1 void(*)(int)
#define F_1I int(*)(int)
#define F_2 void(*)(unsigned int,unsigned int)
#define F_2D int(*)(double,double)
#define F_3 void(*)(unsigned int,unsigned int,int)
#define F_3T void(*)(int,int,char*)
#define F_4 void(*)(int,int,int,int)
#define F_5 void(*)(int,int,int,int,int)

#define CALL_G_METH0(method) CALL_G_METH(method,0,(),RETVOID)
#define CALL_G_METH1(method,arg1) CALL_G_METH(method,1,(arg1),RETVOID)
#define CALL_G_METH1I(method,arg1) CALL_G_METH(method,1I,(arg1),RETINT)
#define CALL_G_METH2(method,arg1,arg2) \
		CALL_G_METH(method,2,((arg1),(arg2)),RETVOID)
#define CALL_G_METH2D(method,arg1,arg2) \
		CALL_G_METH(method,2D,((arg1),(arg2)),RETINT)
#define CALL_G_METH3(method,arg1,arg2,arg3) \
		CALL_G_METH(method,3,((arg1),(arg2),(arg3)),RETVOID)
#define CALL_G_METH3T(method,arg1,arg2,arg3) \
		CALL_G_METH(method,3T,((arg1),(arg2),(arg3)),RETVOID)
#define CALL_G_METH4(method,arg1,arg2,arg3,arg4) \
		CALL_G_METH(method,4,((arg1),(arg2),(arg3),(arg4)),RETVOID)
#define CALL_G_METH5(method,arg1,arg2,arg3,arg4,arg5) \
		CALL_G_METH(method,5,((arg1),(arg2),(arg3),(arg4),(arg5)),RETVOID)

#define CALL_G_METH(method,mult,args,returnval)    (		\
       (term<=0) ? (						\
	 croak("No terminal specified") returnval		\
       ) :							\
       (*(CAT2(F_,mult))my_term_tbl[term].method)args		\
     )

#define init()	CALL_G_METH0(init)
#define reset()	CALL_G_METH0(reset)
#define text()	CALL_G_METH0(text)
#define graphics()	CALL_G_METH0(graphics)
#define linetype(lt)	CALL_G_METH1(linetype,lt)
#define justify_text(mode)	CALL_G_METH1I(justify_text,mode)
#define text_angle(ang)	CALL_G_METH1I(text_angle,ang)
#define scale(xs,ys)	CALL_G_METH2D(scale,xs,ys)
#define move(x,y)	CALL_G_METH2(move,x,y)
#define vector(x,y)	CALL_G_METH2(vector,x,y)
#define put_text(x,y,str)	CALL_G_METH3T(put_text,x,y,str)
#define point(x,y,p)	CALL_G_METH3(point,x,y,p)
#define arrow(sx,sy,ex,ey,head)	CALL_G_METH5(arrow,sx,sy,ex,ey,head)

#define termprop(prop) (my_term_tbl[term].prop)
#define termset(term) my_change_term(term,strlen(term))
int change_term(char*,int);

#ifdef DYNAMIC_PLOTTING			/* Can load plotting DLL later */

UNKNOWN_null()
{
    err(talker,"gnuplot-like plotting environment not loaded yet");
}

static FUNC_PTR change_term_p;

int 
my_change_term(char*s,int l) 
{
    if (!change_term_p)
	UNKNOWN_null();
    term = (*change_term_p)(s,l);
}

#  define change_term(p,l) my_change_term(p,l)
#  define term_tbl (my_term_tbl)

static struct termentry dummy_term_tbl[] = {
    {"unknown", "Unknown terminal type - not a plotting device",
	  100, 100, 1, 1,
	  1, 1, UNKNOWN_null, UNKNOWN_null, UNKNOWN_null, 
	  UNKNOWN_null, UNKNOWN_null, UNKNOWN_null, UNKNOWN_null, UNKNOWN_null, 
	  UNKNOWN_null, UNKNOWN_null, UNKNOWN_null,
     UNKNOWN_null, UNKNOWN_null, UNKNOWN_null},
};
static struct termentry *my_term_tbl = dummy_term_tbl;

/* This function should be called before any graphic code can be used... */
set_term_funcp(FUNC_PTR change_p, struct termentry *term_p)
{
    my_term_tbl = term_p;
    change_term_p = change_p;
}

#else /* !DYNAMIC_PLOTTING */

#  define my_change_term change_term
#  define my_term_tbl term_tbl

#endif /* DYNAMIC_PLOTTING */



#ifdef __cplusplus
  }
#endif
