/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                     PLOT EN HAUTE RESOLUTION                    */
/*                                                                 */
/*                       copyright Babe Cool                       */
/*                                                                 */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

# include "genpari.h"
#include "rect.h"

#ifdef HPPA
#ifndef __GNUC__
typedef char *caddr_t;
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <sys/types.h>
#include <unistd.h>
#ifdef __cplusplus
}
#endif


Colormap PARI_Colormap;		/* The Colormap We Are Using */
XColor  *PARI_Colors;
XColor  *PARI_ExactColors;

static char *PARI_DefaultColors[8] = {
  " ",
  "sienna", /* Courbe impaires dans plothmult */
  "cornsilk",   
  "red",        /* Courbes, ou courbes paires dans plothmult */
  "black",                    
  "grey",                     
  "blue",         /* axes */
  "gainsboro",              
};

void
PARI_ColorSetUp(Display *display, char **Colors, int n)
{
  static int NotInitialised = 1;
  if (NotInitialised)
  {
    int i;
    NotInitialised = 0;
    PARI_Colormap = DefaultColormap(display, 0);
    PARI_Colors = (XColor *) malloc((n + 1) * sizeof(XColor));
    PARI_ExactColors = (XColor *) malloc((n + 1) * sizeof(XColor));
    for (i = 1; i < n; i++)
    {
      XAllocNamedColor(display, PARI_Colormap, Colors[i],
		       &PARI_ExactColors[i], &PARI_Colors[i]);
    }
  }
}

#ifdef SONY

typedef int (*XErrorHandler) (
#if NeedFunctionPrototypes
  Display*	
  XErrorEvent*
#endif
  );

extern XErrorHandler XSetErrorHandler (
#if NeedFunctionPrototypes
  XErrorHandler	
#endif
  );


typedef int (*XIOErrorHandler) ( 
#if NeedFunctionPrototypes
  Display*	
#endif
  );

extern XIOErrorHandler XSetIOErrorHandler (
#if NeedFunctionPrototypes
  XIOErrorHandler
#endif
  );
#endif

int
xerror()
{
  err(ploter1);return 0;
}

GEN
rectdraw(GEN list)
{
  long *ptx,*pty,*numpoints,*numtexts,*xtexts,*ytexts;
  long n,i,j,x0,y0,av=avma;
  long a,b,ro_cnt[ROt_MAX+1],ne;
  char **texts;
  Rect *e;
  RectObj *p1;

  Display *display;
  Window win;
  XSizeHints size_hints;
  XEvent report;
  GC gc;
  XFontStruct *font_info;
  int screen;
  XPoint *points,**lines;
  XSegment *segments;
  XRectangle *rectangles;

  if(fork()) return gnil;

  if(typ(list)!=17) err(rploter3);
  n=lg(list)-1;if(n%3) err(rploter4);
  n=n/3;if(!n) {abort();return gnil;}

  PARI_get_plot();

  if (!(display = XOpenDisplay(NULL))) err(ploter2);
  if (!(font_info = XLoadQueryFont(display, "9x15"))) err(ploter3);
  XSetErrorHandler((XErrorHandler)xerror);
  XSetIOErrorHandler((XIOErrorHandler)xerror);
  PARI_ColorSetUp(display,PARI_DefaultColors, 8);

  screen = DefaultScreen(display);
  ro_cnt[ROt_MV]=ro_cnt[ROt_PT]=ro_cnt[ROt_LN]=ro_cnt[ROt_BX]
    =ro_cnt[ROt_MP]=ro_cnt[ROt_ML]=ro_cnt[ROt_ST]=0;
  for(i=0;i<n;i++)
  {
    if(typ((GEN)list[3*i+1])!=1) err(rploter5);
    ne=itos((GEN)list[3*i+1]);if(!GOODRECT(ne)) err(rploter2);
    e=rectgraph[ne];
    p1=RHead(e);
    while(p1) 
    {
      if(RoType(p1)!=ROt_MP) ro_cnt[RoType(p1)]++;
      else ro_cnt[ROt_PT] += RoMPcnt(p1);
      p1=RoNext(p1);
    }
  }
  points=(XPoint*)malloc(ro_cnt[ROt_PT]*sizeof(XPoint));
  segments=(XSegment*)malloc(ro_cnt[ROt_LN]*sizeof(XSegment));
  rectangles=(XRectangle*)malloc(ro_cnt[ROt_BX]*sizeof(XRectangle));
  lines=(XPoint**)malloc(ro_cnt[ROt_ML]*sizeof(XPoint*));
  numpoints=(long*)malloc(ro_cnt[ROt_ML]*sizeof(long));
  texts=(char**)malloc(ro_cnt[ROt_ST]*sizeof(char*));
  numtexts=(long*)malloc(ro_cnt[ROt_ST]*sizeof(long));
  xtexts=(long*)malloc(ro_cnt[ROt_ST]*sizeof(long));
  ytexts=(long*)malloc(ro_cnt[ROt_ST]*sizeof(long));
  ro_cnt[ROt_PT]=ro_cnt[ROt_LN]=ro_cnt[ROt_BX]=ro_cnt[ROt_ML]=ro_cnt[ROt_ST]=0;
  for(i=0;i<n;i++)
  {
    e=rectgraph[itos((GEN)list[3*i+1])];x0=list[3*i+2];y0=list[3*i+3];
    if((typ((GEN)x0)!=1)||(typ((GEN)y0)!=1)) err(rploter5);
    x0=itos((GEN)x0);y0=itos((GEN)y0);
    p1=RHead(e);
    while(p1)
    {
      switch(RoType(p1))
      {
	case ROt_PT: 
	  points[ro_cnt[ROt_PT]].x=RoPTx(p1)+x0;
	  points[ro_cnt[ROt_PT]].y=RoPTy(p1)+y0;
	  ro_cnt[ROt_PT]++;break;
	case ROt_LN:
	  segments[ro_cnt[ROt_LN]].x1=RoLNx1(p1)+x0;
	  segments[ro_cnt[ROt_LN]].y1=RoLNy1(p1)+y0;
	  segments[ro_cnt[ROt_LN]].x2=RoLNx2(p1)+x0;
	  segments[ro_cnt[ROt_LN]].y2=RoLNy2(p1)+y0;
	  ro_cnt[ROt_LN]++;break;
	case ROt_BX:
	  a=rectangles[ro_cnt[ROt_BX]].x		= RoBXx1(p1)+x0;
	  b=rectangles[ro_cnt[ROt_BX]].y		= RoBXy1(p1)+y0;
	  rectangles[ro_cnt[ROt_BX]].width	= RoBXx2(p1)+x0-a;
	  rectangles[ro_cnt[ROt_BX]].height	= RoBXy2(p1)+y0-b;
	  ro_cnt[ROt_BX]++;break;
	case ROt_MP:
	  ptx=RoMPxs(p1);
	  pty=RoMPys(p1);
	  for(j=0;j<RoMPcnt(p1);j++)
	  {
	    points[ro_cnt[ROt_PT]+j].x=ptx[j]+x0;
	    points[ro_cnt[ROt_PT]+j].y=pty[j]+y0;
	  }
	  ro_cnt[ROt_PT]+=RoMPcnt(p1);break;
	case ROt_ML:
	  ptx=RoMLxs(p1);
	  pty=RoMLys(p1);
	  numpoints[ro_cnt[ROt_ML]]=RoMLcnt(p1);
	  lines[ro_cnt[ROt_ML]]=(XPoint*)malloc(RoMLcnt(p1)*sizeof(XPoint));
	  for(j=0;j<RoMLcnt(p1);j++)
	  {
	    lines[ro_cnt[ROt_ML]][j].x=ptx[j]+x0;
	    lines[ro_cnt[ROt_ML]][j].y=pty[j]+y0;
	  }
	  ro_cnt[ROt_ML]++;break;
	case ROt_ST: 
	  texts[ro_cnt[ROt_ST]]=RoSTs(p1); numtexts[ro_cnt[ROt_ST]]=RoSTl(p1);
	  xtexts[ro_cnt[ROt_ST]]=RoSTx(p1)+x0;ytexts[ro_cnt[ROt_ST]]=RoSTy(p1)+y0;
	  ro_cnt[ROt_ST]++;break;
	default: break;
      }
      p1=RoNext(p1);
    }
  }
  win= XCreateSimpleWindow(display, RootWindow(display, screen), 0, 0, w_width, w_height, 4, BlackPixel(display, screen), WhitePixel(display, screen));
  size_hints.flags = PPosition | PSize;
  size_hints.x = 0;
  size_hints.y = 0;
  size_hints.width = w_width;
  size_hints.height = w_height;
  
  XSetStandardProperties(display, win, "rectplot", NULL, None, NULL, 0, &size_hints);
  XSelectInput(display, win, ExposureMask | ButtonPressMask);
  gc = XCreateGC(display, win, 0, NULL);
  XSetFont(display, gc, font_info->fid);
  XMapWindow(display, win);
  for(;;)
  {
    XSetForeground(display, gc, PARI_Colors[3].pixel);
    XNextEvent(display, &report);
    if (report.type != Expose) break;
    while (XCheckTypedEvent(display, Expose, &report));
    if(ro_cnt[ROt_PT]) XDrawPoints(display, win, gc, points, ro_cnt[ROt_PT], 0);
    if(ro_cnt[ROt_LN]) XDrawSegments(display, win, gc, segments, ro_cnt[ROt_LN]);
    if(ro_cnt[ROt_BX]) XDrawRectangles(display, win, gc, rectangles, ro_cnt[ROt_BX]);
    for(i=0;i<ro_cnt[ROt_ML];i++) XDrawLines(display,win,gc,lines[i],numpoints[i],0);
    for(i=0;i<ro_cnt[ROt_ST];i++) XDrawString(display,win,gc,xtexts[i],ytexts[i],texts[i],numtexts[i]);
  }
  XUnloadFont(display, font_info->fid);
  XFreeGC(display, gc);
  free(points);free(segments);free(rectangles);
  free(numpoints);for(i=0;i<ro_cnt[ROt_ML];i++) free(lines[i]);
  free(lines);free(texts);free(numtexts);free(xtexts);free(ytexts);
  XCloseDisplay(display);
  avma = av;abort();return gnil;
}

void
PARI_get_plot()
{
  Display *display;
  int screen;

  if (pari_plot.init) {
    return;
  }
  if (!(display = XOpenDisplay(NULL))) err(ploter2);

  screen = DefaultScreen(display);
  w_width = DisplayWidth(display, screen) - 40;
  w_height = DisplayHeight(display, screen) - 60; 
  f_height = 15;
  f_width = 9;
  h_unit = 5;
  v_unit = 5;
  pari_plot.init = 1;
  XCloseDisplay(display);
}

long
term_set(char *s)
{
  return -1;
}
