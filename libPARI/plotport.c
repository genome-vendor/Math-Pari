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


PARI_plot pari_plot;

#define DRAWRECT (NUMRECT-2)
#define GRAPHRECT (NUMRECT-1)

double
tickmarks(double a)
{
  long k;
  double m,e;
  
  if(a<=0) return 0.0;
  k=(long)(log(a)/log(10));e=pow(10.0,(double)k);
  m=a/e;
  if(m<=2) return 0.1*e;
  if(m<=5) return 0.2*e;
  return 0.5*e;
}

void
rect_string(long win, long x, long y, char *str)
{
  long i, l = strlen(str);
  GEN vect = cgetg(l+1, 17), xx = cgeti(3), yy = cgeti(3);

  for (i = 0; i < l; i++) {
    vect[i+1] = lgeti(3);
    affsi(str[i], (GEN)vect[i+1]);
  }
  affsi(x, xx);
  affsi(y, yy);
  rectmove(win, xx, yy);
  rectstring(win, vect);
}

#define NBPOINTS 1500

/* 
 * Actually it uses two drawing rectangles: one for strings, another
 * for graph.
 */

/*
 * data is an array of arrays. Its meaning depends on flags.
 * 
 * If flags contain PLOT_PARAMETRIC, the length of array should be be
 * even, and pairs represent curves to plot.
 *
 * If there is no such flag, the first element is an array with
 * x-coordinates, rest with y-coordinates of lines to draw.
 *
 * Additional flags: PLOT_NO_AXE_X, PLOT_NO_AXE_Y, PLOT_NO_FRAME.
 */

void
rectplothrawin(long drawrect, GEN xmin, GEN xmax, GEN ymin, GEN ymax, 
	    GEN data, long prec, uLong flags)
{
  long av = avma,i,j,nc,l, ltype = 0;
  /* double sht,svt; */
  
  rectscale(drawrect, xmin, xmax, ymin, ymax);

  /* sht=tickmarks(rtodbl(xdiff));svt=tickmarks(rtodbl(ydiff)); */

  if (!(flags & PLOT_NO_FRAME)) {
    rectlinetype(drawrect, -2);		/* Frame. */
    rectmove(drawrect, xmin, ymin);
    rectbox(drawrect, xmax, ymax);
  }

  if (!(flags & PLOT_NO_AXE_Y) && gsigne(xmin)*gsigne(xmax)<0) {
    rectlinetype(drawrect, -1);		/* Axes. */
    rectmove(drawrect, gzero, ymin);
    rectline(drawrect, gzero, ymax);
  }
  if(!(flags & PLOT_NO_AXE_X) && gsigne(ymin)*gsigne(ymax)<0) {
    rectlinetype(drawrect, -1);		/* Axes. */
    rectmove(drawrect, xmin, gzero);
    rectline(drawrect, xmax, gzero);
  }

  nc = lg(data);
  for(j = (flags & PLOT_PARAMETRIC ? 1 : 2); j < nc; j++)
    {
      rectpointtype(drawrect, ltype);		/* Graphs. */
      rectlinetype(drawrect, ltype++);		/* Graphs. */
      if (flags & PLOT_PARAMETRIC) {
	if (flags & PLOT_POINTS) {
	  rectpoints(drawrect, (GEN)data[j], (GEN)data[j+1]);
	} else {
	  rectlines(drawrect, (GEN)data[j], (GEN)data[j+1]);
	}
	j++;
      } else {
	if (flags & PLOT_POINTS) {
	  rectpoints(drawrect, (GEN)data[1], (GEN)data[j]);
	} else {
	  rectlines(drawrect, (GEN)data[1], (GEN)data[j]);
	}
      }
    }

  avma = av;
}

GEN
recplothraw(long drawrect, GEN frame, GEN data, long flags, long prec)
{
  if(typ(frame)>18 || typ(frame) <17) err(talker,"not a vector in rploth");
  rectplothrawin(drawrect, (GEN)frame[1], (GEN)frame[2], (GEN)frame[3],
	 (GEN) frame[4], data, prec, flags);
  return gnil;
}


/*
 *  Pure graphing.
 *  If testpoints vanishes, it is set to the width of rectangle.
 *  Returns the vector of the bounding box.
 */

GEN
recplothmultin(long graphrect, entree *ep, GEN a, GEN b, char *ch, long prec, uLong flags, long testpoints)
{
  long av = avma,av2,i,j,sig,nc, offset = 0, points;
  GEN p1,p2,ysml,ybig,xsml,xbig,x,ydiff,xdiff,dx,y,dr,res;

  ysml=cgetr(3);ybig=cgetr(3);
  xsml=cgetr(3);xbig=cgetr(3);
  res = cgetg(5,17);
  res[1] = (long)xsml;
  res[2] = (long)xbig;
  res[3] = (long)ysml;
  res[4] = (long)ybig;

  if (!testpoints) {
    points = 2 * (rectgraph[graphrect]->sizex +
		  rectgraph[graphrect]->sizey + 2);
    testpoints = min(points, NBPOINTS);
  }

  sig=gcmp(b,a); if(!sig) return gnil;
  if(sig<0) {x=a;a=b;b=x;}

  newvalue(ep,cgetr(prec)); x=(GEN)(ep->value);
  gaffect(a,x);p1=lisexpr(ch);
  if (flags & PLOT_PARAMETRIC) {
    dx=gdivgs(gsub(b,a), testpoints - 1);
  } else {
    offset = 1;
    gaffect(a,xsml);gaffect(b,xbig);
    dx=gdivgs(gsub(b,a), testpoints - 1);
  }
  if (flags & PLOT_SINGLE) {
    if(typ(p1)>=17) err(talker,"vector or matrix as input in ploth");
    nc=2;				/* Length+1. */
    gaffect(p1,ysml);gaffect(p1,ybig);
  } else {
    if(typ(p1)!=17) err(talker,"not a row vector in plothmult or ploth2");
    nc=lg(p1);if(nc<=1) {avma=av;return gnil;}
    if (flags & PLOT_PARAMETRIC) {
      gaffect((GEN)p1[1],xsml); gaffect(xsml,xbig);
      gaffect((GEN)p1[2],ysml); gaffect(ysml,ybig);
      for (i = 2; 2*i < nc; i++) {
	if (gcmp(xsml,(GEN)p1[2*i - 1]) > 0) {
	  gaffect((GEN)p1[2*i - 1],xsml);
	}
	if (gcmp(xbig,(GEN)p1[2*i - 1]) < 0) {
	  gaffect((GEN)p1[2*i - 1],xbig);
	}
	if (gcmp(ysml,(GEN)p1[2*i]) > 0) {
	  gaffect((GEN)p1[2*i],ysml);
	}
	if (gcmp(ybig,(GEN)p1[2*i]) < 0) {
	  gaffect((GEN)p1[2*i],ybig);
	}
      }
    } else {
      gaffect(vecmin(p1),ysml); gaffect(vecmax(p1),ybig);
    }
  }

  y = cgetg(nc + offset,17);

  for(i = 1; i < nc + offset; i++)
  {
    p1 = cgetg(testpoints+1,17);
    y[i] = (long)p1;
    for(j=1;j<testpoints+1;j++) p1[j] = lgetr(3);
  }
  av2=avma;
  for(i=1;i<=testpoints;i++)
  {
    p1 = lisexpr(ch);
    if (flags & PLOT_SINGLE) {
      gaffect(p1, (GEN)((GEN)y[2])[i]);
    } else {
      for (j=1; j<nc; j++) {
	gaffect((GEN)p1[j], (GEN)((GEN)y[j + offset])[i]);
      }
    }
    if (!(flags & PLOT_PARAMETRIC)) {
      gaffect(x, (GEN)((GEN)y[1])[i]);
    }
    if (flags & PLOT_PARAMETRIC) {
      for (j = 1; 2*j < nc; j++) {
	if (gcmp(xsml,(GEN)p1[2*j - 1]) > 0) {
	  gaffect((GEN)p1[2*j - 1],xsml);
	}
	if (gcmp(xbig,(GEN)p1[2*j - 1]) < 0) {
	  gaffect((GEN)p1[2*j - 1],xbig);
	}
	if (gcmp(ysml,(GEN)p1[2*j]) > 0) {
	  gaffect((GEN)p1[2*j],ysml);
	}
	if (gcmp(ybig,(GEN)p1[2*j]) < 0) {
	  gaffect((GEN)p1[2*j],ybig);
	}
      }
    } else if (flags & PLOT_SINGLE) {
      if(gcmp(p1, ysml) < 0) gaffect(p1,ysml);
      else if(gcmp(p1, ybig) > 0) gaffect(p1,ybig);
    } else {
      if(gcmp(p2=vecmin(p1), ysml) < 0) gaffect(p2,ysml);
      if(gcmp(p2=vecmax(p1), ybig) > 0) gaffect(p2,ybig);
    }
    gaddz(x,dx,x);avma=av2;
  }

  rectplothrawin(graphrect, xsml, xbig, ysml, ybig, y, prec, flags);

  killvalue(ep);
  avma = (long)res;

  return res;
}

/* 
 * Logic is like that:
 * The window sizes in w_height and w_width are sizes in pixels, and
 * correct pixels are 0..w_width-1.
 *
 * On the other hand, rect functions work with windows in which pixels
 * go 0..width. So we use some translation.
 */
void
recplothmult(long drawrect, long graphrect, entree *ep, GEN a, GEN b, char *ch, long prec, uLong flags)
{
  long av = avma,av2,i,j,sig,is,js,nc, offset = 0;
  GEN p1,p2,ysml,ybig,xsml,xbig,x,ydiff,xdiff,dx,y,dr, res;
  char c1[20], c2[20], c3[20], c4[20];
  
  PARI_get_plot();
  is = w_width - (lmargin + rmargin);
  js = w_height - (tmargin + bmargin);
  initrect(drawrect, w_width - 1, w_height - 1);
  if (graphrect != drawrect) initrect(graphrect, is - 1, js - 1);

  res = recplothmultin(graphrect, ep, a, b, ch, prec, flags, 0);

  xsml = (GEN)res[1];
  xbig = (GEN)res[2];
  ysml = (GEN)res[3];
  ybig = (GEN)res[4];

  /* Convert to floating point, then to double, then to string. */
  p1=cgetr(4); 
  gaffect(ybig,p1); sprintf(c1,"%9.3f",rtodbl(p1));
  gaffect(ysml,p1); sprintf(c2,"%9.3f",rtodbl(p1));
  gaffect(xsml,p1); sprintf(c3,"%9.3f",rtodbl(p1));
  gaffect(xbig,p1); sprintf(c4,"%9.3f",rtodbl(p1));

  rectlinetype(drawrect, -2);		/* Frame. */
  rect_string(drawrect, 0, f_height - 1, c1);
  rect_string(drawrect, 0, w_height - (bmargin + 2 * v_unit), c2);
  rect_string(drawrect, lmargin - (f_width*2), w_height - 1 - bmargin + f_height, c3);
  rect_string(drawrect, w_width - (f_width*10), w_height - 1 - bmargin + f_height, c4);

  dr = cgetg(7,17);
  dr[1] = lgeti(3); affsi(drawrect, (GEN)dr[1]);
  dr[2] = dr[3] = (long)gzero;
  dr[4] = lgeti(3); affsi(graphrect, (GEN)dr[4]);
  dr[5] = lgeti(3); affsi(lmargin, (GEN)dr[5]);
  dr[6] = lgeti(3); affsi(tmargin, (GEN)dr[6]);
  rectdraw(dr);
  avma = av;
  killvalue(ep);
  killrect(drawrect);
  if (graphrect != drawrect) killrect(graphrect);
}

GEN
rectplothraw(long drawrect, long graphrect, GEN listx, GEN listy, int flags)
{
  long av = avma,av2,i,lx,is,js;
  char c1[20],c2[20],c3[20],c4[20];
  GEN p1,xsml,xbig,ysml,ybig,dx,dy,scal,scaly,data,dr,xm,ym;

  if((typ(listx)<17)||(typ(listx)>18)||(typ(listy)<17)||(typ(listy)>18))
    err(ploter4);
  lx=lg(listx);
  if(lg(listy)!=lx) err(ploter5);
  if(lx==1) {return gnil;}

  PARI_get_plot();
  is = w_width - lmargin - rmargin;
  js = w_height - (tmargin + bmargin);

  xsml = vecmin(listx);
  ysml = vecmin(listy);
  xbig = vecmax(listx);
  ybig = vecmax(listy);
/*   av=avma; */
  dx=gdivgs(gsub(xbig,xsml),2);
  dy=gdivgs(gsub(ybig,ysml),2);
  xm = gadd(xsml, dx);
  ym = gadd(ysml, dy);
  if(gcmp0(dx) && gcmp0(dy)) {
    dx=gun;
    dy=gun;
  } else if (gcmp(gmulsg(is,dy),gmulsg(js,dx)) > 0) {
    dx=gdivgs(gmulsg(is,dy),js);
  } else {
    dy=gdivgs(gmulsg(js,dx),is);
  }
  xsml = gsub(xm,dx);
  xbig = gadd(xm,dx);
  ysml = gsub(ym,dy);
  ybig = gadd(ym,dy);

  initrect(drawrect, w_width, w_height);
  if (graphrect != drawrect) initrect(graphrect, is, js);
  av2=avma;
  data = cgetg(3,17);
  data[1] = (long) listx;
  data[2] = (long) listy;
  rectplothrawin(graphrect, xsml, xbig, 
	       ysml, ybig,
	       data, prec, 
	       PLOT_PARAMETRIC | (flags == PL_POINTS ? PLOT_POINTS : 0)
	       );
  avma=av2;

  /* Convert to floating point, then to double, then to string. */
  p1=cgetr(4); 
  gaffect(ybig,p1); sprintf(c1,"%9.3f",rtodbl(p1));
  gaffect(ysml,p1); sprintf(c2,"%9.3f",rtodbl(p1));
  gaffect(xsml,p1); sprintf(c3,"%9.3f",rtodbl(p1));
  gaffect(xbig,p1); sprintf(c4,"%9.3f",rtodbl(p1));

  rect_string(drawrect, 0, f_height - 1, c1);
  rect_string(drawrect, 0, w_height - (bmargin + 2 * v_unit), c2);
  rect_string(drawrect, lmargin - (f_width*2), w_height - bmargin + f_height - 1, c3);
  rect_string(drawrect, w_width - (f_width*10), w_height - bmargin + f_height - 1, c4);

  dr = cgetg(7,17);
  dr[1] = lgeti(3); affsi(drawrect, (GEN)dr[1]);
  dr[2] = dr[3] = (long)gzero;
  dr[4] = lgeti(3); affsi(graphrect, (GEN)dr[4]);
  dr[5] = lgeti(3); affsi(lmargin, (GEN)dr[5]);
  dr[6] = lgeti(3); affsi(tmargin, (GEN)dr[6]);
  rectdraw(dr);
  killrect(drawrect);
  if (graphrect != drawrect) killrect(graphrect);
  avma = av;
  return gun;
}


GEN
plothraw(GEN listx, GEN listy)
{
  return rectplothraw(DRAWRECT, GRAPHRECT, listx, listy, PL_POINTS);
}

GEN
plothrawlines(GEN listx, GEN listy)
{
  return rectplothraw(DRAWRECT, GRAPHRECT, listx, listy, 0);
}

GEN
ploth(entree *ep, GEN a, GEN b, char *ch, long prec)
{
  recplothmult(DRAWRECT, GRAPHRECT, ep,a,b,ch,prec,PLOT_SINGLE);
  return gun;
}


GEN
plothmult(entree *ep, GEN a, GEN b, char *ch, long prec)
{
  recplothmult(DRAWRECT, GRAPHRECT, ep,a,b,ch,prec,0);
  return gun;
}

GEN
ploth2(entree *ep, GEN a, GEN b, char *ch, long prec)
{
  recplothmult(DRAWRECT, GRAPHRECT, ep,a,b,ch,prec,PLOT_PARAMETRIC);
  return gun;
}

GEN
plothsizes()
{
  GEN vect = cgetg(1+6,17);
  int i;

  PARI_get_plot();
  for (i=1;i<7;i++) vect[i] = lgeti(3);
  affsi(w_width, (GEN)vect[1]);
  affsi(w_height, (GEN)vect[2]);
  affsi(h_unit, (GEN)vect[3]);
  affsi(v_unit, (GEN)vect[4]);
  affsi(f_width, (GEN)vect[5]);
  affsi(f_height, (GEN)vect[6]);
  return vect;
}	
