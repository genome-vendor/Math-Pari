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

#include        <suntool/sunview.h>
#include        <suntool/canvas.h>
#include        <suntool/textsw.h>
#include        <suntool/panel.h>

typedef struct spoint {
  int x,y;} SPoint; 
typedef struct ssegment {
  int x1,y1,x2,y2;} SSegment;
typedef struct srectangle {
  int x,y,width,height;} SRectangle;

#define ISCR 1120 /* 1400 en haute resolution */     
#define JSCR 800  /* 1120 en haute resolution */     


GEN
rectdraw(GEN list)
{
  long *e,*p1,*ptx,*pty,*numpoints,*numtexts,*xtexts,*ytexts;
  long n,i,j,x0,y0,av=avma;
  long a,b,c,d,nd[10],ne;
  char **texts;

  Frame ecran;
  Canvas canevas;
  Pixwin *pw;
  Pixfont *font;
  SPoint *points, **lines, *SLine;
  SSegment *segments; 
  SRectangle *rectangles, SRec;

  if(typ(list)!=17) err(rploter3);
  n=lg(list)-1;if(n%3) err(rploter4);
  n=n/3;if(!n) return gnil;
  nd[0]=nd[1]=nd[2]=nd[3]=nd[4]=nd[5]=nd[6]=0;
  for(i=0;i<n;i++)
  {
    if(typ(list[3*i+1])!=1) err(rploter5);
    ne=itos(list[3*i+1]);if((ne<0)||(ne>15)) err(rploter2);
    e=((long**)rectgraph)[ne];
    p1=(long*)e[0];while((long)p1) 
    {
      if(p1[1]!=4) nd[p1[1]]++;
      else nd[1]+=p1[2];
      p1=(long*)p1[0];
    }
  }
  points=(SPoint*)malloc(nd[1]*sizeof(SPoint));
  segments=(SSegment*)malloc(nd[2]*sizeof(SSegment));
  rectangles=(SRectangle*)malloc(nd[3]*sizeof(SRectangle));
  lines=(SPoint**)malloc(nd[5]*sizeof(SPoint*));
  numpoints=(long*)malloc(nd[5]*sizeof(long));
  texts=(char**)malloc(nd[6]*sizeof(char*));
  numtexts=(long*)malloc(nd[6]*sizeof(long));
  xtexts=(long*)malloc(nd[6]*sizeof(long));
  ytexts=(long*)malloc(nd[6]*sizeof(long));
  nd[1]=nd[2]=nd[3]=nd[5]=nd[6]=0;
  for(i=0;i<n;i++)
  {
    e=((long**)rectgraph)[itos(list[3*i+1])];x0=list[3*i+2];y0=list[3*i+3];
    if((typ(x0)!=1)||(typ(y0)!=1)) err(rploter5);
    x0=itos(x0);y0=itos(y0);
    p1=(long*)e[0];
    while((long)p1)
    {
      switch(p1[1])
      {
	case 1: 
	  points[nd[1]].x=p1[2]+x0;
	  points[nd[1]].y=p1[3]+y0;
	  nd[1]++;break;
	case 2:
	  segments[nd[2]].x1=p1[2]+x0;
	  segments[nd[2]].y1=p1[3]+y0;
	  segments[nd[2]].x2=p1[4]+x0;
	  segments[nd[2]].y2=p1[5]+y0;
	  nd[2]++;break;
	case 3:
	  a=rectangles[nd[3]].x=p1[2]+x0;
	  b=rectangles[nd[3]].y=p1[3]+y0;
	  rectangles[nd[3]].width=p1[4]+x0-a;
	  rectangles[nd[3]].height=p1[5]+y0-b;
	  nd[3]++;break;
	case 4:
	  ptx=(long*)p1[3];pty=(long*)p1[4];
	  for(j=0;j<p1[2];j++)
	  {
	    points[nd[1]+j].x=ptx[j]+x0;
	    points[nd[1]+j].y=pty[j]+y0;
	  }
	  nd[1]+=p1[2];break;
	case 5:
	  ptx=(long*)p1[3];pty=(long*)p1[4];
	  numpoints[nd[5]]=p1[2];
	  lines[nd[5]]=(SPoint*)malloc(p1[2]*sizeof(SPoint));
	  for(j=0;j<p1[2];j++)
	  {
	    lines[nd[5]][j].x=ptx[j]+x0;
	    lines[nd[5]][j].y=pty[j]+y0;
	  }
	  nd[5]++;break;
	case 6: 
	  texts[nd[6]]=(char*)p1[3];numtexts[nd[6]]=p1[2];
	  xtexts[nd[6]]=p1[4]+x0;ytexts[nd[6]]=p1[5]+y0;
	  nd[6]++;break;
	default: break;
      }
      p1=(long*)p1[0];
    }
  }
  ecran=window_create(NULL,FRAME,FRAME_LABEL,"rectplot",
                      WIN_ERROR_MSG,"you must be in suntools",0);
  canevas=window_create(ecran,CANVAS,WIN_HEIGHT,JSCR,
                        WIN_WIDTH,ISCR,0);
  window_fit(ecran);pw=canvas_pixwin(canevas);

  font=pw_pfsysopen();
  for(i=0;i<nd[1];i++) pw_put(pw,points[i].x,points[i].y,1);
  for(i=0;i<nd[2];i++) pw_vector(pw,segments[i].x1,segments[i].y1,segments[i].x2,segments[i].y2,PIX_SRC,1);
  for(i=0;i<nd[3];i++) 
  {
    SRec=rectangles[i];a=SRec.x;b=SRec.y;c=a+SRec.width;
    d=b+SRec.height;
    pw_vector(pw,a,b,c,b,PIX_SRC,1);
    pw_vector(pw,c,b,c,d,PIX_SRC,1);
    pw_vector(pw,a,d,c,d,PIX_SRC,1);
    pw_vector(pw,a,b,a,d,PIX_SRC,1);
  }
  for(i=0;i<nd[5];i++) 
  {
    SLine=lines[i];
    for(j=1;j<numpoints[i];j++)
      pw_vector(pw,SLine[j-1].x,SLine[j-1].y,SLine[j].x,SLine[j].y,PIX_SRC,1);
  }
  for(i=0;i<nd[6];i++) 
    for(j=0;texts[i][j];j++)
      pw_char(pw,xtexts[i]+9*j,ytexts[i],PIX_SRC|PIX_DST,font,texts[i][j]);
  window_main_loop(ecran);
  free(points);free(segments);free(rectangles);
  free(numpoints);for(i=0;i<nd[5];i++) free(lines[i]);
  free(lines);free(texts);free(numtexts);free(xtexts);free(ytexts);
  avma = av;return gnil;
}

void
PARI_get_plot()
{
  Display *display;
  int screen;

  if (pari_plot.init) {
    return;
  }

  w_width = ISCR;
  w_height = JSCR; 
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
