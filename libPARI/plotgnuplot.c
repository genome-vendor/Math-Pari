/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                 */
/*                     PLOT EN HAUTE RESOLUTION                    */
/*                                                                 */
/*                       copyright Babe Cool                       */
/*                       (C) Ilya Zakharevich                      */
/*                                                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

# include "genpari.h"
#include "rect.h"
#define croak(str) err(talker,str)
#include "Gnuplot.h"

#ifdef __EMX__
#  define DEF_TERM "pm"
#else
#  define DEF_TERM "X11"
#endif

GEN
rectdraw(GEN list)
{
  long *ptx,*pty;
  long n,i,j,x0,y0,av=avma;
  long a,b,ne, good;
  char **texts;
  int point_type = -1, line_type = 0;
  Rect *e;
  RectObj *p1;

  if(typ(list)!=17) err(rploter3);
  n=lg(list)-1;if(n%3) err(rploter4);
  n=n/3;if(!n) {abort();return gnil;}

  PARI_get_plot();

  graphics();				/* Switch on terminal. */
  linetype(line_type);			/* X does not work otherwise. */
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
	  point(RoPTx(p1)+x0, w_height - 1 - RoPTy(p1) - y0, point_type);
	  break;
	case ROt_LN:
	  move(RoLNx1(p1)+x0, w_height - 1 - RoLNy1(p1) - y0);
	  vector(RoLNx2(p1)+x0, w_height - 1 - RoLNy2(p1) - y0);
	  break;
	case ROt_BX:
	  move(RoBXx1(p1)+x0, w_height - 1 - RoBXy1(p1) - y0);
	  vector(RoBXx2(p1)+x0, w_height - 1 - RoBXy1(p1) - y0);
	  vector(RoBXx2(p1)+x0, w_height - 1 - RoBXy2(p1) - y0);
	  vector(RoBXx1(p1)+x0, w_height - 1 - RoBXy2(p1) - y0);
	  vector(RoBXx1(p1)+x0, w_height - 1 - RoBXy1(p1) - y0);
	  break;
	case ROt_MP:
	  ptx=RoMPxs(p1);
	  pty=RoMPys(p1);
	  for(j=0;j<RoMPcnt(p1);j++)
	  {
	    point(ptx[j]+x0,  w_height - 1 - pty[j] - y0, point_type);
	  }
	  break;
	case ROt_ML:
	  ptx=RoMLxs(p1);
	  pty=RoMLys(p1);
	  j = 0;
	  if (ptx[j]+x0 < 0 || ptx[j]+x0 >= w_width
	      || pty[j] + y0 < 0 || pty[j] + y0 >= w_height) {
	    good = 0;
	  } else {
	    move(ptx[j]+x0, w_height - 1 - pty[j] - y0);
	    good = 1;
	  }
	  for(j=1;j<RoMLcnt(p1);j++)
	  {
	    if (good) {
	      if (ptx[j]+x0 < 0 || ptx[j]+x0 >= w_width
		  || pty[j] + y0 < 0 || pty[j] + y0 >= w_height) {
		good = 0;
	      } else {
		vector(ptx[j]+x0, w_height - 1 - pty[j] - y0);
	      }
	    } else {
	      if (ptx[j]+x0 < 0 || ptx[j]+x0 >= w_width
		  || pty[j] + y0 < 0 || pty[j] + y0 >= w_height) {
	      } else {
		move(ptx[j]+x0, w_height - 1 - pty[j] - y0);
		good = 1;
	      }
	    }
	  }
	  break;
	case ROt_ST:
	  if (RoSTx(p1)+x0 < 0 || RoSTx(p1)+x0+RoSTl(p1)-1 >= w_width
	      || RoSTy(p1) + y0 < 0 || RoSTy(p1) + y0 >= w_height) {
	  } else {
	    put_text(RoSTx(p1)+x0,  w_height - 1 - RoSTy(p1) - y0, RoSTs(p1));
	  }
	  break;
	case ROt_PTT:
	  point_type = RoPTTpen(p1);
	  break;
	case ROt_LNT:
	  linetype(RoLNTpen(p1));
	  break;
	default: break;
      }
      p1=RoNext(p1);
    }
  }

  text();				/* Reset terminal */
  avma = av;return gnil;
}

void
PARI_get_plot()
{
  if (pari_plot.init) {
    return;
  }
  term_set( DEF_TERM );
}


long
term_set(char *s)
{
  if (strlen(s) > PLOT_NAME_LEN) err(talker,"too long name for terminal");
  if (*pari_plot.name && (strcmp(pari_plot.name,s) != 0)) {
	reset();
  }
  strcpy(pari_plot.name,s);
  if (!termset( s )) err(talker,"unknown terminal name");
  init();				/* Init terminal. */

  w_width = termprop(xmax);
  w_height = termprop(ymax);
  f_height = termprop(v_char);
  f_width = termprop(h_char);
  h_unit = termprop(h_tic);
  v_unit = termprop(v_tic);
  pari_plot.init = 1;
  return 1;
}
