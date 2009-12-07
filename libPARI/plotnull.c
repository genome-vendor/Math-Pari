#include "genpari.h"

GEN
rectdraw(GEN list)
{
  return gnil;
}


void
PARI_get_plot()
{
  if (pari_plot.init) {
    return;
  }
  w_width = 2;
  w_height = 2; 
  f_height = 1;
  f_width = 1;
  h_unit = 1;
  v_unit = 1;
  pari_plot.init = 1;
}

long
term_set(char *s)
{
  return -1;
}
