#include "genpari.h"

#ifdef __cplusplus
char *pariversion = "             (C++ PowerPC 601 version)\n";
#else
char *pariversion = "                 (PowerPC 601 version)\n";
#endif

uLong overflow,hiremainder;

long shiftl(uLong x, uLong y)
{
  hiremainder=x>>(BITS_IN_LONG-y);return (x<<y);
}

long shiftlr(uLong x, uLong y)
{
  hiremainder=x<<(BITS_IN_LONG-y);return (x>>y);
}

long mulmodll(uLong a,uLong b,uLong c)
{
  divll(mulll(a,b),c);
  return hiremainder;
}

#ifdef macintosh

#include <Events.h>
     
long
timer(void)
{
  static long oldticks;
  long ticks = TickCount();
  long delay = ticks - oldticks;
  oldticks = ticks;
  return 50 * delay / 3;
}

long
timer2(void)
{
  static long oldticks;
  long ticks = TickCount();
  long delay = ticks - oldticks;
  oldticks = ticks;
  return 50 * delay / 3;
}

#endif