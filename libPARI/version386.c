#include "genpari.h"

#ifdef __cplusplus
char *pariversion = "                    (C++ i386 version)\n";
#else
char *pariversion = "                      (i386 version)\n";
#endif

#include <sys/timeb.h>

long
timer()
{
  static long oldmsec;
  static long oldsec;
  long delay;
  struct timeb t;
  ftime(&t);
  delay = 1000 * (t.time - oldsec) + (t.millitm - oldmsec);
  oldmsec = t.millitm;
  oldsec = t.time;
  return delay;
}

long
timer2()
{
  static long oldmsec;
  static long oldsec;
  long delay;
  struct timeb t;
  ftime(&t);
  delay = 1000 * (t.time - oldsec) + (t.millitm - oldmsec);
  oldmsec = t.millitm;
  oldsec = t.time;
  return delay;
}

