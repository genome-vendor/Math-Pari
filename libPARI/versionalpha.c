#include "genpari.h"

#ifdef __cplusplus
char *pariversion = "                 (C++ Alpha 64-bit version)\n";
#else
char *pariversion = "                   (Alpha 64-bit version)\n";
#endif

long
mulmodll(uLong a, uLong b, uLong c)
{
  divll(mulll(a,b),c);return hiremainder;
}

#include <sys/time.h>
#include <sys/resource.h>

long
timer(void)
{
  static long oldmusec;
  static long oldsec;
  long delay;
  struct rusage r;
  struct timeval t;
  getrusage(0,&r);t=r.ru_utime;
  delay = 1000 * (t.tv_sec - oldsec) + (t.tv_usec - oldmusec) / 1000;
  oldmusec = t.tv_usec;
  oldsec = t.tv_sec;
  return delay;
}

long
timer2(void)
{
  static long oldmusec;
  static long oldsec;
  long delay;
  struct rusage r;
  struct timeval t;
  getrusage(0,&r);t=r.ru_utime;
  delay = 1000 * (t.tv_sec - oldsec) + (t.tv_usec - oldmusec) / 1000;
  oldmusec = t.tv_usec;
  oldsec = t.tv_sec;
  return delay;
}

