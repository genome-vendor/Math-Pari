#include "genpari.h"

#ifdef __cplusplus
char *pariversion = "                   (C++ Sparcv7 version)\n";
#else
char *pariversion = "                     (Sparcv7 version)\n";
#endif

long
mulmodll(uLong a, uLong b, uLong c)
{
  divll(mulll(a,b),c);return hiremainder;
}

#ifdef SOLARIS

#include <sys/times.h>
#include <limits.h>

long
timer(void)
{
  static clock_t old_ticks;
  clock_t delay;
  struct tms t;
  times(&t);
  delay = (t.tms_utime - old_ticks) * (1000 / CLK_TCK);
  old_ticks = t.tms_utime;
  return (long) delay;
}

long
timer2(void)
{
  static clock_t old_ticks;
  clock_t delay;
  struct tms t;
  times(&t);
  delay = (t.tms_utime - old_ticks) * (1000 / CLK_TCK);
  old_ticks = t.tms_utime;
  return (long) delay;
}

#else

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
#endif




