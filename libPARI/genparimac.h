#ifndef __GENPARI__
#define __GENPARI__

#include <setjmp.h>
#include <signal.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <Types.h>
#include <StdLib.h>
#include <stdarg.h>

#ifdef LONG_IS_64BIT
#undef LONG_IS_64BIT
#endif

#ifndef macintosh
#define macintosh
#endif

#include "gencom.h"
#include "erreurs.h"

#if __MWERKS__
#include <SIOUX.h>
#include <Memory.h>
#define malloc(x) NewPtr(x)
#define calloc(x, y) NewPtrClear((x)*(y))
#define free(x) DisposePtr((Ptr)(x))
void *macrealloc(void *p, size_t olds, size_t news);
#define CodeWarrior_Bug
#pragma pointers_in_D0
#endif

#if powerc
#include "genport.h"
#else
#include "gen68k.h"
#endif

#include "mpansi.h"

#if __MWERKS__
#pragma pointers_in_A0
#endif

#endif
