#ifdef __cplusplus
#include        <stdlib.h>
#endif
#include        <stdarg.h>
#include        <setjmp.h>
#include        <signal.h>
#include        <stdio.h>
#include        <string.h>
#include        <math.h>
#include        <ctype.h>

/* uncomment the line below if you do not have the exp2 and log2 functions
   in your math library (or get one that has them), or better put
   * -DNOEXP2 into CC flags in Makefile. */

/* #define NOEXP2 */

#include "gencom.h"
#include "erreurs.h"
#include "gen68k.h"
#include "mpansi.h"

