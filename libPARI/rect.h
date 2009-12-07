
typedef struct RectObj {
  struct RectObj *next;
  long code;
} RectObj;

typedef struct Rect {
  RectObj *head;
  RectObj *tail;
  long sizex;
  long sizey;
  GEN cursorx;
  GEN cursory;
  GEN xscale;
  GEN xshift;
  GEN yscale;
  GEN yshift;
} Rect;

#define ROt_MV 0			/* Move */
#define ROt_PT 1			/* Point */
#define ROt_LN 2			/* Line */
#define ROt_BX 3			/* Box */
#define ROt_MP 4			/* Multiple point */
#define ROt_ML 5			/* Multiple lines */
#define ROt_ST 6			/* String */
#define ROt_PTT 7			/* Point type change */
#define ROt_LNT 8			/* Line type change */
#define ROt_NULL 9			/* To be the start of the chain */

#define ROt_MAX 9			/* Maximal type */

/* The structures below are "subclasses" of RectObj. */

typedef struct RectObj1P {
  struct RectObj *next;
  long code;
  long x;
  long y;
} RectObj1P;

typedef struct RectObj2P {
  struct RectObj *next;
  long code;
  long x1;
  long y1;
  long x2;
  long y2;
} RectObj2P;

typedef struct RectObjMP {
  struct RectObj *next;
  long code;
  long count;
  long *xs;
  long *ys;
} RectObjMP;

typedef struct RectObjST {
  struct RectObj *next;
  long code;
  long length;
  char *s;
  long x;
  long y;
} RectObjST;

typedef struct RectObjPN {
  struct RectObj *next;
  long code;
  long pen;
} RectObjPN;

/* Pointer conversion. */

#define RoMV(rop) ((RectObj1P*)rop)
#define RoPT(rop) ((RectObj1P*)rop)
#define RoLN(rop) ((RectObj2P*)rop)
#define RoBX(rop) ((RectObj2P*)rop)
#define RoMP(rop) ((RectObjMP*)rop)
#define RoML(rop) ((RectObjMP*)rop)
#define RoST(rop) ((RectObjST*)rop)
#define RoPTT(rop) ((RectObjPN*)rop)
#define RoLNT(rop) ((RectObjPN*)rop)
#define RoNULL(rop) ((RectObj*)rop)

/* All the access to the rectangle data _should_ go via these macros! */

#define RHead(rp) ((rp)->head)
#define RTail(rp) ((rp)->tail)
#define RXsize(rp) ((rp)->sizex)
#define RYsize(rp) ((rp)->sizey)
#define RXcursor(rp) ((rp)->cursorx)
#define RYcursor(rp) ((rp)->cursory)
#define RXshift(rp) ((rp)->xshift)
#define RYshift(rp) ((rp)->yshift)
#define RXscale(rp) ((rp)->xscale)
#define RYscale(rp) ((rp)->yscale)


#define RoNext(rop) ((rop)->next)
#define RoType(rop) ((rop)->code)
#define RoMVx(rop) (RoMV(rop)->x)
#define RoMVy(rop) (RoMV(rop)->y)
#define RoPTx(rop) (RoPT(rop)->x)
#define RoPTy(rop) (RoPT(rop)->y)
#define RoLNx1(rop) (RoLN(rop)->x1)
#define RoLNy1(rop) (RoLN(rop)->y1)
#define RoLNx2(rop) (RoLN(rop)->x2)
#define RoLNy2(rop) (RoLN(rop)->y2)
#define RoBXx1(rop) (RoBX(rop)->x1)
#define RoBXy1(rop) (RoBX(rop)->y1)
#define RoBXx2(rop) (RoBX(rop)->x2)
#define RoBXy2(rop) (RoBX(rop)->y2)

#define RoMPcnt(rop) (RoMP(rop)->count)
#define RoMPxs(rop) (RoMP(rop)->xs)
#define RoMPys(rop) (RoMP(rop)->ys)

#define RoMLcnt(rop) (RoML(rop)->count)
#define RoMLxs(rop) (RoML(rop)->xs)
#define RoMLys(rop) (RoML(rop)->ys)

#define RoSTs(rop) (RoST(rop)->s)
#define RoSTl(rop) (RoST(rop)->length)
#define RoSTx(rop) (RoST(rop)->x)
#define RoSTy(rop) (RoST(rop)->y)

#define RoPTTpen(rop) (RoPTT(rop)->pen)
#define RoLNTpen(rop) (RoLNT(rop)->pen)

#define NUMRECT 18
#define GOODRECT(r) (0 <= r && r < NUMRECT)

#define PL_POINTS 1


#define PLOT_PARAMETRIC 1
#define PLOT_SINGLE 2
#define PLOT_NO_AXE_X 8
#define PLOT_NO_AXE_Y 16
#define PLOT_NO_FRAME 32
#define PLOT_POINTS 64

extern  Rect    **rectgraph;
GEN rectlinetype(long,long);
GEN rectpointtype(long,long);
