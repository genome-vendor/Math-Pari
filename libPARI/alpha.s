	.text	

        .set noreorder
        .align 3
        .globl addll
        .ent addll 0
addll:	
        .frame $30,0,$26,0
        .prologue 0
	addq	$16,$17,$0
	cmpult	$0,$16,$1
	stq	$1,overflow
	ret	$31,($26),1
	.end	addll			

        .set noreorder
        .align 3
        .globl addllx
        .ent addllx 0
addllx:	
        .frame $30,0,$26,0
        .prologue 0
	ldq	$2,overflow
	addq	$16,$17,$3
	addq	$2,$3,$0
	cmpeq	$0,$16,$4
	cmpult	$0,$16,$1
	and	$4,$2,$5
	or	$5,$1,$6
	stq	$6,overflow
	ret	$31,($26),1
	.end	addllx
	
        .set noreorder
        .align 3
        .globl subll
        .ent subll 0
subll:	
        .frame $30,0,$26,0
        .prologue 0
	subq	$16,$17,$0
	cmpult	$16,$0,$1
	stq	$1,overflow
	ret	$31,($26),1
	.end	subll			

        .set noreorder
        .align 3
        .globl subllx
        .ent subllx 0
subllx:	
        .frame $30,0,$26,0
        .prologue 0
	ldq	$2,overflow
	subq	$16,$17,$3
	subq	$3,$2,$0
	cmpeq	$0,$16,$4
	cmpult	$16,$0,$1
	and	$4,$2,$5
	or	$5,$1,$6
	stq	$6,overflow
	ret	$31,($26),1
	.end	subllx

        .set noreorder
        .align 3
        .globl shiftl
        .ent shiftl 0
shiftl:	
        .frame $30,0,$26,0
        .prologue 0
	subq	$31,$17,$1
	sll	$16,$17,$0
	srl	$16,$1,$2
	stq	$2,hiremainder
	ret	$31,($26),1
	.end	shiftl

        .set noreorder
        .align 3
        .globl shiftlr
        .ent shiftlr 0
shiftlr:	
        .frame $30,0,$26,0
        .prologue 0
	subq	$31,$17,$1
	srl	$16,$17,$0
	sll	$16,$1,$2
	stq	$2,hiremainder
	ret	$31,($26),1
	.end	shiftlr

        .set noreorder
        .align 3
        .globl bfffo
        .ent bfffo 0
bfffo:	
        .frame $30,0,$26,0
        .prologue 0
	lda	$6,tabshi
	ldiq	$7,60
	and	$16,0xffffffff00000000,$1
	beq	$1,$32
	srl	$16,32,$16
	subq	$7,32,$7
$32:	cmpule	$16,0xffff,$2
	bne	$2,$33
	srl	$16,16,$16
	subq	$7,16,$7
$33:	cmpule	$16,0xff,$3
	bne	$3,$34
	srl	$16,8,$16
	subq	$7,8,$7
$34:	cmpule	$16,0xf,$4
	bne	$4,$35
	srl	$16,4,$16
	subq	$7,4,$7
$35:	s8addq	$16,$6,$5
	ldq	$1,0($5)
	addq	$1,$7,$0
	ret	$31,($26),1
	.end	bfffo
	
        .set noreorder
        .align 3
        .globl mulll
        .ent mulll 0
mulll:	
        .frame $30,0,$26,0
        .prologue 0
	umulh	$16,$17,$1
	mulq	$16,$17,$0
	stq	$1,hiremainder
	ret	$31,($26),1
	.end	mulll

        .set noreorder
        .align 3
        .globl addmul
        .ent addmul 0
addmul:	
        .frame $30,0,$26,0
        .prologue 0
	mulq	$16,$17,$2
	umulh	$16,$17,$1
	ldq	$3,hiremainder
	addq	$2,$3,$0
	cmpult	$0,$2,$4
	addq	$1,$4,$5
	stq	$5,hiremainder
	ret	$31,($26),1
	.end	addmul
		
 # This program is a modification of a file contained in the gmp-1.9
 # library, copyright Free Software Foundation, with permission.
	.globl	err
        .set noreorder
        .align 3
        .globl divll
        .ent divll 0
divll:	
        .frame $30,0,$26,0
        .prologue 0
#define cnt	$2
#define tmp	$3
#define n1	$7
#define n0	$16
#define d	$17
#define qb	$20

	ldq	n1,hiremainder
	ldiq	cnt,16
	cmpule	d,n1,tmp
	bne	tmp,errorhandler
	blt	d,Largedivisor

Loop1:	cmplt	n0,0,tmp
	addq	n1,n1,n1
	bis	n1,tmp,n1
	addq	n0,n0,n0
	cmpule	d,n1,qb
	subq	n1,d,tmp
	cmovne	qb,tmp,n1
	bis	n0,qb,n0
	cmplt	n0,0,tmp
	addq	n1,n1,n1
	bis	n1,tmp,n1
	addq	n0,n0,n0
	cmpule	d,n1,qb
	subq	n1,d,tmp
	cmovne	qb,tmp,n1
	bis	n0,qb,n0
	cmplt	n0,0,tmp
	addq	n1,n1,n1
	bis	n1,tmp,n1
	addq	n0,n0,n0
	cmpule	d,n1,qb
	subq	n1,d,tmp
	cmovne	qb,tmp,n1
	bis	n0,qb,n0
	cmplt	n0,0,tmp
	addq	n1,n1,n1
	bis	n1,tmp,n1
	addq	n0,n0,n0
	cmpule	d,n1,qb
	subq	n1,d,tmp
	cmovne	qb,tmp,n1
	bis	n0,qb,n0
	subq	cnt,1,cnt
	bgt	cnt,Loop1
	stq	n1,hiremainder
	bis	$31,n0,$0
	ret	$31,($26),1

Largedivisor:
	and	n0,1,$4

	srl	n0,1,n0
	sll	n1,63,tmp
	or	tmp,n0,n0
	srl	n1,1,n1

	and	d,1,$6
	srl	d,1,$5
	addq	$5,$6,$5

Loop2:	cmplt	n0,0,tmp
	addq	n1,n1,n1
	bis	n1,tmp,n1
	addq	n0,n0,n0
	cmpule	$5,n1,qb
	subq	n1,$5,tmp
	cmovne	qb,tmp,n1
	bis	n0,qb,n0
	cmplt	n0,0,tmp
	addq	n1,n1,n1
	bis	n1,tmp,n1
	addq	n0,n0,n0
	cmpule	$5,n1,qb
	subq	n1,$5,tmp
	cmovne	qb,tmp,n1
	bis	n0,qb,n0
	cmplt	n0,0,tmp
	addq	n1,n1,n1
	bis	n1,tmp,n1
	addq	n0,n0,n0
	cmpule	$5,n1,qb
	subq	n1,$5,tmp
	cmovne	qb,tmp,n1
	bis	n0,qb,n0
	cmplt	n0,0,tmp
	addq	n1,n1,n1
	bis	n1,tmp,n1
	addq	n0,n0,n0
	cmpule	$5,n1,qb
	subq	n1,$5,tmp
	cmovne	qb,tmp,n1
	bis	n0,qb,n0
	subq	cnt,1,cnt
	bgt	cnt,Loop2

	addq	n1,n1,n1
	addq	$4,n1,n1
	bne	$6,Odd
	stq	n1,hiremainder
	bis	$31,n0,$0
	ret	$31,($26),1

Odd:
	! q' in n0. r' in n1
	addq	n1,n0,n1
	cmpult	n1,n0,tmp	# tmp := carry from addq
	beq	tmp,LLp6
	addq	n0,1,n0
	subq	n1,d,n1
LLp6:	cmpult	n1,d,tmp
	bne	tmp,LLp7
	addq	n0,1,n0
	subq	n1,d,n1
LLp7:
	stq	n1,hiremainder
	bis	$31,n0,$0
	ret	$31,($26),1

errorhandler:
	ldiq	$16,0x2f
	jmp	err
	
	.end	divll

	.align	3
tabshi:	.quad	4,3,2,2,1,1,1,1,0,0,0,0,0,0,0,0
	.globl	hiremainder
	.comm	hiremainder,4
	.globl	overflow
	.comm	overflow,4
