global findHighLow
global updateErrorCache
segment .data
fHigh:		dd	10000000.0
fLow:		dd	-1000000.0
segment .text

;;
;load value from array to FPU
;@param pointer to array
;@param position in array
;;
%macro	arrayToFPU 2
	push	edx
	mov	edx,%1
	fld	dword[edx+4*%2]
	pop	edx
%endmacro

;;
;load value from array
;@param pointer to array
;@param pointer to destination
;@param position in array
;;
%macro	loadArrayPosition 3
	push	edx
	mov	edx,%1
	mov	dword %2, dword[edx+4*%3]
	pop	edx
%endmacro

;;
;macro evalutes expression:
;|first-second| <= error
;eax = 1 if true; otherwise eax = 0
;@param first operand
;@param second operand
;@param error
;;
%macro equalsWithTolerance 3
	fsub	%1,%2
	fabs
	fcomi	%1,%3
	jbe	%%equal
	jmp	%%notEqual
%%equal:
	mov	eax,1
	jmp	%%finish
%%notEqual:
	mov	eax,0
%%finish:
	nop
%endmacro

;;
;macro evalutes kernel
;caching enabled
;@param number of position in dest array
;@param cached kernel array offset in edx
;@param X[second] offset in ebx
;;
%macro getKernel 3
	add	eax,1
	push	eax
	mov	ecx,dword [edx+%2]	; cachedKernel
	imul	eax,4
	add	eax,dword [edx+64]	; const offset
	add	ecx,eax			; cachedKernel[i]
	pop	eax
	fld	dword [ecx]
	fld1
	fchs
	fcomi	ST0,ST1			; cachedKernel == -1?
	fstp	ST0
	fstp	dword [ebx+4*%1]
	jne	%%finish
	lea	esp,[esp-4]		; result
	push	dword [ebx+32]		; number of features
	mov	ecx,dword [ebx+32]	; number of features
	imul	ecx,eax			; x offset
	imul	ecx,4
	add	edi,ecx
	push	edi			; X[i]
	sub	edi,ecx
	push	dword [ebx+%3]		; X[iHigh]
	push	esp			; structure address
	call	computeLinearKernel
	lea	esp,[esp+16]
	push	eax
	mov	eax,[esp+4]
	mov	dword [ebx+4*%1],eax	; high kernel 0
	pop	eax
	push	eax
	fld	dword [ebx+4*%1]
	mov	ecx,dword [edx+%2]
	imul	eax,4
	add	eax,dword [edx+64]	; const offset
	add	ecx,eax
	fstp	dword [ecx]		; save cached kernel
	pop	eax
	lea	esp,[esp+4]
%%finish:
	nop
%endmacro

;;
;finds iHigh and iLow in alpha array
;@param ParallelAsmData structure
;;
findHighLow:
   	push    ebp	;save ebp
	mov     ebp,esp	
	push	ebx
	xor	eax,eax
	mov	edx,[ebp+8]	;structure
	mov	ecx,[edx]	;trainDataSize
	push	dword	0	; temp_esp [ebx+24]
	push	dword	[fHigh] ; fHigh [ebx+20]
	push	dword	[fLow]	; fLow [ebx+16]
	push	dword	-1	; iHigh [ebx+12]
	push	dword	-1	; iLow [ebx+8]
	push	dword	0	; lowFlag [ebx+4]
	push	dword	0	; highFlag [ebx]
	mov	ebx,esp		; stack pointer in ebx
	mov	[ebx+24],esp
	finit
	fld	dword [edx+16]	;stack: cost
	fld	dword [edx+20]	;stack: error,cost
;for(i = 0;i < trainDataSize)
LOOP:	
	mov	eax,[edx]
	sub	eax,ecx
	arrayToFPU	[edx+24],eax ;stack: y(i),error,cost
	arrayToFPU	[edx+28],eax ;stack: alphas(i),y(i),error,cost
	mov	dword [ebx+4],0
	mov	dword [ebx],0
	
;check if error < alpha < cost-error
middle:	
	clc
	fldz				;stack: 0,alphas(i),y(i),error,cost
	fadd	ST0,ST3			;add error to low border
	fcomi	ST0,ST1			;compare low+error > val
	fstp	ST0			;stack: alphas(i),y(i),error,cost
	jb	is			;if ST0 < ST1 - may be in interval
	jmp	alpha0			;second condition
;if alpha > error
is:	
	fld	ST3			;stack: cost,alphas(i),y(i),error,cost
	fsub	ST0,ST3			;cost = cost-error
	fcomi	ST0,ST1			;compare high-error > val
	fstp	ST0			;stack: alphas(i),y(i),error,cost
	jb	alpha0			;if ST0 < ST1 - is not in interval
	mov	eax,1
	mov	dword [ebx],eax		;check high condition
	mov	dword [ebx+4],eax	;check low condition
	jmp	highCheck		;check conditions
	
	
;check if 0 <= alpha <= error
alpha0:	
	;only testAlphaMiddle shoudl pass!
	;jmp	highCheck
	fldz				;stack: 0,alphas(i),y(i),error,cost
	equalsWithTolerance ST0,ST1,ST3 ;check if alphas(i) <= error
	cmp	eax,0
	fstp	ST0			;stack: alphas(i),y(i),error,cost
	je	alphaCost		;alphas(i) > error
	fld1				;stack: 0,alphas(i),y(i),error,cost
	fcomi	ST0,ST2
	fstp	ST0			;stack: alphas(i),y(i),error,cost
	je	alpha0YPos		;y(i) > 0
;if y(i) < 0 -> check iLow
alpha0YNeg:
	mov	eax,1
	mov	dword [ebx+4],eax	;check iLow
	jmp	alphaCost
;if y(i) > 0 -> check iHigh
alpha0YPos:
	mov	eax,1
	mov	dword [ebx],eax		;check iHigh
	jmp	alphaCost
	
;check if cost-error <= alpha <= cost
alphaCost:
	fld	ST3			;stack: cost,alphas(i),y(i),error,cost
	equalsWithTolerance ST0,ST1,ST3 ;check if |alphas(i)-cost| <= error
	cmp	eax,0
	fstp	ST0			;stack: alphas(i),y(i),error,cost
	je	highCheck		;alphas(i) < cost-error
	fld1				;stack: 0,alphas(i),y(i),error,cost
	fcomi	ST0,ST2
	fstp	ST0			;stack: alphas(i),y(i),error,cost
	
	je	alphaCostYPos		;y(i) > 0
;if y(i) < 0 -> check iHigh
alphaCostYNeg:
	mov	eax,1
	mov	dword [ebx],eax		;check iHigh
	jmp	highCheck
;if y(i) > 0 -> check iLow
alphaCostYPos:
	mov	eax,1
	mov	dword [ebx+4],eax	;check iLow
	jmp	highCheck

highCheck:
	mov	eax,dword [ebx]		;check if flag is set
	cmp	eax,0
	je	lowCheck		;highCheck flag = 0
	mov	eax,[edx]
	sub	eax,ecx
	arrayToFPU [edx+32],eax		;stack: error(i),alphas(i),y(i),error,cost
	fld	dword [ebx+20]		;stack: fHigh,error(i),alphas(i),y(i),error,cost
	fcomi	ST0,ST1			;compare fHigh and error(i)
	jnbe	highCheckAssign		;set new value
	fstp	ST0			;stack: error(i),alphas(i),y(i),error,cost
	fstp	ST0			;stack: alphas(i),y(i),error,cost
	jmp	lowCheck		;check fLow	
highCheckAssign:
	mov	eax,[edx]
	sub	eax,ecx
	mov	dword [ebx+12],eax	;iHigh = i
	fstp	ST0			;stack: error(i),alphas(i),y(i),error,cost
	fstp	dword [ebx+20]		;fHigh = ST0 | stack: alphas(i),y(i),error,cost
	
lowCheck:
	mov	eax,dword [ebx+4]	;check if flag is set
	cmp	eax,0
	je	endLoop			;lowCheck flag = 0
	mov	eax,[edx]
	sub	eax,ecx			;get number of iteration
	arrayToFPU [edx+32],eax		;stack: error(i),alphas(i),y(i),error,cost
	fld	dword [ebx+16]		;stack: fLow,error(i),alphas(i),y(i),error,cost
	fcomi	ST0,ST1			;compare fLow and error(i)
	jb	lowCheckAssign		;set new value
	fstp	ST0			;stack: error(i),alphas(i),y(i),error,cost
	fstp	ST0			;stack: alphas(i),y(i),error,cost
	jmp	endLoop			;end iteration	
lowCheckAssign:
	mov	eax,[edx]
	sub	eax,ecx
	mov	dword [ebx+8],eax	;iLow = i
	fstp	ST0			;stack: error(i),alphas(i),y(i),error,cost
	fstp	dword [ebx+16]		;fLow = ST0 | stack: alphas(i),y(i),error,cost
	
endLoop:
	fstp	ST0			;stack: y(i),error,cost
	fstp	ST0			;stack: error,cost
	dec	ecx			;counter--
	jnz	LOOP			;if counter = 0 -> end
	fild	dword [ebx+12]
	fistp	dword [edx+4]		;ptr->iHigh = iHigh
	fild	dword [ebx+8]
	fistp	dword [edx+8]		;ptr->iLow = iLow
end:	
	mov	esp,[ebx+24]
	lea	esp,[esp+28]
	pop	ebx
	pop     ebp
	ret
	
;;
;computes linear kernel
;result in eax
;@param pointer to structure: first,second,size,result
;;
computeLinearKernel:
   	push    ebp			;save ebp
	mov     ebp,esp	
	push	ebx
	push	ecx
	push	edx
	push	eax			;save registers
	push	edi
	push	esi
	xor	eax,eax
	mov	edx,[ebp+8]		;structure
	lea	esp,[esp-24]		;four numbers+sum+temp_esp
	mov	ebx,esp			; stack pointer in ebx
	mov	[ebx+20],esp		; save esp
	mov	esi,[edx]		; first_ptr in esi
	mov	edi,[edx+4]		; second_ptr in edi
	xorps	xmm0,xmm0		; 0 in xmm0
	mov	dword [ebx+16],0	; 0 in temp_sum
	mov	ecx,dword [edx+8]	; get size
	cmp	ecx,12			; check if size >= 12
	jb	leftovers		; less than 12 elements
loop:	
	movups	xmm2,[esi]		;move 4 values from first_ptr etc.
	movups	xmm3,[edi]		;move 4 values from second_ptr etc.
	movups	xmm4,[esi+16]
	movups	xmm5,[edi+16]
	movups	xmm6,[esi+32]
	movups	xmm7,[edi+32]
		
	mulps	xmm2, xmm3		;multiply
	mulps	xmm4, xmm5
	mulps	xmm6, xmm7
	
	addps	xmm0,xmm2
	addps	xmm0,xmm4
	addps	xmm0,xmm6		;sum everything in xmm0
	
	add	esi,48
	add	edi,48
	mov	eax,dword [esi]
	sub	ecx,12			;we computed 12 values
	cmp	ecx,12
	jl	sum			;less then 12 values left
	jmp	loop
sum:
	xor	eax,eax
	movups	[ebx], xmm0		;4 sums to stack
	movss	xmm0, [ebx]
	movss	xmm1, [ebx+4]
	movss	xmm2, [ebx+8]
	movss	xmm3, [ebx+12]

	addss	xmm0, xmm1
	addss	xmm2, xmm3
	addss	xmm0, xmm2		;sum of all elements in xmm0
	movss	[ebx+16],xmm0		;temp sum
	xorps	xmm0,xmm0
	cmp	ecx,0			;finish?
	je	computeKernelEnd
leftovers:
	movss	xmm1, [esi]		;first_ptr
	movss	xmm2, [edi]		;second_ptr
	mulss	xmm1,xmm2		;first_ptr*second_ptr
	addss	xmm0,xmm1		;xmm0 += first_ptr*second_ptr
	add	esi,4
	add	edi,4
	dec	ecx			;--counter
	jnz	leftovers	
computeKernelEnd:
	fld	dword [ebx+16]		;read data from sum
	movss	[ebx+16],xmm0		;read data from leftovers
	fld	dword [ebx+16]		
	fadd	ST0,ST1			;add sum+leftovers
	fstp	dword [edx+12]		;save result
	fstp	ST0
	mov	esp,[ebx+20]		;retrieve esp
	lea	esp,[esp+24]		;free memory
	pop	esi
	pop	edi
	pop	eax
	pop	edx
	pop	ecx
	pop	ebx
	mov	esp,ebp
	pop     ebp
	ret

;;
;Function performs updating cached error.
;@param pointer to ParallelAsmData structure
;;
updateErrorCache:
   	push    ebp			;save ebp
	mov     ebp,esp	
	push	ebx
	push	ecx
	push	edx
	push	eax			; save registers
	xor	eax,eax
	finit
	mov	edx,[ebp+8]		; structure
	lea	esp,[esp-68]		; 4 places for kernel High and 4 for Low
	mov	ebx,esp			; stack pointer in ebx
	mov	[ebx+64],esp		; save esp
		
	mov	eax,dword [edx+56]
	mov	dword [ebx+32],eax	; number of features
	mov	eax,dword [edx+60]
	mov	dword [ebx+36],eax	; number of examples
	mov	eax,dword [edx+4]
	mov	dword [ebx+40],eax	; iHigh
	mov	eax,dword [edx+8]
	mov	dword [ebx+44],eax	; iLow
	
	mov	esi,dword [edx+32]	; esi - error array
	mov	edi,dword [edx+52]	; edi - X array
	mov	ecx,dword [edx]		; counter - trainDataSize
	mov	dword [ebx+56],ecx
	
	mov	ecx,dword [ebx+32]	; number of features
	imul	ecx,dword [ebx+40]	; nof*iHigh
	imul	ecx,4
	add	edi,ecx			; X[iHigh]
	mov	dword [ebx+48],edi	; X[iHigh]	
	sub	edi,ecx
	mov	ecx,dword [ebx+32]	; number of features
	imul	ecx,dword [ebx+44]	; nof*iLow
	imul	ecx,4
	add	edi,ecx			; X[iLow]
	mov	dword [ebx+52],edi	; X[iLow]
	sub	edi,ecx
	mov	ecx,dword [ebx+32]	; number of features
	imul	ecx,dword [edx]		; nof*threadDataSize
	imul	ecx,4
	imul	ecx,dword [edx+12]	; threadID
	add	edi,ecx			; X[threadDataSize*threadID]
	mov	dword [ebx+60],edi	;

	mov	ecx,dword [ebx+56]	; counter - trainDataSize
	cmp	ecx,4
	jl	errorEnd
	
errorLoop:	
	mov	eax,dword [edx]
	sub	eax,dword [ebx+56]	; eax - iteration number
	sub	eax,1
	getKernel 0,48,48
	getKernel 1,48,48
	getKernel 2,48,48
	getKernel 3,48,48
	sub	eax,4			; decrement counter
	getKernel 4,44,52
	getKernel 5,44,52
	getKernel 6,44,52
	getKernel 7,44,52
	sub	eax,3			; decrement counter
errorLoopBody:
	movups	xmm2,[ebx]		; high cache
	movups	xmm4,[ebx+16]		; low cache	
	mov	eax,dword [edx+36]	; load errorUpdateHigh
	mov	dword [ebx],eax		; save errorUpdateHigh 4 times on stack
	mov	dword [ebx+4],eax		
	mov	dword [ebx+8],eax		
	mov	dword [ebx+12],eax		
	movups	xmm1,[ebx]		; xmm1 - errorUpdateHigh
	mov	eax,dword [edx+40]	; load errorUpdateLow
	mov	dword [ebx+16],eax	; save errorUpdateLow 4 times on stack
	mov	dword [ebx+20],eax		
	mov	dword [ebx+24],eax		
	mov	dword [ebx+28],eax	
	movups	xmm3,[ebx+16]		; xmm3 - errorUpdateLow
	movups	xmm0,[esi]		; xmm0 - error
	mulps	xmm1,xmm2
	mulps	xmm3,xmm4
	addps	xmm0,xmm1
	addps	xmm0,xmm3
	movups	[esi],xmm0		; save updated error
	add	esi,16
	mov	ecx,dword [ebx+56]	; counter
	sub	ecx,4
	mov	dword [ebx+56],ecx
	cmp	ecx,0
	jne	errorLoop
errorEnd:	
	mov	esp,[ebx+64]		;retrieve esp
	lea	esp,[esp+68]		;free memory
	pop	eax
	pop	edx
	pop	ecx
	pop	ebx
	pop     ebp
	ret
	
;;
;Function updates alpha values.
;@param pointer to ParallelAsmData structure
;;
updateAlpha:
   	push    ebp			;save ebp
	mov     ebp,esp	
	push	ebx
	push	ecx
	push	edx
	push	eax			; save registers
	finit
	mov	edx,[ebp+8]		; structure
	lea	esp,[esp-40]		; 7 places + esp
	mov	ebx,esp			; stack pointer in ebx
	mov	[ebx+36],esp		; save esp
	
	mov	edi,dword [edx+52]	; edi - X array
	mov	ecx,dword [edx+56]	; number of features
	imul	ecx,dword [edx+4]	; nof*iHigh
	imul	ecx,4
	add	edi,ecx			; X[iHigh]
	mov	dword [ebx+12],edi	; X[iHigh]	
	sub	edi,ecx
	mov	ecx,dword [edx+56]	; number of features
	imul	ecx,dword [edx+8]	; nof*iLow
	imul	ecx,4
	add	edi,ecx			; X[iLow]
	mov	dword [ebx+16],edi	; X[iLow]
	sub	edi,ecx
	
	mov	esi,dword [edx+28]	; alpha array
	mov	ecx,dword [edx+4]
	imul	ecx,4			; iHigh offset
	add	esi,ecx
	mov	eax,dword [esi]
	mov	dword [ebx+28],eax	; alpha(high)
	sub	esi,ecx
	
	mov	ecx,dword [edx+8]
	imul	ecx,4			; iLow offset
	add	esi,ecx
	mov	eax,dword [esi]
	mov	dword [ebx+32],eax	; alpha(low)
	sub	esi,ecx

	;compute eta
	xor	ecx,ecx
	mov	eax,dword [edx+4]	; iHigh
	getKernel 0,48,12		; kernel(high,high)
	getKernel 1,44,16		; kernel(low,low)
	getKernel 2,44,12		; kernel(low,high)
	fld	dword [ebx]
	fld	dword [ebx+4]
	fld	dword [ebx+8]
	fadd	ST0,ST0
	fchs
	fadd	ST0,ST1
	fadd	ST0,ST2
	fstp	dword [ebx+12]		; save eta
	fstp	ST0
	fstp	ST0
	mov	eax,dword [ebx+32]
	mov	dword [ebx+16],eax	; save alphaLowOld
	
	;compute alpha low boundaries
	mov	edi,dword [edx+24]	; y array	
	mov	ecx,dword [edx+8]
	imul	ecx,4			; iLow offset
	add	edi,ecx
	fld	dword [edi]		; y(ilow)
	sub	edi,ecx	
	mov	ecx,dword [edx+4]
	imul	ecx,4			; iHigh offset
	add	edi,ecx
	fld	dword [edi]		; y(iHigh)
	sub	edi,ecx
	fmul	ST0,ST1
	fldz
	fcomi	ST0,ST1
	fstp	ST0
	fstp	ST0
	fstp	ST0
	jb	signPositive		; 0 < y(high)*y(low)?
signNegative:
	fld	dword [ebx+32]
	fld	dword [ebx+28]
	fchs
	fadd	ST0,ST1			; alpha(low) - alpha(high)
	fldz
	fcomi	ST0,ST1
	fstp	ST0			; remove zero
	jbe	signNegAlphaDiffNeg	; 0 < alpha(low) - alpha(high)?
	;alpha diff positive
	fstp	dword [ebx+20]		; save lower bound
	fstp	ST0
	mov	eax, dword [edx+16]	; load cost
	mov	dword [ebx+24],eax	; save upper bound
	jmp	checkEta
signNegAlphaDiffNeg:
	fld	dword [edx+16]
	fadd	ST0,ST1			; alpha(low) - alpha(high) + cost
	fstp	dword [ebx+24]		; save upper bound
	fstp	ST0
	fstp	ST0
	mov	dword [ebx+20],0	; save lower bound
	jmp	checkEta
signPositive:	
	fld	dword [ebx+32]
	fld	dword [ebx+28]
	fadd	ST0,ST1			; alpha(low) + alpha(high)
	fld	dword [edx+16]		
	fcomi	ST0,ST1			; cost > alpha(low) + alpha(high) ?
	jbe	signPosAlphaSumGreater
	;alpha diff lower than cost
	fstp	ST0
	fstp	dword [ebx+28]		; save upper bound
	fstp	ST0
	mov	dword [ebx+32],0	; save lower bound
	jmp	checkEta
signPosAlphaSumGreater:
	fsub	ST1,ST0			; alpha(low)+alpha(high)-C
	fstp	dword [ebx+28]		; save upper bound
	fsub	dword [ebx+32]		; save lower bound
	fstp	ST0
	jmp	checkEta
	
checkEta:
alphaEnd:	
	mov	esp,[ebx+36]		;retrieve esp
	lea	esp,[esp+40]		;free memory
	pop	eax
	pop	edx
	pop	ecx
	pop	ebx
	pop     ebp
	ret