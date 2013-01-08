extern _GLOBAL_OFFSET_TABLE_

global findHighLow:function
segment .data
zeroConst:	dd	0
fHigh:		dd	10000000.0
fLow:		dd	-1000000.0
iHigh:		dd	-1
iLow:		dd	-1
lowFlag		dd	0
highFlag	dd	0
segment .text

%macro	arrayToFPU 2
	push	edx
	mov	edx,%1
	fld	dword[edx+4*%2]
	pop	edx
%endmacro

; ST0 - first operand
; ST1 - second operand
; ST3 - maximal error
; eax = 1 when |ST0-ST1| <= error
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

;	Input: ParallelAsmData struct
;	Output: modifications in struct
;	Finds High and Low alpha
global findHighLow
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
	
	;only testAlpha0 should pass!
	;jmp	alpha0
	;jmp	alphaCost
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
	;jmp	highCheck
	fld	ST3			;stack: cost,alphas(i),y(i),error,cost
	equalsWithTolerance ST0,ST1,ST3 ;check if |alphas(i)-cost| <= error
	;mov	eax,1
	cmp	eax,0
	fstp	ST0			;stack: alphas(i),y(i),error,cost
	je	highCheck		;alphas(i) < cost-error
	fld1				;stack: 0,alphas(i),y(i),error,cost
	fcomi	ST0,ST2
	fstp	ST0			;stack: alphas(i),y(i),error,cost
	
	;fld	ST0
	;fstp	dword [edx+4]
	je	alphaCostYPos		;y(i) > 0
	;jmp	alphaCostYPos
;if y(i) < 0 -> check iHigh
alphaCostYNeg:
	;mov	dword [edx+4],0
	;fstsw	word [edx+4]
	;fild	dword [edx+4]
	;fld	ST4
	;fstp	dword [edx+4]
	mov	eax,1
	mov	dword [ebx],eax		;check iHigh
	jmp	highCheck
;if y(i) > 0 -> check iLow
alphaCostYPos:
	;fst	dword [edx+8]
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
	jnb	highCheckAssign		;set new value
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
	fstp	dword [edx+4]		;ptr->iHigh = iHigh
	fild	dword [ebx+8]
	fstp	dword [edx+8]		;ptr->iLow = iLow
end:	
	mov	esp,[ebx+24]
	lea	esp,[esp+28]
	pop	ebx
	pop     ebp
	ret

;	Input: Pointers to data,size
;	Output: sum in eax
;	Finds High and Low alpha
global computeLinearKernel
computeLinearKernel:
   	push    ebp	;save ebp
	mov     ebp,esp	
	push	ebx
	push	ecx
	push	edx
	push	eax
	xor	eax,eax
	mov	edx,[ebp+8]	;structure
	lea	esp,[esp-24]	;four numbers+sum+temp_esp
	mov	ebx,esp		; stack pointer in ebx
	mov	[ebx+20],esp	;save esp
	mov	esi,[edx]	; first_ptr in esi
	mov	edi,[edx+4]	; second_ptr in edi
	xorps	xmm0,xmm0	; 0 in xmm0
	mov	ecx,dword [edx+8];get size
	cmp	ecx,12		;check if size >= 12
	jb	sum		;less than 12 elements
	;mov	dword [edx+12],2
	;jmp	computeKernelEnd
loop:	
	movups	xmm2,[esi]
	movups	xmm3,[edi]
	movups	xmm4,[esi+16]
	movups	xmm5,[edi+16]
	movups	xmm6,[esi+32]
	movups	xmm7,[edi+32]
		
	mulps	xmm2, xmm3
	mulps	xmm4, xmm5
	mulps	xmm6, xmm7
	
	addps	xmm0,xmm2
	addps	xmm0,xmm4
	addps	xmm0,xmm6
	sub	ecx,12
	cmp	ecx,12
	jl	sum
	add	esi,48
	add	edi,48
	jmp	loop
sum:
	xor	eax,eax
	movups	[ebx], xmm0	;4 sums to stack
	movss	xmm0, [ebx]
	movss	xmm1, [ebx+4]
	movss	xmm2, [ebx+8]
	movss	xmm3, [ebx+12]

	addss	xmm0, xmm1
	addss	xmm2, xmm3
	addss	xmm0, xmm2
	movss	[ebx+16],xmm0
	movss	[edx+12],xmm0	
	cmp	ecx,0
	xorps	xmm0,xmm0
	jmp	computeKernelEnd
leftovers:
	movss	xmm1, [esi]
	movss	xmm2, [edi]
	mulss	xmm1,xmm2
	movss	xmm0,xmm1
	sub	ecx,1
	jnz	leftovers	
computeKernelEnd:
	mov	eax,dword [ebx+16]	;read data from sum
	movss	[ebx+16],xmm0		;read data from leftovers
	add	eax,dword [ebx+16]	;add data
	mov	dword [edx+12],eax	;save result
	mov	esp,[ebx+20]
	lea	esp,[esp+24]
	pop	eax
	pop	edx
	pop	ecx
	pop	ebx
	pop     ebp
	ret

global foo:function
segment .data
str: db "Ala ma kota", 0 ;string
segment .text
;	Input: adress in memory
;	Output: value "5"
;	Writes 'A' letter to adress in memory
global foo
foo:
   	push    ebp
	mov     ebp,esp
    push    ebx
	call    .get_GOT
.get_GOT:
    pop     ebx
    add     ebx,_GLOBAL_OFFSET_TABLE_+$$-.get_GOT wrt ..gotpc

	mov		edx,[ebp+8]						;edx = program parameter 1
	xor 	eax,eax							;eax=0
	mov		al,byte [ebx+str wrt ..gotoff]	;al='A'
	mov 	[edx],eax						;*edx = 'A'
	;mov		ecx,[ebp+12]					;ecx = program parameter 2
	;mov		edx,[ebp+16]					;edx = program parameter 3
	;mov		eax,[ebp+20]
	;movups	xmm0,[ecx]
	;movups	xmm1,[edx]
	;mulps	xmm0,xmm1
	;movups	[eax],xmm0
	mov     eax,5
    mov     esp,ebp
    pop     ebp
    ret
