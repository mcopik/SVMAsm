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
; eax = 1 when |ST0-ST1| < error
%macro equalsWithTolerance 3
	fsub	%1,%2
	fabs
	fcomi	%1,%3
	jbe	@@equal
	jmp	@@notEqual
@@equal:
	mov	eax,1
	jmp	@@finish
@@notEqual:
	mov	eax,0
@@finish:
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

;checks if error < alpha < cost-error
middle:	
	fldz				;stack: 0,alphas(i),y(i),error,cost
	fadd	ST0,ST3			;add error to low border
	fcomi	ST0,ST1			;compare low+error > val
	fstp	ST0			;stack: alphas(i),y(i),error,cost
	jb	is			;if ST0 < ST1 - may be in interval
	jmp	alpha0			;second condition
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
;checks if 0 <= alpha <= error
alpha0:	
	fldz				;stack: 0,alphas(i),y(i),error,cost
	equalsWithTolerance ST0,ST1,ST3 ;check if alphas(i) <= error
	cmp	eax,0
	fstp	ST0			;stack: alphas(i),y(i),error,cost
	je	endLoop			;alphas(i) > error
	fldz				;stack: 0,alphas(i),y(i),error,cost
	fcomi	ST0,ST2
	fstp	ST0			;stack: alphas(i),y(i),error,cost
	jb	alpha0YPos		;y(i) > 0
alpha0YNeg:
	fst	dword [edx+8]
	mov	eax,1
	mov	dword [ebx+4],eax	;check iLow
	jmp	highCheck
alpha0YPos:
	mov	eax,1
	mov	dword [ebx],eax		;check iHigh
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
	fst	dword [edx+4]		;ptr->iHigh = iHigh
	fild	dword [ebx+8]
	fst	dword [edx+8]		;ptr->iLow = iLow
end:	
	mov	esp,[ebx+24]
	lea	esp,[esp+28]
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
