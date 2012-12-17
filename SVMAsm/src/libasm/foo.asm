extern _GLOBAL_OFFSET_TABLE_

global findHighLow:function
segment .data
fHigh:	dd	10000000.0
fLow:	dd	-1000000
zeroConst: dd	0
iHigh:	dd	1
segment .text
%macro	arrayToFPU 2
	push	edx
	mov	edx,%1
	fld	dword[edx+4*%2]
	pop	edx
%endmacro

%macro inInterval 4
	push	edx
	push	ecx
	mov	edx,%1
	add	edx,%3
	cmp	%4,edx
	jg	is
	mov	edx,%2
	sub	edx,%3
	cmp	edx,%4
	jg	is
	jmp	isnt
is:	mov	eax,1
	jmp	intervalFinal
isnt	mov	eax,0
intervalFinal:	pop	edx
	pop	ecx  
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
	xor	ecx,ecx
	add	ecx,1
	xor	eax,eax
	finit
	fld	dword [edx+16]	;load cost	ST4
	fld	dword [edx+20]	;load error	ST3
;for(i = 0;i < trainDataSize)
LOOP:	
	arrayToFPU	[edx+24],0 ;load y(i)	ST2
	arrayToFPU	[edx+28],2 ;load alpha(i) ST1
	
	;checks if error < alpha < cost-error
First:	fild	dword	[zeroConst]	;load zero ST0
	fadd	ST0,ST3		;add error to low border
	fcomi	ST0,ST1		;compare low+error > val
	;fstp	dword [edx+16]
	;fstp	dword [edx+20]
	;fstp	ST0
	;fst	dword [edx+16]
	jb	is		;if ST0 < ST1 - may be in interval
	jmp	isnt
is:	fstp	ST0
	fld	ST3		;cost at stack top
	;
	fsub	ST0,ST3		;cost = cost-error
	fcomi	ST0,ST1		;compare high-error > val
	fstp	ST0
	jb	isnt		;if ST0 < ST1 - is not in interval
	jmp	highCheck
isnt	jmp	Second

Second: jmp endLoop

highCheck:
	arrayToFPU	[edx+32],2;load error(i)
	fld	dword [fHigh]
	fcomi	ST0,ST1
	jnb	highCheckLower
	fstp	ST0
	fstp	ST0
	jmp endLoop
highCheckLower:
	mov	dword [iHigh],eax ;iHigh = i
	fstp	ST0
	fstp	dword [fHigh]	;fHigh = ST0
	
endLoop:
	add	eax,1
	loop	LOOP
	;arrayToFPU	[edx+28],2
	;fst	dword [edx+20]
	fld	dword [fHigh]
	fst	dword [edx+20]
	fild	dword [iHigh]
	fst	dword [edx+16]
	;fstp	ST0
	;fst	dword [edx+16]	
end:	pop	ebx
	mov     esp,ebp
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
