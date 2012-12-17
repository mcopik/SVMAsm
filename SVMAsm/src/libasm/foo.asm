extern _GLOBAL_OFFSET_TABLE_

global findHighLow:function
segment .data
fHigh:	dd	10000000
fLow:	dd	-1000000
zeroConst: dd	0
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
	fld	dword [edx+16]	;load cost
	fld	dword [edx+20]	;load error	
	fld	dword [fHigh]	;fHigh init val
	fld	dword [fLow]	;fLow init val
	arrayToFPU	[edx+28],2
;for(i = 0;i < trainDataSize)
LOOP:	
	add	ebx,1
	loop	LOOP
	sub	ebx,1
	mov	[edx],ebx 
	;fld	dword [edx+16]
	;fld	dword [edx+20]
	;comi	ST0,ST1
	;jne	final
	;jmp	end
	
final:	;fadd	dword [edx+20]
	;fstp	dword [edx+16]
	;mov	[edx],ebx
	;sub	ebx,1
	;mov	ecx, dword [fHigh]	
	push	edx
	push	ecx
	fild	dword	[zeroConst]
	fadd	ST0,ST4		;add error to low border
	fcomi	ST0,ST1		;compare low+error > val
	jb	is		;if ST0 < ST1 - may be in interval
	jmp	isnt
is:	fstp	ST0
	;fst	dword [edx+16]
	fld	ST4		;cost at stack top
	fsub	ST0,ST4		;cost = cost-error
	fcomi	ST0,ST1		;compare high-error > val
	;fst	dword [edx+20]
	jb	isnt		;if ST0 < ST1 - is not in interval
	mov	eax,1
	jmp	intervalFinal
isnt	mov	eax,0
intervalFinal:	pop	ecx
	pop	edx
	;arrayToFPU	[edx+28],2
	;fst	dword [edx+20]
	mov	[edx], eax
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
