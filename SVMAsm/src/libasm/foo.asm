extern _GLOBAL_OFFSET_TABLE_


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
