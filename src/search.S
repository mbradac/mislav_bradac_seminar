.section .text
.global SearchNormal16

# rdi - query_score
# rsi - query_length
# rdx - f_array
# rcx - h_array
# ymm0 - q
# ymm1 - r
# ymm2 - results
# --
# r9 - j
# ymm3 - h
# ymm4 - z
# ymm5 - e
# ymm6 - prev_h
# ymm7 - f
# r10 - query

SearchNormal16:
	xorq	%r9,		%r9
	vpxor	%ymm3,		%ymm3,	%ymm3
	vpxor	%ymm4,		%ymm4,	%ymm4
	vpxor	%ymm5,		%ymm5,	%ymm5
	shlq	$3,		%rsi
	jmp	.btm
.top:	movq	(%rdi,%r9),	%r10
	vmovdqa	32(%rcx,%r9,4),	%ymm6
	vmovdqa	32(%rdx,%r9,4),	%ymm7
	vpaddsw (%r10),		%ymm3,	%ymm3
	vpmaxsw %ymm3,		%ymm4,	%ymm3
	vpmaxsw %ymm3,		%ymm5,	%ymm3
	vpmaxsw %ymm3,		%ymm7,	%ymm3
	vmovdqa %ymm3,		32(%rcx,%r9,4)
	vpmaxsw %ymm2,		%ymm3,	%ymm2
	vpsubsw %ymm1,		%ymm5,	%ymm5
	vpsubsw %ymm1,		%ymm7,	%ymm7
	vpsubsw %ymm0,		%ymm3,	%ymm3
	vpmaxsw %ymm5,		%ymm3,	%ymm5
	vpmaxsw %ymm7,		%ymm3,	%ymm7
	vmovdqa %ymm6,		%ymm3
	vmovdqa %ymm7,		32(%rdx,%r9,4)
	addq	$8,		%r9
.btm:	cmpq	%r9,		%rsi
	jne	.top
	vmovdqa	%ymm2,		%ymm0
	ret
