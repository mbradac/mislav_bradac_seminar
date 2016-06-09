.section .text
.global SearchNormal16

# rdi - query_array
# rsi - query_length
# rdx - score_array
# rcx - f_array
# r8 - h_array
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
# ymm8 - score

SearchNormal16:
	xorq	%r9,		%r9
	vpxor	%ymm3,		%ymm3,	%ymm3
	vpxor	%ymm4,		%ymm4,	%ymm4
	vpxor	%ymm5,		%ymm5,	%ymm5
	shlq	$2,		%rsi
	jmp	.btm
.top:	vmovdqa	32(%r8,%r9,8),	%ymm6
	vmovdqa	32(%rcx,%r9,8),	%ymm7
	xorq	%r10,		%r10 # TODO(maknut ovo)
	movb	(%rdi),		%r10b
	shlq	$5,		%r10
	vmovdqa (%rdx,%r10),	%ymm8
	vpaddsw %ymm3,		%ymm8,	%ymm3
	vpmaxsw %ymm3,		%ymm4,	%ymm3
	vpmaxsw %ymm3,		%ymm5,	%ymm3
	vpmaxsw %ymm3,		%ymm7,	%ymm3
	vmovdqa %ymm3,		32(%r8,%r9,8)
	vpmaxsw %ymm2,		%ymm3,	%ymm2
	vpsubsw %ymm1,		%ymm5,	%ymm5
	vpsubsw %ymm1,		%ymm7,	%ymm7
	vpsubsw %ymm0,		%ymm3,	%ymm3
	vpmaxsw %ymm5,		%ymm3,	%ymm5
	vpmaxsw %ymm7,		%ymm3,	%ymm7
	vmovdqa %ymm6,		%ymm3
	vmovdqa %ymm7,		32(%rcx,%r9,8)
	addq	$4,		%r9
	addq	$1,		%rdi
.btm:	cmpq	%r9,		%rsi
	jne	.top
	vmovdqa	%ymm2,		%ymm0
	ret
