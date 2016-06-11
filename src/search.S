.section .text
.global SearchNormal16

# rdi - query_score
# rsi - query_length
# rdx - f_array
# rcx - h_array
# ymm0 - q
# ymm1 - r
# ymm2 - results
# r8 - z
# --
# r9 - j
# ymm3 - f
# ymm4 - h1
# ymm5 - e1
# ymm6 - prev_h1
# ymm7 - h2
# ymm8 - e2
# ymm9 - prev_h2
# ymm10 - h3
# ymm11 - e3
# ymm12 - prev_h3
# ymm13 - h4
# ymm14 - e4
# ymm15 - prev_h4
# r10 - query

SearchNormal16:
	xorq	%r9,		%r9
	vmovdqa (%r8),		%ymm4
	vmovdqa (%r8),		%ymm5
	vmovdqa (%r8),		%ymm7
	vmovdqa (%r8),		%ymm8
	vmovdqa (%r8),		%ymm10
	vmovdqa (%r8),		%ymm11
	vmovdqa (%r8),		%ymm13
	vmovdqa (%r8),		%ymm14
	shlq	$3,		%rsi
	jmp	.btm

.top:	movq	(%rdi,%r9,4),	%r10
	vmovdqa	32(%rcx,%r9,4),	%ymm6
	vmovdqa	32(%rdx,%r9,4),	%ymm3

	vpaddsw (%r10),		%ymm4,	%ymm4
	vpmaxsw %ymm4,		%ymm5,	%ymm4
	vpmaxsw %ymm4,		%ymm3,	%ymm9
	vpmaxsw %ymm2,		%ymm9,	%ymm2
	vpsubsw %ymm1,		%ymm5,	%ymm5
	vpsubsw %ymm1,		%ymm3,	%ymm3
	vpsubsw %ymm0,		%ymm9,	%ymm4
	vpmaxsw %ymm5,		%ymm4,	%ymm5
	vpmaxsw %ymm3,		%ymm4,	%ymm3

	vpaddsw 32(%r10),	%ymm7,	%ymm7
	vpmaxsw %ymm7,		%ymm8,	%ymm7
	vpmaxsw %ymm7,		%ymm3,	%ymm12
	vpmaxsw %ymm2,		%ymm12,	%ymm2
	vpsubsw %ymm1,		%ymm8,	%ymm8
	vpsubsw %ymm1,		%ymm3,	%ymm3
	vpsubsw %ymm0,		%ymm12,	%ymm7
	vpmaxsw %ymm8,		%ymm7,	%ymm8
	vpmaxsw %ymm3,		%ymm7,	%ymm3

	vpaddsw 64(%r10),	%ymm10,	%ymm10
	vpmaxsw %ymm10,		%ymm11,	%ymm10
	vpmaxsw %ymm10,		%ymm3,	%ymm15
	vpmaxsw %ymm2,		%ymm15,	%ymm2
	vpsubsw %ymm1,		%ymm11,	%ymm11
	vpsubsw %ymm1,		%ymm3,	%ymm3
	vpsubsw %ymm0,		%ymm15,	%ymm10
	vpmaxsw %ymm11,		%ymm10,	%ymm11
	vpmaxsw %ymm3,		%ymm10,	%ymm3

	vpaddsw 96(%r10),	%ymm13,	%ymm13
	vpmaxsw %ymm13,		%ymm14,	%ymm13
	vpmaxsw %ymm13,		%ymm3,	%ymm13
	vpmaxsw %ymm2,		%ymm13,	%ymm2
	vmovdqa %ymm13,		32(%rcx,%r9,4)
	vpsubsw %ymm1,		%ymm14,	%ymm14
	vpsubsw %ymm1,		%ymm3,	%ymm3
	vpsubsw %ymm0,		%ymm13,	%ymm13
	vpmaxsw %ymm14,		%ymm13,	%ymm14
	vpmaxsw %ymm3,		%ymm13,	%ymm3

	vmovdqa %ymm3,		32(%rdx,%r9,4)

	vmovdqa	%ymm6,		%ymm4
	vmovdqa	%ymm9,		%ymm7
	vmovdqa	%ymm12,		%ymm10
	vmovdqa	%ymm15,		%ymm13

	addq	$8,		%r9
.btm:	cmpq	%r9,		%rsi
	jne	.top
	vmovdqa	%ymm2,		%ymm0
	ret
