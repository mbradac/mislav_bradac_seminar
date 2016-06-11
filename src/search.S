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

.macro STEP s, h, prevh, e
	vpaddsw \s(%r10),	\h,	\h
	vpmaxsw \h,		\e,	\h
	vpmaxsw \h,		%ymm3,	\prevh
	vpmaxsw %ymm2,		\prevh,	%ymm2
	vpsubsw %ymm1,		\e,	\e
	vpsubsw %ymm1,		%ymm3,	%ymm3
	vpsubsw %ymm0,		\prevh,	\h
	vpmaxsw \e,		\h,	\e
	vpmaxsw %ymm3,		\h,	%ymm3
.endm

.macro STEP4 s, h, e
	vpaddsw \s(%r10),	\h,	\h
	vpmaxsw \h,		\e,	\h
	vpmaxsw \h,		%ymm3,	\h
	vpmaxsw %ymm2,		\h,	%ymm2
	vmovdqa \h,		32(%rcx,%r9,4)
	vpsubsw %ymm1,		\e,	\e
	vpsubsw %ymm1,		%ymm3,	%ymm3
	vpsubsw %ymm0,		\h,	\h
	vpmaxsw \e,		\h,	\e
	vpmaxsw %ymm3,		\h,	%ymm3
.endm

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
	STEP 0, %ymm4, %ymm9, %ymm5
	STEP 32, %ymm7, %ymm12, %ymm8
	STEP 64, %ymm10, %ymm15, %ymm11
	STEP4 96, %ymm13, %ymm14
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
