.section .text
.global SearchNormal16, SearchMasked16

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
# r11 - real_query_length

.macro INIT_16
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
	movq	%rsi,		%r11
	shrq	$4,		%rsi
	shlq	$4,		%rsi
.endm

.macro STEP_16 s, h, prevh, e
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

.macro STEP4_16 s, h, prevh, e
	vpaddsw \s(%r10),	\h,	\h
	vpmaxsw \h,		\e,	\h
	vpmaxsw \h,		%ymm3,	\h
	vpmaxsw %ymm2,		\h,	%ymm2
	vmovdqa \h,		\prevh(%rcx,%r9,4)
	vpsubsw %ymm1,		\e,	\e
	vpsubsw %ymm1,		%ymm3,	%ymm3
	vpsubsw %ymm0,		\h,	\h
	vpmaxsw \e,		\h,	\e
	vpmaxsw %ymm3,		\h,	%ymm3
.endm

SearchNormal16:
	INIT_16
	jmp	.btms

.tops:	movq	(%rdi,%r9),	%r10
	vmovdqa	32(%rcx,%r9,4),	%ymm6
	vmovdqa	32(%rdx,%r9,4),	%ymm3
	STEP_16 0, %ymm4, %ymm9, %ymm5
	STEP_16 32, %ymm7, %ymm12, %ymm8
	STEP_16 64, %ymm10, %ymm15, %ymm11
	STEP4_16 96, %ymm13, 32, %ymm14
	vmovdqa %ymm3,		32(%rdx,%r9,4)

	movq	8(%rdi,%r9),	%r10
	vmovdqa	64(%rcx,%r9,4),	%ymm4
	vmovdqa	64(%rdx,%r9,4),	%ymm3
	STEP_16 0, %ymm6, %ymm7, %ymm5
	STEP_16 32, %ymm9, %ymm10, %ymm8
	STEP_16 64, %ymm12, %ymm13, %ymm11
	STEP4_16 96, %ymm15, 64, %ymm14
	vmovdqa %ymm3,		64(%rdx,%r9,4)

	addq	$16,		%r9
.btms:	cmpq	%r9,		%rsi
	jne	.tops
	cmpq	%r9,		%r11
	je	.ends

	movq	(%rdi,%r9),	%r10
	vmovdqa	32(%rcx,%r9,4),	%ymm6
	vmovdqa	32(%rdx,%r9,4),	%ymm3
	STEP_16 0, %ymm4, %ymm9, %ymm5
	STEP_16 32, %ymm7, %ymm12, %ymm8
	STEP_16 64, %ymm10, %ymm15, %ymm11
	STEP4_16 96, %ymm13, 32, %ymm14
	vmovdqa %ymm3,		32(%rdx,%r9,4)

.ends:	vmovdqa	%ymm2,		%ymm0
	ret

SearchMasked16:
	INIT_16
	vpaddsw	%ymm2,		%ymm3,	%ymm2
	vpaddsw	%ymm2,		%ymm3,	%ymm2
	vmovdqa	%ymm3,		%ymm15
	jmp	.btmms

.topms:	movq	(%rdi,%r9),	%r10
	vmovdqa	32(%rcx,%r9,4),	%ymm6
	vmovdqa	32(%rdx,%r9,4),	%ymm3
	vpaddsw	%ymm6,		%ymm15,	%ymm6
	vpaddsw	%ymm6,		%ymm15,	%ymm6
	vpaddsw	%ymm3,		%ymm15,	%ymm3
	vpaddsw	%ymm3,		%ymm15,	%ymm3
	STEP_16 0, %ymm4, %ymm9, %ymm5
	STEP_16 32, %ymm7, %ymm12, %ymm8
	vmovdqa	%ymm15,		%ymm7
	STEP_16 64, %ymm10, %ymm15, %ymm11
	STEP4_16 96, %ymm13, 32, %ymm14
	vmovdqa %ymm3,		32(%rdx,%r9,4)

	movq	8(%rdi,%r9),	%r10
	vmovdqa	64(%rcx,%r9,4),	%ymm4
	vmovdqa	64(%rdx,%r9,4),	%ymm3
	vpaddsw	%ymm4,		%ymm7,	%ymm4
	vpaddsw	%ymm4,		%ymm7,	%ymm4
	vpaddsw	%ymm3,		%ymm7,	%ymm3
	vpaddsw	%ymm3,		%ymm7,	%ymm3
	vmovdqa	%ymm7,		%ymm13
	STEP_16 0, %ymm6, %ymm7, %ymm5
	STEP_16 32, %ymm9, %ymm10, %ymm8
	vmovdqa	%ymm13,		%ymm6
	STEP_16 64, %ymm12, %ymm13, %ymm11
	STEP4_16 96, %ymm15, 64, %ymm14
	vmovdqa %ymm3,		64(%rdx,%r9,4)
	vmovdqa	%ymm6,		%ymm15

	addq	$16,		%r9
.btmms:	cmpq	%r9,		%rsi
	jne	.topms
	cmpq	%r9,		%r11
	je	.endms

	movq	(%rdi,%r9),	%r10
	vmovdqa	32(%rcx,%r9,4),	%ymm6
	vmovdqa	32(%rdx,%r9,4),	%ymm3
	vpaddsw	%ymm6,		%ymm15,	%ymm6
	vpaddsw	%ymm6,		%ymm15,	%ymm6
	vpaddsw	%ymm3,		%ymm15,	%ymm3
	vpaddsw	%ymm3,		%ymm15,	%ymm3
	STEP_16 0, %ymm4, %ymm9, %ymm5
	STEP_16 32, %ymm7, %ymm12, %ymm8
	vmovdqa	%ymm15,		%ymm7
	STEP_16 64, %ymm10, %ymm15, %ymm11
	STEP4_16 96, %ymm13, 32, %ymm14
	vmovdqa %ymm3,		32(%rdx,%r9,4)

.endms:	vmovdqa	%ymm2,		%ymm0
	ret
