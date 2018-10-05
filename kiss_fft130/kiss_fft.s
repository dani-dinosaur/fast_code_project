	.file	"kiss_fft.c"
# GNU C (GCC) version 4.8.1 (x86_64-unknown-linux-gnu)
#	compiled by GNU C version 4.8.1, GMP version 5.1.2, MPFR version 3.1.2, MPC version 1.0.1
# GGC heuristics: --param ggc-min-expand=100 --param ggc-min-heapsize=131072
# options passed:  kiss_fft.c --param l1-cache-size=0
# --param l1-cache-line-size=0 --param l2-cache-size=256 -mtune=generic
# -march=x86-64 -O3 -ffast-math -fomit-frame-pointer -fverbose-asm
# options enabled:  -faggressive-loop-optimizations -fassociative-math
# -fasynchronous-unwind-tables -fauto-inc-dec -fbranch-count-reg
# -fcaller-saves -fcombine-stack-adjustments -fcommon -fcompare-elim
# -fcprop-registers -fcrossjumping -fcse-follow-jumps -fcx-limited-range
# -fdefer-pop -fdelete-null-pointer-checks -fdevirtualize -fdwarf2-cfi-asm
# -fearly-inlining -feliminate-unused-debug-types -fexpensive-optimizations
# -ffinite-math-only -fforward-propagate -ffunction-cse -fgcse
# -fgcse-after-reload -fgcse-lm -fgnu-runtime -fguess-branch-probability
# -fhoist-adjacent-loads -fident -fif-conversion -fif-conversion2
# -findirect-inlining -finline -finline-atomics -finline-functions
# -finline-functions-called-once -finline-small-functions -fipa-cp
# -fipa-cp-clone -fipa-profile -fipa-pure-const -fipa-reference -fipa-sra
# -fira-hoist-pressure -fira-share-save-slots -fira-share-spill-slots
# -fivopts -fkeep-static-consts -fleading-underscore -fmerge-constants
# -fmerge-debug-strings -fmove-loop-invariants -fomit-frame-pointer
# -foptimize-register-move -foptimize-sibling-calls -foptimize-strlen
# -fpartial-inlining -fpeephole -fpeephole2 -fpredictive-commoning
# -fprefetch-loop-arrays -freciprocal-math -free -freg-struct-return
# -fregmove -freorder-blocks -freorder-functions -frerun-cse-after-loop
# -fsched-critical-path-heuristic -fsched-dep-count-heuristic
# -fsched-group-heuristic -fsched-interblock -fsched-last-insn-heuristic
# -fsched-rank-heuristic -fsched-spec -fsched-spec-insn-heuristic
# -fsched-stalled-insns-dep -fschedule-insns2 -fshow-column -fshrink-wrap
# -fsplit-ivs-in-unroller -fsplit-wide-types -fstrict-aliasing
# -fstrict-overflow -fstrict-volatile-bitfields -fsync-libcalls
# -fthread-jumps -ftoplevel-reorder -ftree-bit-ccp -ftree-builtin-call-dce
# -ftree-ccp -ftree-ch -ftree-coalesce-vars -ftree-copy-prop
# -ftree-copyrename -ftree-cselim -ftree-dce -ftree-dominator-opts
# -ftree-dse -ftree-forwprop -ftree-fre -ftree-loop-distribute-patterns
# -ftree-loop-if-convert -ftree-loop-im -ftree-loop-ivcanon
# -ftree-loop-optimize -ftree-parallelize-loops= -ftree-partial-pre
# -ftree-phiprop -ftree-pre -ftree-pta -ftree-reassoc -ftree-scev-cprop
# -ftree-sink -ftree-slp-vectorize -ftree-slsr -ftree-sra
# -ftree-switch-conversion -ftree-tail-merge -ftree-ter
# -ftree-vect-loop-version -ftree-vectorize -ftree-vrp -funit-at-a-time
# -funsafe-math-optimizations -funswitch-loops -funwind-tables
# -fvect-cost-model -fverbose-asm -fzero-initialized-in-bss
# -m128bit-long-double -m64 -m80387 -maccumulate-outgoing-args
# -malign-stringops -mfancy-math-387 -mfp-ret-in-387 -mfxsr -mglibc
# -mlong-double-80 -mmmx -mno-sse4 -mpush-args -mred-zone -msse -msse2
# -mtls-direct-seg-refs

	.text
	.p2align 4,,15
	.type	kf_work, @function
kf_work:
.LFB38:
	.cfi_startproc
# BLOCK 2 freq:44 seq:0
# PRED: ENTRY [100.0%]  (FALLTHRU)
	pushq	%r15	#
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14	#
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13	#
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	movq	%rdi, %r13	# Fout, Fout
	pushq	%r12	#
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	movq	%rsi, %r12	# f, f
	pushq	%rbp	#
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movq	%r9, %rbp	# st, st
	pushq	%rbx	#
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$72, %rsp	#,
	.cfi_def_cfa_offset 128
	movl	4(%r8), %eax	# MEM[(int *)factors_11(D) + 4B], m
	movl	(%r8), %r14d	# *factors_11(D), p
	movq	%rdx, 32(%rsp)	# fstride, %sfp
	movl	%ecx, 12(%rsp)	# in_stride, %sfp
	movl	%eax, %edx	# m, D.5049
	movl	%eax, 48(%rsp)	# m, %sfp
	imull	%r14d, %edx	# p, D.5049
	cmpl	$1, %eax	#, m
	movslq	%edx, %rdx	# D.5049, D.5050
	leaq	(%rdi,%rdx,8), %rbx	#, Fout_end
# SUCC: 22 [28.0%]  (CAN_FALLTHRU) 3 [72.0%]  (FALLTHRU,CAN_FALLTHRU)
	je	.L40	#,
# BLOCK 3 freq:32 seq:1
# PRED: 2 [72.0%]  (FALLTHRU,CAN_FALLTHRU)
	movq	32(%rsp), %rdi	# %sfp, fstride
	movslq	12(%rsp), %rsi	# %sfp, D.5050
	movq	%r13, %r15	# Fout, Fout
	movslq	48(%rsp), %rdx	# %sfp, D.5050
	movl	%r14d, 52(%rsp)	# p, %sfp
	movq	%r13, 56(%rsp)	# Fout, %sfp
	movq	%r12, %r13	# f, f
	movq	%r15, %r12	# Fout, Fout
	movq	%rdi, %rax	# fstride, fstride
	salq	$3, %rax	#, D.5050
	imulq	%rax, %rsi	# D.5050, D.5050
	movq	%rax, 40(%rsp)	# D.5050, %sfp
	movslq	%r14d, %rax	# p, D.5050
	imulq	%rdi, %rax	# fstride, D.5050
	leaq	8(%r8), %rdi	#, factors
	movq	%rdi, %r15	# factors, factors
	movq	%rsi, 16(%rsp)	# D.5050, %sfp
	movq	%rax, 24(%rsp)	# D.5050, %sfp
	leaq	0(,%rdx,8), %rax	#, D.5050
# SUCC: 4 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	movq	%rax, %r14	# D.5050, D.5050
# BLOCK 4 freq:352 seq:2
# PRED: 3 [100.0%]  (FALLTHRU,CAN_FALLTHRU) 4 [91.0%]  (DFS_BACK,CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L5:
	movl	12(%rsp), %ecx	# %sfp,
	movq	24(%rsp), %rdx	# %sfp,
	movq	%r13, %rsi	# f,
	movq	%r12, %rdi	# Fout,
	movq	%rbp, %r9	# st,
	movq	%r15, %r8	# factors,
	addq	%r14, %r12	# D.5050, Fout
	call	kf_work	#
	addq	16(%rsp), %r13	# %sfp, f
	cmpq	%rbx, %r12	# Fout_end, Fout
# SUCC: 4 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 5 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	jne	.L5	#,
# BLOCK 5 freq:32 seq:3
# PRED: 4 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	movl	52(%rsp), %r14d	# %sfp, p
	movq	56(%rsp), %r13	# %sfp, Fout
	cmpl	$3, %r14d	#, p
# SUCC: 25 [20.0%]  (CAN_FALLTHRU) 6 [80.0%]  (FALLTHRU,CAN_FALLTHRU)
	je	.L7	#,
# BLOCK 6 freq:35 seq:4
# PRED: 5 [80.0%]  (FALLTHRU,CAN_FALLTHRU) 24 [80.0%]  (CAN_FALLTHRU)
.L43:
# SUCC: 7 [62.5%]  (FALLTHRU,CAN_FALLTHRU) 18 [37.5%]  (CAN_FALLTHRU)
	jle	.L41	#,
# BLOCK 7 freq:22 seq:5
# PRED: 6 [62.5%]  (FALLTHRU,CAN_FALLTHRU)
	cmpl	$4, %r14d	#, p
# SUCC: 13 [40.0%]  (CAN_FALLTHRU) 8 [60.0%]  (FALLTHRU,CAN_FALLTHRU)
	je	.L10	#,
# BLOCK 8 freq:13 seq:6
# PRED: 7 [60.0%]  (FALLTHRU,CAN_FALLTHRU)
	cmpl	$5, %r14d	#, p
	.p2align 4,,3
# SUCC: 9 [66.7%]  (FALLTHRU,CAN_FALLTHRU) 28 [33.3%]  (CAN_FALLTHRU)
	jne	.L6	#,
# BLOCK 9 freq:9 seq:7
# PRED: 8 [66.7%]  (FALLTHRU,CAN_FALLTHRU)
	movl	48(%rsp), %ebx	# %sfp, m
	movq	32(%rsp), %rcx	# %sfp, D.5050
	addq	$264, %rbp	#, twiddles
	movslq	%ebx, %rdx	# m, D.5050
	imulq	%rdx, %rcx	# D.5050, D.5050
	salq	$3, %rdx	#, D.5050
	salq	$3, %rcx	#, D.5050
	leaq	0(%rbp,%rcx), %rax	#, D.5056
	movss	(%rax), %xmm3	# MEM[(struct kiss_fft_cpx *)_231], ya$r
	movss	%xmm3, 12(%rsp)	# ya$r, %sfp
	movss	4(%rax), %xmm3	# MEM[(struct kiss_fft_cpx *)_231 + 4B], ya$i
	addq	%rcx, %rax	# D.5050, D.5056
	leaq	0(%r13,%rdx), %rcx	#, Fout1
	movss	%xmm3, 16(%rsp)	# ya$i, %sfp
	leaq	(%rcx,%rdx), %rsi	#, Fout2
	movss	(%rax), %xmm3	# MEM[(struct kiss_fft_cpx *)_235], yb$r
	movss	%xmm3, 24(%rsp)	# yb$r, %sfp
	leaq	(%rsi,%rdx), %rdi	#, Fout3
	movss	4(%rax), %xmm4	# MEM[(struct kiss_fft_cpx *)_235 + 4B], yb$i
	addq	%rdi, %rdx	# Fout3, Fout4
	testl	%ebx, %ebx	# m
# SUCC: 10 [91.0%]  (FALLTHRU,CAN_FALLTHRU) 12 [9.0%]  (CAN_FALLTHRU)
	jle	.L1	#,
# BLOCK 10 freq:8 seq:8
# PRED: 9 [91.0%]  (FALLTHRU,CAN_FALLTHRU)
	movq	40(%rsp), %r10	# %sfp, D.5050
	movl	%ebx, %eax	# m, D.5058
	movq	%rbp, %r8	# twiddles, ivtmp.221
	subl	$1, %eax	#, D.5058
	movss	16(%rsp), %xmm15	# %sfp, D.5055
	leaq	8(%rdx,%rax,8), %r9	#, D.5056
	xorl	%eax, %eax	# ivtmp.208
# SUCC: 11 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	leaq	(%r10,%r10,2), %r11	#, tmp522
# BLOCK 11 freq:89 seq:9
# PRED: 11 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 10 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
.L20:
	movss	0(%rbp,%rax), %xmm2	# MEM[base: twiddles_227, index: ivtmp.208_181, offset: 0B], D.5055
	addq	$8, %rdx	#, Fout4
	movss	4(%rbp,%rax), %xmm1	# MEM[base: twiddles_227, index: ivtmp.208_181, offset: 4B], D.5055
	addq	$8, %r13	#, Fout
	movaps	%xmm2, %xmm12	# D.5055, D.5055
	addq	$8, %rcx	#, Fout1
	movaps	%xmm1, %xmm3	# D.5055, D.5055
	addq	$8, %rsi	#, Fout2
	movss	-8(%rcx), %xmm0	# MEM[base: Fout1_512, offset: 0B], D.5055
	addq	$8, %rdi	#, Fout3
	movss	-4(%rcx), %xmm13	# MEM[base: Fout1_512, offset: 4B], D.5055
	mulss	%xmm0, %xmm12	# D.5055, D.5055
	mulss	%xmm13, %xmm3	# D.5055, D.5055
	movss	-4(%rsi), %xmm14	# MEM[base: Fout2_513, offset: 4B], D.5055
	mulss	%xmm1, %xmm0	# D.5055, D.5055
	movss	4(%rbp,%rax,2), %xmm1	# MEM[base: twiddles_227, index: ivtmp.208_181, step: 2, offset: 4B], D.5055
	mulss	%xmm2, %xmm13	# D.5055, D.5055
	movss	0(%rbp,%rax,2), %xmm2	# MEM[base: twiddles_227, index: ivtmp.208_181, step: 2, offset: 0B], D.5055
	movaps	%xmm2, %xmm11	# D.5055, D.5055
	subss	%xmm3, %xmm12	# D.5055, D.5055
	movaps	%xmm1, %xmm3	# D.5055, D.5055
	movss	(%r8), %xmm5	# MEM[base: _1024, offset: 0B], D.5055
	addss	%xmm0, %xmm13	# D.5055, D.5055
	movss	-8(%rsi), %xmm0	# MEM[base: Fout2_513, offset: 0B], D.5055
	mulss	%xmm14, %xmm3	# D.5055, D.5055
	movss	-4(%rdx), %xmm8	# MEM[base: Fout4_509, offset: 4B], D.5055
	mulss	%xmm0, %xmm11	# D.5055, D.5055
	movss	-8(%r13), %xmm9	# MEM[base: Fout_514, offset: 0B], scratch$0$r
	mulss	%xmm1, %xmm0	# D.5055, D.5055
	movaps	%xmm5, %xmm1	# D.5055, D.5055
	mulss	%xmm2, %xmm14	# D.5055, D.5055
	movss	-8(%rdi), %xmm2	# MEM[base: Fout3_510, offset: 0B], D.5055
	mulss	%xmm2, %xmm1	# D.5055, D.5055
	movss	-4(%r13), %xmm6	# MEM[base: Fout_514, offset: 4B], scratch$0$i
	subss	%xmm3, %xmm11	# D.5055, D.5055
	movss	4(%r8), %xmm3	# MEM[base: _1024, offset: 4B], D.5055
	movaps	%xmm3, %xmm7	# D.5055, D.5055
	mulss	%xmm3, %xmm2	# D.5055, D.5055
	addss	%xmm0, %xmm14	# D.5055, D.5055
	movss	-4(%rdi), %xmm0	# MEM[base: Fout3_510, offset: 4B], D.5055
	mulss	%xmm0, %xmm7	# D.5055, D.5055
	addq	%r11, %r8	# tmp522, ivtmp.221
	mulss	%xmm5, %xmm0	# D.5055, D.5055
	movss	-8(%rdx), %xmm5	# MEM[base: Fout4_509, offset: 0B], D.5055
	subss	%xmm7, %xmm1	# D.5055, D.5055
	movss	4(%rbp,%rax,4), %xmm7	# MEM[base: twiddles_227, index: ivtmp.208_181, step: 4, offset: 4B], D.5055
	addss	%xmm2, %xmm0	# D.5055, D.5055
	movss	0(%rbp,%rax,4), %xmm2	# MEM[base: twiddles_227, index: ivtmp.208_181, step: 4, offset: 0B], D.5055
	movaps	%xmm2, %xmm3	# D.5055, D.5055
	mulss	%xmm8, %xmm2	# D.5055, D.5055
	movaps	%xmm7, %xmm10	# D.5055, D.5055
	addq	%r10, %rax	# D.5050, ivtmp.208
	mulss	%xmm5, %xmm3	# D.5055, D.5055
	mulss	%xmm8, %xmm10	# D.5055, D.5055
	movaps	%xmm1, %xmm8	# D.5055, D.5055
	mulss	%xmm7, %xmm5	# D.5055, D.5055
	addss	%xmm11, %xmm8	# D.5055, D.5055
	subss	%xmm1, %xmm11	# D.5055, D.5055
	subss	%xmm10, %xmm3	# D.5055, D.5055
	addss	%xmm5, %xmm2	# D.5055, D.5055
	movaps	%xmm0, %xmm5	# D.5055, D.5055
	addss	%xmm14, %xmm5	# D.5055, D.5055
	movaps	%xmm3, %xmm10	# D.5055, D.5055
	subss	%xmm0, %xmm14	# D.5055, D.5055
	movaps	%xmm2, %xmm7	# D.5055, D.5055
	addss	%xmm12, %xmm10	# D.5055, D.5055
	addss	%xmm13, %xmm7	# D.5055, D.5055
	subss	%xmm3, %xmm12	# D.5055, D.5055
	movss	12(%rsp), %xmm3	# %sfp, ya$r
	movaps	%xmm3, %xmm1	# ya$r, D.5055
	movaps	%xmm10, %xmm0	# D.5055, D.5055
	mulss	%xmm10, %xmm1	# D.5055, D.5055
	subss	%xmm2, %xmm13	# D.5055, D.5055
	addss	%xmm9, %xmm0	# scratch$0$r, D.5055
	movss	24(%rsp), %xmm2	# %sfp, yb$r
	addss	%xmm9, %xmm1	# scratch$0$r, D.5055
	mulss	%xmm13, %xmm15	# D.5055, D.5055
	addss	%xmm8, %xmm0	# D.5055, D.5055
	mulss	%xmm4, %xmm13	# yb$i, D.5055
	movss	%xmm0, -8(%r13)	# D.5055, MEM[base: Fout_514, offset: 0B]
	movaps	%xmm7, %xmm0	# D.5055, D.5055
	addss	%xmm6, %xmm0	# scratch$0$i, D.5055
	addss	%xmm5, %xmm0	# D.5055, D.5055
	movss	%xmm0, -4(%r13)	# D.5055, MEM[base: Fout_514, offset: 4B]
	movaps	%xmm2, %xmm0	# yb$r, D.5055
	mulss	%xmm5, %xmm2	# D.5055, D.5055
	mulss	%xmm8, %xmm0	# D.5055, D.5055
	addss	%xmm0, %xmm1	# D.5055, D.5055
	movaps	%xmm3, %xmm0	# ya$r, D.5055
	movaps	%xmm14, %xmm3	# D.5055, D.5055
	mulss	%xmm7, %xmm0	# D.5055, D.5055
	mulss	%xmm4, %xmm3	# yb$i, D.5055
	addss	%xmm6, %xmm0	# scratch$0$i, D.5055
	addss	%xmm15, %xmm3	# D.5055, D.5055
	movss	16(%rsp), %xmm15	# %sfp, D.5055
	mulss	%xmm12, %xmm15	# D.5055, D.5055
	mulss	%xmm4, %xmm12	# yb$i, D.5055
	addss	%xmm2, %xmm0	# D.5055, D.5055
	movaps	%xmm11, %xmm2	# D.5055, D.5055
	mulss	%xmm4, %xmm2	# yb$i, D.5055
	addss	%xmm15, %xmm2	# D.5055, D.5055
	movaps	%xmm1, %xmm15	# D.5055, D.5055
	addss	%xmm3, %xmm1	# D.5055, D.5055
	subss	%xmm3, %xmm15	# D.5055, D.5055
	movss	%xmm15, -8(%rcx)	# D.5055, MEM[base: Fout1_512, offset: 0B]
	movaps	%xmm0, %xmm15	# D.5055, D.5055
	subss	%xmm2, %xmm0	# D.5055, D.5055
	addss	%xmm2, %xmm15	# D.5055, D.5055
	movss	24(%rsp), %xmm2	# %sfp, yb$r
	mulss	%xmm2, %xmm10	# yb$r, D.5055
	mulss	%xmm2, %xmm7	# yb$r, D.5055
	movss	%xmm15, -4(%rcx)	# D.5055, MEM[base: Fout1_512, offset: 4B]
	movss	%xmm0, -4(%rdx)	# D.5055, MEM[base: Fout4_509, offset: 4B]
	addss	%xmm9, %xmm10	# scratch$0$r, D.5055
	movss	12(%rsp), %xmm0	# %sfp, ya$r
	addss	%xmm6, %xmm7	# scratch$0$i, D.5055
	mulss	%xmm0, %xmm8	# ya$r, D.5055
	movss	16(%rsp), %xmm15	# %sfp, ya$i
	mulss	%xmm15, %xmm14	# ya$i, D.5055
	movss	%xmm1, -8(%rdx)	# D.5055, MEM[base: Fout4_509, offset: 0B]
	mulss	%xmm0, %xmm5	# ya$r, D.5055
	mulss	%xmm15, %xmm11	# ya$i, D.5055
	addss	%xmm8, %xmm10	# D.5055, D.5055
	subss	%xmm13, %xmm14	# D.5055, D.5055
	addss	%xmm5, %xmm7	# D.5055, D.5055
	subss	%xmm11, %xmm12	# D.5055, D.5055
	movaps	%xmm10, %xmm0	# D.5055, D.5055
	subss	%xmm14, %xmm10	# D.5055, D.5055
	addss	%xmm14, %xmm0	# D.5055, D.5055
	movss	%xmm0, -8(%rsi)	# D.5055, MEM[base: Fout2_513, offset: 0B]
	movaps	%xmm7, %xmm0	# D.5055, D.5055
	subss	%xmm12, %xmm7	# D.5055, D.5055
	addss	%xmm12, %xmm0	# D.5055, D.5055
	movss	%xmm0, -4(%rsi)	# D.5055, MEM[base: Fout2_513, offset: 4B]
	movss	%xmm10, -8(%rdi)	# D.5055, MEM[base: Fout3_510, offset: 0B]
	movss	%xmm7, -4(%rdi)	# D.5055, MEM[base: Fout3_510, offset: 4B]
	cmpq	%r9, %rdx	# D.5056, Fout4
# SUCC: 11 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 12 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	jne	.L20	#,
# BLOCK 12 freq:17 seq:10
# PRED: 15 [9.0%]  (CAN_FALLTHRU,LOOP_EXIT) 9 [9.0%]  (CAN_FALLTHRU) 11 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
.L1:
	addq	$72, %rsp	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx	#
	.cfi_def_cfa_offset 48
	popq	%rbp	#
	.cfi_def_cfa_offset 40
	popq	%r12	#
	.cfi_def_cfa_offset 32
	popq	%r13	#
	.cfi_def_cfa_offset 24
	popq	%r14	#
	.cfi_def_cfa_offset 16
	popq	%r15	#
	.cfi_def_cfa_offset 8
# SUCC: EXIT [100.0%] 
	ret
# BLOCK 13 freq:9 seq:11
# PRED: 7 [40.0%]  (CAN_FALLTHRU)
.L10:
	.cfi_restore_state
	movslq	48(%rsp), %r9	# %sfp, k
	movq	32(%rsp), %rax	# %sfp, fstride
	leaq	264(%rbp), %rsi	#, tw3
	movl	4(%rbp), %r11d	# st_33(D)->inverse, D.5049
	movq	40(%rsp), %rbp	# %sfp, D.5050
	movq	%rsi, %r8	# tw3, tw3
	movq	%rsi, %rdi	# tw3, tw3
	movq	%rax, %r10	# fstride, D.5050
	leaq	(%rax,%rax,2), %rbx	#, tmp491
	leaq	0(,%r9,8), %rcx	#, D.5050
	salq	$4, %r10	#, D.5050
	salq	$3, %rbx	#, tmp492
	leaq	0(%r13,%rcx), %rdx	#, ivtmp.175
	leaq	(%rdx,%rcx), %rax	#, ivtmp.180
	addq	%rax, %rcx	# ivtmp.180, ivtmp.185
# SUCC: 16 [100.0%] 
	jmp	.L19	#
# BLOCK 14 freq:49 seq:12
# PRED: 16 [50.0%]  (CAN_FALLTHRU)
.L42:
	movaps	%xmm9, %xmm0	# D.5055, D.5055
	subss	%xmm1, %xmm0	# D.5055, D.5055
	addss	%xmm9, %xmm1	# D.5055, D.5055
	movss	%xmm0, (%rdx)	# D.5055, MEM[base: _1048, offset: 0B]
	movaps	%xmm3, %xmm0	# D.5055, D.5055
	addss	%xmm8, %xmm0	# D.5055, D.5055
	subss	%xmm3, %xmm8	# D.5055, D.5055
	movss	%xmm0, 4(%rdx)	# D.5055, MEM[base: _1048, offset: 4B]
	movss	%xmm1, (%rcx)	# D.5055, MEM[base: _1044, offset: 0B]
# SUCC: 15 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	movss	%xmm8, 4(%rcx)	# D.5055, MEM[base: _1044, offset: 4B]
# BLOCK 15 freq:98 seq:13
# PRED: 14 [100.0%]  (FALLTHRU,CAN_FALLTHRU) 17 [100.0%] 
.L18:
	addq	$8, %r13	#, Fout
	addq	$8, %rdx	#, ivtmp.175
	addq	$8, %rax	#, ivtmp.180
	addq	$8, %rcx	#, ivtmp.185
	subq	$1, %r9	#, k
# SUCC: 16 [91.0%]  (FALLTHRU,DFS_BACK,CAN_FALLTHRU) 12 [9.0%]  (CAN_FALLTHRU,LOOP_EXIT)
	je	.L1	#,
# BLOCK 16 freq:98 seq:14
# PRED: 15 [91.0%]  (FALLTHRU,DFS_BACK,CAN_FALLTHRU) 13 [100.0%] 
.L19:
	movss	(%rdx), %xmm2	# MEM[base: _1048, offset: 0B], D.5055
	movss	4(%rdx), %xmm4	# MEM[base: _1048, offset: 4B], D.5055
	movaps	%xmm2, %xmm3	# D.5055, D.5055
	movaps	%xmm4, %xmm5	# D.5055, D.5055
	movss	(%rsi), %xmm1	# MEM[base: tw3_157, offset: 0B], D.5055
	movss	4(%rsi), %xmm0	# MEM[base: tw3_157, offset: 4B], D.5055
	mulss	%xmm1, %xmm3	# D.5055, D.5055
	addq	%rbp, %rsi	# D.5050, tw3
	mulss	%xmm0, %xmm5	# D.5055, D.5055
	movss	4(%rax), %xmm6	# MEM[base: _1046, offset: 4B], D.5055
	mulss	%xmm2, %xmm0	# D.5055, D.5055
	movss	(%rax), %xmm2	# MEM[base: _1046, offset: 0B], D.5055
	mulss	%xmm4, %xmm1	# D.5055, D.5055
	movaps	%xmm6, %xmm7	# D.5055, D.5055
	movss	(%rdi), %xmm4	# MEM[base: tw3_170, offset: 0B], D.5055
	subss	%xmm5, %xmm3	# D.5055, D.5055
	movaps	%xmm2, %xmm5	# D.5055, D.5055
	movss	4(%rcx), %xmm8	# MEM[base: _1044, offset: 4B], D.5055
	addss	%xmm0, %xmm1	# D.5055, D.5055
	movss	4(%rdi), %xmm0	# MEM[base: tw3_170, offset: 4B], D.5055
	mulss	%xmm0, %xmm7	# D.5055, D.5055
	movaps	%xmm8, %xmm9	# D.5055, D.5055
	mulss	%xmm4, %xmm5	# D.5055, D.5055
	addq	%r10, %rdi	# D.5050, tw3
	mulss	%xmm2, %xmm0	# D.5055, D.5055
	mulss	%xmm6, %xmm4	# D.5055, D.5055
	movss	4(%r8), %xmm6	# MEM[base: tw3_183, offset: 4B], D.5055
	mulss	%xmm6, %xmm9	# D.5055, D.5055
	subss	%xmm7, %xmm5	# D.5055, D.5055
	movss	(%rcx), %xmm7	# MEM[base: _1044, offset: 0B], D.5055
	movaps	%xmm7, %xmm2	# D.5055, D.5055
	mulss	%xmm7, %xmm6	# D.5055, D.5055
	addss	%xmm0, %xmm4	# D.5055, D.5055
	movss	(%r8), %xmm0	# MEM[base: tw3_183, offset: 0B], D.5055
	mulss	%xmm0, %xmm2	# D.5055, D.5055
	movss	0(%r13), %xmm7	# MEM[base: Fout_154, offset: 0B], D.5055
	mulss	%xmm8, %xmm0	# D.5055, D.5055
	addq	%rbx, %r8	# tmp492, tw3
	testl	%r11d, %r11d	# D.5049
	subss	%xmm9, %xmm2	# D.5055, D.5055
	movaps	%xmm7, %xmm9	# D.5055, D.5055
	addss	%xmm6, %xmm0	# D.5055, D.5055
	movss	4(%r13), %xmm6	# MEM[base: Fout_154, offset: 4B], D.5055
	movaps	%xmm6, %xmm8	# D.5055, D.5055
	subss	%xmm5, %xmm9	# D.5055, D.5055
	subss	%xmm4, %xmm8	# D.5055, D.5055
	addss	%xmm7, %xmm5	# D.5055, D.5055
	movaps	%xmm2, %xmm7	# D.5055, D.5055
	addss	%xmm6, %xmm4	# D.5055, D.5055
	movaps	%xmm0, %xmm6	# D.5055, D.5055
	addss	%xmm3, %xmm7	# D.5055, D.5055
	addss	%xmm1, %xmm6	# D.5055, D.5055
	movss	%xmm5, 0(%r13)	# D.5055, MEM[base: Fout_154, offset: 0B]
	subss	%xmm2, %xmm3	# D.5055, D.5055
	movss	%xmm4, 4(%r13)	# D.5055, MEM[base: Fout_154, offset: 4B]
	subss	%xmm0, %xmm1	# D.5055, D.5055
	subss	%xmm7, %xmm5	# D.5055, D.5055
	subss	%xmm6, %xmm4	# D.5055, D.5055
	movss	%xmm5, (%rax)	# D.5055, MEM[base: _1046, offset: 0B]
	movss	%xmm4, 4(%rax)	# D.5055, MEM[base: _1046, offset: 4B]
	addss	0(%r13), %xmm7	# MEM[base: Fout_154, offset: 0B], D.5055
	addss	4(%r13), %xmm6	# MEM[base: Fout_154, offset: 4B], D.5055
	movss	%xmm7, 0(%r13)	# D.5055, MEM[base: Fout_154, offset: 0B]
	movss	%xmm6, 4(%r13)	# D.5055, MEM[base: Fout_154, offset: 4B]
# SUCC: 14 [50.0%]  (CAN_FALLTHRU) 17 [50.0%]  (FALLTHRU,CAN_FALLTHRU)
	jne	.L42	#,
# BLOCK 17 freq:49 seq:15
# PRED: 16 [50.0%]  (FALLTHRU,CAN_FALLTHRU)
	movaps	%xmm1, %xmm0	# D.5055, D.5055
	addss	%xmm9, %xmm0	# D.5055, D.5055
	subss	%xmm1, %xmm9	# D.5055, D.5055
	movss	%xmm0, (%rdx)	# D.5055, MEM[base: _1048, offset: 0B]
	movaps	%xmm8, %xmm0	# D.5055, D.5055
	subss	%xmm3, %xmm0	# D.5055, D.5055
	addss	%xmm8, %xmm3	# D.5055, D.5055
	movss	%xmm0, 4(%rdx)	# D.5055, MEM[base: _1048, offset: 4B]
	movss	%xmm9, (%rcx)	# D.5055, MEM[base: _1044, offset: 0B]
	movss	%xmm3, 4(%rcx)	# D.5055, MEM[base: _1044, offset: 4B]
# SUCC: 15 [100.0%] 
	jmp	.L18	#
# BLOCK 18 freq:13 seq:16
# PRED: 6 [37.5%]  (CAN_FALLTHRU)
.L41:
	cmpl	$2, %r14d	#, p
# SUCC: 19 [66.7%]  (FALLTHRU,CAN_FALLTHRU) 28 [33.3%]  (CAN_FALLTHRU)
	jne	.L6	#,
# BLOCK 19 freq:9 seq:17
# PRED: 18 [66.7%]  (FALLTHRU,CAN_FALLTHRU)
	movl	48(%rsp), %edx	# %sfp, m
	movq	40(%rsp), %rcx	# %sfp, D.5050
	addq	$264, %rbp	#, tw1
	movslq	%edx, %rax	# m, D.5050
# SUCC: 20 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	leaq	0(%r13,%rax,8), %rax	#, Fout2
# BLOCK 20 freq:98 seq:18
# PRED: 19 [100.0%]  (FALLTHRU,CAN_FALLTHRU) 20 [91.0%]  (DFS_BACK,CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L14:
	movss	0(%rbp), %xmm0	# MEM[base: tw1_56, offset: 0B], D.5055
	addq	$8, %rax	#, Fout2
	movss	4(%rbp), %xmm3	# MEM[base: tw1_56, offset: 4B], D.5055
	addq	$8, %r13	#, Fout
	movaps	%xmm0, %xmm1	# D.5055, D.5055
	addq	%rcx, %rbp	# D.5050, tw1
	movaps	%xmm3, %xmm5	# D.5055, D.5055
	movss	-8(%rax), %xmm2	# MEM[base: Fout2_54, offset: 0B], D.5055
	movss	-4(%rax), %xmm4	# MEM[base: Fout2_54, offset: 4B], D.5055
	mulss	%xmm2, %xmm1	# D.5055, D.5055
	mulss	%xmm4, %xmm5	# D.5055, D.5055
	mulss	%xmm3, %xmm2	# D.5055, D.5055
	mulss	%xmm4, %xmm0	# D.5055, D.5055
	subss	%xmm5, %xmm1	# D.5055, D.5055
	addss	%xmm2, %xmm0	# D.5055, D.5055
	movss	-8(%r13), %xmm2	# MEM[base: Fout_68, offset: 0B], D.5055
	subss	%xmm1, %xmm2	# D.5055, D.5055
	movss	%xmm2, -8(%rax)	# D.5055, MEM[base: Fout2_54, offset: 0B]
	movss	-4(%r13), %xmm2	# MEM[base: Fout_68, offset: 4B], D.5055
	subss	%xmm0, %xmm2	# D.5055, D.5055
	movss	%xmm2, -4(%rax)	# D.5055, MEM[base: Fout2_54, offset: 4B]
	addss	-8(%r13), %xmm1	# MEM[base: Fout_68, offset: 0B], D.5055
	addss	-4(%r13), %xmm0	# MEM[base: Fout_68, offset: 4B], D.5055
	movss	%xmm1, -8(%r13)	# D.5055, MEM[base: Fout_68, offset: 0B]
	movss	%xmm0, -4(%r13)	# D.5055, MEM[base: Fout_68, offset: 4B]
	subl	$1, %edx	#, m
# SUCC: 20 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 21 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	jne	.L14	#,
# BLOCK 21 freq:9 seq:19
# PRED: 20 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	addq	$72, %rsp	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx	#
	.cfi_def_cfa_offset 48
	popq	%rbp	#
	.cfi_def_cfa_offset 40
	popq	%r12	#
	.cfi_def_cfa_offset 32
	popq	%r13	#
	.cfi_def_cfa_offset 24
	popq	%r14	#
	.cfi_def_cfa_offset 16
	popq	%r15	#
	.cfi_def_cfa_offset 8
# SUCC: EXIT [100.0%] 
	ret
# BLOCK 22 freq:12 seq:20
# PRED: 2 [28.0%]  (CAN_FALLTHRU)
.L40:
	.cfi_restore_state
	movq	32(%rsp), %rax	# %sfp, fstride
	movslq	%ecx, %r15	#, D.5050
	movq	%rdi, %rdx	# Fout, Fout
	salq	$3, %rax	#, D.5050
	imulq	%rax, %r15	# D.5050, D.5050
# SUCC: 23 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	movq	%rax, 40(%rsp)	# D.5050, %sfp
# BLOCK 23 freq:137 seq:21
# PRED: 22 [100.0%]  (FALLTHRU,CAN_FALLTHRU) 23 [91.0%]  (DFS_BACK,CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L4:
	movq	(%r12), %rax	# MEM[base: f_3, offset: 0B], MEM[base: f_3, offset: 0B]
	addq	$8, %rdx	#, Fout
	addq	%r15, %r12	# D.5050, f
	movq	%rax, -8(%rdx)	# MEM[base: f_3, offset: 0B], MEM[base: Fout_1, offset: 0B]
	cmpq	%rbx, %rdx	# Fout_end, Fout
# SUCC: 23 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 24 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	jne	.L4	#,
# BLOCK 24 freq:12 seq:22
# PRED: 23 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	cmpl	$3, %r14d	#, p
# SUCC: 25 [20.0%]  (FALLTHRU,CAN_FALLTHRU) 6 [80.0%]  (CAN_FALLTHRU)
	jne	.L43	#,
# BLOCK 25 freq:9 seq:23
# PRED: 5 [20.0%]  (CAN_FALLTHRU) 24 [20.0%]  (FALLTHRU,CAN_FALLTHRU)
.L7:
	movq	32(%rsp), %rdi	# %sfp, fstride
	movslq	48(%rsp), %rsi	# %sfp, m
	movq	40(%rsp), %r8	# %sfp, D.5050
	movsd	.LC0(%rip), %xmm8	#, tmp595
	movq	%rdi, %rax	# fstride, D.5050
	imulq	%rsi, %rax	# m, D.5050
	leaq	0(,%rsi,8), %rdx	#, D.5050
	movss	268(%rbp,%rax,8), %xmm9	# MEM[(struct kiss_fft_state *)_83 + 4B], epi3$i
	movq	%rdi, %rax	# fstride, D.5050
	salq	$4, %rax	#, D.5050
	addq	$264, %rbp	#, tw2
	movq	%rax, %rdi	# D.5050, D.5050
	leaq	0(%r13,%rdx), %rax	#, ivtmp.146
	movq	%rbp, %rcx	# tw2, tw2
# SUCC: 26 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	addq	%rax, %rdx	# ivtmp.146, ivtmp.151
# BLOCK 26 freq:98 seq:24
# PRED: 25 [100.0%]  (FALLTHRU,CAN_FALLTHRU) 26 [91.0%]  (DFS_BACK,CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L16:
	movss	(%rax), %xmm3	# MEM[base: _388, offset: 0B], D.5055
	addq	$8, %r13	#, Fout
	movss	4(%rax), %xmm4	# MEM[base: _388, offset: 4B], D.5055
	addq	$8, %rdx	#, ivtmp.151
	movaps	%xmm3, %xmm0	# D.5055, D.5055
	addq	$8, %rax	#, ivtmp.146
	movaps	%xmm4, %xmm5	# D.5055, D.5055
	movss	0(%rbp), %xmm1	# MEM[base: tw2_90, offset: 0B], D.5055
	movss	4(%rbp), %xmm2	# MEM[base: tw2_90, offset: 4B], D.5055
	mulss	%xmm1, %xmm0	# D.5055, D.5055
	addq	%r8, %rbp	# D.5050, tw2
	mulss	%xmm2, %xmm5	# D.5055, D.5055
	movss	-4(%rdx), %xmm6	# MEM[base: _386, offset: 4B], D.5055
	mulss	%xmm3, %xmm2	# D.5055, D.5055
	movaps	%xmm6, %xmm7	# D.5055, D.5055
	mulss	%xmm4, %xmm1	# D.5055, D.5055
	movss	4(%rcx), %xmm4	# MEM[base: tw2_103, offset: 4B], D.5055
	mulss	%xmm4, %xmm7	# D.5055, D.5055
	subss	%xmm5, %xmm0	# D.5055, D.5055
	movss	-8(%rdx), %xmm5	# MEM[base: _386, offset: 0B], D.5055
	movaps	%xmm5, %xmm3	# D.5055, D.5055
	mulss	%xmm5, %xmm4	# D.5055, D.5055
	addss	%xmm2, %xmm1	# D.5055, D.5055
	movss	(%rcx), %xmm2	# MEM[base: tw2_103, offset: 0B], D.5055
	mulss	%xmm2, %xmm3	# D.5055, D.5055
	addq	%rdi, %rcx	# D.5050, tw2
	mulss	%xmm6, %xmm2	# D.5055, D.5055
	subss	%xmm7, %xmm3	# D.5055, D.5055
	addss	%xmm4, %xmm2	# D.5055, D.5055
	movss	-8(%r13), %xmm4	# MEM[base: Fout_87, offset: 0B], D.5057
	cvtps2pd	%xmm4, %xmm4	# D.5057, D.5057
	movaps	%xmm3, %xmm7	# D.5055, D.5055
	movaps	%xmm2, %xmm6	# D.5055, D.5055
	addss	%xmm0, %xmm7	# D.5055, D.5055
	addss	%xmm1, %xmm6	# D.5055, D.5055
	subss	%xmm2, %xmm1	# D.5055, D.5055
	subss	%xmm3, %xmm0	# D.5055, D.5055
	unpcklps	%xmm7, %xmm7	# D.5055, D.5055
	unpcklps	%xmm6, %xmm6	# D.5055, D.5055
	mulss	%xmm9, %xmm1	# epi3$i, D.5055
	cvtps2pd	%xmm7, %xmm5	# D.5055, D.5057
	mulss	%xmm9, %xmm0	# epi3$i, D.5055
	mulsd	%xmm8, %xmm5	# tmp595, D.5057
	subsd	%xmm5, %xmm4	# D.5057, D.5057
	unpcklpd	%xmm4, %xmm4	# D.5057
	cvtpd2ps	%xmm4, %xmm5	# D.5057, tmp716
	movss	%xmm5, -8(%rax)	# tmp716, MEM[base: _388, offset: 0B]
	cvtps2pd	%xmm6, %xmm5	# D.5055, D.5057
	movss	-4(%r13), %xmm4	# MEM[base: Fout_87, offset: 4B], D.5057
	mulsd	%xmm8, %xmm5	# tmp595, D.5057
	cvtps2pd	%xmm4, %xmm4	# D.5057, D.5057
	subsd	%xmm5, %xmm4	# D.5057, D.5057
	unpcklpd	%xmm4, %xmm4	# D.5057
	cvtpd2ps	%xmm4, %xmm5	# D.5057, tmp717
	movss	%xmm5, -4(%rax)	# tmp717, MEM[base: _388, offset: 4B]
	addss	-8(%r13), %xmm7	# MEM[base: Fout_87, offset: 0B], D.5055
	addss	-4(%r13), %xmm6	# MEM[base: Fout_87, offset: 4B], D.5055
	movss	%xmm7, -8(%r13)	# D.5055, MEM[base: Fout_87, offset: 0B]
	movss	%xmm6, -4(%r13)	# D.5055, MEM[base: Fout_87, offset: 4B]
	movss	-8(%rax), %xmm2	# MEM[base: _388, offset: 0B], D.5055
	addss	%xmm1, %xmm2	# D.5055, D.5055
	movss	%xmm2, -8(%rdx)	# D.5055, MEM[base: _386, offset: 0B]
	movss	-4(%rax), %xmm2	# MEM[base: _388, offset: 4B], D.5055
	subss	%xmm0, %xmm2	# D.5055, D.5055
	movss	%xmm2, -4(%rdx)	# D.5055, MEM[base: _386, offset: 4B]
	addss	-4(%rax), %xmm0	# MEM[base: _388, offset: 4B], D.5055
	movss	-8(%rax), %xmm2	# MEM[base: _388, offset: 0B], D.5055
	subss	%xmm1, %xmm2	# D.5055, D.5055
	movss	%xmm0, -4(%rax)	# D.5055, MEM[base: _388, offset: 4B]
	movss	%xmm2, -8(%rax)	# D.5055, MEM[base: _388, offset: 0B]
	subq	$1, %rsi	#, m
# SUCC: 26 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 27 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	jne	.L16	#,
# BLOCK 27 freq:9 seq:25
# PRED: 26 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	addq	$72, %rsp	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx	#
	.cfi_def_cfa_offset 48
	popq	%rbp	#
	.cfi_def_cfa_offset 40
	popq	%r12	#
	.cfi_def_cfa_offset 32
	popq	%r13	#
	.cfi_def_cfa_offset 24
	popq	%r14	#
	.cfi_def_cfa_offset 16
	popq	%r15	#
	.cfi_def_cfa_offset 8
# SUCC: EXIT [100.0%] 
	ret
# BLOCK 28 freq:9 seq:26
# PRED: 18 [33.3%]  (CAN_FALLTHRU) 8 [33.3%]  (CAN_FALLTHRU)
.L6:
	.cfi_restore_state
	movslq	%r14d, %rdi	# p, D.5050
	movl	0(%rbp), %ebx	# st_33(D)->nfft, Norig
	leaq	264(%rbp), %r12	#, twiddles
	salq	$3, %rdi	#, D.5050
	call	malloc	#
	movl	48(%rsp), %esi	# %sfp, m
	testl	%esi, %esi	# m
# SUCC: 29 [91.0%]  (FALLTHRU,CAN_FALLTHRU) 39 [9.0%]  (CAN_FALLTHRU)
	jle	.L21	#,
# BLOCK 29 freq:8 seq:27
# PRED: 28 [91.0%]  (FALLTHRU,CAN_FALLTHRU)
	movq	32(%rsp), %rdx	# %sfp, fstride
	movslq	%esi, %rdi	# m, D.5050
	movq	%r13, %r15	# Fout, ivtmp.110
	salq	$3, %rdi	#, D.5050
	movl	$0, 16(%rsp)	#, %sfp
	movl	$0, 12(%rsp)	#, %sfp
	movl	%edx, %ebp	# fstride, D.5054
	movl	%edx, 24(%rsp)	# tmp767, %sfp
	leal	-2(%r14), %edx	#, D.5053
	imull	%esi, %ebp	# m, D.5054
# SUCC: 30 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	leaq	8(%rax,%rdx,8), %r9	#, D.5053
# BLOCK 30 freq:89 seq:28
# PRED: 38 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 29 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
.L22:
	testl	%r14d, %r14d	# p
# SUCC: 31 [91.0%]  (FALLTHRU,CAN_FALLTHRU) 38 [9.0%]  (CAN_FALLTHRU)
	jle	.L26	#,
# BLOCK 31 freq:81 seq:29
# PRED: 30 [91.0%]  (FALLTHRU,CAN_FALLTHRU)
	movq	%r15, %rcx	# ivtmp.110, ivtmp.102
# SUCC: 32 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	xorl	%edx, %edx	# ivtmp.97
# BLOCK 32 freq:900 seq:30
# PRED: 32 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 31 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L30:
	movq	(%rcx), %rsi	# MEM[base: _499, offset: 0B], MEM[base: _499, offset: 0B]
	addq	%rdi, %rcx	# D.5050, ivtmp.102
	movq	%rsi, (%rax,%rdx,8)	# MEM[base: _499, offset: 0B], MEM[base: scratch_378, index: ivtmp.97_500, step: 8, offset: 0B]
	addq	$1, %rdx	#, ivtmp.97
	cmpl	%edx, %r14d	# ivtmp.97, p
# SUCC: 32 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 33 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	jg	.L30	#,
# BLOCK 33 freq:81 seq:31
# PRED: 32 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	movl	16(%rsp), %r8d	# %sfp, ivtmp.92
	movq	(%rax), %r11	# *scratch_378, *scratch_378
	movq	%r15, %rsi	# ivtmp.110, ivtmp.94
# SUCC: 34 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	xorl	%r10d, %r10d	# q1
# BLOCK 34 freq:900 seq:32
# PRED: 37 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 33 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L29:
	cmpl	$1, %r14d	#, p
	movq	%r11, (%rsi)	# *scratch_378, MEM[base: _291, offset: 0B]
# SUCC: 35 [91.0%]  (FALLTHRU,CAN_FALLTHRU) 37 [9.0%]  (CAN_FALLTHRU)
	je	.L27	#,
# BLOCK 35 freq:819 seq:33
# PRED: 34 [91.0%]  (FALLTHRU,CAN_FALLTHRU)
	movss	(%rsi), %xmm3	# MEM[base: _291, offset: 0B], D.5055
	movq	%rax, %rdx	# scratch, ivtmp.81
	movss	4(%rsi), %xmm2	# MEM[base: _291, offset: 4B], D.5055
# SUCC: 36 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	xorl	%ecx, %ecx	# twidx
# BLOCK 36 freq:9100 seq:34
# PRED: 36 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 35 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L28:
	addl	%r8d, %ecx	# ivtmp.92, D.5054
	movss	8(%rdx), %xmm5	# MEM[base: _25, offset: 8B], D.5055
	movl	%ecx, %r13d	# D.5054, twidx
	movss	12(%rdx), %xmm6	# MEM[base: _25, offset: 12B], D.5055
	subl	%ebx, %r13d	# Norig, twidx
	cmpl	%ecx, %ebx	# D.5054, Norig
	cmovle	%r13d, %ecx	# twidx,, twidx
	movaps	%xmm5, %xmm1	# D.5055, D.5055
	movslq	%ecx, %r13	# twidx, D.5050
	movaps	%xmm6, %xmm7	# D.5055, D.5055
	leaq	(%r12,%r13,8), %r13	#, D.5056
	addq	$8, %rdx	#, ivtmp.81
	cmpq	%r9, %rdx	# D.5053, ivtmp.81
	movss	0(%r13), %xmm0	# _409->r, D.5055
	movss	4(%r13), %xmm4	# _409->i, D.5055
	mulss	%xmm0, %xmm1	# D.5055, D.5055
	mulss	%xmm4, %xmm7	# D.5055, D.5055
	mulss	%xmm6, %xmm0	# D.5055, D.5055
	mulss	%xmm5, %xmm4	# D.5055, D.5055
	subss	%xmm7, %xmm1	# D.5055, D.5055
	addss	%xmm4, %xmm0	# D.5055, D.5055
	addss	%xmm1, %xmm3	# D.5055, D.5055
	addss	%xmm0, %xmm2	# D.5055, D.5055
	movss	%xmm3, (%rsi)	# D.5055, MEM[base: _291, offset: 0B]
	movss	%xmm2, 4(%rsi)	# D.5055, MEM[base: _291, offset: 4B]
# SUCC: 36 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 37 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	jne	.L28	#,
# BLOCK 37 freq:900 seq:35
# PRED: 34 [9.0%]  (CAN_FALLTHRU) 36 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
.L27:
	addl	$1, %r10d	#, q1
	addl	%ebp, %r8d	# D.5054, ivtmp.92
	addq	%rdi, %rsi	# D.5050, ivtmp.94
	cmpl	%r14d, %r10d	# p, q1
# SUCC: 34 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 38 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	jne	.L29	#,
# BLOCK 38 freq:89 seq:36
# PRED: 30 [9.0%]  (CAN_FALLTHRU) 37 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
.L26:
	addl	$1, 12(%rsp)	#, %sfp
	movl	24(%rsp), %esi	# %sfp, D.5054
	addq	$8, %r15	#, ivtmp.110
	addl	%esi, 16(%rsp)	# D.5054, %sfp
	movl	48(%rsp), %esi	# %sfp, m
	cmpl	%esi, 12(%rsp)	# m, %sfp
# SUCC: 30 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 39 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	jne	.L22	#,
# BLOCK 39 freq:9 seq:37
# PRED: 28 [9.0%]  (CAN_FALLTHRU) 38 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
.L21:
	addq	$72, %rsp	#,
	.cfi_def_cfa_offset 56
	movq	%rax, %rdi	# scratch,
	popq	%rbx	#
	.cfi_def_cfa_offset 48
	popq	%rbp	#
	.cfi_def_cfa_offset 40
	popq	%r12	#
	.cfi_def_cfa_offset 32
	popq	%r13	#
	.cfi_def_cfa_offset 24
	popq	%r14	#
	.cfi_def_cfa_offset 16
	popq	%r15	#
	.cfi_def_cfa_offset 8
# SUCC: EXIT [100.0%]  (ABNORMAL,SIBCALL)
	jmp	free	#
	.cfi_endproc
.LFE38:
	.size	kf_work, .-kf_work
	.p2align 4,,15
	.globl	kiss_fft_alloc
	.type	kiss_fft_alloc, @function
kiss_fft_alloc:
.LFB40:
	.cfi_startproc
# BLOCK 2 freq:501 seq:0
# PRED: ENTRY [100.0%]  (FALLTHRU)
	pushq	%r15	#
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	leal	-1(%rdi), %eax	#, D.5089
	pushq	%r14	#
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	cltq
	pushq	%r13	#
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12	#
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp	#
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx	#
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	movl	%edi, %ebx	# nfft, nfft
	leaq	272(,%rax,8), %rdi	#, memneeded
	subq	$40, %rsp	#,
	.cfi_def_cfa_offset 96
	testq	%rcx, %rcx	# lenmem
# SUCC: 25 [6.7%]  (CAN_FALLTHRU) 3 [93.3%]  (FALLTHRU,CAN_FALLTHRU)
	je	.L71	#,
# BLOCK 3 freq:468 seq:1
# PRED: 2 [93.3%]  (FALLTHRU,CAN_FALLTHRU)
	testq	%rdx, %rdx	# mem
	movq	%rdx, %rbp	# mem, mem
# SUCC: 4 [85.0%]  (FALLTHRU,CAN_FALLTHRU) 5 [15.0%]  (CAN_FALLTHRU)
	je	.L47	#,
# BLOCK 4 freq:397 seq:2
# PRED: 3 [85.0%]  (FALLTHRU,CAN_FALLTHRU)
	cmpq	(%rcx), %rdi	# *lenmem_13(D), memneeded
	movl	$0, %eax	#, tmp143
# SUCC: 5 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	cmova	%rax, %rbp	# mem,, tmp143, mem
# BLOCK 5 freq:468 seq:3
# PRED: 3 [15.0%]  (CAN_FALLTHRU) 4 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
.L47:
# SUCC: 6 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	movq	%rdi, (%rcx)	# memneeded, *lenmem_13(D)
# BLOCK 6 freq:501 seq:4
# PRED: 5 [100.0%]  (FALLTHRU,CAN_FALLTHRU) 25 [100.0%] 
.L46:
	testq	%rbp, %rbp	# mem
# SUCC: 7 [89.9%]  (FALLTHRU,CAN_FALLTHRU) 23 [10.1%]  (CAN_FALLTHRU)
	je	.L68	#,
# BLOCK 7 freq:451 seq:5
# PRED: 6 [89.9%]  (FALLTHRU,CAN_FALLTHRU)
	cvtsi2sd	%ebx, %xmm3	# nfft, D.5091
	testl	%ebx, %ebx	# nfft
	movl	%ebx, 0(%rbp)	# nfft, MEM[(struct kiss_fft_state *)mem_2].nfft
	movl	%esi, 4(%rbp)	# inverse_fft, MEM[(struct kiss_fft_state *)mem_2].inverse
	movsd	%xmm3, 8(%rsp)	# D.5091, %sfp
# SUCC: 8 [91.0%]  (FALLTHRU,CAN_FALLTHRU) 10 [9.0%]  (CAN_FALLTHRU)
	jle	.L53	#,
# BLOCK 8 freq:410 seq:6
# PRED: 7 [91.0%]  (FALLTHRU,CAN_FALLTHRU)
	movsd	.LC5(%rip), %xmm0	#, D.5091
	xorl	%r12d, %r12d	# i
	movsd	.LC6(%rip), %xmm2	#, D.5091
	testl	%esi, %esi	# inverse_fft
	divsd	%xmm3, %xmm0	# D.5091, D.5091
	movq	%rbp, %r13	# mem, ivtmp.255
	leaq	24(%rsp), %r15	#, tmp140
	leaq	16(%rsp), %r14	#, tmp139
	mulsd	%xmm0, %xmm2	# D.5091, D.5091
	movsd	%xmm2, (%rsp)	# D.5091, %sfp
# SUCC: 9 [50.0%]  (FALLTHRU,CAN_FALLTHRU) 19 [50.0%]  (CAN_FALLTHRU)
	jne	.L55	#,
# BLOCK 9 freq:2278 seq:7
# PRED: 8 [50.0%]  (FALLTHRU,CAN_FALLTHRU) 9 [91.0%]  (DFS_BACK,CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L54:
	cvtsi2sd	%r12d, %xmm0	# i, phase
	movq	%r14, %rsi	# tmp139,
	movq	%r15, %rdi	# tmp140,
	addl	$1, %r12d	#, i
	addq	$8, %r13	#, ivtmp.255
	mulsd	(%rsp), %xmm0	# %sfp, phase
	call	sincos	#
	movsd	24(%rsp), %xmm0	#, D.5088
	movsd	16(%rsp), %xmm4	#, tmp151
	unpcklpd	%xmm0, %xmm0	# D.5088
	cvtpd2ps	%xmm4, %xmm4	# tmp151, tmp151
	cvtpd2ps	%xmm0, %xmm5	# D.5088, tmp152
	movss	%xmm4, 256(%r13)	# tmp151, MEM[base: _62, offset: 264B]
	movss	%xmm5, 260(%r13)	# tmp152, MEM[base: _62, offset: 268B]
	cmpl	%ebx, %r12d	# nfft, i
# SUCC: 9 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 10 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	jne	.L54	#,
# BLOCK 10 freq:451 seq:8
# PRED: 7 [9.0%]  (CAN_FALLTHRU) 9 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT) 20 [100.0%] 
.L53:
	movsd	.LC4(%rip), %xmm0	#, tmp114
	sqrtsd	8(%rsp), %xmm1	# %sfp, D.5091
	leaq	8(%rbp), %rsi	#, facbuf
	movsd	.LC3(%rip), %xmm2	#, tmp112
	andpd	%xmm1, %xmm0	# D.5091, tmp114
	ucomisd	%xmm0, %xmm2	# tmp114, tmp112
# SUCC: 11 [50.0%]  (FALLTHRU,CAN_FALLTHRU) 24 [50.0%]  (CAN_FALLTHRU)
	ja	.L72	#,
# BLOCK 11 freq:452 seq:9
# PRED: 10 [50.0%]  (FALLTHRU,CAN_FALLTHRU) 24 [100.0%] 
.L50:
	movsd	.LC2(%rip), %xmm3	#, D.5091
	movl	$4, %ecx	#, nfft
# SUCC: 12 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	movsd	.LC1(%rip), %xmm2	#, D.5091
# BLOCK 12 freq:5450 seq:10
# PRED: 11 [100.0%]  (FALLTHRU,CAN_FALLTHRU) 18 [100.0%] 
	.p2align 4,,10
	.p2align 3
.L51:
	movl	%ebx, %eax	# nfft, tmp135
	cltd
	idivl	%ecx	# nfft
	testl	%edx, %edx	# D.5089
# SUCC: 15 [50.0%]  (CAN_FALLTHRU) 13 [50.0%]  (FALLTHRU,CAN_FALLTHRU)
	jne	.L60	#,
# BLOCK 13 freq:5000 seq:11
# PRED: 12 [50.0%]  (FALLTHRU,CAN_FALLTHRU) 14 [50.0%]  (CAN_FALLTHRU)
.L73:
	movl	%eax, %ebx	# tmp135, nfft
	movl	%ecx, (%rsi)	# nfft, *facbuf_48
	leaq	8(%rsi), %rax	#, facbuf
	cmpl	$1, %ebx	#, nfft
	movl	%ebx, 4(%rsi)	# nfft, MEM[(int *)facbuf_48 + 4B]
# SUCC: 14 [91.0%]  (FALLTHRU,CAN_FALLTHRU) 23 [9.0%]  (CAN_FALLTHRU,LOOP_EXIT)
	jle	.L68	#,
# BLOCK 14 freq:4550 seq:12
# PRED: 13 [91.0%]  (FALLTHRU,CAN_FALLTHRU)
	movq	%rax, %rsi	# facbuf, facbuf
	movl	%ebx, %eax	# nfft, tmp135
	cltd
	idivl	%ecx	# nfft
	testl	%edx, %edx	# D.5089
# SUCC: 15 [50.0%]  (FALLTHRU,CAN_FALLTHRU) 13 [50.0%]  (CAN_FALLTHRU)
	je	.L73	#,
# BLOCK 15 freq:5000 seq:13
# PRED: 12 [50.0%]  (CAN_FALLTHRU) 14 [50.0%]  (FALLTHRU,CAN_FALLTHRU)
.L60:
	cmpl	$2, %ecx	#, nfft
# SUCC: 21 [33.3%]  (CAN_FALLTHRU) 16 [66.7%]  (FALLTHRU,CAN_FALLTHRU)
	je	.L62	#,
# BLOCK 16 freq:3334 seq:14
# PRED: 15 [66.7%]  (FALLTHRU,CAN_FALLTHRU)
	cmpl	$4, %ecx	#, nfft
# SUCC: 17 [50.0%]  (FALLTHRU,CAN_FALLTHRU) 22 [50.0%]  (CAN_FALLTHRU)
	jne	.L69	#,
# BLOCK 17 freq:1667 seq:15
# PRED: 16 [50.0%]  (FALLTHRU,CAN_FALLTHRU)
	movapd	%xmm2, %xmm0	# D.5091, D.5091
# SUCC: 18 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	movb	$2, %cl	#,
# BLOCK 18 freq:4999 seq:16
# PRED: 17 [100.0%]  (FALLTHRU,CAN_FALLTHRU) 21 [100.0%]  22 [100.0%] 
.L57:
	comisd	%xmm1, %xmm0	# D.5091, D.5091
	cmova	%ebx, %ecx	# nfft,, nfft, nfft
# SUCC: 12 [100.0%] 
	jmp	.L51	#
# BLOCK 19 freq:2278 seq:17
# PRED: 19 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 8 [50.0%]  (CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L55:
	cvtsi2sd	%r12d, %xmm0	# i, phase
	movq	%r14, %rsi	# tmp139,
	movq	%r15, %rdi	# tmp140,
	addl	$1, %r12d	#, i
	addq	$8, %r13	#, ivtmp.262
	mulsd	(%rsp), %xmm0	# %sfp, phase
	xorpd	.LC7(%rip), %xmm0	#, phase
	call	sincos	#
	movsd	24(%rsp), %xmm0	#, D.5088
	movsd	16(%rsp), %xmm6	#, tmp154
	unpcklpd	%xmm0, %xmm0	# D.5088
	cvtpd2ps	%xmm6, %xmm6	# tmp154, tmp154
	cvtpd2ps	%xmm0, %xmm7	# D.5088, tmp155
	movss	%xmm6, 256(%r13)	# tmp154, MEM[base: _57, offset: 264B]
	movss	%xmm7, 260(%r13)	# tmp155, MEM[base: _57, offset: 268B]
	cmpl	%ebx, %r12d	# nfft, i
# SUCC: 19 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 20 [9.0%]  (FALLTHRU)
	jne	.L55	#,
# BLOCK 20 freq:205 seq:18
# PRED: 19 [9.0%]  (FALLTHRU)
# SUCC: 10 [100.0%] 
	jmp	.L53	#
# BLOCK 21 freq:1667 seq:19
# PRED: 15 [33.3%]  (CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L62:
	movapd	%xmm3, %xmm0	# D.5091, D.5091
	movl	$3, %ecx	#, nfft
# SUCC: 18 [100.0%] 
	jmp	.L57	#
# BLOCK 22 freq:1666 seq:20
# PRED: 16 [50.0%]  (CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L69:
	addl	$2, %ecx	#, nfft
	cvtsi2sd	%ecx, %xmm0	# nfft, D.5091
# SUCC: 18 [100.0%] 
	jmp	.L57	#
# BLOCK 23 freq:1002 seq:21
# PRED: 13 [9.0%]  (CAN_FALLTHRU,LOOP_EXIT) 6 [10.1%]  (CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L68:
	addq	$40, %rsp	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	movq	%rbp, %rax	# mem,
	popq	%rbx	#
	.cfi_def_cfa_offset 48
	popq	%rbp	#
	.cfi_def_cfa_offset 40
	popq	%r12	#
	.cfi_def_cfa_offset 32
	popq	%r13	#
	.cfi_def_cfa_offset 24
	popq	%r14	#
	.cfi_def_cfa_offset 16
	popq	%r15	#
	.cfi_def_cfa_offset 8
# SUCC: EXIT [100.0%] 
	ret
# BLOCK 24 freq:226 seq:22
# PRED: 10 [50.0%]  (CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L72:
	.cfi_restore_state
	cvttsd2siq	%xmm1, %rax	# D.5091, tmp117
	cvtsi2sdq	%rax, %xmm1	# tmp117, D.5091
# SUCC: 11 [100.0%] 
	jmp	.L50	#
# BLOCK 25 freq:34 seq:23
# PRED: 2 [6.7%]  (CAN_FALLTHRU)
.L71:
	movl	%esi, (%rsp)	# inverse_fft, %sfp
	call	malloc	#
	movl	(%rsp), %esi	# %sfp, inverse_fft
	movq	%rax, %rbp	#, mem
# SUCC: 6 [100.0%] 
	jmp	.L46	#
	.cfi_endproc
.LFE40:
	.size	kiss_fft_alloc, .-kiss_fft_alloc
	.p2align 4,,15
	.globl	kiss_fft_stride
	.type	kiss_fft_stride, @function
kiss_fft_stride:
.LFB41:
	.cfi_startproc
# BLOCK 2 freq:10000 seq:0
# PRED: ENTRY [100.0%]  (FALLTHRU)
	pushq	%r12	#
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	pushq	%rbp	#
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	movq	%rsi, %rbp	# fin, fin
	pushq	%rbx	#
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	movq	%rdi, %rbx	# st, st
	movq	%rdx, %rdi	# fout, fout
	subq	$16, %rsp	#,
	.cfi_def_cfa_offset 48
	cmpq	%rdx, %rsi	# fout, fin
# SUCC: 4 [10.1%]  (CAN_FALLTHRU) 3 [89.9%]  (FALLTHRU,CAN_FALLTHRU)
	je	.L77	#,
# BLOCK 3 freq:8986 seq:1
# PRED: 2 [89.9%]  (FALLTHRU,CAN_FALLTHRU)
	addq	$16, %rsp	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	movq	%rbx, %r9	# st,
	leaq	8(%rbx), %r8	#,
	popq	%rbx	#
	.cfi_def_cfa_offset 24
	popq	%rbp	#
	.cfi_def_cfa_offset 16
	popq	%r12	#
	.cfi_def_cfa_offset 8
	movl	$1, %edx	#,
# SUCC: EXIT [100.0%]  (ABNORMAL,SIBCALL)
	jmp	kf_work	#
# BLOCK 4 freq:1014 seq:2
# PRED: 2 [10.1%]  (CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L77:
	.cfi_restore_state
	movslq	(%rbx), %rdi	# st_5(D)->nfft, D.5115
	movl	%ecx, 12(%rsp)	# in_stride, %sfp
	salq	$3, %rdi	#, D.5115
	call	malloc	#
	movl	12(%rsp), %ecx	# %sfp, in_stride
	leaq	8(%rbx), %r8	#,
	movq	%rbx, %r9	# st,
	movq	%rbp, %rsi	# fin,
	movq	%rax, %rdi	# tmp75,
	movl	$1, %edx	#,
	movq	%rax, %r12	#, tmp75
	call	kf_work	#
	movslq	(%rbx), %rdx	# st_5(D)->nfft, D.5115
	movq	%rbp, %rdi	# fin,
	movq	%r12, %rsi	# tmp75,
	salq	$3, %rdx	#, D.5115
	call	memcpy	#
	addq	$16, %rsp	#,
	.cfi_def_cfa_offset 32
	movq	%r12, %rdi	# tmp75,
	popq	%rbx	#
	.cfi_def_cfa_offset 24
	popq	%rbp	#
	.cfi_def_cfa_offset 16
	popq	%r12	#
	.cfi_def_cfa_offset 8
# SUCC: EXIT [100.0%]  (ABNORMAL,SIBCALL)
	jmp	free	#
	.cfi_endproc
.LFE41:
	.size	kiss_fft_stride, .-kiss_fft_stride
	.p2align 4,,15
	.globl	kiss_fft
	.type	kiss_fft, @function
kiss_fft:
.LFB42:
	.cfi_startproc
# BLOCK 2 freq:10000 seq:0
# PRED: ENTRY [100.0%]  (FALLTHRU)
	pushq	%r12	#
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	cmpq	%rdx, %rsi	# fout, fin
	pushq	%rbp	#
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	movq	%rsi, %rbp	# fin, fin
	pushq	%rbx	#
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	movq	%rdi, %rbx	# cfg, cfg
	movq	%rdx, %rdi	# fout, fout
# SUCC: 4 [10.1%]  (CAN_FALLTHRU) 3 [89.9%]  (FALLTHRU,CAN_FALLTHRU)
	je	.L81	#,
# BLOCK 3 freq:8986 seq:1
# PRED: 2 [89.9%]  (FALLTHRU,CAN_FALLTHRU)
	movq	%rbx, %r9	# cfg,
	leaq	8(%rbx), %r8	#,
	movl	$1, %ecx	#,
	popq	%rbx	#
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbp	#
	.cfi_def_cfa_offset 16
	popq	%r12	#
	.cfi_def_cfa_offset 8
	movl	$1, %edx	#,
# SUCC: EXIT [100.0%]  (ABNORMAL,SIBCALL)
	jmp	kf_work	#
# BLOCK 4 freq:1014 seq:2
# PRED: 2 [10.1%]  (CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L81:
	.cfi_restore_state
	movslq	(%rbx), %rdi	# cfg_2(D)->nfft, D.5131
	salq	$3, %rdi	#, D.5131
	call	malloc	#
	leaq	8(%rbx), %r8	#,
	movq	%rbx, %r9	# cfg,
	movq	%rbp, %rsi	# fin,
	movq	%rax, %rdi	# tmp74,
	movl	$1, %ecx	#,
	movl	$1, %edx	#,
	movq	%rax, %r12	#, tmp74
	call	kf_work	#
	movslq	(%rbx), %rdx	# cfg_2(D)->nfft, D.5131
	movq	%rbp, %rdi	# fin,
	movq	%r12, %rsi	# tmp74,
	salq	$3, %rdx	#, D.5131
	call	memcpy	#
	popq	%rbx	#
	.cfi_def_cfa_offset 24
	popq	%rbp	#
	.cfi_def_cfa_offset 16
	movq	%r12, %rdi	# tmp74,
	popq	%r12	#
	.cfi_def_cfa_offset 8
# SUCC: EXIT [100.0%]  (ABNORMAL,SIBCALL)
	jmp	free	#
	.cfi_endproc
.LFE42:
	.size	kiss_fft, .-kiss_fft
	.p2align 4,,15
	.globl	kiss_fft_cleanup
	.type	kiss_fft_cleanup, @function
kiss_fft_cleanup:
.LFB43:
	.cfi_startproc
# BLOCK 2 freq:10000 seq:0
# PRED: ENTRY [100.0%]  (FALLTHRU)
# SUCC: EXIT [100.0%] 
	rep; ret
	.cfi_endproc
.LFE43:
	.size	kiss_fft_cleanup, .-kiss_fft_cleanup
	.p2align 4,,15
	.globl	kiss_fft_next_fast_size
	.type	kiss_fft_next_fast_size, @function
kiss_fft_next_fast_size:
.LFB44:
	.cfi_startproc
# BLOCK 2 freq:81 seq:0
# PRED: ENTRY [100.0%]  (FALLTHRU)
	movl	$1431655766, %r8d	#, tmp126
# SUCC: 3 [100.0%]  (FALLTHRU,CAN_FALLTHRU)
	movl	$1717986919, %esi	#, tmp127
# BLOCK 3 freq:900 seq:1
# PRED: 2 [100.0%]  (FALLTHRU,CAN_FALLTHRU) 10 [100.0%] 
	.p2align 4,,10
	.p2align 3
.L91:
	testb	$1, %dil	#, n
	movl	%edi, %ecx	# n, n
# SUCC: 4 [91.0%]  (FALLTHRU,CAN_FALLTHRU) 5 [9.0%]  (CAN_FALLTHRU)
	jne	.L84	#,
# BLOCK 4 freq:9100 seq:2
# PRED: 4 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 3 [91.0%]  (FALLTHRU,CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L85:
	movl	%ecx, %eax	# n, tmp76
	shrl	$31, %eax	#, tmp76
	addl	%eax, %ecx	# tmp76, n
	sarl	%ecx	# n
	testb	$1, %cl	#, n
# SUCC: 4 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 5 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	je	.L85	#,
# BLOCK 5 freq:900 seq:3
# PRED: 3 [9.0%]  (CAN_FALLTHRU) 4 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
.L84:
	movl	%ecx, %eax	# n, tmp128
	imull	%r8d	# tmp126
	movl	%ecx, %eax	# n, tmp83
	sarl	$31, %eax	#, tmp83
	subl	%eax, %edx	# tmp83, tmp80
	leal	(%rdx,%rdx,2), %eax	#, tmp86
	cmpl	%eax, %ecx	# tmp86, n
# SUCC: 6 [91.0%]  (FALLTHRU,CAN_FALLTHRU) 7 [9.0%]  (CAN_FALLTHRU)
	jne	.L86	#,
# BLOCK 6 freq:9100 seq:4
# PRED: 6 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 5 [91.0%]  (FALLTHRU,CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L87:
	movl	%ecx, %eax	# n, tmp129
	sarl	$31, %ecx	#, tmp91
	imull	%r8d	# tmp126
	subl	%ecx, %edx	# tmp91, n
	movl	%edx, %eax	# n, tmp130
	movl	%edx, %ecx	# n, n
	imull	%r8d	# tmp126
	movl	%ecx, %r9d	# n, tmp95
	sarl	$31, %r9d	#, tmp95
	subl	%r9d, %edx	# tmp95, tmp92
	leal	(%rdx,%rdx,2), %edx	#, tmp98
	cmpl	%edx, %ecx	# tmp98, n
# SUCC: 6 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 7 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	je	.L87	#,
# BLOCK 7 freq:900 seq:5
# PRED: 6 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT) 5 [9.0%]  (CAN_FALLTHRU)
.L86:
	movl	%ecx, %eax	# n, tmp131
	imull	%esi	# tmp127
	movl	%ecx, %eax	# n, tmp104
	sarl	$31, %eax	#, tmp104
	sarl	%edx	# tmp100
	subl	%eax, %edx	# tmp104, tmp100
	leal	(%rdx,%rdx,4), %eax	#, tmp107
	cmpl	%eax, %ecx	# tmp107, n
# SUCC: 8 [91.0%]  (FALLTHRU,CAN_FALLTHRU) 9 [9.0%]  (CAN_FALLTHRU)
	jne	.L88	#,
# BLOCK 8 freq:9100 seq:6
# PRED: 8 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 7 [91.0%]  (FALLTHRU,CAN_FALLTHRU)
	.p2align 4,,10
	.p2align 3
.L89:
	movl	%ecx, %eax	# n, tmp132
	sarl	$31, %ecx	#, tmp113
	imull	%esi	# tmp127
	sarl	%edx	# tmp112
	subl	%ecx, %edx	# tmp113, n
	movl	%edx, %eax	# n, tmp133
	movl	%edx, %ecx	# n, n
	imull	%esi	# tmp127
	movl	%ecx, %r9d	# n, tmp118
	sarl	$31, %r9d	#, tmp118
	sarl	%edx	# tmp114
	subl	%r9d, %edx	# tmp118, tmp114
	leal	(%rdx,%rdx,4), %edx	#, tmp121
	cmpl	%edx, %ecx	# tmp121, n
# SUCC: 8 [91.0%]  (DFS_BACK,CAN_FALLTHRU) 9 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT)
	je	.L89	#,
# BLOCK 9 freq:900 seq:7
# PRED: 8 [9.0%]  (FALLTHRU,CAN_FALLTHRU,LOOP_EXIT) 7 [9.0%]  (CAN_FALLTHRU)
.L88:
	cmpl	$1, %ecx	#, n
# SUCC: 11 [9.0%]  (CAN_FALLTHRU,LOOP_EXIT) 10 [91.0%]  (FALLTHRU,CAN_FALLTHRU)
	jle	.L90	#,
# BLOCK 10 freq:819 seq:8
# PRED: 9 [91.0%]  (FALLTHRU,CAN_FALLTHRU)
	addl	$1, %edi	#, n
# SUCC: 3 [100.0%] 
	jmp	.L91	#
# BLOCK 11 freq:81 seq:9
# PRED: 9 [9.0%]  (CAN_FALLTHRU,LOOP_EXIT)
.L90:
	movl	%edi, %eax	# n,
	.p2align 4,,3
# SUCC: EXIT [100.0%] 
	ret
	.cfi_endproc
.LFE44:
	.size	kiss_fft_next_fast_size, .-kiss_fft_next_fast_size
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC0:
	.long	0
	.long	1071644672
	.align 8
.LC1:
	.long	0
	.long	1073741824
	.align 8
.LC2:
	.long	0
	.long	1074266112
	.align 8
.LC3:
	.long	0
	.long	1127219200
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC4:
	.long	4294967295
	.long	2147483647
	.long	0
	.long	0
	.section	.rodata.cst8
	.align 8
.LC5:
	.long	0
	.long	1072693248
	.align 8
.LC6:
	.long	1413754136
	.long	-1072094725
	.section	.rodata.cst16
	.align 16
.LC7:
	.long	0
	.long	-2147483648
	.long	0
	.long	0
	.ident	"GCC: (GNU) 4.8.1"
	.section	.note.GNU-stack,"",@progbits
