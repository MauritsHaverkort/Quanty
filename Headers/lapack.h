#ifdef BLASunderscore
#define FORTRAN_WRAPPER(x) x ## _

#define dlamch FORTRAN_WRAPPER(dlamch)
#define dlamc3 FORTRAN_WRAPPER(dlamc3)
#define slamch FORTRAN_WRAPPER(slamch)
#define slamc3 FORTRAN_WRAPPER(slamc3)
#define cbbcsd FORTRAN_WRAPPER(cbbcsd)
#define cbdsqr FORTRAN_WRAPPER(cbdsqr)
#define cgbbrd FORTRAN_WRAPPER(cgbbrd)
#define cgbcon FORTRAN_WRAPPER(cgbcon)
#define cgbequ FORTRAN_WRAPPER(cgbequ)
#define cgbequb FORTRAN_WRAPPER(cgbequb)
#define cgbrfs FORTRAN_WRAPPER(cgbrfs)
#define cgbsv FORTRAN_WRAPPER(cgbsv)
#define cgbsvx FORTRAN_WRAPPER(cgbsvx)
#define cgbtf2 FORTRAN_WRAPPER(cgbtf2)
#define cgbtrf FORTRAN_WRAPPER(cgbtrf)
#define cgbtrs FORTRAN_WRAPPER(cgbtrs)
#define cgebak FORTRAN_WRAPPER(cgebak)
#define cgebal FORTRAN_WRAPPER(cgebal)
#define cgebd2 FORTRAN_WRAPPER(cgebd2)
#define cgebrd FORTRAN_WRAPPER(cgebrd)
#define cgecon FORTRAN_WRAPPER(cgecon)
#define cgeequ FORTRAN_WRAPPER(cgeequ)
#define cgeequb FORTRAN_WRAPPER(cgeequb)
#define cgees FORTRAN_WRAPPER(cgees)
#define cgeesx FORTRAN_WRAPPER(cgeesx)
#define cgeev FORTRAN_WRAPPER(cgeev)
#define cgeevx FORTRAN_WRAPPER(cgeevx)
#define cgegs FORTRAN_WRAPPER(cgegs)
#define cgegv FORTRAN_WRAPPER(cgegv)
#define cgehd2 FORTRAN_WRAPPER(cgehd2)
#define cgehrd FORTRAN_WRAPPER(cgehrd)
#define cgelq2 FORTRAN_WRAPPER(cgelq2)
#define cgelqf FORTRAN_WRAPPER(cgelqf)
#define cgels FORTRAN_WRAPPER(cgels)
#define cgelsd FORTRAN_WRAPPER(cgelsd)
#define cgelss FORTRAN_WRAPPER(cgelss)
#define cgelsx FORTRAN_WRAPPER(cgelsx)
#define cgelsy FORTRAN_WRAPPER(cgelsy)
#define cgemqrt FORTRAN_WRAPPER(cgemqrt)
#define cgeql2 FORTRAN_WRAPPER(cgeql2)
#define cgeqlf FORTRAN_WRAPPER(cgeqlf)
#define cgeqp3 FORTRAN_WRAPPER(cgeqp3)
#define cgeqpf FORTRAN_WRAPPER(cgeqpf)
#define cgeqr2 FORTRAN_WRAPPER(cgeqr2)
#define cgeqr2p FORTRAN_WRAPPER(cgeqr2p)
#define cgeqrf FORTRAN_WRAPPER(cgeqrf)
#define cgeqrfp FORTRAN_WRAPPER(cgeqrfp)
#define cgeqrt FORTRAN_WRAPPER(cgeqrt)
#define cgeqrt2 FORTRAN_WRAPPER(cgeqrt2)
#define cgeqrt3 FORTRAN_WRAPPER(cgeqrt3)
#define cgerfs FORTRAN_WRAPPER(cgerfs)
#define cgerq2 FORTRAN_WRAPPER(cgerq2)
#define cgerqf FORTRAN_WRAPPER(cgerqf)
#define cgesc2 FORTRAN_WRAPPER(cgesc2)
#define cgesdd FORTRAN_WRAPPER(cgesdd)
#define cgesv FORTRAN_WRAPPER(cgesv)
#define cgesvd FORTRAN_WRAPPER(cgesvd)
#define cgesvx FORTRAN_WRAPPER(cgesvx)
#define cgetc2 FORTRAN_WRAPPER(cgetc2)
#define cgetf2 FORTRAN_WRAPPER(cgetf2)
#define cgetrf FORTRAN_WRAPPER(cgetrf)
#define cgetri FORTRAN_WRAPPER(cgetri)
#define cgetrs FORTRAN_WRAPPER(cgetrs)
#define cggbak FORTRAN_WRAPPER(cggbak)
#define cggbal FORTRAN_WRAPPER(cggbal)
#define cgges FORTRAN_WRAPPER(cgges)
#define cggesx FORTRAN_WRAPPER(cggesx)
#define cggev FORTRAN_WRAPPER(cggev)
#define cggevx FORTRAN_WRAPPER(cggevx)
#define cggglm FORTRAN_WRAPPER(cggglm)
#define cgghrd FORTRAN_WRAPPER(cgghrd)
#define cgglse FORTRAN_WRAPPER(cgglse)
#define cggqrf FORTRAN_WRAPPER(cggqrf)
#define cggrqf FORTRAN_WRAPPER(cggrqf)
#define cggsvd FORTRAN_WRAPPER(cggsvd)
#define cggsvp FORTRAN_WRAPPER(cggsvp)
#define cgtcon FORTRAN_WRAPPER(cgtcon)
#define cgtrfs FORTRAN_WRAPPER(cgtrfs)
#define cgtsv FORTRAN_WRAPPER(cgtsv)
#define cgtsvx FORTRAN_WRAPPER(cgtsvx)
#define cgttrf FORTRAN_WRAPPER(cgttrf)
#define cgttrs FORTRAN_WRAPPER(cgttrs)
#define cgtts2 FORTRAN_WRAPPER(cgtts2)
#define chbev FORTRAN_WRAPPER(chbev)
#define chbevd FORTRAN_WRAPPER(chbevd)
#define chbevx FORTRAN_WRAPPER(chbevx)
#define chbgst FORTRAN_WRAPPER(chbgst)
#define chbgv FORTRAN_WRAPPER(chbgv)
#define chbgvd FORTRAN_WRAPPER(chbgvd)
#define chbgvx FORTRAN_WRAPPER(chbgvx)
#define chbtrd FORTRAN_WRAPPER(chbtrd)
#define checon FORTRAN_WRAPPER(checon)
#define cheequb FORTRAN_WRAPPER(cheequb)
#define cheev FORTRAN_WRAPPER(cheev)
#define cheevd FORTRAN_WRAPPER(cheevd)
#define cheevr FORTRAN_WRAPPER(cheevr)
#define cheevx FORTRAN_WRAPPER(cheevx)
#define chegs2 FORTRAN_WRAPPER(chegs2)
#define chegst FORTRAN_WRAPPER(chegst)
#define chegv FORTRAN_WRAPPER(chegv)
#define chegvd FORTRAN_WRAPPER(chegvd)
#define chegvx FORTRAN_WRAPPER(chegvx)
#define cherfs FORTRAN_WRAPPER(cherfs)
#define chesv FORTRAN_WRAPPER(chesv)
#define chesvx FORTRAN_WRAPPER(chesvx)
#define cheswapr FORTRAN_WRAPPER(cheswapr)
#define chetd2 FORTRAN_WRAPPER(chetd2)
#define chetf2 FORTRAN_WRAPPER(chetf2)
#define chetrd FORTRAN_WRAPPER(chetrd)
#define chetrf FORTRAN_WRAPPER(chetrf)
#define chetri FORTRAN_WRAPPER(chetri)
#define chetri2 FORTRAN_WRAPPER(chetri2)
#define chetri2x FORTRAN_WRAPPER(chetri2x)
#define chetrs FORTRAN_WRAPPER(chetrs)
#define chetrs2 FORTRAN_WRAPPER(chetrs2)
#define chfrk FORTRAN_WRAPPER(chfrk)
#define chgeqz FORTRAN_WRAPPER(chgeqz)
#define chla_transtype FORTRAN_WRAPPER(chla_transtype)
#define chpcon FORTRAN_WRAPPER(chpcon)
#define chpev FORTRAN_WRAPPER(chpev)
#define chpevd FORTRAN_WRAPPER(chpevd)
#define chpevx FORTRAN_WRAPPER(chpevx)
#define chpgst FORTRAN_WRAPPER(chpgst)
#define chpgv FORTRAN_WRAPPER(chpgv)
#define chpgvd FORTRAN_WRAPPER(chpgvd)
#define chpgvx FORTRAN_WRAPPER(chpgvx)
#define chprfs FORTRAN_WRAPPER(chprfs)
#define chpsv FORTRAN_WRAPPER(chpsv)
#define chpsvx FORTRAN_WRAPPER(chpsvx)
#define chptrd FORTRAN_WRAPPER(chptrd)
#define chptrf FORTRAN_WRAPPER(chptrf)
#define chptri FORTRAN_WRAPPER(chptri)
#define chptrs FORTRAN_WRAPPER(chptrs)
#define chsein FORTRAN_WRAPPER(chsein)
#define chseqr FORTRAN_WRAPPER(chseqr)
#define clabrd FORTRAN_WRAPPER(clabrd)
#define clacgv FORTRAN_WRAPPER(clacgv)
#define clacn2 FORTRAN_WRAPPER(clacn2)
#define clacon FORTRAN_WRAPPER(clacon)
#define clacp2 FORTRAN_WRAPPER(clacp2)
#define clacpy FORTRAN_WRAPPER(clacpy)
#define clacrm FORTRAN_WRAPPER(clacrm)
#define clacrt FORTRAN_WRAPPER(clacrt)
#define cladiv FORTRAN_WRAPPER(cladiv)
#define claed0 FORTRAN_WRAPPER(claed0)
#define claed7 FORTRAN_WRAPPER(claed7)
#define claed8 FORTRAN_WRAPPER(claed8)
#define claein FORTRAN_WRAPPER(claein)
#define claesy FORTRAN_WRAPPER(claesy)
#define claev2 FORTRAN_WRAPPER(claev2)
#define clag2z FORTRAN_WRAPPER(clag2z)
#define clags2 FORTRAN_WRAPPER(clags2)
#define clagtm FORTRAN_WRAPPER(clagtm)
#define clahef FORTRAN_WRAPPER(clahef)
#define clahqr FORTRAN_WRAPPER(clahqr)
#define clahr2 FORTRAN_WRAPPER(clahr2)
#define clahrd FORTRAN_WRAPPER(clahrd)
#define claic1 FORTRAN_WRAPPER(claic1)
#define clals0 FORTRAN_WRAPPER(clals0)
#define clalsa FORTRAN_WRAPPER(clalsa)
#define clalsd FORTRAN_WRAPPER(clalsd)
#define clangb FORTRAN_WRAPPER(clangb)
#define clange FORTRAN_WRAPPER(clange)
#define clangt FORTRAN_WRAPPER(clangt)
#define clanhb FORTRAN_WRAPPER(clanhb)
#define clanhe FORTRAN_WRAPPER(clanhe)
#define clanhf FORTRAN_WRAPPER(clanhf)
#define clanhp FORTRAN_WRAPPER(clanhp)
#define clanhs FORTRAN_WRAPPER(clanhs)
#define clanht FORTRAN_WRAPPER(clanht)
#define clansb FORTRAN_WRAPPER(clansb)
#define clansp FORTRAN_WRAPPER(clansp)
#define clansy FORTRAN_WRAPPER(clansy)
#define clantb FORTRAN_WRAPPER(clantb)
#define clantp FORTRAN_WRAPPER(clantp)
#define clantr FORTRAN_WRAPPER(clantr)
#define clapll FORTRAN_WRAPPER(clapll)
#define clapmr FORTRAN_WRAPPER(clapmr)
#define clapmt FORTRAN_WRAPPER(clapmt)
#define claqgb FORTRAN_WRAPPER(claqgb)
#define claqge FORTRAN_WRAPPER(claqge)
#define claqhb FORTRAN_WRAPPER(claqhb)
#define claqhe FORTRAN_WRAPPER(claqhe)
#define claqhp FORTRAN_WRAPPER(claqhp)
#define claqp2 FORTRAN_WRAPPER(claqp2)
#define claqps FORTRAN_WRAPPER(claqps)
#define claqr0 FORTRAN_WRAPPER(claqr0)
#define claqr1 FORTRAN_WRAPPER(claqr1)
#define claqr2 FORTRAN_WRAPPER(claqr2)
#define claqr3 FORTRAN_WRAPPER(claqr3)
#define claqr4 FORTRAN_WRAPPER(claqr4)
#define claqr5 FORTRAN_WRAPPER(claqr5)
#define claqsb FORTRAN_WRAPPER(claqsb)
#define claqsp FORTRAN_WRAPPER(claqsp)
#define claqsy FORTRAN_WRAPPER(claqsy)
#define clar1v FORTRAN_WRAPPER(clar1v)
#define clar2v FORTRAN_WRAPPER(clar2v)
#define clarcm FORTRAN_WRAPPER(clarcm)
#define clarf FORTRAN_WRAPPER(clarf)
#define clarfb FORTRAN_WRAPPER(clarfb)
#define clarfg FORTRAN_WRAPPER(clarfg)
#define clarfgp FORTRAN_WRAPPER(clarfgp)
#define clarft FORTRAN_WRAPPER(clarft)
#define clarfx FORTRAN_WRAPPER(clarfx)
#define clargv FORTRAN_WRAPPER(clargv)
#define clarnv FORTRAN_WRAPPER(clarnv)
#define clarrv FORTRAN_WRAPPER(clarrv)
#define clartg FORTRAN_WRAPPER(clartg)
#define clartv FORTRAN_WRAPPER(clartv)
#define clarz FORTRAN_WRAPPER(clarz)
#define clarzb FORTRAN_WRAPPER(clarzb)
#define clarzt FORTRAN_WRAPPER(clarzt)
#define clascl FORTRAN_WRAPPER(clascl)
#define claset FORTRAN_WRAPPER(claset)
#define clasr FORTRAN_WRAPPER(clasr)
#define classq FORTRAN_WRAPPER(classq)
#define claswp FORTRAN_WRAPPER(claswp)
#define clasyf FORTRAN_WRAPPER(clasyf)
#define clatbs FORTRAN_WRAPPER(clatbs)
#define clatdf FORTRAN_WRAPPER(clatdf)
#define clatps FORTRAN_WRAPPER(clatps)
#define clatrd FORTRAN_WRAPPER(clatrd)
#define clatrs FORTRAN_WRAPPER(clatrs)
#define clatrz FORTRAN_WRAPPER(clatrz)
#define clatzm FORTRAN_WRAPPER(clatzm)
#define clauu2 FORTRAN_WRAPPER(clauu2)
#define clauum FORTRAN_WRAPPER(clauum)
#define cpbcon FORTRAN_WRAPPER(cpbcon)
#define cpbequ FORTRAN_WRAPPER(cpbequ)
#define cpbrfs FORTRAN_WRAPPER(cpbrfs)
#define cpbstf FORTRAN_WRAPPER(cpbstf)
#define cpbsv FORTRAN_WRAPPER(cpbsv)
#define cpbsvx FORTRAN_WRAPPER(cpbsvx)
#define cpbtf2 FORTRAN_WRAPPER(cpbtf2)
#define cpbtrf FORTRAN_WRAPPER(cpbtrf)
#define cpbtrs FORTRAN_WRAPPER(cpbtrs)
#define cpftrf FORTRAN_WRAPPER(cpftrf)
#define cpftri FORTRAN_WRAPPER(cpftri)
#define cpftrs FORTRAN_WRAPPER(cpftrs)
#define cpocon FORTRAN_WRAPPER(cpocon)
#define cpoequ FORTRAN_WRAPPER(cpoequ)
#define cpoequb FORTRAN_WRAPPER(cpoequb)
#define cporfs FORTRAN_WRAPPER(cporfs)
#define cposv FORTRAN_WRAPPER(cposv)
#define cposvx FORTRAN_WRAPPER(cposvx)
#define cpotf2 FORTRAN_WRAPPER(cpotf2)
#define cpotrf FORTRAN_WRAPPER(cpotrf)
#define cpotri FORTRAN_WRAPPER(cpotri)
#define cpotrs FORTRAN_WRAPPER(cpotrs)
#define cppcon FORTRAN_WRAPPER(cppcon)
#define cppequ FORTRAN_WRAPPER(cppequ)
#define cpprfs FORTRAN_WRAPPER(cpprfs)
#define cppsv FORTRAN_WRAPPER(cppsv)
#define cppsvx FORTRAN_WRAPPER(cppsvx)
#define cpptrf FORTRAN_WRAPPER(cpptrf)
#define cpptri FORTRAN_WRAPPER(cpptri)
#define cpptrs FORTRAN_WRAPPER(cpptrs)
#define cpstf2 FORTRAN_WRAPPER(cpstf2)
#define cpstrf FORTRAN_WRAPPER(cpstrf)
#define cptcon FORTRAN_WRAPPER(cptcon)
#define cpteqr FORTRAN_WRAPPER(cpteqr)
#define cptrfs FORTRAN_WRAPPER(cptrfs)
#define cptsv FORTRAN_WRAPPER(cptsv)
#define cptsvx FORTRAN_WRAPPER(cptsvx)
#define cpttrf FORTRAN_WRAPPER(cpttrf)
#define cpttrs FORTRAN_WRAPPER(cpttrs)
#define cptts2 FORTRAN_WRAPPER(cptts2)
#define crot FORTRAN_WRAPPER(crot)
#define cspcon FORTRAN_WRAPPER(cspcon)
#define cspmv FORTRAN_WRAPPER(cspmv)
#define cspr FORTRAN_WRAPPER(cspr)
#define csprfs FORTRAN_WRAPPER(csprfs)
#define cspsv FORTRAN_WRAPPER(cspsv)
#define cspsvx FORTRAN_WRAPPER(cspsvx)
#define csptrf FORTRAN_WRAPPER(csptrf)
#define csptri FORTRAN_WRAPPER(csptri)
#define csptrs FORTRAN_WRAPPER(csptrs)
#define csrscl FORTRAN_WRAPPER(csrscl)
#define cstedc FORTRAN_WRAPPER(cstedc)
#define cstegr FORTRAN_WRAPPER(cstegr)
#define cstein FORTRAN_WRAPPER(cstein)
#define cstemr FORTRAN_WRAPPER(cstemr)
#define csteqr FORTRAN_WRAPPER(csteqr)
#define csycon FORTRAN_WRAPPER(csycon)
#define csyconv FORTRAN_WRAPPER(csyconv)
#define csyequb FORTRAN_WRAPPER(csyequb)
#define csymv FORTRAN_WRAPPER(csymv)
#define csyr FORTRAN_WRAPPER(csyr)
#define csyrfs FORTRAN_WRAPPER(csyrfs)
#define csysv FORTRAN_WRAPPER(csysv)
#define csysvx FORTRAN_WRAPPER(csysvx)
#define csyswapr FORTRAN_WRAPPER(csyswapr)
#define csytf2 FORTRAN_WRAPPER(csytf2)
#define csytrf FORTRAN_WRAPPER(csytrf)
#define csytri FORTRAN_WRAPPER(csytri)
#define csytri2 FORTRAN_WRAPPER(csytri2)
#define csytri2x FORTRAN_WRAPPER(csytri2x)
#define csytrs FORTRAN_WRAPPER(csytrs)
#define csytrs2 FORTRAN_WRAPPER(csytrs2)
#define ctbcon FORTRAN_WRAPPER(ctbcon)
#define ctbrfs FORTRAN_WRAPPER(ctbrfs)
#define ctbtrs FORTRAN_WRAPPER(ctbtrs)
#define ctfsm FORTRAN_WRAPPER(ctfsm)
#define ctftri FORTRAN_WRAPPER(ctftri)
#define ctfttp FORTRAN_WRAPPER(ctfttp)
#define ctfttr FORTRAN_WRAPPER(ctfttr)
#define ctgevc FORTRAN_WRAPPER(ctgevc)
#define ctgex2 FORTRAN_WRAPPER(ctgex2)
#define ctgexc FORTRAN_WRAPPER(ctgexc)
#define ctgsen FORTRAN_WRAPPER(ctgsen)
#define ctgsja FORTRAN_WRAPPER(ctgsja)
#define ctgsna FORTRAN_WRAPPER(ctgsna)
#define ctgsy2 FORTRAN_WRAPPER(ctgsy2)
#define ctgsyl FORTRAN_WRAPPER(ctgsyl)
#define ctpcon FORTRAN_WRAPPER(ctpcon)
#define ctpmqrt FORTRAN_WRAPPER(ctpmqrt)
#define ctpqrt FORTRAN_WRAPPER(ctpqrt)
#define ctpqrt2 FORTRAN_WRAPPER(ctpqrt2)
#define ctprfb FORTRAN_WRAPPER(ctprfb)
#define ctprfs FORTRAN_WRAPPER(ctprfs)
#define ctptri FORTRAN_WRAPPER(ctptri)
#define ctptrs FORTRAN_WRAPPER(ctptrs)
#define ctpttf FORTRAN_WRAPPER(ctpttf)
#define ctpttr FORTRAN_WRAPPER(ctpttr)
#define ctrcon FORTRAN_WRAPPER(ctrcon)
#define ctrevc FORTRAN_WRAPPER(ctrevc)
#define ctrexc FORTRAN_WRAPPER(ctrexc)
#define ctrrfs FORTRAN_WRAPPER(ctrrfs)
#define ctrsen FORTRAN_WRAPPER(ctrsen)
#define ctrsna FORTRAN_WRAPPER(ctrsna)
#define ctrsyl FORTRAN_WRAPPER(ctrsyl)
#define ctrti2 FORTRAN_WRAPPER(ctrti2)
#define ctrtri FORTRAN_WRAPPER(ctrtri)
#define ctrtrs FORTRAN_WRAPPER(ctrtrs)
#define ctrttf FORTRAN_WRAPPER(ctrttf)
#define ctrttp FORTRAN_WRAPPER(ctrttp)
#define ctzrqf FORTRAN_WRAPPER(ctzrqf)
#define ctzrzf FORTRAN_WRAPPER(ctzrzf)
#define cunbdb FORTRAN_WRAPPER(cunbdb)
#define cuncsd FORTRAN_WRAPPER(cuncsd)
#define cung2l FORTRAN_WRAPPER(cung2l)
#define cung2r FORTRAN_WRAPPER(cung2r)
#define cungbr FORTRAN_WRAPPER(cungbr)
#define cunghr FORTRAN_WRAPPER(cunghr)
#define cungl2 FORTRAN_WRAPPER(cungl2)
#define cunglq FORTRAN_WRAPPER(cunglq)
#define cungql FORTRAN_WRAPPER(cungql)
#define cungqr FORTRAN_WRAPPER(cungqr)
#define cungr2 FORTRAN_WRAPPER(cungr2)
#define cungrq FORTRAN_WRAPPER(cungrq)
#define cungtr FORTRAN_WRAPPER(cungtr)
#define cunm2l FORTRAN_WRAPPER(cunm2l)
#define cunm2r FORTRAN_WRAPPER(cunm2r)
#define cunmbr FORTRAN_WRAPPER(cunmbr)
#define cunmhr FORTRAN_WRAPPER(cunmhr)
#define cunml2 FORTRAN_WRAPPER(cunml2)
#define cunmlq FORTRAN_WRAPPER(cunmlq)
#define cunmql FORTRAN_WRAPPER(cunmql)
#define cunmqr FORTRAN_WRAPPER(cunmqr)
#define cunmr2 FORTRAN_WRAPPER(cunmr2)
#define cunmr3 FORTRAN_WRAPPER(cunmr3)
#define cunmrq FORTRAN_WRAPPER(cunmrq)
#define cunmrz FORTRAN_WRAPPER(cunmrz)
#define cunmtr FORTRAN_WRAPPER(cunmtr)
#define cupgtr FORTRAN_WRAPPER(cupgtr)
#define cupmtr FORTRAN_WRAPPER(cupmtr)
#define dbbcsd FORTRAN_WRAPPER(dbbcsd)
#define dbdsdc FORTRAN_WRAPPER(dbdsdc)
#define dbdsqr FORTRAN_WRAPPER(dbdsqr)
#define ddisna FORTRAN_WRAPPER(ddisna)
#define dgbbrd FORTRAN_WRAPPER(dgbbrd)
#define dgbcon FORTRAN_WRAPPER(dgbcon)
#define dgbequ FORTRAN_WRAPPER(dgbequ)
#define dgbequb FORTRAN_WRAPPER(dgbequb)
#define dgbrfs FORTRAN_WRAPPER(dgbrfs)
#define dgbsv FORTRAN_WRAPPER(dgbsv)
#define dgbsvx FORTRAN_WRAPPER(dgbsvx)
#define dgbtf2 FORTRAN_WRAPPER(dgbtf2)
#define dgbtrf FORTRAN_WRAPPER(dgbtrf)
#define dgbtrs FORTRAN_WRAPPER(dgbtrs)
#define dgebak FORTRAN_WRAPPER(dgebak)
#define dgebal FORTRAN_WRAPPER(dgebal)
#define dgebd2 FORTRAN_WRAPPER(dgebd2)
#define dgebrd FORTRAN_WRAPPER(dgebrd)
#define dgecon FORTRAN_WRAPPER(dgecon)
#define dgeequ FORTRAN_WRAPPER(dgeequ)
#define dgeequb FORTRAN_WRAPPER(dgeequb)
#define dgees FORTRAN_WRAPPER(dgees)
#define dgeesx FORTRAN_WRAPPER(dgeesx)
#define dgeev FORTRAN_WRAPPER(dgeev)
#define dgeevx FORTRAN_WRAPPER(dgeevx)
#define dgegs FORTRAN_WRAPPER(dgegs)
#define dgegv FORTRAN_WRAPPER(dgegv)
#define dgehd2 FORTRAN_WRAPPER(dgehd2)
#define dgehrd FORTRAN_WRAPPER(dgehrd)
#define dgejsv FORTRAN_WRAPPER(dgejsv)
#define dgelq2 FORTRAN_WRAPPER(dgelq2)
#define dgelqf FORTRAN_WRAPPER(dgelqf)
#define dgels FORTRAN_WRAPPER(dgels)
#define dgelsd FORTRAN_WRAPPER(dgelsd)
#define dgelss FORTRAN_WRAPPER(dgelss)
#define dgelsx FORTRAN_WRAPPER(dgelsx)
#define dgelsy FORTRAN_WRAPPER(dgelsy)
#define dgemqrt FORTRAN_WRAPPER(dgemqrt)
#define dgeql2 FORTRAN_WRAPPER(dgeql2)
#define dgeqlf FORTRAN_WRAPPER(dgeqlf)
#define dgeqp3 FORTRAN_WRAPPER(dgeqp3)
#define dgeqpf FORTRAN_WRAPPER(dgeqpf)
#define dgeqr2 FORTRAN_WRAPPER(dgeqr2)
#define dgeqr2p FORTRAN_WRAPPER(dgeqr2p)
#define dgeqrf FORTRAN_WRAPPER(dgeqrf)
#define dgeqrfp FORTRAN_WRAPPER(dgeqrfp)
#define dgeqrt FORTRAN_WRAPPER(dgeqrt)
#define dgeqrt2 FORTRAN_WRAPPER(dgeqrt2)
#define dgeqrt3 FORTRAN_WRAPPER(dgeqrt3)
#define dgerfs FORTRAN_WRAPPER(dgerfs)
#define dgerq2 FORTRAN_WRAPPER(dgerq2)
#define dgerqf FORTRAN_WRAPPER(dgerqf)
#define dgesc2 FORTRAN_WRAPPER(dgesc2)
#define dgesdd FORTRAN_WRAPPER(dgesdd)
#define dgesv FORTRAN_WRAPPER(dgesv)
#define dgesvd FORTRAN_WRAPPER(dgesvd)
#define dgesvj FORTRAN_WRAPPER(dgesvj)
#define dgesvx FORTRAN_WRAPPER(dgesvx)
#define dgetc2 FORTRAN_WRAPPER(dgetc2)
#define dgetf2 FORTRAN_WRAPPER(dgetf2)
#define dgetrf FORTRAN_WRAPPER(dgetrf)
#define dgetri FORTRAN_WRAPPER(dgetri)
#define dgetrs FORTRAN_WRAPPER(dgetrs)
#define dggbak FORTRAN_WRAPPER(dggbak)
#define dggbal FORTRAN_WRAPPER(dggbal)
#define dgges FORTRAN_WRAPPER(dgges)
#define dggesx FORTRAN_WRAPPER(dggesx)
#define dggev FORTRAN_WRAPPER(dggev)
#define dggevx FORTRAN_WRAPPER(dggevx)
#define dggglm FORTRAN_WRAPPER(dggglm)
#define dgghrd FORTRAN_WRAPPER(dgghrd)
#define dgglse FORTRAN_WRAPPER(dgglse)
#define dggqrf FORTRAN_WRAPPER(dggqrf)
#define dggrqf FORTRAN_WRAPPER(dggrqf)
#define dggsvd FORTRAN_WRAPPER(dggsvd)
#define dggsvp FORTRAN_WRAPPER(dggsvp)
#define dgsvj0 FORTRAN_WRAPPER(dgsvj0)
#define dgsvj1 FORTRAN_WRAPPER(dgsvj1)
#define dgtcon FORTRAN_WRAPPER(dgtcon)
#define dgtrfs FORTRAN_WRAPPER(dgtrfs)
#define dgtsv FORTRAN_WRAPPER(dgtsv)
#define dgtsvx FORTRAN_WRAPPER(dgtsvx)
#define dgttrf FORTRAN_WRAPPER(dgttrf)
#define dgttrs FORTRAN_WRAPPER(dgttrs)
#define dgtts2 FORTRAN_WRAPPER(dgtts2)
#define dhgeqz FORTRAN_WRAPPER(dhgeqz)
#define dhsein FORTRAN_WRAPPER(dhsein)
#define dhseqr FORTRAN_WRAPPER(dhseqr)
#define disnan FORTRAN_WRAPPER(disnan)
#define dlabad FORTRAN_WRAPPER(dlabad)
#define dlabrd FORTRAN_WRAPPER(dlabrd)
#define dlacn2 FORTRAN_WRAPPER(dlacn2)
#define dlacon FORTRAN_WRAPPER(dlacon)
#define dlacpy FORTRAN_WRAPPER(dlacpy)
#define dladiv FORTRAN_WRAPPER(dladiv)
#define dlae2 FORTRAN_WRAPPER(dlae2)
#define dlaebz FORTRAN_WRAPPER(dlaebz)
#define dlaed0 FORTRAN_WRAPPER(dlaed0)
#define dlaed1 FORTRAN_WRAPPER(dlaed1)
#define dlaed2 FORTRAN_WRAPPER(dlaed2)
#define dlaed3 FORTRAN_WRAPPER(dlaed3)
#define dlaed4 FORTRAN_WRAPPER(dlaed4)
#define dlaed5 FORTRAN_WRAPPER(dlaed5)
#define dlaed6 FORTRAN_WRAPPER(dlaed6)
#define dlaed7 FORTRAN_WRAPPER(dlaed7)
#define dlaed8 FORTRAN_WRAPPER(dlaed8)
#define dlaed9 FORTRAN_WRAPPER(dlaed9)
#define dlaeda FORTRAN_WRAPPER(dlaeda)
#define dlaein FORTRAN_WRAPPER(dlaein)
#define dlaev2 FORTRAN_WRAPPER(dlaev2)
#define dlaexc FORTRAN_WRAPPER(dlaexc)
#define dlag2 FORTRAN_WRAPPER(dlag2)
#define dlag2s FORTRAN_WRAPPER(dlag2s)
#define dlags2 FORTRAN_WRAPPER(dlags2)
#define dlagtf FORTRAN_WRAPPER(dlagtf)
#define dlagtm FORTRAN_WRAPPER(dlagtm)
#define dlagts FORTRAN_WRAPPER(dlagts)
#define dlagv2 FORTRAN_WRAPPER(dlagv2)
#define dlahqr FORTRAN_WRAPPER(dlahqr)
#define dlahr2 FORTRAN_WRAPPER(dlahr2)
#define dlahrd FORTRAN_WRAPPER(dlahrd)
#define dlaic1 FORTRAN_WRAPPER(dlaic1)
#define dlaisnan FORTRAN_WRAPPER(dlaisnan)
#define dlaln2 FORTRAN_WRAPPER(dlaln2)
#define dlals0 FORTRAN_WRAPPER(dlals0)
#define dlalsa FORTRAN_WRAPPER(dlalsa)
#define dlalsd FORTRAN_WRAPPER(dlalsd)
#define dlamrg FORTRAN_WRAPPER(dlamrg)
#define dlaneg FORTRAN_WRAPPER(dlaneg)
#define dlangb FORTRAN_WRAPPER(dlangb)
#define dlange FORTRAN_WRAPPER(dlange)
#define dlangt FORTRAN_WRAPPER(dlangt)
#define dlanhs FORTRAN_WRAPPER(dlanhs)
#define dlansb FORTRAN_WRAPPER(dlansb)
#define dlansf FORTRAN_WRAPPER(dlansf)
#define dlansp FORTRAN_WRAPPER(dlansp)
#define dlanst FORTRAN_WRAPPER(dlanst)
#define dlansy FORTRAN_WRAPPER(dlansy)
#define dlantb FORTRAN_WRAPPER(dlantb)
#define dlantp FORTRAN_WRAPPER(dlantp)
#define dlantr FORTRAN_WRAPPER(dlantr)
#define dlanv2 FORTRAN_WRAPPER(dlanv2)
#define dlapll FORTRAN_WRAPPER(dlapll)
#define dlapmr FORTRAN_WRAPPER(dlapmr)
#define dlapmt FORTRAN_WRAPPER(dlapmt)
#define dlapy2 FORTRAN_WRAPPER(dlapy2)
#define dlapy3 FORTRAN_WRAPPER(dlapy3)
#define dlaqgb FORTRAN_WRAPPER(dlaqgb)
#define dlaqge FORTRAN_WRAPPER(dlaqge)
#define dlaqp2 FORTRAN_WRAPPER(dlaqp2)
#define dlaqps FORTRAN_WRAPPER(dlaqps)
#define dlaqr0 FORTRAN_WRAPPER(dlaqr0)
#define dlaqr1 FORTRAN_WRAPPER(dlaqr1)
#define dlaqr2 FORTRAN_WRAPPER(dlaqr2)
#define dlaqr3 FORTRAN_WRAPPER(dlaqr3)
#define dlaqr4 FORTRAN_WRAPPER(dlaqr4)
#define dlaqr5 FORTRAN_WRAPPER(dlaqr5)
#define dlaqsb FORTRAN_WRAPPER(dlaqsb)
#define dlaqsp FORTRAN_WRAPPER(dlaqsp)
#define dlaqsy FORTRAN_WRAPPER(dlaqsy)
#define dlaqtr FORTRAN_WRAPPER(dlaqtr)
#define dlar1v FORTRAN_WRAPPER(dlar1v)
#define dlar2v FORTRAN_WRAPPER(dlar2v)
#define dlarf FORTRAN_WRAPPER(dlarf)
#define dlarfb FORTRAN_WRAPPER(dlarfb)
#define dlarfg FORTRAN_WRAPPER(dlarfg)
#define dlarfgp FORTRAN_WRAPPER(dlarfgp)
#define dlarft FORTRAN_WRAPPER(dlarft)
#define dlarfx FORTRAN_WRAPPER(dlarfx)
#define dlargv FORTRAN_WRAPPER(dlargv)
#define dlarnv FORTRAN_WRAPPER(dlarnv)
#define dlarra FORTRAN_WRAPPER(dlarra)
#define dlarrb FORTRAN_WRAPPER(dlarrb)
#define dlarrc FORTRAN_WRAPPER(dlarrc)
#define dlarrd FORTRAN_WRAPPER(dlarrd)
#define dlarre FORTRAN_WRAPPER(dlarre)
#define dlarrf FORTRAN_WRAPPER(dlarrf)
#define dlarrj FORTRAN_WRAPPER(dlarrj)
#define dlarrk FORTRAN_WRAPPER(dlarrk)
#define dlarrr FORTRAN_WRAPPER(dlarrr)
#define dlarrv FORTRAN_WRAPPER(dlarrv)
#define dlartg FORTRAN_WRAPPER(dlartg)
#define dlartgp FORTRAN_WRAPPER(dlartgp)
#define dlartgs FORTRAN_WRAPPER(dlartgs)
#define dlartv FORTRAN_WRAPPER(dlartv)
#define dlaruv FORTRAN_WRAPPER(dlaruv)
#define dlarz FORTRAN_WRAPPER(dlarz)
#define dlarzb FORTRAN_WRAPPER(dlarzb)
#define dlarzt FORTRAN_WRAPPER(dlarzt)
#define dlas2 FORTRAN_WRAPPER(dlas2)
#define dlascl FORTRAN_WRAPPER(dlascl)
#define dlasd0 FORTRAN_WRAPPER(dlasd0)
#define dlasd1 FORTRAN_WRAPPER(dlasd1)
#define dlasd2 FORTRAN_WRAPPER(dlasd2)
#define dlasd3 FORTRAN_WRAPPER(dlasd3)
#define dlasd4 FORTRAN_WRAPPER(dlasd4)
#define dlasd5 FORTRAN_WRAPPER(dlasd5)
#define dlasd6 FORTRAN_WRAPPER(dlasd6)
#define dlasd7 FORTRAN_WRAPPER(dlasd7)
#define dlasd8 FORTRAN_WRAPPER(dlasd8)
#define dlasda FORTRAN_WRAPPER(dlasda)
#define dlasdq FORTRAN_WRAPPER(dlasdq)
#define dlasdt FORTRAN_WRAPPER(dlasdt)
#define dlaset FORTRAN_WRAPPER(dlaset)
#define dlasq1 FORTRAN_WRAPPER(dlasq1)
#define dlasq2 FORTRAN_WRAPPER(dlasq2)
#define dlasq3 FORTRAN_WRAPPER(dlasq3)
#define dlasq4 FORTRAN_WRAPPER(dlasq4)
#define dlasq5 FORTRAN_WRAPPER(dlasq5)
#define dlasq6 FORTRAN_WRAPPER(dlasq6)
#define dlasr FORTRAN_WRAPPER(dlasr)
#define dlasrt FORTRAN_WRAPPER(dlasrt)
#define dlassq FORTRAN_WRAPPER(dlassq)
#define dlasv2 FORTRAN_WRAPPER(dlasv2)
#define dlaswp FORTRAN_WRAPPER(dlaswp)
#define dlasy2 FORTRAN_WRAPPER(dlasy2)
#define dlasyf FORTRAN_WRAPPER(dlasyf)
#define dlat2s FORTRAN_WRAPPER(dlat2s)
#define dlatbs FORTRAN_WRAPPER(dlatbs)
#define dlatdf FORTRAN_WRAPPER(dlatdf)
#define dlatps FORTRAN_WRAPPER(dlatps)
#define dlatrd FORTRAN_WRAPPER(dlatrd)
#define dlatrs FORTRAN_WRAPPER(dlatrs)
#define dlatrz FORTRAN_WRAPPER(dlatrz)
#define dlatzm FORTRAN_WRAPPER(dlatzm)
#define dlauu2 FORTRAN_WRAPPER(dlauu2)
#define dlauum FORTRAN_WRAPPER(dlauum)
#define dopgtr FORTRAN_WRAPPER(dopgtr)
#define dopmtr FORTRAN_WRAPPER(dopmtr)
#define dorbdb FORTRAN_WRAPPER(dorbdb)
#define dorcsd FORTRAN_WRAPPER(dorcsd)
#define dorg2l FORTRAN_WRAPPER(dorg2l)
#define dorg2r FORTRAN_WRAPPER(dorg2r)
#define dorgbr FORTRAN_WRAPPER(dorgbr)
#define dorghr FORTRAN_WRAPPER(dorghr)
#define dorgl2 FORTRAN_WRAPPER(dorgl2)
#define dorglq FORTRAN_WRAPPER(dorglq)
#define dorgql FORTRAN_WRAPPER(dorgql)
#define dorgqr FORTRAN_WRAPPER(dorgqr)
#define dorgr2 FORTRAN_WRAPPER(dorgr2)
#define dorgrq FORTRAN_WRAPPER(dorgrq)
#define dorgtr FORTRAN_WRAPPER(dorgtr)
#define dorm2l FORTRAN_WRAPPER(dorm2l)
#define dorm2r FORTRAN_WRAPPER(dorm2r)
#define dormbr FORTRAN_WRAPPER(dormbr)
#define dormhr FORTRAN_WRAPPER(dormhr)
#define dorml2 FORTRAN_WRAPPER(dorml2)
#define dormlq FORTRAN_WRAPPER(dormlq)
#define dormql FORTRAN_WRAPPER(dormql)
#define dormqr FORTRAN_WRAPPER(dormqr)
#define dormr2 FORTRAN_WRAPPER(dormr2)
#define dormr3 FORTRAN_WRAPPER(dormr3)
#define dormrq FORTRAN_WRAPPER(dormrq)
#define dormrz FORTRAN_WRAPPER(dormrz)
#define dormtr FORTRAN_WRAPPER(dormtr)
#define dpbcon FORTRAN_WRAPPER(dpbcon)
#define dpbequ FORTRAN_WRAPPER(dpbequ)
#define dpbrfs FORTRAN_WRAPPER(dpbrfs)
#define dpbstf FORTRAN_WRAPPER(dpbstf)
#define dpbsv FORTRAN_WRAPPER(dpbsv)
#define dpbsvx FORTRAN_WRAPPER(dpbsvx)
#define dpbtf2 FORTRAN_WRAPPER(dpbtf2)
#define dpbtrf FORTRAN_WRAPPER(dpbtrf)
#define dpbtrs FORTRAN_WRAPPER(dpbtrs)
#define dpftrf FORTRAN_WRAPPER(dpftrf)
#define dpftri FORTRAN_WRAPPER(dpftri)
#define dpftrs FORTRAN_WRAPPER(dpftrs)
#define dpocon FORTRAN_WRAPPER(dpocon)
#define dpoequ FORTRAN_WRAPPER(dpoequ)
#define dpoequb FORTRAN_WRAPPER(dpoequb)
#define dporfs FORTRAN_WRAPPER(dporfs)
#define dposv FORTRAN_WRAPPER(dposv)
#define dposvx FORTRAN_WRAPPER(dposvx)
#define dpotf2 FORTRAN_WRAPPER(dpotf2)
#define dpotrf FORTRAN_WRAPPER(dpotrf)
#define dpotri FORTRAN_WRAPPER(dpotri)
#define dpotrs FORTRAN_WRAPPER(dpotrs)
#define dppcon FORTRAN_WRAPPER(dppcon)
#define dppequ FORTRAN_WRAPPER(dppequ)
#define dpprfs FORTRAN_WRAPPER(dpprfs)
#define dppsv FORTRAN_WRAPPER(dppsv)
#define dppsvx FORTRAN_WRAPPER(dppsvx)
#define dpptrf FORTRAN_WRAPPER(dpptrf)
#define dpptri FORTRAN_WRAPPER(dpptri)
#define dpptrs FORTRAN_WRAPPER(dpptrs)
#define dpstf2 FORTRAN_WRAPPER(dpstf2)
#define dpstrf FORTRAN_WRAPPER(dpstrf)
#define dptcon FORTRAN_WRAPPER(dptcon)
#define dpteqr FORTRAN_WRAPPER(dpteqr)
#define dptrfs FORTRAN_WRAPPER(dptrfs)
#define dptsv FORTRAN_WRAPPER(dptsv)
#define dptsvx FORTRAN_WRAPPER(dptsvx)
#define dpttrf FORTRAN_WRAPPER(dpttrf)
#define dpttrs FORTRAN_WRAPPER(dpttrs)
#define dptts2 FORTRAN_WRAPPER(dptts2)
#define drscl FORTRAN_WRAPPER(drscl)
#define dsbev FORTRAN_WRAPPER(dsbev)
#define dsbevd FORTRAN_WRAPPER(dsbevd)
#define dsbevx FORTRAN_WRAPPER(dsbevx)
#define dsbgst FORTRAN_WRAPPER(dsbgst)
#define dsbgv FORTRAN_WRAPPER(dsbgv)
#define dsbgvd FORTRAN_WRAPPER(dsbgvd)
#define dsbgvx FORTRAN_WRAPPER(dsbgvx)
#define dsbtrd FORTRAN_WRAPPER(dsbtrd)
#define dsfrk FORTRAN_WRAPPER(dsfrk)
#define dsgesv FORTRAN_WRAPPER(dsgesv)
#define dspcon FORTRAN_WRAPPER(dspcon)
#define dspev FORTRAN_WRAPPER(dspev)
#define dspevd FORTRAN_WRAPPER(dspevd)
#define dspevx FORTRAN_WRAPPER(dspevx)
#define dspgst FORTRAN_WRAPPER(dspgst)
#define dspgv FORTRAN_WRAPPER(dspgv)
#define dspgvd FORTRAN_WRAPPER(dspgvd)
#define dspgvx FORTRAN_WRAPPER(dspgvx)
#define dsposv FORTRAN_WRAPPER(dsposv)
#define dsprfs FORTRAN_WRAPPER(dsprfs)
#define dspsv FORTRAN_WRAPPER(dspsv)
#define dspsvx FORTRAN_WRAPPER(dspsvx)
#define dsptrd FORTRAN_WRAPPER(dsptrd)
#define dsptrf FORTRAN_WRAPPER(dsptrf)
#define dsptri FORTRAN_WRAPPER(dsptri)
#define dsptrs FORTRAN_WRAPPER(dsptrs)
#define dstebz FORTRAN_WRAPPER(dstebz)
#define dstedc FORTRAN_WRAPPER(dstedc)
#define dstegr FORTRAN_WRAPPER(dstegr)
#define dstein FORTRAN_WRAPPER(dstein)
#define dstemr FORTRAN_WRAPPER(dstemr)
#define dsteqr FORTRAN_WRAPPER(dsteqr)
#define dsterf FORTRAN_WRAPPER(dsterf)
#define dstev FORTRAN_WRAPPER(dstev)
#define dstevd FORTRAN_WRAPPER(dstevd)
#define dstevr FORTRAN_WRAPPER(dstevr)
#define dstevx FORTRAN_WRAPPER(dstevx)
#define dsycon FORTRAN_WRAPPER(dsycon)
#define dsyconv FORTRAN_WRAPPER(dsyconv)
#define dsyequb FORTRAN_WRAPPER(dsyequb)
#define dsyev FORTRAN_WRAPPER(dsyev)
#define dsyevd FORTRAN_WRAPPER(dsyevd)
#define dsyevr FORTRAN_WRAPPER(dsyevr)
#define dsyevx FORTRAN_WRAPPER(dsyevx)
#define dsygs2 FORTRAN_WRAPPER(dsygs2)
#define dsygst FORTRAN_WRAPPER(dsygst)
#define dsygv FORTRAN_WRAPPER(dsygv)
#define dsygvd FORTRAN_WRAPPER(dsygvd)
#define dsygvx FORTRAN_WRAPPER(dsygvx)
#define dsyrfs FORTRAN_WRAPPER(dsyrfs)
#define dsysv FORTRAN_WRAPPER(dsysv)
#define dsysvx FORTRAN_WRAPPER(dsysvx)
#define dsyswapr FORTRAN_WRAPPER(dsyswapr)
#define dsytd2 FORTRAN_WRAPPER(dsytd2)
#define dsytf2 FORTRAN_WRAPPER(dsytf2)
#define dsytrd FORTRAN_WRAPPER(dsytrd)
#define dsytrf FORTRAN_WRAPPER(dsytrf)
#define dsytri FORTRAN_WRAPPER(dsytri)
#define dsytri2 FORTRAN_WRAPPER(dsytri2)
#define dsytri2x FORTRAN_WRAPPER(dsytri2x)
#define dsytrs FORTRAN_WRAPPER(dsytrs)
#define dsytrs2 FORTRAN_WRAPPER(dsytrs2)
#define dtbcon FORTRAN_WRAPPER(dtbcon)
#define dtbrfs FORTRAN_WRAPPER(dtbrfs)
#define dtbtrs FORTRAN_WRAPPER(dtbtrs)
#define dtfsm FORTRAN_WRAPPER(dtfsm)
#define dtftri FORTRAN_WRAPPER(dtftri)
#define dtfttp FORTRAN_WRAPPER(dtfttp)
#define dtfttr FORTRAN_WRAPPER(dtfttr)
#define dtgevc FORTRAN_WRAPPER(dtgevc)
#define dtgex2 FORTRAN_WRAPPER(dtgex2)
#define dtgexc FORTRAN_WRAPPER(dtgexc)
#define dtgsen FORTRAN_WRAPPER(dtgsen)
#define dtgsja FORTRAN_WRAPPER(dtgsja)
#define dtgsna FORTRAN_WRAPPER(dtgsna)
#define dtgsy2 FORTRAN_WRAPPER(dtgsy2)
#define dtgsyl FORTRAN_WRAPPER(dtgsyl)
#define dtpcon FORTRAN_WRAPPER(dtpcon)
#define dtpmqrt FORTRAN_WRAPPER(dtpmqrt)
#define dtpqrt FORTRAN_WRAPPER(dtpqrt)
#define dtpqrt2 FORTRAN_WRAPPER(dtpqrt2)
#define dtprfb FORTRAN_WRAPPER(dtprfb)
#define dtprfs FORTRAN_WRAPPER(dtprfs)
#define dtptri FORTRAN_WRAPPER(dtptri)
#define dtptrs FORTRAN_WRAPPER(dtptrs)
#define dtpttf FORTRAN_WRAPPER(dtpttf)
#define dtpttr FORTRAN_WRAPPER(dtpttr)
#define dtrcon FORTRAN_WRAPPER(dtrcon)
#define dtrevc FORTRAN_WRAPPER(dtrevc)
#define dtrexc FORTRAN_WRAPPER(dtrexc)
#define dtrrfs FORTRAN_WRAPPER(dtrrfs)
#define dtrsen FORTRAN_WRAPPER(dtrsen)
#define dtrsna FORTRAN_WRAPPER(dtrsna)
#define dtrsyl FORTRAN_WRAPPER(dtrsyl)
#define dtrti2 FORTRAN_WRAPPER(dtrti2)
#define dtrtri FORTRAN_WRAPPER(dtrtri)
#define dtrtrs FORTRAN_WRAPPER(dtrtrs)
#define dtrttf FORTRAN_WRAPPER(dtrttf)
#define dtrttp FORTRAN_WRAPPER(dtrttp)
#define dtzrqf FORTRAN_WRAPPER(dtzrqf)
#define dtzrzf FORTRAN_WRAPPER(dtzrzf)
#define dzsum1 FORTRAN_WRAPPER(dzsum1)
#define icmax1 FORTRAN_WRAPPER(icmax1)
#define ieeeck FORTRAN_WRAPPER(ieeeck)
#define ilaclc FORTRAN_WRAPPER(ilaclc)
#define ilaclr FORTRAN_WRAPPER(ilaclr)
#define iladiag FORTRAN_WRAPPER(iladiag)
#define iladlc FORTRAN_WRAPPER(iladlc)
#define iladlr FORTRAN_WRAPPER(iladlr)
#define ilaenv FORTRAN_WRAPPER(ilaenv)
#define ilaprec FORTRAN_WRAPPER(ilaprec)
#define ilaslc FORTRAN_WRAPPER(ilaslc)
#define ilaslr FORTRAN_WRAPPER(ilaslr)
#define ilatrans FORTRAN_WRAPPER(ilatrans)
#define ilauplo FORTRAN_WRAPPER(ilauplo)
#define ilaver FORTRAN_WRAPPER(ilaver)
#define ilazlc FORTRAN_WRAPPER(ilazlc)
#define ilazlr FORTRAN_WRAPPER(ilazlr)
#define iparmq FORTRAN_WRAPPER(iparmq)
#define izmax1 FORTRAN_WRAPPER(izmax1)
#define lsamen FORTRAN_WRAPPER(lsamen)
#define sbbcsd FORTRAN_WRAPPER(sbbcsd)
#define sbdsdc FORTRAN_WRAPPER(sbdsdc)
#define sbdsqr FORTRAN_WRAPPER(sbdsqr)
#define scsum1 FORTRAN_WRAPPER(scsum1)
#define sdisna FORTRAN_WRAPPER(sdisna)
#define sgbbrd FORTRAN_WRAPPER(sgbbrd)
#define sgbcon FORTRAN_WRAPPER(sgbcon)
#define sgbequ FORTRAN_WRAPPER(sgbequ)
#define sgbequb FORTRAN_WRAPPER(sgbequb)
#define sgbrfs FORTRAN_WRAPPER(sgbrfs)
#define sgbsv FORTRAN_WRAPPER(sgbsv)
#define sgbsvx FORTRAN_WRAPPER(sgbsvx)
#define sgbtf2 FORTRAN_WRAPPER(sgbtf2)
#define sgbtrf FORTRAN_WRAPPER(sgbtrf)
#define sgbtrs FORTRAN_WRAPPER(sgbtrs)
#define sgebak FORTRAN_WRAPPER(sgebak)
#define sgebal FORTRAN_WRAPPER(sgebal)
#define sgebd2 FORTRAN_WRAPPER(sgebd2)
#define sgebrd FORTRAN_WRAPPER(sgebrd)
#define sgecon FORTRAN_WRAPPER(sgecon)
#define sgeequ FORTRAN_WRAPPER(sgeequ)
#define sgeequb FORTRAN_WRAPPER(sgeequb)
#define sgees FORTRAN_WRAPPER(sgees)
#define sgeesx FORTRAN_WRAPPER(sgeesx)
#define sgeev FORTRAN_WRAPPER(sgeev)
#define sgeevx FORTRAN_WRAPPER(sgeevx)
#define sgegs FORTRAN_WRAPPER(sgegs)
#define sgegv FORTRAN_WRAPPER(sgegv)
#define sgehd2 FORTRAN_WRAPPER(sgehd2)
#define sgehrd FORTRAN_WRAPPER(sgehrd)
#define sgejsv FORTRAN_WRAPPER(sgejsv)
#define sgelq2 FORTRAN_WRAPPER(sgelq2)
#define sgelqf FORTRAN_WRAPPER(sgelqf)
#define sgels FORTRAN_WRAPPER(sgels)
#define sgelsd FORTRAN_WRAPPER(sgelsd)
#define sgelss FORTRAN_WRAPPER(sgelss)
#define sgelsx FORTRAN_WRAPPER(sgelsx)
#define sgelsy FORTRAN_WRAPPER(sgelsy)
#define sgemqrt FORTRAN_WRAPPER(sgemqrt)
#define sgeql2 FORTRAN_WRAPPER(sgeql2)
#define sgeqlf FORTRAN_WRAPPER(sgeqlf)
#define sgeqp3 FORTRAN_WRAPPER(sgeqp3)
#define sgeqpf FORTRAN_WRAPPER(sgeqpf)
#define sgeqr2 FORTRAN_WRAPPER(sgeqr2)
#define sgeqr2p FORTRAN_WRAPPER(sgeqr2p)
#define sgeqrf FORTRAN_WRAPPER(sgeqrf)
#define sgeqrfp FORTRAN_WRAPPER(sgeqrfp)
#define sgeqrt FORTRAN_WRAPPER(sgeqrt)
#define sgeqrt2 FORTRAN_WRAPPER(sgeqrt2)
#define sgeqrt3 FORTRAN_WRAPPER(sgeqrt3)
#define sgerfs FORTRAN_WRAPPER(sgerfs)
#define sgerq2 FORTRAN_WRAPPER(sgerq2)
#define sgerqf FORTRAN_WRAPPER(sgerqf)
#define sgesc2 FORTRAN_WRAPPER(sgesc2)
#define sgesdd FORTRAN_WRAPPER(sgesdd)
#define sgesv FORTRAN_WRAPPER(sgesv)
#define sgesvd FORTRAN_WRAPPER(sgesvd)
#define sgesvj FORTRAN_WRAPPER(sgesvj)
#define sgesvx FORTRAN_WRAPPER(sgesvx)
#define sgetc2 FORTRAN_WRAPPER(sgetc2)
#define sgetf2 FORTRAN_WRAPPER(sgetf2)
#define sgetrf FORTRAN_WRAPPER(sgetrf)
#define sgetri FORTRAN_WRAPPER(sgetri)
#define sgetrs FORTRAN_WRAPPER(sgetrs)
#define sggbak FORTRAN_WRAPPER(sggbak)
#define sggbal FORTRAN_WRAPPER(sggbal)
#define sgges FORTRAN_WRAPPER(sgges)
#define sggesx FORTRAN_WRAPPER(sggesx)
#define sggev FORTRAN_WRAPPER(sggev)
#define sggevx FORTRAN_WRAPPER(sggevx)
#define sggglm FORTRAN_WRAPPER(sggglm)
#define sgghrd FORTRAN_WRAPPER(sgghrd)
#define sgglse FORTRAN_WRAPPER(sgglse)
#define sggqrf FORTRAN_WRAPPER(sggqrf)
#define sggrqf FORTRAN_WRAPPER(sggrqf)
#define sggsvd FORTRAN_WRAPPER(sggsvd)
#define sggsvp FORTRAN_WRAPPER(sggsvp)
#define sgsvj0 FORTRAN_WRAPPER(sgsvj0)
#define sgsvj1 FORTRAN_WRAPPER(sgsvj1)
#define sgtcon FORTRAN_WRAPPER(sgtcon)
#define sgtrfs FORTRAN_WRAPPER(sgtrfs)
#define sgtsv FORTRAN_WRAPPER(sgtsv)
#define sgtsvx FORTRAN_WRAPPER(sgtsvx)
#define sgttrf FORTRAN_WRAPPER(sgttrf)
#define sgttrs FORTRAN_WRAPPER(sgttrs)
#define sgtts2 FORTRAN_WRAPPER(sgtts2)
#define shgeqz FORTRAN_WRAPPER(shgeqz)
#define shsein FORTRAN_WRAPPER(shsein)
#define shseqr FORTRAN_WRAPPER(shseqr)
#define sisnan FORTRAN_WRAPPER(sisnan)
#define slabad FORTRAN_WRAPPER(slabad)
#define slabrd FORTRAN_WRAPPER(slabrd)
#define slacn2 FORTRAN_WRAPPER(slacn2)
#define slacon FORTRAN_WRAPPER(slacon)
#define slacpy FORTRAN_WRAPPER(slacpy)
#define sladiv FORTRAN_WRAPPER(sladiv)
#define slae2 FORTRAN_WRAPPER(slae2)
#define slaebz FORTRAN_WRAPPER(slaebz)
#define slaed0 FORTRAN_WRAPPER(slaed0)
#define slaed1 FORTRAN_WRAPPER(slaed1)
#define slaed2 FORTRAN_WRAPPER(slaed2)
#define slaed3 FORTRAN_WRAPPER(slaed3)
#define slaed4 FORTRAN_WRAPPER(slaed4)
#define slaed5 FORTRAN_WRAPPER(slaed5)
#define slaed6 FORTRAN_WRAPPER(slaed6)
#define slaed7 FORTRAN_WRAPPER(slaed7)
#define slaed8 FORTRAN_WRAPPER(slaed8)
#define slaed9 FORTRAN_WRAPPER(slaed9)
#define slaeda FORTRAN_WRAPPER(slaeda)
#define slaein FORTRAN_WRAPPER(slaein)
#define slaev2 FORTRAN_WRAPPER(slaev2)
#define slaexc FORTRAN_WRAPPER(slaexc)
#define slag2 FORTRAN_WRAPPER(slag2)
#define slag2d FORTRAN_WRAPPER(slag2d)
#define slags2 FORTRAN_WRAPPER(slags2)
#define slagtf FORTRAN_WRAPPER(slagtf)
#define slagtm FORTRAN_WRAPPER(slagtm)
#define slagts FORTRAN_WRAPPER(slagts)
#define slagv2 FORTRAN_WRAPPER(slagv2)
#define slahqr FORTRAN_WRAPPER(slahqr)
#define slahr2 FORTRAN_WRAPPER(slahr2)
#define slahrd FORTRAN_WRAPPER(slahrd)
#define slaic1 FORTRAN_WRAPPER(slaic1)
#define slaisnan FORTRAN_WRAPPER(slaisnan)
#define slaln2 FORTRAN_WRAPPER(slaln2)
#define slals0 FORTRAN_WRAPPER(slals0)
#define slalsa FORTRAN_WRAPPER(slalsa)
#define slalsd FORTRAN_WRAPPER(slalsd)
#define slamrg FORTRAN_WRAPPER(slamrg)
#define slaneg FORTRAN_WRAPPER(slaneg)
#define slangb FORTRAN_WRAPPER(slangb)
#define slange FORTRAN_WRAPPER(slange)
#define slangt FORTRAN_WRAPPER(slangt)
#define slanhs FORTRAN_WRAPPER(slanhs)
#define slansb FORTRAN_WRAPPER(slansb)
#define slansf FORTRAN_WRAPPER(slansf)
#define slansp FORTRAN_WRAPPER(slansp)
#define slanst FORTRAN_WRAPPER(slanst)
#define slansy FORTRAN_WRAPPER(slansy)
#define slantb FORTRAN_WRAPPER(slantb)
#define slantp FORTRAN_WRAPPER(slantp)
#define slantr FORTRAN_WRAPPER(slantr)
#define slanv2 FORTRAN_WRAPPER(slanv2)
#define slapll FORTRAN_WRAPPER(slapll)
#define slapmr FORTRAN_WRAPPER(slapmr)
#define slapmt FORTRAN_WRAPPER(slapmt)
#define slapy2 FORTRAN_WRAPPER(slapy2)
#define slapy3 FORTRAN_WRAPPER(slapy3)
#define slaqgb FORTRAN_WRAPPER(slaqgb)
#define slaqge FORTRAN_WRAPPER(slaqge)
#define slaqp2 FORTRAN_WRAPPER(slaqp2)
#define slaqps FORTRAN_WRAPPER(slaqps)
#define slaqr0 FORTRAN_WRAPPER(slaqr0)
#define slaqr1 FORTRAN_WRAPPER(slaqr1)
#define slaqr2 FORTRAN_WRAPPER(slaqr2)
#define slaqr3 FORTRAN_WRAPPER(slaqr3)
#define slaqr4 FORTRAN_WRAPPER(slaqr4)
#define slaqr5 FORTRAN_WRAPPER(slaqr5)
#define slaqsb FORTRAN_WRAPPER(slaqsb)
#define slaqsp FORTRAN_WRAPPER(slaqsp)
#define slaqsy FORTRAN_WRAPPER(slaqsy)
#define slaqtr FORTRAN_WRAPPER(slaqtr)
#define slar1v FORTRAN_WRAPPER(slar1v)
#define slar2v FORTRAN_WRAPPER(slar2v)
#define slarf FORTRAN_WRAPPER(slarf)
#define slarfb FORTRAN_WRAPPER(slarfb)
#define slarfg FORTRAN_WRAPPER(slarfg)
#define slarfgp FORTRAN_WRAPPER(slarfgp)
#define slarft FORTRAN_WRAPPER(slarft)
#define slarfx FORTRAN_WRAPPER(slarfx)
#define slargv FORTRAN_WRAPPER(slargv)
#define slarnv FORTRAN_WRAPPER(slarnv)
#define slarra FORTRAN_WRAPPER(slarra)
#define slarrb FORTRAN_WRAPPER(slarrb)
#define slarrc FORTRAN_WRAPPER(slarrc)
#define slarrd FORTRAN_WRAPPER(slarrd)
#define slarre FORTRAN_WRAPPER(slarre)
#define slarrf FORTRAN_WRAPPER(slarrf)
#define slarrj FORTRAN_WRAPPER(slarrj)
#define slarrk FORTRAN_WRAPPER(slarrk)
#define slarrr FORTRAN_WRAPPER(slarrr)
#define slarrv FORTRAN_WRAPPER(slarrv)
#define slartg FORTRAN_WRAPPER(slartg)
#define slartgp FORTRAN_WRAPPER(slartgp)
#define slartgs FORTRAN_WRAPPER(slartgs)
#define slartv FORTRAN_WRAPPER(slartv)
#define slaruv FORTRAN_WRAPPER(slaruv)
#define slarz FORTRAN_WRAPPER(slarz)
#define slarzb FORTRAN_WRAPPER(slarzb)
#define slarzt FORTRAN_WRAPPER(slarzt)
#define slas2 FORTRAN_WRAPPER(slas2)
#define slascl FORTRAN_WRAPPER(slascl)
#define slasd0 FORTRAN_WRAPPER(slasd0)
#define slasd1 FORTRAN_WRAPPER(slasd1)
#define slasd2 FORTRAN_WRAPPER(slasd2)
#define slasd3 FORTRAN_WRAPPER(slasd3)
#define slasd4 FORTRAN_WRAPPER(slasd4)
#define slasd5 FORTRAN_WRAPPER(slasd5)
#define slasd6 FORTRAN_WRAPPER(slasd6)
#define slasd7 FORTRAN_WRAPPER(slasd7)
#define slasd8 FORTRAN_WRAPPER(slasd8)
#define slasda FORTRAN_WRAPPER(slasda)
#define slasdq FORTRAN_WRAPPER(slasdq)
#define slasdt FORTRAN_WRAPPER(slasdt)
#define slaset FORTRAN_WRAPPER(slaset)
#define slasq1 FORTRAN_WRAPPER(slasq1)
#define slasq2 FORTRAN_WRAPPER(slasq2)
#define slasq3 FORTRAN_WRAPPER(slasq3)
#define slasq4 FORTRAN_WRAPPER(slasq4)
#define slasq5 FORTRAN_WRAPPER(slasq5)
#define slasq6 FORTRAN_WRAPPER(slasq6)
#define slasr FORTRAN_WRAPPER(slasr)
#define slasrt FORTRAN_WRAPPER(slasrt)
#define slassq FORTRAN_WRAPPER(slassq)
#define slasv2 FORTRAN_WRAPPER(slasv2)
#define slaswp FORTRAN_WRAPPER(slaswp)
#define slasy2 FORTRAN_WRAPPER(slasy2)
#define slasyf FORTRAN_WRAPPER(slasyf)
#define slatbs FORTRAN_WRAPPER(slatbs)
#define slatdf FORTRAN_WRAPPER(slatdf)
#define slatps FORTRAN_WRAPPER(slatps)
#define slatrd FORTRAN_WRAPPER(slatrd)
#define slatrs FORTRAN_WRAPPER(slatrs)
#define slatrz FORTRAN_WRAPPER(slatrz)
#define slatzm FORTRAN_WRAPPER(slatzm)
#define slauu2 FORTRAN_WRAPPER(slauu2)
#define slauum FORTRAN_WRAPPER(slauum)
#define sopgtr FORTRAN_WRAPPER(sopgtr)
#define sopmtr FORTRAN_WRAPPER(sopmtr)
#define sorbdb FORTRAN_WRAPPER(sorbdb)
#define sorcsd FORTRAN_WRAPPER(sorcsd)
#define sorg2l FORTRAN_WRAPPER(sorg2l)
#define sorg2r FORTRAN_WRAPPER(sorg2r)
#define sorgbr FORTRAN_WRAPPER(sorgbr)
#define sorghr FORTRAN_WRAPPER(sorghr)
#define sorgl2 FORTRAN_WRAPPER(sorgl2)
#define sorglq FORTRAN_WRAPPER(sorglq)
#define sorgql FORTRAN_WRAPPER(sorgql)
#define sorgqr FORTRAN_WRAPPER(sorgqr)
#define sorgr2 FORTRAN_WRAPPER(sorgr2)
#define sorgrq FORTRAN_WRAPPER(sorgrq)
#define sorgtr FORTRAN_WRAPPER(sorgtr)
#define sorm2l FORTRAN_WRAPPER(sorm2l)
#define sorm2r FORTRAN_WRAPPER(sorm2r)
#define sormbr FORTRAN_WRAPPER(sormbr)
#define sormhr FORTRAN_WRAPPER(sormhr)
#define sorml2 FORTRAN_WRAPPER(sorml2)
#define sormlq FORTRAN_WRAPPER(sormlq)
#define sormql FORTRAN_WRAPPER(sormql)
#define sormqr FORTRAN_WRAPPER(sormqr)
#define sormr2 FORTRAN_WRAPPER(sormr2)
#define sormr3 FORTRAN_WRAPPER(sormr3)
#define sormrq FORTRAN_WRAPPER(sormrq)
#define sormrz FORTRAN_WRAPPER(sormrz)
#define sormtr FORTRAN_WRAPPER(sormtr)
#define spbcon FORTRAN_WRAPPER(spbcon)
#define spbequ FORTRAN_WRAPPER(spbequ)
#define spbrfs FORTRAN_WRAPPER(spbrfs)
#define spbstf FORTRAN_WRAPPER(spbstf)
#define spbsv FORTRAN_WRAPPER(spbsv)
#define spbsvx FORTRAN_WRAPPER(spbsvx)
#define spbtf2 FORTRAN_WRAPPER(spbtf2)
#define spbtrf FORTRAN_WRAPPER(spbtrf)
#define spbtrs FORTRAN_WRAPPER(spbtrs)
#define spftrf FORTRAN_WRAPPER(spftrf)
#define spftri FORTRAN_WRAPPER(spftri)
#define spftrs FORTRAN_WRAPPER(spftrs)
#define spocon FORTRAN_WRAPPER(spocon)
#define spoequ FORTRAN_WRAPPER(spoequ)
#define spoequb FORTRAN_WRAPPER(spoequb)
#define sporfs FORTRAN_WRAPPER(sporfs)
#define sposv FORTRAN_WRAPPER(sposv)
#define sposvx FORTRAN_WRAPPER(sposvx)
#define spotf2 FORTRAN_WRAPPER(spotf2)
#define spotrf FORTRAN_WRAPPER(spotrf)
#define spotri FORTRAN_WRAPPER(spotri)
#define spotrs FORTRAN_WRAPPER(spotrs)
#define sppcon FORTRAN_WRAPPER(sppcon)
#define sppequ FORTRAN_WRAPPER(sppequ)
#define spprfs FORTRAN_WRAPPER(spprfs)
#define sppsv FORTRAN_WRAPPER(sppsv)
#define sppsvx FORTRAN_WRAPPER(sppsvx)
#define spptrf FORTRAN_WRAPPER(spptrf)
#define spptri FORTRAN_WRAPPER(spptri)
#define spptrs FORTRAN_WRAPPER(spptrs)
#define spstf2 FORTRAN_WRAPPER(spstf2)
#define spstrf FORTRAN_WRAPPER(spstrf)
#define sptcon FORTRAN_WRAPPER(sptcon)
#define spteqr FORTRAN_WRAPPER(spteqr)
#define sptrfs FORTRAN_WRAPPER(sptrfs)
#define sptsv FORTRAN_WRAPPER(sptsv)
#define sptsvx FORTRAN_WRAPPER(sptsvx)
#define spttrf FORTRAN_WRAPPER(spttrf)
#define spttrs FORTRAN_WRAPPER(spttrs)
#define sptts2 FORTRAN_WRAPPER(sptts2)
#define srscl FORTRAN_WRAPPER(srscl)
#define ssbev FORTRAN_WRAPPER(ssbev)
#define ssbevd FORTRAN_WRAPPER(ssbevd)
#define ssbevx FORTRAN_WRAPPER(ssbevx)
#define ssbgst FORTRAN_WRAPPER(ssbgst)
#define ssbgv FORTRAN_WRAPPER(ssbgv)
#define ssbgvd FORTRAN_WRAPPER(ssbgvd)
#define ssbgvx FORTRAN_WRAPPER(ssbgvx)
#define ssbtrd FORTRAN_WRAPPER(ssbtrd)
#define ssfrk FORTRAN_WRAPPER(ssfrk)
#define sspcon FORTRAN_WRAPPER(sspcon)
#define sspev FORTRAN_WRAPPER(sspev)
#define sspevd FORTRAN_WRAPPER(sspevd)
#define sspevx FORTRAN_WRAPPER(sspevx)
#define sspgst FORTRAN_WRAPPER(sspgst)
#define sspgv FORTRAN_WRAPPER(sspgv)
#define sspgvd FORTRAN_WRAPPER(sspgvd)
#define sspgvx FORTRAN_WRAPPER(sspgvx)
#define ssprfs FORTRAN_WRAPPER(ssprfs)
#define sspsv FORTRAN_WRAPPER(sspsv)
#define sspsvx FORTRAN_WRAPPER(sspsvx)
#define ssptrd FORTRAN_WRAPPER(ssptrd)
#define ssptrf FORTRAN_WRAPPER(ssptrf)
#define ssptri FORTRAN_WRAPPER(ssptri)
#define ssptrs FORTRAN_WRAPPER(ssptrs)
#define sstebz FORTRAN_WRAPPER(sstebz)
#define sstedc FORTRAN_WRAPPER(sstedc)
#define sstegr FORTRAN_WRAPPER(sstegr)
#define sstein FORTRAN_WRAPPER(sstein)
#define sstemr FORTRAN_WRAPPER(sstemr)
#define ssteqr FORTRAN_WRAPPER(ssteqr)
#define ssterf FORTRAN_WRAPPER(ssterf)
#define sstev FORTRAN_WRAPPER(sstev)
#define sstevd FORTRAN_WRAPPER(sstevd)
#define sstevr FORTRAN_WRAPPER(sstevr)
#define sstevx FORTRAN_WRAPPER(sstevx)
#define ssycon FORTRAN_WRAPPER(ssycon)
#define ssyconv FORTRAN_WRAPPER(ssyconv)
#define ssyequb FORTRAN_WRAPPER(ssyequb)
#define ssyev FORTRAN_WRAPPER(ssyev)
#define ssyevd FORTRAN_WRAPPER(ssyevd)
#define ssyevr FORTRAN_WRAPPER(ssyevr)
#define ssyevx FORTRAN_WRAPPER(ssyevx)
#define ssygs2 FORTRAN_WRAPPER(ssygs2)
#define ssygst FORTRAN_WRAPPER(ssygst)
#define ssygv FORTRAN_WRAPPER(ssygv)
#define ssygvd FORTRAN_WRAPPER(ssygvd)
#define ssygvx FORTRAN_WRAPPER(ssygvx)
#define ssyrfs FORTRAN_WRAPPER(ssyrfs)
#define ssysv FORTRAN_WRAPPER(ssysv)
#define ssysvx FORTRAN_WRAPPER(ssysvx)
#define ssyswapr FORTRAN_WRAPPER(ssyswapr)
#define ssytd2 FORTRAN_WRAPPER(ssytd2)
#define ssytf2 FORTRAN_WRAPPER(ssytf2)
#define ssytrd FORTRAN_WRAPPER(ssytrd)
#define ssytrf FORTRAN_WRAPPER(ssytrf)
#define ssytri FORTRAN_WRAPPER(ssytri)
#define ssytri2 FORTRAN_WRAPPER(ssytri2)
#define ssytri2x FORTRAN_WRAPPER(ssytri2x)
#define ssytrs FORTRAN_WRAPPER(ssytrs)
#define ssytrs2 FORTRAN_WRAPPER(ssytrs2)
#define stbcon FORTRAN_WRAPPER(stbcon)
#define stbrfs FORTRAN_WRAPPER(stbrfs)
#define stbtrs FORTRAN_WRAPPER(stbtrs)
#define stfsm FORTRAN_WRAPPER(stfsm)
#define stftri FORTRAN_WRAPPER(stftri)
#define stfttp FORTRAN_WRAPPER(stfttp)
#define stfttr FORTRAN_WRAPPER(stfttr)
#define stgevc FORTRAN_WRAPPER(stgevc)
#define stgex2 FORTRAN_WRAPPER(stgex2)
#define stgexc FORTRAN_WRAPPER(stgexc)
#define stgsen FORTRAN_WRAPPER(stgsen)
#define stgsja FORTRAN_WRAPPER(stgsja)
#define stgsna FORTRAN_WRAPPER(stgsna)
#define stgsy2 FORTRAN_WRAPPER(stgsy2)
#define stgsyl FORTRAN_WRAPPER(stgsyl)
#define stpcon FORTRAN_WRAPPER(stpcon)
#define stpmqrt FORTRAN_WRAPPER(stpmqrt)
#define stpqrt FORTRAN_WRAPPER(stpqrt)
#define stpqrt2 FORTRAN_WRAPPER(stpqrt2)
#define stprfb FORTRAN_WRAPPER(stprfb)
#define stprfs FORTRAN_WRAPPER(stprfs)
#define stptri FORTRAN_WRAPPER(stptri)
#define stptrs FORTRAN_WRAPPER(stptrs)
#define stpttf FORTRAN_WRAPPER(stpttf)
#define stpttr FORTRAN_WRAPPER(stpttr)
#define strcon FORTRAN_WRAPPER(strcon)
#define strevc FORTRAN_WRAPPER(strevc)
#define strexc FORTRAN_WRAPPER(strexc)
#define strrfs FORTRAN_WRAPPER(strrfs)
#define strsen FORTRAN_WRAPPER(strsen)
#define strsna FORTRAN_WRAPPER(strsna)
#define strsyl FORTRAN_WRAPPER(strsyl)
#define strti2 FORTRAN_WRAPPER(strti2)
#define strtri FORTRAN_WRAPPER(strtri)
#define strtrs FORTRAN_WRAPPER(strtrs)
#define strttf FORTRAN_WRAPPER(strttf)
#define strttp FORTRAN_WRAPPER(strttp)
#define stzrqf FORTRAN_WRAPPER(stzrqf)
#define stzrzf FORTRAN_WRAPPER(stzrzf)
#define xerbla_array FORTRAN_WRAPPER(xerbla_array)
#define zbbcsd FORTRAN_WRAPPER(zbbcsd)
#define zbdsqr FORTRAN_WRAPPER(zbdsqr)
#define zcgesv FORTRAN_WRAPPER(zcgesv)
#define zcposv FORTRAN_WRAPPER(zcposv)
#define zdrscl FORTRAN_WRAPPER(zdrscl)
#define zgbbrd FORTRAN_WRAPPER(zgbbrd)
#define zgbcon FORTRAN_WRAPPER(zgbcon)
#define zgbequ FORTRAN_WRAPPER(zgbequ)
#define zgbequb FORTRAN_WRAPPER(zgbequb)
#define zgbrfs FORTRAN_WRAPPER(zgbrfs)
#define zgbsv FORTRAN_WRAPPER(zgbsv)
#define zgbsvx FORTRAN_WRAPPER(zgbsvx)
#define zgbtf2 FORTRAN_WRAPPER(zgbtf2)
#define zgbtrf FORTRAN_WRAPPER(zgbtrf)
#define zgbtrs FORTRAN_WRAPPER(zgbtrs)
#define zgebak FORTRAN_WRAPPER(zgebak)
#define zgebal FORTRAN_WRAPPER(zgebal)
#define zgebd2 FORTRAN_WRAPPER(zgebd2)
#define zgebrd FORTRAN_WRAPPER(zgebrd)
#define zgecon FORTRAN_WRAPPER(zgecon)
#define zgeequ FORTRAN_WRAPPER(zgeequ)
#define zgeequb FORTRAN_WRAPPER(zgeequb)
#define zgees FORTRAN_WRAPPER(zgees)
#define zgeesx FORTRAN_WRAPPER(zgeesx)
#define zgeev FORTRAN_WRAPPER(zgeev)
#define zgeevx FORTRAN_WRAPPER(zgeevx)
#define zgegs FORTRAN_WRAPPER(zgegs)
#define zgegv FORTRAN_WRAPPER(zgegv)
#define zgehd2 FORTRAN_WRAPPER(zgehd2)
#define zgehrd FORTRAN_WRAPPER(zgehrd)
#define zgelq2 FORTRAN_WRAPPER(zgelq2)
#define zgelqf FORTRAN_WRAPPER(zgelqf)
#define zgels FORTRAN_WRAPPER(zgels)
#define zgelsd FORTRAN_WRAPPER(zgelsd)
#define zgelss FORTRAN_WRAPPER(zgelss)
#define zgelsx FORTRAN_WRAPPER(zgelsx)
#define zgelsy FORTRAN_WRAPPER(zgelsy)
#define zgemqrt FORTRAN_WRAPPER(zgemqrt)
#define zgeql2 FORTRAN_WRAPPER(zgeql2)
#define zgeqlf FORTRAN_WRAPPER(zgeqlf)
#define zgeqp3 FORTRAN_WRAPPER(zgeqp3)
#define zgeqpf FORTRAN_WRAPPER(zgeqpf)
#define zgeqr2 FORTRAN_WRAPPER(zgeqr2)
#define zgeqr2p FORTRAN_WRAPPER(zgeqr2p)
#define zgeqrf FORTRAN_WRAPPER(zgeqrf)
#define zgeqrfp FORTRAN_WRAPPER(zgeqrfp)
#define zgeqrt FORTRAN_WRAPPER(zgeqrt)
#define zgeqrt2 FORTRAN_WRAPPER(zgeqrt2)
#define zgeqrt3 FORTRAN_WRAPPER(zgeqrt3)
#define zgerfs FORTRAN_WRAPPER(zgerfs)
#define zgerq2 FORTRAN_WRAPPER(zgerq2)
#define zgerqf FORTRAN_WRAPPER(zgerqf)
#define zgesc2 FORTRAN_WRAPPER(zgesc2)
#define zgesdd FORTRAN_WRAPPER(zgesdd)
#define zgesv FORTRAN_WRAPPER(zgesv)
#define zgesvd FORTRAN_WRAPPER(zgesvd)
#define zgesvx FORTRAN_WRAPPER(zgesvx)
#define zgetc2 FORTRAN_WRAPPER(zgetc2)
#define zgetf2 FORTRAN_WRAPPER(zgetf2)
#define zgetrf FORTRAN_WRAPPER(zgetrf)
#define zgetri FORTRAN_WRAPPER(zgetri)
#define zgetrs FORTRAN_WRAPPER(zgetrs)
#define zggbak FORTRAN_WRAPPER(zggbak)
#define zggbal FORTRAN_WRAPPER(zggbal)
#define zgges FORTRAN_WRAPPER(zgges)
#define zggesx FORTRAN_WRAPPER(zggesx)
#define zggev FORTRAN_WRAPPER(zggev)
#define zggevx FORTRAN_WRAPPER(zggevx)
#define zggglm FORTRAN_WRAPPER(zggglm)
#define zgghrd FORTRAN_WRAPPER(zgghrd)
#define zgglse FORTRAN_WRAPPER(zgglse)
#define zggqrf FORTRAN_WRAPPER(zggqrf)
#define zggrqf FORTRAN_WRAPPER(zggrqf)
#define zggsvd FORTRAN_WRAPPER(zggsvd)
#define zggsvp FORTRAN_WRAPPER(zggsvp)
#define zgtcon FORTRAN_WRAPPER(zgtcon)
#define zgtrfs FORTRAN_WRAPPER(zgtrfs)
#define zgtsv FORTRAN_WRAPPER(zgtsv)
#define zgtsvx FORTRAN_WRAPPER(zgtsvx)
#define zgttrf FORTRAN_WRAPPER(zgttrf)
#define zgttrs FORTRAN_WRAPPER(zgttrs)
#define zgtts2 FORTRAN_WRAPPER(zgtts2)
#define zhbev FORTRAN_WRAPPER(zhbev)
#define zhbevd FORTRAN_WRAPPER(zhbevd)
#define zhbevx FORTRAN_WRAPPER(zhbevx)
#define zhbgst FORTRAN_WRAPPER(zhbgst)
#define zhbgv FORTRAN_WRAPPER(zhbgv)
#define zhbgvd FORTRAN_WRAPPER(zhbgvd)
#define zhbgvx FORTRAN_WRAPPER(zhbgvx)
#define zhbtrd FORTRAN_WRAPPER(zhbtrd)
#define zhecon FORTRAN_WRAPPER(zhecon)
#define zheequb FORTRAN_WRAPPER(zheequb)
#define zheev FORTRAN_WRAPPER(zheev)
#define zheevd FORTRAN_WRAPPER(zheevd)
#define zheevr FORTRAN_WRAPPER(zheevr)
#define zheevx FORTRAN_WRAPPER(zheevx)
#define zhegs2 FORTRAN_WRAPPER(zhegs2)
#define zhegst FORTRAN_WRAPPER(zhegst)
#define zhegv FORTRAN_WRAPPER(zhegv)
#define zhegvd FORTRAN_WRAPPER(zhegvd)
#define zhegvx FORTRAN_WRAPPER(zhegvx)
#define zherfs FORTRAN_WRAPPER(zherfs)
#define zhesv FORTRAN_WRAPPER(zhesv)
#define zhesvx FORTRAN_WRAPPER(zhesvx)
#define zheswapr FORTRAN_WRAPPER(zheswapr)
#define zhetd2 FORTRAN_WRAPPER(zhetd2)
#define zhetf2 FORTRAN_WRAPPER(zhetf2)
#define zhetrd FORTRAN_WRAPPER(zhetrd)
#define zhetrf FORTRAN_WRAPPER(zhetrf)
#define zhetri FORTRAN_WRAPPER(zhetri)
#define zhetri2 FORTRAN_WRAPPER(zhetri2)
#define zhetri2x FORTRAN_WRAPPER(zhetri2x)
#define zhetrs FORTRAN_WRAPPER(zhetrs)
#define zhetrs2 FORTRAN_WRAPPER(zhetrs2)
#define zhfrk FORTRAN_WRAPPER(zhfrk)
#define zhgeqz FORTRAN_WRAPPER(zhgeqz)
#define zhpcon FORTRAN_WRAPPER(zhpcon)
#define zhpev FORTRAN_WRAPPER(zhpev)
#define zhpevd FORTRAN_WRAPPER(zhpevd)
#define zhpevx FORTRAN_WRAPPER(zhpevx)
#define zhpgst FORTRAN_WRAPPER(zhpgst)
#define zhpgv FORTRAN_WRAPPER(zhpgv)
#define zhpgvd FORTRAN_WRAPPER(zhpgvd)
#define zhpgvx FORTRAN_WRAPPER(zhpgvx)
#define zhprfs FORTRAN_WRAPPER(zhprfs)
#define zhpsv FORTRAN_WRAPPER(zhpsv)
#define zhpsvx FORTRAN_WRAPPER(zhpsvx)
#define zhptrd FORTRAN_WRAPPER(zhptrd)
#define zhptrf FORTRAN_WRAPPER(zhptrf)
#define zhptri FORTRAN_WRAPPER(zhptri)
#define zhptrs FORTRAN_WRAPPER(zhptrs)
#define zhsein FORTRAN_WRAPPER(zhsein)
#define zhseqr FORTRAN_WRAPPER(zhseqr)
#define zlabrd FORTRAN_WRAPPER(zlabrd)
#define zlacgv FORTRAN_WRAPPER(zlacgv)
#define zlacn2 FORTRAN_WRAPPER(zlacn2)
#define zlacon FORTRAN_WRAPPER(zlacon)
#define zlacp2 FORTRAN_WRAPPER(zlacp2)
#define zlacpy FORTRAN_WRAPPER(zlacpy)
#define zlacrm FORTRAN_WRAPPER(zlacrm)
#define zlacrt FORTRAN_WRAPPER(zlacrt)
#define zladiv FORTRAN_WRAPPER(zladiv)
#define zlaed0 FORTRAN_WRAPPER(zlaed0)
#define zlaed7 FORTRAN_WRAPPER(zlaed7)
#define zlaed8 FORTRAN_WRAPPER(zlaed8)
#define zlaein FORTRAN_WRAPPER(zlaein)
#define zlaesy FORTRAN_WRAPPER(zlaesy)
#define zlaev2 FORTRAN_WRAPPER(zlaev2)
#define zlag2c FORTRAN_WRAPPER(zlag2c)
#define zlags2 FORTRAN_WRAPPER(zlags2)
#define zlagtm FORTRAN_WRAPPER(zlagtm)
#define zlahef FORTRAN_WRAPPER(zlahef)
#define zlahqr FORTRAN_WRAPPER(zlahqr)
#define zlahr2 FORTRAN_WRAPPER(zlahr2)
#define zlahrd FORTRAN_WRAPPER(zlahrd)
#define zlaic1 FORTRAN_WRAPPER(zlaic1)
#define zlals0 FORTRAN_WRAPPER(zlals0)
#define zlalsa FORTRAN_WRAPPER(zlalsa)
#define zlalsd FORTRAN_WRAPPER(zlalsd)
#define zlangb FORTRAN_WRAPPER(zlangb)
#define zlange FORTRAN_WRAPPER(zlange)
#define zlangt FORTRAN_WRAPPER(zlangt)
#define zlanhb FORTRAN_WRAPPER(zlanhb)
#define zlanhe FORTRAN_WRAPPER(zlanhe)
#define zlanhf FORTRAN_WRAPPER(zlanhf)
#define zlanhp FORTRAN_WRAPPER(zlanhp)
#define zlanhs FORTRAN_WRAPPER(zlanhs)
#define zlanht FORTRAN_WRAPPER(zlanht)
#define zlansb FORTRAN_WRAPPER(zlansb)
#define zlansp FORTRAN_WRAPPER(zlansp)
#define zlansy FORTRAN_WRAPPER(zlansy)
#define zlantb FORTRAN_WRAPPER(zlantb)
#define zlantp FORTRAN_WRAPPER(zlantp)
#define zlantr FORTRAN_WRAPPER(zlantr)
#define zlapll FORTRAN_WRAPPER(zlapll)
#define zlapmr FORTRAN_WRAPPER(zlapmr)
#define zlapmt FORTRAN_WRAPPER(zlapmt)
#define zlaqgb FORTRAN_WRAPPER(zlaqgb)
#define zlaqge FORTRAN_WRAPPER(zlaqge)
#define zlaqhb FORTRAN_WRAPPER(zlaqhb)
#define zlaqhe FORTRAN_WRAPPER(zlaqhe)
#define zlaqhp FORTRAN_WRAPPER(zlaqhp)
#define zlaqp2 FORTRAN_WRAPPER(zlaqp2)
#define zlaqps FORTRAN_WRAPPER(zlaqps)
#define zlaqr0 FORTRAN_WRAPPER(zlaqr0)
#define zlaqr1 FORTRAN_WRAPPER(zlaqr1)
#define zlaqr2 FORTRAN_WRAPPER(zlaqr2)
#define zlaqr3 FORTRAN_WRAPPER(zlaqr3)
#define zlaqr4 FORTRAN_WRAPPER(zlaqr4)
#define zlaqr5 FORTRAN_WRAPPER(zlaqr5)
#define zlaqsb FORTRAN_WRAPPER(zlaqsb)
#define zlaqsp FORTRAN_WRAPPER(zlaqsp)
#define zlaqsy FORTRAN_WRAPPER(zlaqsy)
#define zlar1v FORTRAN_WRAPPER(zlar1v)
#define zlar2v FORTRAN_WRAPPER(zlar2v)
#define zlarcm FORTRAN_WRAPPER(zlarcm)
#define zlarf FORTRAN_WRAPPER(zlarf)
#define zlarfb FORTRAN_WRAPPER(zlarfb)
#define zlarfg FORTRAN_WRAPPER(zlarfg)
#define zlarfgp FORTRAN_WRAPPER(zlarfgp)
#define zlarft FORTRAN_WRAPPER(zlarft)
#define zlarfx FORTRAN_WRAPPER(zlarfx)
#define zlargv FORTRAN_WRAPPER(zlargv)
#define zlarnv FORTRAN_WRAPPER(zlarnv)
#define zlarrv FORTRAN_WRAPPER(zlarrv)
#define zlartg FORTRAN_WRAPPER(zlartg)
#define zlartv FORTRAN_WRAPPER(zlartv)
#define zlarz FORTRAN_WRAPPER(zlarz)
#define zlarzb FORTRAN_WRAPPER(zlarzb)
#define zlarzt FORTRAN_WRAPPER(zlarzt)
#define zlascl FORTRAN_WRAPPER(zlascl)
#define zlaset FORTRAN_WRAPPER(zlaset)
#define zlasr FORTRAN_WRAPPER(zlasr)
#define zlassq FORTRAN_WRAPPER(zlassq)
#define zlaswp FORTRAN_WRAPPER(zlaswp)
#define zlasyf FORTRAN_WRAPPER(zlasyf)
#define zlat2c FORTRAN_WRAPPER(zlat2c)
#define zlatbs FORTRAN_WRAPPER(zlatbs)
#define zlatdf FORTRAN_WRAPPER(zlatdf)
#define zlatps FORTRAN_WRAPPER(zlatps)
#define zlatrd FORTRAN_WRAPPER(zlatrd)
#define zlatrs FORTRAN_WRAPPER(zlatrs)
#define zlatrz FORTRAN_WRAPPER(zlatrz)
#define zlatzm FORTRAN_WRAPPER(zlatzm)
#define zlauu2 FORTRAN_WRAPPER(zlauu2)
#define zlauum FORTRAN_WRAPPER(zlauum)
#define zpbcon FORTRAN_WRAPPER(zpbcon)
#define zpbequ FORTRAN_WRAPPER(zpbequ)
#define zpbrfs FORTRAN_WRAPPER(zpbrfs)
#define zpbstf FORTRAN_WRAPPER(zpbstf)
#define zpbsv FORTRAN_WRAPPER(zpbsv)
#define zpbsvx FORTRAN_WRAPPER(zpbsvx)
#define zpbtf2 FORTRAN_WRAPPER(zpbtf2)
#define zpbtrf FORTRAN_WRAPPER(zpbtrf)
#define zpbtrs FORTRAN_WRAPPER(zpbtrs)
#define zpftrf FORTRAN_WRAPPER(zpftrf)
#define zpftri FORTRAN_WRAPPER(zpftri)
#define zpftrs FORTRAN_WRAPPER(zpftrs)
#define zpocon FORTRAN_WRAPPER(zpocon)
#define zpoequ FORTRAN_WRAPPER(zpoequ)
#define zpoequb FORTRAN_WRAPPER(zpoequb)
#define zporfs FORTRAN_WRAPPER(zporfs)
#define zposv FORTRAN_WRAPPER(zposv)
#define zposvx FORTRAN_WRAPPER(zposvx)
#define zpotf2 FORTRAN_WRAPPER(zpotf2)
#define zpotrf FORTRAN_WRAPPER(zpotrf)
#define zpotri FORTRAN_WRAPPER(zpotri)
#define zpotrs FORTRAN_WRAPPER(zpotrs)
#define zppcon FORTRAN_WRAPPER(zppcon)
#define zppequ FORTRAN_WRAPPER(zppequ)
#define zpprfs FORTRAN_WRAPPER(zpprfs)
#define zppsv FORTRAN_WRAPPER(zppsv)
#define zppsvx FORTRAN_WRAPPER(zppsvx)
#define zpptrf FORTRAN_WRAPPER(zpptrf)
#define zpptri FORTRAN_WRAPPER(zpptri)
#define zpptrs FORTRAN_WRAPPER(zpptrs)
#define zpstf2 FORTRAN_WRAPPER(zpstf2)
#define zpstrf FORTRAN_WRAPPER(zpstrf)
#define zptcon FORTRAN_WRAPPER(zptcon)
#define zpteqr FORTRAN_WRAPPER(zpteqr)
#define zptrfs FORTRAN_WRAPPER(zptrfs)
#define zptsv FORTRAN_WRAPPER(zptsv)
#define zptsvx FORTRAN_WRAPPER(zptsvx)
#define zpttrf FORTRAN_WRAPPER(zpttrf)
#define zpttrs FORTRAN_WRAPPER(zpttrs)
#define zptts2 FORTRAN_WRAPPER(zptts2)
#define zrot FORTRAN_WRAPPER(zrot)
#define zspcon FORTRAN_WRAPPER(zspcon)
#define zspmv FORTRAN_WRAPPER(zspmv)
#define zspr FORTRAN_WRAPPER(zspr)
#define zsprfs FORTRAN_WRAPPER(zsprfs)
#define zspsv FORTRAN_WRAPPER(zspsv)
#define zspsvx FORTRAN_WRAPPER(zspsvx)
#define zsptrf FORTRAN_WRAPPER(zsptrf)
#define zsptri FORTRAN_WRAPPER(zsptri)
#define zsptrs FORTRAN_WRAPPER(zsptrs)
#define zstedc FORTRAN_WRAPPER(zstedc)
#define zstegr FORTRAN_WRAPPER(zstegr)
#define zstein FORTRAN_WRAPPER(zstein)
#define zstemr FORTRAN_WRAPPER(zstemr)
#define zsteqr FORTRAN_WRAPPER(zsteqr)
#define zsycon FORTRAN_WRAPPER(zsycon)
#define zsyconv FORTRAN_WRAPPER(zsyconv)
#define zsyequb FORTRAN_WRAPPER(zsyequb)
#define zsymv FORTRAN_WRAPPER(zsymv)
#define zsyr FORTRAN_WRAPPER(zsyr)
#define zsyrfs FORTRAN_WRAPPER(zsyrfs)
#define zsysv FORTRAN_WRAPPER(zsysv)
#define zsysvx FORTRAN_WRAPPER(zsysvx)
#define zsyswapr FORTRAN_WRAPPER(zsyswapr)
#define zsytf2 FORTRAN_WRAPPER(zsytf2)
#define zsytrf FORTRAN_WRAPPER(zsytrf)
#define zsytri FORTRAN_WRAPPER(zsytri)
#define zsytri2 FORTRAN_WRAPPER(zsytri2)
#define zsytri2x FORTRAN_WRAPPER(zsytri2x)
#define zsytrs FORTRAN_WRAPPER(zsytrs)
#define zsytrs2 FORTRAN_WRAPPER(zsytrs2)
#define ztbcon FORTRAN_WRAPPER(ztbcon)
#define ztbrfs FORTRAN_WRAPPER(ztbrfs)
#define ztbtrs FORTRAN_WRAPPER(ztbtrs)
#define ztfsm FORTRAN_WRAPPER(ztfsm)
#define ztftri FORTRAN_WRAPPER(ztftri)
#define ztfttp FORTRAN_WRAPPER(ztfttp)
#define ztfttr FORTRAN_WRAPPER(ztfttr)
#define ztgevc FORTRAN_WRAPPER(ztgevc)
#define ztgex2 FORTRAN_WRAPPER(ztgex2)
#define ztgexc FORTRAN_WRAPPER(ztgexc)
#define ztgsen FORTRAN_WRAPPER(ztgsen)
#define ztgsja FORTRAN_WRAPPER(ztgsja)
#define ztgsna FORTRAN_WRAPPER(ztgsna)
#define ztgsy2 FORTRAN_WRAPPER(ztgsy2)
#define ztgsyl FORTRAN_WRAPPER(ztgsyl)
#define ztpcon FORTRAN_WRAPPER(ztpcon)
#define ztpmqrt FORTRAN_WRAPPER(ztpmqrt)
#define ztpqrt FORTRAN_WRAPPER(ztpqrt)
#define ztpqrt2 FORTRAN_WRAPPER(ztpqrt2)
#define ztprfb FORTRAN_WRAPPER(ztprfb)
#define ztprfs FORTRAN_WRAPPER(ztprfs)
#define ztptri FORTRAN_WRAPPER(ztptri)
#define ztptrs FORTRAN_WRAPPER(ztptrs)
#define ztpttf FORTRAN_WRAPPER(ztpttf)
#define ztpttr FORTRAN_WRAPPER(ztpttr)
#define ztrcon FORTRAN_WRAPPER(ztrcon)
#define ztrevc FORTRAN_WRAPPER(ztrevc)
#define ztrexc FORTRAN_WRAPPER(ztrexc)
#define ztrrfs FORTRAN_WRAPPER(ztrrfs)
#define ztrsen FORTRAN_WRAPPER(ztrsen)
#define ztrsna FORTRAN_WRAPPER(ztrsna)
#define ztrsyl FORTRAN_WRAPPER(ztrsyl)
#define ztrti2 FORTRAN_WRAPPER(ztrti2)
#define ztrtri FORTRAN_WRAPPER(ztrtri)
#define ztrtrs FORTRAN_WRAPPER(ztrtrs)
#define ztrttf FORTRAN_WRAPPER(ztrttf)
#define ztrttp FORTRAN_WRAPPER(ztrttp)
#define ztzrqf FORTRAN_WRAPPER(ztzrqf)
#define ztzrzf FORTRAN_WRAPPER(ztzrzf)
#define zunbdb FORTRAN_WRAPPER(zunbdb)
#define zuncsd FORTRAN_WRAPPER(zuncsd)
#define zung2l FORTRAN_WRAPPER(zung2l)
#define zung2r FORTRAN_WRAPPER(zung2r)
#define zungbr FORTRAN_WRAPPER(zungbr)
#define zunghr FORTRAN_WRAPPER(zunghr)
#define zungl2 FORTRAN_WRAPPER(zungl2)
#define zunglq FORTRAN_WRAPPER(zunglq)
#define zungql FORTRAN_WRAPPER(zungql)
#define zungqr FORTRAN_WRAPPER(zungqr)
#define zungr2 FORTRAN_WRAPPER(zungr2)
#define zungrq FORTRAN_WRAPPER(zungrq)
#define zungtr FORTRAN_WRAPPER(zungtr)
#define zunm2l FORTRAN_WRAPPER(zunm2l)
#define zunm2r FORTRAN_WRAPPER(zunm2r)
#define zunmbr FORTRAN_WRAPPER(zunmbr)
#define zunmhr FORTRAN_WRAPPER(zunmhr)
#define zunml2 FORTRAN_WRAPPER(zunml2)
#define zunmlq FORTRAN_WRAPPER(zunmlq)
#define zunmql FORTRAN_WRAPPER(zunmql)
#define zunmqr FORTRAN_WRAPPER(zunmqr)
#define zunmr2 FORTRAN_WRAPPER(zunmr2)
#define zunmr3 FORTRAN_WRAPPER(zunmr3)
#define zunmrq FORTRAN_WRAPPER(zunmrq)
#define zunmrz FORTRAN_WRAPPER(zunmrz)
#define zunmtr FORTRAN_WRAPPER(zunmtr)
#define zupgtr FORTRAN_WRAPPER(zupgtr)
#define zupmtr FORTRAN_WRAPPER(zupmtr)

#endif

extern "C"
{
char chla_transtype(uint32_t *trans);
void dbbcsd(const char * const jobu1, const char * const jobu2, const char * const jobv1t, const char * const jobv2t, const char * const trans, uint32_t *m, uint32_t *p, uint32_t *q, double *theta, double *phi, double *u1, const uint32_t * const ldu1, double *u2, const uint32_t * const ldu2, double *v1t, const uint32_t * const ldv1t, double *v2t, const uint32_t * const ldv2t, double *b11d, double *b11e, double *b12d, double *b12e, double *b21d, double *b21e, double *b22d, double *b22e, double *work, const uint32_t * const lwork, int *info);
void dbdsdc(const char * const uplo, const char * const compq, const uint32_t * const n, double *d, double *e, double *u, const uint32_t * const ldu, double *vt, const uint32_t * const ldvt, double *q, uint32_t *iq, double *work, int *iwork, int *info);
void dbdsqr(const char * const uplo, const uint32_t * const n, uint32_t *ncvt, uint32_t *nru, uint32_t *ncc, double *d, double *e, double *vt, const uint32_t * const ldvt, double *u, const uint32_t * const ldu, double *c, const uint32_t * const ldc, double *work, int *info);
void ddisna(const char * const job, uint32_t *m, const uint32_t * const n, double *d, double *sep, int *info);
void dgbbrd(const char * const vect, uint32_t *m, const uint32_t * const n, uint32_t *ncc, uint32_t *kl, uint32_t *ku, double *ab, const uint32_t * const ldab, double *d, double *e, double *q, const uint32_t * const ldq, double *pt, const uint32_t * const ldpt, double *c, const uint32_t * const ldc, double *work, int *info);
void dgbcon(const char * const norm, const uint32_t * const n, uint32_t *kl, uint32_t *ku, double *ab, const uint32_t * const ldab, int *ipiv, double *anorm, double *rcond, double *work, int *iwork, int *info);
void dgbequ(uint32_t *m, const uint32_t * const n, uint32_t *kl, uint32_t *ku, double *ab, const uint32_t * const ldab, double *r, double *c, double *rowcnd, double *colcnd, double *amax, int *info);
void dgbequb(uint32_t *m, const uint32_t * const n, uint32_t *kl, uint32_t *ku, double *ab, const uint32_t * const ldab, double *r, double *c, double *rowcnd, double *colcnd, double *amax, int *info);
void dgbrfs(const char * const trans, const uint32_t * const n, uint32_t *kl, uint32_t *ku, uint32_t *nrhs, double *ab, const uint32_t * const ldab, double *afb, const uint32_t * const ldafb, int *ipiv, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *ferr, double *berr, double *work, int *iwork, int *info);
void dgbsv(const uint32_t * const n, uint32_t *kl, uint32_t *ku, uint32_t *nrhs, double *ab, const uint32_t * const ldab, int *ipiv, double *b, const uint32_t * const ldb, int *info);
void dgbsvx(const char * const fact, const char * const trans, const uint32_t * const n, uint32_t *kl, uint32_t *ku, uint32_t *nrhs, double *ab, const uint32_t * const ldab, double *afb, const uint32_t * const ldafb, int *ipiv, const char * const equed, double *r, double *c, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);
void dgbtf2(uint32_t *m, const uint32_t * const n, uint32_t *kl, uint32_t *ku, double *ab, const uint32_t * const ldab, int *ipiv, int *info);
void dgbtrf(uint32_t *m, const uint32_t * const n, uint32_t *kl, uint32_t *ku, double *ab, const uint32_t * const ldab, int *ipiv, int *info);
void dgbtrs(const char * const trans, const uint32_t * const n, uint32_t *kl, uint32_t *ku, uint32_t *nrhs, double *ab, const uint32_t * const ldab, int *ipiv, double *b, const uint32_t * const ldb, int *info);
void dgebak(const char * const job, const char * const side, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, double *scale, uint32_t *m, double *v, const uint32_t * const ldv, int *info);
void dgebal(const char * const job, const uint32_t * const n, double *a, const uint32_t * const lda, uint32_t *ilo, uint32_t *ihi, double *scale, int *info);
void dgebd2(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *d, double *e, double *tauq, double *taup, double *work, int *info);
void dgebrd(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *d, double *e, double *tauq, double *taup, double *work, const uint32_t * const lwork, int *info);
void dgecon(const char * const norm, const uint32_t * const n, double *a, const uint32_t * const lda, double *anorm, double *rcond, double *work, int *iwork, int *info);
void dgeequ(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax, int *info);
void dgeequb(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax, int *info);
void dgees(const char * const jobvs, const char * const sort, uint32_t (*select)(), const uint32_t * const n, double *a, const uint32_t * const lda, uint32_t *sdim, double *wr, double *wi, double *vs, const uint32_t * const ldvs, double *work, const uint32_t * const lwork, uint32_t *bwork, int *info);
void dgeesx(const char * const jobvs, const char * const sort, uint32_t (*select)(), const char * const sense, const uint32_t * const n, double *a, const uint32_t * const lda, uint32_t *sdim, double *wr, double *wi, double *vs, const uint32_t * const ldvs, double *rconde, double *rcondv, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, uint32_t *bwork, int *info);
void dgeev(const char * const jobvl, const char * const jobvr, const uint32_t * const n, double *a, const uint32_t * const lda, double *wr, double *wi, double *vl, const uint32_t * const ldvl, double *vr, const uint32_t * const ldvr, double *work, const uint32_t * const lwork, int *info);
void dgeevx(const char * const balanc, const char * const jobvl, const char * const jobvr, const char * const sense, const uint32_t * const n, double *a, const uint32_t * const lda, double *wr, double *wi, double *vl, const uint32_t * const ldvl, double *vr, const uint32_t * const ldvr, uint32_t *ilo, uint32_t *ihi, double *scale, double *abnrm, double *rconde, double *rcondv, double *work, const uint32_t * const lwork, int *iwork, int *info);
void dgegs(const char * const jobvsl, const char * const jobvsr, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *alphar, double *alphai, double *beta, double *vsl, const uint32_t * const ldvsl, double *vsr, const uint32_t * const ldvsr, double *work, const uint32_t * const lwork, int *info);
void dgegv(const char * const jobvl, const char * const jobvr, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *alphar, double *alphai, double *beta, double *vl, const uint32_t * const ldvl, double *vr, const uint32_t * const ldvr, double *work, const uint32_t * const lwork, int *info);
void dgehd2(const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, double *a, const uint32_t * const lda, double *tau, double *work, int *info);
void dgehrd(const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, double *a, const uint32_t * const lda, double *tau, double *work, const uint32_t * const lwork, int *info);
void dgejsv(const char * const joba, const char * const jobu, const char * const jobv, const char * const jobr, const char * const jobt, const char * const jobp, uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *sva, double *u, const uint32_t * const ldu, double *v, const uint32_t * const ldv, double *work, const uint32_t * const lwork, int *iwork, int *info);
void dgelq2(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *tau, double *work, int *info);
void dgelqf(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *tau, double *work, const uint32_t * const lwork, int *info);
void dgels(const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *work, const uint32_t * const lwork, int *info);
void dgelsd(uint32_t *m, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *s, double *rcond, uint32_t *rank, double *work, const uint32_t * const lwork, int *iwork, int *info);
void dgelss(uint32_t *m, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *s, double *rcond, uint32_t *rank, double *work, const uint32_t * const lwork, int *info);
void dgelsx(uint32_t *m, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, uint32_t *jpvt, double *rcond, uint32_t *rank, double *work, int *info);
void dgelsy(uint32_t *m, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, uint32_t *jpvt, double *rcond, uint32_t *rank, double *work, const uint32_t * const lwork, int *info);
void dgemqrt(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, uint32_t *nb, double *v, const uint32_t * const ldv, double *t, const uint32_t * const ldt, double *c, const uint32_t * const ldc, double *work, int *info);
void dgeql2(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *tau, double *work, int *info);
void dgeqlf(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *tau, double *work, const uint32_t * const lwork, int *info);
void dgeqp3(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, uint32_t *jpvt, double *tau, double *work, const uint32_t * const lwork, int *info);
void dgeqpf(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, uint32_t *jpvt, double *tau, double *work, int *info);
void dgeqr2(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *tau, double *work, int *info);
void dgeqr2p(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *tau, double *work, int *info);
void dgeqrf(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *tau, double *work, const uint32_t * const lwork, int *info);
void dgeqrfp(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *tau, double *work, const uint32_t * const lwork, int *info);
void dgeqrt(uint32_t *m, const uint32_t * const n, uint32_t *nb, double *a, const uint32_t * const lda, double *t, const uint32_t * const ldt, double *work, int *info);
void dgeqrt2(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *t, const uint32_t * const ldt, int *info);
void dgeqrt3(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *t, const uint32_t * const ldt, int *info);
void dgerfs(const char * const trans, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *af, const uint32_t * const ldaf, int *ipiv, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *ferr, double *berr, double *work, int *iwork, int *info);
void dgerq2(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *tau, double *work, int *info);
void dgerqf(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *tau, double *work, const uint32_t * const lwork, int *info);
void dgesc2(const uint32_t * const n, double *a, const uint32_t * const lda, double *rhs, int *ipiv, int *jpiv, double *scale);
void dgesdd(const char * const jobz, uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *s, double *u, const uint32_t * const ldu, double *vt, const uint32_t * const ldvt, double *work, const uint32_t * const lwork, int *iwork, int *info);
void dgesv(const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, int *ipiv, double *b, const uint32_t * const ldb, int *info);
void dgesvd(const char * const jobu, const char * const jobvt, uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *s, double *u, const uint32_t * const ldu, double *vt, const uint32_t * const ldvt, double *work, const uint32_t * const lwork, int *info);
void dgesvj(const char * const joba, const char * const jobu, const char * const jobv, uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *sva, uint32_t *mv, double *v, const uint32_t * const ldv, double *work, const uint32_t * const lwork, int *info);
void dgesvx(const char * const fact, const char * const trans, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *af, const uint32_t * const ldaf, int *ipiv, const char * const equed, double *r, double *c, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);
void dgetc2(const uint32_t * const n, double *a, const uint32_t * const lda, int *ipiv, int *jpiv, int *info);
void dgetf2(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, int *ipiv, int *info);
void dgetrf(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, int *ipiv, int *info);
void dgetri(const uint32_t * const n, double *a, const uint32_t * const lda, int *ipiv, double *work, const uint32_t * const lwork, int *info);
void dgetrs(const char * const trans, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, int *ipiv, double *b, const uint32_t * const ldb, int *info);
void dggbak(const char * const job, const char * const side, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, double *lscale, double *rscale, uint32_t *m, double *v, const uint32_t * const ldv, int *info);
void dggbal(const char * const job, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, uint32_t *ilo, uint32_t *ihi, double *lscale, double *rscale, double *work, int *info);
void dgges(const char * const jobvsl, const char * const jobvsr, const char * const sort, uint32_t (*selctg)(), const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, uint32_t *sdim, double *alphar, double *alphai, double *beta, double *vsl, const uint32_t * const ldvsl, double *vsr, const uint32_t * const ldvsr, double *work, const uint32_t * const lwork, uint32_t *bwork, int *info);
void dggesx(const char * const jobvsl, const char * const jobvsr, const char * const sort, uint32_t (*selctg)(), const char * const sense, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, uint32_t *sdim, double *alphar, double *alphai, double *beta, double *vsl, const uint32_t * const ldvsl, double *vsr, const uint32_t * const ldvsr, double *rconde, double *rcondv, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, uint32_t *bwork, int *info);
void dggev(const char * const jobvl, const char * const jobvr, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *alphar, double *alphai, double *beta, double *vl, const uint32_t * const ldvl, double *vr, const uint32_t * const ldvr, double *work, const uint32_t * const lwork, int *info);
void dggevx(const char * const balanc, const char * const jobvl, const char * const jobvr, const char * const sense, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *alphar, double *alphai, double *beta, double *vl, const uint32_t * const ldvl, double *vr, const uint32_t * const ldvr, uint32_t *ilo, uint32_t *ihi, double *lscale, double *rscale, double *abnrm, double *bbnrm, double *rconde, double *rcondv, double *work, const uint32_t * const lwork, int *iwork, uint32_t *bwork, int *info);
void dggglm(const uint32_t * const n, uint32_t *m, uint32_t *p, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *d, double *x, double *y, double *work, const uint32_t * const lwork, int *info);
void dgghrd(const char * const compq, const char * const compz, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *q, const uint32_t * const ldq, double *z, const uint32_t * const ldz, int *info);
void dgglse(uint32_t *m, const uint32_t * const n, uint32_t *p, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *c, double *d, double *x, double *work, const uint32_t * const lwork, uint32_t *info);
void dggqrf(const uint32_t * const n, uint32_t *m, uint32_t *p, double *a, const uint32_t * const lda, double *taua, double *b, const uint32_t * const ldb, double *taub, double *work, const uint32_t * const lwork, int *info);
void dggrqf(uint32_t *m, uint32_t *p, const uint32_t * const n, double *a, const uint32_t * const lda, double *taua, double *b, const uint32_t * const ldb, double *taub, double *work, const uint32_t * const lwork, int *info);
void dggsvd(const char * const jobu, const char * const jobv, const char * const jobq, uint32_t *m, const uint32_t * const n, uint32_t *p, uint32_t *k, uint32_t *l, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *alpha, double *beta, double *u, const uint32_t * const ldu, double *v, const uint32_t * const ldv, double *q, const uint32_t * const ldq, double *work, int *iwork, int *info);
void dggsvp(const char * const jobu, const char * const jobv, const char * const jobq, uint32_t *m, uint32_t *p, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *tola, double *tolb, uint32_t *k, uint32_t *l, double *u, const uint32_t * const ldu, double *v, const uint32_t * const ldv, double *q, const uint32_t * const ldq, int *iwork, double *tau, double *work, int *info);
void dgsvj0(const char * const jobv, uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *d, double *sva, uint32_t *mv, double *v, const uint32_t * const ldv, double *eps, double *sfmin, double *tol, uint32_t *nsweep, double *work, const uint32_t * const lwork, int *info);
void dgsvj1(const char * const jobv, uint32_t *m, const uint32_t * const n, uint32_t *n1, double *a, const uint32_t * const lda, double *d, double *sva, uint32_t *mv, double *v, const uint32_t * const ldv, double *eps, double *sfmin, double *tol, uint32_t *nsweep, double *work, const uint32_t * const lwork, int *info);
void dgtcon(const char * const norm, const uint32_t * const n, double *dl, double *d, double *du, double *du2, int *ipiv, double *anorm, double *rcond, double *work, int *iwork, int *info);
void dgtrfs(const char * const trans, const uint32_t * const n, uint32_t *nrhs, double *dl, double *d, double *du, double *dlf, double *df, double *duf, double *du2, int *ipiv, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *ferr, double *berr, double *work, int *iwork, int *info);
void dgtsv(const uint32_t * const n, uint32_t *nrhs, double *dl, double *d, double *du, double *b, const uint32_t * const ldb, int *info);
void dgtsvx(const char * const fact, const char * const trans, const uint32_t * const n, uint32_t *nrhs, double *dl, double *d, double *du, double *dlf, double *df, double *duf, double *du2, int *ipiv, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);
void dgttrf(const uint32_t * const n, double *dl, double *d, double *du, double *du2, int *ipiv, int *info);
void dgttrs(const char * const trans, const uint32_t * const n, uint32_t *nrhs, double *dl, double *d, double *du, double *du2, int *ipiv, double *b, const uint32_t * const ldb, int *info);
void dgtts2(uint32_t *itrans, const uint32_t * const n, uint32_t *nrhs, double *dl, double *d, double *du, double *du2, int *ipiv, double *b, const uint32_t * const ldb);
void dhgeqz(const char * const job, const char * const compq, const char * const compz, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, double *h, const uint32_t * const ldh, double *t, const uint32_t * const ldt, double *alphar, double *alphai, double *beta, double *q, const uint32_t * const ldq, double *z, const uint32_t * const ldz, double *work, const uint32_t * const lwork, int *info);
void dhsein(const char * const side, const char * const eigsrc, const char * const initv, uint32_t *select, const uint32_t * const n, double *h, const uint32_t * const ldh, double *wr, double *wi, double *vl, const uint32_t * const ldvl, double *vr, const uint32_t * const ldvr, uint32_t *mm, uint32_t *m, double *work, uint32_t *ifaill, uint32_t *ifailr, int *info);
void dhseqr(const char * const job, const char * const compz, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, double *h, const uint32_t * const ldh, double *wr, double *wi, double *z, const uint32_t * const ldz, double *work, const uint32_t * const lwork, int *info);
uint32_t disnan(double *din);
void dlabad(double *small, double *large);
void dlabrd(uint32_t *m, const uint32_t * const n, uint32_t *nb, double *a, const uint32_t * const lda, double *d, double *e, double *tauq, double *taup, double *x, const uint32_t * const ldx, double *y, const uint32_t * const ldy);
void dlacn2(const uint32_t * const n, double *v, double *x, uint32_t *isgn, double *est, uint32_t *kase, uint32_t *isave);
void dlacon(const uint32_t * const n, double *v, double *x, uint32_t *isgn, double *est, uint32_t *kase);
void dlacpy(const char * const uplo, uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb);
void dladiv(double *a, double *b, double *c, double *d, double *p, double *q);
void dlae2(double *a, double *b, double *c, double *rt1, double *rt2);
void dlaebz(uint32_t *ijob, uint32_t *nitmax, const uint32_t * const n, uint32_t *mmax, uint32_t *minp, uint32_t *nbmin, double *abstol, double *reltol, double *pivmin, double *d, double *e, double *e2, uint32_t *nval, double *ab, double *c, uint32_t *mout, uint32_t *nab, double *work, int *iwork, int *info);
void dlaed0(uint32_t *icompq, uint32_t *qsiz, const uint32_t * const n, double *d, double *e, double *q, const uint32_t * const ldq, double *qstore, const uint32_t * const ldqs, double *work, int *iwork, int *info);
void dlaed1(const uint32_t * const n, double *d, double *q, const uint32_t * const ldq, uint32_t *indxq, double *rho, uint32_t *cutpnt, double *work, int *iwork, int *info);
void dlaed2(uint32_t *k, const uint32_t * const n, uint32_t *n1, double *d, double *q, const uint32_t * const ldq, uint32_t *indxq, double *rho, double *z, double *dlamda, double *w, double *q2, uint32_t *indx, uint32_t *indxc, uint32_t *indxp, uint32_t *coltyp, int *info);
void dlaed3(uint32_t *k, const uint32_t * const n, uint32_t *n1, double *d, double *q, const uint32_t * const ldq, double *rho, double *dlamda, double *q2, uint32_t *indx, uint32_t *ctot, double *w, double *s, int *info);
void dlaed4(const uint32_t * const n, uint32_t *i, double *d, double *z, double *delta, double *rho, double *dlam, int *info);
void dlaed5(uint32_t *i, double *d, double *z, double *delta, double *rho, double *dlam);
void dlaed6(uint32_t *kniter, uint32_t *orgati, double *rho, double *d, double *z, double *finit, double *tau, int *info);
void dlaed7(uint32_t *icompq, const uint32_t * const n, uint32_t *qsiz, uint32_t *tlvls, uint32_t *curlvl, uint32_t *curpbm, double *d, double *q, const uint32_t * const ldq, uint32_t *indxq, double *rho, uint32_t *cutpnt, double *qstore, uint32_t *qptr, uint32_t *prmptr, uint32_t *perm, uint32_t *givptr, uint32_t *givcol, double *givnum, double *work, int *iwork, int *info);
void dlaed8(uint32_t *icompq, uint32_t *k, const uint32_t * const n, uint32_t *qsiz, double *d, double *q, const uint32_t * const ldq, uint32_t *indxq, double *rho, uint32_t *cutpnt, double *z, double *dlamda, double *q2, const uint32_t * const ldq2, double *w, uint32_t *perm, uint32_t *givptr, uint32_t *givcol, double *givnum, uint32_t *indxp, uint32_t *indx, int *info);
void dlaed9(uint32_t *k, uint32_t *kstart, uint32_t *kstop, const uint32_t * const n, double *d, double *q, const uint32_t * const ldq, double *rho, double *dlamda, double *w, double *s, const uint32_t * const lds, int *info);
void dlaeda(const uint32_t * const n, uint32_t *tlvls, uint32_t *curlvl, uint32_t *curpbm, uint32_t *prmptr, uint32_t *perm, uint32_t *givptr, uint32_t *givcol, double *givnum, double *q, uint32_t *qptr, double *z, double *ztemp, int *info);
void dlaein(uint32_t *rightv, uint32_t *noinit, const uint32_t * const n, double *h, const uint32_t * const ldh, double *wr, double *wi, double *vr, double *vi, double *b, const uint32_t * const ldb, double *work, double *eps3, double *smlnum, double *bignum, int *info);
void dlaev2(double *a, double *b, double *c, double *rt1, double *rt2, double *cs1, double *sn1);
void dlaexc(uint32_t *wantq, const uint32_t * const n, double *t, const uint32_t * const ldt, double *q, const uint32_t * const ldq, uint32_t *j1_, uint32_t *n1, uint32_t *n2, double *work, int *info);
void dlag2(double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *safmin, double *scale1, double *scale2, double *wr1, double *wr2, double *wi);
void dlags2(uint32_t *upper, double *a1, double *a2, double *a3, double *b1, double *b2, double *b3, double *csu, double *snu, double *csv, double *snv, double *csq, double *snq);
void dlagtf(const uint32_t * const n, double *a, double *lambda, double *b, double *c, double *tol, double *d, uint32_t *in, int *info);
void dlagtm(const char * const trans, const uint32_t * const n, uint32_t *nrhs, double *alpha, double *dl, double *d, double *du, double *x, const uint32_t * const ldx, double *beta, double *b, const uint32_t * const ldb);
void dlagts(uint32_t *job, const uint32_t * const n, double *a, double *b, double *c, double *d, uint32_t *in, double *y, double *tol, int *info);
void dlagv2(double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *alphar, double *alphai, double *beta, double *csl, double *snl, double *csr, double *snr);
void dlahqr(uint32_t *wantt, uint32_t *wantz, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, double *h, const uint32_t * const ldh, double *wr, double *wi, uint32_t *iloz, uint32_t *ihiz, double *z, const uint32_t * const ldz, int *info);
void dlahr2(const uint32_t * const n, uint32_t *k, uint32_t *nb, double *a, const uint32_t * const lda, double *tau, double *t, const uint32_t * const ldt, double *y, const uint32_t * const ldy);
void dlahrd(const uint32_t * const n, uint32_t *k, uint32_t *nb, double *a, const uint32_t * const lda, double *tau, double *t, const uint32_t * const ldt, double *y, const uint32_t * const ldy);
void dlaic1(uint32_t *job, uint32_t *j, double *x, double *sest, double *w, double *gamma_, double *sestpr, double *s, double *c);
uint32_t dlaisnan(double *din1, double *din2);
void dlaln2(uint32_t *ltrans, uint32_t *na, uint32_t *nw, double *smin, double *ca, double *a, const uint32_t * const lda, double *d1, double *d2, double *b, const uint32_t * const ldb, double *wr, double *wi, double *x, const uint32_t * const ldx, double *scale, double *xnorm, int *info);
void dlals0(uint32_t *icompq, uint32_t *nl, uint32_t *nr, uint32_t *sqre, uint32_t *nrhs, double *b, const uint32_t * const ldb, double *bx, const uint32_t * const ldbx, uint32_t *perm, uint32_t *givptr, uint32_t *givcol, const uint32_t * const ldgcol, double *givnum, const uint32_t * const ldgnum, double *poles, double *difl, double *difr, double *z, uint32_t *k, double *c, double *s, double *work, int *info);
void dlalsa(uint32_t *icompq, uint32_t *smlsiz, const uint32_t * const n, uint32_t *nrhs, double *b, const uint32_t * const ldb, double *bx, const uint32_t * const ldbx, double *u, const uint32_t * const ldu, double *vt, uint32_t *k, double *difl, double *difr, double *z, double *poles, uint32_t *givptr, uint32_t *givcol, const uint32_t * const ldgcol, uint32_t *perm, double *givnum, double *c, double *s, double *work, int *iwork, int *info);
void dlalsd(const char * const uplo, uint32_t *smlsiz, const uint32_t * const n, uint32_t *nrhs, double *d, double *e, double *b, const uint32_t * const ldb, double *rcond, uint32_t *rank, double *work, int *iwork, int *info);
void dlamrg(uint32_t *n1, uint32_t *n2, double *a, uint32_t *dtrd1, uint32_t *dtrd2, uint32_t *index);
uint32_t dlaneg(const uint32_t * const n, double *d, double *lld, double *sigma, double *pivmin, uint32_t *r);
double dlangb(const char * const norm, const uint32_t * const n, uint32_t *kl, uint32_t *ku, double *ab, const uint32_t * const ldab, double *work);
double dlange(const char * const norm, uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *work);
double dlangt(const char * const norm, const uint32_t * const n, double *dl, double *d, double *du);
double dlanhs(const char * const norm, const uint32_t * const n, double *a, const uint32_t * const lda, double *work);
double dlansb(const char * const norm, const char * const uplo, const uint32_t * const n, uint32_t *k, double *ab, const uint32_t * const ldab, double *work);
double dlansf(const char * const norm, const char * const transr, const char * const uplo, const uint32_t * const n, double *a, double *work);
double dlansp(const char * const norm, const char * const uplo, const uint32_t * const n, double *ap, double *work);
double dlanst(const char * const norm, const uint32_t * const n, double *d, double *e);
double dlansy(const char * const norm, const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *work);
double dlantb(const char * const norm, const char * const uplo, const char * const diag, const uint32_t * const n, uint32_t *k, double *ab, const uint32_t * const ldab, double *work);
double dlantp(const char * const norm, const char * const uplo, const char * const diag, const uint32_t * const n, double *ap, double *work);
double dlantr(const char * const norm, const char * const uplo, const char * const diag, uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *work);
void dlanv2(double *a, double *b, double *c, double *d, double *rt1r, double *rt1i, double *rt2r, double *rt2i, double *cs, double *sn);
void dlapll(const uint32_t * const n, double *x, const uint32_t * const incx, double *y, const uint32_t * const incy, double *ssmin);
void dlapmr(uint32_t *forwrd, uint32_t *m, const uint32_t * const n, double *x, const uint32_t * const ldx, uint32_t *k);
void dlapmt(uint32_t *forwrd, uint32_t *m, const uint32_t * const n, double *x, const uint32_t * const ldx, uint32_t *k);
double dlapy2(double *x, double *y);
double dlapy3(double *x, double *y, double *z);
void dlaqgb(uint32_t *m, const uint32_t * const n, uint32_t *kl, uint32_t *ku, double *ab, const uint32_t * const ldab, double *r, double *c, double *rowcnd, double *colcnd, double *amax, const char * const equed);
void dlaqge(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax, const char * const equed);
void dlaqp2(uint32_t *m, const uint32_t * const n, uint32_t *offset, double *a, const uint32_t * const lda, uint32_t *jpvt, double *tau, double *vn1, double *vn2, double *work);
void dlaqps(uint32_t *m, const uint32_t * const n, uint32_t *offset, uint32_t *nb, uint32_t *kb, double *a, const uint32_t * const lda, uint32_t *jpvt, double *tau, double *vn1, double *vn2, double *auxv, double *f, const uint32_t * const ldf);
void dlaqr0(uint32_t *wantt, uint32_t *wantz, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, double *h, const uint32_t * const ldh, double *wr, double *wi, uint32_t *iloz, uint32_t *ihiz, double *z, const uint32_t * const ldz, double *work, const uint32_t * const lwork, int *info);
void dlaqr1(const uint32_t * const n, double *h, const uint32_t * const ldh, double *sr1, double *si1, double *sr2, double *si2, double *v);
void dlaqr2(uint32_t *wantt, uint32_t *wantz, const uint32_t * const n, uint32_t *ktop, uint32_t *kbot, uint32_t *nw, double *h, const uint32_t * const ldh, uint32_t *iloz, uint32_t *ihiz, double *z, const uint32_t * const ldz, uint32_t *ns, uint32_t *nd, double *sr, double *si, double *v, const uint32_t * const ldv, uint32_t *nh, double *t, const uint32_t * const ldt, uint32_t *nv, double *wv, const uint32_t * const ldwv, double *work, const uint32_t * const lwork);
void dlaqr3(uint32_t *wantt, uint32_t *wantz, const uint32_t * const n, uint32_t *ktop, uint32_t *kbot, uint32_t *nw, double *h, const uint32_t * const ldh, uint32_t *iloz, uint32_t *ihiz, double *z, const uint32_t * const ldz, uint32_t *ns, uint32_t *nd, double *sr, double *si, double *v, const uint32_t * const ldv, uint32_t *nh, double *t, const uint32_t * const ldt, uint32_t *nv, double *wv, const uint32_t * const ldwv, double *work, const uint32_t * const lwork);
void dlaqr4(uint32_t *wantt, uint32_t *wantz, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, double *h, const uint32_t * const ldh, double *wr, double *wi, uint32_t *iloz, uint32_t *ihiz, double *z, const uint32_t * const ldz, double *work, const uint32_t * const lwork, int *info);
void dlaqr5(uint32_t *wantt, uint32_t *wantz, uint32_t *kacc22, const uint32_t * const n, uint32_t *ktop, uint32_t *kbot, uint32_t *nshfts, double *sr, double *si, double *h, const uint32_t * const ldh, uint32_t *iloz, uint32_t *ihiz, double *z, const uint32_t * const ldz, double *v, const uint32_t * const ldv, double *u, const uint32_t * const ldu, uint32_t *nv, double *wv, const uint32_t * const ldwv, uint32_t *nh, double *wh, const uint32_t * const ldwh);
void dlaqsb(const char * const uplo, const uint32_t * const n, uint32_t *kd, double *ab, const uint32_t * const ldab, double *s, double *scond, double *amax, const char * const equed);
void dlaqsp(const char * const uplo, const uint32_t * const n, double *ap, double *s, double *scond, double *amax, const char * const equed);
void dlaqsy(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *s, double *scond, double *amax, const char * const equed);
void dlaqtr(uint32_t *ltran, uint32_t *lreal, const uint32_t * const n, double *t, const uint32_t * const ldt, double *b, double *w, double *scale, double *x, double *work, int *info);
void dlar1v(const uint32_t * const n, uint32_t *b1, uint32_t *bn, double *lambda, double *d, double *l, double *ld, double *lld, double *pivmin, double *gaptol, double *z, uint32_t *wantnc, uint32_t *negcnt, double *ztz, double *mingma, uint32_t *r, uint32_t *isuppz, double *nrminv, double *resid, double *rqcorr, double *work);
void dlar2v(const uint32_t * const n, double *x, double *y, double *z, const uint32_t * const incx, double *c, double *s, const uint32_t * const incc);
void dlarf(const char * const side, uint32_t *m, const uint32_t * const n, double *v, const uint32_t * const incv, double *tau, double *c, const uint32_t * const ldc, double *work);
void dlarfb(const char * const side, const char * const trans, const char * const direct, const char * const storev, uint32_t *m, const uint32_t * const n, uint32_t *k, double *v, const uint32_t * const ldv, double *t, const uint32_t * const ldt, double *c, const uint32_t * const ldc, double *work, const uint32_t * const ldwork);
void dlarfg(const uint32_t * const n, double *alpha, double *x, const uint32_t * const incx, double *tau);
void dlarfgp(const uint32_t * const n, double *alpha, double *x, const uint32_t * const incx, double *tau);
void dlarft(const char * const direct, const char * const storev, const uint32_t * const n, uint32_t *k, double *v, const uint32_t * const ldv, double *tau, double *t, const uint32_t * const ldt);
void dlarfx(const char * const side, uint32_t *m, const uint32_t * const n, double *v, double *tau, double *c, const uint32_t * const ldc, double *work);
void dlargv(const uint32_t * const n, double *x, const uint32_t * const incx, double *y, const uint32_t * const incy, double *c, const uint32_t * const incc);
void dlarnv(uint32_t *idist, uint32_t *iseed, const uint32_t * const n, double *x);
void dlarra(const uint32_t * const n, double *d, double *e, double *e2, double *spltol, double *tnrm, uint32_t *nsplit, uint32_t *isplit, int *info);
void dlarrb(const uint32_t * const n, double *d, double *lld, uint32_t *ifirst, uint32_t *ilast, double *rtol1, double *rtol2, uint32_t *offset, double *w, double *wgap, double *werr, double *work, int *iwork, double *pivmin, double *spdiam, uint32_t *twist, int *info);
void dlarrc(const char * const jobt, const uint32_t * const n, double *vl, double *vu, double *d, double *e, double *pivmin, uint32_t *eigcnt, uint32_t *lcnt, uint32_t *rcnt, int *info);
void dlarrd(const char * const range, const char * const order, const uint32_t * const n, double *vl, double *vu, uint32_t *il, uint32_t *iu, double *gers, double *reltol, double *d, double *e, double *e2, double *pivmin, uint32_t *nsplit, uint32_t *isplit, uint32_t *m, double *w, double *werr, double *wl, double *wu, uint32_t *iblock, uint32_t *indexw, double *work, int *iwork, int *info);
void dlarre(const char * const range, const uint32_t * const n, double *vl, double *vu, uint32_t *il, uint32_t *iu, double *d, double *e, double *e2, double *rtol1, double *rtol2, double *spltol, uint32_t *nsplit, uint32_t *isplit, uint32_t *m, double *w, double *werr, double *wgap, uint32_t *iblock, uint32_t *indexw, double *gers, double *pivmin, double *work, int *iwork, int *info);
void dlarrf(const uint32_t * const n, double *d, double *l, double *ld, uint32_t *clstrt, uint32_t *clend, double *w, double *wgap, double *werr, double *spdiam, double *clgapl, double *clgapr, double *pivmin, double *sigma, double *dplus, double *lplus, double *work, int *info);
void dlarrj(const uint32_t * const n, double *d, double *e2, uint32_t *ifirst, uint32_t *ilast, double *rtol, uint32_t *offset, double *w, double *werr, double *work, int *iwork, double *pivmin, double *spdiam, int *info);
void dlarrk(const uint32_t * const n, uint32_t *iw, double *gl, double *gu, double *d, double *e2, double *pivmin, double *reltol, double *w, double *werr, int *info);
void dlarrr(const uint32_t * const n, double *d, double *e, int *info);
void dlarrv(const uint32_t * const n, double *vl, double *vu, double *d, double *l, double *pivmin, uint32_t *isplit, uint32_t *m, uint32_t *dol, uint32_t *dou, double *minrgp, double *rtol1, double *rtol2, double *w, double *werr, double *wgap, uint32_t *iblock, uint32_t *indexw, double *gers, double *z, const uint32_t * const ldz, uint32_t *isuppz, double *work, int *iwork, int *info);
void dlartg(double *f, double *g, double *cs, double *sn, double *r);
void dlartgp(double *f, double *g, double *cs, double *sn, double *r);
void dlartgs(double *x, double *y, double *sigma, double *cs, double *sn);
void dlartv(const uint32_t * const n, double *x, const uint32_t * const incx, double *y, const uint32_t * const incy, double *c, double *s, const uint32_t * const incc);
void dlaruv(uint32_t *iseed, const uint32_t * const n, double *x);
void dlarz(const char * const side, uint32_t *m, const uint32_t * const n, uint32_t *l, double *v, const uint32_t * const incv, double *tau, double *c, const uint32_t * const ldc, double *work);
void dlarzb(const char * const side, const char * const trans, const char * const direct, const char * const storev, uint32_t *m, const uint32_t * const n, uint32_t *k, uint32_t *l, double *v, const uint32_t * const ldv, double *t, const uint32_t * const ldt, double *c, const uint32_t * const ldc, double *work, const uint32_t * const ldwork);
void dlarzt(const char * const direct, const char * const storev, const uint32_t * const n, uint32_t *k, double *v, const uint32_t * const ldv, double *tau, double *t, const uint32_t * const ldt);
void dlas2(double *f, double *g, double *h, double *ssmin, double *ssmax);
void dlascl(const char * const type, uint32_t *kl, uint32_t *ku, double *cfrom, double *cto, uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, int *info);
void dlasd0(const uint32_t * const n, uint32_t *sqre, double *d, double *e, double *u, const uint32_t * const ldu, double *vt, const uint32_t * const ldvt, uint32_t *smlsiz, int *iwork, double *work, int *info);
void dlasd1(uint32_t *nl, uint32_t *nr, uint32_t *sqre, double *d, double *alpha, double *beta, double *u, const uint32_t * const ldu, double *vt, const uint32_t * const ldvt, uint32_t *idxq, int *iwork, double *work, int *info);
void dlasd2(uint32_t *nl, uint32_t *nr, uint32_t *sqre, uint32_t *k, double *d, double *z, double *alpha, double *beta, double *u, const uint32_t * const ldu, double *vt, const uint32_t * const ldvt, double *dsigma, double *u2, const uint32_t * const ldu2, double *vt2, const uint32_t * const ldvt2, uint32_t *idxp, uint32_t *idx, uint32_t *idxc, uint32_t *idxq, uint32_t *coltyp, int *info);
void dlasd3(uint32_t *nl, uint32_t *nr, uint32_t *sqre, uint32_t *k, double *d, double *q, const uint32_t * const ldq, double *dsigma, double *u, const uint32_t * const ldu, double *u2, const uint32_t * const ldu2, double *vt, const uint32_t * const ldvt, double *vt2, const uint32_t * const ldvt2, uint32_t *idxc, uint32_t *ctot, double *z, int *info);
void dlasd4(const uint32_t * const n, uint32_t *i, double *d, double *z, double *delta, double *rho, double *sigma, double *work, int *info);
void dlasd5(uint32_t *i, double *d, double *z, double *delta, double *rho, double *dsigma, double *work);
void dlasd6(uint32_t *icompq, uint32_t *nl, uint32_t *nr, uint32_t *sqre, double *d, double *vf, double *vl, double *alpha, double *beta, uint32_t *idxq, uint32_t *perm, uint32_t *givptr, uint32_t *givcol, const uint32_t * const ldgcol, double *givnum, const uint32_t * const ldgnum, double *poles, double *difl, double *difr, double *z, uint32_t *k, double *c, double *s, double *work, int *iwork, int *info);
void dlasd7(uint32_t *icompq, uint32_t *nl, uint32_t *nr, uint32_t *sqre, uint32_t *k, double *d, double *z, double *zw, double *vf, double *vfw, double *vl, double *vlw, double *alpha, double *beta, double *dsigma, uint32_t *idx, uint32_t *idxp, uint32_t *idxq, uint32_t *perm, uint32_t *givptr, uint32_t *givcol, const uint32_t * const ldgcol, double *givnum, const uint32_t * const ldgnum, double *c, double *s, int *info);
void dlasd8(uint32_t *icompq, uint32_t *k, double *d, double *z, double *vf, double *vl, double *difl, double *difr, const uint32_t * const lddifr, double *dsigma, double *work, int *info);
void dlasda(uint32_t *icompq, uint32_t *smlsiz, const uint32_t * const n, uint32_t *sqre, double *d, double *e, double *u, const uint32_t * const ldu, double *vt, uint32_t *k, double *difl, double *difr, double *z, double *poles, uint32_t *givptr, uint32_t *givcol, const uint32_t * const ldgcol, uint32_t *perm, double *givnum, double *c, double *s, double *work, int *iwork, int *info);
void dlasdq(const char * const uplo, uint32_t *sqre, const uint32_t * const n, uint32_t *ncvt, uint32_t *nru, uint32_t *ncc, double *d, double *e, double *vt, const uint32_t * const ldvt, double *u, const uint32_t * const ldu, double *c, const uint32_t * const ldc, double *work, int *info);
void dlasdt(const uint32_t * const n, uint32_t *lvl, uint32_t *nd, uint32_t *inode, uint32_t *ndiml, uint32_t *ndimr, uint32_t *msub);
void dlaset(const char * const uplo, uint32_t *m, const uint32_t * const n, double *alpha, double *beta, double *a, const uint32_t * const lda);
void dlasq1(const uint32_t * const n, double *d, double *e, double *work, int *info);
void dlasq2(const uint32_t * const n, double *z, int *info);
void dlasq3(uint32_t *i0, uint32_t *n0, double *z, uint32_t *pp, double *dmin, double *sigma, double *desig, double *qmax, uint32_t *nfail, uint32_t *iter, uint32_t *ndiv, uint32_t *ieee, uint32_t *ttype, double *dmin1, double *dmin2, double *dn, double *dn1, double *dn2, double *g, double *tau);
void dlasq4(uint32_t *i0, uint32_t *n0, double *z, uint32_t *pp, uint32_t *n0in, double *dmin, double *dmin1, double *dmin2, double *dn, double *dn1, double *dn2, double *tau, uint32_t *ttype, double *g);
void dlasq5(uint32_t *i0, uint32_t *n0, double *z, uint32_t *pp, double *tau, double *sigma, double *dmin, double *dmin1, double *dmin2, double *dn, double *dnm1, double *dnm2, uint32_t *ieee, double *eps);
void dlasq6(uint32_t *i0, uint32_t *n0, double *z, uint32_t *pp, double *dmin, double *dmin1, double *dmin2, double *dn, double *dnm1, double *dnm2);
void dlasr(const char * const side, const char * const pivot, const char * const direct, uint32_t *m, const uint32_t * const n, double *c, double *s, double *a, const uint32_t * const lda);
void dlasrt(const char * const id, const uint32_t * const n, double *d, int *info);
void dlassq(const uint32_t * const n, double *x, const uint32_t * const incx, double *scale, double *sumsq);
void dlasv2(double *f, double *g, double *h, double *ssmin, double *ssmax, double *snr, double *csr, double *snl, double *csl);
void dlaswp(const uint32_t * const n, double *a, const uint32_t * const lda, uint32_t *k1, uint32_t *k2, int *ipiv, const uint32_t * const incx);
void dlasy2(uint32_t *ltranl, uint32_t *ltranr, uint32_t *isgn, uint32_t *n1, uint32_t *n2, double *tl, const uint32_t * const ldtl, double *tr, const uint32_t * const ldtr, double *b, const uint32_t * const ldb, double *scale, double *x, const uint32_t * const ldx, double *xnorm, int *info);
void dlasyf(const char * const uplo, const uint32_t * const n, uint32_t *nb, uint32_t *kb, double *a, const uint32_t * const lda, int *ipiv, double *w, const uint32_t * const ldw, int *info);
void dlatbs(const char * const uplo, const char * const trans, const char * const diag, const char * const normin, const uint32_t * const n, uint32_t *kd, double *ab, const uint32_t * const ldab, double *x, double *scale, double *cnorm, int *info);
void dlatdf(uint32_t *ijob, const uint32_t * const n, double *z, const uint32_t * const ldz, double *rhs, double *rdsum, double *rdscal, int *ipiv, int *jpiv);
void dlatps(const char * const uplo, const char * const trans, const char * const diag, const char * const normin, const uint32_t * const n, double *ap, double *x, double *scale, double *cnorm, int *info);
void dlatrd(const char * const uplo, const uint32_t * const n, uint32_t *nb, double *a, const uint32_t * const lda, double *e, double *tau, double *w, const uint32_t * const ldw);
void dlatrs(const char * const uplo, const char * const trans, const char * const diag, const char * const normin, const uint32_t * const n, double *a, const uint32_t * const lda, double *x, double *scale, double *cnorm, int *info);
void dlatrz(uint32_t *m, const uint32_t * const n, uint32_t *l, double *a, const uint32_t * const lda, double *tau, double *work);
void dlatzm(const char * const side, uint32_t *m, const uint32_t * const n, double *v, const uint32_t * const incv, double *tau, double *c1, double *c2, const uint32_t * const ldc, double *work);
void dlauu2(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, int *info);
void dlauum(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, int *info);
void dopgtr(const char * const uplo, const uint32_t * const n, double *ap, double *tau, double *q, const uint32_t * const ldq, double *work, int *info);
void dopmtr(const char * const side, const char * const uplo, const char * const trans, uint32_t *m, const uint32_t * const n, double *ap, double *tau, double *c, const uint32_t * const ldc, double *work, int *info);
void dorbdb(const char * const trans, const char * const signs, uint32_t *m, uint32_t *p, uint32_t *q, double *x11, const uint32_t * const ldx11, double *x12, const uint32_t * const ldx12, double *x21, const uint32_t * const ldx21, double *x22, const uint32_t * const ldx22, double *theta, double *phi, double *taup1, double *taup2, double *tauq1, double *tauq2, double *work, const uint32_t * const lwork, int *info);
void dorcsd(const char * const jobu1, const char * const jobu2, const char * const jobv1t, const char * const jobv2t, const char * const trans, const char * const signs, uint32_t *m, uint32_t *p, uint32_t *q, double *x11, const uint32_t * const ldx11, double *x12, const uint32_t * const ldx12, double *x21, const uint32_t * const ldx21, double *x22, const uint32_t * const ldx22, double *theta, double *u1, const uint32_t * const ldu1, double *u2, const uint32_t * const ldu2, double *v1t, const uint32_t * const ldv1t, double *v2t, const uint32_t * const ldv2t, double *work, const uint32_t * const lwork, int *iwork, int *info);
void dorg2l(uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *work, int *info);
void dorg2r(uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *work, int *info);
void dorgbr(const char * const vect, uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *work, const uint32_t * const lwork, int *info);
void dorghr(const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, double *a, const uint32_t * const lda, double *tau, double *work, const uint32_t * const lwork, int *info);
void dorgl2(uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *work, int *info);
void dorglq(uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *work, const uint32_t * const lwork, int *info);
void dorgql(uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *work, const uint32_t * const lwork, int *info);
void dorgqr(uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *work, const uint32_t * const lwork, int *info);
void dorgr2(uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *work, int *info);
void dorgrq(uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *work, const uint32_t * const lwork, int *info);
void dorgtr(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *tau, double *work, const uint32_t * const lwork, int *info);
void dorm2l(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *c, const uint32_t * const ldc, double *work, int *info);
void dorm2r(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *c, const uint32_t * const ldc, double *work, int *info);
void dormbr(const char * const vect, const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *c, const uint32_t * const ldc, double *work, const uint32_t * const lwork, int *info);
void dormhr(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, double *a, const uint32_t * const lda, double *tau, double *c, const uint32_t * const ldc, double *work, const uint32_t * const lwork, int *info);
void dorml2(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *c, const uint32_t * const ldc, double *work, int *info);
void dormlq(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *c, const uint32_t * const ldc, double *work, const uint32_t * const lwork, int *info);
void dormql(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *c, const uint32_t * const ldc, double *work, const uint32_t * const lwork, int *info);
void dormqr(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *c, const uint32_t * const ldc, double *work, const uint32_t * const lwork, int *info);
void dormr2(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *c, const uint32_t * const ldc, double *work, int *info);
void dormr3(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, uint32_t *l, double *a, const uint32_t * const lda, double *tau, double *c, const uint32_t * const ldc, double *work, int *info);
void dormrq(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, double *a, const uint32_t * const lda, double *tau, double *c, const uint32_t * const ldc, double *work, const uint32_t * const lwork, int *info);
void dormrz(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, uint32_t *l, double *a, const uint32_t * const lda, double *tau, double *c, const uint32_t * const ldc, double *work, const uint32_t * const lwork, int *info);
void dormtr(const char * const side, const char * const uplo, const char * const trans, uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *tau, double *c, const uint32_t * const ldc, double *work, const uint32_t * const lwork, int *info);
void dpbcon(const char * const uplo, const uint32_t * const n, uint32_t *kd, double *ab, const uint32_t * const ldab, double *anorm, double *rcond, double *work, int *iwork, int *info);
void dpbequ(const char * const uplo, const uint32_t * const n, uint32_t *kd, double *ab, const uint32_t * const ldab, double *s, double *scond, double *amax, int *info);
void dpbrfs(const char * const uplo, const uint32_t * const n, uint32_t *kd, uint32_t *nrhs, double *ab, const uint32_t * const ldab, double *afb, const uint32_t * const ldafb, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *ferr, double *berr, double *work, int *iwork, int *info);
void dpbstf(const char * const uplo, const uint32_t * const n, uint32_t *kd, double *ab, const uint32_t * const ldab, int *info);
void dpbsv(const char * const uplo, const uint32_t * const n, uint32_t *kd, uint32_t *nrhs, double *ab, const uint32_t * const ldab, double *b, const uint32_t * const ldb, int *info);
void dpbsvx(const char * const fact, const char * const uplo, const uint32_t * const n, uint32_t *kd, uint32_t *nrhs, double *ab, const uint32_t * const ldab, double *afb, const uint32_t * const ldafb, const char * const equed, double *s, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);
void dpbtf2(const char * const uplo, const uint32_t * const n, uint32_t *kd, double *ab, const uint32_t * const ldab, int *info);
void dpbtrf(const char * const uplo, const uint32_t * const n, uint32_t *kd, double *ab, const uint32_t * const ldab, int *info);
void dpbtrs(const char * const uplo, const uint32_t * const n, uint32_t *kd, uint32_t *nrhs, double *ab, const uint32_t * const ldab, double *b, const uint32_t * const ldb, int *info);
void dpftrf(const char * const transr, const char * const uplo, const uint32_t * const n, double *a, int *info);
void dpftri(const char * const transr, const char * const uplo, const uint32_t * const n, double *a, int *info);
void dpftrs(const char * const transr, const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *a, double *b, const uint32_t * const ldb, int *info);
void dpocon(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *anorm, double *rcond, double *work, int *iwork, int *info);
void dpoequ(const uint32_t * const n, double *a, const uint32_t * const lda, double *s, double *scond, double *amax, int *info);
void dpoequb(const uint32_t * const n, double *a, const uint32_t * const lda, double *s, double *scond, double *amax, int *info);
void dporfs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *af, const uint32_t * const ldaf, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *ferr, double *berr, double *work, int *iwork, int *info);
void dposv(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, int *info);
void dposvx(const char * const fact, const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *af, const uint32_t * const ldaf, const char * const equed, double *s, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);
void dpotf2(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, int *info);
void dpotrf(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, int *info);
void dpotri(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, int *info);
void dpotrs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, int *info);
void dppcon(const char * const uplo, const uint32_t * const n, double *ap, double *anorm, double *rcond, double *work, int *iwork, int *info);
void dppequ(const char * const uplo, const uint32_t * const n, double *ap, double *s, double *scond, double *amax, int *info);
void dpprfs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *ap, double *afp, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *ferr, double *berr, double *work, int *iwork, int *info);
void dppsv(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *ap, double *b, const uint32_t * const ldb, int *info);
void dppsvx(const char * const fact, const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *ap, double *afp, const char * const equed, double *s, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);
void dpptrf(const char * const uplo, const uint32_t * const n, double *ap, int *info);
void dpptri(const char * const uplo, const uint32_t * const n, double *ap, int *info);
void dpptrs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *ap, double *b, const uint32_t * const ldb, int *info);
void dpstf2(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, uint32_t *piv, uint32_t *rank, double *tol, double *work, int *info);
void dpstrf(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, uint32_t *piv, uint32_t *rank, double *tol, double *work, int *info);
void dptcon(const uint32_t * const n, double *d, double *e, double *anorm, double *rcond, double *work, int *info);
void dpteqr(const char * const compz, const uint32_t * const n, double *d, double *e, double *z, const uint32_t * const ldz, double *work, int *info);
void dptrfs(const uint32_t * const n, uint32_t *nrhs, double *d, double *e, double *df, double *ef, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *ferr, double *berr, double *work, int *info);
void dptsv(const uint32_t * const n, uint32_t *nrhs, double *d, double *e, double *b, const uint32_t * const ldb, int *info);
void dptsvx(const char * const fact, const uint32_t * const n, uint32_t *nrhs, double *d, double *e, double *df, double *ef, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *rcond, double *ferr, double *berr, double *work, int *info);
void dpttrf(const uint32_t * const n, double *d, double *e, int *info);
void dpttrs(const uint32_t * const n, uint32_t *nrhs, double *d, double *e, double *b, const uint32_t * const ldb, int *info);
void dptts2(const uint32_t * const n, uint32_t *nrhs, double *d, double *e, double *b, const uint32_t * const ldb);
void drscl(const uint32_t * const n, double *sa, double *sx, const uint32_t * const incx);
void dsbev(const char * const jobz, const char * const uplo, const uint32_t * const n, uint32_t *kd, double *ab, const uint32_t * const ldab, double *w, double *z, const uint32_t * const ldz, double *work, int *info);
void dsbevd(const char * const jobz, const char * const uplo, const uint32_t * const n, uint32_t *kd, double *ab, const uint32_t * const ldab, double *w, double *z, const uint32_t * const ldz, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void dsbevx(const char * const jobz, const char * const range, const char * const uplo, const uint32_t * const n, uint32_t *kd, double *ab, const uint32_t * const ldab, double *q, const uint32_t * const ldq, double *vl, double *vu, uint32_t *il, uint32_t *iu, double *abstol, uint32_t *m, double *w, double *z, const uint32_t * const ldz, double *work, int *iwork, uint32_t *ifail, int *info);
void dsbgst(const char * const vect, const char * const uplo, const uint32_t * const n, uint32_t *ka, uint32_t *kb, double *ab, const uint32_t * const ldab, double *bb, const uint32_t * const ldbb, double *x, const uint32_t * const ldx, double *work, int *info);
void dsbgv(const char * const jobz, const char * const uplo, const uint32_t * const n, uint32_t *ka, uint32_t *kb, double *ab, const uint32_t * const ldab, double *bb, const uint32_t * const ldbb, double *w, double *z, const uint32_t * const ldz, double *work, int *info);
void dsbgvd(const char * const jobz, const char * const uplo, const uint32_t * const n, uint32_t *ka, uint32_t *kb, double *ab, const uint32_t * const ldab, double *bb, const uint32_t * const ldbb, double *w, double *z, const uint32_t * const ldz, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void dsbgvx(const char * const jobz, const char * const range, const char * const uplo, const uint32_t * const n, uint32_t *ka, uint32_t *kb, double *ab, const uint32_t * const ldab, double *bb, const uint32_t * const ldbb, double *q, const uint32_t * const ldq, double *vl, double *vu, uint32_t *il, uint32_t *iu, double *abstol, uint32_t *m, double *w, double *z, const uint32_t * const ldz, double *work, int *iwork, uint32_t *ifail, int *info);
void dsbtrd(const char * const vect, const char * const uplo, const uint32_t * const n, const uint32_t * const kd, double *ab, const uint32_t * const ldab, double *d, double *e, double *q, const uint32_t * const ldq, double *work, int *info);
void dsfrk(const char * const transr, const char * const uplo, const char * const trans, const uint32_t * const n, uint32_t *k, double *alpha, double *a, const uint32_t * const lda, double *beta, double *c);
void dspcon(const char * const uplo, const uint32_t * const n, double *ap, int *ipiv, double *anorm, double *rcond, double *work, int *iwork, int *info);
void dspev(const char * const jobz, const char * const uplo, const uint32_t * const n, double *ap, double *w, double *z, const uint32_t * const ldz, double *work, int *info);
void dspevd(const char * const jobz, const char * const uplo, const uint32_t * const n, double *ap, double *w, double *z, const uint32_t * const ldz, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void dspevx(const char * const jobz, const char * const range, const char * const uplo, const uint32_t * const n, double *ap, double *vl, double *vu, uint32_t *il, uint32_t *iu, double *abstol, uint32_t *m, double *w, double *z, const uint32_t * const ldz, double *work, int *iwork, uint32_t *ifail, int *info);
void dspgst(uint32_t *itype, const char * const uplo, const uint32_t * const n, double *ap, double *bp, int *info);
void dspgv(uint32_t *itype, const char * const jobz, const char * const uplo, const uint32_t * const n, double *ap, double *bp, double *w, double *z, const uint32_t * const ldz, double *work, int *info);
void dspgvd(uint32_t *itype, const char * const jobz, const char * const uplo, const uint32_t * const n, double *ap, double *bp, double *w, double *z, const uint32_t * const ldz, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void dspgvx(uint32_t *itype, const char * const jobz, const char * const range, const char * const uplo, const uint32_t * const n, double *ap, double *bp, double *vl, double *vu, uint32_t *il, uint32_t *iu, double *abstol, uint32_t *m, double *w, double *z, const uint32_t * const ldz, double *work, int *iwork, uint32_t *ifail, int *info);
void dsprfs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *ap, double *afp, int *ipiv, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *ferr, double *berr, double *work, int *iwork, int *info);
void dspsv(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *ap, int *ipiv, double *b, const uint32_t * const ldb, int *info);
void dspsvx(const char * const fact, const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *ap, double *afp, int *ipiv, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);
void dsptrd(const char * const uplo, const uint32_t * const n, double *ap, double *d, double *e, double *tau, int *info);
void dsptrf(const char * const uplo, const uint32_t * const n, double *ap, int *ipiv, int *info);
void dsptri(const char * const uplo, const uint32_t * const n, double *ap, int *ipiv, double *work, int *info);
void dsptrs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *ap, int *ipiv, double *b, const uint32_t * const ldb, int *info);
void dstebz(const char * const range, const char * const order, const uint32_t * const n, double *vl, double *vu, uint32_t *il, uint32_t *iu, double *abstol, double *d, double *e, uint32_t *m, uint32_t *nsplit, double *w, uint32_t *iblock, uint32_t *isplit, double *work, int *iwork, int *info);
void dstedc(const char * const compz, const uint32_t * const n, double *d, double *e, double *z, const uint32_t * const ldz, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void dstegr(const char * const jobz, const char * const range, const uint32_t * const n, double *d, double *e, double *vl, double *vu, uint32_t *il, uint32_t *iu, double *abstol, uint32_t *m, double *w, double *z, const uint32_t * const ldz, uint32_t *isuppz, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void dstein(const uint32_t * const n, double *d, double *e, uint32_t *m, double *w, uint32_t *iblock, uint32_t *isplit, double *z, const uint32_t * const ldz, double *work, int *iwork, uint32_t *ifail, int *info);
void dstemr(const char * const jobz, const char * const range, const uint32_t * const n, double *d, double *e, double *vl, double *vu, uint32_t *il, uint32_t *iu, uint32_t *m, double *w, double *z, const uint32_t * const ldz, uint32_t *nzc, uint32_t *isuppz, uint32_t *tryrac, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void dsteqr(const char * const compz, const uint32_t * const n, double *d, double *e, double *z, const uint32_t * const ldz, double *work, int *info);
void dsterf(const uint32_t * const n, double *d, double *e, int *info);
void dstev(const char * const jobz, const uint32_t * const n, double *d, double *e, double *z, const uint32_t * const ldz, double *work, int *info);
void dstevd(const char * const jobz, const uint32_t * const n, double *d, double *e, double *z, const uint32_t * const ldz, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void dstevr(const char * const jobz, const char * const range, const uint32_t * const n, double *d, double *e, double *vl, double *vu, uint32_t *il, uint32_t *iu, double *abstol, uint32_t *m, double *w, double *z, const uint32_t * const ldz, uint32_t *isuppz, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void dstevx(const char * const jobz, const char * const range, const uint32_t * const n, double *d, double *e, double *vl, double *vu, uint32_t *il, uint32_t *iu, double *abstol, uint32_t *m, double *w, double *z, const uint32_t * const ldz, double *work, int *iwork, uint32_t *ifail, int *info);
void dsycon(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, int *ipiv, double *anorm, double *rcond, double *work, int *iwork, int *info);
void dsyconv(const char * const uplo, const char * const way, const uint32_t * const n, double *a, const uint32_t * const lda, int *ipiv, double *work, int *info);
void dsyequb(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *s, double *scond, double *amax, double *work, int *info);
void dsyev(const char * const jobz, const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *w, double *work, const uint32_t * const lwork, int *info);
void dsyevd(const char * const jobz, const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *w, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void dsyevr(const char * const jobz, const char * const range, const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *vl, double *vu, uint32_t *il, uint32_t *iu, double *abstol, uint32_t *m, double *w, double *z, const uint32_t * const ldz, uint32_t *isuppz, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void dsyevx(const char * const jobz, const char * const range, const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *vl, double *vu, uint32_t *il, uint32_t *iu, double *abstol, uint32_t *m, double *w, double *z, const uint32_t * const ldz, double *work, const uint32_t * const lwork, int *iwork, uint32_t *ifail, int *info);
void dsygs2(uint32_t *itype, const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, int *info);
void dsygst(uint32_t *itype, const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, int *info);
void dsygv(uint32_t *itype, const char * const jobz, const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *w, double *work, const uint32_t * const lwork, int *info);
void dsygvd(uint32_t *itype, const char * const jobz, const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *w, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void dsygvx(uint32_t *itype, const char * const jobz, const char * const range, const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *vl, double *vu, uint32_t *il, uint32_t *iu, double *abstol, uint32_t *m, double *w, double *z, const uint32_t * const ldz, double *work, const uint32_t * const lwork, int *iwork, uint32_t *ifail, int *info);
void dsyrfs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *af, const uint32_t * const ldaf, int *ipiv, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *ferr, double *berr, double *work, int *iwork, int *info);
void dsysv(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, int *ipiv, double *b, const uint32_t * const ldb, double *work, const uint32_t * const lwork, int *info);
void dsysvx(const char * const fact, const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *af, const uint32_t * const ldaf, int *ipiv, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *rcond, double *ferr, double *berr, double *work, const uint32_t * const lwork, int *iwork, int *info);
void dsyswapr(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, uint32_t *i1, uint32_t *i2);
void dsytd2(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *d, double *e, double *tau, int *info);
void dsytf2(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, int *ipiv, int *info);
void dsytrd(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *d, double *e, double *tau, double *work, const uint32_t * const lwork, int *info);
void dsytrf(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, int *ipiv, double *work, const uint32_t * const lwork, int *info);
void dsytri(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, int *ipiv, double *work, int *info);
void dsytri2(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, int *ipiv, double *work, const uint32_t * const lwork, int *info);
void dsytri2x(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, int *ipiv, double *work, uint32_t *nb, int *info);
void dsytrs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, int *ipiv, double *b, const uint32_t * const ldb, int *info);
void dsytrs2(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, int *ipiv, double *b, const uint32_t * const ldb, double *work, int *info);
void dtbcon(const char * const norm, const char * const uplo, const char * const diag, const uint32_t * const n, uint32_t *kd, double *ab, const uint32_t * const ldab, double *rcond, double *work, int *iwork, int *info);
void dtbrfs(const char * const uplo, const char * const trans, const char * const diag, const uint32_t * const n, uint32_t *kd, uint32_t *nrhs, double *ab, const uint32_t * const ldab, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *ferr, double *berr, double *work, int *iwork, int *info);
void dtbtrs(const char * const uplo, const char * const trans, const char * const diag, const uint32_t * const n, uint32_t *kd, uint32_t *nrhs, double *ab, const uint32_t * const ldab, double *b, const uint32_t * const ldb, int *info);
void dtfsm(const char * const transr, const char * const side, const char * const uplo, const char * const trans, const char * const diag, uint32_t *m, const uint32_t * const n, double *alpha, double *a, double *b, const uint32_t * const ldb);
void dtftri(const char * const transr, const char * const uplo, const char * const diag, const uint32_t * const n, double *a, int *info);
void dtfttp(const char * const transr, const char * const uplo, const uint32_t * const n, double *arf, double *ap, int *info);
void dtfttr(const char * const transr, const char * const uplo, const uint32_t * const n, double *arf, double *a, const uint32_t * const lda, int *info);
void dtgevc(const char * const side, const char * const howmny, uint32_t *select, const uint32_t * const n, double *s, const uint32_t * const lds, double *p, const uint32_t * const ldp, double *vl, const uint32_t * const ldvl, double *vr, const uint32_t * const ldvr, uint32_t *mm, uint32_t *m, double *work, int *info);
void dtgex2(uint32_t *wantq, uint32_t *wantz, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *q, const uint32_t * const ldq, double *z, const uint32_t * const ldz, uint32_t *j1_, uint32_t *n1, uint32_t *n2, double *work, const uint32_t * const lwork, int *info);
void dtgexc(uint32_t *wantq, uint32_t *wantz, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *q, const uint32_t * const ldq, double *z, const uint32_t * const ldz, uint32_t *ifst, uint32_t *ilst, double *work, const uint32_t * const lwork, int *info);
void dtgsen(uint32_t *ijob, uint32_t *wantq, uint32_t *wantz, uint32_t *select, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *alphar, double *alphai, double *beta, double *q, const uint32_t * const ldq, double *z, const uint32_t * const ldz, uint32_t *m, double *pl, double *pr, double *dif, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void dtgsja(const char * const jobu, const char * const jobv, const char * const jobq, uint32_t *m, uint32_t *p, const uint32_t * const n, uint32_t *k, uint32_t *l, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *tola, double *tolb, double *alpha, double *beta, double *u, const uint32_t * const ldu, double *v, const uint32_t * const ldv, double *q, const uint32_t * const ldq, double *work, uint32_t *ncycle, int *info);
void dtgsna(const char * const job, const char * const howmny, uint32_t *select, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *vl, const uint32_t * const ldvl, double *vr, const uint32_t * const ldvr, double *s, double *dif, uint32_t *mm, uint32_t *m, double *work, const uint32_t * const lwork, int *iwork, int *info);
void dtgsy2(const char * const trans, uint32_t *ijob, uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *c, const uint32_t * const ldc, double *d, const uint32_t * const ldd, double *e, const uint32_t * const lde, double *f, const uint32_t * const ldf, double *scale, double *rdsum, double *rdscal, int *iwork, uint32_t *pq, int *info);
void dtgsyl(const char * const trans, uint32_t *ijob, uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *c, const uint32_t * const ldc, double *d, const uint32_t * const ldd, double *e, const uint32_t * const lde, double *f, const uint32_t * const ldf, double *scale, double *dif, double *work, const uint32_t * const lwork, int *iwork, int *info);
void dtpcon(const char * const norm, const char * const uplo, const char * const diag, const uint32_t * const n, double *ap, double *rcond, double *work, int *iwork, int *info);
void dtpmqrt(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, uint32_t *l, uint32_t *nb, double *v, const uint32_t * const ldv, double *t, const uint32_t * const ldt, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *work, int *info);
void dtpqrt(uint32_t *m, const uint32_t * const n, uint32_t *l, uint32_t *nb, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *t, const uint32_t * const ldt, double *work, int *info);
void dtpqrt2(uint32_t *m, const uint32_t * const n, uint32_t *l, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *t, const uint32_t * const ldt, int *info);
void dtprfb(const char * const side, const char * const trans, const char * const direct, const char * const storev, uint32_t *m, const uint32_t * const n, uint32_t *k, uint32_t *l, double *v, const uint32_t * const ldv, double *t, const uint32_t * const ldt, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *work, const uint32_t * const ldwork);
void dtprfs(const char * const uplo, const char * const trans, const char * const diag, const uint32_t * const n, uint32_t *nrhs, double *ap, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *ferr, double *berr, double *work, int *iwork, int *info);
void dtptri(const char * const uplo, const char * const diag, const uint32_t * const n, double *ap, int *info);
void dtptrs(const char * const uplo, const char * const trans, const char * const diag, const uint32_t * const n, uint32_t *nrhs, double *ap, double *b, const uint32_t * const ldb, int *info);
void dtpttf(const char * const transr, const char * const uplo, const uint32_t * const n, double *ap, double *arf, int *info);
void dtpttr(const char * const uplo, const uint32_t * const n, double *ap, double *a, const uint32_t * const lda, int *info);
void dtrcon(const char * const norm, const char * const uplo, const char * const diag, const uint32_t * const n, double *a, const uint32_t * const lda, double *rcond, double *work, int *iwork, int *info);
void dtrevc(const char * const side, const char * const howmny, uint32_t *select, const uint32_t * const n, double *t, const uint32_t * const ldt, double *vl, const uint32_t * const ldvl, double *vr, const uint32_t * const ldvr, uint32_t *mm, uint32_t *m, double *work, int *info);
void dtrexc(const char * const compq, const uint32_t * const n, double *t, const uint32_t * const ldt, double *q, const uint32_t * const ldq, uint32_t *ifst, uint32_t *ilst, double *work, int *info);
void dtrrfs(const char * const uplo, const char * const trans, const char * const diag, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *x, const uint32_t * const ldx, double *ferr, double *berr, double *work, int *iwork, int *info);
void dtrsen(const char * const job, const char * const compq, uint32_t *select, const uint32_t * const n, double *t, const uint32_t * const ldt, double *q, const uint32_t * const ldq, double *wr, double *wi, uint32_t *m, double *s, double *sep, double *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void dtrsna(const char * const job, const char * const howmny, uint32_t *select, const uint32_t * const n, double *t, const uint32_t * const ldt, double *vl, const uint32_t * const ldvl, double *vr, const uint32_t * const ldvr, double *s, double *sep, uint32_t *mm, uint32_t *m, double *work, const uint32_t * const ldwork, int *iwork, int *info);
void dtrsyl(const char * const trana, const char * const tranb, uint32_t *isgn, uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, double *c, const uint32_t * const ldc, double *scale, int *info);
void dtrti2(const char * const uplo, const char * const diag, const uint32_t * const n, double *a, const uint32_t * const lda, int *info);
void dtrtri(const char * const uplo, const char * const diag, const uint32_t * const n, double *a, const uint32_t * const lda, int *info);
void dtrtrs(const char * const uplo, const char * const trans, const char * const diag, const uint32_t * const n, uint32_t *nrhs, double *a, const uint32_t * const lda, double *b, const uint32_t * const ldb, int *info);
void dtrttf(const char * const transr, const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *arf, int *info);
void dtrttp(const char * const uplo, const uint32_t * const n, double *a, const uint32_t * const lda, double *ap, int *info);
void dtzrqf(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *tau, int *info);
void dtzrzf(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda, double *tau, double *work, const uint32_t * const lwork, int *info);
double dzsum1(const uint32_t * const n, double *cx, const uint32_t * const incx);
int iladiag(const char * const diag);
int iladlc(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda);
int iladlr(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda);
int ilaenv(uint32_t *ispec, const char * const name, const char * const opts, uint32_t *n1, uint32_t *n2, uint32_t *n3, uint32_t *n4, uint32_t name_len, uint32_t opts_len);
int ilaprec(const char * const prec);
int ilatrans(const char * const trans);
int ilauplo(const char * const uplo);
void ilaver(uint32_t *vers_major, uint32_t *vers_minor, uint32_t *vers_patch);
int ilazlc(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda);
int ilazlr(uint32_t *m, const uint32_t * const n, double *a, const uint32_t * const lda);
int iparmq(uint32_t *ispec, const char * const name, const char * const opts, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, const uint32_t * const lwork);
int izmax1(const uint32_t * const n, double *cx, const uint32_t * const incx);
int lsamen(const uint32_t * const n, const char * const ca, const char * const cb, uint32_t ca_len, uint32_t cb_len);
void slasdt(const uint32_t * const n, uint32_t *lvl, uint32_t *nd, uint32_t *inode, uint32_t *ndiml, uint32_t *ndimr, uint32_t *msub);
void xerbla_array(const char * const srname_array, uint32_t *srname_len, int *info);
void zbbcsd(const char * const jobu1, const char * const jobu2, const char * const jobv1t, const char * const jobv2t, const char * const trans, uint32_t *m, uint32_t *p, uint32_t *q, Complex16 *theta, Complex16 *phi, Complex16 *u1, const uint32_t * const ldu1, Complex16 *u2, const uint32_t * const ldu2, Complex16 *v1t, const uint32_t * const ldv1t, Complex16 *v2t, const uint32_t * const ldv2t, Complex16 *b11d, Complex16 *b11e, Complex16 *b12d, Complex16 *b12e, Complex16 *b21d, Complex16 *b21e, Complex16 *b22d, Complex16 *b22e, double *rwork, const uint32_t * const lrwork, int *info);
void zbdsqr(const char * const uplo, const uint32_t * const n, uint32_t *ncvt, uint32_t *nru, uint32_t *ncc, Complex16 *d, Complex16 *e, Complex16 *vt, const uint32_t * const ldvt, Complex16 *u, const uint32_t * const ldu, Complex16 *c, const uint32_t * const ldc, double *rwork, int *info);
void zdrscl(const uint32_t * const n, Complex16 *sa, Complex16 *sx, const uint32_t * const incx);
void zgbbrd(const char * const vect, uint32_t *m, const uint32_t * const n, uint32_t *ncc, uint32_t *kl, uint32_t *ku, Complex16 *ab, const uint32_t * const ldab, Complex16 *d, Complex16 *e, Complex16 *q, const uint32_t * const ldq, Complex16 *pt, const uint32_t * const ldpt, Complex16 *c, const uint32_t * const ldc, Complex16 *work, double *rwork, int *info);
void zgbcon(const char * const norm, const uint32_t * const n, uint32_t *kl, uint32_t *ku, Complex16 *ab, const uint32_t * const ldab, int *ipiv, Complex16 *anorm, Complex16 *rcond, Complex16 *work, double *rwork, int *info);
void zgbequ(uint32_t *m, const uint32_t * const n, uint32_t *kl, uint32_t *ku, Complex16 *ab, const uint32_t * const ldab, Complex16 *r, Complex16 *c, Complex16 *rowcnd, Complex16 *colcnd, Complex16 *amax, int *info);
void zgbequb(uint32_t *m, const uint32_t * const n, uint32_t *kl, uint32_t *ku, Complex16 *ab, const uint32_t * const ldab, Complex16 *r, Complex16 *c, Complex16 *rowcnd, Complex16 *colcnd, Complex16 *amax, int *info);
void zgbrfs(const char * const trans, const uint32_t * const n, uint32_t *kl, uint32_t *ku, uint32_t *nrhs, Complex16 *ab, const uint32_t * const ldab, Complex16 *afb, const uint32_t * const ldafb, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zgbsv(const uint32_t * const n, uint32_t *kl, uint32_t *ku, uint32_t *nrhs, Complex16 *ab, const uint32_t * const ldab, int *ipiv, Complex16 *b, const uint32_t * const ldb, int *info);
void zgbsvx(const char * const fact, const char * const trans, const uint32_t * const n, uint32_t *kl, uint32_t *ku, uint32_t *nrhs, Complex16 *ab, const uint32_t * const ldab, Complex16 *afb, const uint32_t * const ldafb, int *ipiv, const char * const equed, Complex16 *r, Complex16 *c, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *rcond, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zgbtf2(uint32_t *m, const uint32_t * const n, uint32_t *kl, uint32_t *ku, Complex16 *ab, const uint32_t * const ldab, int *ipiv, int *info);
void zgbtrf(uint32_t *m, const uint32_t * const n, uint32_t *kl, uint32_t *ku, Complex16 *ab, const uint32_t * const ldab, int *ipiv, int *info);
void zgbtrs(const char * const trans, const uint32_t * const n, uint32_t *kl, uint32_t *ku, uint32_t *nrhs, Complex16 *ab, const uint32_t * const ldab, int *ipiv, Complex16 *b, const uint32_t * const ldb, int *info);
void zgebak(const char * const job, const char * const side, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, Complex16 *scale, uint32_t *m, Complex16 *v, const uint32_t * const ldv, int *info);
void zgebal(const char * const job, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, uint32_t *ilo, uint32_t *ihi, Complex16 *scale, int *info);
void zgebd2(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *d, Complex16 *e, Complex16 *tauq, Complex16 *taup, Complex16 *work, int *info);
void zgebrd(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *d, Complex16 *e, Complex16 *tauq, Complex16 *taup, Complex16 *work, const uint32_t * const lwork, int *info);
void zgecon(const char * const norm, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *anorm, Complex16 *rcond, Complex16 *work, double *rwork, int *info);
void zgeequ(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *r, Complex16 *c, Complex16 *rowcnd, Complex16 *colcnd, Complex16 *amax, int *info);
void zgeequb(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *r, Complex16 *c, Complex16 *rowcnd, Complex16 *colcnd, Complex16 *amax, int *info);
void zgees(const char * const jobvs, const char * const sort, uint32_t (*select)(), const uint32_t * const n, Complex16 *a, const uint32_t * const lda, uint32_t *sdim, Complex16 *w, Complex16 *vs, const uint32_t * const ldvs, Complex16 *work, const uint32_t * const lwork, double *rwork, uint32_t *bwork, int *info);
void zgeesx(const char * const jobvs, const char * const sort, uint32_t (*select)(), const char * const sense, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, uint32_t *sdim, Complex16 *w, Complex16 *vs, const uint32_t * const ldvs, Complex16 *rconde, Complex16 *rcondv, Complex16 *work, const uint32_t * const lwork, double *rwork, uint32_t *bwork, int *info);
void zgeev(const char * const jobvl, const char * const jobvr, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *w, Complex16 *vl, const uint32_t * const ldvl, Complex16 *vr, const uint32_t * const ldvr, Complex16 *work, const uint32_t * const lwork, double *rwork, int *info);
void zgeevx(const char * const balanc, const char * const jobvl, const char * const jobvr, const char * const sense, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *w, Complex16 *vl, const uint32_t * const ldvl, Complex16 *vr, const uint32_t * const ldvr, uint32_t *ilo, uint32_t *ihi, Complex16 *scale, Complex16 *abnrm, Complex16 *rconde, Complex16 *rcondv, Complex16 *work, const uint32_t * const lwork, double *rwork, int *info);
void zgegs(const char * const jobvsl, const char * const jobvsr, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *alpha, Complex16 *beta, Complex16 *vsl, const uint32_t * const ldvsl, Complex16 *vsr, const uint32_t * const ldvsr, Complex16 *work, const uint32_t * const lwork, double *rwork, int *info);
void zgegv(const char * const jobvl, const char * const jobvr, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *alpha, Complex16 *beta, Complex16 *vl, const uint32_t * const ldvl, Complex16 *vr, const uint32_t * const ldvr, Complex16 *work, const uint32_t * const lwork, double *rwork, int *info);
void zgehd2(const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, int *info);
void zgehrd(const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zgelq2(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, int *info);
void zgelqf(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zgels(const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *work, const uint32_t * const lwork, int *info);
void zgelsd(uint32_t *m, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *s, Complex16 *rcond, uint32_t *rank, Complex16 *work, const uint32_t * const lwork, double *rwork, int *iwork, int *info);
void zgelss(uint32_t *m, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *s, Complex16 *rcond, uint32_t *rank, Complex16 *work, const uint32_t * const lwork, double *rwork, int *info);
void zgelsx(uint32_t *m, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, uint32_t *jpvt, Complex16 *rcond, uint32_t *rank, Complex16 *work, double *rwork, int *info);
void zgelsy(uint32_t *m, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, uint32_t *jpvt, Complex16 *rcond, uint32_t *rank, Complex16 *work, const uint32_t * const lwork, double *rwork, int *info);
void zgemqrt(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, uint32_t *nb, Complex16 *v, const uint32_t * const ldv, Complex16 *t, const uint32_t * const ldt, Complex16 *c, const uint32_t * const ldc, Complex16 *work, int *info);
void zgeql2(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, int *info);
void zgeqlf(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zgeqp3(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, uint32_t *jpvt, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, double *rwork, int *info);
void zgeqpf(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, uint32_t *jpvt, Complex16 *tau, Complex16 *work, double *rwork, int *info);
void zgeqr2(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, int *info);
void zgeqr2p(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, int *info);
void zgeqrf(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zgeqrfp(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zgeqrt(uint32_t *m, const uint32_t * const n, uint32_t *nb, Complex16 *a, const uint32_t * const lda, Complex16 *t, const uint32_t * const ldt, Complex16 *work, int *info);
void zgeqrt2(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *t, const uint32_t * const ldt, int *info);
void zgeqrt3(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *t, const uint32_t * const ldt, int *info);
void zgerfs(const char * const trans, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *af, const uint32_t * const ldaf, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zgerq2(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, int *info);
void zgerqf(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zgesc2(const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *rhs, int *ipiv, int *jpiv, Complex16 *scale);
void zgesdd(const char * const jobz, uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *s, Complex16 *u, const uint32_t * const ldu, Complex16 *vt, const uint32_t * const ldvt, Complex16 *work, const uint32_t * const lwork, double *rwork, int *iwork, int *info);
void zgesv(const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *b, const uint32_t * const ldb, int *info);
void zgesvd(const char * const jobu, const char * const jobvt, uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, double *s, Complex16 *u, const uint32_t * const ldu, Complex16 *vt, const uint32_t * const ldvt, Complex16 *work, const uint32_t * const lwork, double *rwork, int *info);
void zgesvx(const char * const fact, const char * const trans, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *af, const uint32_t * const ldaf, int *ipiv, const char * const equed, Complex16 *r, Complex16 *c, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *rcond, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zgetc2(const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, int *jpiv, int *info);
void zgetf2(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, int *info);
void zgetrf(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, int *info);
void zgetri(const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *work, const uint32_t * const lwork, int *info);
void zgetrs(const char * const trans, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *b, const uint32_t * const ldb, int *info);
void zggbak(const char * const job, const char * const side, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, Complex16 *lscale, Complex16 *rscale, uint32_t *m, Complex16 *v, const uint32_t * const ldv, int *info);
void zggbal(const char * const job, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, uint32_t *ilo, uint32_t *ihi, Complex16 *lscale, Complex16 *rscale, Complex16 *work, int *info);
void zgges(const char * const jobvsl, const char * const jobvsr, const char * const sort, uint32_t (*selctg)(), const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, uint32_t *sdim, Complex16 *alpha, Complex16 *beta, Complex16 *vsl, const uint32_t * const ldvsl, Complex16 *vsr, const uint32_t * const ldvsr, Complex16 *work, const uint32_t * const lwork, double *rwork, uint32_t *bwork, int *info);
void zggesx(const char * const jobvsl, const char * const jobvsr, const char * const sort, uint32_t (*selctg)(), const char * const sense, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, uint32_t *sdim, Complex16 *alpha, Complex16 *beta, Complex16 *vsl, const uint32_t * const ldvsl, Complex16 *vsr, const uint32_t * const ldvsr, Complex16 *rconde, Complex16 *rcondv, Complex16 *work, const uint32_t * const lwork, double *rwork, int *iwork, const uint32_t * const liwork, uint32_t *bwork, int *info);
void zggev(const char * const jobvl, const char * const jobvr, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *alpha, Complex16 *beta, Complex16 *vl, const uint32_t * const ldvl, Complex16 *vr, const uint32_t * const ldvr, Complex16 *work, const uint32_t * const lwork, double *rwork, int *info);
void zggevx(const char * const balanc, const char * const jobvl, const char * const jobvr, const char * const sense, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *alpha, Complex16 *beta, Complex16 *vl, const uint32_t * const ldvl, Complex16 *vr, const uint32_t * const ldvr, uint32_t *ilo, uint32_t *ihi, Complex16 *lscale, Complex16 *rscale, Complex16 *abnrm, Complex16 *bbnrm, Complex16 *rconde, Complex16 *rcondv, Complex16 *work, const uint32_t * const lwork, double *rwork, int *iwork, uint32_t *bwork, int *info);
void zggglm(const uint32_t * const n, uint32_t *m, uint32_t *p, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *d, Complex16 *x, Complex16 *y, Complex16 *work, const uint32_t * const lwork, int *info);
void zgghrd(const char * const compq, const char * const compz, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *q, const uint32_t * const ldq, Complex16 *z, const uint32_t * const ldz, int *info);
void zgglse(uint32_t *m, const uint32_t * const n, uint32_t *p, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *c, Complex16 *d, Complex16 *x, Complex16 *work, const uint32_t * const lwork, int *info);
void zggqrf(const uint32_t * const n, uint32_t *m, uint32_t *p, Complex16 *a, const uint32_t * const lda, Complex16 *taua, Complex16 *b, const uint32_t * const ldb, Complex16 *taub, Complex16 *work, const uint32_t * const lwork, int *info);
void zggrqf(uint32_t *m, uint32_t *p, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *taua, Complex16 *b, const uint32_t * const ldb, Complex16 *taub, Complex16 *work, const uint32_t * const lwork, int *info);
void zggsvd(const char * const jobu, const char * const jobv, const char * const jobq, uint32_t *m, const uint32_t * const n, uint32_t *p, uint32_t *k, uint32_t *l, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *alpha, Complex16 *beta, Complex16 *u, const uint32_t * const ldu, Complex16 *v, const uint32_t * const ldv, Complex16 *q, const uint32_t * const ldq, Complex16 *work, double *rwork, int *iwork, int *info);
void zggsvp(const char * const jobu, const char * const jobv, const char * const jobq, uint32_t *m, uint32_t *p, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *tola, Complex16 *tolb, uint32_t *k, uint32_t *l, Complex16 *u, const uint32_t * const ldu, Complex16 *v, const uint32_t * const ldv, Complex16 *q, const uint32_t * const ldq, int *iwork, double *rwork, Complex16 *tau, Complex16 *work, int *info);
void zgtcon(const char * const norm, const uint32_t * const n, Complex16 *dl, Complex16 *d, Complex16 *du, Complex16 *du2, int *ipiv, Complex16 *anorm, Complex16 *rcond, Complex16 *work, int *info);
void zgtrfs(const char * const trans, const uint32_t * const n, uint32_t *nrhs, Complex16 *dl, Complex16 *d, Complex16 *du, Complex16 *dlf, Complex16 *df, Complex16 *duf, Complex16 *du2, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zgtsv(const uint32_t * const n, uint32_t *nrhs, Complex16 *dl, Complex16 *d, Complex16 *du, Complex16 *b, const uint32_t * const ldb, int *info);
void zgtsvx(const char * const fact, const char * const trans, const uint32_t * const n, uint32_t *nrhs, Complex16 *dl, Complex16 *d, Complex16 *du, Complex16 *dlf, Complex16 *df, Complex16 *duf, Complex16 *du2, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *rcond, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zgttrf(const uint32_t * const n, Complex16 *dl, Complex16 *d, Complex16 *du, Complex16 *du2, int *ipiv, int *info);
void zgttrs(const char * const trans, const uint32_t * const n, uint32_t *nrhs, Complex16 *dl, Complex16 *d, Complex16 *du, Complex16 *du2, int *ipiv, Complex16 *b, const uint32_t * const ldb, int *info);
void zgtts2(uint32_t *itrans, const uint32_t * const n, uint32_t *nrhs, Complex16 *dl, Complex16 *d, Complex16 *du, Complex16 *du2, int *ipiv, Complex16 *b, const uint32_t * const ldb);
void zhbev(const char * const jobz, const char * const uplo, const uint32_t * const n, uint32_t *kd, Complex16 *ab, const uint32_t * const ldab, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, double *rwork, int *info);
void zhbevd(const char * const jobz, const char * const uplo, const uint32_t * const n, uint32_t *kd, Complex16 *ab, const uint32_t * const ldab, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, const uint32_t * const lwork, double *rwork, const uint32_t * const lrwork, int *iwork, const uint32_t * const liwork, int *info);
void zhbevx(const char * const jobz, const char * const range, const char * const uplo, const uint32_t * const n, uint32_t *kd, Complex16 *ab, const uint32_t * const ldab, Complex16 *q, const uint32_t * const ldq, Complex16 *vl, Complex16 *vu, uint32_t *il, uint32_t *iu, Complex16 *abstol, uint32_t *m, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, double *rwork, int *iwork, uint32_t *ifail, int *info);
void zhbgst(const char * const vect, const char * const uplo, const uint32_t * const n, uint32_t *ka, uint32_t *kb, Complex16 *ab, const uint32_t * const ldab, Complex16 *bb, const uint32_t * const ldbb, Complex16 *x, const uint32_t * const ldx, Complex16 *work, double *rwork, int *info);
void zhbgv(const char * const jobz, const char * const uplo, const uint32_t * const n, uint32_t *ka, uint32_t *kb, Complex16 *ab, const uint32_t * const ldab, Complex16 *bb, const uint32_t * const ldbb, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, double *rwork, int *info);
void zhbgvd(const char * const jobz, const char * const uplo, const uint32_t * const n, uint32_t *ka, uint32_t *kb, Complex16 *ab, const uint32_t * const ldab, Complex16 *bb, const uint32_t * const ldbb, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, const uint32_t * const lwork, double *rwork, const uint32_t * const lrwork, int *iwork, const uint32_t * const liwork, int *info);
void zhbgvx(const char * const jobz, const char * const range, const char * const uplo, const uint32_t * const n, uint32_t *ka, uint32_t *kb, Complex16 *ab, const uint32_t * const ldab, Complex16 *bb, const uint32_t * const ldbb, Complex16 *q, const uint32_t * const ldq, Complex16 *vl, Complex16 *vu, uint32_t *il, uint32_t *iu, Complex16 *abstol, uint32_t *m, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, double *rwork, int *iwork, uint32_t *ifail, int *info);
void zhbtrd(const char * const vect, const char * const uplo, const uint32_t * const n, const uint32_t * const kd, Complex16 *ab, const uint32_t * const ldab, double *d, double *e, Complex16 *q, const uint32_t * const ldq, Complex16 *work, int *info);
void zhecon(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *anorm, Complex16 *rcond, Complex16 *work, int *info);
void zheequb(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *s, Complex16 *scond, Complex16 *amax, Complex16 *work, int *info);
void zheev(const char * const jobz, const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, double *w, Complex16 *work, const uint32_t * const lwork, double *rwork, int *info);
void zheevd(const char * const jobz, const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, double *w, Complex16 *work, const uint32_t * const lwork, double *rwork, const uint32_t * const lrwork, int *iwork, const uint32_t * const liwork, int *info);
void zheevr(const char * const jobz, const char * const range, const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *vl, Complex16 *vu, uint32_t *il, uint32_t *iu, Complex16 *abstol, uint32_t *m, Complex16 *w, Complex16 *z, const uint32_t * const ldz, uint32_t *isuppz, Complex16 *work, const uint32_t * const lwork, double *rwork, const uint32_t * const lrwork, int *iwork, const uint32_t * const liwork, int *info);
void zheevx(const char * const jobz, const char * const range, const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *vl, Complex16 *vu, uint32_t *il, uint32_t *iu, Complex16 *abstol, uint32_t *m, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, const uint32_t * const lwork, double *rwork, int *iwork, uint32_t *ifail, int *info);
void zhegs2(uint32_t *itype, const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, int *info);
void zhegst(uint32_t *itype, const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, int *info);
void zhegv(uint32_t *itype, const char * const jobz, const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, double *w, Complex16 *work, const uint32_t * const lwork, double *rwork, int *info);
void zhegvd(uint32_t *itype, const char * const jobz, const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, double *w, Complex16 *work, const uint32_t * const lwork, double *rwork, const uint32_t * const lrwork, int *iwork, const uint32_t * const liwork, int *info);
void zhegvx(uint32_t *itype, const char * const jobz, const char * const range, const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *vl, Complex16 *vu, uint32_t *il, uint32_t *iu, Complex16 *abstol, uint32_t *m, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, const uint32_t * const lwork, double *rwork, int *iwork, uint32_t *ifail, int *info);
void zherfs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *af, const uint32_t * const ldaf, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zhesv(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *work, const uint32_t * const lwork, int *info);
void zhesvx(const char * const fact, const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *af, const uint32_t * const ldaf, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *rcond, Complex16 *ferr, Complex16 *berr, Complex16 *work, const uint32_t * const lwork, double *rwork, int *info);
void zheswapr(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, uint32_t *i1, uint32_t *i2);
void zhetd2(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *d, Complex16 *e, Complex16 *tau, int *info);
void zhetf2(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, int *info);
void zhetrd(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *d, Complex16 *e, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zhetrf(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *work, const uint32_t * const lwork, int *info);
void zhetri(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *work, int *info);
void zhetri2(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *work, const uint32_t * const lwork, int *info);
void zhetri2x(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *work, uint32_t *nb, int *info);
void zhetrs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *b, const uint32_t * const ldb, int *info);
void zhetrs2(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *work, int *info);
void zhfrk(const char * const transr, const char * const uplo, const char * const trans, const uint32_t * const n, uint32_t *k, Complex16 *alpha, Complex16 *a, const uint32_t * const lda, Complex16 *beta, Complex16 *c);
void zhgeqz(const char * const job, const char * const compq, const char * const compz, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, Complex16 *h, const uint32_t * const ldh, Complex16 *t, const uint32_t * const ldt, Complex16 *alpha, Complex16 *beta, Complex16 *q, const uint32_t * const ldq, Complex16 *z, const uint32_t * const ldz, Complex16 *work, const uint32_t * const lwork, double *rwork, int *info);
void zhpcon(const char * const uplo, const uint32_t * const n, Complex16 *ap, int *ipiv, Complex16 *anorm, Complex16 *rcond, Complex16 *work, int *info);
void zhpev(const char * const jobz, const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, double *rwork, int *info);
void zhpevd(const char * const jobz, const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, const uint32_t * const lwork, double *rwork, const uint32_t * const lrwork, int *iwork, const uint32_t * const liwork, int *info);
void zhpevx(const char * const jobz, const char * const range, const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *vl, Complex16 *vu, uint32_t *il, uint32_t *iu, Complex16 *abstol, uint32_t *m, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, double *rwork, int *iwork, uint32_t *ifail, int *info);
void zhpgst(uint32_t *itype, const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *bp, int *info);
void zhpgv(uint32_t *itype, const char * const jobz, const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *bp, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, double *rwork, int *info);
void zhpgvd(uint32_t *itype, const char * const jobz, const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *bp, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, const uint32_t * const lwork, double *rwork, const uint32_t * const lrwork, int *iwork, const uint32_t * const liwork, int *info);
void zhpgvx(uint32_t *itype, const char * const jobz, const char * const range, const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *bp, Complex16 *vl, Complex16 *vu, uint32_t *il, uint32_t *iu, Complex16 *abstol, uint32_t *m, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, double *rwork, int *iwork, uint32_t *ifail, int *info);
void zhprfs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *ap, Complex16 *afp, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zhpsv(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *ap, int *ipiv, Complex16 *b, const uint32_t * const ldb, int *info);
void zhpsvx(const char * const fact, const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *ap, Complex16 *afp, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *rcond, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zhptrd(const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *d, Complex16 *e, Complex16 *tau, int *info);
void zhptrf(const char * const uplo, const uint32_t * const n, Complex16 *ap, int *ipiv, int *info);
void zhptri(const char * const uplo, const uint32_t * const n, Complex16 *ap, int *ipiv, Complex16 *work, int *info);
void zhptrs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *ap, int *ipiv, Complex16 *b, const uint32_t * const ldb, int *info);
void zhsein(const char * const side, const char * const eigsrc, const char * const initv, uint32_t *select, const uint32_t * const n, Complex16 *h, const uint32_t * const ldh, Complex16 *w, Complex16 *vl, const uint32_t * const ldvl, Complex16 *vr, const uint32_t * const ldvr, uint32_t *mm, uint32_t *m, Complex16 *work, double *rwork, uint32_t *ifaill, uint32_t *ifailr, int *info);
void zhseqr(const char * const job, const char * const compz, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, Complex16 *h, const uint32_t * const ldh, Complex16 *w, Complex16 *z, const uint32_t * const ldz, Complex16 *work, const uint32_t * const lwork, int *info);
void zlabrd(uint32_t *m, const uint32_t * const n, uint32_t *nb, Complex16 *a, const uint32_t * const lda, Complex16 *d, Complex16 *e, Complex16 *tauq, Complex16 *taup, Complex16 *x, const uint32_t * const ldx, Complex16 *y, const uint32_t * const ldy);
void zlacgv(const uint32_t * const n, Complex16 *x, const uint32_t * const incx);
void zlacn2(const uint32_t * const n, Complex16 *v, Complex16 *x, Complex16 *est, uint32_t *kase, uint32_t *isave);
void zlacon(const uint32_t * const n, Complex16 *v, Complex16 *x, Complex16 *est, uint32_t *kase);
void zlacp2(const char * const uplo, uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb);
void zlacpy(const char * const uplo, uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb);
void zlacrm(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *c, const uint32_t * const ldc, double *rwork);
void zlacrt(const uint32_t * const n, Complex16 *cx, const uint32_t * const incx, Complex16 *cy, const uint32_t * const incy, Complex16 *c, Complex16 *s);
void zladiv(Complex16 * retval, Complex16 *x, Complex16 *y);
void zlaed0(uint32_t *qsiz, const uint32_t * const n, Complex16 *d, Complex16 *e, Complex16 *q, const uint32_t * const ldq, Complex16 *qstore, const uint32_t * const ldqs, double *rwork, int *iwork, int *info);
void zlaed7(const uint32_t * const n, uint32_t *cutpnt, uint32_t *qsiz, uint32_t *tlvls, uint32_t *curlvl, uint32_t *curpbm, Complex16 *d, Complex16 *q, const uint32_t * const ldq, Complex16 *rho, uint32_t *indxq, Complex16 *qstore, uint32_t *qptr, uint32_t *prmptr, uint32_t *perm, uint32_t *givptr, uint32_t *givcol, Complex16 *givnum, Complex16 *work, double *rwork, int *iwork, int *info);
void zlaed8(uint32_t *k, const uint32_t * const n, uint32_t *qsiz, Complex16 *q, const uint32_t * const ldq, Complex16 *d, Complex16 *rho, uint32_t *cutpnt, Complex16 *z, Complex16 *dlamda, Complex16 *q2, const uint32_t * const ldq2, Complex16 *w, uint32_t *indxp, uint32_t *indx, uint32_t *indxq, uint32_t *perm, uint32_t *givptr, uint32_t *givcol, Complex16 *givnum, int *info);
void zlaein(uint32_t *rightv, uint32_t *noinit, const uint32_t * const n, Complex16 *h, const uint32_t * const ldh, Complex16 *w, Complex16 *v, Complex16 *b, const uint32_t * const ldb, double *rwork, Complex16 *eps3, Complex16 *smlnum, int *info);
void zlaesy(Complex16 *a, Complex16 *b, Complex16 *c, Complex16 *rt1, Complex16 *rt2, Complex16 *evscal, Complex16 *cs1, Complex16 *sn1);
void zlaev2(Complex16 *a, Complex16 *b, Complex16 *c, Complex16 *rt1, Complex16 *rt2, Complex16 *cs1, Complex16 *sn1);
void zlags2(uint32_t *upper, Complex16 *a1, Complex16 *a2, Complex16 *a3, Complex16 *b1, Complex16 *b2, Complex16 *b3, Complex16 *csu, Complex16 *snu, Complex16 *csv, Complex16 *snv, Complex16 *csq, Complex16 *snq);
void zlagtm(const char * const trans, const uint32_t * const n, uint32_t *nrhs, Complex16 *alpha, Complex16 *dl, Complex16 *d, Complex16 *du, Complex16 *x, const uint32_t * const ldx, Complex16 *beta, Complex16 *b, const uint32_t * const ldb);
void zlahef(const char * const uplo, const uint32_t * const n, uint32_t *nb, uint32_t *kb, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *w, const uint32_t * const ldw, int *info);
void zlahqr(uint32_t *wantt, uint32_t *wantz, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, Complex16 *h, const uint32_t * const ldh, Complex16 *w, uint32_t *iloz, uint32_t *ihiz, Complex16 *z, const uint32_t * const ldz, int *info);
void zlahr2(const uint32_t * const n, uint32_t *k, uint32_t *nb, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *t, const uint32_t * const ldt, Complex16 *y, const uint32_t * const ldy);
void zlahrd(const uint32_t * const n, uint32_t *k, uint32_t *nb, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *t, const uint32_t * const ldt, Complex16 *y, const uint32_t * const ldy);
void zlaic1(uint32_t *job, uint32_t *j, Complex16 *x, Complex16 *sest, Complex16 *w, Complex16 *gamma_, Complex16 *sestpr, Complex16 *s, Complex16 *c);
void zlals0(uint32_t *icompq, uint32_t *nl, uint32_t *nr, uint32_t *sqre, uint32_t *nrhs, Complex16 *b, const uint32_t * const ldb, Complex16 *bx, const uint32_t * const ldbx, uint32_t *perm, uint32_t *givptr, uint32_t *givcol, const uint32_t * const ldgcol, Complex16 *givnum, const uint32_t * const ldgnum, Complex16 *poles, Complex16 *difl, Complex16 *difr, Complex16 *z, uint32_t *k, Complex16 *c, Complex16 *s, double *rwork, int *info);
void zlalsa(uint32_t *icompq, uint32_t *smlsiz, const uint32_t * const n, uint32_t *nrhs, Complex16 *b, const uint32_t * const ldb, Complex16 *bx, const uint32_t * const ldbx, Complex16 *u, const uint32_t * const ldu, Complex16 *vt, uint32_t *k, Complex16 *difl, Complex16 *difr, Complex16 *z, Complex16 *poles, uint32_t *givptr, uint32_t *givcol, const uint32_t * const ldgcol, uint32_t *perm, Complex16 *givnum, Complex16 *c, Complex16 *s, double *rwork, int *iwork, int *info);
void zlalsd(const char * const uplo, uint32_t *smlsiz, const uint32_t * const n, uint32_t *nrhs, Complex16 *d, Complex16 *e, Complex16 *b, const uint32_t * const ldb, Complex16 *rcond, uint32_t *rank, Complex16 *work, double *rwork, int *iwork, int *info);
Complex16 zlangb(const char * const norm, const uint32_t * const n, uint32_t *kl, uint32_t *ku, Complex16 *ab, const uint32_t * const ldab, Complex16 *work);
Complex16 zlange(const char * const norm, uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *work);
Complex16 zlangt(const char * const norm, const uint32_t * const n, Complex16 *dl, Complex16 *d, Complex16 *du);
Complex16 zlanhb(const char * const norm, const char * const uplo, const uint32_t * const n, uint32_t *k, Complex16 *ab, const uint32_t * const ldab, Complex16 *work);
Complex16 zlanhe(const char * const norm, const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *work);
Complex16 zlanhf(const char * const norm, const char * const transr, const char * const uplo, const uint32_t * const n, Complex16 *a, Complex16 *work);
Complex16 zlanhp(const char * const norm, const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *work);
Complex16 zlanhs(const char * const norm, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *work);
Complex16 zlanht(const char * const norm, const uint32_t * const n, Complex16 *d, Complex16 *e);
Complex16 zlansb(const char * const norm, const char * const uplo, const uint32_t * const n, uint32_t *k, Complex16 *ab, const uint32_t * const ldab, Complex16 *work);
Complex16 zlansp(const char * const norm, const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *work);
Complex16 zlansy(const char * const norm, const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *work);
Complex16 zlantb(const char * const norm, const char * const uplo, const char * const diag, const uint32_t * const n, uint32_t *k, Complex16 *ab, const uint32_t * const ldab, Complex16 *work);
Complex16 zlantp(const char * const norm, const char * const uplo, const char * const diag, const uint32_t * const n, Complex16 *ap, Complex16 *work);
Complex16 zlantr(const char * const norm, const char * const uplo, const char * const diag, uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *work);
void zlapll(const uint32_t * const n, Complex16 *x, const uint32_t * const incx, Complex16 *y, const uint32_t * const incy, Complex16 *ssmin);
void zlapmr(uint32_t *forwrd, uint32_t *m, const uint32_t * const n, Complex16 *x, const uint32_t * const ldx, uint32_t *k);
void zlapmt(uint32_t *forwrd, uint32_t *m, const uint32_t * const n, Complex16 *x, const uint32_t * const ldx, uint32_t *k);
void zlaqgb(uint32_t *m, const uint32_t * const n, uint32_t *kl, uint32_t *ku, Complex16 *ab, const uint32_t * const ldab, Complex16 *r, Complex16 *c, Complex16 *rowcnd, Complex16 *colcnd, Complex16 *amax, const char * const equed);
void zlaqge(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *r, Complex16 *c, Complex16 *rowcnd, Complex16 *colcnd, Complex16 *amax, const char * const equed);
void zlaqhb(const char * const uplo, const uint32_t * const n, uint32_t *kd, Complex16 *ab, const uint32_t * const ldab, Complex16 *s, Complex16 *scond, Complex16 *amax, const char * const equed);
void zlaqhe(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *s, Complex16 *scond, Complex16 *amax, const char * const equed);
void zlaqhp(const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *s, Complex16 *scond, Complex16 *amax, const char * const equed);
void zlaqp2(uint32_t *m, const uint32_t * const n, uint32_t *offset, Complex16 *a, const uint32_t * const lda, uint32_t *jpvt, Complex16 *tau, Complex16 *vn1, Complex16 *vn2, Complex16 *work);
void zlaqps(uint32_t *m, const uint32_t * const n, uint32_t *offset, uint32_t *nb, uint32_t *kb, Complex16 *a, const uint32_t * const lda, uint32_t *jpvt, Complex16 *tau, Complex16 *vn1, Complex16 *vn2, Complex16 *auxv, Complex16 *f, const uint32_t * const ldf);
void zlaqr0(uint32_t *wantt, uint32_t *wantz, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, Complex16 *h, const uint32_t * const ldh, Complex16 *w, uint32_t *iloz, uint32_t *ihiz, Complex16 *z, const uint32_t * const ldz, Complex16 *work, const uint32_t * const lwork, int *info);
void zlaqr1(const uint32_t * const n, Complex16 *h, const uint32_t * const ldh, Complex16 *s1, Complex16 *s2, Complex16 *v);
void zlaqr2(uint32_t *wantt, uint32_t *wantz, const uint32_t * const n, uint32_t *ktop, uint32_t *kbot, uint32_t *nw, Complex16 *h, const uint32_t * const ldh, uint32_t *iloz, uint32_t *ihiz, Complex16 *z, const uint32_t * const ldz, uint32_t *ns, uint32_t *nd, Complex16 *sh, Complex16 *v, const uint32_t * const ldv, uint32_t *nh, Complex16 *t, const uint32_t * const ldt, uint32_t *nv, Complex16 *wv, const uint32_t * const ldwv, Complex16 *work, const uint32_t * const lwork);
void zlaqr3(uint32_t *wantt, uint32_t *wantz, const uint32_t * const n, uint32_t *ktop, uint32_t *kbot, uint32_t *nw, Complex16 *h, const uint32_t * const ldh, uint32_t *iloz, uint32_t *ihiz, Complex16 *z, const uint32_t * const ldz, uint32_t *ns, uint32_t *nd, Complex16 *sh, Complex16 *v, const uint32_t * const ldv, uint32_t *nh, Complex16 *t, const uint32_t * const ldt, uint32_t *nv, Complex16 *wv, const uint32_t * const ldwv, Complex16 *work, const uint32_t * const lwork);
void zlaqr4(uint32_t *wantt, uint32_t *wantz, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, Complex16 *h, const uint32_t * const ldh, Complex16 *w, uint32_t *iloz, uint32_t *ihiz, Complex16 *z, const uint32_t * const ldz, Complex16 *work, const uint32_t * const lwork, int *info);
void zlaqr5(uint32_t *wantt, uint32_t *wantz, uint32_t *kacc22, const uint32_t * const n, uint32_t *ktop, uint32_t *kbot, uint32_t *nshfts, Complex16 *s, Complex16 *h, const uint32_t * const ldh, uint32_t *iloz, uint32_t *ihiz, Complex16 *z, const uint32_t * const ldz, Complex16 *v, const uint32_t * const ldv, Complex16 *u, const uint32_t * const ldu, uint32_t *nv, Complex16 *wv, const uint32_t * const ldwv, uint32_t *nh, Complex16 *wh, const uint32_t * const ldwh);
void zlaqsb(const char * const uplo, const uint32_t * const n, uint32_t *kd, Complex16 *ab, const uint32_t * const ldab, Complex16 *s, Complex16 *scond, Complex16 *amax, const char * const equed);
void zlaqsp(const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *s, Complex16 *scond, Complex16 *amax, const char * const equed);
void zlaqsy(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *s, Complex16 *scond, Complex16 *amax, const char * const equed);
void zlar1v(const uint32_t * const n, uint32_t *b1, uint32_t *bn, Complex16 *lambda, Complex16 *d, Complex16 *l, Complex16 *ld, Complex16 *lld, Complex16 *pivmin, Complex16 *gaptol, Complex16 *z, uint32_t *wantnc, uint32_t *negcnt, Complex16 *ztz, Complex16 *mingma, uint32_t *r, uint32_t *isuppz, Complex16 *nrminv, Complex16 *resid, Complex16 *rqcorr, Complex16 *work);
void zlar2v(const uint32_t * const n, Complex16 *x, Complex16 *y, Complex16 *z, const uint32_t * const incx, Complex16 *c, Complex16 *s, const uint32_t * const incc);
void zlarcm(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *c, const uint32_t * const ldc, double *rwork);
void zlarf(const char * const side, uint32_t *m, const uint32_t * const n, Complex16 *v, const uint32_t * const incv, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work);
void zlarfb(const char * const side, const char * const trans, const char * const direct, const char * const storev, uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *v, const uint32_t * const ldv, Complex16 *t, const uint32_t * const ldt, Complex16 *c, const uint32_t * const ldc, Complex16 *work, const uint32_t * const ldwork);
void zlarfg(const uint32_t * const n, Complex16 *alpha, Complex16 *x, const uint32_t * const incx, Complex16 *tau);
void zlarfgp(const uint32_t * const n, Complex16 *alpha, Complex16 *x, const uint32_t * const incx, Complex16 *tau);
void zlarft(const char * const direct, const char * const storev, const uint32_t * const n, uint32_t *k, Complex16 *v, const uint32_t * const ldv, Complex16 *tau, Complex16 *t, const uint32_t * const ldt);
void zlarfx(const char * const side, uint32_t *m, const uint32_t * const n, Complex16 *v, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work);
void zlargv(const uint32_t * const n, Complex16 *x, const uint32_t * const incx, Complex16 *y, const uint32_t * const incy, Complex16 *c, const uint32_t * const incc);
void zlarnv(uint32_t *idist, uint32_t *iseed, const uint32_t * const n, Complex16 *x);
void zlarrv(const uint32_t * const n, Complex16 *vl, Complex16 *vu, Complex16 *d, Complex16 *l, Complex16 *pivmin, uint32_t *isplit, uint32_t *m, uint32_t *dol, uint32_t *dou, Complex16 *minrgp, Complex16 *rtol1, Complex16 *rtol2, Complex16 *w, Complex16 *werr, Complex16 *wgap, uint32_t *iblock, uint32_t *indexw, Complex16 *gers, Complex16 *z, const uint32_t * const ldz, uint32_t *isuppz, Complex16 *work, int *iwork, int *info);
void zlartg(Complex16 *f, Complex16 *g, Complex16 *cs, Complex16 *sn, Complex16 *r);
void zlartv(const uint32_t * const n, Complex16 *x, const uint32_t * const incx, Complex16 *y, const uint32_t * const incy, Complex16 *c, Complex16 *s, const uint32_t * const incc);
void zlarz(const char * const side, uint32_t *m, const uint32_t * const n, uint32_t *l, Complex16 *v, const uint32_t * const incv, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work);
void zlarzb(const char * const side, const char * const trans, const char * const direct, const char * const storev, uint32_t *m, const uint32_t * const n, uint32_t *k, uint32_t *l, Complex16 *v, const uint32_t * const ldv, Complex16 *t, const uint32_t * const ldt, Complex16 *c, const uint32_t * const ldc, Complex16 *work, const uint32_t * const ldwork);
void zlarzt(const char * const direct, const char * const storev, const uint32_t * const n, uint32_t *k, Complex16 *v, const uint32_t * const ldv, Complex16 *tau, Complex16 *t, const uint32_t * const ldt);
void zlascl(const char * const type, uint32_t *kl, uint32_t *ku, Complex16 *cfrom, Complex16 *cto, uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *info);
void zlaset(const char * const uplo, uint32_t *m, const uint32_t * const n, Complex16 *alpha, Complex16 *beta, Complex16 *a, const uint32_t * const lda);
void zlasr(const char * const side, const char * const pivot, const char * const direct, uint32_t *m, const uint32_t * const n, Complex16 *c, Complex16 *s, Complex16 *a, const uint32_t * const lda);
void zlassq(const uint32_t * const n, Complex16 *x, const uint32_t * const incx, Complex16 *scale, Complex16 *sumsq);
void zlaswp(const uint32_t * const n, Complex16 *a, const uint32_t * const lda, uint32_t *k1, uint32_t *k2, int *ipiv, const uint32_t * const incx);
void zlasyf(const char * const uplo, const uint32_t * const n, uint32_t *nb, uint32_t *kb, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *w, const uint32_t * const ldw, int *info);
void zlatbs(const char * const uplo, const char * const trans, const char * const diag, const char * const normin, const uint32_t * const n, uint32_t *kd, Complex16 *ab, const uint32_t * const ldab, Complex16 *x, Complex16 *scale, Complex16 *cnorm, int *info);
void zlatdf(uint32_t *ijob, const uint32_t * const n, Complex16 *z, const uint32_t * const ldz, Complex16 *rhs, Complex16 *rdsum, Complex16 *rdscal, int *ipiv, int *jpiv);
void zlatps(const char * const uplo, const char * const trans, const char * const diag, const char * const normin, const uint32_t * const n, Complex16 *ap, Complex16 *x, Complex16 *scale, Complex16 *cnorm, int *info);
void zlatrd(const char * const uplo, const uint32_t * const n, uint32_t *nb, Complex16 *a, const uint32_t * const lda, Complex16 *e, Complex16 *tau, Complex16 *w, const uint32_t * const ldw);
void zlatrs(const char * const uplo, const char * const trans, const char * const diag, const char * const normin, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *x, Complex16 *scale, Complex16 *cnorm, int *info);
void zlatrz(uint32_t *m, const uint32_t * const n, uint32_t *l, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work);
void zlatzm(const char * const side, uint32_t *m, const uint32_t * const n, Complex16 *v, const uint32_t * const incv, Complex16 *tau, Complex16 *c1, Complex16 *c2, const uint32_t * const ldc, Complex16 *work);
void zlauu2(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *info);
void zlauum(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *info);
void zpbcon(const char * const uplo, const uint32_t * const n, uint32_t *kd, Complex16 *ab, const uint32_t * const ldab, Complex16 *anorm, Complex16 *rcond, Complex16 *work, double *rwork, int *info);
void zpbequ(const char * const uplo, const uint32_t * const n, uint32_t *kd, Complex16 *ab, const uint32_t * const ldab, Complex16 *s, Complex16 *scond, Complex16 *amax, int *info);
void zpbrfs(const char * const uplo, const uint32_t * const n, uint32_t *kd, uint32_t *nrhs, Complex16 *ab, const uint32_t * const ldab, Complex16 *afb, const uint32_t * const ldafb, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zpbstf(const char * const uplo, const uint32_t * const n, uint32_t *kd, Complex16 *ab, const uint32_t * const ldab, int *info);
void zpbsv(const char * const uplo, const uint32_t * const n, uint32_t *kd, uint32_t *nrhs, Complex16 *ab, const uint32_t * const ldab, Complex16 *b, const uint32_t * const ldb, int *info);
void zpbsvx(const char * const fact, const char * const uplo, const uint32_t * const n, uint32_t *kd, uint32_t *nrhs, Complex16 *ab, const uint32_t * const ldab, Complex16 *afb, const uint32_t * const ldafb, const char * const equed, Complex16 *s, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *rcond, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zpbtf2(const char * const uplo, const uint32_t * const n, uint32_t *kd, Complex16 *ab, const uint32_t * const ldab, int *info);
void zpbtrf(const char * const uplo, const uint32_t * const n, uint32_t *kd, Complex16 *ab, const uint32_t * const ldab, int *info);
void zpbtrs(const char * const uplo, const uint32_t * const n, uint32_t *kd, uint32_t *nrhs, Complex16 *ab, const uint32_t * const ldab, Complex16 *b, const uint32_t * const ldb, int *info);
void zpftrf(const char * const transr, const char * const uplo, const uint32_t * const n, Complex16 *a, int *info);
void zpftri(const char * const transr, const char * const uplo, const uint32_t * const n, Complex16 *a, int *info);
void zpftrs(const char * const transr, const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, Complex16 *b, const uint32_t * const ldb, int *info);
void zpocon(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *anorm, Complex16 *rcond, Complex16 *work, double *rwork, int *info);
void zpoequ(const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *s, Complex16 *scond, Complex16 *amax, int *info);
void zpoequb(const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *s, Complex16 *scond, Complex16 *amax, int *info);
void zporfs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *af, const uint32_t * const ldaf, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zposv(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, int *info);
void zposvx(const char * const fact, const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *af, const uint32_t * const ldaf, const char * const equed, Complex16 *s, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *rcond, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zpotf2(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *info);
void zpotrf(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *info);
void zpotri(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *info);
void zpotrs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, int *info);
void zppcon(const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *anorm, Complex16 *rcond, Complex16 *work, double *rwork, int *info);
void zppequ(const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *s, Complex16 *scond, Complex16 *amax, int *info);
void zpprfs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *ap, Complex16 *afp, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zppsv(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *ap, Complex16 *b, const uint32_t * const ldb, int *info);
void zppsvx(const char * const fact, const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *ap, Complex16 *afp, const char * const equed, Complex16 *s, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *rcond, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zpptrf(const char * const uplo, const uint32_t * const n, Complex16 *ap, int *info);
void zpptri(const char * const uplo, const uint32_t * const n, Complex16 *ap, int *info);
void zpptrs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *ap, Complex16 *b, const uint32_t * const ldb, int *info);
void zpstf2(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, uint32_t *piv, uint32_t *rank, Complex16 *tol, Complex16 *work, int *info);
void zpstrf(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, uint32_t *piv, uint32_t *rank, Complex16 *tol, Complex16 *work, int *info);
void zptcon(const uint32_t * const n, Complex16 *d, Complex16 *e, Complex16 *anorm, Complex16 *rcond, double *rwork, int *info);
void zpteqr(const char * const compz, const uint32_t * const n, Complex16 *d, Complex16 *e, Complex16 *z, const uint32_t * const ldz, Complex16 *work, int *info);
void zptrfs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *d, Complex16 *e, Complex16 *df, Complex16 *ef, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zptsv(const uint32_t * const n, uint32_t *nrhs, Complex16 *d, Complex16 *e, Complex16 *b, const uint32_t * const ldb, int *info);
void zptsvx(const char * const fact, const uint32_t * const n, uint32_t *nrhs, Complex16 *d, Complex16 *e, Complex16 *df, Complex16 *ef, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *rcond, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zpttrf(const uint32_t * const n, Complex16 *d, Complex16 *e, int *info);
void zpttrs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *d, Complex16 *e, Complex16 *b, const uint32_t * const ldb, int *info);
void zptts2(uint32_t *iuplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *d, Complex16 *e, Complex16 *b, const uint32_t * const ldb);
void zrot(const uint32_t * const n, Complex16 *cx, const uint32_t * const incx, Complex16 *cy, const uint32_t * const incy, Complex16 *c, Complex16 *s);
void zspcon(const char * const uplo, const uint32_t * const n, Complex16 *ap, int *ipiv, Complex16 *anorm, Complex16 *rcond, Complex16 *work, int *info);
void zspmv(const char * const uplo, const uint32_t * const n, Complex16 *alpha, Complex16 *ap, Complex16 *x, const uint32_t * const incx, Complex16 *beta, Complex16 *y, const uint32_t * const incy);
void zspr(const char * const uplo, const uint32_t * const n, Complex16 *alpha, Complex16 *x, const uint32_t * const incx, Complex16 *ap);
void zsprfs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *ap, Complex16 *afp, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zspsv(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *ap, int *ipiv, Complex16 *b, const uint32_t * const ldb, int *info);
void zspsvx(const char * const fact, const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *ap, Complex16 *afp, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *rcond, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zsptrf(const char * const uplo, const uint32_t * const n, Complex16 *ap, int *ipiv, int *info);
void zsptri(const char * const uplo, const uint32_t * const n, Complex16 *ap, int *ipiv, Complex16 *work, int *info);
void zsptrs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *ap, int *ipiv, Complex16 *b, const uint32_t * const ldb, int *info);
void zstedc(const char * const compz, const uint32_t * const n, Complex16 *d, Complex16 *e, Complex16 *z, const uint32_t * const ldz, Complex16 *work, const uint32_t * const lwork, double *rwork, const uint32_t * const lrwork, int *iwork, const uint32_t * const liwork, int *info);
void zstegr(const char * const jobz, const char * const range, const uint32_t * const n, Complex16 *d, Complex16 *e, Complex16 *vl, Complex16 *vu, uint32_t *il, uint32_t *iu, Complex16 *abstol, uint32_t *m, Complex16 *w, Complex16 *z, const uint32_t * const ldz, uint32_t *isuppz, Complex16 *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void zstein(const uint32_t * const n, Complex16 *d, Complex16 *e, uint32_t *m, Complex16 *w, uint32_t *iblock, uint32_t *isplit, Complex16 *z, const uint32_t * const ldz, Complex16 *work, int *iwork, uint32_t *ifail, int *info);
void zstemr(const char * const jobz, const char * const range, const uint32_t * const n, Complex16 *d, Complex16 *e, Complex16 *vl, Complex16 *vu, uint32_t *il, uint32_t *iu, uint32_t *m, Complex16 *w, Complex16 *z, const uint32_t * const ldz, uint32_t *nzc, uint32_t *isuppz, uint32_t *tryrac, Complex16 *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void zsteqr(const char * const compz, const uint32_t * const n, Complex16 *d, Complex16 *e, Complex16 *z, const uint32_t * const ldz, Complex16 *work, int *info);
void zsycon(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *anorm, Complex16 *rcond, Complex16 *work, int *info);
void zsyconv(const char * const uplo, const char * const way, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *work, int *info);
void zsyequb(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *s, Complex16 *scond, Complex16 *amax, Complex16 *work, int *info);
void zsymv(const char * const uplo, const uint32_t * const n, Complex16 *alpha, Complex16 *a, const uint32_t * const lda, Complex16 *x, const uint32_t * const incx, Complex16 *beta, Complex16 *y, const uint32_t * const incy);
void zsyr(const char * const uplo, const uint32_t * const n, Complex16 *alpha, Complex16 *x, const uint32_t * const incx, Complex16 *a, const uint32_t * const lda);
void zsyrfs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *af, const uint32_t * const ldaf, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void zsysv(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *work, const uint32_t * const lwork, int *info);
void zsysvx(const char * const fact, const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *af, const uint32_t * const ldaf, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *rcond, Complex16 *ferr, Complex16 *berr, Complex16 *work, const uint32_t * const lwork, double *rwork, int *info);
void zsyswapr(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, uint32_t *i1, uint32_t *i2);
void zsytf2(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, int *info);
void zsytrf(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *work, const uint32_t * const lwork, int *info);
void zsytri(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *work, int *info);
void zsytri2(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *work, const uint32_t * const lwork, int *info);
void zsytri2x(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *work, uint32_t *nb, int *info);
void zsytrs(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *b, const uint32_t * const ldb, int *info);
void zsytrs2(const char * const uplo, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, int *ipiv, Complex16 *b, const uint32_t * const ldb, Complex16 *work, int *info);
void ztbcon(const char * const norm, const char * const uplo, const char * const diag, const uint32_t * const n, uint32_t *kd, Complex16 *ab, const uint32_t * const ldab, Complex16 *rcond, Complex16 *work, double *rwork, int *info);
void ztbrfs(const char * const uplo, const char * const trans, const char * const diag, const uint32_t * const n, uint32_t *kd, uint32_t *nrhs, Complex16 *ab, const uint32_t * const ldab, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void ztbtrs(const char * const uplo, const char * const trans, const char * const diag, const uint32_t * const n, uint32_t *kd, uint32_t *nrhs, Complex16 *ab, const uint32_t * const ldab, Complex16 *b, const uint32_t * const ldb, int *info);
void ztfsm(const char * const transr, const char * const side, const char * const uplo, const char * const trans, const char * const diag, uint32_t *m, const uint32_t * const n, Complex16 *alpha, Complex16 *a, Complex16 *b, const uint32_t * const ldb);
void ztftri(const char * const transr, const char * const uplo, const char * const diag, const uint32_t * const n, Complex16 *a, int *info);
void ztfttp(const char * const transr, const char * const uplo, const uint32_t * const n, Complex16 *arf, Complex16 *ap, int *info);
void ztfttr(const char * const transr, const char * const uplo, const uint32_t * const n, Complex16 *arf, Complex16 *a, const uint32_t * const lda, int *info);
void ztgevc(const char * const side, const char * const howmny, uint32_t *select, const uint32_t * const n, Complex16 *s, const uint32_t * const lds, Complex16 *p, const uint32_t * const ldp, Complex16 *vl, const uint32_t * const ldvl, Complex16 *vr, const uint32_t * const ldvr, uint32_t *mm, uint32_t *m, Complex16 *work, double *rwork, int *info);
void ztgex2(uint32_t *wantq, uint32_t *wantz, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *q, const uint32_t * const ldq, Complex16 *z, const uint32_t * const ldz, uint32_t *j1_, int *info);
void ztgexc(uint32_t *wantq, uint32_t *wantz, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *q, const uint32_t * const ldq, Complex16 *z, const uint32_t * const ldz, uint32_t *ifst, uint32_t *ilst, int *info);
void ztgsen(uint32_t *ijob, uint32_t *wantq, uint32_t *wantz, uint32_t *select, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *alpha, Complex16 *beta, Complex16 *q, const uint32_t * const ldq, Complex16 *z, const uint32_t * const ldz, uint32_t *m, Complex16 *pl, Complex16 *pr, Complex16 *dif, Complex16 *work, const uint32_t * const lwork, int *iwork, const uint32_t * const liwork, int *info);
void ztgsja(const char * const jobu, const char * const jobv, const char * const jobq, uint32_t *m, uint32_t *p, const uint32_t * const n, uint32_t *k, uint32_t *l, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *tola, Complex16 *tolb, Complex16 *alpha, Complex16 *beta, Complex16 *u, const uint32_t * const ldu, Complex16 *v, const uint32_t * const ldv, Complex16 *q, const uint32_t * const ldq, Complex16 *work, uint32_t *ncycle, int *info);
void ztgsna(const char * const job, const char * const howmny, uint32_t *select, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *vl, const uint32_t * const ldvl, Complex16 *vr, const uint32_t * const ldvr, Complex16 *s, Complex16 *dif, uint32_t *mm, uint32_t *m, Complex16 *work, const uint32_t * const lwork, int *iwork, int *info);
void ztgsy2(const char * const trans, uint32_t *ijob, uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *c, const uint32_t * const ldc, Complex16 *d, const uint32_t * const ldd, Complex16 *e, const uint32_t * const lde, Complex16 *f, const uint32_t * const ldf, Complex16 *scale, Complex16 *rdsum, Complex16 *rdscal, int *info);
void ztgsyl(const char * const trans, uint32_t *ijob, uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *c, const uint32_t * const ldc, Complex16 *d, const uint32_t * const ldd, Complex16 *e, const uint32_t * const lde, Complex16 *f, const uint32_t * const ldf, Complex16 *scale, Complex16 *dif, Complex16 *work, const uint32_t * const lwork, int *iwork, int *info);
void ztpcon(const char * const norm, const char * const uplo, const char * const diag, const uint32_t * const n, Complex16 *ap, Complex16 *rcond, Complex16 *work, double *rwork, int *info);
void ztpmqrt(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, uint32_t *l, uint32_t *nb, Complex16 *v, const uint32_t * const ldv, Complex16 *t, const uint32_t * const ldt, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *work, int *info);
void ztpqrt(uint32_t *m, const uint32_t * const n, uint32_t *l, uint32_t *nb, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *t, const uint32_t * const ldt, Complex16 *work, int *info);
void ztpqrt2(uint32_t *m, const uint32_t * const n, uint32_t *l, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *t, const uint32_t * const ldt, int *info);
void ztprfb(const char * const side, const char * const trans, const char * const direct, const char * const storev, uint32_t *m, const uint32_t * const n, uint32_t *k, uint32_t *l, Complex16 *v, const uint32_t * const ldv, Complex16 *t, const uint32_t * const ldt, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *work, const uint32_t * const ldwork);
void ztprfs(const char * const uplo, const char * const trans, const char * const diag, const uint32_t * const n, uint32_t *nrhs, Complex16 *ap, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void ztptri(const char * const uplo, const char * const diag, const uint32_t * const n, Complex16 *ap, int *info);
void ztptrs(const char * const uplo, const char * const trans, const char * const diag, const uint32_t * const n, uint32_t *nrhs, Complex16 *ap, Complex16 *b, const uint32_t * const ldb, int *info);
void ztpttf(const char * const transr, const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *arf, int *info);
void ztpttr(const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *a, const uint32_t * const lda, int *info);
void ztrcon(const char * const norm, const char * const uplo, const char * const diag, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *rcond, Complex16 *work, double *rwork, int *info);
void ztrevc(const char * const side, const char * const howmny, uint32_t *select, const uint32_t * const n, Complex16 *t, const uint32_t * const ldt, Complex16 *vl, const uint32_t * const ldvl, Complex16 *vr, const uint32_t * const ldvr, uint32_t *mm, uint32_t *m, Complex16 *work, double *rwork, int *info);
void ztrexc(const char * const compq, const uint32_t * const n, Complex16 *t, const uint32_t * const ldt, Complex16 *q, const uint32_t * const ldq, uint32_t *ifst, uint32_t *ilst, int *info);
void ztrrfs(const char * const uplo, const char * const trans, const char * const diag, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *x, const uint32_t * const ldx, Complex16 *ferr, Complex16 *berr, Complex16 *work, double *rwork, int *info);
void ztrsen(const char * const job, const char * const compq, uint32_t *select, const uint32_t * const n, Complex16 *t, const uint32_t * const ldt, Complex16 *q, const uint32_t * const ldq, Complex16 *w, uint32_t *m, Complex16 *s, Complex16 *sep, Complex16 *work, const uint32_t * const lwork, int *info);
void ztrsna(const char * const job, const char * const howmny, uint32_t *select, const uint32_t * const n, Complex16 *t, const uint32_t * const ldt, Complex16 *vl, const uint32_t * const ldvl, Complex16 *vr, const uint32_t * const ldvr, Complex16 *s, Complex16 *sep, uint32_t *mm, uint32_t *m, Complex16 *work, const uint32_t * const ldwork, double *rwork, int *info);
void ztrsyl(const char * const trana, const char * const tranb, uint32_t *isgn, uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, Complex16 *c, const uint32_t * const ldc, Complex16 *scale, int *info);
void ztrti2(const char * const uplo, const char * const diag, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *info);
void ztrtri(const char * const uplo, const char * const diag, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, int *info);
void ztrtrs(const char * const uplo, const char * const trans, const char * const diag, const uint32_t * const n, uint32_t *nrhs, Complex16 *a, const uint32_t * const lda, Complex16 *b, const uint32_t * const ldb, int *info);
void ztrttf(const char * const transr, const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *arf, int *info);
void ztrttp(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *ap, int *info);
void ztzrqf(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *tau, int *info);
void ztzrzf(uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zunbdb(const char * const trans, const char * const signs, uint32_t *m, uint32_t *p, uint32_t *q, Complex16 *x11, const uint32_t * const ldx11, Complex16 *x12, const uint32_t * const ldx12, Complex16 *x21, const uint32_t * const ldx21, Complex16 *x22, const uint32_t * const ldx22, Complex16 *theta, Complex16 *phi, Complex16 *taup1, Complex16 *taup2, Complex16 *tauq1, Complex16 *tauq2, Complex16 *work, const uint32_t * const lwork, int *info);
void zuncsd(const char * const jobu1, const char * const jobu2, const char * const jobv1t, const char * const jobv2t, const char * const trans, const char * const signs, uint32_t *m, uint32_t *p, uint32_t *q, Complex16 *x11, const uint32_t * const ldx11, Complex16 *x12, const uint32_t * const ldx12, Complex16 *x21, const uint32_t * const ldx21, Complex16 *x22, const uint32_t * const ldx22, Complex16 *theta, Complex16 *u1, const uint32_t * const ldu1, Complex16 *u2, const uint32_t * const ldu2, Complex16 *v1t, const uint32_t * const ldv1t, Complex16 *v2t, const uint32_t * const ldv2t, Complex16 *work, const uint32_t * const lwork, double *rwork, const uint32_t * const lrwork, int *iwork, int *info);
void zung2l(uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, int *info);
void zung2r(uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, int *info);
void zungbr(const char * const vect, uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zunghr(const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zungl2(uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, int *info);
void zunglq(uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zungql(uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zungqr(uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zungr2(uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, int *info);
void zungrq(uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zungtr(const char * const uplo, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *work, const uint32_t * const lwork, int *info);
void zunm2l(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work, int *info);
void zunm2r(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work, int *info);
void zunmbr(const char * const vect, const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work, const uint32_t * const lwork, int *info);
void zunmhr(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *ilo, uint32_t *ihi, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work, const uint32_t * const lwork, int *info);
void zunml2(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work, int *info);
void zunmlq(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work, const uint32_t * const lwork, int *info);
void zunmql(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work, const uint32_t * const lwork, int *info);
void zunmqr(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work, const uint32_t * const lwork, int *info);
void zunmr2(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work, int *info);
void zunmr3(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, uint32_t *l, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work, int *info);
void zunmrq(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work, const uint32_t * const lwork, int *info);
void zunmrz(const char * const side, const char * const trans, uint32_t *m, const uint32_t * const n, uint32_t *k, uint32_t *l, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work, const uint32_t * const lwork, int *info);
void zunmtr(const char * const side, const char * const uplo, const char * const trans, uint32_t *m, const uint32_t * const n, Complex16 *a, const uint32_t * const lda, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work, const uint32_t * const lwork, int *info);
void zupgtr(const char * const uplo, const uint32_t * const n, Complex16 *ap, Complex16 *tau, Complex16 *q, const uint32_t * const ldq, Complex16 *work, int *info);
void zupmtr(const char * const side, const char * const uplo, const char * const trans, uint32_t *m, const uint32_t * const n, Complex16 *ap, Complex16 *tau, Complex16 *c, const uint32_t * const ldc, Complex16 *work, int *info);
}
