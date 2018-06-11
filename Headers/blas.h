
#ifdef BLASunderscore
#define FORTRAN_WRAPPER(x) x ## _

#define caxpy FORTRAN_WRAPPER(caxpy)
#define ccopy FORTRAN_WRAPPER(ccopy)
#define cdotc FORTRAN_WRAPPER(cdotc)
#define cdotu FORTRAN_WRAPPER(cdotu)
#define cgbmv FORTRAN_WRAPPER(cgbmv)
#define cgemm FORTRAN_WRAPPER(cgemm)
#define cgemv FORTRAN_WRAPPER(cgemv)
#define cgerc FORTRAN_WRAPPER(cgerc)
#define cgeru FORTRAN_WRAPPER(cgeru)
#define chbmv FORTRAN_WRAPPER(chbmv)
#define chemm FORTRAN_WRAPPER(chemm)
#define chemv FORTRAN_WRAPPER(chemv)
#define cher FORTRAN_WRAPPER(cher)
#define cher2 FORTRAN_WRAPPER(cher2)
#define cher2k FORTRAN_WRAPPER(cher2k)
#define cherk FORTRAN_WRAPPER(cherk)
#define chpmv FORTRAN_WRAPPER(chpmv)
#define chpr FORTRAN_WRAPPER(chpr)
#define chpr2 FORTRAN_WRAPPER(chpr2)
#define crotg FORTRAN_WRAPPER(crotg)
#define cscal FORTRAN_WRAPPER(cscal)
#define csrot FORTRAN_WRAPPER(csrot)
#define csscal FORTRAN_WRAPPER(csscal)
#define cswap FORTRAN_WRAPPER(cswap)
#define csymm FORTRAN_WRAPPER(csymm)
#define csyr2k FORTRAN_WRAPPER(csyr2k)
#define csyrk FORTRAN_WRAPPER(csyrk)
#define ctbmv FORTRAN_WRAPPER(ctbmv)
#define ctbsv FORTRAN_WRAPPER(ctbsv)
#define ctpmv FORTRAN_WRAPPER(ctpmv)
#define ctpsv FORTRAN_WRAPPER(ctpsv)
#define ctrmm FORTRAN_WRAPPER(ctrmm)
#define ctrmv FORTRAN_WRAPPER(ctrmv)
#define ctrsm FORTRAN_WRAPPER(ctrsm)
#define ctrsv FORTRAN_WRAPPER(ctrsv)
#define dasum FORTRAN_WRAPPER(dasum)
#define daxpy FORTRAN_WRAPPER(daxpy)
#define dcabs1 FORTRAN_WRAPPER(dcabs1)
#define dcopy FORTRAN_WRAPPER(dcopy)
#define ddot FORTRAN_WRAPPER(ddot)
#define dgbmv FORTRAN_WRAPPER(dgbmv)
#define dgemm FORTRAN_WRAPPER(dgemm)
#define dgemv FORTRAN_WRAPPER(dgemv)
#define dger FORTRAN_WRAPPER(dger)
#define dnrm2 FORTRAN_WRAPPER(dnrm2)
#define drot FORTRAN_WRAPPER(drot)
#define drotg FORTRAN_WRAPPER(drotg)
#define drotm FORTRAN_WRAPPER(drotm)
#define drotmg FORTRAN_WRAPPER(drotmg)
#define dsbmv FORTRAN_WRAPPER(dsbmv)
#define dscal FORTRAN_WRAPPER(dscal)
#define dsdot FORTRAN_WRAPPER(dsdot)
#define dspmv FORTRAN_WRAPPER(dspmv)
#define dspr FORTRAN_WRAPPER(dspr)
#define dspr2 FORTRAN_WRAPPER(dspr2)
#define dswap FORTRAN_WRAPPER(dswap)
#define dsymm FORTRAN_WRAPPER(dsymm)
#define dsymv FORTRAN_WRAPPER(dsymv)
#define dsyr FORTRAN_WRAPPER(dsyr)
#define dsyr2 FORTRAN_WRAPPER(dsyr2)
#define dsyr2k FORTRAN_WRAPPER(dsyr2k)
#define dsyrk FORTRAN_WRAPPER(dsyrk)
#define dtbmv FORTRAN_WRAPPER(dtbmv)
#define dtbsv FORTRAN_WRAPPER(dtbsv)
#define dtpmv FORTRAN_WRAPPER(dtpmv)
#define dtpsv FORTRAN_WRAPPER(dtpsv)
#define dtrmm FORTRAN_WRAPPER(dtrmm)
#define dtrmv FORTRAN_WRAPPER(dtrmv)
#define dtrsm FORTRAN_WRAPPER(dtrsm)
#define dtrsv FORTRAN_WRAPPER(dtrsv)
#define dzasum FORTRAN_WRAPPER(dzasum)
#define dznrm2 FORTRAN_WRAPPER(dznrm2)
#define icamax FORTRAN_WRAPPER(icamax)
#define idamax FORTRAN_WRAPPER(idamax)
#define isamax FORTRAN_WRAPPER(isamax)
#define izamax FORTRAN_WRAPPER(izamax)
#define lsame FORTRAN_WRAPPER(lsame)
#define sasum FORTRAN_WRAPPER(sasum)
#define saxpy FORTRAN_WRAPPER(saxpy)
#define scasum FORTRAN_WRAPPER(scasum)
#define scnrm2 FORTRAN_WRAPPER(scnrm2)
#define scopy FORTRAN_WRAPPER(scopy)
#define sdot FORTRAN_WRAPPER(sdot)
#define sdsdot FORTRAN_WRAPPER(sdsdot)
#define sgbmv FORTRAN_WRAPPER(sgbmv)
#define sgemm FORTRAN_WRAPPER(sgemm)
#define sgemv FORTRAN_WRAPPER(sgemv)
#define sger FORTRAN_WRAPPER(sger)
#define snrm2 FORTRAN_WRAPPER(snrm2)
#define srot FORTRAN_WRAPPER(srot)
#define srotg FORTRAN_WRAPPER(srotg)
#define srotm FORTRAN_WRAPPER(srotm)
#define srotmg FORTRAN_WRAPPER(srotmg)
#define ssbmv FORTRAN_WRAPPER(ssbmv)
#define sscal FORTRAN_WRAPPER(sscal)
#define sspmv FORTRAN_WRAPPER(sspmv)
#define sspr FORTRAN_WRAPPER(sspr)
#define sspr2 FORTRAN_WRAPPER(sspr2)
#define sswap FORTRAN_WRAPPER(sswap)
#define ssymm FORTRAN_WRAPPER(ssymm)
#define ssymv FORTRAN_WRAPPER(ssymv)
#define ssyr FORTRAN_WRAPPER(ssyr)
#define ssyr2 FORTRAN_WRAPPER(ssyr2)
#define ssyr2k FORTRAN_WRAPPER(ssyr2k)
#define ssyrk FORTRAN_WRAPPER(ssyrk)
#define stbmv FORTRAN_WRAPPER(stbmv)
#define stbsv FORTRAN_WRAPPER(stbsv)
#define stpmv FORTRAN_WRAPPER(stpmv)
#define stpsv FORTRAN_WRAPPER(stpsv)
#define strmm FORTRAN_WRAPPER(strmm)
#define strmv FORTRAN_WRAPPER(strmv)
#define strsm FORTRAN_WRAPPER(strsm)
#define strsv FORTRAN_WRAPPER(strsv)
#define xerbla FORTRAN_WRAPPER(xerbla)
#define zaxpy FORTRAN_WRAPPER(zaxpy)
#define zcopy FORTRAN_WRAPPER(zcopy)
#define zdotc FORTRAN_WRAPPER(zdotc)
#define zdotu FORTRAN_WRAPPER(zdotu)
#define zdrot FORTRAN_WRAPPER(zdrot)
#define zdscal FORTRAN_WRAPPER(zdscal)
#define zgbmv FORTRAN_WRAPPER(zgbmv)
#define zgemm FORTRAN_WRAPPER(zgemm)
#define zgemv FORTRAN_WRAPPER(zgemv)
#define zgerc FORTRAN_WRAPPER(zgerc)
#define zgeru FORTRAN_WRAPPER(zgeru)
#define zhbmv FORTRAN_WRAPPER(zhbmv)
#define zhemm FORTRAN_WRAPPER(zhemm)
#define zhemv FORTRAN_WRAPPER(zhemv)
#define zher FORTRAN_WRAPPER(zher)
#define zher2 FORTRAN_WRAPPER(zher2)
#define zher2k FORTRAN_WRAPPER(zher2k)
#define zherk FORTRAN_WRAPPER(zherk)
#define zhpmv FORTRAN_WRAPPER(zhpmv)
#define zhpr FORTRAN_WRAPPER(zhpr)
#define zhpr2 FORTRAN_WRAPPER(zhpr2)
#define zrotg FORTRAN_WRAPPER(zrotg)
#define zscal FORTRAN_WRAPPER(zscal)
#define zswap FORTRAN_WRAPPER(zswap)
#define zsymm FORTRAN_WRAPPER(zsymm)
#define zsyr2k FORTRAN_WRAPPER(zsyr2k)
#define zsyrk FORTRAN_WRAPPER(zsyrk)
#define ztbmv FORTRAN_WRAPPER(ztbmv)
#define ztbsv FORTRAN_WRAPPER(ztbsv)
#define ztpmv FORTRAN_WRAPPER(ztpmv)
#define ztpsv FORTRAN_WRAPPER(ztpsv)
#define ztrmm FORTRAN_WRAPPER(ztrmm)
#define ztrmv FORTRAN_WRAPPER(ztrmv)
#define ztrsm FORTRAN_WRAPPER(ztrsm)
#define ztrsv FORTRAN_WRAPPER(ztrsv)

#endif

double dasum(const uint32_t * const n,double *dx,const uint32_t * const incx);
void daxpy(const uint32_t * const n,double *da,double *dx,const uint32_t * const incx,double *dy,const uint32_t * const incy);
double dcabs1(double *z);
void dcopy(const uint32_t * const n,double *dx,const uint32_t * const incx,double *dy,const uint32_t * const incy);
double ddot(const uint32_t * const n,const double *const dx,const uint32_t * const incx,const double *const dy,const uint32_t *const incy);
void dgbmv(const char * const trans,const uint32_t * const m,const uint32_t * const n,uint32_t *kl,uint32_t *ku,double *alpha,double *a,const uint32_t * const lda,double *x,const uint32_t * const incx,double *beta,double *y,const uint32_t * const incy);
void dgemm(const char * const transa,const char * const transb,const uint32_t * const m,const uint32_t * const n,const uint32_t * const k,double *alpha,double *a,const uint32_t * const lda,double *b,const uint32_t * const ldb,double *beta,double *c,const uint32_t * const ldc);
void dgemv(const char * const trans,const uint32_t * const m,const uint32_t * const n,double *alpha,double *a,const uint32_t * const lda,double *x,const uint32_t * const incx,double *beta,double *y,const uint32_t * const incy);
void dger(const uint32_t * const m,const uint32_t * const n,double *alpha,double *x,const uint32_t * const incx,double *y,const uint32_t * const incy,double *a,const uint32_t * const lda);
double dnrm2(const uint32_t * const n,double *x,const uint32_t * const incx);
void drot(const uint32_t * const n,double *dx,const uint32_t * const incx,double *dy,const uint32_t * const incy,double *c,double *s);
void drotg(double *da,double *db,double *c,double *s);
void drotm(const uint32_t * const n,double *dx,const uint32_t * const incx,double *dy,const uint32_t * const incy,double *dparam);
void drotmg(double *dd1,double *dd2,double *dx1,double *dy1,double *dparam);
void dsbmv(const char * const uplo,const uint32_t * const n,uint32_t *k,double *alpha,double *a,const uint32_t * const lda,double *x,const uint32_t * const incx,double *beta,double *y,const uint32_t * const incy);
void dscal(const uint32_t * const n,double *da,double *dx,const uint32_t * const incx);
void dspmv(const char * const uplo,const uint32_t * const n,double *alpha,double *ap,double *x,const uint32_t * const incx,double *beta,double *y,const uint32_t * const incy);
void dspr(const char * const uplo,const uint32_t * const n,double *alpha,double *x,const uint32_t * const incx,double *ap);
void dspr2(const char * const uplo,const uint32_t * const n,double *alpha,double *x,const uint32_t * const incx,double *y,const uint32_t * const incy,double *ap);
void dswap(const uint32_t * const n,double *dx,const uint32_t * const incx,double *dy,const uint32_t * const incy);
void dsymm(const char * const side,const char * const uplo,const uint32_t * const m,const uint32_t * const n,double *alpha,double *a,const uint32_t * const lda,double *b,const uint32_t * const ldb,double *beta,double *c,const uint32_t * const ldc);
void dsymv(const char * const uplo,const uint32_t * const n, const double * const alpha,const double * const a,const uint32_t * const lda,double *x,const uint32_t * const incx, const double * const beta,double *y,const uint32_t * const incy);
void dsyr(const char * const uplo,const uint32_t * const n,double *alpha,double *x,const uint32_t * const incx,double *a,const uint32_t * const lda);
void dsyr2(const char * const uplo,const uint32_t * const n,double *alpha,double *x,const uint32_t * const incx,double *y,const uint32_t * const incy,double *a,const uint32_t * const lda);
void dsyr2k(const char * const uplo,const char * const trans,const uint32_t * const n,uint32_t *k,double *alpha,double *a,const uint32_t * const lda,double *b,const uint32_t * const ldb,double *beta,double *c,const uint32_t * const ldc);
void dsyrk(const char * const uplo,const char * const trans,const uint32_t * const n,uint32_t *k,double *alpha,double *a,const uint32_t * const lda,double *beta,double *c,const uint32_t * const ldc);
void dtbmv(const char * const uplo,const char * const trans,const char * const diag,const uint32_t * const n,uint32_t *k,double *a,const uint32_t * const lda,double *x,const uint32_t * const incx);
void dtbsv(const char * const uplo,const char * const trans,const char * const diag,const uint32_t * const n,uint32_t *k,double *a,const uint32_t * const lda,double *x,const uint32_t * const incx);
void dtpmv(const char * const uplo,const char * const trans,const char * const diag,const uint32_t * const n,double *ap,double *x,const uint32_t * const incx);
void dtpsv(const char * const uplo,const char * const trans,const char * const diag,const uint32_t * const n,double *ap,double *x,const uint32_t * const incx);
void dtrmm(const char * const side,const char * const uplo,const char * const transa,const char * const diag,const uint32_t * const m,const uint32_t * const n,double *alpha,double *a,const uint32_t * const lda,double *b,const uint32_t * const ldb);
void dtrmv(const char * const uplo,const char * const trans,const char * const diag,const uint32_t * const n,double *a,const uint32_t * const lda,double *x,const uint32_t * const incx);
void dtrsm(const char * const side,const char * const uplo,const char * const transa,const char * const diag,const uint32_t * const m,const uint32_t * const n,double *alpha,double *a,const uint32_t * const lda,double *b,const uint32_t * const ldb);
void dtrsv(const char * const uplo,const char * const trans,const char * const diag,const uint32_t * const n,double *a,const uint32_t * const lda,double *x,const uint32_t * const incx);
double dzasum(const uint32_t * const n,double *zx,const uint32_t * const incx);
double dznrm2(const uint32_t * const n,double *x,const uint32_t * const incx);
int idamax(const uint32_t * const n,double *dx,const uint32_t * const incx);
int izamax(const uint32_t * const n,double *zx,const uint32_t * const incx);
int lsame(const char * const ca,const char * const cb);
void xerbla(const char * const srname,uint32_t *info);
void zaxpy(const uint32_t * const n,Complex16 *za,Complex16 *zx,const uint32_t * const incx,Complex16 *zy,const uint32_t * const incy);
void zcopy(const uint32_t * const n,Complex16 *zx,const uint32_t * const incx,Complex16 *zy,const uint32_t * const incy);
void zdotc(Complex16 * retval, const uint32_t * const n,Complex16 *zx,const uint32_t * const incx,Complex16 *zy,const uint32_t * const incy);
void zdotu(Complex16 * retval, const uint32_t * const n,Complex16 *zx,const uint32_t * const incx,Complex16 *zy,const uint32_t * const incy);
void zdrot(const uint32_t * const n,Complex16 *cx,const uint32_t * const incx,Complex16 *cy,const uint32_t * const incy,Complex16 *c,Complex16 *s);
void zdscal(const uint32_t * const n,Complex16 *da,Complex16 *zx,const uint32_t * const incx);
void zgbmv(const char * const trans,const uint32_t * const m,const uint32_t * const n,uint32_t *kl,uint32_t *ku,Complex16 *alpha,Complex16 *a,const uint32_t * const lda,Complex16 *x,const uint32_t * const incx,Complex16 *beta,Complex16 *y,const uint32_t * const incy);
void zgemm(const char * const transa,const char * const transb,const uint32_t * const m,const uint32_t * const n,const uint32_t * const k,Complex16 *alpha,Complex16 *a,const uint32_t * const lda,Complex16 *b,const uint32_t * const ldb,Complex16 *beta,Complex16 *c,const uint32_t * const ldc);
void zgemv(const char * const trans,const uint32_t * const m,const uint32_t * const n,Complex16 *alpha,Complex16 *a,const uint32_t * const lda,Complex16 *x,const uint32_t * const incx,Complex16 *beta,Complex16 *y,const uint32_t * const incy);
void zgerc(const uint32_t * const m,const uint32_t * const n,Complex16 *alpha,Complex16 *x,const uint32_t * const incx,Complex16 *y,const uint32_t * const incy,Complex16 *a,const uint32_t * const lda);
void zgeru(const uint32_t * const m,const uint32_t * const n,Complex16 *alpha,Complex16 *x,const uint32_t * const incx,Complex16 *y,const uint32_t * const incy,Complex16 *a,const uint32_t * const lda);
void zhbmv(const char * const uplo,const uint32_t * const n,uint32_t *k,Complex16 *alpha,Complex16 *a,const uint32_t * const lda,Complex16 *x,const uint32_t * const incx,Complex16 *beta,Complex16 *y,const uint32_t * const incy);
void zhemm(const char * const side,const char * const uplo,const uint32_t * const m,const uint32_t * const n,Complex16 *alpha,Complex16 *a,const uint32_t * const lda,Complex16 *b,const uint32_t * const ldb,Complex16 *beta,Complex16 *c,const uint32_t * const ldc);
void zhemv(const char * const uplo,const uint32_t * const n,Complex16 *alpha,Complex16 *a,const uint32_t * const lda,Complex16 *x,const uint32_t * const incx,Complex16 *beta,Complex16 *y,const uint32_t * const incy);
void zher(const char * const uplo,const uint32_t * const n,Complex16 *alpha,Complex16 *x,const uint32_t * const incx,Complex16 *a,const uint32_t * const lda);
void zher2(const char * const uplo,const uint32_t * const n,Complex16 *alpha,Complex16 *x,const uint32_t * const incx,Complex16 *y,const uint32_t * const incy,Complex16 *a,const uint32_t * const lda);
void zher2k(const char * const uplo,const char * const trans,const uint32_t * const n,uint32_t *k,Complex16 *alpha,Complex16 *a,const uint32_t * const lda,Complex16 *b,const uint32_t * const ldb,Complex16 *beta,Complex16 *c,const uint32_t * const ldc);
void zherk(const char * const uplo,const char * const trans,const uint32_t * const n,uint32_t *k,Complex16 *alpha,Complex16 *a,const uint32_t * const lda,Complex16 *beta,Complex16 *c,const uint32_t * const ldc);
void zhpmv(const char * const uplo,const uint32_t * const n,Complex16 *alpha,Complex16 *ap,Complex16 *x,const uint32_t * const incx,Complex16 *beta,Complex16 *y,const uint32_t * const incy);
void zhpr(const char * const uplo,const uint32_t * const n,Complex16 *alpha,Complex16 *x,const uint32_t * const incx,Complex16 *ap);
void zhpr2(const char * const uplo,const uint32_t * const n,Complex16 *alpha,Complex16 *x,const uint32_t * const incx,Complex16 *y,const uint32_t * const incy,Complex16 *ap);
void zrotg(Complex16 *ca,Complex16 *cb,Complex16 *c,Complex16 *s);
void zscal(const uint32_t * const n,Complex16 *za,Complex16 *zx,const uint32_t * const incx);
void zswap(const uint32_t * const n,Complex16 *zx,const uint32_t * const incx,Complex16 *zy,const uint32_t * const incy);
void zsymm(const char * const side,const char * const uplo,const uint32_t * const m,const uint32_t * const n,Complex16 *alpha,Complex16 *a,const uint32_t * const lda,Complex16 *b,const uint32_t * const ldb,Complex16 *beta,Complex16 *c,const uint32_t * const ldc);
void zsyr2k(const char * const uplo,const char * const trans,const uint32_t * const n,uint32_t *k,Complex16 *alpha,Complex16 *a,const uint32_t * const lda,Complex16 *b,const uint32_t * const ldb,Complex16 *beta,Complex16 *c,const uint32_t * const ldc);
void zsyrk(const char * const uplo,const char * const trans,const uint32_t * const n,uint32_t *k,Complex16 *alpha,Complex16 *a,const uint32_t * const lda,Complex16 *beta,Complex16 *c,const uint32_t * const ldc);
void ztbmv(const char * const uplo,const char * const trans,const char * const diag,const uint32_t * const n,uint32_t *k,Complex16 *a,const uint32_t * const lda,Complex16 *x,const uint32_t * const incx);
void ztbsv(const char * const uplo,const char * const trans,const char * const diag,const uint32_t * const n,uint32_t *k,Complex16 *a,const uint32_t * const lda,Complex16 *x,const uint32_t * const incx);
void ztpmv(const char * const uplo,const char * const trans,const char * const diag,const uint32_t * const n,Complex16 *ap,Complex16 *x,const uint32_t * const incx);
void ztpsv(const char * const uplo,const char * const trans,const char * const diag,const uint32_t * const n,Complex16 *ap,Complex16 *x,const uint32_t * const incx);
void ztrmm(const char * const side,const char * const uplo,const char * const transa,const char * const diag,const uint32_t * const m,const uint32_t * const n,Complex16 *alpha,Complex16 *a,const uint32_t * const lda,Complex16 *b,const uint32_t * const ldb);
void ztrmv(const char * const uplo,const char * const trans,const char * const diag,const uint32_t * const n,Complex16 *a,const uint32_t * const lda,Complex16 *x,const uint32_t * const incx);
void ztrsm(const char * const side,const char * const uplo,const char * const transa,const char * const diag,const uint32_t * const m,const uint32_t * const n,Complex16 *alpha,Complex16 *a,const uint32_t * const lda,Complex16 *b,const uint32_t * const ldb);
void ztrsv(const char * const uplo,const char * const trans,const char * const diag,const uint32_t * const n,Complex16 *a,const uint32_t * const lda,Complex16 *x,const uint32_t * const incx);
