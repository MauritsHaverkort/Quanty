/*  
 *  u = M v     (see MatrixVector files)
 *  u = v^* M v (see MatrixVector files)
 *  u = v
 *  u = alpha v
 *  alpha = u v
 *  u = u + alpha v
 *  
 *
 */

void RealVectorCopy(double *const u, double *const v, uint32_t N)
{
  uint32_t i;
  for (i = 0; i < N; i++)
  {
    u[i] = v[i];
  }
}

void ComplexVectorCopy(double *const uRe, double *const uIm, double *const vRe, double *const vIm, uint32_t N)
{
  uint32_t i;
  for (i = 0; i < N; i++)
  {
    uRe[i] = vRe[i];
    uIm[i] = vIm[i];
  }
}

void RealVectorMultiply(double *const u, const uint32_t N, const double alphaRe)
{
  uint32_t i;
  for (i = 0; i < N; i++)
  {
    u[i] *= alphaRe;
  }
}

void ComplexVectorMultiply(double *const uRe, double *const uIm, const uint32_t N, const double alphaRe,
                           const double alphaIm)
{
  uint32_t i;
  double temp;
  for (i = 0; i < N; i++)
  {
    temp = uRe[i];
    uRe[i] = alphaRe * uRe[i] - alphaIm * uIm[i];
    uIm[i] = alphaRe * uIm[i] + alphaIm * temp;
  }
}

/*
 *
 */
void RealVectorDotVectorMahan(const double *const u, const double *const v, uint32_t N, double *const a)
{
  uint32_t i;
  volatile double c, y, t;
  *a = 0;
  c = 0.0;
  y = 0.0;
  t = 0.0;
  for (i = 0; i < N; i++)
  {
    y = u[i] * v[i] - c;
    t = *a + y;
    c = (t - *a);
    c = c - y;
    *a = t;
  }
}


void RealVectorDotVector(const double *const u, const double *const v, uint32_t N, double *const a)
{
  int const one=1;
  *a = ddot(&N, u, &one, v, &one);
}


void ComplexVectorDotVector(const double *const uRe, const double *const uIm, const double *const vRe,
                            const double *const vIm, const uint32_t N, double *const aRe, double *const aIm)
{
  /* Warning u is complex conjugated */
  uint32_t i;
  *aRe = 0;
  *aIm = 0;
  for (i = 0; i < N; i++)
  {
    *aRe += uRe[i] * vRe[i] + uIm[i] * vIm[i];
    *aIm += uRe[i] * vIm[i] - uIm[i] * vRe[i];
  }
}

void RealVectorAdd(double *const u, double *const v, uint32_t N, double alphaRe)
{
  uint32_t i;
  for (i = 0; i < N; i++)
  {
    u[i] += alphaRe * v[i];
  }
}

void ComplexVectorAdd(double *const uRe, double *const uIm, double *const vRe, double *const vIm, const uint32_t N,
                      const double alphaRe, const double alphaIm)
{
  uint32_t i;
  for (i = 0; i < N; i++)
  {
    uRe[i] += alphaRe * vRe[i] - alphaIm * vIm[i];
    uIm[i] += alphaRe * vIm[i] + alphaIm * vRe[i];
  }
}
