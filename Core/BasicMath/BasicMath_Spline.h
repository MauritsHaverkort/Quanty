
int InitCubicSpline(CubicSplineType *S);
void FreeCubicSpline(CubicSplineType *S);
void CreateCubicSpline(CubicSplineType *S, double yp1, double ypn);
void CreateCubicSpline(double *x, double *y, uint32_t n, double yp1, double yp2, double *y2);
double EvaluateCubicSpline(CubicSplineType *S, double x);
double EvaluateCubicSpline(double *xa, double *ya, double *y2a, uint32_t n, double x);
