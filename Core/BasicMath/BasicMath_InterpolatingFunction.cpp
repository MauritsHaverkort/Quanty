int InitInterpolatingFunctionAsCubicSpline(InterpolatingFunctionType *fun, uint32_t N)
{
  fun->type='C';
  fun->f = (void*)malloc(sizeof(CubicSplineType));
  ((CubicSplineType*)(fun->f))->N=N;
  InitCubicSpline((CubicSplineType*)(fun->f));
  return 0;
}


void FreeInterpolatingFunction(InterpolatingFunctionType *fun)
{
  switch (fun->type)
  {
    case 'C':
    {
      FreeCubicSpline((CubicSplineType*)fun->f);
      fun->type='\0';
      break;
    }
    default:
    {
      printf("ERROR in Free Interpolating Function. Unrecognized interpolation type \"%c\"",fun->type);
    }
  }
}

double EvaluateInterpolatingFunction(InterpolatingFunctionType *fun, double x)
{
  switch (fun->type)
  {
    case 'C':
    {
      return EvaluateCubicSpline((CubicSplineType*)fun->f, x);
    }
    default:
    {
      printf("ERROR in Evaluate Interpolating Function. Unrecognized interpolation type \"%c\"",fun->type);
      return NAN;
    }
  }
}
