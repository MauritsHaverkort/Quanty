#!/bin/sh
# The DFT/HF add-on and the atomic code it builds on compute in double precision only
# (Doc/Claude/QuantyDFT_StageA, design section 0.1): no float and no single-precision
# math functions. Run from the repository root, as CI does. Lists the offending lines
# and exits 1 when it finds any.
if grep -rnwE "float|sqrtf|expf|powf|logf|fabsf|erff|sinf|cosf|tanf|atanf" Core/BasicMath AtomicOperators DFT Martin; then
  echo "single precision found in the atomic code"
  exit 1
fi
echo "passed"
