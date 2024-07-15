
void FreeCompactTensor(CompactTensorType * const M);
void FreeCompactTensorIncludingPointer(CompactTensorType * const M);
uint32_t InitCompactTensor(CompactTensorType * const M);
int CompactTensorToString(CompactTensorType * const M, DynamicString *Ad);
int CompactTensorToString(CompactTensorType * const M, char **A);
void MakeCompactComplexTensorReal(CompactTensorType * M);
int CompactTensorCopy(CompactTensorType * A, CompactTensorType * copyA);
int CompactTensorChop(CompactTensorType * A, double epsilon);





