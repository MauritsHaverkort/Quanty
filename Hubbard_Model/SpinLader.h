uint32_t CreateSpinLaderSpinOperatorsSz(SpinLaderSettings Set, OperatorType * S);
uint32_t CreateSpinLaderSpinOperatorsSplus(SpinLaderSettings Set, OperatorType * S);
uint32_t CreateSpinLaderSpinOperatorsSminus(SpinLaderSettings Set, OperatorType * S);
uint32_t SetValuesOfUnnededOperatorsToZero(OperatorType * H);
uint32_t CreateSpinLaderHamiltonian(SpinLaderSettings Set, OperatorType * H);
uint32_t CreateSpinLaderStartWaveFunction(const SpinLaderSettings Set, WaveFunctionType * psi);
uint32_t CreateSpinLaderSpinFlipOperator(SpinLaderSettings Set, OperatorType * O);
int main(int argc, const char *argv[]);
