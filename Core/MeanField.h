/* This file was automatically generated.  Do not edit! */
uint32_t MeanFieldGroundstate(const OperatorType *const O,CompactMatrixType *const Rho,const double kbT,double *Etot);
uint32_t CalculateDensityMatrixFromOneElectronEigensystem(const CompactMatrixType *const fun,double *val,CompactMatrixType *Rho,const double kbT,const double Nelec,double *Etot);
uint32_t MeanFieldOperator(const OperatorType *const O,OperatorType *const MO,const CompactMatrixType Rho);
uint32_t FreeOperator(OperatorType *const H);
uint32_t QDetNotOrderedFcFaBcBa(const uint16_t A,const uint16_t B,const uint32_t NFermion);
uint32_t NormalOrderOperator(const OperatorType O,OperatorType *const Oo,uint32_t(*QDetNotOrdered)(uint16_t,uint16_t,uint32_t));
uint32_t InitOperatorCopySettings(OperatorType *const O,const OperatorType Oset);
uint32_t MeanFieldOperatorOneElectronMatrix(const OperatorType *const O,CompactMatrixType *const MO,const CompactMatrixType Rho);
uint32_t OperatorCrAnLength2ToOneElectronCompactMatrix(const OperatorType *const O,CompactMatrixType *const Oone);
uint32_t OperatorToOneElectronMatrix(const OperatorType *const O,MatrixType *const Oone);
uint32_t RealOperatorAddLader(OperatorType *const O,const uint16_t *const oadd,const uint32_t length,const double VRe);
uint32_t ComplexOperatorAddLader(OperatorType *const O,const uint16_t *const oadd,const uint32_t length,const double VRe,const double VIm);
uint32_t MakeRealOperatorFixedLengthComplex(OperatorType *const O,const uint32_t length);
uint32_t InitOperatorFixedLength(OperatorType *const O,const uint32_t length,const uint32_t NBitsKey);
uint32_t InitOperator(OperatorType *const O);
uint32_t CompactMatrixToOneParticleOperator(const CompactMatrixType *const Omat,OperatorType *const O);
