
uint32_t BADoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType *O,BADoubleTriDiagonalMatrixType Mup,BADoubleTriDiagonalMatrixType Mdn,double U,double mu);
uint32_t DoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType *O,DoubleTriDiagonalMatrixType Mup,DoubleTriDiagonalMatrixType Mdn,double U,double mu);
uint32_t SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType *O,SecondDoubleTriDiagonalMatrixType Mup,SecondDoubleTriDiagonalMatrixType Mdn,double U,double mu);
uint32_t AndersonMatrixToAndersonImpurityOperator(OperatorType *O,AndersonMatrixType Mup,AndersonMatrixType Mdn,double U,double mu);
uint32_t QDetNotOrderedFcFaBcBa(const uint16_t A,const uint16_t B,const uint32_t NFermion);
uint32_t RealOperatorAddLaderNormalOrder(OperatorType *const O,const uint16_t *const oadd,const uint32_t length,const double VRe,uint32_t(*QDetNotOrdered)(uint16_t,uint16_t,uint32_t));
uint32_t InitOperatorFixedLength(OperatorType *const O,const uint32_t length,const uint32_t NBitsKey);
uint32_t InitNamedOperator(OperatorType *const O,char *name);
uint32_t TriDiagonalMatrixToAndersonImpurityOperator(OperatorType *O,TriDiagonalMatrixType Mup,TriDiagonalMatrixType Mdn,double U,double mu);
