printf("====    written by Maurits W. Haverkort                  ====\n");
uint32_t InitializeTiming(TimeType * T);
uint32_t EnlargeTimingStructure(TimeType * T);
void FreeTiming(TimeType * T);
uint32_t BeginTiming(const char *s, TimeType * T);
uint32_t EndTiming(const char *s, TimeType * T);
uint32_t StreamPrintTiming(TimeType T, FILE * stream);
uint32_t FilePrintTiming(TimeType T, const char *const flname, const char *const mode);
uint32_t PrintTiming(TimeType T);
