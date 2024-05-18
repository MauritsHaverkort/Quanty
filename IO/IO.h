#include "BitMap.h"
#include "DresdenFPLO.h"
#include "FilePrint.h"
#include "FileRead.h"
#include "Print.h"
#include "PrintMathematicaExpression.h"
#include "Read.h"
#include "StreamPrint.h"
#include "StreamRead.h"
#include "StringPrint.h"
#include "StuttgartNMTO.h"
#include "postscript.h"
#include "Wannier90.h"

int EncodeStringBase64(const DynamicString *in, DynamicString *out);
void RemoveSpacesAndNewlines(char *str);
void RemoveDoubleSpacesAndNewlines(char *str);
void RemoveCommentLines(char *str, char comment_char);
