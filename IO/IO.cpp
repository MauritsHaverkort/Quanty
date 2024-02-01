#include "BitMap.cpp"
#include "DresdenFPLO.cpp"
#include "FilePrint.cpp"
#include "FileRead.cpp"
#include "Print.cpp"
#include "PrintMathematicaExpression.cpp"
#include "Read.cpp"
#include "StreamPrint.cpp"
#include "StreamRead.cpp"
#include "StringPrint.cpp"
#include "StuttgartNMTO.cpp"
//#include "postscript.cpp"
#include "Wannier90.cpp"

int EncodeStringBase64(const DynamicString *in, DynamicString *out)
{
  int i,j;
  uint32_t a,b,c, abc;
  out->length = 4 * ((in->length + 2) / 3);
  out->allocated = out->length+1;
  out->str = (char *)malloc(out->allocated*sizeof(char));
  out->str[out->length] = '\0';
  if (out->str == NULL) return 1;
  j=0;
  for (i = 0; i < in->length;)
  {
    a = i < in->length ? (unsigned char)in->str[i++] : 0;
    b = i < in->length ? (unsigned char)in->str[i++] : 0;
    c = i < in->length ? (unsigned char)in->str[i++] : 0;
    abc = (a << 0x10) + (b << 0x08) + c;
    out->str[j++] = BASE64_CHARS[(abc >> 3 * 6) & 0x3F];
    out->str[j++] = BASE64_CHARS[(abc >> 2 * 6) & 0x3F];
    out->str[j++] = BASE64_CHARS[(abc >> 1 * 6) & 0x3F];
    out->str[j++] = BASE64_CHARS[(abc >> 0 * 6) & 0x3F];
  }
  switch(in->length % 3)
  {
    case 1:
      // One leftover byte -> two padding characters
      out->str[out->length - 2] = '=';
    case 2:
      // Two leftover bytes -> one padding character
      out->str[out->length - 1] = '=';
  }
  return 0;
}


class ExpectationValues
{
  //  Constructors:
public:
  ExpectationValues(const vector<WaveFunctionType> &psiList, OperatorType hamil, const vector<OperatorType> &opList, const map<string, vector<double> > &extraLists, bool hasHamiltonian, double zero =1E-10, bool sortByEnergy = false, bool calcSigma = false)
  {
    this->hasHamiltonian = hasHamiltonian;
    for (uint i = 0; i < (uint) psiList.size(); ++i)
    {
      values.push_back(vector<double>());
    }
    vector<double> energies;
    if (hasHamiltonian)
    {
      headers.push_back("E");
      energies = expectationValuesSingleOperator(psiList, hamil, zero);
      for (uint i = 0; i < (uint) values.size(); ++i)
      {
        values[i].push_back(energies[i]);
      }
    }
    for (uint j = 0; j < (uint) opList.size(); ++j)
    {
      headers.push_back(opList[j].Name);
      vector<double> tmp = expectationValuesSingleOperator(psiList, opList[j], zero);
      for (uint i = 0; i < (uint) values.size(); ++i)
      {
        values[i].push_back(tmp[i]);
      }
    }
    for (map<string, vector<double> >::const_iterator iter = extraLists.begin(); iter != extraLists.end(); ++iter)
    {
      headers.push_back(iter->first);
      if (iter->second.size() != psiList.size())
        cout << "Warning: Size " << iter->second.size() << " of List with name " << iter->first << " is not equal to " << psiList.size() << ", the number of Wavefunctions." << endl;
      for (uint i = 0; i < (uint) iter->second.size(); ++i)
        values[i].push_back(iter->second[i]);
    }
    if (hasHamiltonian && calcSigma)
    {
      headers.push_back("std-dev");
      vector<double> tmp = sigmaValues(psiList, hamil, energies);
      for (uint i = 0; i < (uint) values.size(); ++i)
        values[i].push_back(tmp[i]);
    }
    if (hasHamiltonian && sortByEnergy)
      sort(values.begin(), values.end(), matrixSorterByFirstValue);
  }
  
  //  Member Variables:
  vector<string> headers;
  vector<vector<double> > values;
  bool hasHamiltonian;
  
  //  Methods:
public:
  
  void print(uint colWidth = 8, uint printHeaderEvery = 0)
  {
    for (uint i = 0; i < (uint) values.size(); ++i)
    {
      if (i == 0 || (printHeaderEvery > 0 && i % printHeaderEvery == 0))
      {
        printHeader(colWidth);
      }
      printLine(i, colWidth);
    }
  }
  
private:
  
  void printHeader(int colWidth)
  {
    cout << endl;
    cout << blankFill(colWidth);
    cout << " ";
    for (uint i = 0; i < (uint) headers.size(); ++i)
    {
      string name(headers[i]);
      if (name.length() > colWidth - 1)
        name = name.substr(0, colWidth - 1);
      cout << name;
      cout << blankFill(colWidth, name.length());
    }
    cout << endl;
  }
  
  void printLine(uint i, uint colWidth)
  {
    string num = numberToStringMaxWidth(i + 1, colWidth);
    cout << num;
    cout << blankFill(colWidth - (int)num.length());
    for (uint j = 0; j < (uint) values[i].size(); ++j)
    {
      num = numberToStringMaxWidth(values[i][j], colWidth - 1);
      cout << num;
      cout << blankFill(colWidth, num.length());
    }
    cout << endl;
  }
};
