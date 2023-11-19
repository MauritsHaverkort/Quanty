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
#include "postscript.cpp"
#include "Wannier90.cpp"


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
