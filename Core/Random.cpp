#include "Random.h"

double mwhdrand()
{
  return (double)rand()/((double)(RAND_MAX)+1);
}

void mwhsrand(unsigned int seed)
{
  srand(seed);
  return;
}
