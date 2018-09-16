/*
 * InterpolatingFunction.cpp
 *
 *  Created on: Sep 8, 2018
 *      Author: brass
 */

namespace Interpolation {

InterpolatingFunction::InterpolatingFunction() {
	// TODO Auto-generated constructor stub

}

InterpolatingFunction::~InterpolatingFunction() {
	// TODO Auto-generated destructor stub
}

double InterpolatingFunction::operator()(double t)
{
	printf("error: call of virtual InterpolatingFunction::operator()\n");
	return NAN;
}
double InterpolatingFunction::operator()(double t, unsigned int n)
{
	printf("error: call of virtual InterpolatingFunction::operator()\n");
	return NAN;
}
std::vector<double> InterpolatingFunction::operator()(std::vector<double> t)
{
	printf("error: call of virtual InterpolatingFunction::operator()\n");
	return std::vector<double>();
}

void InterpolatingFunction::Create(std::vector<double>, std::vector<double>)
{

}

std::vector<double> InterpolatingFunction::GetX()
{
	return x;
}

std::vector<double> InterpolatingFunction::GetY()
{
	return y;
}

} /* namespace Interpolation */
