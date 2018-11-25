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
std::vector<double> InterpolatingFunction::operator()(const std::vector<double> &t)
{
	printf("error: call of virtual InterpolatingFunction::operator()\n");
	return std::vector<double>();
}
std::vector<double> InterpolatingFunction::operator()(const std::vector<double> &t, unsigned int k)
{
	printf("error: call of virtual InterpolatingFunction::operator()\n");
	return std::vector<double>();
}

void InterpolatingFunction::Create(const std::vector<double> &u, const std::vector<double> &v)
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

double InterpolatingFunction::GetX(uint32_t i)
{
	return x[i];
}

double InterpolatingFunction::GetY(uint32_t i)
{
	return y[i];
}

} /* namespace Interpolation */
