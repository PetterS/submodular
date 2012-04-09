
#include <iostream>
#include <sstream>
#include "submodular.h"

void test()
{
	std::cout << "Submodular library!" << std::endl;
}


PseudoBoolean::PseudoBoolean()
{
}

void PseudoBoolean::add_monomial(double c, int i, int j, int k, int l)
{
    if (l>=0) {
        Petter::PseudoBoolean<double>::add_monomial(i,j,k,l,c);
    }
    else if (k>=0) {
        Petter::PseudoBoolean<double>::add_monomial(i,j,k,c);
    }
    else if (j>=0) {
        Petter::PseudoBoolean<double>::add_monomial(i,j,c);
    }
    else {
        Petter::PseudoBoolean<double>::add_monomial(i,c);
    }
}

int PseudoBoolean::nvars()
{
    return Petter::PseudoBoolean<double>::nvars();
}

double PseudoBoolean::eval(signed char* x, int nx)
{
    std::vector<Petter::label> xvec(x, x+nx);
    return Petter::PseudoBoolean<double>::eval(xvec);
}

char* PseudoBoolean::__str__()
{
    static char my_str[1024];
    std::stringstream sout;
    Petter::PseudoBoolean<double>::print_helper(sout);
    std::strncpy(my_str, sout.str().c_str(), 1023);
    for (size_t i=0;i<1024;++i) {
        if (my_str[i] == '\n') my_str[i] = ' ';
    }
    return my_str;
}

double PseudoBoolean::minimize(Method method, signed char* x, int nx) 
{
    size_t n = nvars();
    std::vector<Petter::label> xvec(n,-1);
    double bound;
    int labeled;

    if (method == GRD) {
        bound = Petter::PseudoBoolean<double>::minimize(xvec, labeled, Petter::GRD, 0);
    }
    else if (method == HOCR) {
        bound = Petter::PseudoBoolean<double>::minimize(xvec, Petter::HOCR);
    }
    else {
        throw std::runtime_error("Unknown method");
    }
    
    for (size_t i=0; i<n && i<nx && i<xvec.size(); ++i) {
        x[i] = xvec[i];
    }
    
    return bound;
}

