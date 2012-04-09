//#ifndef SUBMODULAR_ASD234
//#define SUBMODULAR_ASD234

void test();

#include <vector>
#include <PseudoBoolean.h>

enum Method {
    GRD,
    GRD_gen,
    HOCR
};

class PseudoBoolean : public Petter::PseudoBoolean<double>
{
    public:
        PseudoBoolean();
        void add_monomial(double c, int i, int j=-1, int k=-1, int l=-1);
        int nvars();
        double eval(signed char* x, int nx);
        
        char* __str__();
        
        double minimize(Method method, signed char* x, int nx);
};

//#endif
