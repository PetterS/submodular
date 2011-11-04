#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif

#define INSTANTIATE(c) \
	template class Petter::PseudoBoolean<c>; \
	template class Petter::SymmetricPseudoBoolean<c>; 

INSTANTIATE(double);
INSTANTIATE(int);
