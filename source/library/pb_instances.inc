#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif

#define INSTANTIATE(c) \
	template class Petter::PseudoBoolean<c>; \
	template class Petter::SymmetricPseudoBoolean<c>; \
	template class Petter::GeneratorPseudoBoolean<c>; \
	template class Petter::Generators<c>; 

INSTANTIATE(double);
INSTANTIATE(int);
