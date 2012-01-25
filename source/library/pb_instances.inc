#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif

#define INSTANTIATE(c) \
	template class Petter::PseudoBoolean<c>; \
	template class Petter::SymmetricPseudoBoolean<c>; \
	template class Petter::GeneratorPseudoBoolean<c>; \
	template class Petter::Generators<c>; \
	template c     Petter::branch_and_bound<c>(const Petter::PseudoBoolean<c>& f, vector<label>& x, Petter::BBInfo* bbinfo);

INSTANTIATE(double);
INSTANTIATE(int);
