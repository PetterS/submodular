#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif

#define INSTANTIATE(c) \
	template class Petter::BipartiteVertexPacking<c>; \
	template class Petter::VertexPacking<c>;

INSTANTIATE(double);
INSTANTIATE(int);
