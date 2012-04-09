
import numpy as np;
import submodular
from submodular import test
from submodular import PseudoBoolean

test()


# f(x) = x0 - 2x1 - x2  - x0x1 + 4x1x2
f = PseudoBoolean()
f.add_monomial(1.0, 0)
f.add_monomial(-2.0, 1)
f.add_monomial(-1.0, 2)
f.add_monomial(-1.0, 0,1)
f.add_monomial(4.0, 1,2)

print f
print 'f(0,1,1) = ', f.eval(np.array([0,1,1], dtype=np.int8) )

method = submodular.GRD
n = f.nvars()
x = np.zeros( n, dtype=np.int8)
bound = f.minimize(method, x)

print "bound=", bound
print "x=", x
