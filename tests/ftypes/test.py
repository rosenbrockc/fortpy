from celib.symmetry_module import make_primitive, get_spacegroup
from numpy import array, reshape
a = array([[.5,.5,0],[0,.5,.5],[.5,0,.5]])
b = array([1])
c = array([[0],[0],[0]])
result = get_spacegroup(a, b, c, True, 1e-10)
print(result.sg_op.flags)
for i in range(7):
    print(result.sg_op[:,:,i])
#result = make_primitive(a, b, c, True, 1e-10)

