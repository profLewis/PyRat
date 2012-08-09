import numpy
from mayavi import mlab

def test_triangular_mesh():
    """An example of a cone, ie a non-regular mesh defined by its
        triangles.
    """
    n = 8
    t = numpy.linspace(-numpy.pi, numpy.pi, n)
    z = numpy.exp(1j*t)
    x = z.real.copy()
    y = z.imag.copy()
    z = numpy.zeros_like(x)

    triangles = [(0, i, i+1) for i in range(1, n)]
    x = numpy.r_[0, x]
    y = numpy.r_[0, y]
    z = numpy.r_[1, z]
    t = numpy.r_[0, t]

    return mlab.triangular_mesh(x, y, z, triangles, scalars=t)

@mlab.show
def do():
  test_triangular_mesh()

do()
