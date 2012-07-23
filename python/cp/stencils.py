"""
Spatial differentiation operator stencils for Closest Point Method
codes
"""
from numpy import array as a
from numpy import isscalar

def Laplacian_2nd(dim):
    """
    Second-order approximation to the Laplacian (the 5 pt stencil in
    2D)
    """
    if (dim == 2):
        def f(dx):
            """Weights for second-order approx to Laplacian"""
            if (isscalar(dx)):
                dx2 = dx*dx
                return [-4/dx2, 1/dx2, 1/dx2, 1/dx2, 1/dx2]
            else:
                dx2 = dx[0]*dx[0]
                dy2 = dx[1]*dx[1]
                return [-2/dx2-2/dy2, 1/dx2, 1/dx2, 1/dy2, 1/dy2]
        # a function of dx for the weights
        DiffWeightsFcn = f
        # the points in the stencil
        DiffStencil = [ a([ 0,  0]), \
                        a([ 1,  0]), \
                        a([-1,  0]), \
                        a([ 0,  1]), \
                        a([ 0, -1]) ]
        # stencil is typically a hypercross: longest arm of it is used
        # in calculating e.g., bandwidths
        DiffLongestArm = 1
    elif (dim == 3):
        def f(dx):
            if (isscalar(dx)):
                dx2 = dx*dx
                return [-6/dx2, 1/dx2, 1/dx2, 1/dx2, 1/dx2, 1/dx2, 1/dx2]
            else:
                dx2 = dx[0]*dx[0]
                dy2 = dx[1]*dx[1]
                dz2 = dx[2]*dx[2]
                return [-2/dx2-2/dy2-2/dz2, 1/dx2, 1/dx2, 1/dy2, 1/dy2, 1/dz2, 1/dz2]
        DiffWeightsFcn = f
        DiffStencil = [ a([ 0,  0,  0]), \
                        a([ 1,  0,  0]), \
                        a([-1,  0,  0]), \
                        a([ 0,  1,  0]), \
                        a([ 0, -1,  0]), \
                        a([ 0,  0,  1]), \
                        a([ 0,  0, -1]) ]
        DiffLongestArm = 1
    else:
        raise NameError('dim ' + str(dim) + ' not implemented')
    return (DiffWeightsFcn, DiffStencil, DiffLongestArm)


def Laplacian_4th(dim):
    """
    Forth-order approximation to the Laplacian
    """
    if (dim == 2):
        def f(dx):
            d2 = dx*dx
            #return [ -5.0/d2, \
            #        (-1.0/12.0)/d2, (4.0/3.0)/d2, (4.0/3.0)/d2, (-1.0/12.0)/d2, \
            #        (-1.0/12.0)/d2, (4.0/3.0)/d2, (4.0/3.0)/d2, (-1.0/12.0)/d2 ]
            print "** WARNING: hardcoded f96 stuff **"
            # TODO: this won't work if dx is integer
            #fl = type(dx)
            from numpy import float96
            fl = float96
            # f1o12: float of 1/12
            f1o12 = fl(1)/fl(12)
            f4o3 = fl(4)/fl(3)

            return [ -fl(5)/d2, \
                     -f1o12/d2, f4o3/d2, f4o3/d2, -f1o12/d2, \
                     -f1o12/d2, f4o3/d2, f4o3/d2, -f1o12/d2 ]
        DiffWeightsFcn = f
        DiffStencil = [ a([ 0,   0]), \
                        a([-2,   0]), \
                        a([-1,   0]), \
                        a([ 1,   0]), \
                        a([ 2,   0]), \
                        a([ 0,  -2]), \
                        a([ 0,  -1]), \
                        a([ 0,   1]), \
                        a([ 0,   2]) ]
        DiffLongestArm = 2
    elif (dim == 3):
        print 'WARNING: not tested yet since copy-paste from matlab code'
        def f(dx):
            d2 = dx*dx
            return [ (-15.0/2.0)/d2, \
		     (-1.0/12.0)/d2, (4.0/3.0)/d2, (4.0/3.0)/d2, (-1.0/12.0)/d2, \
		     (-1.0/12.0)/d2, (4.0/3.0)/d2, (4.0/3.0)/d2, (-1.0/12.0)/d2, \
		     (-1.0/12.0)/d2, (4.0/3.0)/d2, (4.0/3.0)/d2, (-1.0/12.0)/d2 ]
        DiffWeightsFcn = f
        DiffStencil = [ a([ 0,  0,  0]), \
		        a([-2,  0,  0]), \
		        a([-1,  0,  0]), \
		        a([ 1,  0,  0]), \
		        a([ 2,  0,  0]), \
		        a([ 0, -2,  0]), \
		        a([ 0, -1,  0]), \
		        a([ 0,  1,  0]), \
		        a([ 0,  2,  0]), \
		        a([ 0,  0, -2]), \
		        a([ 0,  0, -1]), \
		        a([ 0,  0,  1]), \
		        a([ 0,  0,  2]) ]
        DiffLongestArm = 2
    else:
        raise NameError('dim ' + str(dim) + ' not implemented')
    return (DiffWeightsFcn, DiffStencil, DiffLongestArm)



def Biharmonic_2nd(dim):
    """
    WARNING: PROBABLY YOU DON'T WANT THIS FOR THE CLOSEST POINT
    METHOD.  It will be inconsistent.  Instead square the Laplacian
    matrix like in [Macdonald&Ruuth 2009].
    """
    if dim==2:
        def f(dx):
            d4 = dx*dx*dx*dx
            return [ -20/d4,  -1/d4,  8/d4,  8/d4,  -1/d4,  -1/d4,  8/d4,  8/d4,  -1/d4, \
                     -2/d4,  -2/d4,  -2/d4,  -2/d4 ]
        DiffWeightsFcn = f
        DiffStencil = [ a([ 0,  0]), \
		        a([-2,  0]), \
		        a([-1,  0]), \
		        a([ 1,  0]), \
		        a([ 2,  0]), \
		        a([ 0, -2]), \
		        a([ 0, -1]), \
		        a([ 0,  1]), \
		        a([ 0,  2]), \
                        a([ 1,  1]), \
                        a([ 1, -1]), \
                        a([-1,  1]), \
                        a([-1, -1]) ]
        DiffLongestArm = 2
    else:
        raise NameError('dim ' + str(dim) + ' not implemented')
    return (DiffWeightsFcn, DiffStencil, DiffLongestArm)