import unittest
from utils import deriv_num
from math import sin, pi

class TestStringMethods(unittest.TestCase):

    def test_upper(self):

        f = lambda x : x*x

        d = deriv_num(f, 1, 0.000001)
        self.assertAlmostEqual(d, 2.0, places=4)

        f = lambda x: sin(x)

        d = deriv_num(f, pi/2)
        self.assertAlmostEqual(d, 0.0, places=5)

        d = deriv_num(f, 0.0, 1e-6)
        self.assertAlmostEqual(d, 1.0, places=5)        

        f = lambda x,y: x*x + y*y + x*y

        yconst = 2
        d= deriv_num(lambda x: f(x, yconst), 2)
        self.assertAlmostEqual(d, 6.0, places=5)   
        print d


if __name__ == '__main__':
    unittest.main()