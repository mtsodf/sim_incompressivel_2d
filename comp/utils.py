def deriv_num(f, x, dx=None):
    if dx is None:
        dx = 1e-6*x
        if dx == 0:
            print "Problemas com a derivada numerica"
            dx = 0.1

    return (f(x+dx)-f(x))/dx



