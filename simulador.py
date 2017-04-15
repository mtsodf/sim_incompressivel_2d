import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid


def plot(vetor, nx, ny):
    F = plt.figure(1, (5.5, 3.5))
    my_cmap = plt.cm.get_cmap('seismic')
    grid = ImageGrid(F, 111,  # similar to subplot(111)
                     nrows_ncols=(1, 1),
                     axes_pad=0.1,
                     add_all=True,
                     label_mode="L",
                     cbar_location="right",
                     cbar_mode="single",
                     )



    g = np.zeros([ny, nx])

    ind = 0
    for j in range(ny):
        for i in range(nx):
            g[j, i] = vetor[ind]
            ind += 1

    im1 = g

    vmin, vmax = g.min(), g.max()
    ax = grid[0]
    im = ax.imshow(im1, origin="lower", vmin=vmin, vmax=vmax, cmap=my_cmap, interpolation="nearest")
    grid.cbar_axes[0].colorbar(im)

    plt.draw()
    plt.show()



def trans(k1,k2,d1,d2,visco,a1, esp):
    return 2*k1*k2*a1*esp/(visco*(d1*k2 + d2*k1))

class Modelo(object):
    """docstring for Modelo"""
    def __init__(self, Leitor):
        super(Modelo, self).__init__()
        self.nx, self.ny = Leitor.read_integers("DIM")
        self.kx = Leitor.read_float("KX")
        self.ky = Leitor.read_float("KY")
        self.dx = Leitor.read_float("DX")
        self.dy = Leitor.read_float("DY")
        self.visco = Leitor.read_float("VISCO")[0]
        self.esp = Leitor.read_float("ESPESSURA")[0]
        self.pocos = Leitor.read_pocos()


        self.numcels = self.nx*self.ny


    def print_modelo(self):
        print "Dimensoes: %5d, %5d" % (self.nx, self.ny)
        print "Viscosidade: %10.2f" % self.visco



    def twod_ind(self, i , j):
        return i+j*self.nx

    def ind_2d(self, ind):
        i = ind % j
        j = ind/self.nx
        return i, j


    def verificar_compatibilidade(self):
        massa = 0

        for poco in self.pocos:
            massa += poco[-1]


        return massa == 0.0



    def calcTrans(self, i, j, dir):
        visco = self.visco
        esp = self.esp
        if dir == "X+":
            k1 = self.kx[self.twod_ind(i,j)]
            k2 = self.kx[self.twod_ind(i+1,j)]
            d1 = self.dx[i]
            d2 = self.dx[i+1]
            a1 = self.dy[j]
            return trans(k1,k2,d1,d2,visco,a1,esp)

        if dir == "X-":
            k1 = self.kx[self.twod_ind(i,j)]
            k2 = self.kx[self.twod_ind(i-1,j)]
            d1 = self.dx[i]
            d2 = self.dx[i-1]
            a1 = self.dy[j]
            return trans(k1,k2,d1,d2,visco,a1,esp)

        if dir == "Y+":
            k1 = self.ky[self.twod_ind(i,j)]
            k2 = self.ky[self.twod_ind(i,j+1)]
            d1 = self.dy[j]
            d2 = self.dy[j+1]
            a1 = self.dx[i]
            return trans(k1,k2,d1,d2,visco,a1,esp)

        if dir == "Y-":
            k1 = self.ky[self.twod_ind(i,j)]
            k2 = self.ky[self.twod_ind(i,j-1)]
            d1 = self.dy[j]
            d2 = self.dy[j-1]
            a1 = self.dx[i]
            return trans(k1,k2,d1,d2,visco,a1,esp)

    def build_matrix(self):
        self.matriz=np.zeros([self.nx*self.ny, self.nx*self.ny])

        for j in range(self.ny):
            for i in range(self.nx):
                ind  = self.twod_ind(i,j)
                if i > 0:
                    t = self.calcTrans(i, j, "X-")
                    self.matriz[ind, ind] += t
                    ind2 = self.twod_ind(i-1, j)
                    self.matriz[ind, ind2] -= t

                if i < self.nx-1:
                    t = self.calcTrans(i, j, "X+")
                    self.matriz[ind, ind] += t
                    ind2 = self.twod_ind(i+1, j)
                    self.matriz[ind, ind2] -= t

                if j > 0:
                    t = self.calcTrans(i, j, "Y-")
                    self.matriz[ind, ind] += t
                    ind2 = self.twod_ind(i, j-1)
                    self.matriz[ind, ind2] -= t

                if j < self.nx-1:
                    t = self.calcTrans(i, j, "Y+")
                    self.matriz[ind, ind] += t
                    ind2 = self.twod_ind(i, j+1)
                    self.matriz[ind, ind2] -= t      
       

    def tem_poco(self, i, j):
        for poco in self.pocos:
            if(i == poco[0] and j == poco[1]):
                return True

        return False


    def set_zero_matriz(self):
        for j in range(self.ny):
            for i in range(self.nx):
                if not self.tem_poco(i, j):
                    ind = self.twod_ind(i, j)
                    for k in range(self.numcels):
                        self.matriz[ind, k] = 0

                    self.matriz[ind, ind] = 1

                    return

        print "Todas as celulas possuem pocos..."




    def calc_ld(self):
        ld = np.zeros(self.nx*self.ny)

        for i,j,q in self.pocos:
            ind = self.twod_ind(i, j)
            ld[ind] = q

        self.ld = ld



    def solve(self):
        self.build_matrix()
        self.set_zero_matriz()
        self.calc_ld()
        x = np.linalg.solve(self.matriz, self.ld)

        x = x - np.mean(x)
        self.x = x



class LeitorDeArquivo(object):
    """docstring for LeitorDeArquivo"""
    def __init__(self, filename):
        super(LeitorDeArquivo, self).__init__()
        self.filename = filename

        #self.nx, self.ny = self.read_integers("DIM")

        #print self.nx, self.ny

    def read_pocos(self):
        l = self.read_keyword("POCOS", True)
        pocos = []
        for x in l:
            i, j, vazao = x
            i = int(i)
            j = int(j)
            vazao = float(vazao)

            pocos.append([i, j, vazao])
        return pocos


    def read_integers(self, keyword):
        l = self.read_keyword(keyword)
        values = []

        for x in l:
            if "*" in x:
                qtd, value = x.split("*")
                qtd = int(qtd)
                value = int(value)
                values.extend(qtd*[value])
            else:
                values.append(int(x))

        return values

    def read_float(self, keyword):
        l = self.read_keyword(keyword)
        values = []

        for x in l:
            if "*" in x:
                qtd, value = x.split("*")
                qtd = int(qtd)
                value = float(value)
                values.extend(qtd*[value])
            else:
                values.append(float(x))

        return values

    def read_keyword(self, keyword, multiple_lines=False):
        f = open(self.filename)

        for line in f:
            if keyword in line:
                line = next(f)
                if multiple_lines:
                    r = []
                    while len(line.strip()) > 0:
                        r.append(line.split())
                        line = next(f)
                    f.close()
                    return r
                f.close()
                return line.split()


Leitor = LeitorDeArquivo(sys.argv[1])

modelo = Modelo(Leitor)

print modelo.verificar_compatibilidade()

modelo.print_modelo()

modelo.build_matrix()

modelo.calc_ld()

modelo.set_zero_matriz()

modelo.solve()


plot(modelo.x, modelo.nx, modelo.ny)

