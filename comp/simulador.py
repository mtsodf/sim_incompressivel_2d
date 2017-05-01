import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from math import sqrt, pow, pi, log



def trans(k1,k2,d1,d2,area):
    denominador = d1/k1
    denominador += d2/k2

    return (d1+d2) * area / denominador

class Fluido(object):

    def __init__(self, visref, compf, rhoref, pref, vispr):
        self.visref = visref
        self.rhoref = rhoref
        self.compf = compf
        self.pref = pref
        self.vispr = vispr

    def viscosidade(self, pres):
        return self.visref +  self.vispr*(pres - self.pref)

    def densidade(self, pres):
        return self.rhoref*(1 + self.compf*(pres - self.pref))

    #Derivada da densidade pela pressao
    def d_densidade(self, pres=None):
        return self.rhoref*self.compf

    def mobilidade_massica(self, pres):
        return self.densidade(pres)/self.viscosidade(pres)

    #Derivada da mobilidade massica
    def d_mobilidade_massica(self, pres):
        des = self.densidade(pres)
        visco = self.viscosidade(pres)

        r = des*self.vispr - visco*self.d_densidade(pres)

        r /= visco*visco

        return r

    def __str__(self):
        s =  "Dados do Fluido\n"
        s += "\t Viscosidade Ref = %f cp\n" % self.visref
        s += "\t Densidade   Ref = %f kg/m3\n" % self.rhoref
        s += "\t Compressibilidade = %f cm2/kgf\n" % self.compf


        return s


class Poco(object):

    def __init__(self, nx, ny, nz, radius, pf):
        self.nx= nx
        self.ny = ny
        self.nz = nz
        self.pf = pf
        self.radius = radius

    def set_wi(self, kx, ky, dx, dy, h):
        kfrac = ky/kx
        req = 0.28*sqrt(sqrt(kfrac)*dx*dx+sqrt(1/kfrac)*dy*dy)
        req /= (pow(kfrac, 0.25) + pow(1/kfrac, 0.25))

        self.wi = 2*pi*sqrt(kx*ky)*h/log(req/self.radius)

class Modelo(object):
    """docstring for Modelo"""
    def __init__(self, Leitor):
        super(Modelo, self).__init__()

        #Lendo Caracteristicas da Malha
        self.nx, self.ny, self.nz = Leitor.read_integers("DIM")
        self.numcels = self.nx*self.ny*self.nz

        self.kx = Leitor.read_float("KX", self.numcels)
        self.ky = Leitor.read_float("KY", self.numcels)
        self.kz = Leitor.read_float("KZ", self.numcels)
        self.dx = Leitor.read_float("DX", self.nx)
        self.dy = Leitor.read_float("DY", self.ny)
        self.dz = Leitor.read_float("DZ", self.nz)
        self.poro = Leitor.read_float("POR", self.numcels)
        self.pref = Leitor.read_float("PREF")[0]
        self.pres = []
        self.pres.append(Leitor.read_float("PRES", self.numcels))


        #Lendo Tempo
        self.num_tsteps = Leitor.read_integers("TIMESTEPS")[0]
        self.delta_t = Leitor.read_float("DELTAT")[0]


        #Lendo Caracteristicas do Fluido
        visref = Leitor.read_float("VISREF")[0]
        rhoref = Leitor.read_float("RHOREF")[0]
        compf = Leitor.read_float("COMPF")[0]
        vispr = Leitor.read_float("VISPR")[0]

        self.fluido = Fluido(visref, compf, rhoref, self.pref, vispr)

        self.compr = Leitor.read_float("COMPR")[0]

        self.pocos = Leitor.read_pocos()


        self.calc_wis()


    def calc_wis(self):
        for poco in self.pocos:
            i,j,k = poco.nx, poco.ny, poco.nz
            ind = self.twod_ind(i,j,k)
            kx = self.kx[ind]
            ky = self.ky[ind]
            dx = self.dx[i]
            dy = self.dy[j]
            h = self.dz[k]

            poco.set_wi(kx,ky,dx,dy,h)


    def print_modelo(self):
        print "Dimensoes: %5d, %5d %5d" % (self.nx, self.ny, self.nz)
        print self.fluido


    def twod_ind(self, i, j, k):
        return i+j*self.nx+k*(self.nx*self.ny)


    def ind_2d(self, ind):
        k = ind / (self.nx * self.ny)
        ind = ind % (self.nx * self.ny)
        j = ind/self.nx
        i = ind%self.nx
        return i, j, k


    def verificar_compatibilidade(self):
        massa = 0

        for poco in self.pocos:
            massa += poco[-1]


        return massa == 0.0


    def porosidade(self, i, j, k, p):
        poro = self.poro[self.twod_ind(i, j, k)]

        return poro*(1+self.compr*(p-self.pref))

    def d_porosidade(self,i,j,k,pres=None):
        poro = self.poro[self.twod_ind(i, j, k)]
        return poro * self.compr


    def termo_acumulacao(self, i, j, k):

        pass


    def volume(self, i, j, k):
        return self.dx[i]*self.dy[j]*self.dz[k]

    def calc_jacob(self, pres1, pres0, dt):
        A = np.zeros((self.numcels, self.numcels))

        for ind in xrange(self.numcels):

            #Derivada do termo de acumulacao

            i, j, k = self.ind_2d(ind)
            vol = self.volume(i,j,k)

            p1 = pres1[ind]
            p0 = pres0[ind]

            porosidade = self.porosidade(i,j,k,p1)

            

            d_acum = (vol*dt)* (self.d_porosidade(i,j,k) * self.fluido.densidade(p1) + porosidade*self.fluido.d_densidade(p1))

            A[ind,ind] += d_acum



            #TODO colocar termo gravitacional
            if(i > 0):
                indviz = self.twod_ind(i-1,j,k)
                pviz = pres1[indviz]

                transmissibilidade = self.calcTrans(i, j, k, p1, pviz, "I-")

                trans += transmissibilidade*(p1 - pviz)
                #Derivada em relacao a p(n+1,i)
                A[ind, ind] += transmissibilidade*(self.fluido.d_mobilidade_massica(p1)*(p1-pviz) 
                                                    + self.fluido.mobilidade_massica(p1))

                #TODO verificar se a mobilidade eh calculada na celula
                A[ind, indviz] += transmissibilidade*(-self.fluido.mobilidade_massica(p1))


            if(i < self.nx -1):
                indviz = self.twod_ind(i+1,j,k)
                pviz = pres1[indviz]

                transmissibilidade = self.calcTrans(i, j, k, p1, pviz, "I+")

                #Derivada em relacao a p(n+1,i)
                A[ind, ind] += transmissibilidade*(self.fluido.d_mobilidade_massica(p1)*(p1-pviz) 
                                                    + self.fluido.mobilidade_massica(p1))

                #TODO verificar se a mobilidade eh calculada na celula
                A[ind, indviz] += transmissibilidade*(-self.fluido.mobilidade_massica(p1))

            if(j > 0):
                indviz = self.twod_ind(i,j-1,k)
                pviz = pres1[indviz]

                transmissibilidade = self.calcTrans(i, j, k, p1, pviz, "J-")

                #Derivada em relacao a p(n+1,i)
                A[ind, ind] += transmissibilidade*(self.fluido.d_mobilidade_massica(p1)*(p1-pviz) 
                                                    + self.fluido.mobilidade_massica(p1))

                #TODO verificar se a mobilidade eh calculada na celula
                A[ind, indviz] += transmissibilidade*(-self.fluido.mobilidade_massica(p1))

            if(j<self.ny-1):
                indviz = self.twod_ind(i,j+1,k)
                pviz = pres1[indviz]

                transmissibilidade = self.calcTrans(i, j, k, p1, pviz, "J+")

                #Derivada em relacao a p(n+1,i)
                A[ind, ind] += transmissibilidade*(self.fluido.d_mobilidade_massica(p1)*(p1-pviz) 
                                                    + self.fluido.mobilidade_massica(p1))

                #TODO verificar se a mobilidade eh calculada na celula
                A[ind, indviz] += transmissibilidade*(-self.fluido.mobilidade_massica(p1))

            if(k > 0):
                indviz = self.twod_ind(i,j,k-1)
                pviz = pres1[indviz]

                transmissibilidade = self.calcTrans(i, j, k, p1, pviz, "K-")

                #Derivada em relacao a p(n+1,i)
                A[ind, ind] += transmissibilidade*(self.fluido.d_mobilidade_massica(p1)*(p1-pviz) 
                                                    + self.fluido.mobilidade_massica(p1))

                #TODO verificar se a mobilidade eh calculada na celula
                A[ind, indviz] += transmissibilidade*(-self.fluido.mobilidade_massica(p1))
            
            if(k<self.nz-1):
                indviz = self.twod_ind(i,j,k+1)
                pviz = pres1[indviz]

                transmissibilidade = self.calcTrans(i, j, k, p1, pviz, "K+")

                #Derivada em relacao a p(n+1,i)
                A[ind, ind] += transmissibilidade*(self.fluido.d_mobilidade_massica(p1)*(p1-pviz) 
                                                    + self.fluido.mobilidade_massica(p1))

                #TODO verificar se a mobilidade eh calculada na celula
                A[ind, indviz] += transmissibilidade*(-self.fluido.mobilidade_massica(p1))

        print A

    def calc_residuo(self, pres1, pres0, dt):

        r = []

        for ind in xrange(self.numcels):
            i, j, k = self.ind_2d(ind)
            

            #Calculo do termo de acumulacao
            vol = self.volume(i, j, k)

            p1 = pres1[ind]
            p0 = pres0[ind]

            por1 = self.porosidade(i, j, k, p1)
            por0 = self.porosidade(i, j, k, p0)

            rho1 = self.fluido.densidade(p1)
            rho0 = self.fluido.densidade(p0)

            acum = vol * (rho1*por1 - rho0*por0) / dt

            #Transmissibilidades
            trans = 0


            #TODO colocar o termo gravitacional
            if(i > 0):
                pviz = pres1[self.twod_ind(i-1,j,k)]
                transmissibilidade = self.calcTrans(i, j, k, p1, pviz, "I-")
                trans += transmissibilidade*(p1 - pres1[self.twod_ind(i-1,j,k)])
            if(i < self.nx -1):
                pviz = pres1[self.twod_ind(i+1,j,k)]
                transmissibilidade = self.calcTrans(i, j, k, p1, pviz, "I+")
                trans += transmissibilidade*(p1 - pres1[self.twod_ind(i+1,j,k)])
            if(j > 0):
                pviz = pres1[self.twod_ind(i,j-1,k)]
                transmissibilidade = self.calcTrans(i, j, k, p1, pviz, "J-")
                trans += transmissibilidade*(p1 - pres1[self.twod_ind(i,j-1,k)])
            if(j<self.ny-1):
                pviz = pres1[self.twod_ind(i,j+1,k)]
                transmissibilidade = self.calcTrans(i, j, k, p1, pviz, "J+")
                trans += transmissibilidade*(p1 - pres1[self.twod_ind(i,j+1,k)])
            if(k > 0):
                pviz = pres1[self.twod_ind(i,j,k-1)]
                transmissibilidade = self.calcTrans(i, j, k, p1, pviz, "K-")
                trans += transmissibilidade*(p1 - pres1[self.twod_ind(i,j,k-1)])
            if(k<self.nz-1):
                pviz = pres1[self.twod_ind(i,j,k+1)]
                transmissibilidade = self.calcTrans(i, j, k, p1, pviz, "K+")
                trans += transmissibilidade*(p1 - pres1[self.twod_ind(i,j,k+1)])


            
            r.append(trans + acum)

        #TODO adicionar o termo dos pocos
        for poco in self.pocos:
            i,j,k=poco.nx, poco.ny, poco.nz

            ind = self.twod_ind(i,j,k)

            r[ind] -= self.fluido.mobilidade_massica(p1)*poco.wi*(p1-poco.pf)



        return r

    def simular(self):

       print self.calc_residuo(self.pres[0], self.pres[0], 10)

       self.calc_jacob(self.pres[0], self.pres[0], 10)

    def vizinhos(self, i, j, k):
        v = []
        if(i > 0):
            v.append((i-1,j,k))
        if(i < self.nx -1):
            v.append((i+1,j,k))
        if(j > 0):
            v.append((i,j-1,k))
        if(j<self.ny-1):
            v.append((i,j+1,k))
        if(k > 0):
            v.append((i,j,k-1))
        if(k<self.nz-1):
            v.append((i,j,k+1))

        return v


    def calcTrans(self, i, j, k, p1, pviz, dir):

        k1=None; k2=None; d1=None; d2=None; area=None

        if dir == "I+":
            k1 = self.kx[self.twod_ind(i,j,k)]
            k2 = self.kx[self.twod_ind(i+1,j,k)]
            d1 = self.dx[i]
            d2 = self.dx[i+1]
            area = self.dy[j]*self.dz[k]

        elif dir == "I-":
            k1 = self.kx[self.twod_ind(i,j,k)]
            k2 = self.kx[self.twod_ind(i-1,j,k)]
            d1 = self.dx[i]
            d2 = self.dx[i-1]
            area = self.dy[j]*self.dz[k]

        elif dir == "J-":
            k1 = self.ky[self.twod_ind(i,j,k)]
            k2 = self.ky[self.twod_ind(i,j-1,k)]
            d1 = self.dy[j]
            d2 = self.dy[j-1]
            area = self.dx[i]*self.dz[k]

        elif dir == "J+":
            k1 = self.ky[self.twod_ind(i,j,k)]
            k2 = self.ky[self.twod_ind(i,j+1,k)]
            d1 = self.dy[j]
            d2 = self.dy[j+1]
            area = self.dx[i]*self.dz[k]

        elif dir == "K-":
            k1 = self.kz[self.twod_ind(i,j,k)]
            k2 = self.kz[self.twod_ind(i,j,k-1)]
            d1 = self.dz[k]
            d2 = self.dz[k-1]
            area = self.dx[i]*self.dy[j]

        elif dir == "K+":
            k1 = self.kz[self.twod_ind(i,j,k)]
            k2 = self.kz[self.twod_ind(i,j-1,k)]
            d1 = self.dz[k]
            d2 = self.dz[k-1]
            area = self.dx[i]*self.dy[j]



        #TODO verificar onde eh calculada a mobilidade.
        mob1 = self.fluido.mobilidade_massica(p1)
        mobviz = self.fluido.mobilidade_massica(pviz)

        #TODO verificar se eh media harmonica
        mob = 2/(1/mob1 + 1/mobviz)
 
        return trans(k1,k2,d1,d2,area)* mob



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
            i, j, k, radius, pf = x
            i = int(i)
            j = int(j)
            k = int(k)
            radius = float(radius)
            pf = float(pf)

            pocos.append(Poco(i, j, k, radius, pf))

        return pocos


    def read_integers(self, keyword, qtd_expected=None):
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

        if qtd_expected and len(values)!=qtd_expected:
            raise Exception("Quantidade esperada para %s eh %d. Quantidade lida %d" % (keyword, qtd_expected, len(values)))

        return values

    def read_float(self, keyword, qtd_expected=None):
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

        if qtd_expected and len(values)!=qtd_expected:
            raise Exception("Quantidade esperada para %s eh %d. Quantidade lida %d" % (keyword, qtd_expected, len(values)))

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

        raise Exception("Nao foi possivel encontrar a palavra %s" % keyword)


Leitor = LeitorDeArquivo(sys.argv[1])

modelo = Modelo(Leitor)


modelo.print_modelo()

modelo.simular()

