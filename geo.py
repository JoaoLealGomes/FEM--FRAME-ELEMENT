import numpy as np
import matplotlib.pyplot as plt

class geometria():
    def __init__(self, coord, mat, ele, apoio) -> None:
        self.coordenadas = coord
        self.materiais = mat
        self.elementos = ele
        self.apoio = apoio
        self.graus_lib = 3
        self.N_nos = self.coordenadas.shape[0]
        self.N_elem = self.elementos.shape[0]
        self.ordem()
        self.kll()

    def ordem(self):
        fixos = []
        livres = []
        self.N_cond = 0
        for i in range(self.N_nos):
            for j in range(1,4):
                if self.apoio[i, j] == 1:
                    self.N_cond += 1
                    fixos.append((i * self.graus_lib) + j - 1)
                else:
                    livres.append((i * self.graus_lib) + j - 1)

        self.ordem = np.append(livres, fixos)

    def kll(self):
        self.dim = self.graus_lib * self.N_nos - self.N_cond


    def plot(self):
        nx = np.transpose(self.coordenadas[:self.N_nos,:1])
        ny = np.transpose(self.coordenadas[:self.N_nos,1:])
        plt.scatter(nx, ny, color='red')  # plotagem dos nó
        for i in range(len(nx)):
            for j in range(len(nx[i])):
                plt.text(nx[i][j] - 0.45, ny[i][j] - 0.15, j)
        for i in range(self.N_elem):
            x = np.array([self.coordenadas[self.elementos[i,0],0],self.coordenadas[self.elementos[i,1],0]])
            y = np.array([self.coordenadas[self.elementos[i,0],1],self.coordenadas[self.elementos[i,1],1]])            
            plt.plot(x, y)       # plotagem da linha
        plt.title("Geometria")
        plt.show()



        pass