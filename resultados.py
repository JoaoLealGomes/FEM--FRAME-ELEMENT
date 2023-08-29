import math
import numpy as np
from solver import frame_solver as fs

class results(fs):
    def __init__(self, geometria, cargas, nodais_eq):
        self.geometria = geometria
        self.cargas = cargas
        self.qnodais = nodais_eq
        super().__init__(self.geometria, self.cargas)

    def reacao_apoio(self):
        self.fp = np.matmul(self.Keg_ord[self.dim:, :self.dim], self.dl)
        self.reacao = self.fp - self.F_ord[self.dim:]


    def __momento(self, L, x):
        return np.array(([0, 1.5*x, L*(-2 +6*x)/8, 0, -1.5*x, L*(2 + 6*x)/8]))


    def __cortante(self, L):
        return np.array(([0, 1.5, 0.75*L, 0, -1.5, 0.75*L]))
    

    def __resultado_final(self):
        self.M_Fletor = np.copy(self.M)
        self.Q_cortante = np.copy(self.Q)

        for i in range(0, self.qnodais.shape[0], 2):
            for elemento in self.elementos:
                if elemento[0] == self.qnodais[i][0] and elemento[1] == self.qnodais[i + 1][0]:
                    index = (np.where(np.all(self.elementos == elemento, axis=1) == True)[0][0])*2
                    self.M_Fletor[index] = self.M_Fletor[index] + self.qnodais[i][2]
                    self.M_Fletor[index + 1] = self.M_Fletor[index + 1] - self.qnodais[i + 1][2]
        
        for elemento in self.elementos:
            no_inic = elemento[0]
            no_final = elemento[1]
            xni = self.coordenadas[no_inic,0]
            yni = self.coordenadas[no_inic,1]
            xnf = self.coordenadas[no_final,0]
            ynf = self.coordenadas[no_final,1]
            sinal = -1 if xnf - xni == 0 and ynf - yni != 0 else 1
            index = (np.where(np.all(self.elementos == elemento, axis=1) == True)[0][0])*2
            self.M_Fletor[index] = sinal*self.M_Fletor[index]
            self.M_Fletor[index + 1] = sinal*self.M_Fletor[index + 1]


        for i in range(0, self.qnodais.shape[0], 2):
            for elemento in self.elementos:
                if elemento[0] == self.qnodais[i][0] and elemento[1] == self.qnodais[i + 1][0]:
                    index = (np.where(np.all(self.elementos == elemento, axis=1) == True)[0][0])*2
                    self.Q_cortante[index] = self.Q_cortante[index] - self.qnodais[i][1]
                    self.Q_cortante[index + 1] = self.Q_cortante[index + 1] + self.qnodais[i + 1][1]
                

    def esforcos(self):
        self.N = np.zeros((self.N_elem, 1)) 
        self.M = np.array([])
        self.Q = np.array([])
        for i in range(self.N_elem):
            
            no_inic = self.elementos[i,0]
            no_final = self.elementos[i,1]
            
            xni = self.coordenadas[no_inic,0]
            yni = self.coordenadas[no_inic,1]
            xnf = self.coordenadas[no_final,0]
            ynf = self.coordenadas[no_final,1]
            
            L = ((xnf-xni)**2 + (ynf-yni)**2)**(1/2)
            E = self.materiais[self.elementos[i,2],0]
            I = self.materiais[self.elementos[i,2],2]
            A = self.materiais[self.elementos[i,2],1]
            
            dx, dy = (xnf - xni),  (ynf - yni)

            if dx == 0:
                if dy < 0:
                    theta = -math.pi/2.0
                else:
                    theta = math.pi/2.0
            else:
                theta = math.atan(dy/dx)       
            
            c, s = math.cos(theta), math.sin(theta)

            T = np.array(([[c, s, 0, 0, 0, 0],
                        [-s, c, 0, 0, 0, 0],
                        [0, 0, 1, 0, 0, 0],
                        [0, 0, 0, c, s, 0],
                        [0, 0, 0, -s, c, 0],
                        [0, 0, 0, 0, 0, 1]]))
            posicao = [3 * no_inic, 3 * no_inic + 1, 3 * no_inic + 2, 3 * no_final, 3 * no_final + 1, 3 * no_final + 2]

            U = np.array([[ self.d[posicao[0]]],
                            [self.d[posicao[1]]],
                            [self.d[posicao[2]]],
                            [self.d[posicao[3]]],
                            [self.d[posicao[4]]],
                            [self.d[posicao[5]]]])
            
            ue = np.dot(T, U)

            self.N[i, 0] = E*A/L*(ue[3,0] - ue[0,0]) 
            
            self.M = np.append(self.M, self.__momento(L, -1)*4*E*I/(L**2) @ ue)
            self.M = np.append(self.M, self.__momento(L, 1)*4*E*I/(L**2) @ ue)
            self.Q = np.append(self.Q, self.__cortante(L)*8*E*I/(L**3) @ ue)
            self.Q = np.append(self.Q, self.__cortante(L)*8*E*I/(L**3) @ ue)

        self.__resultado_final()
