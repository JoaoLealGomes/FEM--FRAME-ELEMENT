import math
import numpy as np
from scipy import linalg

class frame_solver():
    def __init__(self, geometria, cargas):
        self.cargas = cargas
        self.coordenadas = geometria.coordenadas    
        self.elementos = geometria.elementos
        self.materiais = geometria.materiais
        self.N_nos = geometria.N_nos
        self.N_elem = geometria.N_elem
        self.graus_lib = geometria.graus_lib
        self.ordem = geometria.ordem
        self.dim = geometria.dim
        

    def __K_eg(self):
        keg = np.zeros((3 * self.N_nos,3 * self.N_nos))
        
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
            
            k_ele = np.array([[E*A/L, 0, 0, -E*A/L, 0, 0],
                            [0, 12*E*I/(L**3), 6*E*I/(L**2), 0, -12*E*I/(L**3), 6*E*I/(L**2)],
                            [0, 6*E*I/(L**2), 4*E*I/L, 0, -6*E*I/(L**2), 2*E*I/L],
                            [-E*A/L, 0, 0, E*A/L, 0, 0],
                            [0, -12*E*I/(L**3), -6*E*I/(L**2), 0, 12*E*I/(L**3), -6*E*I/(L**2)],
                            [0, 6*E*I/(L**2), 2*E*I/L, 0, -6*E*I/(L**2), 4*E*I/L]])   
 
            k_ele = np.dot(np.dot(np.transpose(T), k_ele), T) # T^T * K_ele * T

            posi = [3 * no_inic, 3 * no_inic + 1, 3 * no_inic + 2, 3 * no_final, 3 * no_final + 1, 3 * no_final + 2]
            
            for j in range(6):
                for k in range(6):
                    keg[posi[j],posi[k]] += k_ele[j,k]
         
        self.Keg = keg
        # self.Keg[7][7] += 35000 # Mola Vertical
        # self.Keg[8][8] += 20000 # Mola Rotacional
    
    def __reordenacao(self):
        keg_ord = np.zeros((3*self.N_nos,3*self.N_nos))
        for i in range(3*self.N_nos):
            for j in range(3*self.N_nos):
                keg_ord[i,j] = self.Keg[self.ordem[i],self.ordem[j]]
        
        self.Keg_ord = keg_ord
    
    def __K_kll(self):
        self.Keg_kll = self.Keg_ord[:self.dim,:self.dim]
    

    def __vetor_forcas(self):
        # Vetor de Forças
        self.F = np.zeros(self.graus_lib * self.N_nos)
        for carga in self.cargas:
            no = carga[0]
            for j in range(self.graus_lib):
                self.F[int(self.graus_lib*no + j)] = carga[j + 1]

    def __F_ord(self):
        self.F_ord = np.zeros(3*self.N_nos)
        for i in range(3*self.N_nos):
            self.F_ord[i] = self.F[self.ordem[i]]
    
    def __FKll(self):
        self.F_kll = self.F_ord[:self.dim]

    
    def __deslocamentos(self):
        self.dl = linalg.solve(self.Keg_kll, self.F_kll)
        self.d = np.zeros((3*self.N_nos))
        for i in range(self.dim):
            self.d[self.ordem[i]] = self.dl[i]

    def solve(self):
        self.__K_eg()
        self.__reordenacao()
        self.__K_kll()
        self.__vetor_forcas()
        self.__F_ord()
        self.__FKll()
        self.__deslocamentos()


                