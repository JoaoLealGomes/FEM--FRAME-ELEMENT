import numpy as np


class carregamento():
    def __init__(self, coord): # [nó inical, nó final, valor inicial, valor final, axis]
        self.coordenadas = coord
        self.cargas = np.empty((0,4))
        self.nodais_eq = np.empty((0,3))
        
    def carga_pontual(self, no, eixo, magnitude):
        if eixo == 'x':
            direcao = np.array([[0, 1, 0, 0]])*magnitude
        elif eixo == 'y':
            direcao = np.array([[0, 0, 1, 0]])*magnitude
        elif eixo == 'm':
            direcao = np.array([[0, 0, 0, 1]])*magnitude
        else:
            raise NameError('O eixo de aplicação da carga deve ser eixo x, y ou m(momento)!')
        direcao[0][0] = no
        self.cargas = np.append(self.cargas, direcao, axis=0)
        
        
    def carga_distribuida(self, no_i, no_f, eixo, qi, qf=None):
        xni = self.coordenadas[no_i,0]
        yni = self.coordenadas[no_i,1]
        xnf = self.coordenadas[no_f,0]
        ynf = self.coordenadas[no_f,1]
              
        L = ((xnf-xni)**2 + (ynf-yni)**2)**(1/2)

        q = (qi + qf)/2 if qf != None else qi
        qx = qf - qi if qf != None else 0

        if eixo == 'y':
            carga = np.array([[no_i, 0, (L/2)*(q - qx/5), (q*L/6 - qx*L/60)*(L/2)]])
            self.cargas = np.append(self.cargas, carga, axis=0)

            neq = np.array([[no_i, (L/2)*(q - qx/5), (q*L/6 - qx*L/60)*(L/2)]])
            self.nodais_eq = np.append(self.nodais_eq, neq, axis=0)

            carga = np.array([[no_f, 0, (L/2)*(q + qx/5), (-q*L/6 - qx*L/60)*(L/2)]])
            self.cargas = np.append(self.cargas, carga, axis=0)

            neq = np.array([[no_f, (L/2)*(q + qx/5), (-q*L/6 - qx*L/60)*(L/2)]])
            self.nodais_eq = np.append(self.nodais_eq, neq, axis=0)
        elif eixo == 'x':
            carga = np.array([[no_i, (L/2)*(q - qx/5), 0, (q*L/6 - qx*L/60)*(L/2)]])
            self.cargas = np.append(self.cargas, carga, axis=0)

            neq = np.array([[no_i, (L/2)*(q - qx/5), (q*L/6 - qx*L/60)*(L/2)]])
            self.nodais_eq = np.append(self.nodais_eq, neq, axis=0)
            
            carga = np.array([[no_f, (L/2)*(q + qx/5), 0, (-q*L/6 - qx*L/60)*(L/2)]])
            self.cargas = np.append(self.cargas, carga, axis=0)

            neq = np.array([[no_f, (L/2)*(q + qx/5), (-q*L/6 - qx*L/60)*(L/2)]])
            self.nodais_eq = np.append(self.nodais_eq, neq, axis=0)
        else:
            raise NameError('O eixo de aplicação da carga deve ser eixo x ou y!')

    def repeticao(self):
        indices = []
        for i in range(len(self.cargas)):
            indices.append(self.cargas[i][0])

        indices = list(set(indices))
        indices = np.array(indices)
        
        cargas_finais = []
        for i in range(indices.shape[0]):
            soma = np.array([0., 0., 0., 0.])
            for j in range(len(self.cargas)):
                if self.cargas[j][0] == indices[i]:
                    soma += self.cargas[j]
            soma[0] = indices[i]
            cargas_finais.append(soma)
        cargas_finais = np.array(cargas_finais)
        self.cargas = cargas_finais

