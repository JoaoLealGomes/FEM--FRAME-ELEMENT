from geo import geometria
from load import carregamento
from resultados import results
import numpy as np


# Coordenadas dos nós [x,y]
a = 4.8
b = 3.2
coordenadas = np.array([[0, 0],
                        [a, b],
                        [a, 0],
                        [2*a, b],
                        [2*a + a/4, b],
                        [2*a, b/2],
                        [2*a, 0]])

# Vetor de materiais [Elas., Area, Inercia]
materiais = np.array([[2e8, 5e-03, 1.28e-04]])

# Matriz de elementos
elementos = np.array([[0, 1, 0],
                        [1, 2, 0],
                        [1, 3, 0],
                        [3, 4, 0],
                        [3, 5, 0],
                        [5, 6, 0]])

# Vinculos [Nó, Horizontal, Vertical, Rotação]
apoio = np.array([[0, 1, 1, 0],
                  [1, 0, 0, 0],
                  [2, 1, 1, 1],
                  [3, 0, 0, 0],
                  [4, 0, 0, 0],
                  [5, 0, 0, 0],
                  [6, 1, 1, 0]])


portico = geometria(coordenadas, materiais, elementos, apoio)
#portico.plot()
cargas = carregamento(coordenadas)

P = 32
q = 12
H = 30

cargas.carga_pontual(1, 'y', -P)
cargas.carga_distribuida(1, 2, 'x', -q)
cargas.carga_distribuida(1, 3, 'y', 0, -q)
cargas.carga_pontual(4, 'y', -P)
cargas.carga_pontual(5, 'x', H)

cargas.repeticao()


solver = results(portico, cargas.cargas, cargas.nodais_eq)
solver.solve()
solver.reacao_apoio()
solver.esforcos()

print(f'Deslocamentos: \n {solver.d}')
print(f'Reações de Apoio: \n {solver.reacao}')
print(f'Momento Fletor: \n {solver.M_Fletor}')
print(f'Esforço Cortante:\n {solver.Q_cortante}')
print(f'Esforço Normal: \n {solver.N}')

print(f'Momentos Nodais: \n {solver.M}')
print(f'Esforço Cortante Nodal:\n {solver.Q}')

