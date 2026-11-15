import sys

from branch_and_bound.BABbase import ler_instancia_tsplib, calcular_custo_rota, branch_and_bound_tsplib
from heuristicas_aproximacoes.christofides import aproximacao_christofides
#from trivial import construir_solucao_trivial


caminho = sys.argv[1]
nome, coords, matriz = ler_instancia_tsplib(caminho)
n = len(coords)

print(f"Instancia: {nome} | Vertices: {n}\n")

heuristicas = {
    #"trivial": construir_solucao_trivial,
    "christofides": aproximacao_christofides
}

for nomeH, heur in heuristicas.items():
    rota = heur(matriz, n)
    custoInicial = calcular_custo_rota(rota, matriz)
    melhorRota, melhorCusto = branch_and_bound_tsplib(matriz, rota)
    print(f"{nomeH}: rota={custoInicial} | bnb={melhorCusto}")