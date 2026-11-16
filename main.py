import sys

from leitor.leitor import ler_instancia_tsplib
from branch_and_bound.BABbase import calcular_custo_rota, branch_and_bound_tsplib
from heuristicas_aproximacoes.christofides import aproximacao_christofides
#from trivial import construir_solucao_trivial

from heuristicas_aproximacoes.twice_around_the_tree import aproximacao_twice_around_the_tree

from heuristicas_aproximacoes.convex_hull import heuristica_convex_hull

from heuristicas_aproximacoes.gulosa import heuristica_gulosa


caminho = sys.argv[1]
nome, coords, matriz = ler_instancia_tsplib(caminho)
n = len(coords)

print(f"Instancia: {nome} | Vertices: {n}\n")

heuristicas = {
    #"trivial": construir_solucao_trivial,
    #"christofides": aproximacao_christofides,
    #"twice_around_the_tree": aproximacao_twice_around_the_tree
    #"convex_hull": lambda m, n: heuristica_convex_hull(m, n, coords),
    "gulosa": heuristica_gulosa
}

for nomeH, heur in heuristicas.items():
    rota = heur(matriz, n)
    custoInicial = calcular_custo_rota(rota, matriz)
    melhorRota, melhorCusto = branch_and_bound_tsplib(matriz, rota)

    print(f"{nomeH}:")
    print(f" rota inicial : {rota} | custo = {custoInicial}")
    print(f" melhor rota : {melhorRota} | custo = {melhorCusto}")