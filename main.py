import sys

from analise.analise import Contador, EstatisticasBAB
from leitor.leitor import ler_instancia_tsplib
from branch_and_bound.BABbase import branch_and_bound_tsplib
from heuristicas_aproximacoes.christofides import aproximacao_christofides
from heuristicas_aproximacoes.trivial import construir_solucao_trivial
from heuristicas_aproximacoes.twice_around_the_tree import aproximacao_twice_around_the_tree
from heuristicas_aproximacoes.convex_hull import heuristica_convex_hull
from heuristicas_aproximacoes.gulosa import heuristica_gulosa


caminho = sys.argv[1]
nome, coords, matriz = ler_instancia_tsplib(caminho)
n = len(coords)

# lista que acumula todas as linhas de saida
linhas = []
linhas.append(f"Instancia: {nome} | Vertices: {n}")
linhas.append("")

heuristicas = {
    "trivial": lambda m, n: construir_solucao_trivial(n),
    "christofides": aproximacao_christofides,
    "twice_around_the_tree": aproximacao_twice_around_the_tree,
    "convex_hull": lambda m, n: heuristica_convex_hull(m, n, coords),
    "gulosa": heuristica_gulosa
}

for nomeH, heur in heuristicas.items():
    # instanciar metricas por heuristica (para comparar warm-start)
    contador = Contador()
    estatisticas = EstatisticasBAB()

    # 1. medir a heuristica
    with contador.medir(f"heuristica_{nomeH}"):
        rota = heur(matriz, n)

    # 2. medir o B&B (injeta contador e estatisticas)
    with contador.medir("bnb_total"):
        melhorRota, melhorCusto, custoInicial = branch_and_bound_tsplib(
            matriz, rota, contador=contador, estatisticas=estatisticas
        )

    # 3. acumular resultado + metricas
    linhas.append(f"{nomeH}:")
    linhas.append(f" rota inicial : {rota} | custo = {custoInicial}")
    linhas.append(f" melhor rota  : {melhorRota} | custo = {melhorCusto}")
    linhas.append(f" {estatisticas.relatorio()}")
    linhas.append(f" tempo heuristica = {contador.tempos.get(f'heuristica_{nomeH}', 0.0):.4f}s")
    linhas.append(f" tempo bnb        = {contador.tempos.get('bnb_total', 0.0):.4f}s")
    linhas.append(
        f" tempo bound      = {contador.tempos.get('bound_caminho', 0.0):.4f}s"
        f" | chamadas = {contador.chamadas.get('bound_caminho', 0)}"
    )
    linhas.append("")

# escrever tudo no arquivo .txt com nome baseado na instancia
nomeArquivoSaida = f"{nome}.txt"
with open(nomeArquivoSaida, "w") as arquivo:
    arquivo.write("\n".join(linhas))

print(f"Resultados salvos em {nomeArquivoSaida}")