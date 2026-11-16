import sys
import csv

from analise.tempo import Contador, EstatisticasBAB
from analise.geometria import extrair_features
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

# extrair features UMA VEZ por instancia (nao depende da heuristica)
features = extrair_features(coords, matriz, n)

print(f"Instancia: {nome} | Vertices: {n}")
print(f"features: {features}\n")

heuristicas = {
    #"trivial": lambda m, n: construir_solucao_trivial(n),
    "christofides": aproximacao_christofides,
    #"twice_around_the_tree": aproximacao_twice_around_the_tree,
    #"convex_hull": lambda m, n: heuristica_convex_hull(m, n, coords),
    #"gulosa": heuristica_gulosa,
}

resultados = []

for nomeH, heur in heuristicas.items():
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

    # 3. montar dicionario com TODAS as metricas + features
    resultado = {
        "instancia": nome,
        "n": n,
        "heuristica": nomeH,
        "custoInicial": custoInicial,
        "melhorCusto": melhorCusto,
        "tempoHeuristica": contador.tempos.get(f"heuristica_{nomeH}", 0.0),
        "tempoBAB": contador.tempos.get("bnb_total", 0.0),
        "tempoBound": contador.tempos.get("bound_caminho", 0.0),
        "chamadasBound": contador.chamadas.get("bound_caminho", 0),
        "nosGerados": estatisticas.nosGerados,
        "nosExpandidos": estatisticas.nosExpandidos,
        "nosPodados": estatisticas.nosPodados,
        "nosFolha": estatisticas.nosFolha,
        "nosPodaAntecipada": estatisticas.nosPodaAntecipada,
        **features,
    }
    resultados.append(resultado)

    print(f"{nomeH}: inicial={custoInicial} | melhor={melhorCusto}")
    print(f"  {estatisticas.relatorio()}")
    print(f"  tempo heur={resultado['tempoHeuristica']:.4f}s | tempo bnb={resultado['tempoBAB']:.4f}s")
    print()

# salvar CSV com nome baseado na instancia (sem .tsp)
nomeBase = nome.replace(".tsp", "")
nomeCsv = f"{nomeBase}.csv"
with open(nomeCsv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=resultados[0].keys())
    writer.writeheader()
    writer.writerows(resultados)

print(f"Resultados salvos em {nomeCsv}")