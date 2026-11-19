import sys
import csv
import os

from analise.tempo import Contador, EstatisticasBAB
from analise.geometria import extrair_features
from leitor.leitor import ler_instancia_tsplib, ler_explicit
from branch_and_bound.BABbase import branch_and_bound_tsplib
from heuristicas_aproximacoes.christofides import aproximacao_christofides
from heuristicas_aproximacoes.trivial import construir_solucao_trivial
from heuristicas_aproximacoes.twice_around_the_tree import aproximacao_twice_around_the_tree
from heuristicas_aproximacoes.convex_hull import heuristica_convex_hull
from heuristicas_aproximacoes.gulosa import heuristica_gulosa


TEMPO_LIMITE = 30 * 60  # segundos (30 min por execucao do B&B)

caminho = sys.argv[1]
#nome, coords, matriz = ler_instancia_tsplib(caminho)
#n = len(coords)

nome, coords, matriz = ler_explicit(caminho)
n = len(matriz) - 1

# extrair features uma vez por instancia
features = extrair_features(coords, matriz, n)

print(f"Instancia: {nome} | Vertices: {n}")
print(f"features: {features}\n")

heuristicas = {
    "trivial": lambda m, n: construir_solucao_trivial(n),
    "christofides": aproximacao_christofides,
    "twice_around_the_tree": aproximacao_twice_around_the_tree,
    #"convex_hull": lambda m, n: heuristica_convex_hull(m, n, coords),
    "gulosa": heuristica_gulosa,
}

resultados = []

for nomeH, heur in heuristicas.items():
    contador = Contador()
    estatisticas = EstatisticasBAB()

    # 1. medir a heuristica
    with contador.medir(f"heuristica_{nomeH}"):
        rota = heur(matriz, n)

    # 2. medir o B&B (acrescenta contador e estatisticas)
    with contador.medir("bnb_total"):
        melhorRota, melhorCusto, custoInicial, terminou = branch_and_bound_tsplib(
            matriz, rota, contador=contador, estatisticas=estatisticas,
            tempoLimite=TEMPO_LIMITE
        )

    # gap de aproximacao: erro da heuristica em relacao ao otimo (melhorCusto do B&B)
    # se abortou por timeout, o otimo eh desconhecido: gap = NA
    if terminou:
        gap = (custoInicial - melhorCusto) / melhorCusto if melhorCusto > 0 else 0.0
    else:
        gap = "NA"

    # 3. montar dicionario com todas as metricas e features
    resultado = {
    "instancia": nome,
    "n": n,
    "heuristica": nomeH,
    "rotaHeuristica": rota,
    "rotaOtima": melhorRota if terminou else "NA",
    "custoInicial": custoInicial,
    "melhorCusto": melhorCusto if terminou else "NA",
    "gap": gap,
    "terminou": terminou,
    "tempoHeuristica": contador.tempos.get(f"heuristica_{nomeH}", 0.0),
    "tempoBNB": contador.tempos.get("bnb_total", 0.0),
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

    if terminou:
        print(f"{nomeH}: inicial={custoInicial} | melhor={melhorCusto} | gap={gap * 100:.2f}%")
    else:
        print(f"{nomeH}: inicial={custoInicial} | melhor=NA | gap=NA (timeout de {TEMPO_LIMITE}s)")
    print(f"  {estatisticas.relatorio()}")
    print(f"  tempo heur={resultado['tempoHeuristica']:.4f}s | tempo bnb={resultado['tempoBNB']:.4f}s")
    print()

# salvar CSV em <nome>_resultados/ (EX.: vinicius_instancias -> vinicius_resultados)
pastaInstancia = os.path.basename(os.path.dirname(os.path.abspath(caminho)))
nomePessoa = pastaInstancia.removesuffix("_instancias")
pastaResultados = f"{nomePessoa}_resultados"
os.makedirs(pastaResultados, exist_ok=True)

# nome do CSV baseado na instancia (sem .tsp)
nomeBase = nome.replace(".tsp", "")
nomeCsv = os.path.join(pastaResultados, f"{nomeBase}_{'_'.join(heuristicas.keys())}.csv")
with open(nomeCsv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=resultados[0].keys())
    writer.writeheader()
    writer.writerows(resultados)

print(f"Resultados salvos em {nomeCsv}")