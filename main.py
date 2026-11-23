import sys
import csv
import os
from analise.tempo import Contador, EstatisticasBAB
from analise.geometria import extrair_features
from branch_and_bound.BABbase import branch_and_bound_tsplib
from leitor.leitor import ler_instancia_tsplib
from heuristicas_aproximacoes.christofides import aproximacao_christofides
from heuristicas_aproximacoes.convex_hull import heuristica_convex_hull
from heuristicas_aproximacoes.gulosa import heuristica_gulosa
from heuristicas_aproximacoes.trivial import construir_solucao_trivial
from heuristicas_aproximacoes.twice_around_the_tree import aproximacao_twice_around_the_tree
from analise.memoria import MedidorMemoria


# segundos (30 min por execucao do B&B)
TEMPO_LIMITE = 30 * 60  

caminho = sys.argv[1]
nome, coords, matriz = ler_instancia_tsplib(caminho)
n = len(matriz) - 1

# extrair features uma vez por instancia
features = extrair_features(coords, matriz, n)

print(f"Instancia: {nome} | Vertices: {n}")
print(f"Features: {features}\n")

# cada heuristica fornece a solucao inicial do B&B
heuristicas = {
    "trivial": lambda m, n: construir_solucao_trivial(n),
    "gulosa": heuristica_gulosa,
    #"convex_hull": lambda m, n: heuristica_convex_hull(m, n, coords),
    "twice_around_the_tree": aproximacao_twice_around_the_tree,
    "christofides": aproximacao_christofides,
}

resultados = []

for nomeH, heur in heuristicas.items():
    contador = Contador()
    estatisticas = EstatisticasBAB()
    medidorMem = MedidorMemoria()

    # 1. executar e medir a heuristica
    with contador.medir(f"heuristica_{nomeH}"):
        with medidorMem.medir(f"heuristica_{nomeH}"):
            rota = heur(matriz, n)

    # 2. executar o B&B usando a rota como limite superior inicial
    with contador.medir("bnb_total"):
        melhorRota, melhorCusto, custoInicial, terminou = branch_and_bound_tsplib(
            matriz,
            rota,
            contador=contador,
            estatisticas=estatisticas,
            tempoLimite=TEMPO_LIMITE,
        )

    # 3. calcular as metricas
    if terminou:
        gap = (custoInicial - melhorCusto) / melhorCusto if melhorCusto > 0 else 0.0
    else:
        gap = "NA"

    melhoriaPercentual = ((custoInicial - melhorCusto) / custoInicial * 100 if custoInicial > 0 else 0.0)

    resultado = {
        "instancia": nome,
        "n": n,
        "heuristica": nomeH,
        "custoInicial": custoInicial,
        "melhorCusto": melhorCusto,
        "rotaInicial": rota,    
        "melhorRota": melhorRota,
        "gap": gap,
        "melhoriaPercentual": melhoriaPercentual,
        "terminou": terminou,
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
        "picoMemHeuristica_MB": medidorMem.picos.get(f"heuristica_{nomeH}", 0.0)
    }

    resultados.append(resultado)

    # 4. exibir resultados
    status = "otimo encontrado" if terminou else "timeout"

    print(
        f"{nomeH}: inicial={custoInicial} | "
        f"melhor={melhorCusto} | "
        f"melhoria={melhoriaPercentual:.2f}% | "
        f"status={status}"
    )

    if terminou:
        print(f"Gap em relacao ao otimo: {gap * 100:.2f}%")

    print(f"{estatisticas.relatorio()}")
    print(
        f"tempo heur={resultado['tempoHeuristica']:.4f}s | "
        f"tempo B&B={resultado['tempoBAB']:.4f}s | "
        f"tempo bound={resultado['tempoBound']:.4f}s"
    )
    print()

# 5. salvar resultados
pastaInstancia = os.path.basename(os.path.dirname(os.path.abspath(caminho)))
nomePessoa = pastaInstancia.removesuffix("_instancias")
pastaResultados = f"{nomePessoa}_resultados"

os.makedirs(pastaResultados, exist_ok=True)

nomeBase = nome.replace(".tsp", "")
nomeHeuristicas = "_".join(heuristicas.keys())
nomeCsv = os.path.join(pastaResultados, f"{nomeBase}_{nomeHeuristicas}.csv")

with open(nomeCsv, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=resultados[0].keys())
    writer.writeheader()
    writer.writerows(resultados)

print(f"Resultados salvos em {nomeCsv}")
