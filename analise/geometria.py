import math
import networkx as nx

from heuristicas_aproximacoes.convex_hull import fecho_convexo_monotone_chain


def razao_de_aspecto(coordenadas):
    """
    Razão de aspecto do bounding box dos pontos.
    Retorna largura / altura (>= 1); 1.0 se altura for zero.
    """
    xs = [p[0] for p in coordenadas.values()]
    ys = [p[1] for p in coordenadas.values()]
    largura = max(xs) - min(xs)
    altura = max(ys) - min(ys)
    if altura == 0:
        return 1.0
    return largura / altura


def proporcao_fecho(coordenadas):
    """
    Proporção de pontos que pertencem ao fecho convexo (entre 0 e 1).
    """
    n = len(coordenadas)
    if n == 0:
        return 0.0
    fecho = fecho_convexo_monotone_chain(coordenadas)
    return len(fecho) / n


def entropia_espacial(coordenadas, bins=8):
    """
    Entropia de Shannon normalizada da distribuição espacial (histograma 2D).
    Quanto maior, mais espalhados; menor, mais aglomerados.
    """
    if len(coordenadas) <= 1:
        return 0.0

    xs = [p[0] for p in coordenadas.values()]
    ys = [p[1] for p in coordenadas.values()]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    largura = (xmax - xmin) if xmax > xmin else 1.0
    altura = (ymax - ymin) if ymax > ymin else 1.0

    histograma = [[0] * bins for _ in range(bins)]
    for x, y in coordenadas.values():
        ix = min(int((x - xmin) / largura * bins), bins - 1)
        iy = min(int((y - ymin) / altura * bins), bins - 1)
        histograma[ix][iy] += 1

    n = len(coordenadas)
    entropia = 0.0
    for linha in histograma:
        for contagem in linha:
            if contagem > 0:
                p = contagem / n
                entropia -= p * math.log(p)

    maxEntropia = math.log(bins * bins)
    return entropia / maxEntropia if maxEntropia > 0 else 0.0


def grau_maximo_agm(matrizDistancias, n):
    """
    Grau maximo de um no na Arvore Geradora Minima (presenca de hubs).
    """
    if n <= 1:
        return 0

    G = nx.Graph()
    G.add_nodes_from(range(1, n + 1))
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            G.add_edge(i, j, weight=matrizDistancias[i][j])

    AGM = nx.minimum_spanning_tree(G, weight="weight")
    return max(dict(AGM.degree()).values())


def cv_vizinho_mais_proximo(matrizDistancias, n):
    """
    Coeficiente de variação (desvio / média) das distancias ao vizinho
    mais proximo. Mede heterogeneidade local.
    """
    if n <= 1:
        return 0.0

    distanciasNN = []
    for i in range(1, n + 1):
        menor = min(matrizDistancias[i][j] for j in range(1, n + 1) if j != i)
        distanciasNN.append(menor)

    media = sum(distanciasNN) / n
    if media == 0:
        return 0.0

    variancia = sum((d - media) ** 2 for d in distanciasNN) / n
    desvio = math.sqrt(variancia)
    return desvio / media


def extrair_features(coordenadas, matrizDistancias, n):
    """
    Extrai todas as features geométricas/topológicas de uma instância.
    Retorna dicionario com as cinco features.
    """
    return {
        "aspectRatio": razao_de_aspecto(coordenadas),
        "propFecho": proporcao_fecho(coordenadas),
        "entropia": entropia_espacial(coordenadas),
        "grauMaxMST": grau_maximo_agm(matrizDistancias, n),
        "cvNN": cv_vizinho_mais_proximo(matrizDistancias, n),
    }