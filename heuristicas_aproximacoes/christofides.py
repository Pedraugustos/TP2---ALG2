"""
O código deve ser da seguinte forma: nome_de_funcao (sneak case), nomeDeVariavel (camel case), NomeDeEstruturas (camel case com inicial maiúscula)
Após a declaração de cada função, dê uma breve descrição do que ela faz com texto bem formatado (acentuação, pontuação, etc) delimitado por """""".
Também dentro das função, ao realizar passos do algoritmo, dê uma breve descrição comentada com #, sem acentuação e iniciando com letra minúscula, a fim de evitar poluição textual, sobre o que o passo faz.
Aqui vai um exemplo de função implementada da forma correta, mas para o algoritmo de Christofides: 
"""
import networkx as nx


def aproximacao_christofides(matrizDistancias, n):
    """
    Algoritmo aproximado de Christofides para o TSP com métrica.
    Gera uma rota inicial (limitante superior) a partir da matriz de
    distâncias.

    Entrada:
        matrizDistancias: matriz (n+1) x (n+1) de distâncias (índice 0
            ignorado).
        n: número de vértices.

    Saída:
        cicloH: ciclo hamiltoniano aproximado (lista de n+1 vértices com
            o primeiro repetido no fim).
    """

    # construir grafo completo nao direcionado a partir da matriz
    G = nx.Graph()
    for i in range(1, n + 1):
        G.add_node(i)
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            G.add_edge(i, j, weight=matrizDistancias[i][j])

    # arvore geradora minima
    AGM = nx.minimum_spanning_tree(G, weight="weight")

    # conjunto de vertices de grau impar na AGM
    impares = [v for v in AGM.nodes() if AGM.degree(v) % 2 == 1]

    # subgrafo induzido pelos vertices impares
    H = G.subgraph(impares)

    # emparelhamento perfeito de custo minimo
    EP = nx.min_weight_matching(H, weight="weight")

    # combinar a AGM com o emparelhamento
    multigrafo = nx.MultiGraph()
    multigrafo.add_nodes_from(G.nodes())
    multigrafo.add_edges_from(AGM.edges(data=True))
    multigrafo.add_edges_from((u, v, G[u][v]) for u, v in EP)

    # circuito euleriano
    circuitoE = list(nx.eulerian_circuit(multigrafo))

    # transformar o circuito euleriano em um ciclo hamiltoniano
    visitados = set()
    cicloH = []

    for u, v in circuitoE:
        if u not in visitados:
            cicloH.append(u)
            visitados.add(u)

    # retornar ao vertice inicial
    cicloH.append(cicloH[0])

    return cicloH