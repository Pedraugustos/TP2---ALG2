import networkx as nx


def aproximacao_twice_around_the_tree(matrizDistancias, n):
    """
    Algoritmo aproximado de Twice-Around-the-Tree para o TSP com métrica.

    Entrada:
        matrizDistancias: matriz (n+1) x (n+1) de distâncias (índice 0 ignorado).
        n: número de vértices (vértices de 1 até n).

    Saída:
        cicloH: ciclo hamiltoniano aproximado, com o primeiro vértice
        repetido no fim para fechar o ciclo.
    """

    # caso trivial: um unico vertice
    if n == 1:
        return [1, 1]

    # construir o grafo completo nao direcionado com os vertices de 1 a n
    G = nx.Graph()
    G.add_nodes_from(range(1, n + 1))
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            G.add_edge(i, j, weight=matrizDistancias[i][j])

    # arvore geradora minima
    AGM = nx.minimum_spanning_tree(G, weight="weight")

    # duplicar as arestas da AGM
    multigrafo = nx.MultiGraph()
    multigrafo.add_nodes_from(G.nodes())
    multigrafo.add_edges_from(AGM.edges(data=True))
    multigrafo.add_edges_from(AGM.edges(data=True))

    # circuito euleriano partindo do vertice 1
    circuitoE = list(nx.eulerian_circuit(multigrafo, source=1))

    # transformar o circuito euleriano em um ciclo hamiltoniano (atalhos)
    visitados = set()
    cicloH = []

    for u, v in circuitoE:
        if u not in visitados:
            cicloH.append(u)
            visitados.add(u)

    # retornar ao vertice inicial
    cicloH.append(cicloH[0])

    return cicloH