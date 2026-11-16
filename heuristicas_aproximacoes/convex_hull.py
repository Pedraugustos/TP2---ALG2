def produto_vetorial(o, a, b):
    """
    Calcula o produto vetorial (componente z) dos vetores OA e OB.

    Entrada:
        o, a, b: pontos no formato (x, y).

    Saída:
        valor: positivo se O -> A -> B faz uma curva à esquerda (anti-horário),
        negativo se faz uma curva à direita e zero se os pontos são colineares.
    """

    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def fecho_convexo_monotone_chain(coordenadas):
    """
    Calcula o fecho convexo dos pontos pelo algoritmo Monotone Chain (Andrew),
    em O(n log n), usando apenas Python puro.

    Pontos colineares sobre as arestas do fecho são descartados; eles serão
    tratados depois como vértices internos e inseridos pela heurística.

    Entrada:
        coordenadas: dicionário {id: (x, y)}.

    Saída:
        fecho: lista de ids dos vértices do fecho em sentido anti-horário,
        sem repetir o primeiro no final.
    """

    # ordenar os ids por (x, y), com o id como desempate
    ids = sorted(coordenadas, key=lambda v: (coordenadas[v][0], coordenadas[v][1], v))

    # caso degenerado: menos de tres pontos
    if len(ids) <= 2:
        return ids

    # cadeia inferior
    inferior = []
    for v in ids:
        while len(inferior) >= 2 and produto_vetorial(
            coordenadas[inferior[-2]], coordenadas[inferior[-1]], coordenadas[v]
        ) <= 0:
            inferior.pop()
        inferior.append(v)

    # cadeia superior
    superior = []
    for v in reversed(ids):
        while len(superior) >= 2 and produto_vetorial(
            coordenadas[superior[-2]], coordenadas[superior[-1]], coordenadas[v]
        ) <= 0:
            superior.pop()
        superior.append(v)

    # juntar as cadeias removendo os extremos repetidos
    fecho = inferior[:-1] + superior[:-1]

    # remover ids repetidos (pontos com coordenadas iguais) mantendo a ordem
    vistos = set()
    fechoSemRepeticao = []
    for v in fecho:
        if v not in vistos:
            fechoSemRepeticao.append(v)
            vistos.add(v)

    return fechoSemRepeticao


def heuristica_convex_hull(matrizDistancias, n, coordenadas):
    """
    Heurística de Inserção Baseada em Fecho Convexo (Convex Hull) para o TSP.

    Monta uma sub-rota inicial com os vértices do fecho convexo e, em seguida,
    insere gulosamente os vértices internos: a cada passo, escolhe o par
    (vértice k, aresta (i, j) da sub-rota) de menor custo de inserção
    d(i, k) + d(k, j) - d(i, j), até que todos os vértices estejam no ciclo.

    Entrada:
        matrizDistancias: matriz (n+1) x (n+1) de distâncias (índice 0 ignorado).
        n: número de vértices (rotulados de 1 até n).
        coordenadas: dicionário {id: (x, y)} devolvido pelo leitor.

    Saída:
        cicloH: ciclo hamiltoniano começando no vértice 1, com o 1 repetido
        no final para fechar o ciclo.
    """

    # caso trivial: um unico vertice
    if n == 1:
        return [1, 1]

    # sub-rota inicial com os vertices do fecho convexo
    subRota = fecho_convexo_monotone_chain(coordenadas)

    # vertices internos que ainda precisam ser inseridos
    naoInseridos = set(range(1, n + 1)) - set(subRota)

    # insercao gulosa ate todos os vertices entrarem no ciclo
    while naoInseridos:
        melhorCusto = float("inf")
        melhorVertice = None
        melhorPosicao = None
        tamanho = len(subRota)

        # avaliar todo par (vertice interno k, aresta (i, j) da sub-rota)
        for k in naoInseridos:
            for posicao in range(tamanho):
                i = subRota[posicao]
                j = subRota[(posicao + 1) % tamanho]
                custoInsercao = (
                    matrizDistancias[i][k]
                    + matrizDistancias[k][j]
                    - matrizDistancias[i][j]
                )
                if custoInsercao < melhorCusto:
                    melhorCusto = custoInsercao
                    melhorVertice = k
                    melhorPosicao = posicao

        # inserir k entre i e j (logo apos a posicao de i)
        subRota.insert(melhorPosicao + 1, melhorVertice)
        naoInseridos.remove(melhorVertice)

    # rotacionar o ciclo para comecar no vertice 1
    indiceInicio = subRota.index(1)
    cicloH = subRota[indiceInicio:] + subRota[:indiceInicio]

    # fechar o ciclo repetindo o vertice inicial
    cicloH.append(cicloH[0])

    return cicloH