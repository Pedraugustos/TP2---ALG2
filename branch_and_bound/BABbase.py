import heapq

def calcular_custo_rota(rota, matrizDistancias):
    """
    Calcula o custo total de uma rota fechada (ciclo hamiltoniano).

    Entrada:
        rota: lista de vértices na ordem de visita.
        matrizDistancias: matriz n x n de distâncias.

    Saída:
        custo: soma das distâncias ao longo do ciclo fechado.
    """

    # somar arestas consecutivas
    custo = 0
    for i in range(len(rota) - 1):
        custo += matrizDistancias[rota[i]][rota[i + 1]]

    # fechar o ciclo voltando ao inicio
    custo += matrizDistancias[rota[-1]][rota[0]]
    return custo


def bound_caminho(matrizDistancias, n, caminho):
    """
    Calcula um limitante inferior, considerando as arestas já
    fixadas noh caminho e restringindo as candidatas de cada vértice.

    Regras por vértice:
      - primeiro do caminho: aresta fixada + menor aresta incidente que não
        seja a fixada e não vai para vértice interno (saturado);
      - último do caminho: mesma ideia, mas para a saída;
      - internos: já saturados (contribução noh custo fixo);
      - não visitados: as duas menores arestas incidentes permitidas.

    O bound e o custo fixo do caminho mais metade da soma das
    arestas faltantes (divisão por 2 corrige a contagem dupla).

    Entrada:
        matrizDistancias: matriz (n+1) x (n+1) de distâncias.
        n: número de vértices.
        caminho: lista de vértices já fixados (caminho parcial simples).

    Saída:
        bound: limitante inferior (float) ou infinito se inviável.
    """

    INF = float("inf")
    visitados = set(caminho)
    custoFixo = 0.0
    grau = {v: 0 for v in range(1, n + 1)}

    # contabilizar arestas ja fixadas noh caminho e seus graus
    for a, b in zip(caminho, caminho[1:]):
        grau[a] += 1
        grau[b] += 1
        custoFixo += matrizDistancias[a][b]

    # identificar as duas pontas do caminho
    pontas = set()
    if len(caminho) >= 1:
        pontas.add(caminho[0])
        pontas.add(caminho[-1])

    # somar as arestas minimas que faltam em cada vertice
    soma = 0.0
    for v in range(1, n + 1):
        faltam = 2 - grau[v]
        if faltam <= 0:
            continue

        # coletar arestas permitidas incidentes a v
        candidatos = []
        for u in range(1, n + 1):
            if u == v:
                continue
            if matrizDistancias[v][u] == INF:
                continue
            # proibir aresta para vertice visitado que nao seja ponta
            # (vertices internos estao saturados e nao podem receber mais)
            if u in visitados and u not in pontas:
                continue
            # proibir aresta entre as duas pontas antes do fechamento final
            # (evita fechar ciclo com menos de n vertices)
            if v in pontas and u in pontas and len(caminho) > 1 and len(caminho) < n:
                continue
            candidatos.append(matrizDistancias[v][u])

        # se nao ha arestas suficientes, caminho eh inviavel
        candidatos.sort()
        if len(candidatos) < faltam:
            return INF

        # somar as "faltam" menores arestas
        soma += sum(candidatos[:faltam])

    # dividir por 2 para corrigir contagem dupla das arestas
    return custoFixo + soma / 2.0


def branch_and_bound_tsplib(matrizDistancias, rotaInicial, contador=None, estatisticas=None):
    """
    Resolve o TSP por Branch and Bound.
    Usa a rota inicial recebida como limitante superior e o limitante
    de caminho (duas menores arestas por vértice, respeitando arestas
    já fixadas) como limitante inferior. Exploração é best-first: uma
    fila de prioridade (heap) expande sempre o nó de menor bound.
    Válido para instâncias métricas.

    Entrada:
        matrizDistancias: matriz (n+1) x (n+1) de distâncias euclidianas.
        rotaInicial: lista de n+1 vértices (primeiro repetido no fim) com
            uma rota fechada válida, usada como limitante superior inicial.
        contador: objeto Contador opcional; se fornecido, mede o tempo de
            bound_caminho (rotulo "bound_caminho") e do BAB (rotulo
            "bnb_total" e "copia_caminho").
        estatisticas: objeto EstatisticasBAB opcional; se fornecido, contabiliza
            nosGerados, nosExpandidos, nosPodados, nosFolha e
            nosPodaAntecipada.

    Saída:
        melhorRota: melhor ciclo hamiltoniano encontrado (fechado).
        melhorCusto: custo do melhor ciclo.
    """

    n = len(matrizDistancias) - 1
    INF = float("inf")
    inicio = 1

    # rota inicial obrigatoria como limitante superior
    melhorRota = rotaInicial[:]
    custoInicial = calcular_custo_rota(rotaInicial, matrizDistancias)
    melhorCusto = custoInicial

    # limitante inferior inicial do caminho raiz [inicio]
    caminhoRaiz = [inicio]
    if contador is not None:
        with contador.medir("bound_caminho"):
            boundRaiz = bound_caminho(matrizDistancias, n, caminhoRaiz)
    else:
        boundRaiz = bound_caminho(matrizDistancias, n, caminhoRaiz)

    # fila de prioridade (min-heap) de subproblemas:
    # (bound, nivel, contadorDesempate, custoFixo, caminhoParcial)
    # o contador de desempate evita comparar listas quando bounds empatam
    contadorDesempate = 0
    filaPrioridade = []
    
    heapq.heappush(filaPrioridade, (boundRaiz, len(caminhoRaiz), contadorDesempate, 0.0, caminhoRaiz))
    contadorDesempate += 1

    # raiz conta como noh gerado (foi inserida na fila)
    if estatisticas is not None:
        estatisticas.nosGerados += 1

    # expansao best-first: sempre retirar o noh de menor bound
    while filaPrioridade:
        boundAtual, nivelAtual, _, custoFixo, caminho = heapq.heappop(filaPrioridade)

        # noh retirado da fila: conta como expandido
        if estatisticas is not None:
            estatisticas.nosExpandidos += 1

        # poda: descartar noh se ja nao pode superar o melhor conhecido
        if boundAtual >= melhorCusto:
            if estatisticas is not None:
                estatisticas.nosPodados += 1
            continue

        # folha: caminho completo, fechar ciclo e comparar
        if nivelAtual == n:
            if estatisticas is not None:
                estatisticas.nosFolha += 1
            custoRota = custoFixo + matrizDistancias[caminho[-1]][inicio]
            if custoRota < melhorCusto:
                melhorCusto = custoRota
                melhorRota = caminho[:] + [inicio]
            continue

        # noh interno: ramificar para todos os proximos vertices nao visitados
        ultimo = caminho[-1]
        verticesRestantes = [v for v in range(1, n + 1) if v not in caminho]

        for proximo in verticesRestantes:
            # aresta inexistente: ignorar (defensivo; matriz e completa)
            if matrizDistancias[ultimo][proximo] == INF:
                continue

            novoCaminho = caminho + [proximo]
            novoCustoFixo = custoFixo + matrizDistancias[ultimo][proximo]

            if contador is not None:
                with contador.medir("bound_caminho"):
                    novoBound = bound_caminho(matrizDistancias, n, novoCaminho)
            else:
                novoBound = bound_caminho(matrizDistancias, n, novoCaminho)

            # so empilhar se o bound do filho ainda permite bater o melhor
            if novoBound < melhorCusto:
                heapq.heappush(
                    filaPrioridade,
                    (novoBound, len(novoCaminho), contadorDesempate, novoCustoFixo, novoCaminho),
                )
                contadorDesempate += 1
                if estatisticas is not None:
                    estatisticas.nosGerados += 1
            else:
                # poda antecipada: filho nao entra na fila
                if estatisticas is not None:
                    estatisticas.nosPodaAntecipada += 1

    return melhorRota, melhorCusto, custoInicial