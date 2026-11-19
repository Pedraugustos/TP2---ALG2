import heapq
import time

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

'''
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
'''
def precalcular_vizinhos_ordenados(matrizDistancias, n):
    """
    Ordena os vizinhos de cada vértice uma única vez,
    do menor para o maior custo de aresta.

    Essa estrutura será reutilizada em todos os cálculos do bound.
    """
    vizinhosOrdenados = [[] for _ in range(n + 1)]

    for v in range(1, n + 1):
        vizinhos = [u for u in range(1, n + 1) if u != v]

        '''
        OTIMIZAÇÃO
        Ordenação de vizinhos mais baratos realizada apenas uma vez por vértice.
        '''
        vizinhos.sort(key=lambda u: matrizDistancias[v][u])

        vizinhosOrdenados[v] = vizinhos

    return vizinhosOrdenados


def bound_caminho(matrizDistancias, n, caminho, vizinhosOrdenados=None):
    """
    Calcula o limitante inferior considerando as arestas já fixadas.

    Reutiliza vizinhos previamente ordenados e seleciona somente
    as arestas permitidas necessárias, sem criar e ordenar listas
    de candidatos a cada chamada.
    """

    INF = float("inf")

    # se a funcao for chamada isoladamente, sem a estrutura pre-calculada,
    # ela ainda funciona. no B&B, passaremos a estrutura ja pronta.
    if vizinhosOrdenados is None:
        vizinhosOrdenados = precalcular_vizinhos_ordenados(
            matrizDistancias, n
        )

    '''
    OTIMIZAÇÃO:
    Para a otimização funcionar, substituimos set por lista.
    '''
    visitado = [False] * (n + 1)
    grau = [0] * (n + 1)

    for v in caminho:
        visitado[v] = True

    custoFixo = 0.0

    # contabilizar as arestas ja fixadas e seus graus.
    for a, b in zip(caminho, caminho[1:]):
        grau[a] += 1
        grau[b] += 1
        custoFixo += matrizDistancias[a][b]

    # as pontas do caminho sao os unicos vertices visitados
    # que ainda podem receber uma aresta adicional.
    primeiraPonta = caminho[0]
    ultimaPonta = caminho[-1]

    soma = 0.0

    for v in range(1, n + 1):

        faltam = 2 - grau[v]

        if faltam <= 0:
            continue

        '''
        OTIMIZACAO:
        Seleção direta das menores arestas permitidas, sem criar ou ordenar uma nova lista de candidatos.
        '''
        escolhidas = 0
        somaVertice = 0.0

        for u in vizinhosOrdenados[v]:

            distancia = matrizDistancias[v][u]

            if distancia == INF:
                continue

            # vertices internos do caminho ja estao saturados.
            if visitado[u] and u != primeiraPonta and u != ultimaPonta:
                continue

            # nao permitir fechar o ciclo antes de visitar todos os vertices.
            if v in (primeiraPonta, ultimaPonta) and u in (primeiraPonta, ultimaPonta) and len(caminho) > 1 and len(caminho) < n:
                continue

            somaVertice += distancia
            escolhidas += 1

            # ja encontramos todas as arestas necessarias.
            if escolhidas == faltam:
                break

        # se nao existem arestas suficientes, o caminho e inviavel.
        if escolhidas < faltam:
            return INF

        soma += somaVertice

    # divisao por 2 corrige a contagem dupla das arestas faltantes.
    return custoFixo + soma / 2.0

def branch_and_bound_tsplib(matrizDistancias, rotaInicial, contador=None, estatisticas=None, tempoLimite=None):
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
        tempoLimite: limite de tempo em segundos (opcional). Se estourar, a
            busca e interrompida e o ultimo valor retornado e False.

    Saída:
        melhorRota: melhor ciclo hamiltoniano encontrado (fechado).
        melhorCusto: custo do melhor ciclo.
    """

    n = len(matrizDistancias) - 1
    INF = float("inf")
    inicio = 1

    # controle de timeout: terminou=False indica busca abortada (otimo nao provado)
    prazo = None if tempoLimite is None else time.perf_counter() + tempoLimite
    terminou = True

    '''
    OTIMIZACAO:
    Pre-calculo dos vizinhos ordenados, realizado apenas uma vez.
    '''
    vizinhosOrdenados = precalcular_vizinhos_ordenados(matrizDistancias, n)

    # rota inicial obrigatoria como limitante superior
    melhorRota = rotaInicial[:]
    custoInicial = calcular_custo_rota(rotaInicial, matrizDistancias)
    melhorCusto = custoInicial

    # limitante inferior da raiz
    caminhoRaiz = [inicio]
    if contador is not None:
        with contador.medir("bound_caminho"):
            '''OTIMIZAÇÃO: vizinhosOrdenados já salva os melhores vizinhos'''
            boundRaiz = bound_caminho(matrizDistancias, n, caminhoRaiz, vizinhosOrdenados)
    else:
        '''OTIMIZAÇÃO: vizinhosOrdenados já salva os melhores vizinhos'''
        boundRaiz = bound_caminho(matrizDistancias, n, caminhoRaiz, vizinhosOrdenados)

    # fila de prioridade (min-heap) de subproblemas:
    # (bound, nivel, contadorDesempate, custoFixo, caminhoParcial)
    # o contador de desempate evita comparar listas quando bounds empatam
    contadorDesempate = 0
    filaPrioridade = []

    heapq.heappush(
        filaPrioridade,
        (boundRaiz, len(caminhoRaiz), contadorDesempate, 0.0, caminhoRaiz)
    )
    contadorDesempate += 1

    if estatisticas is not None:
        estatisticas.nosGerados += 1

    # expansao best-first: sempre retirar o noh de menor bound
    while filaPrioridade:
        # timeout: interromper a busca mantendo as estatisticas acumuladas
        if prazo is not None and time.perf_counter() > prazo:
            terminou = False
            break

        boundAtual, nivelAtual, _, custoFixo, caminho = heapq.heappop(filaPrioridade)

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

        ultimo = caminho[-1]
        verticesRestantes = [v for v in range(1, n + 1) if v not in caminho]

        for proximo in verticesRestantes:

            if prazo is not None and time.perf_counter() > prazo:
                terminou = False
                break

            if matrizDistancias[ultimo][proximo] == INF:
                continue

            novoCaminho = caminho + [proximo]
            novoCustoFixo = custoFixo + matrizDistancias[ultimo][proximo]

            if contador is not None:
                with contador.medir("bound_caminho"):
                    '''OTIMIZAÇÃO: vizinhosOrdenados já salva os melhores vizinhos'''
                    novoBound = bound_caminho(matrizDistancias, n, novoCaminho, vizinhosOrdenados)
            else:
                '''OTIMIZAÇÃO: vizinhosOrdenados já salva os melhores vizinhos'''
                novoBound = bound_caminho(matrizDistancias, n, novoCaminho, vizinhosOrdenados)

            # empilhar apenas se o bound do filho ainda permite bater o melhor
            if novoBound < melhorCusto:
                heapq.heappush(
                    filaPrioridade,
                    (novoBound, len(novoCaminho), contadorDesempate, novoCustoFixo, novoCaminho)
                )
                contadorDesempate += 1
                if estatisticas is not None:
                    estatisticas.nosGerados += 1

            # poda antecipada: filho nao entra na fila
            elif estatisticas is not None:
                estatisticas.nosPodaAntecipada += 1

        if not terminou:
            break

    return melhorRota, melhorCusto, custoInicial, terminou
'''
def branch_and_bound_tsplib(matrizDistancias, rotaInicial, contador=None, estatisticas=None, tempoLimite=None):
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
        tempoLimite: limite de tempo em segundos (opcional). Se estourar, a
            busca e interrompida e o ultimo valor retornado e False.

    Saída:
        melhorRota: melhor ciclo hamiltoniano encontrado (fechado).
        melhorCusto: custo do melhor ciclo.
    """

    n = len(matrizDistancias) - 1
    INF = float("inf")
    inicio = 1

    # controle de timeout: terminou=False indica busca abortada (otimo nao provado)
    prazo = None if tempoLimite is None else time.perf_counter() + tempoLimite
    terminou = True

    # rota inicial obrigatoria como limitante superior
    melhorRota = rotaInicial[:]
    custoInicial = calcular_custo_rota(rotaInicial, matrizDistancias)
    melhorCusto = custoInicial

    """
    OTIMIZAÇÃO:
    Ordenar os vizinhos de cada vertice uma unica vez, do menor para o maior custo de deslocamento.
    """
    vizinhosOrdenados = []

    for i in range(n + 1):
        vizinhos = [j for j in range(1, n + 1) if j != i]

        vizinhos.sort(key=lambda j: matrizDistancias[i][j])

        vizinhosOrdenados.append(vizinhos)

    # limitante inferior da raiz
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

    heapq.heappush(
        filaPrioridade,
        (boundRaiz, len(caminhoRaiz), contadorDesempate, 0.0, caminhoRaiz)
    )
    contadorDesempate += 1

    if estatisticas is not None:
        estatisticas.nosGerados += 1

    # expansao best-first: sempre retirar o noh de menor bound
    while filaPrioridade:
        # timeout: interromper a busca mantendo as estatisticas acumuladas
        if prazo is not None and time.perf_counter() > prazo:
            terminou = False
            break

        boundAtual, nivelAtual, _, custoFixo, caminho = heapq.heappop(filaPrioridade)

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

        """
        OTIMIZAÇÃO
        Percorrer os vizinhos do vértice atual em ordem de menor peso.
        """
        ultimo = caminho[-1]

        for proximo in vizinhosOrdenados[ultimo]:

            # ignorar vertices ja visitados
            if proximo in caminho:
                continue

            if prazo is not None and time.perf_counter() > prazo:
                terminou = False
                break

            if matrizDistancias[ultimo][proximo] == INF:
                continue

            novoCaminho = caminho + [proximo]
            novoCustoFixo = custoFixo + matrizDistancias[ultimo][proximo]

            if contador is not None:
                with contador.medir("bound_caminho"):
                    novoBound = bound_caminho(matrizDistancias, n, novoCaminho)
            else:
                novoBound = bound_caminho(matrizDistancias, n, novoCaminho)

            # empilhar apenas se o bound do filho ainda permite bater o melhor
            if novoBound < melhorCusto:
                heapq.heappush(
                    filaPrioridade,
                    (novoBound, len(novoCaminho), contadorDesempate, novoCustoFixo, novoCaminho)
                )
                contadorDesempate += 1
                if estatisticas is not None:
                    estatisticas.nosGerados += 1

            # poda antecipada: filho nao entra na fila
            elif estatisticas is not None:
                estatisticas.nosPodaAntecipada += 1

        if not terminou:
            break

    return melhorRota, melhorCusto, custoInicial, terminou
'''