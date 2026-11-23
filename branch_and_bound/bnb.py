import time

def calcular_custo_rota(rota, matrizDistancias):
    """
    Calcula o custo total de uma rota, considerando um ciclo.

    Entrada:
        rota: lista de vértices na ordem de visita, com o primeiro repetido no fim.
        matrizDistancias: matriz (n+1) x (n+1) de distâncias.

    Saída:
        custo: soma das distâncias do ciclo.
    """

    # somar as arestas consecutivas, incluindo o fechamento do ciclo
    return sum(matrizDistancias[rota[i]][rota[i + 1]] for i in range(len(rota) - 1))


def validar_rota_inicial(rotaInicial, n):
    """
    Verifica se a rota fornecida como limitante superior é um ciclo hamiltoniano
    válido e retorna sua representação fechada.

    Entrada:
        rotaInicial: ciclo fornecido pela heurística de inicialização.
        n: número de vértices da instância.

    Saída:
        rota: ciclo válido com n+1 posições.
    """

    # remover fechamento duplicado, caso exista
    rota = list(rotaInicial)
    if len(rota) == n + 1 and rota[0] == rota[-1]:
        rota = rota[:-1]

    # validar quantidade e unicidade dos vertices
    if len(rota) != n or len(set(rota)) != n:
        raise ValueError("rotaInicial deve conter cada vertice exatamente uma vez")

    if set(rota) != set(range(1, n + 1)):
        raise ValueError("rotaInicial deve conter os vertices de 1 a n")

    # fechar o ciclo sem alterar a ordem produzida pela heuristica
    return rota + [rota[0]]


def construir_vizinhos_ordenados(matrizDistancias, n):
    """
    Pré-calcula os vizinhos de cada vértice em ordem crescente de distância.

    A ordenação é feita uma única vez, evitando ordenar candidatos em todos
    os nós da árvore de busca.

    Entrada:
        matrizDistancias: matriz (n+1) x (n+1) de distâncias.
        n: número de vértices.

    Saída:
        vizinhosOrdenados: lista indexada pelo vértice, com vizinhos ordenados.
    """

    vizinhosOrdenados = [[] for _ in range(n + 1)]

    # ordenar uma vez os vizinhos de cada vertice
    for v in range(1, n + 1):
        vizinhosOrdenados[v] = sorted((u for u in range(1, n + 1) if u != v), key=lambda u: matrizDistancias[v][u])

    return vizinhosOrdenados


def bound_grau(matrizDistancias, n, caminho, custoFixo, visitados, vizinhosOrdenados):
    """
    Calcula um limitante inferior baseado nos graus que ainda precisam ser
    completados no caminho parcial.

    Cada vértice deve ter grau 2 no ciclo final. O método soma, para cada
    vértice, os custos das menores arestas ainda admissíveis necessárias
    para completar seu grau. Como cada aresta restante pode ser contabilizada
    por suas duas extremidades, divide-se essa soma por dois.


    Entrada:
        matrizDistancias: matriz (n+1) x (n+1) de distâncias.
        n: número de vértices.
        caminho: caminho parcial começando no vértice inicial.
        custoFixo: custo das arestas já escolhidas no caminho.
        visitados: conjunto dos vértices já presentes no caminho.
        vizinhosOrdenados: vizinhos pré-ordenados por distância.

    Saída:
        limiteInferior: custo mínimo estimado para completar o ciclo.
    """

    INF = float("inf")
    inicio = caminho[0]
    ultimo = caminho[-1]
    tamanhoCaminho = len(caminho)
    somaArestas = 0.0

    # calcular quantas arestas ainda faltam para cada vertice atingir grau 2
    for v in range(1, n + 1):
        if v in visitados:
            if tamanhoCaminho == 1:
                # na raiz, o vertice inicial ainda precisa de duas arestas
                faltam = 2
            elif v == inicio or v == ultimo:
                # as duas pontas do caminho precisam de uma aresta adicional
                faltam = 1
            else:
                # vertices internos ja possuem grau 2
                faltam = 0
        else:
            # vertices ainda nao visitados precisam de duas arestas
            faltam = 2

        if faltam == 0:
            continue

        encontrados = 0

        # aproveitar a ordenacao global e parar assim que encontrar as arestas necessarias
        for u in vizinhosOrdenados[v]:
            if matrizDistancias[v][u] == INF:
                continue

            if u in visitados:
                # vertices internos nao podem receber novas arestas
                if u != inicio and u != ultimo:
                    continue

                # nao permitir fechar o ciclo antes de visitar todos os vertices
                if tamanhoCaminho > 1 and v in (inicio, ultimo):
                    continue

            somaArestas += matrizDistancias[v][u]
            encontrados += 1

            if encontrados == faltam:
                break

        if encontrados < faltam:
            return INF

    # cada aresta ainda nao escolhida pode aparecer nas duas extremidades
    return custoFixo + somaArestas / 2.0


def _incrementar_estatistica(estatisticas, atributo, quantidade=1):
    """
    Atualiza uma métrica caso o objeto de estatísticas possua o atributo.

    Entrada:
        estatisticas: objeto opcional de métricas.
        atributo: nome do contador a atualizar.
        quantidade: incremento aplicado.

    Saída:
        Nenhuma.
    """

    # manter compatibilidade com diferentes objetos de estatisticas
    if estatisticas is not None and hasattr(estatisticas, atributo):
        setattr(estatisticas, atributo, getattr(estatisticas, atributo) + quantidade)


def branch_and_bound_tsplib(matrizDistancias, rotaInicial, contador=None, estatisticas=None, tempoLimite=None):
    """
    Resolve o TSP por Branch and Bound, utilizando a rota fornecida
    pela heurística do grupo como limitante superior inicial.

    A busca utiliza DFS, ordenação de candidatos por distância e um limitante
    inferior baseado nos graus ainda não completados. As arestas não são
    removidas: isso preserva a possibilidade de certificar a otimalidade.

    Entrada:
        matrizDistancias: matriz (n+1) x (n+1) de distâncias.
        rotaInicial: ciclo hamiltoniano produzido por uma heurística do grupo.
        contador: objeto opcional para medir subrotinas.
        estatisticas: objeto opcional para contar nós gerados, expandidos e podados.
        tempoLimite: limite de execução em segundos (None = sem limite).

    Saída:
        melhorRota: melhor ciclo encontrado.
        melhorCusto: custo da melhor rota encontrada.
        custoInicial: custo da rota heurística usada como warm-start.
        terminou: True se a otimalidade foi provada; False se houve timeout.
    """

    n = len(matrizDistancias) - 1
    if n < 3:
        raise ValueError("A instancia deve possuir pelo menos tres vertices.")

    # validar a solucao heuristica que sera usada como incumbente inicial
    rotaApertada = validar_rota_inicial(rotaInicial, n)

    melhorRota = rotaApertada[:]
    custoInicial = calcular_custo_rota(melhorRota, matrizDistancias)
    melhorCusto = custoInicial

    # preparar vizinhos ordenados uma unica vez
    vizinhosOrdenados = construir_vizinhos_ordenados(matrizDistancias, n)

    # preparar o controle de tempo e o estado inicial da busca
    prazo = None if tempoLimite is None else time.perf_counter() + tempoLimite
    inicio = melhorRota[0]
    caminho = [inicio]
    visitados = {inicio}
    nosProcessados = 0
    terminou = True

    _incrementar_estatistica(estatisticas, "nosGerados")

    def calcular_bound(caminhoAtual, custoAtual, visitadosAtual):
        # medir o tempo gasto no calculo do limitante inferior, se solicitado
        if contador is not None:
            with contador.medir("bound_caminho"):
                return bound_grau(matrizDistancias, n, caminhoAtual, custoAtual, visitadosAtual, vizinhosOrdenados)

        return bound_grau(matrizDistancias, n, caminhoAtual, custoAtual, visitadosAtual, vizinhosOrdenados)

    # verificar se a raiz pode ser descartada pelo limitante inferior
    boundRaiz = calcular_bound(caminho, 0.0, visitados)
    if boundRaiz >= melhorCusto:
        return melhorRota, melhorCusto, custoInicial, True

    def busca(custoAtual):
        """
        Explora recursivamente os caminhos parciais e atualiza a melhor rota.

        Entrada:
            custoAtual: custo das arestas fixadas no caminho corrente.

        Saída:
            True se a busca terminou; False se o tempo limite foi atingido.
        """

        nonlocal melhorRota, melhorCusto, nosProcessados

        nosProcessados += 1

        # verificar o prazo periodicamente para reduzir o custo do controle
        if prazo is not None and (nosProcessados & 127) == 0:
            if time.perf_counter() >= prazo:
                return False

        ultimo = caminho[-1]

        # ordenar implicitamente os filhos usando os vizinhos pre-ordenados
        for proximo in vizinhosOrdenados[ultimo]:
            if proximo in visitados:
                continue

            novoCusto = custoAtual + matrizDistancias[ultimo][proximo]

            # aplicar uma poda simples antes de calcular o limitante inferior
            if novoCusto >= melhorCusto:
                _incrementar_estatistica(estatisticas, "nosPodaAntecipada")
                _incrementar_estatistica(estatisticas, "nosPodados")
                continue

            caminho.append(proximo)
            visitados.add(proximo)
            _incrementar_estatistica(estatisticas, "nosGerados")

            if len(caminho) == n:
                # fechar o ciclo e atualizar o incumbente, se houver melhoria
                custoRota = novoCusto + matrizDistancias[proximo][inicio]
                _incrementar_estatistica(estatisticas, "nosFolha")

                if custoRota < melhorCusto:
                    melhorCusto = custoRota
                    melhorRota = caminho[:] + [inicio]

                visitados.remove(proximo)
                caminho.pop()
                continue

            # calcular o limitante inferior do novo caminho parcial
            novoBound = calcular_bound(caminho, novoCusto, visitados)

            if novoBound < melhorCusto:
                _incrementar_estatistica(estatisticas, "nosExpandidos")

                # continuar a busca em profundidade
                if not busca(novoCusto):
                    return False
            else:
                _incrementar_estatistica(estatisticas, "nosPodados")
                _incrementar_estatistica(estatisticas, "nosPodaAntecipada")

            # desfazer a escolha antes de testar o proximo filho
            visitados.remove(proximo)
            caminho.pop()

        return True

    # iniciar a exploracao a partir da raiz
    terminou = busca(0.0)

    return melhorRota, melhorCusto, custoInicial, terminou