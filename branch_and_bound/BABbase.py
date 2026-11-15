import math
import heapq

def ler_instancia_tsplib(caminhoArquivo):
    """
    Lê uma instância TSPLIB no formato EUC_2D e retorna os dados básicos.

    Entrada:
        caminhoArquivo: string com o caminho para o arquivo .tsp.

    Saída:
        nomeInstancia: nome declarado na instância.
        coordenadas: dicionário {idVertice: (x, y)} com ids iniciando em 1.
        matrizDistancias: matriz (lista de listas) n x n com distâncias
                          euclidianas arredondadas para inteiro (padrão TSPLIB).
    """

    # abrir arquivo e ler todas as linhas
    with open(caminhoArquivo, "r") as arquivo:
        linhas = arquivo.readlines()

    nomeInstancia = ""
    coordenadas = {} # dicionario do tipo {id: (x, y)}
    lendoSecao = False # flag para leitura

    # varrer linhas extraindo cabecalho e coordenadas
    for linha in linhas:
        linhaLimpa = linha.strip()

        # quatro tipos de linha
        # 1. cabecalho que declara o nome da instancia
        if linhaLimpa.startswith("NAME"):
            nomeInstancia = linhaLimpa.split(":")[-1].strip()
        # 2. inicio de coordenadas de vertices
        elif linhaLimpa.startswith("NODE_COORD_SECTION"):
            lendoSecao = True
        # 3. fim de arquivo
        elif linhaLimpa.startswith("EOF"):
            break
        # 4. coordenadas de vertices
        elif lendoSecao and linhaLimpa:
            partes = linhaLimpa.split()

            if len(partes) >= 3:
                idVertice = int(partes[0])
                x = float(partes[1])
                y = float(partes[2])
                coordenadas[idVertice] = (x, y)

    # construir matriz de distancias euclidianas arredondadas
    n = len(coordenadas)
    # linha e coluna 0 nao sao analisadas
    matrizDistancias = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        xi, yi = coordenadas[i]
        for j in range(1, n + 1):
            if i == j:
                matrizDistancias[i][j] = 0
            else:
                xj, yj = coordenadas[j]
                d = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2)
                # arredondamento TSPLIB (nint)
                # [.0, .5) = 0
                # [.5, .0) = 1
                matrizDistancias[i][j] = int(d + 0.5)

    return nomeInstancia, coordenadas, matrizDistancias


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
    fixadas no caminho e restringindo as candidatas de cada vértice.

    Regras por vértice:
      - primeiro do caminho: aresta fixada + menor aresta incidente que não
        seja a fixada e não va para vertice interno (saturado);
      - último do caminho: mesma ideia, mas para a saída;
      - internos: já saturados (contribução no custo fixo);
      - não visitados: as duas menores arestas incidentes permitidas.

    O bound e o custo fixo do caminho mais metade da soma das arestas
    faltantes (a divisão por 2 corrige a contagem dupla).

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

    # contabilizar arestas ja fixadas no caminho e seus graus
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

        # se nao ha arestas suficientes, caminho e inviavel
        candidatos.sort()
        if len(candidatos) < faltam:
            return INF

        # somar as "faltam" menores arestas
        soma += sum(candidatos[:faltam])

    # dividir por 2 para corrigir contagem dupla das arestas
    return custoFixo + soma / 2.0


def bound_caminho(matrizDistancias, n, caminho):
    """
    Calcula um limitante inferior, considerando as arestas já
    fixadas no caminho e restringindo as candidatas de cada vértice.

    Regras por vértice:
      - primeiro do caminho: aresta fixada + menor aresta incidente que não
        seja a fixada e não va para vertice interno (saturado);
      - último do caminho: mesma ideia, mas para a saída;
      - internos: já saturados (contribução no custo fixo);
      - não visitados: as duas menores arestas incidentes permitidas.

    O bound e o custo fixo do caminho mais metade da soma das arestas
    faltantes (a divisão por 2 corrige a contagem dupla).

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

    # contabilizar arestas ja fixadas no caminho e seus graus
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

        # se nao ha arestas suficientes, caminho e inviavel
        candidatos.sort()
        if len(candidatos) < faltam:
            return INF

        # somar as "faltam" menores arestas
        soma += sum(candidatos[:faltam])

    # dividir por 2 para corrigir contagem dupla das arestas
    return custoFixo + soma / 2.0


def branch_and_bound_tsplib(matrizDistancias, rotaInicial):
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
            Obrigatória. Para cold start, passe a solução trivial
            1 -> 2 -> ... -> n -> 1.

    Saída:
        melhorRota: melhor ciclo hamiltoniano encontrado (fechado).
        melhorCusto: custo do melhor ciclo.
    """

    n = len(matrizDistancias) - 1
    INF = float("inf")
    inicio = 1

    # rota inicial obrigatoria como limitante superior
    melhorRota = rotaInicial[:]
    melhorCusto = calcular_custo_rota(rotaInicial, matrizDistancias)

    # limitante inferior inicial do caminho raiz [inicio]
    caminhoRaiz = [inicio]
    boundRaiz = bound_caminho(matrizDistancias, n, caminhoRaiz)

    # fila de prioridade (min-heap) de subproblemas:
    # (bound, nivel, custoFixo, caminhoParcial)
    # o contador de desempate evita comparar listas quando bounds empatam
    contadorDesempate = 0
    filaPrioridade = []
    heapq.heappush(
        filaPrioridade,
        (boundRaiz, len(caminhoRaiz), contadorDesempate, 0.0, caminhoRaiz),
    )
    contadorDesempate += 1

    # expansao best-first: sempre retirar o no de menor bound
    while filaPrioridade:
        boundAtual, nivelAtual, _, custoFixo, caminho = heapq.heappop(filaPrioridade)

        # poda: descartar no se ja nao pode superar o melhor conhecido
        if boundAtual >= melhorCusto:
            continue

        # folha: caminho completo, fechar ciclo e comparar
        if nivelAtual == n:
            custoRota = custoFixo + matrizDistancias[caminho[-1]][inicio]
            if custoRota < melhorCusto:
                melhorCusto = custoRota
                melhorRota = caminho[:] + [inicio]
            continue

        # no interno: ramificar para todos os proximos vertices nao visitados
        ultimo = caminho[-1]
        verticesRestantes = [v for v in range(1, n + 1) if v not in caminho]

        for proximo in verticesRestantes:
            # aresta inexistente: ignorar (defensivo; matriz e completa)
            if matrizDistancias[ultimo][proximo] == INF:
                continue

            novoCaminho = caminho + [proximo]
            novoCustoFixo = custoFixo + matrizDistancias[ultimo][proximo]
            novoBound = bound_caminho(matrizDistancias, n, novoCaminho)

            # so empilhar se o bound do filho ainda permite bater o melhor
            if novoBound < melhorCusto:
                heapq.heappush(
                    filaPrioridade,
                    (novoBound, len(novoCaminho), contadorDesempate, novoCustoFixo, novoCaminho),
                )
                contadorDesempate += 1

    return melhorRota, melhorCusto


def resolver_instancia_tsplib(caminhoArquivo, rotaInicial):
    """
    Função de conveniência: lê uma instância TSPLIB EUC_2D, executa o
    branch and bound e imprime/exibe o resultado.

    Entrada:
        caminhoArquivo: caminho para o arquivo .tsp.
        rotaInicial: lista de n+1 vértices (primeiro repetido no fim) com
            uma rota fechada válida, usada como limitante superior inicial.
            Obrigatória.

    Saída:
        melhorRota, melhorCusto: solução ótima encontrada.
    """

    # ler instancia
    nomeInstancia, coordenadas, matrizDistancias = ler_instancia_tsplib(caminhoArquivo)
    n = len(coordenadas)
    print(f"Instância: {nomeInstancia} | Vértices: {n}")

    # executar branch and bound com a rota inicial recebida
    melhorRota, melhorCusto = branch_and_bound_tsplib(matrizDistancias, rotaInicial)

    # exibir resultado
    print(f"Melhor rota: {melhorRota}")
    print(f"Melhor custo: {melhorCusto}")

    return melhorRota, melhorCusto