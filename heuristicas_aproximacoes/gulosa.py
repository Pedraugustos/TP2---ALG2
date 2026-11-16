def heuristica_gulosa(matrizDistancias, n):
    """
    Heurística Gulosa do Vizinho Mais Próximo (Nearest Neighbor) para o TSP.

    Parte do vértice 1 e, a cada passo, vai para o vértice não visitado mais
    próximo do vértice atual, até visitar todos; por fim, retorna ao vértice 1.
    Em caso de empate, escolhe o candidato de menor identificador, o que torna
    o algoritmo determinístico.

    Entrada:
        matrizDistancias: matriz (n+1) x (n+1) de distâncias (índice 0 ignorado).
        n: número de vértices (rotulados de 1 até n).

    Saída:
        cicloH: ciclo hamiltoniano começando no vértice 1, com o 1 repetido
        no final para fechar o ciclo.
    """

    # iniciar a rota pelo vertice 1
    verticeAtual = 1
    cicloH = [verticeAtual]
    visitados = {verticeAtual}

    # repetir ate todos os vertices estarem na rota
    while len(cicloH) < n:
        melhorDistancia = float("inf")
        proximoVertice = None

        # buscar o vizinho nao visitado mais proximo (ordem crescente de id
        # com comparacao estrita garante desempate pelo menor identificador)
        for v in range(1, n + 1):
            if v in visitados:
                continue
            if matrizDistancias[verticeAtual][v] < melhorDistancia:
                melhorDistancia = matrizDistancias[verticeAtual][v]
                proximoVertice = v

        # aresta inexistente para todos os restantes: rota inviavel
        if proximoVertice is None:
            raise ValueError("nao ha aresta do vertice atual para nenhum vertice restante")

        # adicionar o vertice escolhido e atualizar o vertice atual
        cicloH.append(proximoVertice)
        visitados.add(proximoVertice)
        verticeAtual = proximoVertice

    # fechar o ciclo retornando ao vertice inicial
    cicloH.append(cicloH[0])

    return cicloH