def construir_solucao_trivial(n):
    """
    Constrói a solução ingênua (trivial) 1 -> 2 -> ... -> n -> 1.

    Entrada:
        n: número de vértices da instância.

    Saída:
        rota: lista [1, 2, ..., n, 1].
    """

    # rota sequencial com retorno ao inicio
    rota = list(range(1, n + 1)) + [1]
    return rota
