import math


def _distancia_euc_2d(xi, yi, xj, yj):
    """
    Distância euclidiana 2D no padrão TSPLIB (arredondamento nint).
    """

    d = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2)
    # arredondamento TSPLIB (nint)
    # [.0, .5) = 0
    # [.5, .0) = 1
    return int(d + 0.5)


def _converter_para_radianos(coord):
    """
    Converte coordenada GEO (DDD.MM) da TSPLIB para radianos.

    A parte inteira é o grau; a parte fracionária (em centésimos de grau)
    representa os minutos: 38.24 = 38 graus e 24 minutos.
    """

    graus = int(coord)
    minutos = coord - graus
    # 1 minuto = 1/60 grau; 5*minutos/3 e a conversao de DDD.MM para graus decimais
    return math.pi * (graus + 5.0 * minutos / 3.0) / 180.0


def _distancia_geo(xi, yi, xj, yj):
    """
    Distância geográfica no padrão TSPLIB (fórmula oficial para GEO).
    """

    # converter latitudes (x) e longitudes (y) para radianos
    lat_i = _converter_para_radianos(xi)
    lon_i = _converter_para_radianos(yi)
    lat_j = _converter_para_radianos(xj)
    lon_j = _converter_para_radianos(yj)

    # formula da TSPLIB para GEO
    q1 = math.cos(lon_i - lon_j)
    q2 = math.cos(lat_i - lat_j)
    q3 = math.cos(lat_i + lat_j)
    # 6378.388 km eh o raio da terra
    d = 6378.388 * math.acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0

    # truncamento (int) e nao arredondamento
    return int(d)


def _distancia_att(xi, yi, xj, yj):
    """
    Distância quase euclidiana no padrão TSPLIB (formato ATT).

    A TSPLIB define ATT como:
        d = ceil( sqrt( (dx^2 + dy^2) / 10 ) )   se o valor for >= 1
        d = 0                                    se o valor for < 1
    """

    dx = xi - xj
    dy = yi - yj
    r = math.sqrt((dx * dx + dy * dy) / 10.0)

    # TSPLIB: arredonda para cima se >= 1; senao trunca
    if r >= 1.0:
        return int(math.ceil(r))
    return int(r)


def _ler_matriz_explicit(numeros, n, formato):
    """
    Constroi a matriz (n+1) x (n+1) simétrica a partir da lista de
    números da seção EDGE_WEIGHT_SECTION, conforme o formato.

    Formatos suportados:
      - LOWER_DIAG_ROW: triângulo inferior com diagonal, linha por linha
      - UPPER_DIAG_ROW: triângulo superior com diagonal, linha por linha
      - FULL_MATRIX: matriz completa n x n, linha por linha

    Entrada:
        numeros: lista de inteiros (todos os valores da seção).
        n: número de vértices.
        formato: string com o EDGE_WEIGHT_FORMAT.

    Saída:
        matrizDistancias: matriz (n+1) x (n+1) simétrica.
    """

    matriz = [[0] * (n + 1) for _ in range(n + 1)]

    if formato == "LOWER_DIAG_ROW":
        # ordem: (1,1), (2,1),(2,2), (3,1),(3,2),(3,3), ...
        k = 0
        for i in range(1, n + 1):
            for j in range(1, i + 1):
                valor = numeros[k]
                k += 1
                matriz[i][j] = valor
                matriz[j][i] = valor
        return matriz

    if formato == "UPPER_ROW":
        # ordem: (1,2),(1,3)..(1,n), (2,3)..(2,n), ...
        # (a diagonal principal nao e lida; assume-se 0)
        k = 0
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                valor = numeros[k]
                k += 1
                matriz[i][j] = valor
                matriz[j][i] = valor  # espelha para o triangulo inferior
        return matriz

    if formato == "UPPER_DIAG_ROW":
        # ordem: (1,1),(1,2),...,(1,n), (2,2),...,(2,n), ...
        k = 0
        for i in range(1, n + 1):
            for j in range(i, n + 1):
                valor = numeros[k]
                k += 1
                matriz[i][j] = valor
                matriz[j][i] = valor
        return matriz

    if formato == "FULL_MATRIX":
        # n x n, linha por linha
        k = 0
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                matriz[i][j] = numeros[k]
                k += 1
        return matriz

    raise ValueError(f"EDGE_WEIGHT_FORMAT nao suportado: {formato}")


def ler_instancia_tsplib(caminhoArquivo):
    """
    Lê uma instância TSPLIB (EUC_2D, GEO, ATT ou EXPLICIT) e retorna os
    dados básicos.

    Para EXPLICIT, suporta LOWER_DIAG_ROW, UPPER_DIAG_ROW e FULL_MATRIX.

    Entrada:
        caminhoArquivo: string com o caminho para o arquivo .tsp.

    Saída:
        nomeInstancia: nome declarado na instância.
        coordenadas: dicionário {idVertice: (x, y)} (vazio se EXPLICIT).
        matrizDistancias: matriz (n+1) x (n+1) de distâncias.
    """

    # abrir arquivo e ler todas as linhas
    with open(caminhoArquivo, "r") as arquivo:
        linhas = arquivo.readlines()

    nomeInstancia = ""
    tipoPeso = "EUC_2D"      # padrao TSPLIB quando o cabecalho nao declara
    formatoPeso = None       # so usado se tipoPeso == EXPLICIT
    dimensao = None          # lido do cabecalho DIMENSION
    coordenadas = {}         # dicionario do tipo {id: (x, y)}
    lendoCoordenadas = False # flag para a secao NODE_COORD_SECTION
    lendoPesos = False       # flag para a secao EDGE_WEIGHT_SECTION
    numerosPeso = []         # valores da secao EDGE_WEIGHT_SECTION

    # varrer linhas
    for linha in linhas:
        linhaLimpa = linha.strip()
        if not linhaLimpa:
            continue

        # cabecalhos
        if linhaLimpa.startswith("NAME"):
            nomeInstancia = linhaLimpa.split(":")[-1].strip()
        elif linhaLimpa.startswith("EDGE_WEIGHT_TYPE"):
            tipoPeso = linhaLimpa.split(":")[-1].strip()
        elif linhaLimpa.startswith("EDGE_WEIGHT_FORMAT"):
            formatoPeso = linhaLimpa.split(":")[-1].strip()
        elif linhaLimpa.startswith("DIMENSION"):
            dimensao = int(linhaLimpa.split(":")[-1].strip())
        elif linhaLimpa.startswith("NODE_COORD_SECTION"):
            lendoCoordenadas = True
            lendoPesos = False
        elif linhaLimpa.startswith("EDGE_WEIGHT_SECTION"):
            lendoPesos = True
            lendoCoordenadas = False
        elif linhaLimpa.startswith("EOF"):
            break
        elif linhaLimpa.startswith("DISPLAY_DATA_SECTION"):
            # coordenadas de visualizacao: nao usadas para calculo
            lendoCoordenadas = False
            lendoPesos = False
        elif lendoCoordenadas:
            partes = linhaLimpa.split()
            if len(partes) >= 3:
                idVertice = int(partes[0])
                x = float(partes[1])
                y = float(partes[2])
                coordenadas[idVertice] = (x, y)
        elif lendoPesos:
            # cada linha pode ter varios numeros
            for token in linhaLimpa.split():
                numerosPeso.append(int(token))

    # caso 1: EXPLICIT (matriz direta na secao EDGE_WEIGHT_SECTION)
    if tipoPeso == "EXPLICIT":
        if formatoPeso is None:
            raise ValueError("EDGE_WEIGHT_TYPE EXPLICIT sem EDGE_WEIGHT_FORMAT")
        n = dimensao if dimensao is not None else len(coordenadas)
        if n is None:
            raise ValueError("DIMENSION ausente e sem coordenadas para inferir n")
        matrizDistancias = _ler_matriz_explicit(numerosPeso, n, formatoPeso)
        return nomeInstancia, coordenadas, matrizDistancias

    # caso 2: coordenadas (EUC_2D, GEO, ATT) -> calcula a matriz
    if tipoPeso == "EUC_2D":
        funcaoDistancia = _distancia_euc_2d
    elif tipoPeso == "GEO":
        funcaoDistancia = _distancia_geo
    elif tipoPeso == "ATT":
        funcaoDistancia = _distancia_att
    else:
        raise ValueError(f"EDGE_WEIGHT_TYPE não suportado: {tipoPeso}")

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
                matrizDistancias[i][j] = funcaoDistancia(xi, yi, xj, yj)

    return nomeInstancia, coordenadas, matrizDistancias