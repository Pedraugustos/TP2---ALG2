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

    A parte inteira e o grau; a parte fracionaria (em centesimos de grau)
    representa os minutos: 38.24 -> 38 graus e 24 minutos.
    """

    graus = int(coord)
    minutos = coord - graus
    # 1 minuto = 1/60 grau; 5*minutos/3 e a conversao de DDD.MM para graus decimais
    return math.pi * (graus + 5.0 * minutos / 3.0) / 180.0


def _distancia_geo(xi, yi, xj, yj):
    """
    Distância geografica no padrão TSPLIB (formula oficial para GEO).
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
    d = 6378.388 * math.acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0

    # truncamento (int) e nao arredondamento
    return int(d)


def ler_instancia_tsplib(caminhoArquivo):
    """
    Lê uma instância TSPLIB (EUC_2D, GEO ou ATT) e retorna os dados básicos.

    Entrada:
        caminhoArquivo: string com o caminho para o arquivo .tsp.

    Saída:
        nomeInstancia: nome declarado na instância.
        coordenadas: dicionário {idVertice: (x, y)} com ids iniciando em 1.
        matrizDistancias: matriz (n+1) x (n+1) de distâncias no padrão
                          TSPLIB do tipo detectado (EDGE_WEIGHT_TYPE).
    """

    # abrir arquivo e ler todas as linhas
    with open(caminhoArquivo, "r") as arquivo:
        linhas = arquivo.readlines()

    nomeInstancia = ""
    tipoPeso = "EUC_2D" # padrao TSPLIB quando o cabecalho nao declara
    coordenadas = {} # dicionario do tipo {id: (x, y)}
    lendoSecao = False # flag para leitura

    # varrer linhas extraindo cabecalho e coordenadas
    for linha in linhas:
        linhaLimpa = linha.strip()

        # tipos de linha do cabecalho e das coordenadas
        # 1. cabecalho que declara o nome da instancia
        if linhaLimpa.startswith("NAME"):
            nomeInstancia = linhaLimpa.split(":")[-1].strip()
        # 2. cabecalho que declara o tipo de peso das arestas
        elif linhaLimpa.startswith("EDGE_WEIGHT_TYPE"):
            tipoPeso = linhaLimpa.split(":")[-1].strip()
        # 3. inicio de coordenadas de vertices
        elif linhaLimpa.startswith("NODE_COORD_SECTION"):
            lendoSecao = True
        # 4. fim de arquivo
        elif linhaLimpa.startswith("EOF"):
            break
        # 5. coordenadas de vertices
        elif lendoSecao and linhaLimpa:
            partes = linhaLimpa.split()

            if len(partes) >= 3:
                idVertice = int(partes[0])
                x = float(partes[1])
                y = float(partes[2])
                coordenadas[idVertice] = (x, y)

    # escolher a funcao de distancia conforme o tipo de peso
    if tipoPeso == "EUC_2D":
        funcaoDistancia = _distancia_euc_2d
    elif tipoPeso == "GEO":
        funcaoDistancia = _distancia_geo
    else:
        raise ValueError(f"EDGE_WEIGHT_TYPE não suportado: {tipoPeso}")

    # construir matriz de distancias arredondadas
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

def ler_explicit(caminhoArquivo):
    """
    Leitor específico para instâncias EXPLICIT (TSPLIB).

    A instância utiliza:
        - EDGE_WEIGHT_TYPE: EXPLICIT
        - EDGE_WEIGHT_FORMAT: LOWER_DIAG_ROW

    Retorna:
        nomeInstancia: nome declarado na instância.
        coordenadas: dicionário vazio (a instância não possui coordenadas).
        matrizDistancias: matriz (n+1) x (n+1) simétrica de distâncias.
    """

    # abrir arquivo e ler todas as linhas
    with open(caminhoArquivo, "r") as arquivo:
        linhas = arquivo.readlines()

    nomeInstancia = ""
    dimensao = None
    lendoPesos = False
    numerosPeso = []

    # varrer linhas
    for linha in linhas:
        linhaLimpa = linha.strip()

        if not linhaLimpa:
            continue

        # cabecalhos
        if linhaLimpa.startswith("NAME"):
            nomeInstancia = linhaLimpa.split(":")[-1].strip()

        elif linhaLimpa.startswith("DIMENSION"):
            dimensao = int(linhaLimpa.split(":")[-1].strip())

        # inicio da secao de pesos
        elif linhaLimpa.startswith("EDGE_WEIGHT_SECTION"):
            lendoPesos = True

        # fim do arquivo
        elif linhaLimpa.startswith("EOF"):
            break

        # leitura dos pesos (cada linha pode conter varios numeros)
        elif lendoPesos:
            for valor in linhaLimpa.split():
                numerosPeso.append(int(valor))

    # verificar se a dimensao foi encontrada
    if dimensao is None:
        raise ValueError("DIMENSION nao encontrada no arquivo.")

    n = dimensao

    # LOWER_DIAG_ROW possui n*(n+1)/2 valores
    quantidadeEsperada = n * (n + 1) // 2

    if len(numerosPeso) != quantidadeEsperada:
        raise ValueError(
            f"Quantidade incorreta de pesos: esperados {quantidadeEsperada}, encontrados {len(numerosPeso)}."
        )

    # construir matriz simetrica com o triangulo inferior
    # linha e coluna 0 nao sao analisadas
    matrizDistancias = [[0] * (n + 1) for _ in range(n + 1)]

    k = 0
    for i in range(1, n + 1):
        for j in range(1, i + 1):
            valor = numerosPeso[k]
            k += 1

            matrizDistancias[i][j] = valor
            matrizDistancias[j][i] = valor

    # a instancia EXPLICIT nao possui coordenadas
    coordenadas = {}

    return nomeInstancia, coordenadas, matrizDistancias