import time
from contextlib import contextmanager


class Contador:
    """
    Acumula tempo e contagem de chamadas por função.
    Uso:
        contador = Contador()
        with contador.medir("bound_caminho"):
            bound_caminho(...)
        print(contador.relatorio())
    """

    def __init__(self):
        self.tempos = {}    # tempo acumulado (s)
        self.chamadas = {}  # numero de chamadas

    @contextmanager
    def medir(self, rotulo):
        inicio = time.perf_counter()
        try:
            yield
        finally:
            duracao = time.perf_counter() - inicio
            self.tempos[rotulo] = self.tempos.get(rotulo, 0.0) + duracao
            self.chamadas[rotulo] = self.chamadas.get(rotulo, 0) + 1

    def relatorio(self):
        linhas = []
        for rotulo in self.tempos:
            linhas.append(
                f"{rotulo:<20} | tempo={self.tempos[rotulo]:>10.6f}s | chamadas={self.chamadas[rotulo]:>10}"
            )
        return "\n".join(linhas)


class EstatisticasBAB:
    """
    Contadores do Branch and Bound para a analise de warm-start.
    """

    def __init__(self):
        self.nosExpandidos = 0   # nohs retirados da fila e processados
        self.nosPodados = 0      # nohs descartados pela poda (bound >= melhor)
        self.nosGerados = 0      # filhos inseridos na fila
        self.nosFolha = 0        # nohs folha (caminho completo)
        self.nosPodaAntecipada = 0

    def relatorio(self):
        return (
            f"nosGerados={self.nosGerados} | nosExpandidos={self.nosExpandidos} | "
            f"nosPodados={self.nosPodados} | nosFolha={self.nosFolha}"
            f"nosPodaAntecipada={self.nosPodaAntecipada}"
        )