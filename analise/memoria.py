import tracemalloc
import gc

class MedidorMemoria:
    """
    Mede memória alocada por um bloco de código, usando tracemalloc.

    Uso:
        medidor = MedidorMemoria()
        with medidor.medir("minha_rotina"):
            rotina()
        print(medidor.relatorio())
    """

    def __init__(self):
        self.picos = {}

    def medir(self, rotulo):
        # implementacao como context manager
        return _ContextoMedicaoMemoria(self, rotulo)

    def _registrar(self, rotulo, pico_bytes):
        # converte para MB para facilitar a leitura
        self.picos[rotulo] = pico_bytes / (1024 * 1024)

    def relatorio(self):
        linhas = []
        for rotulo, pico_mb in self.picos.items():
            linhas.append(f"{rotulo:<20} | pico de memoria = {pico_mb:>8.2f} MB")
        return "\n".join(linhas)


class _ContextoMedicaoMemoria:
    def __init__(self, medidor, rotulo):
        self.medidor = medidor
        self.rotulo = rotulo
        self.inicio = None

    def __enter__(self):
        # limpar o lixo acumulado antes de comecar a medir memoria
        gc.collect()
        # inicia o rastreamento de memoria
        tracemalloc.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # pega o pico de memoria desde o inicio
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        # registra no medidor principal
        self.medidor._registrar(self.rotulo, peak)
        # nao suprime excecoes
        return False