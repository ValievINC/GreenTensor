import matplotlib.pyplot as plt
import math
from typing import Optional
from .LensCalculator import LensCalculator

class LensPlotCreator:
    """Класс для создания графиков"""

    @staticmethod
    def create_plots(lensCalc: LensCalculator) -> tuple[plt.Figure, Optional[plt.Figure]]:
        """
        Создание графиков для заданных расчётов
        
        :param lensCalc: Rалькулятор линзы с выполненными расчётами
        """
        fig_line = plt.figure()
        plt.plot(lensCalc.Tetay, lensCalc.DN_NORM, color='blue', linestyle='-', linewidth=2, label='Green_tensor')
        plt.grid(True)
        plt.legend()
        
        fig_polar, ax = plt.subplots(figsize=(4, 4), subplot_kw={'projection': 'polar'})
        teta_corrected = [angle - math.pi for angle in lensCalc.Teta]
        ax.plot(teta_corrected, lensCalc.DN_NORM, color='blue', linestyle='-', linewidth=1, label='Green_tensor')
        ax.legend(loc='upper right')
        return fig_line, fig_polar