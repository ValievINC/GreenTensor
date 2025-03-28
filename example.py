from GreenTensor.Lens import Lens
from GreenTensor.LensCalculator import LensCalculator
from GreenTensor.LensPlotCreator import LensPlotCreator
import matplotlib.pyplot as plt

# Радиус линзы (коэффициент умножения на PI)
rad = 4
# Количество слоёв линзы (последний слой - воздух)
n = 6
# Нормированные радиусы слоев
a = [0.1, 0.2, 0.3, 0.4, 0.5, 1]
# Диэлектрическую проницаемость материала
eps = [10000000000000000000000000000000000000000, 10000000000000000000000000000000000000000, 10000000000000000000000000000000000000000, 10000000000000000000000000000000000000000, 10000000000000000000000000000000000000000, 1]
# Магнитную проницаемость материала
miy = [1, 1, 1, 1, 1, 1]

# Создаём экземпляр класса Lens
lens = Lens(rad, n, a, eps, miy)

# Делаем все необходимые расчёты с помощью LensCalculator
lensCalc = LensCalculator(lens)

# Создаём графики
fig_line, fig_polar = LensPlotCreator.create_plots(lensCalc)
plt.show()