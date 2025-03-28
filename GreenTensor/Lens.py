import math


class Lens:
    """
    Класс для представления линзы

    Атрибуты
    --------
    Radius : float
        Радиус линзы (коэффициент умножения pi)
    Layers_count : int
        Число слоев линзы (последний слой - воздух)
    Accuracy : int
        Точность расчетов
    Norm_radii : list[float]
        Нормированные радиусы слоев
    Dielectric_constants : list[float]
        Диэлектрическая проницаемость материала слоев
    Magnetic_permeabilities : list[float]
        Магнитная проницаемость материала слоев
    """
    def __init__(
        self,
        radiusRatio: int,
        layers_count: int,
        norm_radii: list[float],
        dielectric_constants: list[float],
        magnetic_permeabilities: list[float]
    ):
        """
        Инициализация линзы с валидацией параметров

        :param radius: Радиус линзы (коэффициент умножения pi)
        :param layers_count: Число слоев линзы (последний слой - воздух)
        :param norm_radii: Нормированные радиусы слоев
        :param dielectric_constants: Диэлектрическая проницаемость материала слоев
        :param magnetic_permeabilities: Магнитная проницаемость материала слоев
        """
        self._validate_inputs(
            radiusRatio,
            layers_count,
            norm_radii,
            dielectric_constants,
            magnetic_permeabilities
        )
        
        self.Radius: float = radiusRatio * math.pi
        self.Layers_count: int = layers_count
        self.Accuracy: int = math.ceil(self.Radius*2)
        self.Norm_radii: list[float] = norm_radii
        self.Dielectric_constants: list[float] = dielectric_constants
        self.Magnetic_permeabilities: list[float] = magnetic_permeabilities

    @staticmethod
    def _validate_inputs(
        radius: int,
        layers_count: int,
        norm_radii: list[float],
        dielectric_constants: list[float],
        magnetic_permeabilities: list[float]
    ) -> None:
        """Валидация всех входных параметров"""
        
        # Проверка типов данных
        if not isinstance(radius, int):
            raise TypeError("radius должен быть целым числом")
            
        if not isinstance(layers_count, int):
            raise TypeError("layers_count должен быть целым числом")
            
        if not all(isinstance(x, (float, int)) for x in norm_radii):
            raise TypeError("norm_radii должен содержать только числа")
            
        if not all(isinstance(x, (float, int)) for x in dielectric_constants):
            raise TypeError("dielectric_constants должен содержать только числа")
            
        if not all(isinstance(x, (float, int)) for x in magnetic_permeabilities):
            raise TypeError("magnetic_permeabilities должен содержать только числа")

        if radius <= 0:
            raise ValueError("radius должен быть положительным числом")
            
        if layers_count <= 0:
            raise ValueError("layers_count должен быть положительным числом")
            
        if len(norm_radii) != layers_count:
            raise ValueError("Количество norm_radii должно соответствовать layers_count")
            
        if len(dielectric_constants) != layers_count:
            raise ValueError("Количество dielectric_constants должно соответствовать layers_count")
            
        if len(magnetic_permeabilities) != layers_count:
            raise ValueError("Количество magnetic_permeabilities должно соответствовать layers_count")

        if not all(0 < r <= 1 for r in norm_radii):
            raise ValueError("Нормированные радиусы должны быть в диапазоне (0, 1]")
            
        if any(e <= 0 for e in dielectric_constants):
            raise ValueError("Диэлектрическая проницаемость должна быть положительной")
            
        if any(m <= 0 for m in magnetic_permeabilities):
            raise ValueError("Магнитная проницаемость должна быть положительной")

        if norm_radii[-1] != 1:
            raise ValueError("Последний слой (воздух) должен иметь радиус 1")
        
    def __str__(self) -> str:
        return (
            f"Lens Parameters:\n"
            f"Radius: {self.Radius:.2f}\n"
            f"Accuracy: {self.Accuracy}\n"
            f"Layers: {self.Layers_count}\n"
            f"Norm Radii: {self.Norm_radii}\n"
            f"Dielectric Constants: {self.Dielectric_constants}\n"
            f"Magnetic Permeabilities: {self.Magnetic_permeabilities}"
        )