from typing import Final
import math

TETA_START: Final[float] = 0.01
TETA_STOP: Final[int] = 360
STEP: Final[float] = math.pi/180
STEPS: Final[int] = int(((abs(TETA_STOP) - abs(TETA_START))*(math.pi/180)) / STEP)