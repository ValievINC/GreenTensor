from .Lens import Lens
import cmath
import scipy
import math
from .Constants import STEPS, STEP

class LensCalculator:
    """
    Класс для представления калькулятора линзы

    Атрибуты
    --------
    Alfa : list[float]
        Углы отклонения для диэлектрических постоянных линзы.
    Beta : list[float]
        Углы отклонения для магнитных проницаемостей линзы.
    Etta : list[float]
        Волновые сопротивления для каждого слоя линзы.
    K : list[list[float]]
        Матрица коэффициентов среды k.
    J : list[list[float]]
        Значения функции Бесселя первого рода.
    Jder : list[list[float]]
        Производные функции Бесселя первого рода.
    N : list[list[float]]
        Значения функции Неймана.
    Nder : list[list[float]]
        Производные функции Неймана.
    C : list[list[float]]
        ***
    Cder : list[list[float]]
        Производные ***
    S : list[list[float]]
        ***
    Sder : list[list[float]]
        Производные ***
    Z : list[list[float]]
        Импедансы.
    Y : list[list[float]]
        Адмитансы.
    MJ : list[float]
        Значения функции Бесселя для радиуса линзы.
    MJder : list[float]
        Производные функции Бесселя для радиуса линзы.
    MH : list[float]
        Значения функции Ханкеля первого рода для радиуса линзы.
    MHder : list[float]
        Производные функции Ханкеля первого рода для радиуса линзы.
    Mn : list[float]
        ***
    Nn : list[float]
        ***
    Teta : list[float]
        Углы наблюдения в радианах.
    Cos_Teta : list[float]
        Косинусы углов наблюдения Teta.
    Pii : list[list[float]]
        Коэффициенты Pii, используемые в расчетах поляризационных характеристик.
    Tay : list[list[float]]
        Коэффициенты Tay, используемые в расчетах поляризационных характеристик.
    E_teta : list[list[float]]
        Электрическое поле на углах Teta.
    P_teta : list[list[float]]
        Поляризационное поле на углах Teta.
    Tetay : list[float]
        Углы наблюдения, нормализованные относительно количества шагов расчета.
    P_teta_max : float
        Максимальное значение поляризационного поля P_teta.
    DN_NORM : list[float]
        Нормированное значение диаграммы направленности.
    """
    def __init__(self, lens: Lens):
        """
        Инициализация калькулятора со всеми подсчётами для заданной линзы

        :param lens: Рассматриваемая линза
        """
        self.Alfa: list[float] = self.__get_Alpha(lens)
        self.Beta: list[float] = self.__get_Beta(lens)
        self.Etta: list[float] = self.__get_Etta(lens)
        self.K: list[list[float]] = self.__get_K(lens)
        self.J: list[list[float]] = self.__get_J(lens)
        self.Jder: list[list[float]] = self.__get_Jder(lens)
        self.N: list[list[float]] = self.__get_N(lens)
        self.Nder: list[list[float]] = self.__get_Nder(lens)
        self.C: list[list[float]] = self.__get_C(lens)
        self.Cder: list[list[float]] = self.__get_Cder(lens)
        self.S: list[list[float]] = self.__get_S(lens)
        self.Sder: list[list[float]] = self.__get_Sder(lens)
        self.Z: list[list[float]]
        self.Y: list[list[float]]
        self.Z, self.Y = self.__get_ZY(lens)
        self.MJ: list[float]
        self.MJder: list[float]
        self.MH: list[float]
        self.MHder: list[float]
        self.MJ, self.MJder, self.MH, self.MHder = self.__get_mLists(lens)
        self.Mn: list[float]
        self.Nn: list[float]
        self.Mn, self.Nn = self.__get_Mn_Nn(lens)
        self.Teta: list[float] = self.__get_Teta(lens)
        self.Cos_Teta: list[float] = self.__get_Cos_Teta(lens)
        self.Pii: list[list[float]]
        self.Tay: list[list[float]]
        self.Pii, self.Tay = self.__get_Pii_Tay(lens)
        self.E_teta: list[list[float]]
        self.P_teta: list[list[float]]
        self.E_teta, self.P_teta = self.__get_EP_teta(lens)
        self.Tetay: list[float] = self.__get_Tetay(lens)
        self.P_teta_max: float = self.__get_P_teta_max()
        self.DN_NORM: list[float] = self.__get_DN_NORM(lens)


    def __get_Alpha(self, lens: Lens) -> list[float]:
        return [math.atan(dc.imag / dc.real) for dc in lens.Dielectric_constants]

    def __get_Beta(self, lens: Lens) -> list[float]:
        return [math.atan(mp.imag / mp.real) for mp in lens.Magnetic_permeabilities]
    
    def __get_Etta(self, lens: Lens) -> list[float]:
        return [math.sqrt(
            math.fabs(lens.Dielectric_constants[i]) * math.fabs(lens.Magnetic_permeabilities[i]))
            for i in range(lens.Layers_count)]

    def __get_K(self, lens: Lens) -> list[list[float]]:
        result = [[0] * lens.Layers_count for i in range(lens.Layers_count)]
        j = 0
        for i in range (lens.Layers_count):
            result[i][j] = lens.Radius * lens.Norm_radii[i] * self.Etta[j]
            if j < lens.Layers_count - 1:
                j = j + 1
                result[i][j] = lens.Radius * lens.Norm_radii[i] * self.Etta[j]
        return result

    def __jfunc(self, i: int, j1: int, j2: int) -> float:
        nu = i + 1
        return (scipy.special.jv(nu + 0.5, self.K[j1][j2])) * (math.sqrt(self.K[j1][j2] * math.pi/2))

    def __get_J(self, lens: Lens) -> list[float]:
        result = [0 * lens.Layers_count for i in range(lens.Accuracy)]
        for i in range(lens.Accuracy):
            result[i] = self.__jfunc(i, 0, 0)
        return result

    def __jderfunc(self, i: int, j1: int, j2: int, tie: bool) -> float:
        nu = i + 1
        if not tie:
            jder = ((nu / (2 * nu + 1)) *  (scipy.special.jv(nu - 0.5, self.K[j1][j2]) * math.sqrt(self.K[j1][j2] * math.pi/2)) - \
            ((nu + 1) / (2 * nu + 1)) *  (scipy.special.jv(nu + 1.5, self.K[j1][j2]) * math.sqrt(self.K[j1][j2] * math.pi/2)) + \
            (self.J[i] / self.K[j1][j2]))
        else:
            jder = ((nu / (2 * nu + 1)) * ((scipy.special.jv(nu - 0.5,self.K[j1][j2]) * (math.sqrt(self.K[j1][j2] * math.pi/2))) / self.K[j1][j2])) * self.K[j1][j2] - \
            (((nu + 1) / (2 * nu + 1)) * ((scipy.special.jv(nu + 1.5,self.K[j1][j2])) * (math.sqrt(self.K[j1][j2] * math.pi/2))) / self.K[j1][j2]) * self.K[j1][j2] + \
            ((scipy.special.jv(nu + 0.5,self.K[j1][j2])) * (math.sqrt(self.K[j1][j2] * math.pi/2))) / self.K[j1][j2]
        return jder
    
    def __get_Jder(self, lens: Lens) -> list[float]:
        result = [0 * lens.Layers_count for i in range(lens.Accuracy)]
        for i in range(lens.Accuracy):
            result[i] = self.__jderfunc(i, 0, 0, False)
        return result
    
    def __nfunc(self, i: int, j1: int, j2: int) -> float:
        return scipy.special.yv((i+1) + 0.5, self.K[j1][j2]) * math.sqrt(self.K[j1][j2]* math.pi/2)
    
    def __get_N(self, lens: Lens) -> list[float]:
        result = [0 * lens.Layers_count for i in range(lens.Accuracy)]
        for i in range(lens.Accuracy):
            result[i] = self.__nfunc(i, 0, 0)
        return result
    
    def __nderfunc(self, i: int, j1: int, j2: int, tie: bool) -> float:
        nu = i + 1
        if not tie:
            nder = ((nu / (2 * nu + 1)) *  (scipy.special.yv(nu - 0.5, self.K[j1][j2]) * math.sqrt(self.K[j1][j2] * math.pi/2)) - \
            ((nu + 1) / (2 * nu + 1)) *  (scipy.special.yv(nu + 1.5, self.K[j1][j2]) * math.sqrt(self.K[j1][j2] * math.pi/2)) + \
            (self.__nfunc(i, j1, j2) / self.K[j1][j2]))
        else:
            nder = (((nu/(2 * nu + 1)) * (((scipy.special.yv(nu - 0.5, self.K[j1][j2])) * (math.sqrt(self.K[j1][j2] * math.pi/2))) / self.K[j1][j2])) * self.K[j1][j2] - \
            (((nu + 1) / (2 * nu + 1)) * ((scipy.special.yv(nu + 1.5, self.K[j1][j2])) * (math.sqrt(self.K[j1][j2] * math.pi/2))) / self.K[j1][j2]) * self.K[j1][j2] + \
            (self.__nfunc(i, j1, j2)) / self.K[j1][j2])
        return nder

    def __get_Nder(self, lens: Lens) -> list[float]:
        result = [0 * lens.Layers_count for i in range(lens.Accuracy)]
        for i in range(lens.Accuracy):
            result[i] = self.__nderfunc(i, 0, 0, False)
        return result
            
    def __get_C(self, lens: Lens) -> list[list[float]]:
        result = [[0] * (len(self.Etta)-1) for i in range(lens.Accuracy)]
        for i in range(lens.Accuracy - 1):
            for j in range(len(self.Etta) - 1):
                result[i][j] = (self.__jfunc(i, (j+1), (j+1)) * self.__nderfunc(i, j, (j+1), True)) \
                    - (self.__nfunc(i, (j+1), (j+1)) * self.__jderfunc(i, j, (j+1), True))
        return result
    
    def __get_Cder(self, lens: Lens) -> list[list[float]]:
        result = [[0] * (len(self.Etta)-1) for i in range(lens.Accuracy)]
        for i in range(lens.Accuracy - 1):
            for j in range(len(self.Etta) - 1):
                result[i][j] = self.__jderfunc(i, (j+1), (j+1), True) * self.__nderfunc(i, j, (j+1), True) \
                    - self.__nderfunc(i, (j+1), (j+1), True) * self.__jderfunc(i, j, (j+1), True)
        return result
    
    def __get_S(self, lens: Lens) -> list[list[float]]:
        result = [[0] * (len(self.Etta)-1) for i in range(lens.Accuracy)]
        for i in range(lens.Accuracy - 1):
            for j in range(len(self.Etta) - 1):
                result[i][j] = self.__nfunc(i,(j+1),(j+1)) * self.__jfunc(i,j,(j+1)) \
                    - self.__jfunc(i,(j+1),(j+1)) * self.__nfunc(i,j,(j+1))
        return result
    
    def __get_Sder(self, lens: Lens) -> list[list[float]]:
        result = [[0] * (len(self.Etta)-1) for i in range(lens.Accuracy)]
        for i in range(lens.Accuracy - 1):
            for j in range(len(self.Etta) - 1):
                result[i][j] = self.__nderfunc(i,(j+1),(j+1), True) * self.__jfunc(i, j, (j+1)) \
                    - self.__jderfunc(i,(j+1),(j+1), True) * self.__nfunc(i, j, (j+1))
        return result
    
    def __get_ZY(self, lens: Lens) \
            -> tuple[list[list[float]], list[list[float]]]:
        dc = lens.Dielectric_constants.copy()
        alpha = self.Alfa.copy()
        if dc[len(dc)-1] != (len(dc)-1):
            alpha.append(0)
            dc.append(len(dc))

        z = [[0] * (len(lens.Norm_radii)) for i in range(lens.Accuracy)]
        y = [[0] * (len(lens.Norm_radii)) for i in range(lens.Accuracy)]

        for i in range(lens.Accuracy - 1):
            for h in range(len(lens.Norm_radii)):
                if h == 0:
                    z[i][h] = (cmath.sqrt((cmath.exp(alpha[1] * 1j) * abs(dc[1])) / ((cmath.exp(alpha[0] * 1j) * abs(dc[0]))))) * ((self.Jder[i])/(self.J[i]))
                    y[i][h] = (cmath.sqrt((cmath.exp(alpha[0] * 1j) * abs(dc[0])) / ((cmath.exp(alpha[1] * 1j) * abs(dc[1]))))) * ((self.Jder[i])/(self.J[i]))
                else:
                    if h == (len(lens.Norm_radii) - 1):
                        z[i][h] = (cmath.sqrt((cmath.exp(alpha[h+1] * 1j) * abs(dc[h+1])) / ((cmath.exp(alpha[h] * 1j) * abs(dc[h]))))) * \
                                (self.Cder[i][h-1] + z[i][h-1] * self.Sder[i][h-1]) / (self.C[i][h-1] + z[i][h-1] * self.S[i][h-1]) / 2
                        y[i][h] = (cmath.sqrt((cmath.exp(alpha[h] * 1j) * abs(dc[h])) / ((cmath.exp(alpha[h+1] * 1j) * abs(dc[h+1]))))) * \
                                (self.Cder[i][h-1] + y[i][h-1] * self.Sder[i][h-1]) / (self.C[i][h-1] + y[i][h-1] * self.S[i][h-1]) * 2
                    else:
                        z[i][h] = (cmath.sqrt((cmath.exp(alpha[h+1] * 1j) * abs(dc[h+1])) / ((cmath.exp(alpha[h] * 1j) * abs(dc[h]))))) * \
                                (self.Cder[i][h-1] + z[i][h-1] * self.Sder[i][h-1]) / (self.C[i][h-1] + z[i][h-1] * self.S[i][h-1])
                        y[i][h] = (cmath.sqrt((cmath.exp(alpha[h] * 1j) * abs(dc[h])) / ((cmath.exp(alpha[h+1] * 1j) * abs(dc[h+1]))))) * \
                                (self.Cder[i][h-1] + y[i][h-1] * self.Sder[i][h-1]) / (self.C[i][h-1] + y[i][h-1] * self.S[i][h-1])   
        return z, y
        
    def __hfunc(self, i: int, rad: float) -> float:
        nu = i + 1
        return (scipy.special.hankel1(nu + 0.5, rad)) * (math.sqrt(rad * math.pi/2))

    def __hderfunc(self, i: int, rad: float) -> float:
        nu = i + 1
        return ((nu / (2 * nu + 1)) * (((scipy.special.hankel1(nu - 0.5, rad) * (cmath.sqrt(rad * math.pi/2))) / rad)) * rad - \
        (((nu + 1) / (2 * nu + 1)) * ((scipy.special.hankel1(nu + 1.5, rad)) * (cmath.sqrt(rad * math.pi/2))) / rad) * rad + \
        ((scipy.special.hankel1(nu + 0.5, rad)) * (cmath.sqrt(rad * math.pi/2))) / rad)

    def __get_mLists(self, lens: Lens) \
            -> tuple[list[float], list[float], list[float], list[float]]:        
        mJ = [0 * lens.Layers_count for i in range(lens.Accuracy)]
        mJder = [0 * lens.Layers_count for i in range(lens.Accuracy)]
        mH = [0 * lens.Layers_count for i in range(lens.Accuracy)]
        mHder = [0 * lens.Layers_count for i in range(lens.Accuracy)]
        
        original_k00 = self.K[0][0]
        self.K[0][0] = lens.Radius
        for i in range(lens.Accuracy):
            mJ[i] = self.__jfunc(i, 0, 0)
            mJder[i] = self.__jderfunc(i, 0, 0, True)
            mH[i] = self.__hfunc(i, lens.Radius)
            mHder[i] = self.__hderfunc(i, lens.Radius)
        self.K[0][0] = original_k00

        return mJ, mJder, mH, mHder
    
    def __get_Mn_Nn(self, lens: Lens) -> tuple[list[float], list[float]]:
        mn = [0 * lens.Layers_count for i in range(lens.Accuracy)]
        nn = [0 * lens.Layers_count for i in range(lens.Accuracy)]
        h = len(lens.Norm_radii) - 1
        for i in range(lens.Accuracy):
            mn[i] = (self.Z[i][h] * self.MJ[i] - self.MJder[i]) / (self.Z[i][h] * self.MH[i] - self.MHder[i])
            mn[i] = mn[i].real - mn[i].imag * 1j
            nn[i] = (self.Y[i][h] * self.MJ[i] - self.MJder[i]) / (self.Y[i][h] * self.MH[i] - self.MHder[i])
            nn[i] = nn[i].real - nn[i].imag * 1j
        return mn, nn
    
    def __get_Teta(self, lens: Lens) -> list[float]:
        teta_start = 0.01
        teta = [0 * lens.Layers_count for i in range(STEPS)]
        for i in range(STEPS):
            if i == 0:
                teta[i] = (teta_start)*(math.pi/180)
            else:
                teta[i] = teta[i-1] + STEP
        return teta
    
    def __get_Cos_Teta(self, lens: Lens) -> list[float]:
        cos_teta = [0 * lens.Layers_count for i in range(STEPS)]
        for i in range(STEPS):
            cos_teta[i] = math.cos(self.Teta[i])
        return cos_teta
    
    def __get_Pii_Tay(self, lens: Lens) -> tuple[list[list[float]], list[list[float]]]:
        M = [0 * lens.Layers_count for i in range(STEPS)]
        Lm0=[0 * lens.Layers_count for i in range(STEPS)]
        Lm1=[0 * lens.Layers_count for i in range(STEPS)]
        Lm2=[0 * lens.Layers_count for i in range(STEPS)]
        pii = [[0] * ((2*STEPS+1)) for i in range(STEPS+1)]
        tay = [[0] * ((2*STEPS+1)) for i in range(STEPS+1)]
        counter=0

        for i in range(lens.Accuracy):
            counter = counter + 1
            M = scipy.special.lpmv(0, counter, self.Cos_Teta)
            Lm0 = M
            M = scipy.special.lpmv(1, counter, self.Cos_Teta)
            Lm1 = M

            if counter<2:
                Lm2 = 0
            else:
                M = scipy.special.lpmv(2, counter, self.Cos_Teta)
                Lm2 = M

            for z in range(len(self.Teta)):
                if (self.Teta[z] > 0) & (self.Teta[z] < math.pi):
                    pii[i][z] = ((1)*Lm1[z])/(math.sin(self.Teta[z]))
                else:
                    if (self.Teta[z] > math.pi) & (self.Teta[z] < 2*math.pi):
                        pii[i][z] = ((-1)*Lm1[z])/(math.sin(self.Teta[z]))

            for z in range(len(self.Teta)):
                if counter<2:
                    tay[i][z] = (1/2)*(-counter*(counter+1)*Lm0[z])
                else:
                    tay[i][z] = (1/2)*(Lm2[z]-counter*(counter+1)*Lm0[z])
                    
        return pii, tay
    
    def __get_EP_teta(self, lens: Lens) -> tuple[list[list[float]], list[list[float]]]:
        E_teta= [[0] * (len(self.Teta)) for i in range(lens.Accuracy)]
        P_teta= [[0] * (len(self.Teta)) for i in range(lens.Accuracy)]
        y=0

        for z in range(len(self.Teta)):
            for p in range(lens.Accuracy):
                y=p+1
                E_teta[p][z]=((((2*y+1)/(y*(y+1)))*((-1)**y))*(self.Tay[p][z]*self.Mn[p]-self.Pii[p][z]*self.Nn[p]))

            for p in range(lens.Accuracy):
                y=p+1
                P_teta[0][z]=(P_teta[0][z]+E_teta[p][z])
                
            for p in range(lens.Accuracy):
                P_teta[0][z]=abs(P_teta[0][z])
            
        return E_teta, P_teta
    
    def __get_Tetay(self, lens: Lens) -> list[float]:
        tetay = [0 * lens.Layers_count for i in range(STEPS)]
        for i in range(len(self.Teta)):
            tetay[i] = self.Teta[i]*(STEPS/math.pi) 
        tetay.reverse()
        return tetay
    
    def __get_P_teta_max(self) -> float:
        P_teta_max = 0
        for i in range(len(self.Teta)):
            if self.P_teta[0][i] > P_teta_max:
                P_teta_max = self.P_teta[0][i]
        return P_teta_max
    
    def __get_DN_NORM(self, lens: Lens) -> list[float]:
        dn_norm = [0 * lens.Layers_count for i in range(len(self.Teta))]
        for i in range(len(self.Teta)):
            dn_norm[i] = 20*math.log10(self.P_teta[0][i]/self.P_teta_max)
        return dn_norm