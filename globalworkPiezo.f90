!	Глобальные рабочие переменные для модуля PiezoMultylayered
	
    Module WorkGlobalsPiezoMultilayered
    use GlobalsPiezoMultilayered
	implicit none
	
		integer									:: mPsi, n0, m0, nm0, nP					!	Граница числа проекторов ([-M; M]), полные объёмы дискретизации (N0 = N + 1, M0 = 2 M + 1), их произведение, число полюсов
		integer, allocatable					:: indCh (:, :)                         	!	Матрица соответствия индексов
		real (dp)								:: ld, ud, pr, ub, eps, t1, t2, t3, t4		!	Для Dinn5~, Star5, Halfc
		real (dp)								:: b, EnuO, kappa, ka, kappa2, e0, epsHc  	!	Постоянная для накладки, -, волновое число с умножением на a и квадратом, погрешность вычисления полюсов
		real (dp), parameter					:: pi = dacos (-1._dp)
		real (dp), dimension (:), allocatable	:: vPk, vSk, alfaM, alfM2, GalM	!	Скорости продольных и поперечных волн, \alpha_m, \alpha_m^2 и 1 / (\kappa^2 - \alpha_m^2)
		complex (dp)							:: eC, epsC, C (2)              !	Комплексные значения погрешности, коэффициенты в представлении v
        complex (dp), allocatable				:: cinf (:)     				!	i^n
        complex (dp), allocatable				:: poles (:)            		!	Полюса
	
    end module WorkGlobalsPiezoMultilayered
