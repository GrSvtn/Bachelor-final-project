!	Глобальные переменные для модуля PiezoMultylayered
!	
!	Точность eps также может быть задана через процедуру init
	
    Module GlobalsPiezoMultilayered
    use iso_fortran_env, only: dp => REAL64
	implicit none
	
		!	Входные параметры
	    integer									:: nL, nPhi							!	Количество слоёв, граница числа членов разложения функции q ([0; N])
		real (dp)								:: omega							!	Частота колебаний
		real (dp)								:: a, hO, EO, nuO, rhoO, d31, Ve	!	Длина, толщина, модуль Юнга, коэффициента Пуассона, плотность и пьезокерамическая постоянная накладки и напряжение
		real (dp), dimension (:), allocatable	:: h, E, nu, rhok					!	Толщины, модули Юнга, коэффициенты Пуассона и плотности слоёв
		
		!	Выходные параметры
		complex (dp), allocatable	:: qc (:)	!	Вектор коэффициентов в разложении q
	
    end module GlobalsPiezoMultilayered
