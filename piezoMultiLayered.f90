!	Модуль поиска коэффициентов разложения функции контактной нагрузки q взаимодействия многослойного плоского изотропного волновода
!		и пьезокерамического актуатора, горизонтальные смещения которого v также отыскиваются модулем
	
	module PiezoMultiLayered
	use WorkGlobalsPiezoMultilayered
	use functions, only: SLAE, BesselJ, ones
    use subroutines, only: HalfcHandler
	use Global_Mult
	implicit none
	private
	public dispcur, qCoeffs, alloc, allocW, Q, vFunc, uv

	contains
	
	
	!===============================================		Коэффициенты в разложении q		===============================================
	
	subroutine qCoeffs
	
		integer										:: i, j
		integer, allocatable						:: oddInd (:)
		complex (dp)								:: kapc, Ec, Cvec (2), cin (0:3), Cmatr (2, 2)
		complex (dp), allocatable					:: Atemp (:)
		complex (dp), dimension (:, :), allocatable	:: Amatr, f, qk
		
		call init
        if (.not. allocated (qc)) allocate (qc (n0))
		allocate (oddInd (mPsi), Atemp (nm0), qk (2, n0), f (2, m0), Amatr (n0, m0))
		
		oddInd = [(i, i = 2 - mod (mPsi, 2), m0, 2)]
		cin = [complex (dp)	:: 1, ci, -1, -ci]
        cinf = [(cin (mod (i, 4)), i = 0, nPhi)]
		f = 2 * sin (ka) * spread (GalM, 1, 2)
        f (:, oddInd) = - f (:, oddInd)
        f (1, :) = kappa * f (1, :)
        f (2, :) = ci * alfaM * f (2, :)
		
		call dinn5 (Fanm, t1, t2, t3, t4, ld, ud, eps, pr, ub, nm0, Atemp)
		Amatr = reshape (Atemp, shape (Amatr))
        
		Amatr (:, oddind) = - Amatr (:, oddind)
		deallocate (Atemp, oddInd)
		
		do i = 1, n0
			Amatr (i, :) = Amatr (i, :) * cinf (i)
		end do
		Amatr = 2 * a * Amatr
        
        f = transpose (f)
		qk = transpose (SLAE (transpose (Amatr), f, m0, n0, 2))	!	Затуп с матрицей A - она должна быть транспонированной
		deallocate (f, Amatr)
        
        kapc = kappa
        Cmatr = reshape ([((Q (i * kapc, qk (j, :)), i = 1, -1, -2), j = 1, 2)], shape (Cmatr))
        Cmatr = matmul (reshape ([1, 1, 1, -1], shape (Cmatr)), Cmatr)
        Cmatr = b * Cmatr * exp (ci * ka) / (4 * kappa)
        Cmatr (1, :) = Cmatr (1, :) / sin (ka)
        Cmatr (2, :) = Cmatr (2, :) / cos (ka)
        
        Ec = e0 / (kappa * cos (ka))
        Cvec = - Ec * Cmatr (:, 2)
        Cmatr = Cmatr - ones (2)
        
        C = SLAE (Cmatr, Cvec, 2)
        C = C + [complex :: 0._dp, Ec]
        
		qc = matmul (C, qk)
				
	end subroutine qCoeffs
	
	
	!===============================================		Подынтегральное выражение элементов матрицы A		===============================================
	
	subroutine Fanm (alpha, Amatr, nm)
	
		integer			:: nm, i, j
		complex (dp)	:: alpha, ma, bj (0:nPhi), Gala (-mPsi:mPsi), Amatr (nm), alf (0:1, -mPsi:mPsi)
		
		if (abs (alpha) .lt. 1e-6) then
            Amatr = 0._dp
			return
		end if
		
		alf = alpha
        alf (1, :) = alfaM
        if (any (abs (alpha - alfaM) .lt. epsHc)) pause 'close to \alpha_m'
        Gala = 1 / (alpha * alpha - alfM2)
		
		ma = sin (a * alpha) * (Kalf (alpha) - b / (kappa2 - alpha * alpha))
		bj = BesselJ (0, a * alpha, n0)
        Amatr = ma * [((alf (mod (i, 2), j) * bj (i) * Gala (j), i = 0, nPhi), j = -mPsi, mPsi)]
        
        
	end subroutine Fanm
	
	
	!===============================================		Матрица Грина		===============================================
	
	complex (dp) function Kalf (alfa)
    
        complex (dp)				:: K
        complex (dp), intent (in)	:: alfa
        
		al = alfa
        call MultiK_Uni
		call form_K2 (K)
		
		Kalf = K
		
	end function Kalf
	
	
	!===============================================		Q (\alpha)		===============================================
	
	complex (dp) function Q (alpha, qcur, mcur)
    
        complex (dp), dimension (n0)        :: coefs, mults
        complex (dp), intent (in)           :: alpha
        complex (dp), optional, intent (in)	:: qcur (n0), mcur (n0)
        
        coefs = qc
        if (present (qcur)) coefs = qcur
        
        mults = cinf
        if (present (mcur)) mults = mcur
        
        Q = a * pi * sum (coefs * mults * BesselJ (0, a * alpha, n0))
        
	end function Q
	
	
	!===============================================		Значение v в наборе точек		===============================================
	
	subroutine vFunc (x, v)
	
		integer 				    :: i, j, n
        real (dp)               	:: t0, x (:)
		complex (dp), allocatable	:: v (:), csc (:, :), csf (:, :)
		
        n = size (x)
        allocate (v (n), csc (n, 2), csf (n, 0:nPhi))
        t0 = kappa / 2
		call dinn5 (Fv, t0, t0, t0, kappa + 1, ld, ud, eps, pr, ub, n, v)
        
        v = a * b * v / (2 * pi) + dot_product (C, [cos (x * kappa), sin (x * kappa)])
        deallocate (csc, csf)
        
        
		contains
    
		subroutine Fv (alpha, res, m)
        
			integer			:: m
            complex (dp)	:: alpha, res (m)
            
            csc = reshape ([cos (x * alpha), -ci * sin (x * alpha)], shape (csc))
            csf = reshape ([(csc (:, 2 - mod (i, 2)), i = 0, nPhi)], shape (csf))
            
            res = [(Q (alpha, mcur = cinf * csf (i, :)), i = 1, n)] / (kappa2 - alpha * alpha)
    
		end subroutine Fv
		
	end subroutine vFunc
	
	
	!===============================================		Значение v в наборе точек		===============================================
	
	subroutine uv (x, u)
	
		integer 				    :: i, j, n
        real (dp)               	:: x (:)
		complex (dp), allocatable	:: u (:)
		
        n = size (x)
		call dinn5 (Fuv, t1, t2, t3, t4, ld, ud, eps, pr, ub, n, u)
        
        u = u / (2 * pi)
        
        
		contains
    
		subroutine Fuv (alpha, res, m)
        
			integer			:: m
            complex (dp)	:: alpha, res (m)
            
            res = Kalf (alpha) * (Q (alpha) * exp (-ci * alpha * x) + Q (-alpha) * exp (ci * alpha * x))
    
		end subroutine Fuv
		
	end subroutine uv
	
	
	!===============================================		Дисперсионные кривые		===============================================
	
	subroutine dispcur
	use consoleMonitoring, only: initialize, retention
    
		integer		:: i, j, k, nom, file, nMax
        real (dp)	:: do_, o1, o2, omc
        
        call init (key = 1)
        !b = 0
        
        open (newunit = file, file = 'data/dc/poles.dat', action = 'write')
        
        k = 0
        nMax = 15
		o1 = 0
        o2 = 10
        do_ = 1e-3_dp
        
        nom = nint ((o2 - o1) / do_)
        call initialize ('Dispersion curves', o1, o2)
        do i = 0, nom
            omc = o1 + i * do_
            call pinit (omc)
			write (file, '(f12.6, <nMax-1> (", ", f12.6))') poles % re, [(0._dp, j = 1, nMax - nP)]
            deallocate (poles)
            call retention (omc, k, 10)
        end do
        
        close (file)
        
	end subroutine dispcur
	
	
	!===============================================		Инициализация		===============================================
	
	subroutine init (eps_, key)
	
		integer								:: i, j
        integer, optional					:: key
        real (dp), allocatable				:: mu (:)
		real (dp), intent (in), optional	:: eps_
		
        if (.not. allocated (cinf)) call allocW
        
		ld = 1e-2_dp	!	0.01	!	Dinn5 ~
        pr = 1e-2_dp	!	0.01
        ub = 15e2_dp	!	1500
        ud = 0._dp		!	0
    	eps = 1e-3_dp	!	0.001
		if (present (eps_)) eps = eps_
		
		epsHc = 1e-8_dp			!	0.00000001
		epsC = dcmplx (eps, 0)	!	(0.001; 0)
        eC = epsC * 1e-2_dp		!	(0.00001; 0)
		
		n0 = nPhi + 1
		m0 = 2 * mPsi + 1
		nm0 = n0 * m0
		
		alfaM = pi * [(i, i = -mPsi, mPsi)] / a
		alfM2 = alfaM * alfaM
	
		EnuO = (1 - nuO * nuO) / EO
		b = EnuO / hO
		e0 = (1	+ nuO) * d31 * Ve / hO
		
        allocate (mu (nL))
        mu = E / (2 * (1 + nu))
		vSk = sqrt (mu / rhok)
		vPk = vSk * sqrt (2 * (1 - nu) / (1 - 2 * nu))
		deallocate (mu)
        
		Ns = nL	!	Global_Mult
		iK = 2
		iY = 0
		iN = 1
		Vp = vPk
		Vs = vSk
		rho = rhok
		hs = h
		z = 0
        
		if (.not. present (key)) call pinit (omega)
		
	end subroutine init
	
	
	!===============================================		Поиск полюсов		===============================================
	
	subroutine pinit (om)
    
		integer         			:: i, nMax, k, file
        real (dp)					:: om, ome, do_
		!complex (dp), allocatable   :: t_ (:)!, poles (:)
        
		!k = 3
  !      ome = om
		!do_ = 1e-2_dp	!	0.01
        call oUpd (om)! + k * do_)
  !      
		nMax = 1d2	!	100
		!allocate (t_ (nMax))
		!call alf0 (t_, nMax)
		!
		!nMax = nP
  !      if (allocated (poles)) deallocate (poles)
		!allocate (poles (nMax))
		!poles = t_ (:nMax)
		!deallocate (t_)
		!
		!do i = k - 1, 0, -1
  !      om = ome
		!	call oUpd (om + i * do_)
		!	call alf1 (poles)
		!end do
        
        call HalfcHandler (poleF, 0._dp, 1 + 1.5 * om, eps, epsHc, nMax, poles, nP)
		
        t1 = minval (abs (poles)) / 2
        t4 = maxval (abs (poles)) + 1
        t2 = t1
        t3 = t1
        
	end subroutine pinit
	
	
	!===============================================		Обновление значения частоты		===============================================
	
	subroutine oUpd (om)
    
		real (dp)	:: om
		
		w = om
		omega = om
        
        kappa2 = EnuO * rhoO * omega * omega
		kappa = sqrt (kappa2)
        ka = a * kappa
		GalM = 1._dp / (kappa2 - alfM2)

	end subroutine oUpd
	
	
	!===============================================		Начальный этап трассировки полюсов		===============================================
	
	subroutine alf0 (t_, nPm)
	
		integer         :: i, j, k, nA, nPm
		complex (dp)    :: tc, t_ (nPm)
		
		nPm = 1
		nA = 1 + nint (1.5 * omega / eps)
		
        k = 1
        t_ = 0._dp
        do
			tc = k * epsC
			call crootw2 (tc, eC, epsHc, 1000, -7, poleF)
            
            k = k + 1
            if ((tc % re .ge. eps) .and. (abs (tc % im) .lt. eps)) exit
        end do
        t_ (1) % re = tc % re
        
		do1: do i = k, nA
			tc = i * epsC
			call crootw2 (tc, eC, epsHc, 1000, -7, poleF)
			
			if ((tc % re .lt. eps) .or. (abs (tc % im) .gt. eps)) cycle do1
			do j = 1, nPm
				if (abs (abs (t_ (j) % re) - abs (tc % re)) .le. 1e-2_dp) cycle do1
			end do
			
			nPm = nPm + 1
			t_ (nPm) % re = tc % re
		end do do1
		
		nP = nPm
		
	end subroutine alf0
	
	
	!===============================================		Трассировка полюсов		===============================================
	
	subroutine alf1 (alfa)
	
		integer     				:: i, j, nPm
		complex (dp)                :: tc, t_ (nP + 1)
		complex (dp), allocatable	:: alfa (:)
		
		nPm = nP + 1
		nP = 0
        
		t_ = [complex (dp) :: epsC, alfa]
		do2: do i = 1, nPm
			tc = t_ (i)
			call crootw2 (tc, eC, epsHc, 1000, -7, poleF)
			
			if ((tc .ne. tc) .or. (tc % re .lt. epsHc) .or. (abs (tc % im) .gt. eps)) cycle do2
			do j = 1, nP
				if (abs (abs (t_ (j) % re) - abs (tc % re)) .lt. epsHc) cycle do2
			end do
			
			nP = nP + 1
			t_ (nP) % re = tc % re
		end do do2
		
		if (nP .ne. nPm - 1) then
			deallocate (alfa)
			allocate (alfa (nP))
        end if
		alfa = t_ (:nP)
                
	end subroutine alf1
	
	
	!===============================================		Функция, нули которой отражают полюса		===============================================
	
	complex (dp) function poleF (alfa)
    
        complex (dp), intent (in)	:: alfa

		poleF = 1 / (Kalf (alfa) - b / (kappa2 - alfa * alfa))
		    
	end function poleF
	
	
	!===============================================		Выделение и очищение фрагмента рабочей памяти		===============================================
	
	subroutine alloc
	
		if (allocated (h)) then
			deallocate (h, E, nu, rhok, vPk, vSk)
        else
			allocate (h (nL), E (nL), nu (nL), rhok (nL), vPk (nL), vSk (nL))
		end if
	
	end subroutine alloc
	
	subroutine allocW (i)
    
		integer, optional	:: i
	
		if (allocated (cinf)) then
			deallocate (cinf, alfaM, alfM2, GalM)
        else
            mPsi = ceiling (real (nPhi) / 2)
            if (present (i)) mPsi = mPsi + i
			allocate (cinf (n0), alfaM (-mPsi:mPsi), alfM2 (-mPsi:mPsi), GalM (-mPsi:mPsi))
		end if
	
	end subroutine allocW
	
	end module PiezoMultiLayered
