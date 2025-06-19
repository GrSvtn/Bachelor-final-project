    program labtask8
	
	    call start
	
	end program labtask8

	
    subroutine start
    
        logical             :: key
        integer             :: ch
        character (len = 2) :: inpt
        
        call wr

        read (*, '(a)') inpt
        ch = 4
        if (trim (inpt) .ne. '') read (inpt, *) ch
        
        print *, 'Do you want nPhi to be choosen recommendedly? Enter - yes, any key - no'
        
        read (*, '(a)') inpt
        key = trim (inpt) .eq. ''
        
        call coefficients (ch, key)
        
    end subroutine start
    
    
    subroutine wr
	
1		format (x, a)

		print '(a)', 'Choose subtask (4 - by default):'
        print 1, '(1) Dispersion curves on [0, 10]'
        print 1, '(2) Coefficients for q with even nPhi >= n'
        print 1, '(3) Coefficients for q with nPhi = n, + 22, 2'
        print 1, '(4) Function v on [-a, a]'
    
	end subroutine wr
	
	
	subroutine coefficients (ch, key)
	use GlobalsPiezoMultilayered
	use PiezoMultiLayered, only: dispcur, qCoeffs, alloc, allocW, vFunc, uv
    use functions, only: int2str
	implicit none
	
        logical                     :: key
	    integer			    		:: i, u, nx, ch
        real (dp)				    :: acopy
        real (dp), allocatable  	:: x (:)
        complex (dp), allocatable   :: v (:)
        
        open (newunit = u, file = 'data/input.dat', action = 'read')	!	Открываю входной файл
		read (u, *) a, hO, EO, nuO, rhoO, d31, Ve						!	Считываю основные параметры
        read (u, *) nL, nPhi, omega
		
        if (key) nPhi = nint (a * omega) + 11   !   При v_S = 1
        nPhi = nPhi + mod (nPhi, 2)	!	Необходимость квадратности системы (n = 2 m)
        
		call alloc	!	Выделение памяти на глобальные параметры
		do i = 1, nL
			read (u, *) h (i), E (i), nu (i), rhok (i)
		end do
		
		close (u)
        
        nx = 100	!	Сколько точек придётся на (0, a]
        allocate (x (2 * nx + 1))
        
        call defx (x, a, nx)
        
        print *, 'nPhi = ', nPhi
        select case (ch)
            case (1)
            
                call dispcur
            
            case (2)
            
                call qCoeffs
            
    			call wrfc ('q', qc)	!	Записывает вычисленные значения в файлы
				deallocate (qc)
            
            case (3)
                
                u = nPhi
                do i = u, u + 22, 2
                    nPhi = i
					call allocW
					call qCoeffs
            
    				call wrfc ('qPhi/q' // int2str (i), qc)
					print *, 'nPhi = ', int2str (i)
                    
					call allocW
				end do
                deallocate (qc)
                
            case (4)
                
				call allocW
				call qCoeffs
                call vFunc (x, v)
            
				open (newunit = u, file = 'data/v/vf.dat')
    			write (u, '(f12.6, <2 * nx> (", ", f12.6))') abs (v)
                close (u)
                
                acopy = real (nint (a + 0.9))
                call defx (x, acopy, nx)
                call uv (x, v)
                
                open (newunit = u, file = 'data/u/uv.dat')
    			write (u, '(f12.6, <2 * nx> (", ", f12.6))') abs (v)
                close (u)
                
                deallocate (v)
				call allocW
			
        end select
		
        deallocate (x)
        call alloc
            
		
		contains
	
		subroutine wrfc (name, A)	!	Записывает комплексные значения в два файла
		implicit none
	
			integer			:: n, u1, u2
			character (*)	:: name
			complex (dp)	:: A (:)
		
			open (newunit = u1, file = 'data/qvCoefs/' // name // 'Re.dat', action = 'write')
			open (newunit = u2, file = 'data/qvCoefs/' // name // 'Im.dat', action = 'write')
		
			n = size (A)
			write (u1, '(f12.6, <n - 1> (", ", f12.6))') A % re
			write (u2, '(f12.6, <n - 1> (", ", f12.6))') A % im
		
			close (u1)
			close (u2)
	
		end subroutine wrfc
	
    end subroutine coefficients
    
    subroutine defx (x, leng, nx)
	use GlobalsPiezoMultilayered
    
        integer     :: i, nx
        real (dp)   :: dx, leng, x (2 * nx + 1)
        
        dx = leng / nx
        x = [(i * dx - leng, i = 0, 2 * nx)]
        
    end subroutine defx
