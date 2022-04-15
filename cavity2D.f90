! gfortran -o post cavity2D.f90 
! gfortran -o post -Ofast cavity2D.f90
! ./post
! JESUS, Anna and GRION, Livia
! Last modification: Mar, 24

	program cavity2D
	
	implicit none
	
  	double precision, parameter :: Lx = 1.d0, Ly = 1.d0				! Comprimento do dominio na direcao x e y
  	integer, parameter :: Nx = 20, Ny = 20 							! Numero de celulas na direcao x e y
  	double precision, parameter :: dx = Lx/float(Nx), dy = Ly/float(Ny)	! Comprimento do volume finito em x e y
	double precision, parameter :: u0 = 1.d0                              ! Velocidade na tampa              
  	double precision, parameter :: vkinem = 1E-2 					! Viscosidade cinematica (mi/rho)   
	double precision, parameter :: rho = 1E0						! Massa especifica
  	double precision, parameter :: w = 1.94d0                             ! Omega do metodo SOR 
  	double precision dx2, dy2, beta								! auxiliares
   	integer, parameter :: nt = 1000								! Numero de passos de tempo
     integer, parameter :: NPRIN = 20								! Exibe a solucao                
 	
 	! Variaveis
	double precision, dimension(Nx+1,Ny) :: ustar, unew, un, deltau
	double precision, dimension(Nx,Ny+1) :: vstar, vnew, vn, deltav
	double precision, dimension(Nx,Ny) :: pnew, rhs, div_unew, erro
	double precision, dimension(Nx,Ny) :: CONVx, CONVy, VISCx, VISCy
	double precision dt, dt1, dt2, t, t1, t2, tfinal 
	double precision tol, divu, eps, dtn
	double precision errou, errov, errop, mu, mv 
	integer i, j, itera, it 
	character*100 arq

 	dx2 = dx*dx
	dy2 = dy*dy
	beta = dx2/dy2

 	itera = 0													! Contador SOR	
 	it = 0													! Contador evolucao temporal
	
	! Condicao inicial
	un = 0.d0 
	vn = 0.d0  
	
	unew = un
	vnew = vn
	pnew = 0.d0

	
	VISCx = 0.d0
	CONVx = 0.d0
	VISCy = 0.d0
	CONVy = 0.d0
	
	!Inicializacao
	errou = 1E0
	errov = 1E0
	erro = 0.d0 
	t = 0.d0	
	
	tfinal = 30.d0
	tol = 1E-8 
	dt = 1E-3 

		
	! Condicao para estacionaridade  norma(du/dt + dv/dt) < eps
	eps = 1E-8
	dtn = 1.d0

 	call cpu_time(t1)

	do while (dtn.ge.eps) ! Estacionario
	!do while (t.le.tfinal)  
	
	! Passo 1: Preditor
	! u^* = u^n + dt*( (mu/rho)*Lu^n - C(u^n) + (1/rho)*f^n)		
		do j =  2,Ny-1
    			do i = 2,Nx-1                 
    
     			VISCx(i+1,j) = (un(i+1,j+1) - 2.d0*un(i+1,j) + un(i+1,j-1))/dy2 			&
                			+ (un(i+2,j) - 2.d0*un(i+1,j) + un(i,j))/dx2  
                			
     			VISCy(i,j+1) = (vn(i,j+2) - 2.d0*vn(i,j+1) + vn(i,j))/dy2 + 			&
                			+ (vn(i+1,j+1) - 2.d0*vn(i,j+1) + vn(i-1,j+1))/dx2  
            
            
				CONVx(i+1,j) = (un(i+1,j+1) + un(i+1,j))*(vn(i+1,j+1) + vn(i,j+1))/dy/4.d0 &
              				- (un(i+1,j-1) + un(i+1,j))*(vn(i+1,j) + vn(i,j))/dy/4.d0  	 	& 
              				+ (un(i+2,j) + un(i+1,j))*(un(i+2,j) + un(i+1,j))/dx/4.d0  	 	& 
              				- (un(i,j) + un(i+1,j))*( un(i,j) + un(i+1,j))/dx/4.d0 
           
         
    				CONVy(i,j+1) = (vn(i,j+2) + vn(i,j+1))*(vn(i,j+2) + vn(i,j+1))/dy/4.d0  	& 
              				- (vn(i,j) + vn(i,j+1))*(vn(i,j) + vn(i,j+1))/dy/4.d0 		  	&
              				+ (un(i+1,j+1) + un(i+1,j))*(vn(i+1,j+1) + vn(i,j+1))/dx/4.d0 	&
              				- (un(i,j) + un(i,j+1))*(vn(i,j+1) + vn(i-1,j+1))/dx/4.d0 

    
				ustar(i+1,j) = un(i+1,j) + dt*( -CONVx(i+1,j) + vkinem*VISCx(i+1,j) )				
				vstar(i,j+1) = vn(i,j+1) + dt*( -CONVy(i,j+1) + vkinem*VISCy(i,j+1) )     
          
          
			end  do    
		end do

 		ustar(2,2:Ny-1) = 0.d0
		ustar(Nx,2:Ny-1) = 0.d0
     
		vstar(2:Nx-1,2) = 0.d0
		vstar(2:Nx-1,Ny) = 0.d0
 
 		! Passo 2 : Lado direito equação de Poisson para a pressao
     	! D G(p^{n+1}) = (rho/dt)*D(u^*) 
		do j = 2,Ny-1 
			do i = 2,Nx-1
				rhs(i,j) = rho/dt*((ustar(i+1,j) - ustar(i,j))/dx + (vstar(i,j+1) - vstar(i,j))/dy) 
			end do
		end do

		errop = 1E0 
		itera = 0 
 
		! Passo 3: Resolve o sistema linear (Equação de Poisson)
		! D G(p^{n+1}) = (rho/dt)*D(u^*)
		do while (errop.ge.tol)
  			
			pnew(1,2:Ny-1) = pnew(2,2:Ny-1) 
			pnew(Nx,2:Ny-1) = pnew(Nx-1,2:Ny-1) 
			pnew(2:Nx-1,1) = pnew(2:Nx-1,2) 
			pnew(2:Nx-1,Ny) = pnew(2:Nx-1,Ny-1) 
			
  			call sor_method(pnew, rhs, Nx, Ny, dx2, beta, w)  		! Metodo SOR 
  		     call calc_erro(pnew, rhs, Nx, Ny, dx2, dy2, erro)				! Residuo
 			errop = norm2(erro)									! norma do residuo
  			itera = itera + 1  	
			pnew = pnew  -  pnew(2,2);					! Normalização da pressao
  			
  			!write(*,*) itera, errop
  		  			   
		end do
		  			
		! Passo 4: Correção da velocidade
		! u^{n+1} = u^* - (dt/rho)G p^{n+1}
		do j = 2,Ny-1
    			do i = 2,Nx-1
        
        			unew(i+1,j) = ustar(i+1,j) - dt/rho*(pnew(i+1,j) - pnew(i,j))/dx      
        			vnew(i,j+1) = vstar(i,j+1) - dt/rho*(pnew(i,j+1) - pnew(i,j))/dy     
        
       			! Calcula a incompressibilidade D u^{n+1} = 0
        			div_unew(i,j) = ((unew(i+1,j) - unew(i,j))/dx + (vnew(i,j+1) - vnew(i,j))/dy)  
        		  
			end do       
		end do
     
         	! Passo 5: Atualiza o campo de velocidade no contorno em t_{n+1}
		unew(3:Nx-1,1) = - unew(3:Nx-1,2)  
		unew(3:Nx-1,Ny) = 2.d0*u0 - unew(3:Nx-1,Ny-1)   
		vnew(1,3:Ny-1) = - vnew(2,3:Ny-1) 
		vnew(Nx,3:Ny-1) = - vnew(Nx-1,3:Ny-1)  
    
     	! Encontra o maximo da velocidade, em módulo
    		!mu = maxval(abs(unew))
    		!mv = maxval(abs(vnew))
    		
    		mu = 2.d0
    		mv = 1.d0
    		
		dt1 = min(dx/mu, dy/mv)									! Condição CFL
		dt2 = 0.5d0/vkinem/(1.d0/dx2 + 1.d0/dy2) 
		dt = 0.4d0*min(dt1,dt2) 									! Constante 0.2< tal< 0.6 => tal = 0.4
         
          ! Delta da velocidade u^{n+1} - u^n
          deltau = un - unew;
         	deltav = vn - vnew;
		
		errou = norm2(deltau) 
		errov = norm2(deltav) 
		divu = norm2(div_unew)

		! dtn ~ du/dt + dv/dt
		dtn =  (errou + errov)/dt;
		
		! if (mod(it,NPRIN).eq.0) then
		!  	print  '("t=",(ES11.4E2), " |u|=",(ES11.4E2), " |v|=",(ES11.4E2), &
		! 	 	" |p|=",(ES11.4E2),  " |divu|=",(ES11.4E2),  " |e|=", (ES11.4E2))', &
		!   	 		 t, errou, errov, errop, divu, dtn 
	     ! 	write(*,*) itera, errop	 		 
          !	end if  
         	 
         	! Atualiza o campo de velocidade  
		un = unew;
		vn = vnew;       
      
      	! Atualiza o instante de tempo   
     	it = it + 1;        
     	t = t + dt;          
            
	        
	end do !endwhile

	call cpu_time(t2)

	print  '("CPU TIME FORTRAN  = ", (ES11.4E2))', t2 - t1
	print ' ("CPU TIME FORTRAN  = ", (100f18.12))', t2 - t1
	print  '("t=",(ES11.4E2), " |u|=",(ES11.4E2), " |v|=",(ES11.4E2), &
		   	 	" |p|=",(ES11.4E2),  " |divu|=",(ES11.4E2),  " |e|=", (ES11.4E2))', &
		   	 		 t, errou, errov, errop, divu, dtn 
   	write(*,*) itera, errop
    ! U =  0.5d0*(un(2:Nx+1,1:Ny) +  un(1:Nx,1:Ny))
    ! V =  0.5d0*(vn(1:Nx,2:Ny+1) +  vn(1:Nx,1:Ny))
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! pos processamento 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(12,file='U20.dat')
		do j = 1, Ny 
			do i = 1, Nx+1 
			write(12,'(100f18.12)') unew(i,j)
		end do
		write(12,*)
	end do 
	close(12)
	
	open(12,file='V20.dat')
		do j = 1, Ny+1 
			do i = 1, Nx  
			write(12,'(100f18.12)') vnew(i,j) 
		end do
		write(12,*)
	end do 
	close(12)
	
	open(12,file='P20.dat')
		do j = 1, Ny 
			do i = 1, Nx 
			write(12,'(100f18.12)') pnew(i,j)
		end do
		write(12,*)
	end do 
	close(12)

   
end program  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine sor_method(f, tf, imax, jmax, dx2, beta, w)

	implicit none
	integer i, j, imax, jmax
  	double precision dx2, dy2, beta, w
  	double precision, dimension(imax, jmax) :: f,tf

		do j = 2, jmax-1
  			do i = 2, imax-1
				f(i,j) = (1.d0 - w) * f(i,j) + w * (-dx2*tf(i,j) + f(i-1,j) + f(i+1,j) + &
             				beta*(f(i,j-1)+f(i,j+1))) / (2.d0*(beta + 1.d0))
  			end do
		end do
  
		return
	
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine calc_erro(f, tf, imax, jmax, dx2, dy2, erro)

	implicit none
	integer i, j, imax, jmax
  	double precision dx2, dy2
  	double precision, dimension(imax, jmax) :: f, tf, erro

		erro = 0.d0
		do j = 2, jmax-1
  			do i = 2, imax-1
    				erro(i,j)  =  abs(tf(i,j) - (f(i-1,j) - 2.d0*f(i,j) + f(i+1,j))/dx2  - & 
    								(f(i,j-1) -2.d0*f(i,j) + f(i,j+1))/dy2) 
  			end do
		end do

		
		return  
	end
