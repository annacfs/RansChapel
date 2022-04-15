//chpl  cavity2D.chpl 
// chpl --fast cavity2D.chpl   
// ./cavity
// JESUS, Anna and GRION, Livia
// Last modification: Mar, 25
 
 	use LinearAlgebra;
	use Time;
	var runtime: Timer;
	
	param Lx = 1.0, Ly = 1.0 : real;								// Comprimento do dominio na direcao x e y
	param Nx = 20, Ny = 20 : int;									// Numero de celulas na direcao x e y
	param dx = Lx/Nx, dy = Ly/Ny : real; 							// Comprimento do volume finito em x e y
	param u0 = 1.0 : real;                 							// Velocidade na tampa
	param vkinem = 1E-2 : real;             						// Viscosidade cinematica (mi/rho)    
	param rho = 1.0 : real;										// Massa especifica
	param w  = 1.940 : real;               							// Omega do metodo SOR
	param p : normType;											// defaut: norma euclidiana 
	param dx2, dy2, beta : real;									// auxiliares
	param nt = 1000: int;										// Numero de passos de tempo
 	param NPRIN = 20: int;										// Exibe a solucao 
 
 	// Variaveis
	var ustar, unew, un, deltau : [1..Nx+1,1..Ny] real;  
	var vstar, vnew, vn, deltav : [1..Nx,1..Ny+1] real;  
	var pnew, rhs, div_unew, erro : [1..Nx,1..Ny] real;  
 	var CONVx, VISCx, CONVy, VISCy : [1..Nx,1..Ny] real;   
 	
	var dt, dt1, dt2, t, t1, t2, tfinal: real;
	var tol, divu, eps, dtn: real;
	var errou, errov, errop, mu, mv: real;  
	var i, j, itera, it : int;
	
	dx2 = dx*dx;
	dy2 = dy*dy;
	beta = dx2/dy2;
	itera = 0;								 				// Contador SOR
	it = 0;									 				// Contador da evolucao temporal

	// Condicao inicial
	un = 0.0;
	vn = 0.0;

	unew = un;
	vnew = vn;
	pnew = 0.0; 
	
	VISCx = 0.0;
	CONVx = 0.0;
	VISCy = 0.0;
	CONVy = 0.0;
	
	// Inicializacao 
	errou = 1E0;
	errov = 1E0;
	erro = 0.0;
	t = 0.0;
	
	tfinal = 30.0;
	tol = 1E-8;
	dt = 1E-3; 

	
	// Condicao para estacionaridade  norma(du/dt + dv/dt) < eps
	eps = 1E-8; 								
	dtn = 1.0;

	runtime.start();

	serial{
	//while(t <= tfinal) do{
	while(dtn >= eps) do{ //Estacionario
	
		// Passo 1: Preditor
		// u^* = u^n + dt*( (mu/rho)*Lu^n - C(u^n) + (1/rho)*f^n)
		for i in 2..Nx-1 do{
			for j in 2..Ny-1 do{
                
				VISCx[i+1,j] = (un[i+1,j+1] - 2.0*un[i+1,j] + un[i+1,j-1])/dy2  
                			+ (un[i+2,j] - 2.0*un[i+1,j] + un[i,j])/dx2;
     
				VISCy[i,j+1] = (vn[i,j+2] - 2.0*vn[i,j+1] + vn[i,j])/dy2 +  
                			+ (vn[i+1,j+1] - 2.0*vn[i,j+1] + vn[i-1,j+1])/dx2; 
            
            
				CONVx[i+1,j] = (un[i+1,j+1] + un[i+1,j])*(vn[i+1,j+1] + vn[i,j+1])/dy/4.0   
              				- (un[i+1,j-1] + un[i+1,j])*(vn[i+1,j] + vn[i,j])/dy/4.0   
              				+ (un[i+2,j] + un[i+1,j])*(un[i+2,j] + un[i+1,j])/dx/4.0  
              				- (un[i,j] + un[i+1,j])*( un[i,j] + un[i+1,j])/dx/4.0;
           
         
				CONVy[i,j+1] = (vn[i,j+2] + vn[i,j+1])*(vn[i,j+2] + vn[i,j+1])/dy/4.0  
              				- (vn[i,j] + vn[i,j+1])*(vn[i,j] + vn[i,j+1])/dy/4.0  
              				+ (un[i+1,j+1] + un[i+1,j])*(vn[i+1,j+1] + vn[i,j+1])/dx/4.0 
              				- (un[i,j] + un[i,j+1])*(vn[i,j+1] + vn[i-1,j+1])/dx/4.0;

				ustar[i+1,j] = un[i+1,j] + dt*( -CONVx[i+1,j] + vkinem*VISCx[i+1,j] );	
				vstar[i,j+1] = vn[i,j+1] + dt*( -CONVy[i,j+1] + vkinem*VISCy[i,j+1] );     
          
			}
		}

		ustar[2,2..Ny-1] = 0.0;
		ustar[Nx,2..Ny-1] = 0.0;
     
		vstar[2..Nx-1,2] = 0.0;
		vstar[2..Nx-1,Ny] = 0.0;
     
     	// Passo 2 : Lado direito equação de Poisson para a pressao
     	// D G(p^{n+1}) = (rho/dt)*D(u^*)     	
		for i in 2..Nx-1 do{
			for j in 2..Ny-1 do{	
				rhs[i,j] = rho/dt*((ustar[i+1,j] - ustar[i,j])/dx + (vstar[i,j+1] - vstar[i,j])/dy);
			}
		}

		errop = 1E0;
		itera = 0;

		// Passo 3: Resolve o sistema linear (Equação de Poisson)
		// D G(p^{n+1}) = (rho/dt)*D(u^*)
		while (errop  >= tol) do{    
     
		 	pnew[1,2..Ny-1] = pnew[2,2..Ny-1];
		 	pnew[2..Nx-1,1] = pnew[2..Nx-1,2];
			pnew[Nx,2..Ny-1] = pnew[Nx-1,2..Ny-1];
			pnew[2..Nx-1,Ny] = pnew[2..Nx-1,Ny-1];
         
			pnew  = sor_method(pnew, rhs, Nx, Ny, dx2, beta, w);		// metodo SOR 			 		 
			erro  = calc_erro(pnew, rhs, Nx, Ny, dx2, dy2);  			// residuo
			errop = norm(erro, p); 								// norma do residuo
	   		itera = itera + 1;		   	
		   	pnew = pnew -  pnew[2,2];					// normalização da pressao
			
			//writef('%i %6.12dr\n', itera, errop); 
			//writef("%i   %12.8er \n",itera, errop);
   
		}
 		

 
 		// Passo 4: Correção da velocidade
 		// u^{n+1} = u^* - (dt/rho)G p^{n+1}
		for i in 2..Nx-1 do{
			for j in 2..Ny-1 do{
        
				unew[i+1,j] = ustar[i+1,j] - dt/rho*(pnew[i+1,j] - pnew[i,j])/dx;      
				vnew[i,j+1] = vstar[i,j+1] - dt/rho*(pnew[i,j+1] - pnew[i,j])/dy;     
        
        			// Calcula a incompressibilidade D u^{n+1} = 0
				div_unew[i,j] = ((unew[i+1,j] - unew[i,j])/dx + (vnew[i,j+1] - vnew[i,j])/dy);    

   			}
   		}
   
		// Passo 5: Atualiza o campo de velocidade no contorno em t_{n+1}
		unew[3..Nx-1,1] = - unew[3..Nx-1,2];
		unew[3..Nx-1,Ny] = 2.0*u0 - unew[3..Nx-1,Ny-1];
		vnew[1,3..Ny-1] = - vnew[2,3..Ny-1];
		vnew[Nx,3..Ny-1] = - vnew[Nx-1,3..Ny-1];  
    
    		// Encontra o maximo da velocidade, em módulo	
		//mu = max reduce(abs(unew));
		//mv = max reduce(abs(vnew)); 
    		mu = 2.0;
    		mv = 1.0;
    		 
 		dt1 = min(dx/mu, dy/mv); 								// Condição CFL
		dt2 = 0.5/vkinem/(1/dx2 + 1/dy2);						
		dt = 0.4*min(dt1,dt2);         							// Constante 0.2 < tal< 0.6 => tal = 0.4
        
        	// Delta da velocidade u^{n+1} - u^n
         	deltau = un - unew;
         	deltav = vn - vnew;
         
		errou = norm(deltau,p);
		errov = norm(deltav,p);
		divu =  norm(div_unew,p);


		dtn =  (errou + errov)/dt;								// dtn ~ du/dt + dv/dt			 
              
		/*	if (mod(it,NPRIN) == 0) then {
					writef("t= %5.4Er |u|= %5.4Er |v|= %5.4Er |p|= %5.4Er |divu|= %5.4Er |e|= %5.4Er\n", 
							t, errou, errov, errop, divu, dtn);  
					writef("%i   %12.8er \n",itera, errop);           

		} */
		
		// Atualiza o campo de velocidade           
          un = unew;
		vn = vnew;  
	 	it = it + 1;        
	 	t = t + dt;          
  

	} // endwhile
}
	runtime.stop();
	writef(" CPU TIME CHAPEL  = %5.4Er \n", runtime.elapsed());
	writef(" CPU TIME CHAPEL  = %10.16dr \n",runtime.elapsed());
	writef("--------------------------------------------------------------------------------------------- \n"); 
	writef("t= %5.4Er |u|= %5.4Er |v|= %5.4Er |p|= %5.4Er |divu|= %5.4Er |e|= %5.4Er\n", 
							t, errou, errov, errop, divu, dtn);   
	writef("%i   %12.8er \n",itera, errop);						                 
	writef("--------------------------------------------------------------------------------------------- \n"); 

	//U =  0.5*(un(2..Nx+1,1..Ny) +  un(1..Nx,..Ny))
	//V =  0.5*(vn(1..Nx,2..Ny+1) +  vn(1..Nx,1..Ny))
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// pos processamento 
// PS: o loop está invertido para facilitar a leitura no pos processamento
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
use IO;                            
const u = openwriter("U20.out");       
		for j in 1..Ny  do{  
			for i in 1..Nx+1  do{             
        			u.writef(" %10.12dr \n", unew[i,j]);
        		}
        u.writef("\n ");
		}
u.close();

const v = openwriter("V20.out"); 
		for j in 1..Ny+1  do{ 
			for i in 1..Nx  do{              
        			v.writef(" %10.12dr \n", vnew[i,j]);
        }
              v.writef("\n ");
	}
v.close();

const pp = openwriter("P20.out");      
		for j in 1..Ny  do{  
			for i in 1..Nx  do{             
        		pp.writef(" %10.12dr \n", pnew[i,j]);
        }
        		pp.writef("\n ");
	}
pp.close();
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	proc sor_method(f, tf, imax: int, jmax: int,  dx2: real, beta: real,  w: real){ 

		for i in 2..imax-1  do{
     		for j in 2..jmax-1  do{
          
				f[i,j] = (1  - w) * f[i,j] + w * (-dx2*tf[i,j] + f[i-1,j] + f[i+1,j] +  
						beta*(f[i,j-1]+f[i,j+1])) / (2.0*(beta + 1.0)) ; 
	 
		}
	}

		return f;     
       
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	proc calc_erro(f, tf, imax: int, jmax: int, dx2: real, dy2: real){

	var errol : [1..Nx,1..Ny] real; 
	errol = 0.0;

	for i in 2..imax-1  do{
		for j in 2..jmax-1 do{
		 
      		errol[i,j] =   abs(tf[i,j] -(f[i-1,j] - 2.0*f[i,j] + f[i+1,j])/dx2 - 
              			(f[i,j-1] - 2.0*f[i,j] + f[i,j+1])/dy2 );      
 		}
  	}
     
		return errol;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
