module am_toy_models
    !
    use am_helpers
    use am_constants
    !
    implicit none
    !
    
 
 
!clear;clc;
! 
!global hbar 
!global m 
!hbar = 1;
!m = 1;
! 
!% input parameters
!nplanewaves = 13 ; % number of fourier components (should be an odd number)
!nkpt = 100; % number of kpoints
! 
!% Describe how to setup the Hamiltonian
!Mkl = @(h,l) 0*(abs(h-l)==1) + 0.05*(abs(h-l)==2) + 0.1*(abs(h-l)==3);
!Ehk = @(h,k) hbar^2/(2*m)*norm(k+h)^2;
!Hamiltonian = @(h,l,k) (h==l)*(Ehk(h,k)) + (h~=l)*Mkl(h,l);
! 
!% if mod(nfourier,2) == 0; nfourier = nfourier+1; end 
!planewave = ceil(-nplanewaves/2):1:floor(nplanewaves/2);
!kpt = linspace(-1/2,1/2,nkpt);
!for k = 1:nkpt
! 
!    % construct the Hamiltonian
!    for j = 1:nplanewaves
!    for i = 1:nplanewaves
!        H(i,j) = Hamiltonian(planewave(i),planewave(j),kpt(k));
!    end 
!    end
!    
!    % Solve for electron energies (eigenvalues)
!    d = eig(H);
!    
!    % Save result : E(band,kpoint)
!    E(:,k) = d(:);
!end
! 
!% shift energy scale to put lowest band at zero.
!E = E-min(min(E));
!kpt=mod(kpt-1/2,1)+1/2;
!plot(kpt,E,'.')
!ylim([0 1])
! 
!%%
!% plot density of states
!Edos = reshape(E,1,[]);
!Edos = sort(Edos);
!dE = 0.005;
!Ei = [0:dE:max(Edos)];
!plot(Ei,hist(Edos,Ei)/dE/nkpt)
!axis([0 1 0 30])
 
 
 


    
    
end module
    