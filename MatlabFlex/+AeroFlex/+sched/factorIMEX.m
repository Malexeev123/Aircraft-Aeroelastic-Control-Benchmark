function fac = factorIMEX(L, dt)
%FACTORIMEX LU factorization used by ROMIntegrator and SimRunner.
gamma = (2 - sqrt(2))/2;
A = speye(size(L)) - gamma*dt*L;
[Lfac,Ufac,piv] = lu(A,'vector');
fac = struct('gamma',gamma,'delta',-2*sqrt(2)/3, ...
             'Lfac',Lfac,'Ufac',Ufac,'piv',piv);
end
