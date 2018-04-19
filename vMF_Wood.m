function X=vMF_Wood(kappa,n,d)
% VMFB Generates a random sample of size N from the von Mises-Fisher distribution
%
% W1=vMF_Wood(kappa,d,n);
%
% kappa         distribution parameter, concentration parameter 
% N             sample size ( where N is a positive integer)
% d             dimension of the space, for a 2D sphere d=3
% Authors: Gy.Terdik, B.Wainwright
%
% Copyright 2018 Gy.Terdik
%
% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this list 
%    of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice, this list 
%    of conditions and the following disclaimer in the documentation and/or other materials 
%    provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its contributors may be used 
%    to endorse or promote products derived from this software without specific prior written 
%    permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS 
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY 
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
% IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
% OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%
%
% The methods implemented here are based on
% "Simulation of the von mises fisher distribution"
% by Andrew T.A. Wood (1994)
% Communications in Statistics - Simulation and Computation, 23:1, 157-164,
% DOI: 10.1080/03610919408813161
%% Wood 

% d=3;  dim of the space / m=3 corresponds to the unit sphere

    kappa1=kappa;
    kappa=abs(kappa);
    %%%%%%%%%%%%%%%%
    %%%  Step 0  %%%
    %%%%%%%%%%%%%%%%
    b  = (-2*kappa +sqrt(4*(kappa^2)+(d-1)^2) )/(d-1);
    x0 = (1-b)/(1+b);
    c  = kappa*x0 + (d-1)*log(1-x0^2);
    A  = (d-1)/2; %shape parameter for Beta distribution
    % n  = min(2^12,2*n);

    X=[];
    while(length(X) < n)
        %%%%%%%%%%%%%%%%
        %%%  Step 1  %%%
        %%%%%%%%%%%%%%%%
        U = rand(n,1); %random n-by-1 vector from uniform(0,1)
        Z = betarnd(A,A,n,1); %random n-by-1 vector from Beta(A,A)
        W = (1-(1+b)*Z)./(1-(1-b)*Z);
       %%%%%%%%%%%%%%%%
        %%Steps 2 & 3%%%
        %%%%%%%%%%%%%%%%
        Step2 = kappa*W+(d-1)*log(1-x0*W)-c;
        X=[X;W( Step2 > log(U) )]; %Good Indices
    end
    X=sign(kappa1)*X(1:n);
 
