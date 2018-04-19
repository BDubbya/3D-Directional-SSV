function X=vMF_Inv(kappa,n)
% VMFB Generates a random sample of size N from the von Mises-Fisher distribution
%
% W1=vMFb(kappa,N)
%
% kappa         distribution parameter, concentration parameter (kappa >= 0)
% N             sample size ( where N is a positive integer)
%
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
% dim of the space / m=3 corresponds to the unit sphere

%% Rubinstein 81, p.39, Fisher 87, p.59

     kappaS=sign(kappa);
     kappa=abs(kappa);
%      kappa1=exp(-2*kappa);
     U = rand(n,1); %random n-by-1 vector from uniform(0,1)
     X=log(2*U*sinh(kappa)+exp(-kappa))/kappa;
%      cos(2*asin(sqrt(-log(U*(1-kappa1)+kappa1)/(2*kappa)))); %

 X=kappaS*X(1:n);
end
