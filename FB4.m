function X=FB4(kappa,gamm,n)
% FB4       Generates a random sample of size n from the Fisher-Bingham_4 distribution
%           as a subfunction that is called by Random_FB4,Random_FB5,
%           Random_FB6
%
% W2 = FB4(kappa,gamm,n);
%
% kappa         distribution parameter (-inf < kappa < inf <)
% gamm          distribution parameter (-inf < gamma < inf)
% n             sample size ( where N is a positive integer)
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
% The methods implemented here are based on
% "The simulation of spherical distributions in the Fisher-Bingham family"
% by Andrew T.A. Wood (1987)
% Communications in Statistics - Simulation and Computation, 16:3, 885-898,
% DOI: 10.1080/03610918708812624
%% 
if kappa==0 && gamm==0;
    U = rand(n,1);
    X = 2*U-1; % cos(theta)
    return
elseif kappa==0;
    X= Watson_LW(gamm,n); % Watson !!!!!!!! which one??
    return
else gamm==0;
    X=vMF_Wood(kappa,n,3);% when gamma==0 FB4 is the von Mises-Fisher distribution
    return
end

%% 

cFB4= @(kappa4,gamm1)(integral(@(u)(exp(kappa4*u+gamm1*u.^2)),-1,1)); % normalizing constant
% ckg=cFB4(kappa-gamm,0);
%% 
if gamm < 0 % if gamma<0, then FB4 is a truncated normal
    if kappa <= -2*gamm
        X=[];
        while length(X)< n
            a0= normcdf(1,-kappa/gamm/2,sqrt(-1/gamm/2))-...
                normcdf(-1,-kappa/gamm/2,sqrt(-1/gamm/2)); % acceptance ratio
            X=normrnd(-kappa/gamm/2,sqrt(-1/gamm/2),[ceil(n/a0),1]);
            X=X(X>=-1 & X<=1);
        end
    elseif kappa > -2*gamm
        X=[];  
        a2 = exp(gamm)*cFB4(kappa,gamm)/cFB4(kappa+2*gamm,0); % acceptance ratio
        while length(X)< n
            x2 = vMF_Wood(kappa+2*gamm,3,ceil(n/a2));
            x_U= rand(ceil(n/a2),1);
            x2 = x2(x_U<=exp(gamm)*exp(-2*gamm*x2+gamm*x2.^2));
            X  = [X;x2 ] ;
        end
    end
elseif gamm > 0 % corresponds to a positive quadratic term
    p1 = cFB4(kappa+gamm,0)/(cFB4(kappa+gamm,0)+cFB4(kappa-gamm,0)); % mixing proportion
    a3= p1*(1+exp(-2*gamm))*cFB4(kappa,gamm)/cFB4(kappa+gamm,0); % acceptance ratio
    X=[];
    while length(X) < n
        x3_minus  = vMF_Wood(kappa-gamm,3,ceil(n/a3)); % marginal vMF distribution
        x3_plus  = vMF_Wood(kappa+gamm,3,ceil(n/a3)); % marginal vMF distribution
        x3_U = rand(length(x3_plus),1);
        x3_mix = x3_plus .* (x3_U <= p1)+x3_minus .* (1-(x3_U <= p1));
        %mix corresponds to the mixed density 
        gV=((1+exp(-2*gamm))*exp(gamm*(x3_mix.^2)))./(exp(gamm*x3_mix)+exp(-gamm*x3_mix));
        U=rand(ceil(n/a3),1);
        x3_mix = x3_mix(U<=gV);
        X = [X;x3_mix];
    end
end
X=X(1:n);
end
