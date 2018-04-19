function X = Watson_LW(gamm,n)
% WATSON_LW Generates a random sample of size N from the Watson 
%           distribution (FB4 family)
%
% Y = Watson_LW(gamm,N)
%
% gamm          distribution parameter (-inf < gamma < inf)
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
% "The simulation of spherical distributions in the Fisher-Bingham family"
% by Andrew T.A. Wood (1987)
% Communications in Statistics - Simulation and Computation, 16:3, 885-898,
% DOI: 10.1080/03610918708812624


X=[];

%% Li Wong
if gamm > 0; % bipolar
    ro=4*gamm/(2*gamm+3+sqrt((2*gamm+3)^2-16*gamm));
    r=(3*ro/2/gamm)^3*exp(-3+2*gamm/ro);
    while length(X)< n
        U1= rand(n,1);  %ceil(N/a2)
        U2= rand(n,1);
        S= U1.^2./(1-ro*(1-U1.^2)); % 
        W=gamm*S;
        V=r*U2.^2./((1-ro*S).^3);
        X1 = sqrt(  S(V <= exp(2*W))); %
        V=rand(length(X1),1);
        X1 =X1.*((V<1/2)-(V>=1/2));
        X=[X; X1]; % cos(theta)
    end
else  % gridle
    b=exp(2*gamm)-1;
    
    while length(X)< n
        U1= rand(n,1);  %ceil(N/a2)
        U2= rand(n,1);
        V=log(1+U1*b)/gamm;
        ksi=2*pi*U2;
        c=cos(ksi);
        S1=V.*c.^2;         S2=V-S1;
        Sind=logical((S1<=1).*(S2<=1));
        S1=S1(Sind);          S2=S2(Sind);%[S1 S2]
        d=sqrt(V);
        X1=d.*c; X2=d.*sin(ksi);
        Z=[X1(Sind)';X2(Sind)'];
        
        X=[X; Z(:)];  % cos(theta)
    end
end
X=X(1:n);

