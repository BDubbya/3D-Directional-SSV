function X = Watson_Fish(gamm,n)
% WATSON_FISH   Generates a random sample of size n from the 
%               Watson-Fisher distribution (in the FB4 family)
%
% X = Watson_Fish(gamm,n)
%
% gamm          distribution parameter (-inf < gamma < inf)
% n             sample size ( where n is a positive integer)
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

%% density

X=[];

if gamm > 0; % bipolar
     C=1/(exp(gamm)-1);
    while length(X)< n
        U1 = rand(n,1); 
        U2 = rand(n,1);
        X0 = log(U1/C+1)/gamm;
        X1 = X0(U2 <= exp(gamm*(X0.^2-X0)));
        V = rand(length(X1),1);
        X1 = X1.*((V<1/2)-(V>=1/2));
        X = [X; X1]; %  cos(theta)
    end
else
     C1 = sqrt(-gamm);
     C2 = atan(C1);
    while length(X)< n
        U1 = rand(n,1);
        U2 = rand(n,1);
        X0 = tan(C2*U1)/C1;
        X1 = X0(U2 <= (1-gamm*X0.^2).*exp(gamm*X0.^2));
        V = rand(length(X1),1);
        X1 = X1.*((V<1/2)-(V>=1/2));
        X = [X; X1]; %  cos(theta)
    end
end
X = X(1:n);
