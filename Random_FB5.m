function Y = Random_FB5(kappa,bet,varargin)
% RANDOM_FB5 Generates a random sample of size n from the Fisher-Bingham_5 
%            distribution
%
% Y = Random_FB5(kappa,bet,Mu,n)
%%
% Examples of correct usage:
%
% n=2^12; 
% kappa = 4.2;
% bet = 4.5 ; 
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
%
% Examples of correct function construction: 
% Y = Random_FB5(kappa,bet,Mu,Psi,n)
% Y = Random_FB5(kappa,bet,Mu,Psi)
% Y = Random_FB5(kappa,bet,Mu)
% Y = Random_FB5(kappa,bet)
%%
%
% Required Inputs
% kappa       distribution parameter (kappa >= 0)
% bet         distribution parameter (beta >= 0)
%
% Optional Inputs
% Mu          distribution parameter (for rotating the North pole to Mu)
% n           sample size ( where n is a positive integer)
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
%% Verify/Default parameter values

% Required parameter characteristics
p = inputParser;
dbl = {'double'};
noneg = {'nonnegative'};
class = {'numeric'};
pos_int = {'integer','positive'};
real = {'real'};
noinf = {'finite'};

% Check characteristics for the required and optional parameters
addRequired(p,'kappa',@(x)validateattributes(x,dbl,real));
addRequired(p,'bet',@(x)validateattributes(x,dbl,noneg));
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'Psi',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'n',[],@(x)validateattributes(x,class,pos_int));

p.parse(kappa,bet,varargin{:});
Mu = p.Results.Mu;
Psi = p.Results.Psi;

n = p.Results.n;

% Default Values
if isempty(Mu)
    Mu = [0 0 1];
end
if isempty(n)
    n = 2^15;
end


%% density
sign_kappa = sign(kappa);
kappa = abs(kappa); % kappa should be non-negative
X = [];

%%

if kappa==0 && bet==0
    Y=Random_Uni_Norm(n);
    return
elseif  bet==0
    Y = Random_vMF_Wood(kappa,Mu,n);
    return
else kappa==0 
   Y = Random_FB4_Beta_GyT(bet,Mu,Psi,n);
    return
end 
%% 
    cFB = integral(@(u)(besseli(0,bet*(1-u.^2)).*...
         exp(kappa*u)),-1,1);% normalizing constant
    cFB0 = integral(@(u)(exp(kappa*u+(-bet)*u.^2)),-1,1); % normalizing constant
    cFB1 = integral(@(u)(exp(kappa*u+(bet)*u.^2)),-1,1); % normalizing constant
    p2 = cFB0/(cFB0+exp(-2*bet)*cFB1); % mixing proportion
    a4 = 2*cFB/(exp(bet)*cFB0+exp(-bet)*cFB1); % acceptance ratio
    
    while length(X) < n
        y3_minus  =FB4(kappa,-bet,ceil(n/a4)); % marginal FB4 distribution
        y3_plus  = FB4(kappa,bet,ceil(n/a4)); % marginal FB4 distribution
        y3_Unif = rand(ceil(n/a4),1);
        y3_mix = y3_minus .* (y3_Unif <= p2)+y3_plus .* (1-(y3_Unif <= p2));
        U = rand(ceil(n/a4),1);
        gFB = besseli(0,bet*(1-U.^2))./cosh(bet*(1-U.^2)); % target dist. eval at u,g(u;kappa,bet,0)
        y3_mix = y3_mix(U<=gFB);
        X = [X;y3_mix];
    end

X = X(1:n);

%% psi vMF; for  bet==0 Ulrich applies
if bet == 0
    V = randn(2,n);
    V = V./kron(ones(2,1),sqrt(sum(V.^2,1)));
    Y = sign_kappa*[kron(ones(1,2),(1-X.^2).^(1/2)).*V' X];
else
    psi = vMF_Circ_Phi(bet*(1-X.^2))+pi*(rand(length(X),1)>1/2);
    Y = sign_kappa*[cos(psi).*sqrt(1-X.^2),sin(psi).*sqrt(1-X.^2),X];
end

%%

Np=[0,0,1]; % z-axis (North Pole)
Mu=Mu/norm(Mu); % should be unit vector
if Mu(3)~= 1
    Ux=cross(Np,Mu);%axis of rotation
    Ux= Ux/norm(Ux);
    thetaX=acos(Mu(3));
    Rg= rotationVectorToMatrix(Ux*thetaX);
    Y=Y*Rg;
end
