function Y = Random_FB4_Wood(kappa,gamm,varargin)
% RANDOM_FB4_WOOD Generates a random sample of size n from the 
%                 Fisher-Bingham_4 distribution
%
% Y = Random_FB4_Wood(kappa,gamm,Mu,n)
%
%%
% Examples of correct usage:
%
% n=2^12; 
% kappa = 4.2;
% gamm = -3.5; 
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
%
% Examples of correct function construction: 
% Y = Random_FB5_Kent(kappa,gamm,Mu,n)
% Y = Random_FB5_Kent(kappa,gamm,Mu)
% Y = Random_FB5_Kent(kappa,gamm)
%%
%
% Required Inputs
% kappa       distribution parameter (kappa \in {Real Numbers})
% gamm        distribution parameter (gamma \in {Real Numbers})
%
% Optional Inputs
% Mu          distribution parameter (for rotating the North pole to Mu)
% n           sample size (where n is a positive integer)
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
class = {'numeric'};
pos_int = {'integer','positive'};
real = {'real'};
noinf = {'finite'};

% Check characteristics for the required and optional parameters
addRequired(p,'kappa',@(x)validateattributes(x,dbl,real));
addRequired(p,'gamm',@(x)validateattributes(x,dbl,real));
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'n',[],@(x)validateattributes(x,class,pos_int));

p.parse(kappa,gamm,varargin{:});
Mu = p.Results.Mu;
n = p.Results.n;

% Default Values
if isempty(Mu)
    Mu = [0 0 1];
end
if isempty(n)
    n = 2^12;
end

%%
% cFB = integral(@(u)(exp(kappa*u+gamm*u.^2)),-1,1); % normalizing constant
% p2 = 1/2; % mixing proportion
a4 = 1; % acceptance ratio
sign_kappa = sign(kappa)+(kappa==0);
kappa = abs(kappa); % kappa should be non-negative

%%
if gamm==0; 
   Y = Random_vMF_Wood(kappa,Mu,n);
   return % when gamma==0 FB4 is the von Mises-Fisher distribution
elseif kapa==0 
    Y = Random_Watsom_Fish(gamm,Mu,n);
    return
else
    X = FB4_Wood(kappa,gamm,ceil(n)); % marginal FB4 distribution
end

%% psi iid uniform
psi = 2*pi*rand(length(X),1);
Xs = sqrt(1-X.^2);%sin(acos(X));
Y = sign_kappa*[cos(psi).*Xs,sin(psi).*Xs,X];

%% rotation (orient the North Pole with mean vector (Mu))
Np=[0,0,1]; % z-axis (North Pole)
Mu=Mu/norm(Mu); % should be unit vector
if norm(Mu-Np) ~= 0
    Y=Rotation_Mu(Y,Mu);
end
