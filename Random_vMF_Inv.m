function Y = Random_vMF_Inv(kappa,varargin)
% RANDOM_VMF_INV Generates a random sample of size n from the 
%                von Mises-Fisher distribution
%
% Y = Random_vMF_Inv(kappa,Mu,n)
%%
% Examples of correct usage:
%
% n=2^12; 
% kappa = 4.2;
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
%
% Examples of correct function construction: 
% Y = Random_vMF_Inv(kappa,Mu,n)
% Y = Random_vMF_Inv(kappa,Mu)
% Y = Random_vMF_Inv(kappa)
%%
%
% Required Inputs
% kappa       distribution parameter (kappa \in {Real Numbers})
%
% Optional Inputs
% Mu          distribution parameter, (for rotating the North pole to Mu)
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
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'n',[],@(x)validateattributes(x,class,pos_int));

p.parse(kappa,varargin{:});
Mu = p.Results.Mu;
n = p.Results.n;

% Default Values
if isempty(Mu)
    Mu = [0 0 1];
end
if isempty(n)
    n = 2^15;
end

%% density
% dim of the space / m=3 corresponds to the unit sphere
if kappa == 0;
    Y = Random_Uni_Norm(n); % kappa=0 corresponds to uniform dist. on the sphere
    return
else
    X = vMF_Inv(kappa,n); % cos(theta)
end

%%
psi = 2*pi*rand(n,1);
sX = sqrt(1-X.^2); 
Y = [cos(psi).*sX,sin(psi).*sX,X];

%% rotation (orient mean vector (mu) with the North Pole)
Np=[0,0,1]; % z-axis (North Pole)
Mu=Mu/norm(Mu); % should be unit vector
if norm(Mu-Np) ~= 0
    Y=Rotation_Mu(Y,Mu);
end
return
