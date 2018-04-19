function Y= Random_Y3_2Compl_square(varargin)
% RANDOM_Y3_2COMPL_SQUARE Generates a random sample of size n from the 
%                         density characterized by the complex valued 
%                         spherical harmonic of degree l=3 and order m=2
%
% Y = Random_Y3_2compl_square(Mu,n)
%%
% Examples of correct usage:
%
% n=2^12; 
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
%
% Examples of correct function construction: 
% Y = Random_Y3_2compl_square(Mu,n)
% Y = Random_Y3_2compl_square(Mu)
% Y = Random_Y3_2compl_square
%%
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
noinf = {'finite'};

% Check characteristics for the required and optional parameters
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'Psi',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'n',[],@(x)validateattributes(x,class,pos_int));

p.parse(varargin{:});
Mu = p.Results.Mu;
n = p.Results.n;

% Default Values
if isempty(Mu)
    Mu = [0 0 1];
end
if isempty(n)
    n = 2^12;
end

%% spherical harmonics parameters: fix degree and order
ld=3; % degree of spherical harmonic
mo=2; % order of spherical harmonic


%% envelope parameters
% X~e(x) = Beta(alpha=2,beta=1.75) - Envelope Distribution
% U~Unif(0,1)
% Reject X if U > f(X)/e(X) - Rejection Criteria

M = 1.074; % acceptance ratio/expect M trials per accepted X (Ulrich p.160)
alpha = 3.08; % shape parameter for Beta distribution (envelope)
beta = 2.5249; % shape parameter for Beta distribution (envelope)
n= floor(M*n); % total number of trials necessary

%%  Step 1: calculate x3 (from (x1,x2,x3))

X2 = betarnd(alpha,beta,n,1);  %random n-by-1 vector of size n from X~e(x)

% Squared, fully normalized associated legendre function of degree ld and 
% orders mo=0,1,...,ld
glm = (legendre(ld,X2,'norm')).^2; % (g_l^m(u))^2
X1 = 2*glm(mo+1,:);

% envelope density function
g_X2 = betapdf(X2,alpha,beta); 

% random n-by-1 vector from uniform(0,1)
U1 = rand(n,1); 
%accepted observations
X_accept = X2(U1<=X1'./g_X2); 

%%  Step 2
% Positive or Negative w/ Prob(1/2)
U2 = rand(length(X_accept),1);
W=(ones(length(X_accept),1)-2*(U2<0.5)).*X_accept;

%%  Step 3
% calculate x1 and x2
V = randn(length(W),2);
V = V./kron(ones(1,2),sqrt(sum(V.^2,2)));
Y = [ kron(ones(1,2),(1-W.^2).^(1/2)).*V W];

%% rotation (orient the North Pole with mean vector (mu))
Mu = Mu/norm(Mu); % should be unit vector
Np=[0,0,1];
if norm(Mu-Np) ~= 0
    Y=Rotation_Mu(Y,Mu);
end
