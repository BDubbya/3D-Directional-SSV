function Y = Random_FB6(kappa,bet,gamm,varargin)
% RANDOM_FB6 Generates a random sample of size n from the Fisher-Bingham_6 
%            distribution
%
% Y = Random_FB6(kappa,bet,gamm,Mu,Psi,n)
%%
% Examples of correct usage:
%
% n=2^12; 
% kappa = 4.2;
% bet = 4.5 ;
% gamm = -3.5; 
% Psi= pi/2;
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
%
% Examples of correct function construction: 
% Y = Random_FB6(kappa,bet,gamm,Mu,Psi,n)
% Y = Random_FB6(kappa,bet,gamm,Mu,Psi)
% Y = Random_FB6(kappa,bet,gamm,Mu)
% Y = Random_FB6(kappa,bet,gamm)
%%
%
% Required Inputs
% kappa       distribution parameter (kappa \in {Real Numbers})
% bet         distribution parameter (beta >= 0)
% gamm        distribution parameter (gamma \in {Real Numbers})
%
% Optional Inputs
% Mu          distribution parameter (for rotating the North pole to Mu)
% Psi         distribution parameter, (-infinity < Psi < infinity,
%                               for rotation to the new frame of reference)
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
noneg = {'nonnegative'};
class = {'numeric'};
pos_int = {'integer','positive'};
real = {'real'};
noinf = {'finite'};

% Check characteristics for the required and optional parameters
addRequired(p,'kappa',@(x)validateattributes(x,dbl,real));
addRequired(p,'bet',@(x)validateattributes(x,dbl,noneg));
addRequired(p,'gamm',@(x)validateattributes(x,dbl,real));
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'Psi',[],@(x)validateattributes(x,dbl,noinf));
addOptional(p,'n',[],@(x)validateattributes(x,class,pos_int));

p.parse(kappa,bet,gamm,varargin{:});
Mu = p.Results.Mu;
Psi = p.Results.Psi;
n = p.Results.n;

% Default Values
if isempty(Mu)
    Mu = [0 0 1];
end
if isempty(Psi)
    Psi = 0;
end
if isempty(n)
    n = 2^15;
end
%% 
%% 
if kappa==0 && bet==0 && gamm==0
    Y=Random_Uni_Norm(n);
    return
elseif kappa==0 && gamm==0
    Y = Random_FB4_Beta_GyT(bet,Mu,Psi,n);
    return
elseif  kappa==0 && bet==0
    Y = Random_Watsom_Fish(gamm,Mu,n);
    return
elseif  gamm==0
    Y = Random_FB5_GyT(kappa,bet,Mu,Psi,n);
    return
else  bet==0
    Y = Random_FB4(kappa,gamm,Mu,n);
    return
end


%% Bingham kappa=0,gamm=bet
if and(kappa==0,gamm==bet) == 1
    Y = Random_Watson_LW(-bet,[0 1 0],n);
    K = [0 1 0; -1 0 0; 0 0 0];
    R = eye(3)+sin(Psi)*K+(1-cos(Psi))*K*K;% roatation around the North pole
    Y = Y*R;
    
    % rotation (orient the North Pole with mean vector (mu))
    Mu = Mu/norm(Mu); % should be unit vector
    Np = [0,0,1]; % z-axis (North Pole)
    if norm(Mu-Np) ~= 0
        Y=Rotation_Mu(Y,Mu);
    end
    return
end

%%
sign_kappa=sign(kappa)+(kappa==0);
kappa=abs(kappa); 
cFB = integral2(@(u,phi)(exp(bet*(1-u.^2).*cos(2*phi)+kappa*u+gamm*u.^2)),-1,1,0,2*pi);% norm. constant
C_minus = integral(@(u)(exp(kappa*u+(gamm-bet)*u.^2)),-1,1);
C_plus = integral(@(u)(exp(kappa*u+(gamm+bet)*u.^2)),-1,1);
C = exp(bet)*C_minus+exp(-bet)*C_plus;
p1 = exp(bet)*C_minus/C;
% a6 = C/2/cFB;
%% 

X = [];
while length(X)<n
    U = rand(n,1);
    X1_minus = FB4(kappa,gamm-bet,ceil(n)); % marginal FB4 dist %FB4_Watson_gN_F
    X1_plus = FB4(kappa,gamm+bet,ceil(n)); % marginal FB4 dist
    X1 = (U<p1).*X1_minus+(1-(U<p1)).*X1_plus;
    X2 = X1(U <= besseli(0,bet*(1-X1.^2))./cosh(bet*(1-X1.^2)));
    X = [X;X2];
end
X=X(1:n);

%% phi calculating the x1 and x2 components;
Phi = vMF_Circ_Phi(bet*(1-X.^2)); %vMF
% U = randi(4,n,1);
% Phi = acos(PhivMF)/2;
% Phi = Phi.*((U==1)+(U==3))+((U==2)+(U==3))*pi-Phi.*((U==2)+(U==4))+(U==4)*2*pi; % doing 2 pi
cPhi = cos(Phi+Psi);
sPhi = sin(Phi+Psi);

Xs = sqrt(1-X.^2); 
Y = sign_kappa*[cPhi.*Xs,sPhi.*Xs,X];
  
%% rotation (orient the North Pole with mean vector (Mu))
Np=[0,0,1]; % z-axis (North Pole)
Mu=Mu/norm(Mu); % should be unit vector
if norm(Mu-Np) ~= 0
    Y=Rotation_Mu(Y,Mu);
end
return
