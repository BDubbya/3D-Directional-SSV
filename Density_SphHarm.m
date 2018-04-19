function r2 = Density_SphHarm(L,m,varargin)
% DENSITY_SPHHARM   Calculates and plots the density of the complex valued 
%                   and real valued spherical harmonics, for different 
%                   degrees and orders.
%
% r2 = Density_SphHarm(L,m,R,Mu,Psi,resolution);
%
% Required Inputs
% L             degree of the spherical harmonic
% m             order of the spherical harmonic
%
% Optional Inputs
% R             logical indicator (true=real/false=complex mod square spherical harmonics)
% Mu            distribution parameter, mean vector (Mu)
% Psi           distribution parameter, (-infinity < Psi < infinity)
% resolution    plot resolution parameter
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
% if abs m <= L then (L,m) only otherwise all

%% Verify/Default parameter values

% Required parameter characteristics
p = inputParser;
class = {'numeric'};
pos_int = {'integer','positive'};
m_ord = {'integer','>=',-L,'<=',L};
yes_no = {'logical'};
bin = {'binary'};
dbl = {'double'};
vec3D = {'size',[1,3],'finite'};
longPsi = {'>=',0,'<=',2*pi,'finite'};
resInt = {'>',16,'integer'};

% Check characteristics for the required and optional parameters
addRequired(p,'L',@(x)validateattributes(x,class,pos_int));
addRequired(p,'m',@(x)validateattributes(x,class,m_ord));
addOptional(p,'R',[],@(x)validateattributes(x,yes_no,bin));
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,vec3D,'Mu'));
addOptional(p,'Psi',[],@(x)validateattributes(x,dbl,longPsi,'Psi'));
addOptional(p,'resolution',[],@(x)validateattributes(x,class,resInt,'resolution'));

p.parse(L,m,varargin{:});
R = p.Results.R;
Mu = p.Results.Mu;
Psi = p.Results.Psi;
resolution = p.Results.resolution;

% Default Values
if isempty(R)
    R = true;
end
if isempty(Mu)
    Mu = [0 0 1];
end
if isempty(Psi)
    Psi = 0;
end
if isempty(resolution)
    resolution = 100;
end

%% density
delta = pi/resolution; % step size
theta = 0:delta:pi; % colatitude
phi = 0:2*delta:2*pi; % longitude
[phi,theta] = meshgrid(phi,theta);
x2 = sin(theta).*cos(phi);
y2 = sin(theta).*sin(phi);
z2 = cos(theta);

P_LM = legendre(L,cos(theta(:,1)));

P_LM = P_LM(abs(m)+1,:)'; % Legendre polynomials
P_LM = repmat(P_LM, [1, size(theta, 1)]);
N_LM = sqrt((2*L+1)/4/pi*factorial(L-abs(m))/factorial(L+abs(m))); % norm. const

% base spherical harmonic functions square
if R ~= 1
    r2= 2*N_LM^2*(P_LM ).^2;
else
    if m>=0
        r2 = 2*N_LM^2*(P_LM .* cos(m*(phi-Psi))).^2;
    else
        r2 = 2*N_LM^2*(P_LM .* sin(abs(m)*(phi-Psi))).^2;
    end
end
% map to sphere surface
%     abra(phi1,theta1,r2);

%% plot
x = (1+r2).* x2;
y = (1+r2).*y2;
z = (1+r2).*z2;
h = surf(x,y,z,r2+1);

%% rotation (orient mean vector (Mu) with the North Pole)
Np=[0,0,1]; % z-axis (North Pole)
Mu=Mu/norm(Mu); % should be unit vector
if Mu(3)~= 1
    Ux=cross(Np,Mu);%axis of rotation
    Ux= Ux/norm(Ux);
    thetaX=rad2deg(acos(Mu(3)));
    rotate(h,Ux,thetaX );
end
%
%% display options
view(3)
set(h,'LineStyle','none');
colormap(jet(1024));
colorbar
camlight left
camlight right
lighting phong
A=max(1.5,max(1+r2(:)));
axis([-A A -A A -A A])
grid on

title(['L = ', num2str(L), ', M = ', num2str(m)])%% plot
end
