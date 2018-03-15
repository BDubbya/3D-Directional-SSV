function gx = Density_FB5(kappa,bet,varargin)
% DENSITY_FB5 Calculates and plots Fisher-Bingham_5 density
%
% gx = Density_FB5(kappa,bet,Mu,Psi,resolution)
%%
% Examples of correct usage:
%
% kappa = 4.2;
% bet = 4.5 ; 
% Psi= pi/2;
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
% resolution=100;
%
% Examples of correct function construction: 
% gx = Density_FB5(kappa,bet,Mu,Psi,resolution) 
% gx = Density_FB5(kappa,bet,Mu,Psi);
% gx = Density_FB5(kappa,bet,Mu);
% gx = Density_FB5(kappa,bet);
%%
%
% Required Inputs
% kappa       distribution parameter (kappa \in {Real Numbers})
% bet         distribution parameter (beta >= 0)
%
% Optional Inputs
% Mu          distribution parameter, (for rotating the North pole to Mu)
% Psi         distribution parameter, (-infinity < Psi < infinity)           
% resolution  plot resolution parameter
%
% Authors: Gy.Terdik, B.Wainwright
%
% Copyright 2018 Gy.Terdik
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

%% Verify/Default parameter values

% Required parameter characteristics
p = inputParser;
dbl = {'double'};
noneg = {'nonnegative'};
class = {'numeric'};
vec3D = {'size',[1,3],'finite'};
longPsi = {'>=',0,'<=',2*pi,'finite'};
resInt = {'>',16,'integer'};
real = {'real'};

% Check characteristics for the required and optional parameters
addRequired(p,'kappa',@(x)validateattributes(x,dbl,real));
addRequired(p,'bet',@(x)validateattributes(x,dbl,noneg));
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,vec3D,'Mu'));
addOptional(p,'Psi',[],@(x)validateattributes(x,dbl,longPsi,'Psi'));
addOptional(p,'resolution',[],@(x)validateattributes(x,class,resInt,'resolution'));

p.parse(kappa,bet,varargin{:});
Mu = p.Results.Mu;
Psi = p.Results.Psi;
resolution = p.Results.resolution;

% Default Values
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
g5=@(u,phi)(exp(bet*(1-u.^2).*cos(2*phi)+ kappa*u));
cFB=integral2(@(u,phi)(exp(bet*(1-u.^2).*cos(2*phi)+ kappa*u)),-1, 1 ,0,2*pi); % norm. constant
delta = pi/resolution; % step size
theta = 0:delta:pi; % colatitude
phi = 0:2*delta: 2*pi; % longitude

%% plot
[theta,phi] = meshgrid(theta,phi);
phi1=phi-Psi;
gx=g5(cos(theta),phi1)/cFB;
x2 = sin(theta).*cos(phi);
y2 = sin(theta).*sin(phi);
z2 = cos(theta);

r = gx;
x = (1+r).*x2;
y = (1+r).*y2;
z = (1+r).*z2;
h = surf(x,y,z,r+1); 

%% rotation (orient mean vector (mu) with the North Pole)
Np=[0,0,1]; % z-axis (North Pole)
Mu=Mu/norm(Mu); % should be unit vector
if Mu(3)~= 1
    Ux=cross(Np,Mu);%axis of rotation
    Ux= Ux/norm(Ux);
    thetaX=rad2deg(acos(Mu(3)));
    rotate(h,Ux,thetaX )
end

%% display options
view(3)
set(h,'LineStyle','none')
colormap(jet(1024));
colorbar
camlight left
camlight right
lighting phong
A=max(1.5,max(1+gx(:)));
axis([-A A -A A -A A])
grid on

title(['FB5; \kappa =',num2str(kappa), '; \beta=',num2str(bet), '; ,\psi=',num2str(Psi)])
xlabel('x-axis')
ylabel('y-axis')