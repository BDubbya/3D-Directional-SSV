function gx = Density_FB4_Gy(kappa,gamm,varargin)
% DENSITY_FB4 Calculates and plots Fisher-Bingham_4 density
%
% gx = Density_FB4(kappa,gamm,Mu,resolution)
%%
% Examples of correct usage:
%
% kappa = 4.2;
% gamm = -3.5; 
%
% Mu=[0 0 1];% Mu=[1 -1 1] ;%  Mu=[0 -1 1]; 
% resolution=100;
%
% Examples of correct function construction: 
% gx = Density_FB4(kappa,gamm,Mu,resolution); 
% gx = Density_FB4(kappa,gamm,Mu);
% gx = Density_FB4(kappa,gamm);
%%
%
% Required Inputs
% kappa       distribution parameter (kappa \in {Real Numbers})
% gamm        distribution parameter (gamma \in {Real Numbers})
%
% Optional Inputs
% Mu          distribution parameter,( for rotating the North pole to Mu)          
% resolution  plot resolution parameter
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
vec3D = {'size',[1,3],'finite'};
resInt = {'>',16,'integer'};
real = {'real'};

% Check characteristics for the required and optional parameters
addRequired(p,'kappa',@(x)validateattributes(x,dbl,real,'kappa'));
addRequired(p,'gamm',@(x)validateattributes(x,dbl,real,'gamm'));
addOptional(p,'Mu',[],@(x)validateattributes(x,dbl,vec3D,'Mu'));
addOptional(p,'resolution',[],@(x)validateattributes(x,class,resInt,'resolution'));

p.parse(kappa,gamm,varargin{:});
Mu = p.Results.Mu;
resolution = p.Results.resolution;

% Default Values
if isempty(Mu)
    Mu = [0 0 1]; % North Pole
end
if isempty(resolution)
    resolution = 100;
end

%% density
cFB = integral(@(u)( exp(kappa*u+gamm*u.^2)),-1,1); % normalization constant
delta = pi/resolution; % step size
theta = pi:-delta:0; % altitude
phi = -pi:2*delta:pi; % azimuth
gx = exp(kappa*cos(theta')+gamm*(cos(theta').^2))/cFB/2/pi;

%% plot
[~,theta] = meshgrid(phi,theta);
gx = repmat(gx, [1, size(theta, 1)]);
[x2,y2,z2] = sphere(resolution);
r = gx;
x = (1+r).* x2;
y = (1+r).*y2;
z = (1+r).*z2;
h = surf(x,y,z,r+1);

%% rotation (orient mean vector (mu) with the North Pole)
Np = [0,0,1]; % z-axis (North Pole)
Mu = Mu/norm(Mu); % should be unit vector
if Mu(3)~= 1
    Ux=cross(Np,Mu); % axis of rotation
    Ux= Ux/norm(Ux);
    thetaX=rad2deg(acos(Mu(3)));
    rotate(h,Ux,thetaX )
end

%% display options
view(3)
set(h,'LineStyle','none');
colormap(jet(1024));
colorbar
camlight left
camlight right
lighting phong
A=max(1.5,max(1+gx(:)));
axis([-A A -A A -A A])
grid on

title(['FB4; \kappa =',num2str(kappa), '; \gamma =',num2str(gamm)])
