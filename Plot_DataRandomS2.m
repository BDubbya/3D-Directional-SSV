function h=Plot_DataRandomS2(Y)
% PLOT_DATARANDOMS2 Plots data points on the unit sphere (S2)
%
% h=Plot_DataRandomS2(Y)
%%
% Examples of correct usage:
%
% n=2^12; 
% Y = Random_Uni_Inv(n);
% h=Plot_DataRandomS2(Y);
%%
%
% Required Inputs
% Y     n-by-3 matrix/dataframe, where each row is a datapoint in standard 
%       cartesian coordinate format, and N is the sample size
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

% Check characteristics for the required and optional parameters
addRequired(p,'Y',@(x)validateattributes(x,dbl,{'ncols',3},'Y'));

%%
Rp=Y'; % convert form to (1:3,:) 
% figure    
resolution = 50; % graphical parameter
c = ones(resolution); 
[x,y,z] = sphere(resolution);
h2 = surf(x,y,z,1*c);
set(h2,'LineStyle','none');
colormap([0 1 0])
axis tight; hold on
% plot3(x1,y1,z1,'.','MarkerSize',15,'color',[0,0,0]) %nonuniform
h=plot3(Rp(1,:)',Rp(2,:)',Rp(3,:)','.','MarkerSize',4,'color',[0,0,0]); %nonuniform
hold off
view(3);
alpha(h2,.9);
