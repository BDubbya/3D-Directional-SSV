function histY_n=Hist2DSphere(Y,varargin)
% HIST2DSPHERE Finds histogram counts on the HEALPix tesselated unit sphere. 
%              Data points (in standard cartesian format) are sorted, 
%              grouped and counted by which HEALPix pixel they are in.
%
% histY_n=Hist2DSphere(Y,nPix);
%%
% Examples of correct usage:
%
% nSide=2^2;
% nPix= nSide2nPix(nSide);
% kappa = 4.2;
% bet = 4.5 ;
% gamm = -3.5; 
% Psi= pi/2;
%
% Mu=[0 0 1];
% n = 2^12
%
% Y =Random_FB4(kappa,gamm,Mu,n);
%
%
% histY_n = Hist2DSphere(Y,nPix); 
% histY_n = Hist2DSphere(Y);
%%
%
% Required Inputs
% Y      n-by-3 dataframe, where each row is a datapoint in standard
%        cartesian coordinate format, and n is the sample size
%
% Optional Inputs
% nPix   number of HEALPix pixels (resolution parameter), should be of the 
%        form nPix = 3*4^k, for k = 1,2,...
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
addOptional(p,'nPix',[],@(x)validateattributes(x,dbl,{'scalar'},'nPix'));
 
p.parse(Y,varargin{:});
nPix = p.Results.nPix;

% Default Values
if isempty(nPix)
    nPix = 768;
end
if mod(log2(nPix/3),2)~= 0
    errormessage = 'Error:\nnPix must be of the form nPix=3*4^k \n for k=1,2,...';
    error('something:anything',errormessage)
end

%% Histogram Count Code
% normalized:  mean(histY)=1!!
nSide=nPix2nSide(nPix); % HEALPix tesselation parameter
Ypix=vec2pix(nSide,num2cell(Y,2)); % pixel numbers ,'nest',false 
Ysort=sort(Ypix); % Ysort(1:10) Ysort(end-10:end)
histY=zeros(nPix,1);
DYsort=diff(Ysort);
yd=find(DYsort~=0);
histY(Ysort(1)) = yd(1); % first one
ind1=cumsum(DYsort(DYsort~=0))+Ysort(1); % ind1(end-3:end)
Li=length(ind1)-1;
histY(ind1(1:Li))=diff(yd)'; % histY(1:3); histY(end-3:end)
histY(Ysort(end)) =length(Ysort)-yd(end); % last one
histY_n=histY*nPix/length(Y)/4/pi ; %;  integral of histogram is 1;
