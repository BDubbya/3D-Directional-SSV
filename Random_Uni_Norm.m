function Y = Random_Uni_Norm(varargin)
% RANDOM_UNI_NORM Generates a random sample of size n from the Uniform 
% distribution in S2
%
% Y = Random_Uni_Norm(n)
%%
% Examples of correct usage:
%
% n=2^12; 
%
% Examples of correct function construction: 
% Y = Random_Uni_Norm(n)
% Y = Random_Uni_Norm
%%
%
% Optional Inputs
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
class = {'numeric'};
pos_int = {'integer','positive'};

% Check characteristics for the required and optional parameters
addOptional(p,'n',[],@(x)validateattributes(x,class,pos_int));

p.parse(varargin{:});
n = p.Results.n;

% Default Values
if isempty(n)
    n = 2^15;
end

%% density
R1 = normrnd(0,1,3,n); % X1_i,X2_i,X3_i iid Normal(0,1) for i=1,2,...,n
R1norm = sqrt(sum(R1.^2,1));
R1norm = kron(ones(3,1),R1norm);
Y = R1'./R1norm'; % normalize vectors, so all are of length 1
return
