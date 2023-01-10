function [P,e] = EigDecomp(obj)
%
% computes decomposition of Hamiltonian obj.H into distinct eigenspaces
% H = sum_k e_k P_k where P_k are projectors onto eigenspaces and e_k
% is the corresponding eigenvalue
%
% input:  obj -- spinnet object
% output: P   -- {P1,P_2,...} list of projectors
%         e   -- [e1,e2,....] vector of eigenvalues

% SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University\
% SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Sophie M Shermer <lw1660@gmail.com>, Swansea University\
% SPDX-License-Identifier: AGPL-3.0-or-later

N = obj.N;
% compute eigendecomposition
[V,E] = eig(obj.H);
e = diag(E);

% identify distinct eigenspaces by finding differences between
% eigenvalues that are greater than a certain tolerance

ind = [0 find(diff(e)>1e-10)' N];
Ns  = length(ind)-1;
for k=1:Ns
    Vk = V(:,ind(k)+1:ind(k+1));
    P{k} = Vk*Vk';
end
e = e(ind(1:end-1)+1);

% sanity check -- check that norm(sum_k e_k P_i - obj.H)<tol
tol = 1e-10;  % this value may be too small for very large systems
if norm(Op(P,e,N)-obj.H)>tol
    error('something is wrong here')
end

% utility routine
function H = Op(Pi,e,N)
H = zeros(N);
for k=1:length(e)
    H = H + e(k)*Pi{k};
end