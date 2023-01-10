function err_act = ComputeError(obj,ispin,jspin,t,TITLE)
% function err_act = ComputeError(obj,ispin,jspin,t,TITLE)

% SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
% SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 SM Shermer <lw1660@gmail.com>, Swansea University
% SPDX-FileCopyrightText: Copyright (C) 2011-2019 Edmond Jonckheere, University of Southern California
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

% sanity check
if norm(Op(P,e,N)-obj.H)>1e-10
    error('something is wrong here')
end

% calculate input/output state projection onto eigenspaces
Pij = cell2mat (arrayfun (@(k) P{k}(ispin,jspin), (1:Ns)', 'UniformOutput', false));

% eliminate dark subspaces
s   = sign(Pij);
nonzero_proj = find(s);
ee = e(nonzero_proj);
ss = s(nonzero_proj);
Pij= Pij(nonzero_proj);

% calculate maximum transition probability
max_p = obj.prob;
max_p = max_p(ispin,jspin);

% calculate actual transfer error
err_act = arrayfun(@(t)1-abs(exp(-1i*t*ee.')*Pij).^2/max_p, t);
[err1,imin1] = min(err_act);
tmin1 = t(imin1);
% optional verification
% u1 = abs(exp(-1i*tmin1*ee.')*Pij).^2;

% calculate proxy error
%phase = (ee-ee(1))/pi*t + (ss-ss(1))/2*ones(size(t));
%err_phs = sum(min(mod(phase,2),mod(-phase,2)));
%[err2,imin2] = min(err_phs);
%tmin2 = t(imin2);
% u2 = abs(exp(-1i*tmin2*ee.')*Pij).^2;

figure(1)
%subplot(2,1,1)
semilogy(t,err_act,tmin1,err1,'rp'), grid on
xlabel('time (1/J)'), ylabel('error 1 -- p(t)/p_{max}')
if exist('TITLE','var')
    title(TITLE)
end
%subplot(2,1,2),
%semilogy(t,error2,tmin2,err2,'rp'), grid on
%xlabel('time'), ylabel('proxy error')

function H = Op(Pi,e,N)
H = zeros(N);
for k=1:length(e)
    H = H + e(k)*Pi{k};
end
