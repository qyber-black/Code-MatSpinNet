function PlotScaling_BiasRing9(ring9)
% plots minimum transfer times vs infidelity for 1-8, 2-7, 3-6, 4-5 transfer
% input:  ring9 -- data structure with subfields which can be named
%         anything and there can be any number but each field must have 
%         the following subfields
%         omg0     -- reference transition
%         theta    -- vector to be approximated
%         p        -- matrix with columns corresponding to numerators
%         q        -- vector of denominators of sim. diophan. approx.
%         parity_ok-- vector of flags indicating if p satisfies parity
%         prob     -- list of matrices of transition probabilities at time
%                     tf = 2*q*pi/omg0
%          
% output: Plot of shortest transfer times found vs infidelity 
%
% This version only works for ring9 at the moment.

% SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
% SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 SM Shermer <lw1660@gmail.com>, Swansea University
% SPDX-FileCopyrightText: Copyright (C) 2011-2019 Edmond Jonckheere, University of Southern California
% SPDX-License-Identifier: AGPL-3.0-or-later

function e = GetFidError(y,m,n)
    e = 1-cellfun(@(x)x(m,n),y);
end
function b = Min(a)
    if isempty(a)
        b = NaN;
    else
        b = min(a);
    end
end

tf = cell2mat(struct2cell(...
    structfun(@(y)2*y.q*pi/y.omg0, ring9,'UniformOutput',false))');

fid18 = cell2mat(struct2cell( ...
    structfun(@(y)GetFidError(y.prob,1,8), ring9,'UniformOutput',false))');
fid27 = cell2mat(struct2cell( ...
    structfun(@(y)GetFidError(y.prob,2,7), ring9,'UniformOutput',false))');
fid36 = cell2mat(struct2cell( ...
    structfun(@(y)GetFidError(y.prob,3,6), ring9,'UniformOutput',false))');
fid45 = cell2mat(struct2cell( ...
    structfun(@(y)GetFidError(y.prob,4,5), ring9,'UniformOutput',false))');

TGT =[0.1 0.06 0.03 0.01 0.006 0.003 0.0013]; 
T18 = arrayfun(@(x)Min(tf(find(fid18<x))),TGT);
T27 = arrayfun(@(x)Min(tf(find(fid27<x))),TGT);
T36 = arrayfun(@(x)Min(tf(find(fid36<x))),TGT);
T45 = arrayfun(@(x)Min(tf(find(fid45<x))),TGT);

lin18 = polyfit(log10(TGT),log10(T18),1);
lin27 = polyfit(log10(TGT),log10(T27),1);
lin36 = polyfit(log10(TGT),log10(T36),1);
lin45 = polyfit(log10(TGT),log10(T45),1);

clf
loglog(TGT,T18,'bo','MarkerFaceColor','b'), hold on, grid on
loglog(TGT,10.^polyval(lin18,log10(TGT)),'b:','LineWidth',1)
loglog(TGT,T27,'gs','MarkerFaceColor','g')
loglog(TGT,10.^polyval(lin27,log10(TGT)),'g:','LineWidth',1)
loglog(TGT,T36,'r>','MarkerFaceColor','r')
loglog(TGT,10.^polyval(lin36,log10(TGT)),'r:','LineWidth',1)
loglog(TGT,T45,'cp','MarkerFaceColor','c')
loglog(TGT,10.^polyval(lin45,log10(TGT)),'c:','LineWidth',1)

xlabel('Target Infidelity')
ylabel('Transfer Time (1/J)')
legend('p_{18}','','p_{27}','','p_{36}','','p_{45}','')

end
