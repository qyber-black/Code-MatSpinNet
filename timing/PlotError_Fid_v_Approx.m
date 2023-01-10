function [err,fid_err,tf] = PlotError_Fid_v_Approx(data,m,n,NOFIT)
% function [err,fid_err,tf] = PlotError_Fid_v_Approx(data,m,n,NOFIT)

% SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
% SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 SM Shermer <lw1660@gmail.com>, Swansea University
% SPDX-FileCopyrightText: Copyright (C) 2011-2019 Edmond Jonckheere, University of Southern California
% SPDX-License-Identifier: AGPL-3.0-or-later

err = cell2mat(struct2cell(...
    structfun(@(y)max(abs(y.err)), data,'UniformOutput',false))');

fid_err = cell2mat(struct2cell( ...
    structfun(@(y)GetFidError(y.prob,m,n), data,'UniformOutput',false))');

tf = cell2mat(struct2cell(...
    structfun(@(y)2*y.q*pi/y.omg0, data,'UniformOutput',false))');

% sort error -- otherwise plots are messed up
[err, ind] = sort(err); fid_err = fid_err(ind); tf = tf(ind);

% linear fit and plot results
if exist('NOFIT','var')
    plot(err,fid_err,'.')
else
    fit1 = polyfit(err(err<0.15),fid_err(err<0.15),1);
    plot(err,fid_err,'.',err,polyval(fit1,err))
end
grid on
xlabel('Diophanine approximation error max_k |\theta_k q - p_k|')
ylabel('Transfer fidelity error')

function e = GetFidError(y,m,n)
    e = 1-cellfun(@(x)x(m,n),y);
end

end
