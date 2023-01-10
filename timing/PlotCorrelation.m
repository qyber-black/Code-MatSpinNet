function PlotCorrelation(ring9,Bias)
% function PlotCorrelation(ring9,Bias)

% SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
% SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 SM Shermer <lw1660@gmail.com>, Swansea University
% SPDX-FileCopyrightText: Copyright (C) 2011-2019 Edmond Jonckheere, University of Southern California
% SPDX-License-Identifier: AGPL-3.0-or-later

[err,fid18,tf]=PlotError_Fid_v_Approx(ring9,1,8,'NOFIT'), hold on
[err,fid27,tf]=PlotError_Fid_v_Approx(ring9,2,7,'NOFIT')
[err,fid36,tf]=PlotError_Fid_v_Approx(ring9,3,6,'NOFIT')
[err,fid45,tf]=PlotError_Fid_v_Approx(ring9,4,5,'NOFIT')
set(gca,'XScale','log','YScale','log')

Corr.p18 = corr([log(err);log(fid18)]')
Corr.p27 = corr([log(err);log(fid27)]')
Corr.p36 = corr([log(err);log(fid36)]')
Corr.p45 = corr([log(err);log(fid45)]')

STR1 = fields(Corr);
STR2 = struct2cell(structfun(@(x)x(1,2),Corr,'UniformOutput',false));
legend(cellfun(@(s,f)sprintf('%s Corr %0.3f',s,f),STR1,STR2,'UniformOutput',false))

title(sprintf('Ring N=9 with Bias=%.0f at n=9',Bias))
