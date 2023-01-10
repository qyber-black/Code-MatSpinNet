function create_csv(mat_file_name)
% reads mat file mat_file_name which must be a matspinnet controller file
% outputs csv file mat_file_name.csv one row for each controller and 
% columns as follows:
% fidelity error, transfer time, controller (bias vector), init time, init guess bias (vector)
%
% Run this from the repo root or it may fail.

% SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>, Swansea University
% SPDX-License-Identifier: AGPL-3.0-or-later

load(mat_file_name)
if ~exist('results','var')
  error('no results variable containing controllers')
end

for k=1:length(results) 
  tmp = results{k}; 
  out(k,:) = [tmp.err, tmp.time, tmp.bias', tmp.init_time, tmp.init_bias']; 
end

csvwrite(sprintf('%s.csv',mat_file_name(1:end-4)),out);