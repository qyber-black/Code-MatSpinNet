function KDE = kdeH (S, D, Range, H_N, varargin);
% Match a Hamiltonian KDE to exerpimental data
%
%    KDE = kdeH (S, D, Range, H_N, varargin);
%
%    S        - QSN object
%    D        - Experimental data
%    Range    - Data range
%    H_N      - Number of samples for H reconstruction
%    varargin - H Resamples (omit to not resample)
%
%    KDE     - Hamiltonian KDE

% SPDX-FileCopyrightText: Copyright (C) 2011-2019 Frank C Langbein <frank@langbein.org>, Cardiff University
% SPDX-License-Identifier: AGPL-3.0-or-later

% Setup sampling
persistent Base;
if isempty(Base)
  Base = [ 2   3   5   7   11  13  17  19  23  29  31  37  41  43  47  53  59  61  67  71 73  79  83  89  97  101   103   107   109   113   127   131   137   139   149   151   157   163   167   173 ];
end

% Sample Hamiltonian space
dim = size(Range, 1);
if dim > size(Base, 2)
  error ('Too many unknowns for H sampling as currently implemented');
end
% FIXME: ld, anisotropic sampling here
L = Range(1,2) - Range(1,1) + 1;
seed = 0;
for l = Range(1,1):Range(1,2);
  %ds = Range(2,1) + i4_to_hammersley_sequence(1, H_N, [0], seed, [1], [2]) * (Range(2,2) - Range(2,1));
  ds = Range(2,1) + (0:1/(H_N-1):1) * (Range(2,2) - Range(2,1));
  Hs(1:2,1+(l-Range(1,1)) * H_N:(l+1-Range(1,1)) * H_N)  = [ l * ones(1,H_N); ds];
  seed = seed + H_N;
end
H_N = H_N * L;

% Compute probabilities of Hamiltonian samples
W = zeros (1, H_N);
parfor n = 1:H_N
  W(n) = prob_binomial_log (Hs(:,n), D);
end
% FIXME: normalisation not needed?  W = W / sum(W);

% Match PDF
bw = (Range(2,2) - Range(2,1)) / H_N * L * 0.6;
KDE = kde(Hs,bw,W,'Gaussian');

% Adjust bandwidth
KDE = ksize(KDE,'hall');

% Resample
if size(varargin,2) == 1
    H_Resamples = varargin{1};
else
    H_Resamples = [];
end

if ~isempty(H_Resamples)

  % Resample according to KDE
  % FIXME: ld, anisotropic sampling here
  Hs = sample (KDE, H_Resamples);

  % Compute probabilities of Hamiltonian re-samples
  W = zeros (1, H_Resamples);
  parfor n = 1:H_Resamples
    W(n) = prob_binomial_log (Hs(:,n), D);
  end
  % FIXME: normalisation not needed?  W = W / sum(W);

  % Match PDF
  bw = (Range(:,2) - Range(:,1)) / (H_Resamples.^(1/dim));
  KDE = kde(Hs,bw,W,'Gaussian');

  % Adjust bandwidth
  KDE = ksize(KDE,'Hall');

end