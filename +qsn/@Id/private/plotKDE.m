function plotKDE (KDE, varargin)
% Plot KDE
%   KDE      - KDE
%   varargin - if set to 1, plot 3D contour, 2 plot mesh, other nothing

% SPDX-FileCopyrightText: Copyright (C) 2011-2019 Frank C Langbein <frank@langbein.org>, Cardiff University
% SPDX-License-Identifier: AGPL-3.0-or-later

if size(varargin, 2) == 1
  mode = varargin{1};
else
  mode = 0;
end

figure(1);
clf;
[H,X,Y] = hist (KDE);
H(isinf(H)) = -1;
H(isnan(H)) = -1;
if mode == 1
  contour3 (X, Y, H',30);
  hold on;
elseif mode == 2
  mesh (X, Y, H');
  hold on;
end
[c,h] = contourf (X, Y, H',30);
colorbar; shading flat
hold off;
view(2);

figure(2);
clf; 
plot (KDE, '*bWS-g'); axis square, axis tight
