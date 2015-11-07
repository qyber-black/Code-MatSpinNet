function setup(a)
  % Setup QSN code
  %  setup('build') - compile code
  %  setup()        - set paths (or any other argument)
  if exist('a','var') && strcmp(a,'build')
    % Compile KDE code
    cd '@kde/mex'
    makemex
    cd '../..'
  end
  addpath '.' './GA' 
end
