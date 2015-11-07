function [F,dF] = dexpma(M,dM,v)
  % [F,dF] = dexpma(M,dM,v)
  %
  % Gradient helper for eal_static_st
  N   = length(M);
  AM  = [M zeros(N);dM M];
  if exist('v','var')
    PSI = expv(AM,[v;v]);
    F   = PSI(1:N);
    dF  = PSI(N+1:end);
  else
    PSI = expm(AM);
    F   = PSI(1:N,1:N);
    dF  = PSI(N+1:2*N,1:N);
  end
end