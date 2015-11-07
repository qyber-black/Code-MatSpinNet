function [x,p,q,theta,OK,err,info] = example_attainable()
  % QSN package example
  %
  % Attainability

  % Biased ring for 1-5 transition
  N = 13;
  target = 6;
  B = zeros(1,N);
  B(1,3) = 0;  % Bias
  net = qsn.QSN('ring',N,'XX', B);

  [x,p,q,theta,OK,err,info] = net.find_attainable (1,target,1)

end

%q13 = [118 33 98 ? 334]