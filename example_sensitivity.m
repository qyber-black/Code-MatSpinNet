function example_sensitivity(C,Tm)
  % QSN package example
  %
  % Sensitivity to pertubations

  figure(3);

  T = 0:0.5:Tm+5;
  T = Tm-5:0.5:Tm+5;
  N = 11;
  init = 1;
  target = 5;

  pert = 0.01;

  %C = [ 10.4516   11.3318   11.6275   11.3318   10.4516  229.1179   67.8795   55.7833   55.7833   67.8795  229.1179 ];

  for r = 1:100
    J = 1 * (1 + pert * (rand(1,N) * 2 - 1));
    chain = qsn.QSN('ring', J, 'XX', C);

    P = chain.prob(T);

    for in = 1:N
      if in == init
        p = arrayfun(@(t) P{t}(in,target),1:size(P,2));
        plot (T,p,'b');
      elseif in ~= target
        p = arrayfun(@(t) P{t}(in,target),1:size(P,2));
        plot (T,p,'m');
      end
      hold on;
    end

  end

  J = 1 * (1 + pert * ones(1,N));
  chain = qsn.QSN('ring', J, 'XX', C);
  P = chain.prob(T);
  p = arrayfun(@(t) P{t}(init,target),1:size(P,2));
  plot (T,p,'r');

  J = 1 * (1 - pert * ones(1,N));
  chain = qsn.QSN('ring', J, 'XX', C);
  P = chain.prob(T);
  p = arrayfun(@(t) P{t}(init,target),1:size(P,2));
  plot (T,p,'g');

  hold off;

end