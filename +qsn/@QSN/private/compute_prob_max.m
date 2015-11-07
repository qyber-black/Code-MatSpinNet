function compute_prob_max (obj)
  % compute_prob_max (obj) - Compute maximum transition probabilities
  %
  %   obj    - Quantum spin network object
  %

  % Eigenanalysis
  [V,E] = eig (obj.H);
  if isa(E,'sym')
    E = eval(E);
  end
  E = diag(E);

  % Identify eigenspaces by finding differences between eigenvalues
  % that are greater than a certain tolerance
  ind = [0 find(diff(E) > obj.eps)' obj.N];

  % Compute maximum probabilities
  p = zeros(size(obj.H,1));
  for k=1:length(ind)-1
    Vk = V(:,ind(k)+1:ind(k+1));
    p = p + abs(Vk*Vk');
  end
  p = p.^2;

  obj.p_max = p;
end