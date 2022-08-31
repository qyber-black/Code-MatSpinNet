function [dpdd,M,lambda] = bias_sensitivity(obj,results,info,Smu)
  % [dpdd,M,lambda] = obj.bias_sensitivity(results,info,Smu)
  %
  % This function computes the partial derivative relative of the probability of
  % transmission of the excitation from |in> to |out> in an amount of time tf in a spin network of
  % Hamiltonian H perturbed as H+Smu, where Smu is a matrix reflecting the perturbation of the couplings.
  %
  % On input
  %
  % obj     : QSN object
  % results : results cell array from bias_control
  % info    : info cell array from bias_control
  % Smu     : structure of the perturbation
  %            1-n    : pertubation on diagonal
  %            n+1-2n : pertubation on ring couplings
  %            'all'  : 1-2n pertubations above
  %            'couplings': coupling pertubations only
  %            'biases': bias pertubations only
  %            matrix : Smu matrix
  %
  % On output
  %
  % dpdd   : partial derivative of probability of transfer relative to Smu for all results
  % M      : the "sign divided difference matrix" for all results
  % lambda : the eigenvalues of the Hamiltonian for all results

  if ~exist('Smu', 'var') || strcmp(Smu,'all') || isempty(Smu)
    Smu = arrayfun(@(mu) GetSmu(obj.N,mu), [1:2*obj.N], 'UniformOutput', false);
  elseif strcmp(Smu, 'couplings')
    Smu = arrayfun(@(mu) GetSmu(obj.N,mu), [obj.N+1:2*obj.N], 'UniformOutput', false);
  elseif strcmp(Smu, 'biases')
    Smu = arrayfun(@(mu) GetSmu(obj.N,mu), [1:obj.N], 'UniformOutput', false);
  elseif ismatrix(Smu) && size(Smu,1) == size(Smu,2) && size(Smu,1) == obj.N
    Smu = { Smu };
  elseif isinteger(Smu) && Smu > 0 && Smu <= 2*obj.N
    Smu = { GetSmu(obj.N,Smu) };
  else
    error('Smu not recognised');
  end

  for run = 1:size(results,2)

    H = obj.H + diag(results{run}.bias);
    [V,e] = eig(H);
    lambda{run} = diag(e);

    S = cellfun(@(x) results{run}.time * V' * x * V, Smu, 'UniformOutput', false);
    SINC = sinc((lambda{run} * ones(1,obj.N) - ones(obj.N,1) * lambda{run}') * results{run}.time / (2*pi));
    SUML = (lambda{run} * ones(1,obj.N) + ones(obj.N,1) * lambda{run}') * results{run}.time / 2;
    %outV = out_vec'*V;
    outV = V(info.args.out,:);
    %Vin  = V'*in_vec;
    Vin = V(info.args.in,:)';
    inout= 2 * (outV' .* conj(Vin));
    for j = 1:obj.N
      SIN  = sin(SUML-lambda{run}(j) * results{run}.time);
      M{run,j} = SIN .* SINC;
      for s = 1:length(S)
        b(s,j) = outV * (M{run,j} .* S{s}) * Vin;
      end
    end
    dpdd{run} = b*inout;

    % Ultra slow version
    %n = size(e,2);
    %for k = 1:n
    %  for l = 1:n
    %    w(k,l) = lambda{run}(k) - lambda{run}(l);
    %  end
    %end
    %for s = 1:size(Smu,2)
    %  sum(s) = 0;
    %  for k = 1:n
    %    for l = 1:n
    %      VSV = (V(:,k) * V(:,k)') * Smu{s} * (V(:,l) * V(:,l)');
    %      s_kl = VSV(info.args.out,info.args.in) * sinc(1/2*results{run}.time*w(k,l));
    %      sum_j = 0;
    %      for j = 1:n
    %        VV = V(:,j) * V(:,j)';
    %        sum_j = sum_j + VV(info.args.in,info.args.out) * sin(1/2*results{run}.time * (w(k,j) + w(l,j)));
    %      end
    %      sum(s) = sum(s) + s_kl * sum_j;
    %    end
    %  end
    %end
    %sum = 2 * results{run}.time * sum;
    %dpdd{run} = sum';

  end

  function Smu = GetSmu(n,mu)
    Smu = zeros(n,n);
    if 0 < mu && mu <= n
      Smu(mu,mu) = 1;
    elseif n < mu && mu < 2*n
      mu = mu - n;
      Smu(mu,mu+1) = 1;
      Smu(mu+1,mu) = 1;
    elseif mu == 2*n
      mu = mu - n;
      Smu(1,n) = 1;
      Smu(n,1) = 1;
    else
      error('Illegal mu');
    end
  end

end