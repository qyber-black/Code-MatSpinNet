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

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Sophie M Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-FileCopyrightText: Copyright (C) 2022 Sean Patrick O'Neil <seanonei@usc.edu>, US Army
  % SPDX-License-Identifier: AGPL-3.0-or-later

  if ~exist('Smu', 'var') || strcmp(Smu,'all') || isempty(Smu)
    Smu = arrayfun(@(mu) GetSmu(obj.N,mu), 1:2*obj.N, 'UniformOutput', false);
  elseif strcmp(Smu, 'couplings')
    Smu = arrayfun(@(mu) GetSmu(obj.N,mu), obj.N+1:2*obj.N, 'UniformOutput', false);
  elseif strcmp(Smu, 'biases')
    Smu = arrayfun(@(mu) GetSmu(obj.N,mu), 1:obj.N, 'UniformOutput', false);
  elseif ismatrix(Smu) && size(Smu,1) == size(Smu,2) && size(Smu,1) == obj.N
    Smu = { Smu };
  elseif isinteger(Smu) && Smu > 0 && Smu <= 2*obj.N
    Smu = { GetSmu(obj.N,Smu) };
  else
    error('Smu not recognised');
  end

  if isempty(info.args.readout)
  % if no readout time window or other constraints are specified then we 
  % evaluate the sensitivity at time results{run}.time

    for run = 1:size(results,2)
    % efficient implementation
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
          end % for s
      end % for j
      dpdd{run} = b*inout;
  
%     original version (more intuitive but *very* slow)
%     n = size(e,2);
%     for k = 1:n
%       for l = 1:n
%         w(k,l) = lambda{run}(k) - lambda{run}(l);
%       end
%     end
%     for s = 1:size(Smu,2)
%       sum(s) = 0;
%       for k = 1:n
%         for l = 1:n
%           VSV = (V(:,k) * V(:,k)') * Smu{s} * (V(:,l) * V(:,l)');
%           s_kl = VSV(info.args.out,info.args.in) * sinc(1/2*results{run}.time*w(k,l));
%           sum_j = 0;
%           for j = 1:n
%             VV = V(:,j) * V(:,j)';
%             sum_j = sum_j + VV(info.args.in,info.args.out) * sin(1/2*results{run}.time * (w(k,j) + w(l,j)));
%           end
%           sum(s) = sum(s) + s_kl * sum_j;
%         end
%       end
%     end
%     sum = 2 * results{run}.time * sum;
%     dpdd{run} = sum';
 
    end % for run

  else % if
  % if time window specified then calculate the bias sensitivity as follow
    H   = obj.H;
    N   = max(size(H));
    I   = eye(N);
    in  = I(:,info.args.in);
    out = I(:,info.args.out);
    for run = 1:length(results)
      sens2 = zeros(length(results),size(Smu,2));
      T  = results{run}.time;
      D  = results{run}.bias;
      dt = info.args.readout(1);
      Hd = H + diag(D);
      [V,lambda] = eig(Hd);
    
      l = diag(lambda);
      o = l*ones(1,N)-ones(N,1)*l';
      io = o.^(-1); io(find(o==0))=0;

%   intuitive version (very slow)
%     for q = 1:size(Smu,2)
%       for k = 1:N
%         for ell = 1:N
%           for p=1:N
%             if l(k) == l(p) && l(ell)==l(p)
%               sens2(run,q) = sens2(run,q) + 0;
%             elseif l(k) == l(ell)
%               sens2(run,q) = sens2(run,q) - 2*(out'*P{k}*Smu{q}*P{k}*in)*(in'*P{p}*out)*( (1/(o(k,p)))*( (T+dt/2)*cos(o(k,p)*(T+dt/2)) - (T-dt/2)*cos(o(k,p)*(T-dt/2))) - (1/(o(k,p)^2))*( sin(o(k,p)*(T+dt/2)) - sin(o(k,p)*(T-dt/2)) ) );
%             else
%               sens2(run,q) = sens2(run,q) + 2*(out'*P{k}*Smu{q}*P{ell}*in)*(in'*P{p}*out)*(1/o(k,ell))*( (T+dt/2)*sinc(o(p,ell)*(T+dt/2)/pi) - (T-dt/2)*sinc(o(p,ell)*(T-dt/2)/pi) - (T+dt/2)*sinc(o(k,p)*(T+dt/2)/pi)+(T-dt/2)*sinc(o(k,p)*(T-dt/2)/pi));
%             end % if
%           end % for p
%         end % for ell 
%       end % for k
%     end % for q
%  dpdd{run} = (sens2(run,:)/dt)';
%
%   much faster alternative
    P = arrayfun(@(n)V(:,n)*V(:,n)', [1:N],'UniformOutput',false);
    Pin    = cell2mat(arrayfun(@(k)P{k}*in,     [1:N],'UniformOutput',false));
    Pout   = cell2mat(arrayfun(@(k)out'*P{k},   [1:N],'UniformOutput',false).');
    Pio    = arrayfun(@(k)in'*P{k}*out,         [1:N]);
    
    T1 = T-dt/2;  T2 = T+dt/2;
    C = io.*(T2*cos(o*T2)-T1*cos(o*T1)) - io.^2.*(sin(o*T2)-sin(o*T1));
    D = T2*sinc(o*T2/pi)-T1*sinc(o*T1/pi);

    for q = 1:size(Smu,2)
      PSP    = Pout*Smu{q}*Pin;
      for k = 1:N
        for ell = 1:N
          for p=1:N
            if l(k) == l(p) && l(ell)==l(p)
               sens2(run,q) = sens2(run,q);
            elseif l(k) == l(ell)
               sens2(run,q) = sens2(run,q) - 2*PSP(k,k)*Pio(p)*C(k,p); 
            else
              sens2(run,q) = sens2(run,q) + 2*PSP(k,ell)*Pio(p)*io(k,ell)*(D(p,ell)-D(k,p));
            end % if
          end % for p
        end % for ell
      end % for k
    end % for q
   dpdd{run} = (sens2(run,:)/dt)';
  end % for run
end % else part

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
  end % subroutine

end % main function