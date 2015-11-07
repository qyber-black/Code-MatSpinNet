function out = ComputeTime(obj,ispin,jspin,sfact,runs)

% function ComputeTime(obj,ispin,jspin,sfact,runs) estimates time required
% to achieve transfer from spin i to spin j in spin network obj.
%
% using simulateneous diaphantine approximation and LLL algorith with
%     x = [rand(1,length(theta))*sfact(1), sfact(2)];
% attemps 500 runs to find p with correct parity.
% if none of the runs succeeds tries combinations of p.
%

% ring07 = qsn.QSN('ring',7,'XX',zeros(1,7));
% out2 = qsn.ComputeTime(1,5,[1,1e-6]) should output something like
% out2 =
%
%              omg0: 2.4450
%             theta: [2x1 double]
%                 p: [2x1 double]
%                 q: 745
%                tf: 1.9145e+03
%    err_dio_approx: 1.5085e-04
%           err_act: 0.0126
%       err_act_rel: 0.0307
%
% omg0 is the reference transition, theta = 2*omg/omg0 excluding omg0, p and
% q are the simultaneous diophantine approximation found.  If it does not
% find a primary approximation so that p satisfy the parity constraints,
% then it will try to combine solutions to satisfy the parity constaints and
% return instead a list of p2 and q2 with the corresponding tf and errors
% for each of the combined solutions or subconvergents as Edmond called them.
%
% tf is the corresponding transfer time, err_dio_approx is the diophantine
% approximation error, err_act is the actual transfer error, i.e., the maximum
% transfer probability minus the actually transfer probability achieved.
% err_act_rel the relative error.
%
% It needs the routine EigDecomp.m which should be included as a method in qsn
% maybe because it's used in various contexts.  Note that EigDecomp has a
% tolerance that may need to be adjust based on the size of the network to avoid
% that it erroneously combines subspaces that should be distinct, which may
% happen if the tolerance tol is too large.  On the other hand if tol is too
% small it will take eigenvectors that should be combined into a single eigenspace
% to be separate and then all the bounds are useless, so I would say when in doubt,
% a slightly larger tol is better.
%
% It's not perfect and not very robust but it's a start, I think, and with a bit
% of effort you could turn it into a qsn method but would need some work because
% at the moment the scaling factor sfact determines the magnitude of the error
% but only roughly.  Ideally, you'd want a routine that calcuates the minimum time
% needed for the absolute or relative transfer error to be below a given threshold.
% Also, selecting subconvergents is an issue -- at the moment I'm just returning
% them all for the user to decide with is the best tradeoff between transfer time
% and transfer fidelity.

% I'd suggest you compare my routine to your is_attainable routine.  There is a lot
% of overlap.  Some differences
%
% 1) selection of reference transition
% 2) I repeat the simultaneous diophantine approx 500 times
% 3) I calculate the actual transfer time and transfer error
% 4) I've added the combination of subconvergents to satisfy parity constraints if
% I cannot find a primary approximation that works.
%
% The way I select the reference transition isn't necessarily better than yours.
% At the moment the choice is random.  The case where there are only two subspaces
% and a single transition is actually trivial -- just compute the diophantine approx
% of the single transition frequency.  The way the subconvergents are handled could
% be improved.

if ~exist('sfact','var')
    sfact = [1 1e-10];
end
if ~exist('runs','var')
    runs   = 500;
end

[P,e] = EigDecomp(obj);
Ns = length(e);

% calculate input/output state projection onto eigenspaces
Pij = cell2mat (arrayfun (@(k) P{k}(ispin,jspin), (1:Ns)', 'UniformOutput', false));

% eliminate dark subspaces
s  = sign(Pij);
nonzero_proj = find(s);
ee = e(nonzero_proj);
ss = s(nonzero_proj);
Pij= Pij(nonzero_proj);

if length(ss)<3,
    error('we have a problem')
end

% sort subspaces
[sss,ind] = sort(ss,1,'descend');
ees  = ee(ind);
Pijs = Pij(ind);
omg  = diff(ees);
par  = mod(diff(sss)/2,2);

% select reference transition
if sss(2)-sss(1)==0,
     % first transition
    omg0  = omg(1);
    theta = 2*omg(2:end)/omg0;
    parity= par(2:end);
else % last transition
    omg0  = omg(end);
    theta = 2*omg(1:end-1)/omg0;
    parity= par(1:end-1);
end

done = 0; k = 1;
while ~done
    x = [rand(1,length(theta))*sfact(1), sfact(2)];
   [p(:,k),q(k)] = sim_dio_approx(theta,x,1);
   % ensure q is positive
   p(:,k) = p(:,k)*sign(q(k));
   q(k)   = q(k)*sign(q(k));
   % check parity
   parity_ok = all(mod(p(:,k),2) == parity);
   done      = (q(k)>0 & parity_ok) || k == runs;
   % only increment if q~=0 and not done
   k = k + (q(k)~=0)*(~done);
end

out.omg0  = omg0;
out.theta = theta;

if parity_ok % terminated with correct parity
    out.p = p(:,end);
    out.q = q(end);
    out.tf    = 2*q(end)*pi/omg0;
    out.err_dio_approx = norm(p(:,end)/q(end)-theta);
    max_p = obj.prob();
    act_p = abs(expm(-i*out.tf*obj.H)).^2;
    out.err_act = (max_p(ispin,jspin) - act_p(ispin,jspin));
    out.err_act_rel = out.err_act/max_p(ispin,jspin);
else
    p2 = []; q2 = [];
    for k = 1:runs
	for l = k+1:runs
	    parity_ok = all(mod(p(:,k)+p(:,l),2) == parity);
    	    if parity_ok
		p2 = p(:,k)+p(:,l);
        	q2 = q(:,k)+q(:,l);
	    end
	end
    end
    out.p2 = p2;
    out.q2 = q2;
    out.tf = 2*q2*pi/omg0;
    max_p = obj.prob();
    for k=1:length(q2)
        out.err_dio_approx(k) = norm(p2(:,k)/q2(k) -theta);
        act_p = abs(expm(-i*out.tf(k)*obj.H)).^2;
        out.err_act(k) = (max_p(ispin,jspin) - act_p(ispin,jspin));
	out.err_act_rel(k) = out.err_act(k)/max_p(ispin,jspin);
    end
end

