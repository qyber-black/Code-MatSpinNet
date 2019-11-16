function solution = reverse2_mutation_random(solution, problem,~,~)
%MUTACAO REVERSE2

pos = randi(problem.n_j);
% se j� estiver no m�ximo
if (sum(solution.p)==problem.p_max)
    %faz muta��o de redu��o na posi��o
    solution.p(pos) = max(solution.p(pos)-1,0);
% se j� estiver no m�nimo
elseif (sum(solution.p)==problem.p_min)
    %faz muta��o de acr�scimo na posi��o
    solution.p(pos) = min(solution.p(pos)+1,1);
% se estiver no meio termo
else
    %faz um bit-flop simples
    solution.p(pos) = mod(solution.p(pos)-1,2);
end

pos = randi(problem.n_k);
% se j� estiver no m�ximo
if (sum(solution.q)==problem.q_max)
    %faz muta��o de redu��o na posi��o
    solution.q(pos) = max(solution.q(pos)-1,0);
% se j� estiver no m�nimo
elseif (sum(solution.q)==problem.q_min)
    %faz muta��o de acr�scimo na posi��o
    solution.q(pos) = min(solution.q(pos)+1,1);
% se estiver no meio termo
else
    %faz um bit-flop simples
    solution.q(pos) = mod(solution.q(pos)-1,2);
end


end

