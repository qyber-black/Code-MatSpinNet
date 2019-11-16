function [ ind ] = reverse2_generate_jayaraman( problem, ~ )
%GERAINDIVIDUO

% Maximum quantity of clients i
maxIteracao = 100;

%armazena as melhores solu��es
solFinalOtima = Inf;

poolSolucoes = inf(1,round(maxIteracao*0.05));
poolEscolhidosJOtimo = zeros(round(maxIteracao*0.05),problem.n_j);
poolEscolhidosKOtimo = zeros(round(maxIteracao*0.05),problem.n_k);

for iteracao=1:maxIteracao,
    clear escolhidosJ;
    clear escolhidosK;
    
    %seleciona uma quantidade de dep�sitos maxJ aleatoriamente, criando o vetor
    %escolhidosJ que cont�m os maxJ escolhidos (marcado como 1)
    i=0;
    escolhidosJ = zeros(1,problem.n_j);
    while(1)
        pos=mod(round(100*rand()),problem.n_j)+1;
        if escolhidosJ(1,pos)==0
            escolhidosJ(1,pos)=1;
            i=i+1;
        else
            continue;
        end
        if(i==problem.p_max)
            break;
        end
    end
    
    %seleciona uma quantidade de fabricas maxJ aleatoriamente, criando o vetor
    %escolhidosK que cont�m os maxK escolhidos (marcado como 1)
    i=0;
    escolhidosK = zeros(1,problem.n_k);
    while(1)
        pos=mod(round(100*rand()),problem.n_k)+1;
        if escolhidosK(1,pos)==0
            escolhidosK(1,pos)=1;
            i=i+1;
        else
            continue
        end
        if(i==problem.q_max)
            break;
        end
    end
    
    ind.p = escolhidosJ;
    ind.q = escolhidosK;
    
    [fx,ind] = reverse2_evaluate(ind, problem);
    
    escolhidosJ = ind.p;
    escolhidosK = ind.q;
    
    %coloca solu��o no pool de solucoes
    [valor, index] = max(poolSolucoes);
    if (fx < valor)
        poolSolucoes(index) = fx;
        poolEscolhidosJOtimo(index,:) = escolhidosJ;
        poolEscolhidosKOtimo(index,:) = escolhidosK;
    end
    
    if (solFinalOtima > fx)
        solFinalOtima = fx;
        escolhidosJOtimo = escolhidosJ;
        escolhidosKOtimo = escolhidosK;
    end
end


%iniciando heuristic concentration. nesse caso, pegamos a melhor solu��o
%encontrada e adicionamos P_max + 2 - n�meros de locais na solu��o. Como o
%n�mero de locais na solu��o � igual a P_max, vamos adicionar 2 locais k e j ao
%conjunto de solu��es, rodar novalmente o simplex e comparar com a melhor solu��o
%obtida. Os locais adicionados s�o os dois que mais apareceram nas outras
%solu��es do pool.

%fprintf('inicia heuristic concentration\n');
for it=1:size(poolEscolhidosJOtimo,1)
poolEscolhidosJOtimo(it,:) =(poolEscolhidosJOtimo(it,:) - escolhidosJOtimo);
end
[~,index] = sort(sum(poolEscolhidosJOtimo,1),'descend');
escolhidosJOtimoHC = escolhidosJOtimo;
escolhidosJOtimoHC(index(1)) = 1; escolhidosJOtimoHC(index(2)) = 1;

for it=1:size(poolEscolhidosKOtimo,1)
poolEscolhidosKOtimo(it,:) = (poolEscolhidosKOtimo(it,:) - escolhidosKOtimo);
end
[~,index] = sort(sum(poolEscolhidosKOtimo,1),'descend');
escolhidosKOtimoHC = escolhidosKOtimo;
escolhidosKOtimoHC(index(1)) = 1; escolhidosKOtimoHC(index(2)) = 1;

indHC.p = escolhidosJOtimoHC;
indHC.q = escolhidosKOtimoHC;

[fx,indHC] = reverse2_evaluate(indHC, problem);


    %verifica se a solu��o respeita maxJ e maxD
    if ( size(indHC.p(indHC.p~=0),2) > (problem.p_max) || size(indHC.q(indHC.q~=0),2) > (problem.q_max)) 
        solFinal = Inf;
    else
        solFinal = fx;
   end

if (solFinalOtima > solFinal)
    %fprintf('Melhor solucao encontrada com a HC!\n');
    ind.p = escolhidosJOtimoHC;
    ind.q = escolhidosKOtimoHC;
    ind.x = indHC.x;
end

%inicia expans�o heuristica
modifica = 1;
while (modifica > 0)
    melhor = zeros(3,1);
    melhorF = inf;    
    modifica = 0;
    
    for i=1:problem.n_j
        if (ind.p(i) == 0)
            novo.p = ind.p;
            novo.q = ind.q;
            novo.p(i) = 1;
            [fx,novo] = reverse2_evaluate(novo, problem); 
            if(solFinalOtima > fx && melhorF > fx) %se a solu��o encontrada for melhor que todas as solu��es j� encontradas
                %fprintf('Melhor solucao encontrada com a HE! Testando se � poss�vel coloc�-la\n');
                novo.p;
                if (size(novo.p(novo.p~=0),2) <= (problem.p_max) && size(novo.q(novo.q~=0),2) <= (problem.q_max)) %se for uma solu��o v�lida
                    melhor = novo;
                    melhorF = fx;
                    modifica = 1;
                end
            end
        end
    end
    
    for i=1:problem.n_k
        if (ind.q(i) == 0)
            novo.p = ind.p;
            novo.q = ind.q;
            novo.q(i) = 1;
            [fx,novo] = reverse2_evaluate(novo, problem); 
            if(solFinalOtima > fx && melhorF > fx) %se a solu��o encontrada for melhor que todas as solu��es j� encontradas
                if (size(novo.p(novo.p~=0),2) <= (problem.p_max) && size(novo.q(novo.q~=0),2) <= (problem.q_max)) %se for uma solu��o v�lida
                    melhor = novo;
                    melhorF = fx;
                    modifica = 1;
                end
            end
        end
    end    

    
    %se a solu��o foi modificada
    if (modifica > 0)
        %fprintf('solu��o melhor encontrada pela expans�o heuristica\n');
        solFinalOtima = melhorF;
        ind = melhor;
    end
%toc;
end

