function [ fx ] = realfunction_evaluate_foxholes( x )
%FOXHOLES funcao

a = [-32, -16, 0, 16, 32, -32, -16, 0 ,16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32;
	 -32, -32, -32, -32, -32, -16, -16, -16, -16, -16, 0, 0, 0, 0, 0, 16, 16, 16, 16, 16, 32, 32, 32, 32, 32];

sum = 0;
for i=1:25
    sum = sum + 1/(i + (x(1)-a(1,i))^6 + (x(2)-a(2,i))^6 );
end


fx = (1/(1/500+sum))-0.9980038378;


end

