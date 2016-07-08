function [ T, N ] = gillespie_choice(rates)

% set the time delay T, and the index N of the next reaction
% according to Gillespie's procedure.
% The input is a vector of reaction rates.
%
% Two random numbers are used. 
% F. Nedelec, 2012

S = sum(rates);

T = -log(rand) / S;

N = 1;
X = rand * S - rates(1);

while X > 0
    N = N + 1;
    X = X - rates(N);
end


end