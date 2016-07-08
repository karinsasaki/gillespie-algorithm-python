function [ ht, hn ] = francois_gillespie(stoch, rates, initial, tmax)

% function [ ht, hn ] = gillespie(stoch, rates, initial, tmax)
%
% Simulate a system of reactions defined in a stoichiometry matrix for a time `tmax`,
% with a set of rates and the initial abundance of each species.
% The system of reaction is simulated stochastically with Gillespie's method.
%
% `stoch` is a stoichiometry matrix:
%    its number of lines should be the number of reaction
%    its number of column should be the number of species
% the vector `rates` should provide the rate of each reaction specified in `stoch`
% the vector `initial` specifies the initial aboundance of each species.
%
% The output is a vector of time `ht`, and a matrix `hn` containing the number of each species,
% for every time point.
%
%
% Example: to simulate a simple decay:
%   A -> B  with rate kA
%
% [ ht, hn ] = gillespie([-1], [kA], [nA], 100);
% plot(ht, hn);
%
% Example: to simulate an equilibrium:
%   A -> B  with rate kA
%   B -> A  with rate kB
%
% [ ht, hn ] = gillespie([ -1, 1; 1, -1 ], [kA, kB], [nA, nB], 100);
% plot(ht, hn);
%
%
% Copyright F. Nedelec, 09.2014




[ nreac, nspec ] = size(stoch);

if length(rates) ~= nreac
    error('the rate vector does not match the number of reactions');
end
rates = reshape(rates, nreac, 1);

if length(initial) ~= nspec
    error('the initial state vector does not match the number of species');
end

t = 0;
n = reshape(initial, 1, nspec);
vo = ones(nreac, 1);
stoch_in = -stoch .* ( stoch < 0 );

est = ceil( 2 * sum(initial) * sum(rates) * tmax );
ht = zeros(est, 1);
hn = zeros(est, nspec);


tix = 1;
ht(tix) = t;
hn(tix, :) = n;

while t < tmax;

    pp = prod(power(vo * n, stoch_in), 2) .* rates;
    [dt, ix] = gillespie_choice(pp);
    t = t + dt;
    n = n + stoch(ix, :);

    tix = tix + 1;
    ht(tix) = t;
    hn(tix, :) = n;
    
end

ht = ht(1:tix);
hn = hn(1:tix, :);

end