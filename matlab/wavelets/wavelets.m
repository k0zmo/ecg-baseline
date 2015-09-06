%% Pobranie danych
ecg_data = rdsamp('mitdb/101', 'maxt', ':40');
t = ecg_data(:,1);
x = ecg_data(:,2);

%% Transformata prosta
level = 7;
y = x(:)'; % wymus postac [a b c] (row vector)
c = [];
l = zeros(1, level+2); % 1 dodatkowo na caly sygnal i 1 na aproksymacje
l(end) = length(y);

for i = 1:level
    [y,d] = dwt1(y,'db4');
    c = [d c]; % dopisz info o detalach
    l(level+2-i) = length(d); % zapisz info o dlugosci detali
end

% Zapisz ostatnia aproksymacje
c = [y c];
l(1) = length(y);

%% Filtracja
% usun wspolczynniki z ostatniej aproksymacji
c(1:l(1)) = zeros(l(1),1);

%% Transformata odwrotna
a = c(1:l(1));
to = l(1);

for i = 1:level
    from = to + 1;
    to = from + l(i+1) - 1;
    a = idwt1(a, c(from:to), 'db4', l(i+2));
    size(a)
end

plot(linspace(0,1,length(a)),a,'linewidth',1), hold on
plot(linspace(0,1,length(x)),x,'r');
