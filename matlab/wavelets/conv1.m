function [y] = conv1(f,g,shape)
lf = length(f);
lg = length(g);

if ~strcmp(shape,'full')
    ly = max(lf-max(0,lg-1),0);
    y = zeros(1,ly);

    for k=1:ly
        for j=1:lg
            y(k) = y(k) + f(k-j+lg) * g(j);
        end
    end
else
    ly = lf+lg-1;
    y = zeros(1,ly);

    for k=1:ly
        for j=1:lg
            if k-j+1 >= 1 && k-j+1 <= lf 
                y(k) = y(k) + f(k-j+1) * g(j);
            end
        end
    end
end