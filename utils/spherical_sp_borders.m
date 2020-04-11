

function B = spherical_sp_borders(S)

[h,w] = size(S);

B = zeros(h,w);
for i=1:h-1
    for j=1:w-1
        label = S(i,j);
        if (label ~= S(i,j+1))
            B(i,j) = 1;
        end
        if (label ~= S(i+1,j))
            B(i,j) = 1;
        end
        if (label ~= S(i+1,j+1))
            B(i,j) = 1;
        end     
    end
end
j = w;
for i=1:h
    if (B(i,j-1))
        B(i,j) = 1;
    end
end