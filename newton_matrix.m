function x = newton_matrix(x,F,J,N)

k = 1;

while(k <= N)
    y = -J(x)\F(x);
    x = x + y;
    k = k + 1;
end
