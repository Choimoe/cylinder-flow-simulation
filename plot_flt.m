hold on
axis equal
X = load("X.txt");
Y = load("Y.txt");
M = 101;
N = 101;

for i = 1 : M
    plot(X(i, :), Y(i, :), '-b');
end

for j = 2 : N - 1
    plot(X(:, j), Y(:, j), '-b');
end


r = 1; theta = 0 : pi / 100 : 2 * pi;
x = r * cos(theta); y = r * sin(theta);
plot(x, y, '-r')
r = 20; 
theta = 0 : pi / 100 : 2 * pi;
x = r * cos(theta); y = r * sin(theta);
plot(x, y, '-r')

hold on