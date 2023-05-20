data_path = "";
data_time = "3000";

PSI3000 = readmatrix([data_path, "PSI", data_time, ".dat"]);
U3000 = readmatrix([data_path, "U", data_time, ".dat"]);
V3000 = readmatrix([data_path, "V", data_time, ".dat"]);

X = readmatrix([data_path, "X.txt"]);
Y = readmatrix([data_path, "Y.txt"]);

figure

quiver(X, Y, U3000, V3000, 'r');

M = sqrt(U3000.^2+V3000.^2);
colorlist = jet;
Mdown = min(M(:));
Mup = max(M(:));
Mlist = linspace(Mdown, Mup, 256);

scaler1 = Mup./M;
U3000 = U3000.*scaler1;
V3000 = V3000.*scaler1;

scaler2 = 0.05;
U3000 = U3000*scaler2;
V3000 = V3000*scaler2;

figure

[m, n] = size(X);
for i = 1 : m
    for j = 1 : n
        Mtemp = abs(M(i, j) - Mlist);
        index = find(Mtemp == min(Mtemp));
        colorarrow = colorlist(index,:);
        q = quiver(X(i, j), Y(i, j), U3000(i, j), V3000(i, j), 'MaxHeadSize', 100);
        q.LineWidth = 1;
        q.Color = colorarrow;
        hold on
    end
end

hc = colorbar;

colormap(jet)

hc.TickLabels = linspace(Mdown, Mup, 11);