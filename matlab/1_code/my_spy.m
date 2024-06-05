function my_spy(S)
    %S = bucky();
    %S(S~=0) = log10(abs(S(S~=0)));
    %S(S < -20) = -20 * ones(size(S(S < -20)));
    %S(S > 0) = log10(S(S>0));
    %S(S < 0) = log10(-S(S<0));
    [m, n] = size(S);
    %[X, Y] = meshgrid(1:m, 1:n);
    %S = (X + Y) .* S;
    nonzeroInd = find(S);
    [x, y] = ind2sub([m n], nonzeroInd);

    figure();
    patch(y, x, log10(abs(S(nonzeroInd))), ...
               'Marker', 's', 'MarkerFaceColor', 'flat', 'MarkerSize', 3, ...
               'EdgeColor', 'none', 'FaceColor', 'none');
    set(gca, 'XLim', [0, n + 1], 'YLim', [0, m + 1], 'YDir', 'reverse', ...
        'PlotBoxAspectRatio', [n + 1, m + 1, 1]);
    colormap(jet(20));
    colorbar();
end