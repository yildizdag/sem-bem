%% B-SPLINE BASIS FUNCTIONS
% TÜBİTAK 3501 - İP-1.1
% Variation of B-spline basis functions along the parametric coordinate

clear;
clc;
close all;

%% ===============================================================
%  PARAMETERS
% ===============================================================

% Polynomial degree
p = 3;

% Open non-uniform knot vector
Xi = [0 0 0 0 ...
      0.20 0.40 0.60 0.80 ...
      1 1 1 1];

% Number of B-spline basis functions
nBasis = length(Xi) - p - 1;

% Parametric coordinate for plotting
xi = linspace(Xi(1), Xi(end), 1500);


%% ===============================================================
%  EVALUATE B-SPLINE BASIS FUNCTIONS
% ===============================================================

N = zeros(nBasis, length(xi));

for i = 1:nBasis
    for k = 1:length(xi)
        N(i,k) = BSplineBasis(i, p, xi(k), Xi);
    end
end


%% ===============================================================
%  PLOT
% ===============================================================

figure('Color','w', ...
       'Position',[100 100 800 480]);

hold on;
box on;

% ---------------------------------------------------------------
% Plot all B-spline basis functions
% ---------------------------------------------------------------

for i = 1:nBasis
    plot(xi, N(i,:), ...
        'LineWidth', 1.8,...
        'Color','k');
end


% ---------------------------------------------------------------
% Show internal knots
% ---------------------------------------------------------------

internalKnots = unique(Xi(p+2:end-p-1));

% for k = 1:length(internalKnots)
% 
%     xline(internalKnots(k), ':', ...
%         'LineWidth', 0.001);
% 
% end


% ---------------------------------------------------------------
% Add labels directly to the basis functions
% ---------------------------------------------------------------

for i = 1:nBasis

    % Maximum value and corresponding location
    [Nmax, idx] = max(N(i,:));

    xText = xi(idx);
    yText = Nmax;

    % Slightly different offsets to avoid overlapping labels
    if i == 1
        yOffset = -0.06;
        xOffset = 0.045;
    elseif i == 8
        yOffset = -0.06;
        xOffset = -0.045;
    else
        yOffset = 0.03;
        xOffset = 0;
    end

    text(xText + xOffset, yText + yOffset, ...
        sprintf('$N_{%d,%d}$', i, p), ...
        'Interpreter','latex', ...
        'FontSize',18, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom');

end


% ---------------------------------------------------------------
% Axis labels
% ---------------------------------------------------------------

xlabel('$\xi$', ...
    'Interpreter','latex', ...
    'FontSize',18);

ylabel('$N_{i,p}(\xi)$', ...
    'Interpreter','latex', ...
    'FontSize',18);


% ---------------------------------------------------------------
% Axis limits
% ---------------------------------------------------------------

xlim([Xi(1) Xi(end)]);
ylim([0 1.12]);


% ---------------------------------------------------------------
% Axis appearance
% ---------------------------------------------------------------

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',18, ...
    'LineWidth',1, ...
    'TickLabelInterpreter','latex');

grid off;


%% ===============================================================
%  OPTIONAL EXPORT
% ===============================================================

% Vector PDF:
% exportgraphics(gcf, ...
%     'BSpline_Basis_Functions.pdf', ...
%     'ContentType','vector');

% High-resolution PNG:
% exportgraphics(gcf, ...
%     'BSpline_Basis_Functions.png', ...
%     'Resolution',600);


%% ===============================================================
%  B-SPLINE BASIS FUNCTION
%  Cox-de Boor recursive definition
% ===============================================================

function N = BSplineBasis(i, p, xi, Xi)

    % -----------------------------------------------------------
    % Zeroth-degree basis function
    % -----------------------------------------------------------

    if p == 0

        if Xi(i) <= xi && xi < Xi(i+1)

            N = 1;

        % Special treatment at the right end of the parameter space
        elseif xi == Xi(end) && ...
               Xi(i+1) == Xi(end) && ...
               Xi(i) < Xi(end)

            N = 1;

        else

            N = 0;

        end

        return;

    end


    % -----------------------------------------------------------
    % First term of Cox-de Boor recursion
    % -----------------------------------------------------------

    denominator1 = Xi(i+p) - Xi(i);

    if denominator1 == 0

        term1 = 0;

    else

        term1 = ...
            (xi - Xi(i)) / denominator1 * ...
            BSplineBasis(i, p-1, xi, Xi);

    end


    % -----------------------------------------------------------
    % Second term of Cox-de Boor recursion
    % -----------------------------------------------------------

    denominator2 = Xi(i+p+1) - Xi(i+1);

    if denominator2 == 0

        term2 = 0;

    else

        term2 = ...
            (Xi(i+p+1) - xi) / denominator2 * ...
            BSplineBasis(i+1, p-1, xi, Xi);

    end


    % -----------------------------------------------------------
    % Basis function value
    % -----------------------------------------------------------

    N = term1 + term2;

end