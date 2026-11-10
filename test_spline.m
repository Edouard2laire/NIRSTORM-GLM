% 1. Paramètres des B-splines
t_start = 0;
t_end = 25;
num_basis = 10;   % Nombre de fonctions de base souhaité
order = 4;        % Ordre 4 = Spline Cubique (de degré 3)

% 2. Calcul du nombre de nœuds internes et du vecteur de nœuds
% Formule : Nombre de nœuds = Nombre de bases + Ordre
num_knots = num_basis + order;
num_inner_knots = num_basis - order; 

% Génération des nœuds internes équidistants
inner_knots = linspace(t_start, t_end, num_inner_knots + 2);
inner_knots = inner_knots(2:end-1); % On garde uniquement l'intérieur

% Construction du vecteur de nœuds clampé (Knot Vector)
knots = [repmat(t_start, 1, order), inner_knots, repmat(t_end, 1, order)];

% 3. Grille d'évaluation fine de 0 à 25s
t_eval = linspace(t_start, t_end, 500);

% 4. Calcul de la matrice des fonctions de base
% Chaque colonne de 'B' représentera une fonction de base
B = spcol(knots, order, t_eval);

% 5. Affichage des 10 fonctions de base
figure;
plot(t_eval, B, 'LineWidth', 2);
grid on;
xlim([t_start, t_end]);
ylim([-0.05, 1.05]);
xlabel('Temps (s)', 'FontSize', 12);
ylabel('Amplitude B_i(t)', 'FontSize', 12);
title(sprintf('%d Fonctions de Base B-splines Cubiques (0 à 25s)', num_basis), 'FontSize', 14);

% Ajout d'une légende propre
labels = arrayfun(@(x) sprintf('B_{%d}(t)', x), 1:num_basis, 'UniformOutput', false);
legend(labels, 'Location', 'eastoutside');
