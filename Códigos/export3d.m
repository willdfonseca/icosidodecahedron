% stl export


%% Cleaning Service
clear all; clc; close all


%% Gera a malha do icosidodecaedro e exporta para STL
[V, F] = plot_icosidodecaedro(1, 'edge', []);  % Exemplo sem círculos

% Inicia a lista de faces trianguladas com as faces já trianguladas
facesTriAll = F.tri;  

% Triangula as faces pentagonais (fan triangulation)
for i = 1:length(F.pent)
    poly = F.pent{i};
    if ~isempty(poly) && numel(poly) >= 3
        % Fan triangulation: fixa o primeiro vértice e cria triângulos
        for j = 2:(length(poly)-1)
            facesTriAll(end+1,:) = [poly(1), poly(j), poly(j+1)];
        end
    end
end

% Cria o objeto triangulation
TR = triangulation(facesTriAll, V);

% Salva em STL
stlwrite(TR, 'icosidodecaedro_1.stl');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gera a malha do icosaedro

clear all; clc; close all

[V, F] = plot_polyhedron(1, 'edge', 'icosaedro');

% Converte para objeto triangulation
TR = triangulation(F, V);

% Salva em STL
stlwrite(TR, 'icosaedro.stl');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gera a malha do dodecaedro

clear all; clc; close all

[V, Fcell] = plot_polyhedron(2, 'edge', 'dodecaedro');

% Converte as faces (polígonos) em triângulos
facesTri = [];
for i = 1:length(Fcell)
    if ~isempty(Fcell{i})
        poly = Fcell{i};
        % Fan triangulation: fixa o primeiro vértice e cria triângulos
        for j = 2:(length(poly)-1)
            facesTri(end+1,:) = [poly(1), poly(j), poly(j+1)]; 
        end
    end
end

% Converte em objeto triangulation
TR = triangulation(facesTri, V);

% Salva em STL
stlwrite(TR, 'dodecaedro.stl');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gera a malha do icosidodecaedro
clear all; clc; close all

[V, Fstruct] = plot_polyhedron(1, 'edge', 'icosidodecaedro');

% Inicia a lista de faces trianguladas com as faces triangulares já existentes
facesTri = Fstruct.tri;  

% Triangula as faces pentagonais (fan triangulation)
for i = 1:length(Fstruct.pent)
    poly = Fstruct.pent{i};
    if ~isempty(poly)
        for j = 2:(length(poly)-1)
            facesTri(end+1,:) = [poly(1), poly(j), poly(j+1)]; 
        end
    end
end

% Converte em objeto triangulation
TR = triangulation(facesTri, V);

% Salva em STL
stlwrite(TR, 'icosidodecaedro.stl');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%