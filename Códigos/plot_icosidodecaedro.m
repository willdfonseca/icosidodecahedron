function [V, F] = plot_icosidodecaedro(value, mode, circleRadius)
% plot_icosidodecaedro(value, mode, circleRadius)
%
% Plota o icosidodecaedro (sólido de Arquimedes com 20 faces triangulares e 12 pentagonais)
% de acordo com um parâmetro geométrico e, ao final, exibe:
%   - o comprimento da aresta (dos pentágonos, que é igual ao das outras faces),
%   - a área de um pentágono regular com essa aresta,
%   - o diâmetro aproximado (distância máxima entre quaisquer dois vértices),
%   - o raio do incírculo (maior círculo inscrito) em um pentágono,
%   - a área de um triângulo equilátero (face triangular) e
%   - a área total do icosidodecaedro.
%
% Parâmetros:
%   value : valor numérico que define o tamanho do sólido.
%           No modo 'edge', value é o comprimento da aresta.
%           No modo 'diameter', value é o diâmetro aproximado.
%
%   mode  : string que pode ser 'edge' (padrão) ou 'diameter'.
%
%   circleRadius (opcional): se fornecido (e não vazio), em cada face pentagonal
%           será plotado um círculo de raio especificado, em cor cinza.
%
% Exemplos:
%   plot_icosidodecaedro(1);                % assume aresta = 1 (modo 'edge')
%   plot_icosidodecaedro(2, 'diameter');      % diâmetro aproximado = 2
%   plot_icosidodecaedro(1, 'edge', 0.05);    % plota também círculos de raio 0.05
%
% Autor:
% Data: 15/02/2025

if nargin < 1
    value = 1;
end
if nargin < 2
    mode = 'edge';
end
if nargin < 3
    circleRadius = [];
end

%% 1. Definindo o icosaedro padrão

phi = (1+sqrt(5))/2;  % razão áurea

Vicos = [...
    -1,  phi,  0;
     1,  phi,  0;
    -1, -phi,  0;
     1, -phi,  0;
     0, -1,  phi;
     0,  1,  phi;
     0, -1, -phi;
     0,  1, -phi;
     phi,  0, -1;
     phi,  0,  1;
    -phi,  0, -1;
    -phi,  0,  1];
 
Ficos = [...
    1, 12, 6;
    1, 6, 2;
    1, 2, 8;
    1, 8, 11;
    1, 11, 12;
    2, 6, 10;
    6, 12, 5;
    12, 11, 3;
    11, 8, 7;
    8, 2, 9;
    4, 10, 5;
    4, 5, 3;
    4, 3, 7;
    4, 7, 9;
    4, 9, 10;
    5, 10, 6;
    3, 5, 12;
    7, 3, 11;
    9, 7, 8;
    10, 9, 2];

%% 2. Obter as arestas únicas e calcular os pontos médios
% Cada face triangular gera 3 arestas
arestas = [];
for i = 1:size(Ficos,1)
    face = Ficos(i,:);
    e1 = sort(face([1 2]));
    e2 = sort(face([2 3]));
    e3 = sort(face([3 1]));
    arestas = [arestas; e1; e2; e3]; %#ok<AGROW>
end
arestas = unique(arestas, 'rows');
numArestas = size(arestas,1);

% Cada aresta gera um vértice do icosidodecaedro (ponto médio)
Vrect = zeros(numArestas, 3);
edgeMap = containers.Map('KeyType','char','ValueType','int32');
for i = 1:numArestas
    idx = arestas(i,:);
    mid = ( Vicos(idx(1),:) + Vicos(idx(2),:) )/2;
    Vrect(i,:) = mid;
    key = sprintf('%d_%d', idx(1), idx(2));
    edgeMap(key) = i;
end

%% 3. Reescalando o sólido de acordo com o modo escolhido

switch lower(mode)
    case 'edge'
        % Usa uma aresta de uma face triangular para determinar o comprimento atual
        faceEx = Ficos(1,:);
        e1 = sort(faceEx([1 2]));
        e2 = sort(faceEx([2 3]));
        key1 = sprintf('%d_%d', e1(1), e1(2));
        key2 = sprintf('%d_%d', e2(1), e2(2));
        v1 = Vrect(edgeMap(key1),:);
        v2 = Vrect(edgeMap(key2),:);
        L_atual = norm(v1 - v2);
        escala = value / L_atual;
    case 'diameter'
        % Diâmetro aproximado: distância máxima entre quaisquer dois vértices
        D_atual = max(pdist(Vrect));
        escala = value / D_atual;
    otherwise
        error('Modo desconhecido. Use "edge" ou "diameter".');
end

Vrect = Vrect * escala;
Vicos = Vicos * escala;

%% 4. Construindo as faces do icosidodecaedro

% (a) Faces triangulares: cada face do icosaedro gera uma face triangular
triFaces = zeros(size(Ficos,1),3);
for i = 1:size(Ficos,1)
    face = Ficos(i,:);
    e1 = sort(face([1 2]));
    e2 = sort(face([2 3]));
    e3 = sort(face([3 1]));
    key1 = sprintf('%d_%d', e1(1), e1(2));
    key2 = sprintf('%d_%d', e2(1), e2(2));
    key3 = sprintf('%d_%d', e3(1), e3(2));
    triFaces(i,:) = [ edgeMap(key1), edgeMap(key2), edgeMap(key3) ];
end

% (b) Faces pentagonais: cada vértice do icosaedro gera um pentágono
numPent = size(Vicos,1);  % 12 vértices => 12 pentágonos
pentFaces = cell(numPent,1);
for v = 1:numPent
    inds = [];
    for k = 1:size(arestas,1)
        if any(arestas(k,:) == v)
            inds(end+1) = k; %#ok<AGROW>
        end
    end
    if isempty(inds)
        warning('Vértice %d não foi encontrado em nenhuma aresta!', v);
        continue;
    end
    
    pts = Vrect(inds,:);
    centro = mean(pts,1);
    ref = pts(1,:) - centro;
    n = Vicos(v,:) / norm(Vicos(v,:));
    b2 = cross(n, ref);
    if norm(b2) < 1e-6
        B = null(ref');
        if isempty(B)
            warning('null(ref'') retornou vazio para o vértice %d.', v);
            continue;
        else
            b2 = B(:,1)';
        end
    end
    angles = zeros(length(inds),1);
    for j = 1:length(inds)
        d = pts(j,:) - centro;
        angles(j) = atan2(dot(d,b2), dot(d,ref));
    end
    [~, ordem] = sort(angles);
    pentFaces{v} = inds(ordem);
end

%% 5. Plotando o sólido

figure;
hold on;
% Faces triangulares
for i = 1:size(triFaces,1)
    patch('Vertices', Vrect, 'Faces', triFaces(i,:), ...
          'FaceColor', [1 0.8 0.8], 'EdgeColor', 'k', 'FaceAlpha', 0.9);
end
% Faces pentagonais
for i = 1:length(pentFaces)
    if ~isempty(pentFaces{i})
        patch('Vertices', Vrect, 'Faces', pentFaces{i}, ...
              'FaceColor', [0.8 1 0.8], 'EdgeColor', 'k', 'FaceAlpha', 0.9);
    end
end

% Se o parâmetro circleRadius foi informado, plota um círculo em cada face pentagonal
if ~isempty(circleRadius)
    for i = 1:length(pentFaces)
        if ~isempty(pentFaces{i}) && numel(pentFaces{i}) >= 3
            indices = pentFaces{i};
            V_face = Vrect(indices, :);
            center_face = mean(V_face,1);
            % Calcula o normal da face usando os três primeiros vértices
            n_face = cross(V_face(2,:) - V_face(1,:), V_face(3,:) - V_face(1,:));
            if norm(n_face) < eps
                n_face = [0,0,1];
            else
                n_face = n_face / norm(n_face);
            end
            % Vetor de referência no plano da face: do centro até o 1º vértice
            r_face = V_face(1,:) - center_face;
            if norm(r_face) < eps
                r_face = [1,0,0];
            else
                r_face = r_face / norm(r_face);
            end
            % Vetor perpendicular no plano
            b_face = cross(n_face, r_face);
            if norm(b_face) < eps
                b_face = [0,1,0];
            else
                b_face = b_face / norm(b_face);
            end
            theta = linspace(0, 2*pi, 100);
            circle_pts = center_face + circleRadius * (cos(theta)'*r_face + sin(theta)'*b_face);
            plot3(circle_pts(:,1), circle_pts(:,2), circle_pts(:,3), 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
        end
    end
end

axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
if strcmpi(mode,'edge')
    title(['Icosidodecaedro com aresta = ' num2str(value)]);
else
    title(['Icosidodecaedro com diâmetro \approx ' num2str(value)]);
end
grid on;
view(3);
hold off;

%% 6. Cálculos finais e exibição dos resultados

% Seleção de uma aresta a partir da primeira face triangular
edge_length = norm(Vrect(triFaces(1,1),:) - Vrect(triFaces(1,2),:));

% Área de um pentágono regular com aresta "edge_length"
area_pentagon = (1/4)*sqrt(5*(5+2*sqrt(5)))*edge_length^2;

% Área de um triângulo equilátero com aresta "edge_length"
area_triangle = (sqrt(3)/4)*edge_length^2;

% Diâmetro aproximado: maior distância entre quaisquer dois vértices
D_approx = max(pdist(Vrect));

% Raio do incírculo em um pentágono regular:
r_incircle = edge_length/(2*tan(pi/5));

% Área total do icosidodecaedro:
area_total = 20*area_triangle + 12*area_pentagon;

fprintf('\nResultados:\n');
fprintf('Comprimento da aresta (dos pentágonos): %.4f\n', edge_length);
fprintf('Área de um pentágono: %.4f\n', area_pentagon);
fprintf('Área de um triângulo: %.4f\n', area_triangle);
fprintf('Área total do icosidodecaedro: %.4f\n', area_total);
fprintf('Diâmetro aproximado do sólido: %.4f\n', D_approx);
fprintf('Raio do círculo inscrito no pentágono: %.4f\n', r_incircle);

%% 7. Retornando os dados da malha
% V são os vértices utilizados na plotagem (Vrect)
% F é uma estrutura contendo as faces triangulares e pentagonais
V = Vrect;
F.tri = triFaces;
F.pent = pentFaces;

end
