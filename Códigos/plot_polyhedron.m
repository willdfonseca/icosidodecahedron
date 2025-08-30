function [V, F] = plot_polyhedron(value, mode, solid, circleRadius)
% plot_polyhedron(value, mode, solid, circleRadius)
%
% Plota um dos sólidos platônicos/arquimédicos a partir de um parâmetro geométrico:
%   - 'icosaedro'         : 20 faces triangulares (sólido dual do dodecaedro);
%   - 'dodecaedro'        : 12 faces pentagonais (dual do icosaedro);
%   - 'icosidodecaedro'   : 20 faces triangulares e 12 faces pentagonais.
%
% Parâmetros:
%   value : valor numérico que define o tamanho do sólido.
%           No modo 'edge', value é o comprimento da aresta.
%           No modo 'diameter', value é o diâmetro aproximado.
%
%   mode  : string que pode ser 'edge' (padrão) ou 'diameter'.
%
%   solid : string que pode ser 'icosaedro', 'dodecaedro' ou 'icosidodecaedro'
%           (padrão: 'icosidodecaedro').
%
%   circleRadius (opcional): se fornecido e não vazio, para os sólidos com faces
%           (triangulares ou pentagonais) será plotado, em cada face, um círculo de
%           raio especificado, em cor cinza.
%
% Exemplos:
%   plot_polyhedron(1, 'edge', 'icosaedro');
%   plot_polyhedron(2, 'diameter', 'dodecaedro', 0.1);
%   plot_polyhedron(1, 'edge', 'icosidodecaedro', 0.05);
%
% Autor: [Seu Nome]
% Data: [Data]

% Valores padrão para os argumentos
if nargin < 1, value = 1; end
if nargin < 2, mode = 'edge'; end
if nargin < 3, solid = 'icosidodecaedro'; end
if nargin < 4, circleRadius = []; end

%% 1. Definição comum: icosaedro padrão
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

%% 2. Seleciona qual sólido será gerado
switch lower(solid)
    
    case 'icosaedro'
        %% 2A. Icosaedro
        switch lower(mode)
            case 'edge'
                s_edge = norm(Vicos(Ficos(1,1),:) - Vicos(Ficos(1,2),:));
                escala = value / s_edge;
            case 'diameter'
                D_actual = max(pdist(Vicos));
                escala = value / D_actual;
            otherwise
                error('Modo desconhecido. Use "edge" ou "diameter".');
        end
        
        Vicos_scaled = Vicos * escala;
        
        figure; hold on;
        % Plotando cada face triangular
        for i = 1:size(Ficos,1)
            patch('Vertices', Vicos_scaled, 'Faces', Ficos(i,:), ...
                  'FaceColor', [1 0.8 0.8], 'EdgeColor', 'k', 'FaceAlpha', 0.9);
        end
        
        % Se circleRadius foi informado, plota um círculo em cada face triangular
        if ~isempty(circleRadius)
            for i = 1:size(Ficos,1)
                % Extrai os vértices da face
                tri = Vicos_scaled(Ficos(i,:), :);
                A = tri(1,:); B = tri(2,:); C = tri(3,:);
                % Calcula os comprimentos dos lados
                a = norm(B - C);
                b = norm(C - A);
                c = norm(A - B);
                % Calcula o incentro (centro do círculo inscrito)
                incenter = (a*A + b*B + c*C) / (a + b + c);
                % Calcula o vetor normal da face
                n_face = cross(B - A, C - A);
                if norm(n_face) < eps
                    n_face = [0, 0, 1];
                else
                    n_face = n_face / norm(n_face);
                end
                % Define um vetor no plano da face (por exemplo, de incenter até A)
                r_face = A - incenter;
                if norm(r_face) < eps
                    r_face = [1, 0, 0];
                else
                    r_face = r_face / norm(r_face);
                end
                % Vetor perpendicular no plano
                b_face = cross(n_face, r_face);
                if norm(b_face) < eps
                    b_face = [0, 1, 0];
                else
                    b_face = b_face / norm(b_face);
                end
                theta = linspace(0, 2*pi, 100);
                circle_pts = incenter + circleRadius * (cos(theta)' * r_face + sin(theta)' * b_face);
                plot3(circle_pts(:,1), circle_pts(:,2), circle_pts(:,3), 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
            end
        end
        
        axis equal; grid on; view(3);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        title(sprintf('Icosaedro com %s = %g', mode, value));
        hold off;
        
        edge_length = norm(Vicos_scaled(Ficos(1,1),:) - Vicos_scaled(Ficos(1,2),:));
        area_triangle = (sqrt(3)/4)*edge_length^2;
        D_approx = max(pdist(Vicos_scaled));
        
        fprintf('\nResultados para Icosaedro:\n');
        fprintf('Comprimento da aresta: %.4f\n', edge_length);
        fprintf('Área de um triângulo: %.4f\n', area_triangle);
        fprintf('Diâmetro aproximado: %.4f\n', D_approx);
        
        % Retorna os dados da malha
        V = Vicos_scaled;
        F = Ficos;
        
    case 'dodecaedro'
        %% 2B. Dodecaedro
        % Vértices do dodecaedro: centróides das 20 faces do icosaedro
        numFacesIco = size(Ficos,1);
        Vdod = zeros(numFacesIco, 3);
        for i = 1:numFacesIco
            Vdod(i,:) = mean(Vicos(Ficos(i,:),:), 1);
        end
        
        % Faces pentagonais: cada um dos 12 vértices do icosaedro gera uma face
        D_faces = cell(12,1);
        for v = 1:12
            faceIdx = find(any(Ficos == v, 2));
            if numel(faceIdx) < 3
                continue;
            end
            pts = Vdod(faceIdx, :);
            centro = mean(pts, 1);
            
            % Projeção para 2D via SVD (plano de melhor ajuste)
            [~, ~, V_svd] = svd(pts - centro, 'econ');
            proj = (pts - centro) * V_svd(:,1:2);
            
            % Obtém a ordem dos pontos usando convhull
            k = convhull(proj(:,1), proj(:,2));
            k = k(1:end-1);
            D_faces{v} = faceIdx(k);
        end
        
        % Escalonamento
        nonEmpty = find(~cellfun(@isempty, D_faces), 1);
        if isempty(nonEmpty)
            error('Nenhuma face foi encontrada para o dodecaedro.');
        end
        faceEx = D_faces{nonEmpty};
        if numel(faceEx) < 2
            error('Face para escalonamento insuficiente.');
        end
        s_edge = norm(Vdod(faceEx(1),:) - Vdod(faceEx(2),:));
        switch lower(mode)
            case 'edge'
                escala = value / s_edge;
            case 'diameter'
                D_actual = max(pdist(Vdod));
                escala = value / D_actual;
            otherwise
                error('Modo desconhecido. Use "edge" ou "diameter".');
        end
        
        Vdod = Vdod * escala;
        
        %% Plotando o dodecaedro
        figure; hold on;
        for v = 1:12
            if ~isempty(D_faces{v})
                patch('Vertices', Vdod, 'Faces', D_faces{v}', ...
                      'FaceColor', [0.8 1 0.8], 'EdgeColor', 'k', 'FaceAlpha', 0.9);
            end
        end
        axis equal; grid on; view(3);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        title(sprintf('Dodecaedro com %s = %g', mode, value));
        
        % Se circleRadius foi informado, plota um círculo em cada face pentagonal
        if ~isempty(circleRadius)
            for v = 1:12
                if ~isempty(D_faces{v}) && numel(D_faces{v}) >= 3
                    indices = D_faces{v};
                    V_face = Vdod(indices, :);
                    center_face = mean(V_face,1);
                    % Calcula o normal da face (usando os três primeiros vértices)
                    n_face = cross(V_face(2,:) - V_face(1,:), V_face(3,:) - V_face(1,:));
                    if norm(n_face) < eps
                        n_face = [0, 0, 1];
                    else
                        n_face = n_face / norm(n_face);
                    end
                    % Vetor de referência no plano da face
                    r_face = V_face(1,:) - center_face;
                    if norm(r_face) < eps
                        r_face = [1, 0, 0];
                    else
                        r_face = r_face / norm(r_face);
                    end
                    % Vetor perpendicular no plano
                    b_face = cross(n_face, r_face);
                    if norm(b_face) < eps
                        b_face = [0, 1, 0];
                    else
                        b_face = b_face / norm(b_face);
                    end
                    theta = linspace(0, 2*pi, 100);
                    circle_pts = center_face + circleRadius * (cos(theta)'*r_face + sin(theta)'*b_face);
                    plot3(circle_pts(:,1), circle_pts(:,2), circle_pts(:,3), 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
                end
            end
        end
        
        hold off;
        
        if ~isempty(faceEx) && numel(faceEx) >= 2
            edge_length = norm(Vdod(faceEx(1),:) - Vdod(faceEx(2),:));
            area_pentagon = (1/4)*sqrt(5*(5+2*sqrt(5)))*edge_length^2;
            D_approx = max(pdist(Vdod));
            r_incircle = edge_length/(2*tan(pi/5));
            area_total = 12 * area_pentagon;
            fprintf('\nResultados para Dodecaedro:\n');
            fprintf('Comprimento da aresta: %.4f\n', edge_length);
            fprintf('Área de um pentágono: %.4f\n', area_pentagon);
            fprintf('Diâmetro aproximado: %.4f\n', D_approx);
            fprintf('Raio do círculo inscrito no pentágono: %.4f\n', r_incircle);
            fprintf('Área total do dodecaedro: %.4f\n', area_total);
        else
            fprintf('\nNão foi possível calcular as propriedades do dodecaedro.\n');
        end
        
        % Retorna os dados da malha
        V = Vdod;
        F = D_faces;
        
    case 'icosidodecaedro'
        %% 2C. Icosidodecaedro
        % (a) Pontos médios das arestas do icosaedro
        arestas = [];
        for i = 1:size(Ficos,1)
            face = Ficos(i,:);
            e1 = sort(face([1 2]));
            e2 = sort(face([2 3]));
            e3 = sort(face([3 1]));
            arestas = [arestas; e1; e2; e3];  %#ok<AGROW>
        end
        arestas = unique(arestas, 'rows');
        numArestas = size(arestas,1);
        
        Vrect = zeros(numArestas, 3);
        edgeMap = containers.Map('KeyType','char','ValueType','int32');
        for i = 1:numArestas
            idx = arestas(i,:);
            mid = (Vicos(idx(1),:) + Vicos(idx(2),:))/2;
            Vrect(i,:) = mid;
            key = sprintf('%d_%d', idx(1), idx(2));
            edgeMap(key) = i;
        end
        
        % Escalonamento
        faceEx = Ficos(1,:);
        e1 = sort(faceEx([1 2]));
        e2 = sort(faceEx([2 3]));
        key1 = sprintf('%d_%d', e1(1), e1(2));
        key2 = sprintf('%d_%d', e2(1), e2(2));
        L_atual = norm(Vrect(edgeMap(key1),:) - Vrect(edgeMap(key2),:));
        switch lower(mode)
            case 'edge'
                escala = value / L_atual;
            case 'diameter'
                D_atual = max(pdist(Vrect));
                escala = value / D_atual;
            otherwise
                error('Modo desconhecido. Use "edge" ou "diameter".');
        end
        
        Vrect = Vrect * escala;
        Vicos_scaled = Vicos * escala;
        
        % (a) Faces triangulares
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
        
        % (b) Faces pentagonais: cada vértice do icosaedro gera um pentágono.
        numPent = size(Vicos,1);
        pentFaces = cell(numPent,1);
        for v = 1:numPent
            inds = [];
            for k = 1:size(arestas,1)
                if any(arestas(k,:) == v)
                    inds(end+1) = k;  %#ok<AGROW>
                end
            end
            pts = Vrect(inds,:);
            centro = mean(pts,1);
            ref = pts(1,:) - centro;
            n = Vicos_scaled(v,:) / norm(Vicos_scaled(v,:));
            b2 = cross(n, ref);
            if norm(b2) < 1e-6
                B = null(ref');
                if isempty(B)
                    b2 = [0 0 1];
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
        
        %% Plotando o icosidodecaedro
        figure; hold on;
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
        axis equal; grid on; view(3);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        title(sprintf('Icosidodecaedro com %s = %g', mode, value));
        
        % Se circleRadius foi informado, plota um círculo em cada face pentagonal
        if ~isempty(circleRadius)
            for v = 1:numPent
                if ~isempty(pentFaces{v}) && numel(pentFaces{v}) >= 3
                    indices = pentFaces{v};
                    V_face = Vrect(indices, :);
                    center_face = mean(V_face,1);
                    n_face = cross(V_face(2,:) - V_face(1,:), V_face(3,:) - V_face(1,:));
                    if norm(n_face) < eps
                        n_face = [0, 0, 1];
                    else
                        n_face = n_face / norm(n_face);
                    end
                    r_face = V_face(1,:) - center_face;
                    if norm(r_face) < eps
                        r_face = [1, 0, 0];
                    else
                        r_face = r_face / norm(r_face);
                    end
                    b_face = cross(n_face, r_face);
                    if norm(b_face) < eps
                        b_face = [0, 1, 0];
                    else
                        b_face = b_face / norm(b_face);
                    end
                    theta = linspace(0, 2*pi, 100);
                    circle_pts = center_face + circleRadius * (cos(theta)'*r_face + sin(theta)'*b_face);
                    plot3(circle_pts(:,1), circle_pts(:,2), circle_pts(:,3), 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
                end
            end
        end
        
        hold off;
        
        edge_length = norm(Vrect(triFaces(1,1),:) - Vrect(triFaces(1,2),:));
        area_pentagon = (1/4)*sqrt(5*(5+2*sqrt(5)))*edge_length^2;
        area_triangle = (sqrt(3)/4)*edge_length^2;
        D_approx = max(pdist(Vrect));
        r_incircle = edge_length/(2*tan(pi/5));
        area_total = 20*area_triangle + 12*area_pentagon;
        
        fprintf('\nResultados para Icosidodecaedro:\n');
        fprintf('Comprimento da aresta (dos pentágonos): %.4f\n', edge_length);
        fprintf('Área de um pentágono: %.4f\n', area_pentagon);
        fprintf('Área de um triângulo: %.4f\n', area_triangle);
        fprintf('Área total do icosidodecaedro: %.4f\n', area_total);
        fprintf('Diâmetro aproximado: %.4f\n', D_approx);
        fprintf('Raio do círculo inscrito no pentágono: %.4f\n', r_incircle);
        
        % Retorna os dados da malha
        V = Vrect;
        % Para o icosidodecaedro, as faces são de dois tipos:
        % - Triangulares (armazenadas em triFaces)
        % - Pentagonais (armazenadas em pentFaces)
        % Retornamos F como uma estrutura contendo ambos.
        F.tri = triFaces;
        F.pent = pentFaces;
        
    otherwise
        error('Sólido desconhecido. Use "icosaedro", "dodecaedro" ou "icosidodecaedro".');
end

end
