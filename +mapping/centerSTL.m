function SceneModel = centerSTL(filename)
    % CENTERSTL Centers an STL model at the origin.
    %   SceneModel = CENTERSTL(filename) returns a structure compatible
    %   with SceneModel for visualization in siteviewer.

    % Load the STL file using an alternative reader
    [vertices, faces, normals, name] = stlReadAlternative(filename);

    % Debug: Check the types and sizes of the outputs
    disp('Vertices size:');
    disp(size(vertices));
    disp('Vertices type:');
    disp(class(vertices)); % Check the type of vertices
    disp('Faces size:');
    disp(size(faces));
    disp('Faces content:');
    disp(faces); % Display faces content
    disp('Normals size:');
    disp(size(normals));

    % Check if vertices is numeric
    if ~isnumeric(vertices)
        error('Vertices must be a numeric array.');
    end

    % Calculate the centroid of the vertices
    centroid = mean(vertices, 1);

    % Center the vertices at the origin
    centered_vertices = vertices - centroid;

    % Prepare SceneModel for siteviewer
    SceneModel.Vertices = centered_vertices; % Centered vertices
    SceneModel.Faces = faces;                % Faces (unchanged)
    SceneModel.Normals = normals;            % Normals (unchanged)
    SceneModel.Name = name;                  % Name of the model

    % Optional: Visualize the centered model
    figure;
    patch('Faces', faces, 'Vertices', centered_vertices, 'FaceColor', 'r', 'FaceAlpha', 0.5);
    title('Centered Model');
    axis equal;
end
