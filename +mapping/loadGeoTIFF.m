% Based off of https://www.mathworks.com/help/comm/ug/visualize-coverage-maps-over-lunar-terrain-using-ray-tracing.html#mw_rtc_VisualizeCoverageMapsUsingRayTracingOverLunarTerrainExample_M_ADBCAD52
% Made-to-Work-With https://pgda.gsfc.nasa.gov/products/90, specifically with https://pgda.gsfc.nasa.gov/data/LOLA_20mpp/LDEC_60S_240MPP_ADJ.TIF
function tri = loadGeoTIFF(filename, origin)
    % Returns a Triangulation Object for use with Raytracing (and SiteViewer)
    % Read in file from `filename` root directory. isfile(filename) should evaluate to true
    % `origin` controls the origin of the map. Place the Base Station here.
    %           - of format [originlat, originlon, originht]
    [Z, R] = readgeoraster(filename,OutputType="double");   % Read in the file...
    if isprop(R, "GeographicCRS")                           % Geographic CRS stores the reference spheroid for the Moon.
        moonSpheroid = R.GeographicCRS.Spheroid;            
    elseif isprop(R, "ProjectedCRS")
        moonSpheroid = R.ProjectedCRS.GeographicCRS.Spheroid;
    else
        error("GeographicCRS.Spheroid not found in file! Please inspect the contents of `R` variable.");
    end
    % Specify the originlat, lon, height:
    originlat = origin(1); originlon = origin(2); originht = origin(3);

    % To use raytracing prop model with a triangulation object, must use a
    %   triangulation object referenced to Cartesian coords
    %%%[gridlat, gridlon] = geographicGrid(R); % Extract a grid of geographic coordiantes from the reference object for terrain data
    % Manually create the geographic grid using meshgrid:
    % Generate X, Y grid using the World Limits and Raster Size
    [xWorld, yWorld] = meshgrid(linspace(R.XWorldLimits(1), R.XWorldLimits(2), R.RasterSize(2)), ...
                                linspace(R.YWorldLimits(1), R.YWorldLimits(2), R.RasterSize(1)));
    
    % Convert the projected coordinates (X, Y) to geographic coordinates (lat, lon)
    [gridlat, gridlon] = projinv(R.ProjectedCRS, xWorld, yWorld);

    % Convert the coords to Cartesian coordinates in an east-north-up (ENU)
    %   local coordinate system. Specify origin of the ENU system using the
    %   provided origin in function args
    [x, y, z] = geodetic2enu(gridlat, gridlon, Z, originlat, originlon, originht, moonSpheroid);

    % Convert Cartesian coords to a patch with triangular faces
    [F,V] = surf2patch(x, y, z, "triangles"); % reps patch with two vars: v = vertices, f = which vertices connect to form a triangular face
    F = F(:, [1 3 2]);                        % prep to visualize object by reordering vertices
    tri = triangulation(F, V);                % then, create triangulation object from patch
end