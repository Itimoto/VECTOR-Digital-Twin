function [aut] = generateArbitraryAntenna(fc, varargin)
    % Generates an Antenna object with arbitrary properties    
    % inputs:   fc - Carrier Frequency for antenna
    %           
    %           (If loading from HFSS)
    %           filename - filename of HFSS antenna pattern, located in
    %               models/aut/ (e.g. 'custompattern-matlabex.csv')
    %           freqvector - frequency range of element. [lowerBound
    %               higherBound] in Hz
    %
    %           (If loading from custom equation)
    %           thetaphipattern - callback in form of @(theta, phi) that
    %               returns a complex value representing the magnitude and
    %               phase of the element (theta, phi - in degrees)

    % Parse Inputs
    p = inputParser;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    addRequired(p, 'fc', validScalarPosNum);
    % Pregenerated Pattern from HFSS
    addOptional(p, 'filename', []);                 % Filename to import HFSS patterns in from
    addOptional(p, 'freqvector', [fc*0.99 fc*1.01]); % Frequency vector describing range for element pattern
    % Pattern generated through (Theta, Phi) (function)
    addOptional(p, 'thetaphipattern', []);          % Complex Callback function through which theta & phi may be accessed
    parse(p, fc, varargin{:});

    % Generate a Custom Antenna Element for the Phased Array
    if ~isempty(p.Results.filename)
        % If we specified a preloaded antenna pattern:
        [pattern_phitheta, phi, theta] = antenna.helperPatternImport(fullfile('./models/aut', filename));
        %%% TODO ^^^ Have this work with any hfss-loaded pattern
        aut = phased.CustomAntennaElement('PatternCoordinateSystem','phi-theta',...
                                        'PhiAngles', phi, ...
                                        'ThetaAngles', theta, ...
                                        'MagnitudePattern', pattern_phitheta, ...
                                        'PhasePattern', zeros(size(pattern_phitheta))); % Assume phase 0 for now (TODO: INCORPORATE PHASE!)

        aut.SpecifyPolarizationPattern = false;     % TODO - Incorporate polarization, if time allows
        aut.FrequencyVector = p.Results.freqvector; % Specify frequency range for antenna

    elseif ~isempty(p.Results.thetaphipattern)
        % We'll generate an antenna with a dipole pattern graphing based on
        %   Theta and Phi
        theta = 0:180;
        phi = 0:360;
        [theta_mg, phi_mg] = meshgrid(theta, phi); % Set up meshgrid for pulling epattern
        aut = phased.CustomAntennaElement('ThetaAngles',theta, ...
                                              'PhiAngles', phi);
        % Pull complex epattern - doesn't account for polarization :(
        epattern = p.Results.thetaphipattern(theta_mg, phi_mg);
        magpattern = mag2db(abs(epattern))';     aut.MagnitudePattern = magpattern;
        phapattern = rad2deg(angle(epattern))';  aut.PhasePattern = phapattern;
      
        aut.SpecifyPolarizationPattern = false;     % TODO - Incorporate polarization, if time allows
        aut.FrequencyVector = p.Results.freqvector; % Specify frequency range for antenna

    else
        % We've got nothing. Notify the user, but through in a stock antenna.
        helper.dlog("WARNING - NO PATTERN FILE OR THETA-PHI PATTERN SPECIFIED");
        helper.dlog("PROCEEDING WITH STOCK ANTENNA");

        stockant = design(monopole, fc);
        % This is way overcomplicating it, but the main reason why we're
        % doing this is to continue working with the `phased.CustomAntennaElement` class
        az = -180:180;
        el = -90:90;

        if isprop(stockant, 'Length')
            D = stockant.Length;
        elseif isprop(stockant, 'Height')
            D = stockant.Height;
        else
            helper.dlog("ANTENNA APERTURE SIZE NOT FOUND. APPROXIMATING...");
            D = 200; % lmao (okay in the interest of Far-Field approx)
        end

        radius = (2 * (10*D)^2) / (physconst('LightSpeed') / fc); % Place in Far-Field, and then some (just in case)
        
        [az_mg, el_mg] = meshgrid(deg2rad(az), deg2rad(el));
        [x, y, z] = sph2cart(az_mg, el_mg, radius); 
        pts = [x(:)'; y(:)'; z(:)']; % Transpose to 3xN (expect length 3 x (361*180=65341)
        %[az_mg, el_mg, radius] = cart2sph(x, y, z); % Run it back, just in case.
        %az = rad2deg(az_mg(1, :)); el = rad2deg(el_mg(:, 1));

        % Pull the E field:
        [E, ~] = EHfields(stockant, fc, pts);
        % Extract [X, Y, Z] components of the Electric Field Vector:
        Ex = reshape(E(1, :), [length(el), length(az)]); % E_x component, size Az x El
        Ey = reshape(E(2, :), [length(el), length(az)]); % E_y component
        Ez = reshape(E(3, :), [length(el), length(az)]); % E_z component

        % Convert [X, Y, Z] components to Az (around-the-horizon), El (horizon-to-z-axis)
        %E_r = Ex .* cos(el_mg) .* cos(az_mg) + Ey .* cos(el_mg) .* sin(az_mg) + Ez .* sin(el_mg); <-- Should be 0 in an ideal medium,for a propagating wave
        %E_theta = -Ex .* sin(el_mg) .* cos(az_mg) - Ey .* sin(el_mg) .* sin(az_mg) + Ez .* cos(el_mg);
        %E_phi = -Ex .* sin(az_mg) + Ey .* cos(az_mg); % ^^ Backups, for future reference.
        E_az = -Ex .* sin(az_mg) + Ey .* cos(az_mg);
        E_el = Ex .* sin(el_mg) .* cos(az_mg) + Ey .* sin(el_mg) .* sin(az_mg) - Ez .* cos(el_mg);

        % Compute the Horizontal & Vertical Components:
        E_horiz = E_az;  % Cross Pol; along the X-Y plane, when viewed from the side
        E_vert = E_el; % Co-pol; along the Z-axis, when viewed from the side

        aut = phased.CustomAntennaElement('PatternCoordinateSystem','az-el',...
                                            'AzimuthAngles', az, ...
                                            'ElevationAngles', el, ...
                                            'SpecifyPolarizationPattern', true, ...
                                            'HorizontalMagnitudePattern', mag2db(abs(E_horiz)), ...
                                            'HorizontalPhasePattern', rad2deg(angle(E_horiz)), ...
                                            'VerticalMagnitudePattern', mag2db(abs(E_vert)), ...
                                            'VerticalPhasePattern', rad2deg(angle(E_horiz)));
        pattern(aut, fc);
    end
end