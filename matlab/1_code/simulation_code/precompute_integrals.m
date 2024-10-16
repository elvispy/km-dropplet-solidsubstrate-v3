function [M, angles] = precompute_integrals(angles, N)
    %% Documentation
    % precomputes the integral int_{cos(angles(ii))}^{cos(angles(ii+1))}
    % P_l(u)/u**3 du
    % for l less than N, including zero. All angles should be in the
    % interval 0, pi/2 inclusive (south pole based). 
    
    % Inputs:
    %   - angles: Vector of angles to be selected. If scalar, will assume
    %   uniformly spaced angles from 0 to pi, and discard values
    %   accordingly. If Nan, will assume N+1 uniformly spaced values
    %   - N: Number of harmonic modes in the legendre decomposition (plus
    %   zero)
    
    % Outputs:
    %   - M: JxN matrix, such that M(jj, l) is the integral from agle(jj)
    %   to angles(jj+1) of order P_(l-1)
    %   - angles: vector of size J+1, representing the angles to be
    %   integrated
    
    % Example:
    % If you input precompute_integrals([pi, 3*pi/4], 3)
    % it will output precisely what wolfram alpha says:
    % https://www.wolframalpha.com/input?i=int_%7B-1%7D%5E%7Bcos%283pi%2F4%29%7D+legendreP%282%2C+u%29%2Fu%5E3+du+%29
    
    
    if isnan(angles)
        angles = linspace(0, pi, N+1);
    elseif isscalar(angles)
        angles = linspace(0, pi, angles);
    else
        if min(angles) < pi*0.5
            %warning("The angles vector has entries outside of the interval [pi/2, pi]. Rescaling values.");
        end
    end
    
    angles = angles(angles < pi * 0.49); angles = angles(angles >= 0);
    angles = [0, (angles(1:(end-1)) + angles(2:end))/2];
    M = zeros(length(angles)-1, N+1);
    for ii = 1:(length(angles)-1)
        M(ii, 2:end) = integral(@(u) collectPl(N, u) ./(u.^3), ...
            cos(angles(ii+1)), cos(angles(ii)), 'RelTol', 1e-5, ...
            'AbsTol', 1e-5, 'ArrayValued', 1);
        M(ii, 1) = integral(@(x) 1./x.^3, cos(angles(ii+1)), cos(angles(ii)));
    end
end