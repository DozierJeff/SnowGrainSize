function [effectiveRadius] = grainSizeFromAbsorption(wavelength,reflectance,cosZ)
% effectiveRadius = grainSizeFromAbsorption(wavelength,reflectance,name/value pairs)
%estimates snow grain size based on the spectral scaled absorption features
%around 1030 nm and 1260 nm along with the cosine of the illumination angle
%
%Input
% wavelength - vector, in nm
% reflectance - vector of same size, dimensionless
% cosZ - cosine illumination angle, scalar or vector (with a vector, can
%   help estimate uncertainty in the retrieval)

%Output
% effectiveRadius - matrix of estimates, a row for each absorption feature
%   (if both calculated), a column for each cosine of illumination
% fit assumes clean snow and fSCA=1
%
%based on Nolin & Dozier, RSE 2000, DOI: 10.1016/S0034-4257(00)00111-5

% thin plate spline functions loaded on first pass so stored as persistent
% variable
persistent Fun

p = inputParser;
addRequired(p,'wavelength',@(x) isvector(x) && isnumeric(x) && all(x>0))
addRequired(p,'reflectance',@(x) isvector(x) && isnumeric(x) && all(x>=0))
addRequired(p,'cosZ',@(x) isnumeric(x) && all(x>0) && all(x<=1))

parse(p,wavelength,reflectance,cosZ)

assert(length(wavelength)==length(reflectance),...
    'vectors of wavelength and reflectance must be same size')

% load thin plate interpolations, computed for radii from 30 to 2000 um
% and cosine illumination angles 0.05 to 1
if isempty(Fun)
    Fun = load('thinPlate.mat');
end
F1030r = Fun.F1030r;
F1260r = Fun.F1260r;

% scaled absorption features from the data
[scaledAbs,exitflag] = snowAbsorptionFeature(wavelength(:),reflectance(:));
if any(scaledAbs<=0) || exitflag<0
    effectiveRadius = NaN;
    return
elseif isnumeric(exitflag) && isscalar(exitflag) && exitflag==1
    featureCovered = true(1,2);
else
    featureCovered = exitflag;
end

effectiveRadius = zeros(length(scaledAbs),length(cosZ));
for k=1:length(scaledAbs)
    
    % which thin plate function to use?
    if length(scaledAbs)==2
        if k==1
            Finterp = F1030r;
        else
            Finterp = F1260r;
        end
    elseif featureCovered(1)
        Finterp = F1030r;
    else
        Finterp = F1260r;
    end
    
    % loop through cosines
    for n=1:length(cosZ)
        effectiveRadius(k,n) = Finterp(scaledAbs(k),cosZ(n));
    end
end
% interpolation returns square roots
effectiveRadius = effectiveRadius.^2;
end