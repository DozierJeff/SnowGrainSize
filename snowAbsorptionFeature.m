function [svec,varargout] = snowAbsorptionFeature(wavelength,refl,varargin)
% [svec] = snowAbsorptionFeature(wavelength,refl [,'plot',true/false])
% [svec,exitflag] = snowAbsorptionFeature(wavelength,refl,____)
% [svec,exitflag,intgR] = snowAbsorptionFeature(wavelength,refl,____)
%nowAbsorptionFeature -- calculates snow absorption features at 1030 nm and
%1260 nm from wavelength and reflectance
%
%Input
%   wavelength - vector in nm (can fix to allow other units, but nm now)
%   refl - reflectance corresponding to wavelength
%Optional input, name-value pair for logical variable
%   'plot' - plot results, default false
%
%Output (both features if the wavelengths go out to at least 1328 nm, if
%not, just the feature for the 1030 nm feature is returned
%   svec - vector of scaled integrals of the two features
%Optional output, in order
%   exitflag - 1 if both absorption bands calculated, otherwise the
%       feature-covered logical vector that indicates which feature
%       calculated
%   intgR - band integrated reflectance across each of the two features, if
%       wavelength ranges covers, otherwise just the feature calculated
%
%based on Nolin & Dozier, RSE 2000, DOI: 10.1016/S0034-4257(00)00111-5

p = inputParser;
addRequired(p,'wavelength',@(x) isvector(x) && all(x>0))
addRequired(p,'refl',@(x) isvector(x) && all(x>=0))
addParameter(p,'plot',false,@(x) isscalar(x) && (isnumeric(x) || islogical(x)))
parse(p,wavelength,refl,varargin{:});
doPlot = logical(p.Results.plot);
wavelength = wavelength(:);

%wavelength ranges that cover the two absorption bands
waveRange = [970 1090; 1158 1328];
% trim the input wavelengths to cover the wave range
allWave = cell(1,2);
featureCovered = false(1,2);
featureWave = [1030 1260];
for k=1:length(featureCovered)
    k1 = find(wavelength<=waveRange(k,1),1,'last');
    k2 = find(wavelength>=waveRange(k,2),1,'first');
    if ~isempty(k1) && ~isempty(k2)
        featureCovered(k) = true;
        allWave{k} = wavelength(k1:k2);
    else
        warning('wavelengths go from %g to %g nm, not sufficient to cover the %g feature',...
            min(wavelength),max(wavelength),featureWave(k))
    end
end
assert(any(featureCovered),'wavelengths cover neither feature')

% adjust allWave in case just one feature covered
if ~all(featureCovered)
    if featureCovered(1)
        allWave = allWave(1);
        waveRange = waveRange(1,:);
    else
        allWave = allWave(2);
        waveRange = waveRange(2,:);
    end
end

% interpolating function that covers the input reflectance
Fr = pchip(wavelength,refl(:));

% scaled differences across each band
waveRange = waveRange.';
svec = zeros(1,length(allWave));
intg = zeros(size(svec));
for k=1:length(allWave)
    % functions for the continuum just linear across each band
    x = waveRange(:,k);
    y = fnval(Fr,x);
    a = -((x(2)*y(1)-x(1)*y(2))/(x(1)-x(2)));
    b = -((y(2)-y(1))/(x(1)-x(2)));
    contValue = a+b*allWave{k};
    curveValue = fnval(Fr,allWave{k});
    sdiff = contValue-curveValue;
    % integral of the scaled differences
    Fpp = pchip(allWave{k},sdiff./contValue);
    svec(k) = diff(fnval(fnint(Fpp),[waveRange(1,k) waveRange(2,k)]))/...
        (waveRange(2,k)-waveRange(1,k));
    % integral under the reflectance itself (to get other measures)
    intg(k) = diff(fnval(fnint(Fr),[waveRange(1,k) waveRange(2,k)]))/...
        (waveRange(2,k)-waveRange(1,k));
    if doPlot
        plot(allWave{k},[curveValue contValue],'k','lineWidth',1)
        hold on;
    end
end

if nargout>1
    if all(featureCovered)
        varargout{1} = 1;
    else
        varargout{1} = featureCovered;
    end
    if nargout>2
        varargout{2} = intg;
    end
end

end