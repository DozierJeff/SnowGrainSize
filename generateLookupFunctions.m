function [S,varargout] = generateLookupFunctions(varargin)
% generate thin plate lookup functions for radius as function of absorption
% features and cosine illumination for clean deep snow
%
% call with any argument except zero or false to plot spectra
%
%Output
% S - structure with sfit functions to calculate sqrt radius from 1030 nm
%   and 1260 nm features
%Optional output, table(s)
% Tbl1 - table of radii, cosZ, and the two features
% Tbl2 - table of wavelength, cosZ, wavelength, and modeled reflectances

p = inputParser;
addOptional(p,'plot',[],@(x) isscalar(x) && (isnumeric(x) || islogical(x)));
parse(p,varargin{:})
if isempty(p.Results.plot)
    doPlot = false;
else
    doPlot = logical(p.Results.plot);
end

% ranges of variables
r = round(linspace(sqrt(30),sqrt(2000),25).^2);
cz = linspace(.05,1,15);
w = 900:5:1400;
[radius,cosZ,wavelength] = ndgrid(r,cz,w);
P = struct('radius',radius(:),'cosZ',cosZ(:),'wavelength',wavelength(:),...
    'waveUnit','nm','lookup',true);
refl = cleanSnowModel(P);
if doPlot
    figure;
    plot(w,reshape(refl,[],length(w))','color',[.75 .75 .75],'linewidth',.25)
    hold on;
    y = [min(refl(:)) max(refl(:))];
    plot([1030 1030],y,'k--','linewidth',1)
    plot([1260 1260],y,'k--','linewidth',1)
    xlabel('wavelength, nm')
    ylabel('spectral reflectance')
    axis padded
end

% absorption features
Tbl = table;
for k=1:length(r)
    for n=1:length(cz)
        t = radius==r(k) & cosZ==cz(n);
        [s,d] = snowAbsorptionFeature(w,refl(t));
        thisTbl = table(r(k),cz(n),s,d,'VariableNames',...
            {'radius','cosZ','scaledAbs','intgR'});
        Tbl = [Tbl;thisTbl]; %#ok<AGROW>
    end
end

% fit thin plated splines to values, using sqrt of radius
% because generally snow reflectance vs sqrt(radius) is closer to linear
% than snow reflectance vs radius
% The F... functions calculate radius as functions of the absorption
% features and cosine illumination
ft = 'thinplateinterp';
F1030r = fit([Tbl.scaledAbs(:,1) Tbl.cosZ],sqrt(Tbl.radius),ft,'normalize','on');
F1260r = fit([Tbl.scaledAbs(:,2) Tbl.cosZ],sqrt(Tbl.radius),ft,'normalize','on');

% return as members of a structure, then save with the '-struct' option
S.F1030r = F1030r;
S.F1260r = F1260r;

%optionally return the table of absorption metrics used to generate and the
%table of reflectances
%(can be used to test retrieval)
if nargout>1
    varargout{1} = Tbl;
    if nargout>2
        varargout{2} = table(radius(:),cosZ(:),wavelength(:),refl,...
            'VariableNames',{'radius','cosZ','wavelength','reflectance'});
    end
end
end

% function to return reflectance of clean, deep snow
% input is clean snow prescription
function [refl] = cleanSnowModel(P)

% Mie scattering properties
if P.lookup
    M = lookupMie('ice',P.radius,'um',P.wavelength,P.waveUnit);
else
    error('right now, ''lookup'' must be true, but Mie calculation routine could go here')
end

% snow reflectance
refl = twostream(P.cosZ,M.omega,M.g,'delta-Eddington');

end