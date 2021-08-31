function [Tbl] = testWithModel(Tr,nr,nc)
% [Tbl] = testWithModel(Tr,nr,nc)
%run grainSizeFromAbsorption with modeled data
%Tr is the 2nd optional output from generateLookupFunction

% unique radii cosZ wavelength
ru = unique(Tr.radius);
cu = unique(Tr.cosZ);
w = unique(Tr.wavelength);

% random selections from each
Nr = sort(randperm(length(ru)-5,nr));
firstPass = true;
Nc = randperm(length(cu),nc);
while firstPass || min(cu(Nc))<.21 || max(cu(Nc))>.8
    Nc = randperm(length(cu)-5,nc)+2;
    firstPass = false;
end
Nc = sort(Nc);

Tbl = table;
for k=Nr
    for n=Nc
        t = Tr.radius==ru(k) & Tr.cosZ==cu(n);
        testCos = cu(n)+[-.2 -.1 0 .1 .2];
        er = grainSizeFromAbsorption(w,Tr.reflectance(t),testCos);
        radiusError = -sqrt(ru(k))+sqrt(er);
        thisT1 = table(1030,ru(k),cu(n),round(er(1,:)),testCos,radiusError(1,:),...
            'VariableNames',...
            {'absorption','radius','cosZ','solvedRadii','cosines','rSqrtError'});
        thisT2 = table(1260,ru(k),cu(n),round(er(2,:)),testCos,radiusError(2,:),...
            'VariableNames',...
            {'absorption','radius','cosZ','solvedRadii','cosines','rSqrtError'});
        Tbl = [Tbl;thisT1;thisT2]; %#ok<AGROW>
    end
end