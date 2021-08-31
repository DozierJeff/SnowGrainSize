function [solvedRadius] = testWithData(Tbl)
%demo with data from spectra sent to Ann
%call with that table
muS = [0.67347 0.40604 0.51543 0.37078 0.77011 0.83557];
solvedRadius = table;
for k=1:6
    r = Tbl{:,k+1};
    er = grainSizeFromAbsorption(Tbl.wavelength,r,muS(k));
    vn = Tbl.Properties.VariableNames{k+1};
    spid = str2double(replace(vn,'ID',''));
    thisTbl = table(spid,er(1),er(2),'VariableNames',...
        {'supPixID','radius1030','radius1260'});
    solvedRadius = [solvedRadius;thisTbl]; %#ok<AGROW>
end
end

