function cpIndex = getClosingParanthesisIndex(str, opIndex)

nCp = 0;
for iChar=opIndex+1:length(str)
    if strcmp(str(iChar), ')')
        if nCp==0
           cpIndex = iChar;
           return;
        else
            nCp = nCp - 1;
        end
    elseif strcmp(str(iChar), '(')
        nCp = nCp + 1;
    end   
end