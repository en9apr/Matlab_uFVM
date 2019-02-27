function F = cfdConvertCelltoMatrix(f)

maxLength = 1;

for i=1:length(f)
    maxLength = max(maxLength,length(f{i}));
end

F = zeros(length(f),maxLength);

for i=1:length(f)
   if length(f{i})<maxLength
       F(i,1:length(f{i})) = f{i};
       F(i,length(f{i})+1:end) = f{i}(end);
   else
       F(i,:) = f{i};
   end       
end