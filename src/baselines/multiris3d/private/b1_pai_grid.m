function grid=b1_pai_grid(search)
lo=search.lowerBound; hi=search.upperBound; d=search.gridStep;
[x,y,z]=ndgrid(lo(1):d:hi(1),lo(2):d:hi(2),lo(3):d:hi(3)); grid=[x(:).';y(:).';z(:).'];
end
