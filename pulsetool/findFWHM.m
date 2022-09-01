function fwhm = findFWHM(offsets, Mxy)


maxval = max(Mxy(:));

greaterThanHalf = Mxy > (maxval/2);

indicesGreaterThanHalf = find(greaterThanHalf==1);
firstIndex = indicesGreaterThanHalf(1);
lastIndex = indicesGreaterThanHalf(end);


x2 = offests(firstIndex);
y2 = Mxy(firstIndex);

x1 = offests(firstIndex-1);
y1 = Mxy(firstIndex-1);

m = (y2-y1)/(x2-x1);
b = y2 - m*x2;





1;
