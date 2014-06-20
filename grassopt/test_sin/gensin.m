T = 400;
t = 1:T;
p0 = 2*cos((2*pi*t)/3) + cos((2*pi*t)/7) + cos((2*pi*t)/10);

for i=1:10
  vn = sprintf('p%d', i);
  eval([vn ' =  p0 + randn(1, T) * 0.02 * i; ']);
  save([vn '.txt'], vn, '-ascii')
  
  vn = sprintf('s%d', i);
  fid = fopen([vn '.txt'], 'wt');
  fprintf(fid, '1 1 7 6\n 394\n 7\n');
  fclose(fid);
end  
  
  
  
