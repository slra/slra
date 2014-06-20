T = 400;
t = 1:T;
p0 = 2*cos((2*pi*t)/3) + cos((2*pi*t)/7) + cos((2*pi*t)/10);
p0pow = rms(p0);


for i=1:2
  vn = sprintf('p%d', i+5);
  sig = (0.01* p0pow) * (10^(i-1));
  pn = ((-1).^(1:T)) * sig;
  SNR = 20 * log10(p0pow / rms(pn))
  eval([vn ' =  p0 + pn; ']);
  save([vn '.txt'], vn, '-ascii')
  
%  figure;
%  eval(['plot(' vn ');']);
  
  vn = sprintf('s%d', i+5);
  fid = fopen([vn '.txt'], 'wt');
  fprintf(fid, '1 1 7 6\n 394\n 7\n');
  fclose(fid);
  
end  
  
  
  
