function gendata(ells, test_examples)
  fid = fopen('table_desc.txt', 'wt');
  fprintf(fid, 'name\tT\tq\tm\tell\n');

  for i=1:length(ells)
    % Save info to the table of the test examples
    name = test_examples{i}; 
    eval([name ';']);
    T = size(y,1);
%    T1 = T; T = round(0.7 * (T1));  y = y(1:T,:); u = u(1:T,:);    
    
    m = size(u,2);
    q = size(u,2) + size(y,2);
    ell = ells(i);
    fprintf(fid, '{$%s$}\t%d\t%d\t%d\t%d\n', name, T, q, m, ell);
  
    % Save a test example
    vn = sprintf('p%d', i);
    eval([vn ' =  [u(:) ; y(:)];']);
    save([vn '.txt'], vn, '-ascii')


    vn = sprintf('s%d', i);
    fid2 = fopen([vn '.txt'], 'wt');

    Sm = (ell + 1) * q;
    fprintf(fid2, '%d %d %d %d\n', 1, q, Sm, ell*q + m); 
    fprintf(fid2, '%d\n', T-ell); 
    fprintf(fid2, '%d ', ones(1,q) * (ell+1)); 
    fclose(fid2);
    
    phi = eye(Sm);
    perm = [];
    for j=1:(ell+1)
      perm = [perm (j:(ell+1):Sm)];
    end
    
    phi = phi(:,perm);
    vn = sprintf('phi%d', i);
    save([vn '.txt'], 'phi', '-ascii')
  end

  fclose(fid);
end  

