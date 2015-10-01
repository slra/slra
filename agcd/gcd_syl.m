function [ph, info] = gcd_syl(p, w, d, opt)
if ~exist('opt', 'var')
  opt = struct();
end
p = reshape(p,length(p), 1);
bfp = cell2mat(p);
bfn = cellfun(@length, p) - 1;
bfell = bfn - d;
Lb = bfell(2) + 1 + sum(bfell(3:end) + 1 + bfn(1));
  Kb = sum(bfn(2:end) + (bfell(1) + 1)); 
  q1 = [zeros(Lb-1,1); -p{1}; zeros(Lb-1,1)];
  q2 = [zeros(bfell(1),1)];
  for k=2:length(p) 
    q2 = [q2; p{k}; zeros(bfell(1),1)];
  end
Lb = bfell(2) + 1 + sum(bfell(3:end) + 1 + bfn(1));
  Kb = sum(bfn(2:end) + (bfell(1) + 1)); 
  q1 = [zeros(Lb-1,1); -p{1}; zeros(Lb-1,1)];
  q2 = [zeros(bfell(1),1)];
  for k=2:length(p) 
    q2 = [q2; p{k}; zeros(bfell(1),1)];
  end
  s.m = [Lb; bfell(1)+1];
  s.n = [Kb];
  phiels = cell(length(p)-1);
  phiels{end} = eye(bfell(1)+bfell(2)+2);
  for k=3:length(p)
    phiels{length(p)+1-k} = [eye(bfell(k)+1)  zeros(bfell(k)+1, bfn(1))];
  end
  s.phi = fliplr(blkdiag(phiels{:}));
if (~exist(w, 'var') || isempty(w))
  w = cellfun(@(x){ones(size(x))}, p);
end
w1 = [Inf * ones(Kb-bfn(1)-1,1); w{1}; Inf * ones(Kb-bfn(1)-1,1)];
w2 = [Inf * ones(bfell(1),1)];
for k=2:length(bfn) 
  w2 = [w2; w{k}; Inf * ones(bfell(1),1)];
end
s.w = [w1;w2];
  opt.reggamma = 10e4;
  [ph, info] = slra([q1;q2], s, size(s.phi,1)-1, opt);
  ph = ph(s.w~=Inf);
  ph = mat2cell(ph, bfn+1, [1]); 
  ph{1} = -ph{1};
end