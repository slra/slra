function [gini] = gini_full_3( p1, p2, p3, d )
  bfn = [length(p1) length(p2) length(p3)] - 1;

  Sd = [multmat(p2, bfn(1)-d) (-multmat(p1, bfn(2)-d)) zeros(bfn(1)+bfn(2)+1-d, bfn(3)-d+1);...
        multmat(p3, bfn(1)-d) zeros(bfn(1)+bfn(3)+1-d, bfn(2)-d+1) (-multmat(p1, bfn(3)-d));...
        zeros(bfn(2)+bfn(3)+1-d, bfn(1)-d+1) multmat(p3, bfn(2)-d) (-multmat(p2, bfn(3)-d))];
  [U,S,V] = svd(Sd);
  gini = V(:,end);
end