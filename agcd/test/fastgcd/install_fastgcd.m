% The following two lines work only in MATLAB 2014b..
%F = webread('http://www.unilim.fr/pages_perso/paola.boito/fast_gcd_wls.m')
%T = textscan(F, '%s', 'Delimiter', {'\n'})

F = fopen('download/fast_gcd_wls.m');
if F == -1
  error('Please download the file http://www.unilim.fr/pages_perso/paola.boito/fast_gcd_wls.m to "download" subdirectory and rerun the file')
end
T = textscan(F, '%s', 'Delimiter', {'\n'})
fclose(F);

F_strs = T{1,1};
fun_heads = (~cellfun('isempty', strfind(F_strs, '----')));
fun_str_num = find(fun_heads);
fun_str_num = fun_str_num(1:2:end);

fun_names = cell(length(fun_str_num)+1,1);

for (i=1:length(fun_str_num))
  SS =  textscan( F_strs{fun_str_num(i)+1}, '%% function %s');
  SS = SS{1,1};
  fun_names{i+1} = SS{1};
end    
fun_begins = [1; fun_str_num];
fun_ends = [fun_str_num-1;  size(F_strs,1)];

fun_names{1} = 'fast_gcd_wls';

for i=1:length(fun_names)
  fid = fopen([fun_names{i} '.m'],'w');
  if (strcmp(fun_names{i}, 'c_f_newton_iter'))
    [~,ind] = max(~cellfun('isempty', strfind(F_strs(fun_begins(i):fun_ends(i)), 'function [')));
    ind = ind + fun_begins(i) - 1;
    F_strs(ind) = strrep(F_strs(ind), '[z,res]', '[z,res,q1,q2,num_iter]');
  end    
  
  fprintf(fid, '%s', strjoin(F_strs(fun_begins(i):fun_ends(i))','\n'));  
  fclose(fid);  
end
