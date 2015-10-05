F = webread('http://www.unilim.fr/pages_perso/paola.boito/fast_gcd_wls.m')
T = textscan(F, '%s', 'Delimiter', {'\n'})

F_strs = T{1,1};
fun_heads = (~cellfun('isempty', strfind(F_strs, '----')))
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
    
    
end
