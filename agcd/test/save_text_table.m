function save_text_table( filename, names, X )
 
  fid = fopen(filename, 'w');
  if (~isempty(names))
    names = ['', names];
    fprintf(fid, '\t%s', names{:});
    fprintf(fid, '\n');
  end
  fclose(fid);
  save(filename, 'X', '-ascii','-append');
end

