function data = load_raw(fn)
  find = fopen(fn);
  if (find == -1)
      data = [];
      return;
  end
  data = textscan(find, '%f');
  fclose(find);
  data = data{1};
end