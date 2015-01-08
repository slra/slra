%% Parameters
opt.format = 'latex';
opt.outputDir = '.';
opt.stylesheet = 'mxdom2doxygen.xsl';
opt.evalCode = false;
opt.showCode = false;

%% List of publish commands
publish_files = {'slra' 'slra_mex_obj', 'reg_slra', 'slra_ext'};

for el = 1:length(publish_files)
  publish(['../' publish_files{el} '.m'], opt);
end  