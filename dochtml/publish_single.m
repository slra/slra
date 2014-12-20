function publish_single(name)
  publish(name, struct('format', 'latex', 'outputDir', '.', 'stylesheet','mxdom2doxygen.xsl', 'evalCode', false)); exit;
end