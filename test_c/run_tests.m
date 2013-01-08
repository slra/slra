function [fields, res] = run_tests(testnos, methods, opt)
%	methods = {'slra_mex'  'slra_mex_chp'  'slra_grass'  'slra_fmincon' 'slra_reg'};
	fields = {'iter'  'fmin'  'time' 'switches' 'fit'};

	for k = 1:length(fields)
	  res{k} = zeros(length(testnos), length(methods));
	end

	for i = 1:length(testnos)
		fprintf('Test #%d: ', testnos(i));
		for j = 1:length(methods) 
		  fprintf('%s, ', methods{j});
		  warning('off', 'all')
		  eval(['info1 = run_test(@' methods{j} ', testnos(i), opt);']);
		  warning('on', 'all')
		  for k = 1:length(fields)
		    if (isfield(info1, fields{k}))
		      eval(['res{k}(i,j) = info1.' fields{k} ';']);
		    end;  
		  end;
		end
		fprintf('\n');
	end  
end




