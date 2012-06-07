function print_fig(file_name)
xlabel('x'), ylabel('y'), title('t'), set(gca, 'fontsize', 25)
eval(sprintf('print -depsc %s.eps', file_name));
