function decades_equal(hAxes,xLimits,yLimits)

  if (nargin < 2) || isempty(xLimits)
    xLimits = get(hAxes,'XLim');
  end
  if (nargin < 3) || isempty(yLimits)
    yLimits = get(hAxes,'YLim');
  end

  logScale = diff(yLimits)/diff(xLimits);
  powerScale = diff(log10(yLimits))/diff(log10(xLimits));

  set(hAxes,'Xlim',xLimits,...
            'YLim',yLimits,...
            'DataAspectRatio',[1 logScale/powerScale 1]);
 
 set(hAxes, 'Position', get(gca, 'OuterPosition') - ...
   get(hAxes, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);           
    set(hAxes, 'LooseInset', [0,0,0,0]);
%  m = get(hAxes,'PlotBoxAspectRatio')
%  set(hAxes,'PlotBoxAspectRatioMode', 'manual')
%  set(hAxes,'PlotBoxAspectRatio', [1 1 1]);
            
            

end
