%***************************************************************************
%                            RECONS.M
%***************************************************************************
% PURPOSE: Reconstructs the time domain vector of the analysed signal.
%***************************************************************************
% CALL:    [vrecn]=recons(trec,freq,damp,ampl,phas)
%***************************************************************************
% INPUT:   trec   -- sampled time sequence of simulated signal for plotting
%          freq   -- estimated frequency 
%          damp   -- estimated damping factor
%          ampl   -- estimated amplitude
%          phas   -- estimated phase
% OUTPUT:  vrecn  -- reconstructed time vector
%***************************************************************************

  function [vrecn]=recoons(trec,freq,damp,ampl,phas)

%-----------reconstruction for time vector---------------
vrecn=[zeros(size(trec))];
%vrecn2=[zeros(size(trec))];
  for i=1:length(freq),
      peak=exp(-damp(i)*trec+sqrt(-1)*(2*pi*freq(i)*trec+phas(i)*pi/180));
      peak=ampl(i)*peak;
      %peak2=ampl(i).*exp(-damp(i)*trec).*cos(2*pi*freq(i)*trec+phas(i)*pi/180);
      vrecn=vrecn+peak;
      %vrecn2=vrecn2+peak2;
  end
%keyboard