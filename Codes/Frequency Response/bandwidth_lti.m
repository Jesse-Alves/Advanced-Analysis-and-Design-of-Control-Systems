function wb=bandwidth_lti(sys)
  
  [mag, W] = sigma(sys);
  
  if max(mag)<1/sqrt(2)
    display('Resposta em frequência abaixo de -3db ')
    return   
  endif
  
  
  [aux,i]=min(abs(mag-1/sqrt(2)));
  
  w=linspace(W(i-1),W(i+1),1000);
  
  [mag, W] = sigma(sys,w);
  
  [aux,i]=min(abs(mag-1/sqrt(2)));
  
  wb=w(i);
end function
