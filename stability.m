function s = stability(eigs) 
%prende in input gli autovalori eigs e restituisce:
%stabile --> s = 1
%instabile --> s = 0
%con parte immaginaria --> s = -1
count = 0;
l = length(eigs);
for j = 1:l %itero per tutta la lunghezza dell'input
   re = real(eigs(j)); %prendo la parte reale dell'autovalore j
  % im = imag(eigs(j)); %prendo la parte immaginaria dell'autovalore j
  % if im ~= 0 %se la parte immaginaria non è nulla, allora esco dal ciclo
      % s = -1;
   %    return %esce dalla funzione
   %else      %se la parte immaginaria è nulla, allora continuo
       if re < 0
           count = count + 1;
       end %end dell'if re < 0
  % end %end dell'if im != 0
  % qui ci capita solo se la funzione non è uscita
  if count == l %se tutti gli autovalori sono stabili
      s = 1;
  else
      s = 0;
  end  
end % end del for

end

