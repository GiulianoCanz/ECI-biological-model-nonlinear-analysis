function H = p_hill_function(var,c_c,hill_coeff,ch_fold) %pointwise hill function
Hm = 1./(1+(var./c_c).^(hill_coeff));
H = Hm + ch_fold.*(1- Hm) ;
%H = (1./((1+var./c_c).^hill_coeff)) + ch_fold.*((var./c_c).^hill_coeff)./(1+(var./c_c).^hill_coeff);
end
