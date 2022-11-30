
function H = hill_function(var,c_c,hill_coeff,ch_fold)
Hm = 1/(1+(var/c_c)^(hill_coeff));
H = Hm + ch_fold*(1- Hm) ;
end

