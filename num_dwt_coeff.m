function n = num_dwt_coeff(N,tmpscale)

while tmpscale > 0
    n = 2*round(N/2);
    N = n/2;
    tmpscale = tmpscale-1;
end