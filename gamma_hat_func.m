function gamma_hat = gamma_hat_func(S,s0)
k = 1;
s = diag(diag(S'*S));
gamma_hat = (eye(size(s0,1)) - diag(diag(s0./s)));
end