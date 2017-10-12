function dy = dyn_class_1(t,x0,Inp)
Cs = Inp.Cs;Cb = Inp.Cb;W = Inp.W;s0 = Inp.s0;b0 = Inp.b0;
% if(t>Inp.tf*0.1)
%     W=zeros(size(W,1),size(W,2));
% end
XN = x0(1:numel(x0)/2);
XNd = x0(numel(x0)/2+1:end);

n = size(Cb,2); % number of nodes
N = reshape(XN,3,n);
Nd = reshape(XNd,3,n);

[N,Nd] = barlengthcorrect(N,Nd,Cb,b0);

Bd = Nd*Cb';
B = N*Cb';
S = N*Cs';
l_hat_2pow = 1.41^-2 * eye(size(B,2));
gamma_hat = gamma_hat_func(S,s0);
m = 1;
m_hat = m*eye(size(B,2));
lambda_hat = 0.5*l_hat_2pow*B'*(N*Cs'*gamma_hat*Cs-W)*Cb' - (1/12)*l_hat_2pow*m_hat*Bd'*Bd;
K = Cs'*gamma_hat*Cs - Cb'*lambda_hat*Cb;
Minv = [m_hat/12 zeros(size(m_hat,1),size(m_hat,2));zeros(size(m_hat,1),size(m_hat,2)) m_hat];
Ndd = (W-N*K)*Minv ;
dy = [reshape(Nd,numel(Nd),1);reshape(Ndd,numel(Ndd),1)];
end