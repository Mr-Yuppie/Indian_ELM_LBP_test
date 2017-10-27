%


function Psi = PCA_Train(X, PC)

[N M] = size(X);

Sigma = (X * X')/M;
  
[Psi, S, V] = svd(Sigma);% ½µÎ¬

Psi  = Psi(:, 1: PC);