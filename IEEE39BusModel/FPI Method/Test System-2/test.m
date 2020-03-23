% tic
% n = 200;
% A = 500;
% a = zeros(n);
% for i = 1:n
%     a(i) = max(abs(eig(rand(A))));
% end
% toc

tic
ticBytes(gcp);
n = 1024;
A = zeros(n);
parfor (i = 1:n)
    A(i,:) = (1:n) .* sin(i*2*pi/1024);
end
tocBytes(gcp)
toc