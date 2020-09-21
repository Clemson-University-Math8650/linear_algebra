clear all;
close all;

n = [1000; 2000; 4000];
t_dense_mat = [6.7772; 72.5228; 833.359];
t_tridiag_mat = [0.044787; 0.180129; 0.73261];
t_symm_mat = [6.19006; 49.4631; 397.558];

fg1 = figure(1);
plot(n,t_dense_mat, n,t_tridiag_mat, n,t_symm_mat);
xlabel("n");
ylabel("t (s)");
legend("dense matrix", "tridiagonal matrix", "symmetric matrix", 'Location', 'NW');
axis tight;
saveas(fg1, "compare_mat_performance.pdf");