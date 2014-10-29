
% This is the spareness equation
sfunc = @(r,N) (1 - ( (sum(r./N)^2)/sum(r.^2/N) )) / (1-1/N);

mean = 25;
std = 50;
for ii = 1:400
    norm_nums = random('norm', mean, std, [500,1]);
    norm_nums(norm_nums < 0) = 0; %half-rectify
    if (ii==1)
        figure; 
        hist(norm_nums,50);
    end
    non_silent = find(norm_nums > 0);
    n_silence = min(length(non_silent), ii);
    norm_nums(non_silent(1:n_silence)) = 0; % shut off activity to see how it affects sparseness
    percent_silent(ii) = 100 * sum(norm_nums == 0) / 500;
    s_norm(ii) = sfunc(norm_nums, length(norm_nums));
    %disp(sprintf('Mean:%d  STD:%d    S:%f', mean, std(ii), s_norm(ii)));
end
figure;
plot(percent_silent, s_norm);
set(gca, 'TickDir', 'out', 'FontSize', 14);
xlabel('Percent NonFiring MCs', 'FontSize', 18);
ylabel('Sparseness Measure', 'FontSize', 18);


uniform_nums = rand(1000,1);
poiss_nums = random('poiss', 10, [1000,1]);
figure;
hist(poiss_nums);
all_ones = ones(1000,1);
all_z = zeros(1000,1);


s_uniform = sfunc(uniform_nums, length(uniform_nums))
s_norm = sfunc(norm_nums, length(norm_nums))
s_poiss = sfunc(poiss_nums, length(poiss_nums))
s_ones = sfunc(all_ones, length(all_ones))
s_zeros = sfunc(all_z, length(all_z))