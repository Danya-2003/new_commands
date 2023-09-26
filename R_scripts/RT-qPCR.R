RQ_sample = E^(cq_cntrl - cq_sample)

RQ_samples_hk = c(RQ_sample_hk_1, RQ_sample_hk_2)

NE_sample_exp_i = RQ_sample/exp(mean(log(RQ_samples_hk)))


NE_sample_exp_1 = 28
NE_sample_exp_2 = 26

NE_sample_cntrl_1 = 16
NE_sample_cntrl_2 = 17

SD_sample_exp = 3
SD_sample_cntrl = 2


NE_sample_exp = c(NE_sample_exp_1, NE_sample_exp_2)
NE_sample_cntrl = c(NE_sample_cntrl_1, NE_sample_cntrl_2)

count_exp = length(NE_sample_exp)
count_cntrl = length(NE_sample_cntrl)

v = count_exp + count_cntrl - 2

t = abs(mean(NE_sample_exp)-mean(NE_sample_cntrl))/
  (((((count_exp-1)*SD_sample_exp^2) + ((count_cntrl-1)*SD_sample_cntrl^2)) / v )^(0.5)*
     (1/count_exp+1/count_cntrl)^(0.5))

A = gamma((v+1)/2) / ((v*pi)^(0.5) * gamma(v/2)) *
  (1+((t^2)/v))^(-(v+1)/2)
A
