exp_steps(7:19,1) = 88

exp_steps(:,2) = exp_steps(:,1)/pixperm
exp_steps(:,3) = exp_steps(:,1)/pixperdeg
exp_steps(:,4) = [0:length(exp_steps)-1]/48

for i=2:length(exp_steps)
    dp(i)=exp_steps(i,1)-exp_steps(i-1,1)
end