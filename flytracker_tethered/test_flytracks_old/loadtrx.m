x = zeros(15,200);
for i = 1:200
    flstrng = sprintf('fly%03d', i+6);
    tmp = load(flstrng);
    x(:,i) = tmp.xh;
end