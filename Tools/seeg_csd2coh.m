function coh = seeg_csd2coh(csd)

for r = 1:size(csd,2)
for c = 1:size(csd,3)
    Cxy = abs(csd(:,r,c)) .^ abs(2 ./ (csd(:,r,r).*csd(:,c,c)) );
    coh(r,c)    = mean(Cxy);
end
end
