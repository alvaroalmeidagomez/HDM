function [KM,KPM] = orderPermutations(korder,d)

KM=nchoosek(1:d,korder);
KPM=nchoosek(1:d,korder+1);

end

