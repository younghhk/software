function [value] = quantile_values(Q_ordered, Z_ordered)
value =  log(sum(Z_ordered.* Q_ordered .* (Z_ordered >= 0) +  Z_ordered .* (Q_ordered-1) .* (Z_ordered < 0)));
end
