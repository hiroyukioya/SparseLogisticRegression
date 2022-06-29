function [out]=soft_threshold(x, delta)

out=sign(x).*max(0, abs(x)-delta);