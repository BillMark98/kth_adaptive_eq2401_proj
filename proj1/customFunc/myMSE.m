function mse = myMSE(x,y, normFactor)

% calculate the mse between two vectors, normalised with a factor
if (length(x) ~= length(y))
    error("length of x and y unequal");
end
m = length(x);
diff_err = (x - y).^2;
mse = sum(diff_err) / m / normFactor;
