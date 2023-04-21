function [r] = measure(p)
 v = [1:length(p)];
        n = 1;
        c = cumsum([0,p(:).']);%cumulative sum of the vector
        c = c/c(end); % make sur the cumulative is 1
        k=rand(1,n);
        [~,i] = histc(k,c);
        r = v(i);% map to v values


end

