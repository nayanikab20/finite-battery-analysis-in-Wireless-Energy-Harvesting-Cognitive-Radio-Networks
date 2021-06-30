function [X] = nkg_sq2( b, m, x, sampleSize)
    
    m = m;
    w = b;

    for i = 1:length(x)
      
    %pdfx(i)= x(i)^(m-1)*exp(-x(i)/(w/m)) / ( (w/m)^m *gamma(m) ) ;
    
    %nakagami
    %pdfx(i)=((2*m^m)/(gamma(m)*w^m))*x(i)^(2*m-1)*exp(-((m/w)*x(i)^2)) ;
    
    %cdfx(i) = igamma(m, x(i)/(w/m)) / gamma(m);
    cdfx(i) = gammainc(x(i)/(w/m), m) / gamma(m);
    end
    
    cdf =cdfx; % when using cdf formula
    %cdf = cumsum(pdfx); % when using pdf formulae
    %cdf = cdf-min(cdf);
    cdf = cdf/max(cdf);
    [cdf, i] = unique(cdf);
    x = x(i);
    rnd = rand(sampleSize, 1);
    X = interp1(cdf, x, rnd, 'pchip', 0);
    
    
end
