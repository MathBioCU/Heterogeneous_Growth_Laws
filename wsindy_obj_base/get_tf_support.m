function m = get_tf_support(phi,N,tauhat,k_x)
%     ks = ((0:N-1)-floor(N/2)).^2;
%     errs = zeros(floor(N/2)-1,1);
%     m = 1;
%     errs(1) = 1;
%     check = 0;
%     while and(check == 0,m<=(N-3)/2)
%         m = m + 1;
%         x = -1:(1/m):1;
%         phi_grid = [0 phi(x(2:end-1)) 0 zeros(1,N-2*m-1)];
%         phi_fft = abs(fft(phi_grid));
%         phi_fft = phi_fft/(sum(phi_fft)/2);
%         errs(m) = sum(phi_fft(1:k_x))-tauhat;
%         errs(m)
%         check = errs(m)>0;
%     end

    ks = ((0:N-1)-floor(N/2)).^2;
    errs = zeros(floor((N-3)/2)-1,1);
    for m=2:floor((N-3)/2)
        x = -1+(1/m):(1/m):1-(1/m);
        phi_grid = [0 phi(x) zeros(1,N-2*m)];
        phi_fft = abs(fft(phi_grid));
        phi_fft = fftshift(phi_fft);
        phi_fft = phi_fft/sum(phi_fft);
        errs(m-1) = abs((k_x/tauhat)^2-sum(phi_fft.*ks));
    end
    [~,m]=min(errs);
    m=m+1;

    if isempty(m)
        m=0;
    end

end