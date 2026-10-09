% sigma = 0.1;
% x = linspace(-1,5,1000);
% h = mean(diff(x));
% f = sin(x) + randn(size(x))*sigma;
% plot(f)
% 
% sig = estimate_sigma(f,h);
% 
% disp(abs(sig-sigma)/sigma)

% % another way to do this:
% 
%      k=6; % choose filter that has first k moments zero
%     C = fdcoeffF(k,0,-k:k);
%     filter = C(:,end);
%     filter = filter/norm(filter,2);
%     sig = rms(conv(f,filter,'valid'));

% function sig = estimate_sigma(f)
% 
%     Ih = (f(4:end-1)-f(2:end-3))/2;
%     I2h = (f(5:end)-f(1:end-4))/4;    
%     sig = sqrt(8/5)*rms(Ih-I2h);
% 
%     if sig<0.01
%         I4h = (f(9:end)-f(1:end-8))/8;
%         Ih_1 = 4/3*(Ih - 1/4*I2h);
%         I2h_1 = 4/3*(I2h(3:end-2) - 1/4*I4h);
%         sig = sqrt(576/714)*rms(Ih_1(3:end-2)-I2h_1);
%     end 
%     
% end

function sig = estimate_sigma(f,k,dim)

    C = fdcoeffF(k,0,-k-2:k+2);
    filter = C(:,end);
    filter = filter/norm(filter,2);
    if dim>1
        sig = rms(reshape(convn(permute(f,[dim 1:dim-1 dim+1:length(size(f))]),filter(:),'valid'),[],1));
    else
        sig = rms(conv(f,filter(:),'valid'));
    end

end

function sig = estimate_sigma(f,t)

    k = 6;

    [M,d] = size(f);
    if max(abs(diff(t)-mean(diff(t))))>10^-14
        vec = zeros(M-2*k-5,d);
        for j=k+3:length(f)-k-2
            C = fdcoeffF(k,t(j),t(j-k-2:j+k+2));
            filter = C(:,end);
            filter = filter/norm(filter,2);
            vec(j,:) = f(j-k-2:j+k+2,:)'*filter;
        end
        sig = rms(vec(:));
    else
        C = fdcoeffF(k,0,-k-2:k+2);
        filter = C(:,end);
        filter = filter/norm(filter,2);
        sig = rms(reshape(conv2(filter(:),1,f,'valid'),[],1));
    end
end