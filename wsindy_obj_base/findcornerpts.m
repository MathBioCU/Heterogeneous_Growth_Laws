% Tf = 40001; ss = 40;
% xobs = x(1:ss:Tf,:); 
% tobs = t(1:ss:Tf);
% 
% for nn=1:2
%     [corner,~] = findcornerpts(xobs(:,nn),tobs);
%     disp(corner(2))
% end

function [corner,sig_est] = findcornerpts(xobs,t,toggle_plot)

    if ~exist('toggle_plot','var')
        toggle_plot=1;
    end

    if isempty(t)
        t = 1:length(xobs);
    end

    xobs = xobs(:);
    xobs = [xobs;flipud(xobs)];
    t = t(:);
    t = [t;t(end)+mean(diff(t))+t-t(1)];

    T = length(t);
    wn = ((0:T-1)-floor(T/2))'*(2*pi)/range(t);
    xx = wn(1:ceil(end/2));
    NN = length(xx);
    Ufft = mean(reshape(abs(fftshift(fft(xobs))),T,[]),2) /sqrt(2*NN);
    Ufft = Ufft(1:ceil(T/2),:);
    [~,Umax] = max(Ufft);
    cornerdat = cumsum(abs(Ufft(1:Umax)));
    tstarind1 = getcorner(cornerdat,xx(1:Umax),toggle_plot);
%     tstarind2 = getcorner(log(abs(Ufft(1:Umax))),xx(1:Umax));
%     tstarind = floor((tstarind1+tstarind2)/2);
    tstarind = tstarind1;
%     tstarind = tstarind2;
    tstar = -xx(tstarind)/2;
    corner = [tstar max(ceil((Umax-tstarind+1)/2),1)];
    sig_est = rms(Ufft(1:min(corner(:,2)),:));
%     plot(1:Umax,cornerdat,tstarind,cornerdat(tstarind),'d','markersize',20)
end