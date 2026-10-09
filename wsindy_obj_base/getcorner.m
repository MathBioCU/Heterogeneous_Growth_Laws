function kx = getcorner(Ufft,xx,toggle_plot)
    if ~exist('toggle_plot','var')
        toggle_plot=0;
    end

    xx = xx(:);
    Ufft = Ufft(:);
    NN = length(Ufft);

    meth=1;

    if meth==1 %tilted graph
        [~,Umax] = max(Ufft);
        y = -flipud(xx(Umax:end));
        U = cumsum(flipud(Ufft(Umax:end)));
        U = U/max(U)*length(U);

        m = diff(U([1 end]))/diff(y([1 end]));
        Rtheta = [[1 m];[-m 1]]/sqrt(1+m^2);
        U_tilt = Rtheta(2,:)*[y(:)';U(:)'];
        [~,kx] = min(U_tilt);
        Ufft = flipud(cumsum(flipud(Ufft)));
        kx = NN - kx +1;
    elseif meth==8 %tilted graph, log
        [~,Umax] = max(Ufft);
        y = -flipud(xx(Umax:end));
        U = cumsum(flipud(log(Ufft(Umax:end))));
        U = U/max(U)*length(U);

        m = diff(U([1 end]))/diff(y([1 end]));
        Rtheta = [[1 m];[-m 1]]/sqrt(1+m^2);
        U_tilt = Rtheta(2,:)*[y(:)';U(:)'];
        [~,kx] = min(U_tilt);
        Ufft = flipud(cumsum(flipud(log(Ufft))));
        kx = NN - kx +1;
    elseif meth==2 %approximate maximum curvature
        [~,Umax] = max(Ufft);
        Umax=1;
        U = flipud(Ufft(Umax:end));
        U = U/sum(U)*length(U);
        tf = testfcn(U,'subinds',1,'meth','direct','param',floor(length(U)/20));
        kappa = abs(tf.test(U,1))./(1+U(tf.rads+1:end-tf.rads).^2).^(3/2);
        [~,kx] = max(kappa);
        % figure(100)
        %     findpeaks(kappa,'MinPeakDistance',floor(length(U)/10))
        %     hold on
        %     plot(kx,kappa(kx),'rx','markersize',15)
        %     hold off
        %     drawnow
        Ufft = flipud(cumsum(flipud(Ufft)));
        kx = NN-kx+1-tf.rads;
    elseif meth==3 %piecewise linear, pointwise relative l2-norm
        [~,Umax] = max(Ufft);
        y = -flipud(xx(Umax:end));
        U = cumsum(flipud(Ufft(Umax:end)));
        U = U/max(U)*length(U);

        errs = zeros(length(U),1);
        for k=1:length(U)
            [L1,L2,m1,m2,b1,b2,Ufft_av1,Ufft_av2]=build_lines(U,y,k);
            inds1 = find(Ufft_av1); inds2 = find(Ufft_av2);
            errs(k) = sqrt(sum(((L1(inds1)-Ufft_av1(inds1))./Ufft_av1(inds1)).^2) + sum(((L2(inds2)-Ufft_av2(inds2))./Ufft_av2(inds2)).^2)); % relative l2
%             plot([[L1;L2] [Ufft_av1;Ufft_av2]]);title(errs(k)); drawnow
        end
        [~,kx] = min(errs);
        kx = NN - kx + 1;
        Ufft = flipud(cumsum(flipud(Ufft)));
    elseif meth==4 %piecewise linear, pointwise relative l2-norm, cumsum(log)
        [~,Umax] = max(Ufft);
        y = -flipud(xx(Umax:end));
        U = cumsum(flipud(log(Ufft(Umax:end))));
        U = U/max(U)*length(U);

        errs = zeros(length(U),1);
        for k=1:length(U)
            [L1,L2,m1,m2,b1,b2,Ufft_av1,Ufft_av2]=build_lines(U,y,k);
            inds1 = find(Ufft_av1); inds2 = find(Ufft_av2);
            errs(k) = sqrt(sum(((L1(inds1)-Ufft_av1(inds1))./Ufft_av1(inds1)).^2) + sum(((L2(inds2)-Ufft_av2(inds2))./Ufft_av2(inds2)).^2)); % relative l2
%             plot([[L1;L2] [Ufft_av1;Ufft_av2]]);title(errs(k)); drawnow
        end
        [~,kx] = min(errs);
        kx = NN - kx + 1;
        Ufft = flipud(cumsum(flipud(log(Ufft))));
    elseif meth==5 %piecewise linear, total least squares
        [~,Umax] = max(Ufft);
        y = -flipud(xx(Umax:end));
        U = cumsum(flipud(Ufft(Umax:end)));
        U = U/max(U)*length(U);
        errs = zeros(length(U),1);
        for k=1:length(U)
            [L1,L2,m1,m2,b1,b2,Ufft_av1,Ufft_av2]=build_lines(U,y,k);
            Rtheta1 = [[1 m1];[-m1 1]]/sqrt(1+m1^2);
            Rtheta2 = [[1 m2];[-m2 1]]/sqrt(1+m2^2);    
            L1_tilt = Rtheta1(2,:)*[y(1:k)'+b1/m1;L1(:)'];
            L2_tilt = Rtheta2(2,:)*[y(k:length(U))'+b2/m2;L2(:)'];
            U1_tilt = Rtheta1(2,:)*[y(1:k)'+b1/m1;Ufft_av1(:)'];
            U2_tilt = Rtheta2(2,:)*[y(k:length(U))'+b2/m2;Ufft_av2(:)'];
            errs(k) = norm([U1_tilt(:);U2_tilt(:)],2);
            %%% regularize as close to kx from other method
    %         errs(k) = norm([U1_tilt(:);U2_tilt(:)],2)+(k-kx)^2/(NN-kx-1)^2;
        end
        [~,kx] = min(errs);
        kx = NN - kx + 1;
        Ufft = flipud(cumsum(flipud(Ufft)));
    elseif meth==6
        U = Ufft.^2; U = U/sum(U);
%         km = 5; N = 1000;
%         [a,b] = kmeans(arrayfun(@(x) find(Ufft>=x,1), (1:N)/N)',km);
%         kx = floor(sum(arrayfun(@(i)b(i)*length(find(a==i))/length(a),1:km)));
%         semilogy(S); drawnow
        kx = ceil(sqrt(sum((0:length(U)-1).^2.*U(:)')));
        Ufft = flipud(cumsum(flipud(Ufft)));
    elseif meth==7
        U = Ufft;
        pad = ceil(length(U)/50);
        U = [min(U)*ones(pad,1);U];
        phi = @(x) (1-x.^2).^3;
        locz = [];pz = [];
        m = length(U);
        numpeaks = 3;
        locs = ones(numpeaks,1);
        mtol = pad;
        while and(length(locs)<=numpeaks,m>pad)
            m=ceil(m/2);
            z=conv(U,phi((-m:m)/m),'same');
            [~,locs,~,p]=findpeaks(log(z));
            [~,I] = sort(p,'descend');    
            if length(I)<numpeaks
                locs = [locs;zeros(numpeaks-length(I),1)];
                p = [p;zeros(numpeaks-length(I),1)];
                locz = [locz,locs-pad+1];
                pz = [pz,p];
            else
                locz = [locz,locs(I(1:numpeaks))-pad+1];
                pz = [pz,p(I(1:numpeaks))];
            end
        end
        if m<mtol
            [~,Umax] = max(Ufft);
            y = -flipud(xx(Umax:end));
            U = cumsum(flipud(Ufft(Umax:end)));
            U = U/max(U)*length(U);
    
            m = diff(U([1 end]))/diff(y([1 end]));
            Rtheta = [[1 m];[-m 1]]/sqrt(1+m^2);
            U_tilt = Rtheta(2,:)*[y(:)';U(:)'];
            [~,kx] = min(U_tilt);
            kx = NN - kx +1;
        else
            pz = pz/sum(pz(:));
            locz = locz.*pz;
            kx = floor(sum(locz(:)));
        end
        Ufft = flipud(cumsum(flipud(Ufft)));
    end

    if toggle_plot

        [L1,L2,m1,m2,b1,b2,Ufft_av1,Ufft_av2]=build_lines(Ufft,xx,kx);
        Rtheta1 = [[1 m1];[-m1 1]]/sqrt(1+m1^2);
        Rtheta2 = [[1 m2];[-m2 1]]/sqrt(1+m2^2);
        L1_tilt = Rtheta1(2,:)*[xx(1:kx)'+b1/m1;L1(:)'];
        L2_tilt = Rtheta2(2,:)*[xx(kx:NN)'+b2/m2;L2(:)'];
        U1_tilt = Rtheta1(2,:)*[xx(1:kx)'+b1/m1;Ufft_av1(:)'];
        U2_tilt = Rtheta2(2,:)*[xx(kx:NN)'+b2/m2;Ufft_av2(:)'];    
        
        figure(toggle_plot)
        plot(1:NN,Ufft,'k',...
            1:kx,L1,'g',kx:NN,L2,'g',...
            ... % 1:kx,L1_tilt,'r',kx:NN,L2_tilt,'r',...
            ... % 1:kx,U1_tilt,'k--',kx:NN,U2_tilt,'k--',...
            kx,Ufft(kx),'d',...
                'markersize',20)
        % xlim([1 NN])
        % ylim([0 NN])
        % axis equal
        drawnow
    end

end

%%% other cost functions
%         errs(k) = sqrt(sum(((L1-Ufft_av1)./Ufft_av1).^2) + sum(((L2-Ufft_av2)./Ufft_av2).^2)); % relative l2   
%         errs(k) = m1^2/(1+m1^2)*norm(1:k)^2 + 1/(1+m1^2)*norm(Ufft_av1-b1)^2 + m2^2/(1+m2^2)*norm(k:NN)^2 + 1/(1+m2^2)*norm(Ufft_av2-b2)^2; % relative l2 
%         errs(k) = norm(L1-Ufft_av1,2)/norm(Ufft_av1,2)+norm(L2-Ufft_av2,2)/norm(Ufft_av2,2); % relative l2 
%         errs(k) = norm(L1-Ufft_av1,inf)/norm(Ufft_av1,inf)+norm(L2-Ufft_av2,inf)/norm(Ufft_av2,inf); % relative l2 

%%% get kx from a local minimum in errs
%     [~,kx,~,P] = findpeaks(-errs,1:NN,'minpeakdistance',10);
%     if isempty(kx)
%         [~,kx] = min(errs);
%     else
%         [~,mp] = max(P);
%         kx = kx(mp);
% %     kx = kx(ceil((end+1)/2));
% %     kx = kx(end);
%     end
%     subplot(2,1,2)
%     findpeaks(-errs,'minpeakdistance',15)
%     xlim([1 NN])

%%% dynamic plot
% if toggle_plot > 1
% 
%     if mod(k,toggle_plot)==0
%         subplot(2,1,1)
%         plot(1:NN,Ufft,'k',...
%             1:k,L1,'g',k:NN,L2,'g',...
%             k,Ufft(k),'d','markersize',20)
%         xlim([1 NN])
%         subplot(2,1,2)
%         plot(1:k,L1_tilt,'r',k:NN,L2_tilt,'r',...
%             1:k,U1_tilt,'k--',k:NN,U2_tilt,'k--','markersize',20)
%         xlim([1 NN])
%         drawnow
%     end
% 
% end

function [L1,L2,m1,m2,b1,b2,Ufft_av1,Ufft_av2]=build_lines(Ufft,xx,k)

   NN = length(Ufft);
   subinds1 = 1:k;
   subinds2 = k:NN;
   Ufft_av1 = Ufft(subinds1);
   Ufft_av2 = Ufft(subinds2);
   
   [m1,b1,L1]=lin_regress(Ufft_av1,xx(subinds1));
   [m2,b2,L2]=lin_regress(Ufft_av2,xx(subinds2));

end

function [m,b,L]=lin_regress(U,x)

   %%% endpoints fixed, continuous piecewise linear
   m = (U(end)-U(1))/(x(end)-x(1));
   b = U(1)-m*x(1);
   L = U(1)+m*(x-x(1));

   %%% least-squares, endpoints not fixed, discontinuous piecewise linear
%    A = [0*x(:)+1 x(:)];
%    c = A \ U(:);
%    b = c(1); m = c(2); L = A*c;

end