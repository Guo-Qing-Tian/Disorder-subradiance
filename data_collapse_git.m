clear
% upload datas for \alpha=0.4 and \varphi=0.5\pi
datas=load('Gma_D05_alpha04.mat');datas=datas.rho_E_mean_std; 
indx_W=1:15;
gma=datas(indx_W,:);
xx=100:100:2000; % decay rates for different N
WW=10.^(linspace(-2,-0.5,20)).'; % decay rates for different W
gma=gma';
Gamma_sub = gma(1:1:end, :);
L_sub = xx(1:1:end);

W_cs=WW(indx_W);
Gamma_sub=Gamma_sub(:,:);

nn_indx=5:1:length(L_sub);

% calculate characteristic scale
xis=zeros(length(nn_indx),length(W_cs));
for nn=1:length(nn_indx)
    for j = 1:length(W_cs)
        G_col = Gamma_sub(1:nn_indx(nn), j);
        yy=G_col.';
        ns=L_sub(1:nn_indx(nn));
        yy=yy/sum(yy);
        j0=sum(ns.*yy);
        xis(nn,j)=sqrt(sum((ns-j0).^2.*yy));
    end
end

L   = xx(nn_indx);
xi1=xis';
Ws=repmat((W_cs),[1,length(L)]);

% optimized by Nelder-Mead ——
init = [mean(W_cs),1];
opts = optimset( ...
    'Display',    'iter', ...
    'TolX',       1e-6, ...
    'TolFun',     1e-6 ...
);

cost = @(p) compute_Pb_Wc_nu(Ws,xi1,L,p(1),p(2),1);
[p_opt, Cmin] = fminsearch(cost, init, opts);
Wc_opt = p_opt(1);
nu_opt   = p_opt(2);

fprintf('nu = %.4f, C_total = %.3e\n', ...
        nu_opt, Cmin);

X = (Ws-Wc_opt) .* (L .^ (1/nu_opt)); % \tilde{W}
XX=X(:);
Y = xi1 .* (L .^ (-1)); % \xi*N_{max}^{-1}
YY=Y(:);

% upload datas for \alpha=0.5 and \varphi=0.4\pi
datas=load('Gma_D04_alpha05.mat');datas=datas.rho_E_mean_std;
gma=datas(indx_W,:);
xx=100:100:2000;

gma=gma';
Gamma_sub = gma(1:1:end, :);
L_sub = xx(1:1:end);

Gamma_sub=Gamma_sub(:,:);
nn_indx=5:1:length(L_sub);


% Extract the characteristic scale from decay rates
xis=zeros(length(nn_indx),length(W_cs));
for nn=1:length(nn_indx)
    for j = 1:length(W_cs)
        G_col = Gamma_sub(1:nn_indx(nn), j);

        yy=G_col.';
        ns=L_sub(1:nn_indx(nn));
        yy=yy/sum(yy);
        j0=sum(ns.*yy);
        xis(nn,j)=sqrt(sum((ns-j0).^2.*yy));
    end
end
xi2=xis';
X2 = (Ws-Wc_opt) .* (L .^ (1/nu_opt));
Y2 = xi2 .* (L .^ (-1));

% upload datas for \alpha=0.6 and \varphi=0.3\pi
datas=load('Gma_D03_alpha06.mat');datas=datas.rho_E_mean_std;
gma=datas(indx_W,:);
xx=100:100:2000;

gma=gma';
Gamma_sub = gma(1:1:end, :);
L_sub = xx(1:1:end);

Gamma_sub=Gamma_sub(:,:);
nn_indx=5:1:length(L_sub);


% Extract the characteristic scale from decay rates
xis=zeros(length(nn_indx),length(W_cs));
for nn=1:length(nn_indx)
    for j = 1:length(W_cs)
        G_col = Gamma_sub(1:nn_indx(nn), j);

        yy=G_col.';
        ns=L_sub(1:nn_indx(nn));
        yy=yy/sum(yy);
        j0=sum(ns.*yy);
        xis(nn,j)=sqrt(sum((ns-j0).^2.*yy));
    end
end
xi3=xis';

W = Ws;
X3 = (Ws-Wc_opt) .* (L .^ (1/nu_opt));
Y3 = xi3 .* (L .^ (-1));
%%
% Find the relative, non-universal constant c. We have set c=1 for the datas corresponding to \alpha=0.4 and \varphi=0.5\pi.
XX1=X(:,1:2:end);YY1=Y(:,1:2:end);XX1=XX1(:);YY1=YY1(:);[XX1,indx]=sort(XX1);YY1=YY1(indx);
XX2=X2(:,1:2:end);YY2=Y2(:,1:2:end);XX2=XX2(:);YY2=YY2(:);[XX2,indx]=sort(XX2);YY2=YY2(indx);
XX3=X3(:,1:2:end);YY3=Y3(:,1:2:end);XX3=XX3(:);YY3=YY3(:);[XX3,indx]=sort(XX3);YY3=YY3(indx);
opts = struct();
opts.b_min = 0.02;
opts.b_max = 5.0;
opts.nGrid = 100;
opts.interpMethod = 'pchip';
opts.errorType = 'MSE';
opts.verbose = false;

[c1, err1, info] = fitScalingFactorBounded(XX1, YY1, XX2, YY2, opts);
[c2, err2, info] = fitScalingFactorBounded(XX1, YY1, XX3, YY3, opts);
%%
% data collapse
fontsize=35;boxwid=1.5;

color=slanCM(5,6);
ind=0;
for mm=1:3:size(X,2)
    ind=ind+1;
    scatter(X(:,mm), Y(:,mm),475,'o','filled','MarkerFaceColor',color(ind,:),'MarkerEdgeColor','none')
    hold on
end

ind=0;
for mm=1:3:size(X2,2)
    ind=ind+1;
    scatter(c1*X2(:,mm), Y2(:,mm),475,'s','filled','MarkerFaceColor',color(ind,:),'MarkerEdgeColor','none')
    hold on
end

ind=0;
for mm=1:3:size(X3,2)
    ind=ind+1;
    scatter(c2*X3(:,mm), Y3(:,mm),475,'d','filled','MarkerFaceColor',color(ind,:),'MarkerEdgeColor','none')
    hold on
end


xs=linspace(0.17,15,500);
p=3.5;
plot(xs,0.265./(1+(xs/5.16).^p).^(1.722/p))

ax=gca;
ax.XScale='log';
ax.YScale='log';

ax.YLim=[0.03 0.32];
ax.XLim=[0.56 19];
pbaspect([1.4 1 1])
ax.YAxis.Exponent = 0;
ax.FontSize=fontsize;%label字体大小
ax.LineWidth = boxwid;%边框厚度
ax.YTick=0.05:0.05:0.3;
ax.XLabel.String = '$(W - W_c)\,L^{1/\nu}$';
ax.XLabel.FontSize = fontsize;
ax.XLabel.Interpreter ='latex';
ax.YLabel.String = '$\xi_{\alpha}\,L^{-1}$';
ax.YLabel.FontSize = fontsize;
ax.YLabel.Interpreter ='latex';
ax.TickLabelInterpreter='latex';
ax.Box='on';
%%
function Pb = compute_Pb_Wc_nu(W, Mdata, L_list, Wc, nu, q)
    kappa=1;
    if nargin < 6
        q = 1;
    end

    [N, M] = size(W);
    residuals = [];

    % Measurement for the fits
    for p = 1:M
        Lp = L_list(p);
        x_p = (W(:,p) - Wc) .* (Lp).^(1/nu);
        y_p = Mdata(:,p)/Lp^(kappa);

        Ep = @(x) interp1(x_p, y_p, x, 'pchip', NaN);

        for j = 1:M
            if j == p
                continue;
            end

            Lj = L_list(j);
            x_j = (W(:,j) - Wc) .* (Lj).^(1/nu);
            y_j = Mdata(:,j)/Lj^(kappa);

            for i = 1:N
                xi = x_j(i);
                yi = y_j(i);
                if xi >= min(x_p) && xi <= max(x_p)
                    y_interp = Ep(xi);
                    if ~isnan(y_interp)
                        residuals(end+1) = abs(yi - y_interp)^q;
                    end
                end
            end
        end
    end

    Nover = length(residuals);
    if Nover == 0
        Pb = NaN;
    else
        Pb = (sum(residuals) / Nover)^(1/q);
        if Wc<0
            Pb = Pb + 10;
        end
    end
end

function [b_opt, fval, info] = fitScalingFactorBounded(x_ref, y_ref, x_target, y_target, opts)

    if nargin < 5
        opts = struct();
    end
    if ~isfield(opts, 'b_min'), opts.b_min = 0.1; end
    if ~isfield(opts, 'b_max'), opts.b_max = 10; end
    if ~isfield(opts, 'nGrid'), opts.nGrid = 50; end
    if ~isfield(opts, 'interpMethod'), opts.interpMethod = 'linear'; end
    if ~isfield(opts, 'errorType'), opts.errorType = 'MSE'; end
    if ~isfield(opts, 'weighting'), opts.weighting = 'uniform'; end
    if ~isfield(opts, 'verbose'), opts.verbose = false; end

    if ~isvector(x_ref) || ~isvector(y_ref) || ~isvector(x_target) || ~isvector(y_target)
        error('Input x_ref, y_ref, x_target, y_target must be vectors');
    end
    x_ref = x_ref(:);
    y_ref = y_ref(:);
    x_target = x_target(:);
    y_target = y_target(:);

    [x_ref, idx0] = sort(x_ref);
    y_ref = y_ref(idx0);
    mask = isfinite(x_ref) & isfinite(y_ref);
    x_ref = x_ref(mask); y_ref = y_ref(mask);

    [x_target, idx1] = sort(x_target);
    y_target = y_target(idx1);
    mask2 = isfinite(x_target) & isfinite(y_target);
    x_target = x_target(mask2); y_target = y_target(mask2);

    function err = cost_b(b)
        x_scaled = b * x_target;
        xmin = min(x_ref); xmax = max(x_ref);
        maskOvlp = (x_scaled >= xmin) & (x_scaled <= xmax);
        if sum(maskOvlp) < 3

            err = Inf;
            return;
        end
        xs = x_scaled(maskOvlp);
        ys = y_target(maskOvlp);

        y_ref_interp = interp1(x_ref, y_ref, xs, opts.interpMethod);

        nanMask = ~isfinite(y_ref_interp);
        if any(nanMask)
            xs(nanMask) = []; ys(nanMask) = []; y_ref_interp(nanMask) = [];
            if numel(xs) < 3
                err = Inf;
                return;
            end
        end

        switch opts.errorType
            case 'MSE'
                dif = ys - y_ref_interp;
                if strcmp(opts.weighting, 'uniform')
                    err = mean(dif.^2);
                elseif strcmp(opts.weighting, 'density')

                    err = mean(dif.^2);
                else

                    if isa(opts.weighting, 'function_handle')
                        w = opts.weighting(xs);
                        err = sum(w .* (dif.^2)) / sum(w);
                    else
                        err = mean(dif.^2);
                    end
                end
            case 'RMSE'
                err = sqrt(mean((ys - y_ref_interp).^2));
            case 'corr'

                C = corrcoef(ys, y_ref_interp);
                if numel(C) < 4
                    err = 1; 
                else
                    r = C(1,2);
                    err = 1 - r;
                end
            otherwise
                error('Unknow errorType: %s', opts.errorType);
        end
    end


    b_grid = linspace(opts.b_min, opts.b_max, opts.nGrid);
    cost_grid = nan(size(b_grid));
    for ii = 1:length(b_grid)
        cost_grid(ii) = cost_b(b_grid(ii));
    end

    if opts.verbose
        figure;
        plot(b_grid, cost_grid, 'o-');
        xlabel('b'); ylabel('Cost'); title('Coarse search of b');
        drawnow;
    end


    [~, idx_min] = min(cost_grid);

    b_lo = opts.b_min;
    b_hi = opts.b_max;

    if idx_min > 1 && idx_min < length(b_grid)

        b_lo = b_grid(max(idx_min-1,1));
        b_hi = b_grid(min(idx_min+1,length(b_grid)));

        delta = b_grid(2) - b_grid(1);
        b_lo = max(opts.b_min, b_lo - delta);
        b_hi = min(opts.b_max, b_hi + delta);
    end


    options_fmin = optimset('TolX',1e-6, 'Display','off');
    try
        [b_opt, fval] = fminbnd(@cost_b, b_lo, b_hi, options_fmin);
    catch

        [b_opt, fval] = fminbnd(@cost_b, opts.b_min, opts.b_max, options_fmin);
    end


    if (abs(b_opt - opts.b_min) < 1e-3*opts.b_min) || (abs(b_opt - opts.b_max) < 1e-3*opts.b_max)
        if opts.verbose
            warning('Optimal b near boundary; consider expanding search range.');
        end
    end


    info = struct();
    info.b_grid = b_grid;
    info.cost_grid = cost_grid;
    info.opts = opts;
    info.bracket = [b_lo, b_hi];
    info.cost_at_opt = fval;

    if opts.verbose
        figure;
        plot(x_ref, y_ref, 'o-', 'DisplayName', 'ref');
        hold on;
        plot(b_opt*x_target, y_target, 's-', 'DisplayName', sprintf('target scaled b=%.3f', b_opt));
        legend;
        xlabel('Scaled variable'); ylabel('Scaled observable');
        title(sprintf('Final data collapse with b=%.4f, cost=%.3e', b_opt, fval));
        drawnow;
    end
end

