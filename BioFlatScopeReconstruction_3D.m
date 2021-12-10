%% Reconstruction code
% Set to working directory
clear;
addpath('./fns/');
protoPSFDir = '.\ExampleData\3D_recon\';

%%
minDpth = 1; %2.6
maxDpth = 30; %3.2

%% Params
gamm_TV2D = 0.001;
gamm_sp = 0.005;
mu_V = 1;
mu_Z = 1;
mu_W = 1;
mu_T = 1;

strSparseSet = {'grpSparse','depthSparse','fullSparse'};
strSparse = 'fullSparse';
switch strSparse
    case strSparseSet{1}
        disp('Executing group sparsity on depth');
    case strSparseSet{2}
        disp('Executing depth sparsity');
    case strSparseSet{3}
        disp('Executing full sparsity');
end

maxIters = 20;

resid_tol = 1.5;
tau_inc = 1.5; %1.1;
tau_dec = 1.5; %1.1;

params = struct('gamm_TV2D',gamm_TV2D,'gamm_sp', gamm_sp,...
    'mu_V', mu_V, 'mu_Z', mu_Z, ...
    'mu_W', mu_W, 'mu_T', mu_T, ...
    'maxIters', maxIters, ...
    'resid_tol', resid_tol, ...
    'tau_inc', tau_inc, ...
    'tau_dec', tau_dec);

%% Get capture file

[testFileList,testCapDir] = deal('RawCap_beads_nonscatter.mat',...
    '.\ExampleData\3D_recon\');
frameImg = importdata([testCapDir, filesep, testFileList]);

%% Get PSF and test image file names and directories

pfload = matfile([protoPSFDir, filesep, 'psf2.60-psf4.00align_magfs.mat']);

drng =  pfload.drng(1,minDpth:maxDpth);
magfs = pfload.mags(minDpth:maxDpth,1);

idx = drng >= 2.6 & drng <= 3.2; 
drng = drng(idx);
psfs = pfload.psfs(:,:,minDpth:maxDpth);
magfs = magfs(idx);

[Ny,Nx,Nz] = size(psfs);

invZ2 = 1./drng.^2;
invZ2 = Nz*invZ2/sum(invZ2);
wght = zeros(1,1,Nz);
wght(:) = invZ2;
% no depth preference
% wght(:) = 1;

dpth_scl = 1./magfs;
dpth_scl = dpth_scl/max(dpth_scl);

dpth_scl_temp = zeros(1,1,Nz);
dpth_scl_temp(:) = dpth_scl(:);
dpth_scl = dpth_scl_temp;

%% Out dir
outDir = [testCapDir,'/reconADMM_magf_wt_2DTV_',strSparse,filesep];
mkdir(outDir)

%% Filter operators
fftshift2 = @(x) fftshift(fftshift(x,1),2);
ifftshift2 = @(x) ifftshift(ifftshift(x,1),2);

Fx2 = @(x) fft2(fftshift2(x));
FiltX2 = @(H,x) real(ifftshift2(ifft2(H.*Fx2(x))));

Fxn = @(x) fftn(fftshift(x));
FiltXn = @(H,x) real(ifftshift(ifftn(H.*Fxn(x))));

%% Process PSFs
py = 5; px = 5;

% padding fns
pad2d = @(x) padarray(x,[py,px],0,'both');

% Smooth the boundary of psf
w_psfCut = 100; %10
kg = fspecial('gaussian',w_psfCut*[1,1],w_psfCut/10); %2);
crpSmth = zeros(Ny,Nx);
crpSmth(w_psfCut+1:end-w_psfCut,w_psfCut+1:end-w_psfCut) = 1;
crpSmth = imfilter(crpSmth,kg,'same');
psfs = bsxfun(@times, psfs, crpSmth);

% Normalized psf with sum val - mono
sms = sum(sum(psfs,1),2);
psfs = bsxfun(@rdivide, psfs,sms);
psfs = 500*psfs;
% psfs = bsxfun(@times, psfs, sms(:,:,1));
psfs = pad2d(psfs);

% % Normalized psf with norm val and pad - mono
% nrms = sqrt(sum(sum(psfs.^2,1),2));
% psfs = bsxfun(@rdivide, psfs,nrms);
% psfs = pad2d(psfs);

% Filter in fourier domain
% Hd = Fx2(psfs);
Hd = bsxfun(@times, Fx2(psfs), dpth_scl);
Hd_conj = conj(Hd);
% Forward and adjoint of filter
Hdfor = @(x) FiltX2(Hd,x);
Hdadj = @(x) FiltX2(Hd_conj,x);
% H'H in fourier domain
HdtHd = abs(Hd.*Hd_conj);

%% Summation operator
Sfor = @(x) sum(x,3); % Summation forward
Sadj = @(x) repmat(x,[1,1,Nz]); % Summation adjoint

ss = zeros(size(psfs));
ss(1,1,:) = 1;
StS = real(fftn(ss));

%% Gradient function-- 2D gradient
% x-grad filter
kx = zeros(size(psfs));
kx(1,1,:) = 1;
kx(1,2,:) = -1;
Kx = fft2(kx);

% y-grad filter
ky = zeros(size(psfs));
ky(1,1,:) = 1;
ky(2,1,:) = -1;
Ky = fft2(ky);

% Psi - Forward gradient (finite difference)
Psi = @(x) deal(FiltX2(Kx,x), FiltX2(Ky,x));

% Adjoint of gradient
PsiT = @(P1,P2) FiltX2(conj(Kx),P1) + FiltX2(conj(Ky),P2);

% Psi'Psi in fourier domain
lapl = zeros(size(psfs));    %Compute laplacian in closed form. This is the kernal to compute Psi'Psi
lapl(1,1,:) = 4;
lapl(1,2,:) = -1;
lapl(2,1,:) = -1;
lapl(1,end,:) = -1;
lapl(end,1,:) = -1;
PsiTPsi = real(fft2(lapl));   %Compute power spectrum of laplacian

%% Loops for test images and PSFs
frame = 0;
%rect = [1401,692,850,850]; % [px,py] = [5,5]
while frame < 1 %hasFrame(vr)

    frame = frame+1;
    ff = 1; %frame;

    outDirFrame = [outDir, filesep, sprintf('frame_psfSumNorm_%03d/',ff)];
    mkdir(outDirFrame);

    %       frameImg = readFrame(vr);
    Ichs = im2double(frameImg(:,:,1));
    b = pad2d(Ichs);
    Stb = Sadj(b);

    %% Initialization
    Xt_prv = zeros(size(psfs));

    rhoV_prv = 0;
    [rhoZdx_prv,rhoZdy_prv] = deal(0);
    rhoW_prv = 0;
    rhoT_prv = 0;

    fdataFid = [];
    fregPen = [];
    fobj = [];
    fobj_alt = [];

    primal_res_V = [];
    dual_res_V = [];
    primal_res_Z = [];
    dual_res_Z = [];
    primal_res_W = [];
    dual_res_W = [];
    primal_res_T = [];
    dual_res_T = [];

    HX_prv = Hdfor(Xt_prv);
    [dx_prv, dy_prv] = Psi(Xt_prv);


    % fh1 = figure; title 'recon';
    % fh2 = figure; title 'primal & dual residues';
    fh4 = figure; title 'obj fn';
    drawnow;


    %% Iterate
    iter = 0;
    while iter < maxIters
        tic;

        iter = iter + 1;

        %% Primal update
        vFilt_mult = 1./(StS + mu_V);
        numr = Stb + mu_V*HX_prv + rhoV_prv;
        V_nxt = FiltXn(vFilt_mult,numr);

        [Zdx_nxt, Zdy_nxt] = shrinkIso2DTV(dx_prv + rhoZdx_prv/mu_Z,dy_prv + rhoZdy_prv/mu_Z, gamm_TV2D/mu_Z);

        W_nxt = max(Xt_prv + rhoW_prv/mu_W,0);

        switch strSparse
            case strSparseSet{1}
                T_nxt = proximal_21_grpSparse(Xt_prv + rhoT_prv/mu_T, gamm_sp/mu_T, wght);
            case strSparseSet{2}
                T_nxt = proximal_12(Xt_prv + rhoT_prv/mu_T, gamm_sp/mu_T, wght);
            case strSparseSet{3}
                T_nxt = proximal_L1(Xt_prv + rhoT_prv/mu_T, gamm_sp/mu_T, wght);
        end

        xFilt_mult = 1./(mu_V*HdtHd + mu_Z*PsiTPsi + mu_W + mu_T);
        numr = Hdadj(mu_V*V_nxt - rhoV_prv) ...
            + PsiT(mu_Z*Zdx_nxt - rhoZdx_prv, mu_Z*Zdy_nxt - rhoZdy_prv) ...
            + mu_W*W_nxt - rhoW_prv + mu_T*T_nxt - rhoT_prv;
        Xt_nxt = FiltX2(xFilt_mult, numr);

        %% Next derivatives
        HX_nxt = Hdfor(Xt_nxt);
        [dx_nxt, dy_nxt] = Psi(Xt_nxt);

        %% Dual update
        rpV = HX_nxt - V_nxt;
        rhoV_nxt = rhoV_prv + mu_V*rpV;

        rpZdx = dx_nxt - Zdx_nxt;
        rpZdy = dy_nxt - Zdy_nxt;
        rhoZdx_nxt = rhoZdx_prv + mu_Z*rpZdx;
        rhoZdy_nxt = rhoZdy_prv + mu_Z*rpZdy;

        rpW = Xt_nxt - W_nxt;
        rhoW_nxt = rhoW_prv + mu_W*rpW;

        rpT = Xt_nxt - T_nxt;
        rhoT_nxt = rhoT_prv + mu_T*rpT;

        %% Objective fn value
        switch strSparse
            case strSparseSet{1}
                penSparse = @(x) sum(squeeze(sqrt(sum(sum(x.^2,1),2)).*wght));
            case strSparseSet{2}
                penSparse = @(x) norm(sum(bsxfun(@times,abs(x),wght),3), 'fro')^2;
            case strSparseSet{3}
                penSparse = @(x) sum(squeeze(sum(sum(abs(x),1),2).*wght));
        end

        fdataFid(iter) = 0.5*(norm(b-Sfor(HX_nxt),'fro')^2);
        fregPen(iter) = gamm_sp * penSparse(Xt_nxt) ...
            + gamm_TV2D * ( sum(abs(dx_nxt(:))) + sum(abs(dy_nxt(:))) );
        fobj(iter) = fdataFid(iter) + fregPen(iter);
        fobj_alt(iter) = 0.5*(norm(b-Sfor(V_nxt),'fro')^2) ...
            + gamm_sp *  penSparse(T_nxt) ...
            + gamm_TV2D * ( sum(abs(Zdx_nxt(:))) + sum(abs(Zdy_nxt(:))) );

        % Residuals
        primal_res_V(iter) = norm(rpV(:));
        dual_res_V(iter) = mu_V*norm(HX_nxt(:)-HX_prv(:));

        primal_res_Z(iter) = sqrt(norm(rpZdx(:))^2 + norm(rpZdy(:))^2);
        dual_res_Z(iter) = mu_Z*sqrt( norm(dx_nxt(:)-dx_prv(:))^2 ...
            + norm(dy_nxt(:)-dy_prv(:))^2 );

        primal_res_W(iter) = norm(rpW(:));
        dual_res_W(iter) = mu_W*norm(Xt_nxt(:)-Xt_prv(:));

        primal_res_T(iter) = norm(rpT(:));
        dual_res_T(iter) = mu_T*norm(Xt_nxt(:)-Xt_prv(:));

        % Update mu (augmented penalties)
        [mu_V, muV_update] = ...
            penaltyUpdater(mu_V,primal_res_V(iter),dual_res_V(iter),resid_tol,tau_inc,tau_dec);
        [mu_Z, muZ_update] = ...
            penaltyUpdater(mu_Z,primal_res_Z(iter),dual_res_Z(iter),resid_tol,tau_inc,tau_dec);
        [mu_W, muW_update] = ...
            penaltyUpdater(mu_W,primal_res_W(iter),dual_res_W(iter),resid_tol,tau_inc,tau_dec);
        [mu_T, muT_update] = ...
            penaltyUpdater(mu_T,primal_res_T(iter),dual_res_T(iter),resid_tol,tau_inc,tau_dec);

        %     if mu_V < 1
        %         mu_V = 1;
        %     end

        %     if muV_update, disp('mu_V was updated'); end
        %     if muU_update, disp('mu_U was updated'); end
        %     if muZ_update, disp('mu_Z was updated'); end
        %     if muW_update, disp('mu_W was updated'); end
        %     if muT_update, disp('mu_T was updated'); end
        fprintf('mu_V: %0.2f, mu_Z: %0.2f, mu_W: %0.2f, mu_T: %0.2f.\n',...
            mu_V, mu_Z, mu_W, mu_T);

        %% Update previous estimate as the current
        Xt_prv = Xt_nxt;
        HX_prv = HX_nxt;
        [dx_prv, dy_prv] = deal(dx_nxt, dy_nxt);

        rhoV_prv = rhoV_nxt;
        [rhoZdx_prv,rhoZdy_prv] = deal(rhoZdx_nxt,rhoZdy_nxt);
        rhoW_prv = rhoW_nxt;
        rhoT_prv = rhoT_nxt;

        toc,
        %% Display
        if ~mod(iter,1)
            Xt_tmp = Xt_nxt;
            %         figure(fh1), imshow(max(Xt_tmp/max(Xt_tmp(:)),[],3)); title 'recon';

%             figure(fh2), title 'primal & dual residues';
%             subplot(4,1,1),semilogy(1:iter,[primal_res_V',dual_res_V'])
%             subplot(4,1,2),semilogy(1:iter,[primal_res_Z',dual_res_Z'])
%             subplot(4,1,3),semilogy(1:iter,[primal_res_W',dual_res_W'])
%             subplot(4,1,4),semilogy(1:iter,[primal_res_T',dual_res_T'])

            figure(fh4), semilogy(1:iter,[fdataFid',fregPen',fobj',fobj_alt']); title 'obj fn';
            legend('DataFid', 'Penalty', 'Orig obj', 'Alt obj');
            drawnow;
        end

        if ~mod(iter,5)
            Xt_Stack = Xt_nxt;
            maxVal = max(Xt_Stack(:));
            Xt_Stack = Xt_Stack/maxVal; %Scale through depth
            implay(Xt_Stack);
            save(sprintf('%s/%s_iters%03d_%gTV2D_%g%s_%gf-%gf_test.mat',outDirFrame,'Hydra_0520_4',iter,gamm_TV2D,gamm_sp,strSparse,drng(1),drng(end)), 'Xt_Stack','maxVal','drng','gamm_TV2D','gamm_sp','strSparse','iter','params','-v7.3');
        end
    end

    %% cropping
    Xt_Stack = single(Xt_nxt);

    maxVal = max(Xt_Stack(:));
    Xt_Stack = Xt_Stack/maxVal; %Scale through depth
    implay(Xt_Stack);

    %% Saves etc.

    save(sprintf('%s/%s_iters%03d_%gTV2D_%g%s_%gf-%gf.mat',outDirFrame,'hydra_0511_1',maxIters,gamm_TV2D,gamm_sp,strSparse,drng(1),drng(end)), 'Xt_Stack','maxVal','drng','gamm_TV2D','gamm_sp','strSparse','maxIters','params','-v7.3');
end