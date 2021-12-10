%% 2D Reconstruction code for Bio-FlatScope


% Set to working directory
protoDir = '.\ExampleData\2D_recon\';

psfExt ='.mat'; % '.mat', '.tiff';
capExt ='.tiff'; % '.mat', '.tiff', '.avi'; % For avi, select only one file

H = 2048; W = 3072; %Imaging SourceMos
binFlag = 0;

%% Params
gamm_l2 = 5000; % Typical 5000-20000 knob for regularization
maxIters = 1; % Don't need to change

%% Get PSF and test image file names and directories
% %{      
switch psfExt
    case '.mat'
        [psfFileList, psfDir] = uigetfile([protoDir, '.\*.mat'], 'Select PSF','Multiselect','on');
        if ~iscell(psfFileList)
            psfFileList = {psfFileList};
        end
    case '.tiff'
        [psfFileList, psfDir] = uigetfile([protoDir, '.\*.tiff'], 'Select PSF','Multiselect','on');
        if ~iscell(psfFileList)
            psfFileList = {psfFileList};
        end
end
%}
switch capExt
    case '.mat'
        [testFileList,testCapDir] = uigetfile([protoDir, '.\*.mat'], 'Select Test','Multiselect','on'); 
        if ~iscell(testFileList)
            testFileList = {testFileList};
        end
        imageCount = length(testFileList);
    case '.tiff'
        [testFileList,testCapDir] = uigetfile([protoDir, '.\*.tiff'], 'Select Test','Multiselect','on');
        if ~iscell(testFileList)
            testFileList = {testFileList};
        end
        imageCount = length(testFileList);
    case '.avi'
        [testFileList,testCapDir] = uigetfile([protoDir, '\*.avi'], 'Select Test');
        if ~iscell(testFileList)
            testFileList = {testFileList};
        end
        vr = VideoReader(fullfile(testCapDir,testFileList{1}));
        k = 1;
        mov = [];
        while hasFrame(vr)
              mov(:,:,:,k) = im2double(readFrame(vr));
            k = k+1;
        end
        imageCount = k-1;
end
outDir = [testCapDir,'/recon/'];
mkdir(outDir)
%% Binning if selected
if binFlag
      H = H/2;
      W = W/2;    
end
%% Filter operators
Fx = @(x) fft2(fftshift(x));
FiltX = @(H,x) real(ifftshift(ifft2(H.*Fx(x))));

%% Loops for test images and PSFs
for tt = 1:imageCount % Loop for multiple files
      Xt_Stack = [];
      
      if strcmp(capExt,'.avi')
            testfile = squeeze(mov(:,:,1,tt));
            testName = sprintf('%sframe%03d',vr.Name,tt);
      else
      testFile = testFileList{tt};
      [~,testName] = fileparts(testFile);
      end

      for pp = 1:length(psfFileList) % Loop for multiple PSFs
            psfFile1 = psfFileList{1};
            psfFile = psfFileList{pp};

%% PSF loading
            switch psfExt
                case '.mat'
                    Ichs = matfile(fullfile(psfDir,psfFile));
                    Ichs = Ichs.avgCap;
                case '.tiff'
                    Ichs = imread(fullfile(psfDir,psfFile));
            end
            psf = im2double(Ichs);

            % Smooth the boundary of psf
            w_psfCut = 100; %10
            kg = fspecial('gaussian',w_psfCut*[1,1],w_psfCut/10); %2);
            crpSmth = zeros(size(psf));
            crpSmth(w_psfCut+1:end-w_psfCut,w_psfCut+1:end-w_psfCut) = 1;
            crpSmth = imfilter(crpSmth,kg,'same');
            psf = bsxfun(@times, psf, crpSmth);

            % computing padding and pad PSF
            [Ny,Nx] = size(psf(:,:,1));
            py = ceil(Ny/2); px = ceil(Nx/2);
            psf = padarray(psf,[py,px],0,'both');

            % Filter in fourier domain
            Hs = Fx(psf);   %Compute 3D spectrum
            Hs_conj = conj(Hs);
            % Forward and adjoint of filter
            Hfor = @(x) FiltX(Hs,x);
            Hadj = @(x) FiltX(Hs_conj,x);
            % H'H in fourier domain
            HtH = abs(Hs.*Hs_conj);

%% Test Image loading
            switch capExt
                  case '.mat'
                        Ichs = matfile(fullfile(testCapDir,testFile));
                        Ichs = Ichs.avgCap;
                  case '.tiff'
                        Ichs = im2double(imread(fullfile(testCapDir,testFile)));
                  case '.avi'
                        Ichs = testfile;
            end

            Ichs = im2double(Ichs); 
            b = padarray(Ichs,[py,px],0,'both');

%% Reconstruction
            fprintf('PSF%02dImage%03d ',pp,tt)
            tic;
            xFilt_mult = 1./(HtH + gamm_l2);
            numerator = Hadj(b);
            Xt_nxt = FiltX(xFilt_mult,numerator);
            toc,
                        
%% Remove bright spot and crop
            Xt_nxt(2049,3073) = 0; % Remove bright spot from deconvolution
            Xt_nxt = Xt_nxt(1449:2648,2323:3822); % Crop center region (1500x1200 px)
            Xt_Stack(:,:,pp) = single(Xt_nxt);
      end
      Xt_Stack = Xt_Stack/max(Xt_Stack(:)); %Scale through stack of PSFs
      Xt_Stacku8 = im2uint8(Xt_Stack);
%% Saves etc.
      for yy = 1:length(psfFileList)
            imwrite(Xt_Stacku8(:,:,yy),sprintf('%s/%s_%05d_%s.png',outDir,testName,gamm_l2,psfFileList{yy}));   
      end
      save(sprintf('%s/%s_%05d_%s-%s.mat',outDir,testName,gamm_l2,psfFile1,psfFile), 'Xt_Stacku8');
end
