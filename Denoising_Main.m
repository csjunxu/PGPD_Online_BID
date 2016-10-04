clear;
addpath('model');
GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_MeanImage\';
GT_fpath = fullfile(GT_Original_image_dir, '*.png');
TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_NoisyImage\';
TT_fpath = fullfile(TT_Original_image_dir, '*.png');
GT_im_dir  = dir(GT_fpath);
TT_im_dir  = dir(TT_fpath);
im_num = length(TT_im_dir);


par.ps = 7; % patch size
par.nlsp = 5;  % number of non-local patches
par.step = 3; % the step of two neighbor patches
par.Win = min(2*par.ps,15); % size of window around the patch
par.IteNum = 1;
delta = 0.07;
par.delta = delta;
% par.mode=2;
par.nSig=10;
par.cls_num = 32;

for eta = 1%0.9:0.05:1.1
    par.eta=eta;
%     for  lambda = 0.07:0.01:0.15
%         par.lambda = lambda;
        for c1 = 0.05:0.05:0.3
            par.c1 = c1*2*sqrt(2);
        PSNR = [];
        SSIM = [];
        CCPSNR = [];
        CCSSIM = [];
        for i = 1 : im_num
            IMin = im2double(imread(fullfile(TT_Original_image_dir,TT_im_dir(i).name) ));
            IM_GT = im2double(imread(fullfile(GT_Original_image_dir, GT_im_dir(i).name)));
            S = regexp(TT_im_dir(i).name, '\.', 'split');
            IMname = S{1};
            [h,w,ch] = size(IMin);
            fprintf('%s: \n',TT_im_dir(i).name);
            CCPSNR = [CCPSNR csnr( IMin*255,IM_GT*255, 0, 0 )];
            CCSSIM = [CCSSIM cal_ssim( IMin*255, IM_GT*255, 0, 0 )];
            fprintf('The initial PSNR = %2.4f, SSIM = %2.4f. \n', CCPSNR(end), CCSSIM(end));
            %% 3 channels
            IMout = zeros(size(IMin));
            for cc = 1:ch
                IMin_cc = IMin(:,:,cc);
                IM_GT_cc = IM_GT(:,:,cc);
                fprintf('Channel %d: The initial PSNR = %2.4f, SSIM = %2.4f. \n', cc, csnr( IMin_cc*255,IM_GT_cc*255, 0, 0 ), cal_ssim( IMin_cc*255, IM_GT_cc*255, 0, 0 ));
                % read clean image
                par.I = IM_GT_cc;
                par.nim =   IMin_cc;
                % PGPD denoising
                [IMout_cc,par]  =  Denoising(par);
                fprintf('Channel %d: The final PSNR = %2.4f, SSIM = %2.4f. \n', cc, csnr( IMout_cc*255, IM_GT_cc*255, 0, 0 ), cal_ssim( IMout_cc*255, IM_GT_cc*255, 0, 0 ));
                IMout(:,:,cc) = IMout_cc;
            end
            %% output
            PSNR = [PSNR csnr( IMout*255, IM_GT*255, 0, 0 )];
            SSIM = [SSIM cal_ssim( IMout*255, IM_GT*255, 0, 0 )];
            %% output
            %             imwrite(IMout, ['../cc_Results/Real_PGPD/OnlinePGPD_3Chs_BID_' IMname '.png']);
        end
        mPSNR = mean(PSNR);
        mSSIM = mean(SSIM);
        mCCPSNR = mean(CCPSNR);
        mCCSSIM = mean(CCSSIM);
        save(['C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\OnlinePGPD_3Chs_BID_' num2str(eta) '_' num2str(c1) '.mat'],'PSNR','mPSNR','SSIM','mSSIM','CCPSNR','mCCPSNR','CCSSIM','mCCSSIM');
    end
end