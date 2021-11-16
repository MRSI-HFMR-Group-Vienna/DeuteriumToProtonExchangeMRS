clear all
close all

% % % % % % % % % % % % % % % % % % % % % % 

subj='SUBJ1_S2';

n_blocks=7;

cv_thresh=12;
fwhm_thresh=0.1;
snr_thresh=5;
phase_thresh=40;

source='QualityAndOutlier_Clip';

% % % % % % % % % % % % % % % % % % % % % % 

main=sprintf('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/MRSI_temporal_stability/DEU/NewBasis_%s',subj);
stat=sprintf('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/MRSI_temporal_stability/DEU/NewBasis_%s/QC',subj);
%main=sprintf('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/MRSI_temporal_stability/DEU/ExtendedBasis_%s',subj);
%stat=sprintf('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/MRSI_temporal_stability/DEU/ExtendedBasis_%s/QualityControlMask',subj);


for counter=1:5
    
    switch counter
        case 1
            meta='Cr+PCr_amp_map';
            meta_short='Cr+PCr';
        case 2
            meta='GPC+PCh_amp_map';
            meta_short='GPC+PCh';
        case 3
            meta='NAA+NAAG_amp_map';
            meta_short='NAA+NAAG';
        case 4
            meta='PCr+Cr_amp_map';
            meta_short='Cr+PCr';
        case 5
            meta='PCh+GPC_amp_map';
            meta_short='GPC+PCh';
    end
    
    cd(main);
    
    path=sprintf('%s/BLOCK1_withDeuBasis/maps/%s/%s.mnc',main,source,meta);
    if isfile(path) 
    else
        continue
    end

    command=char(['mkdir QC']);
    system (command);
    cd(stat);
    for i=1:n_blocks
        command=char([sprintf('mkdir BL%d',i)]);
        system (command);

        system(sprintf('cp %s/BLOCK%d_withDeuBasis/maps/%s/%s.mnc %s/BL%d/%s.mnc',main,i,source,meta,stat,i,meta_short));

        cd(sprintf('BL%d',i));

        system(sprintf('mnc2nii %s.mnc %s.nii',meta_short,meta_short));  
        eval(sprintf('META%d=load_nifti(''%s.nii'');',i,meta_short));

        command='rm -rf *.nii';
        system (command); 

        cd('..');
    end
    
    CV_nii=META1;

    [n_h,n_i,n_j]=size(META1.vol);
    
    META_full=zeros(n_h,n_i,n_j,n_blocks);
    
    for i=1:n_blocks
        eval(sprintf('META_full(:,:,:,%d)=META%d.vol(:,:,:);',i,i));
        eval(sprintf('clear META%d',i));
    end

    tester=zeros(n_h,n_i,n_j);

    for h=1:n_h
    for i=1:n_i
    for j=1:n_j

        for a=1:n_blocks

            if META_full(h,i,j,a) > 0
                tester(h,i,j)=tester(h,i,j)+1;
            end

        end

        if tester(h,i,j) == n_blocks
            tester(h,i,j)=1;
        else
            tester(h,i,j)=0;
        end

    end
    end
    end 

    CV=zeros(n_h,n_i,n_j);  

    for h=1:n_h
        for i=1:n_i
            for j=1:n_j

                if tester(h,i,j) == 1

                    std_bl = std(META_full(h,i,j,:));
                    mean_bl = mean(META_full(h,i,j,:));
                    cv_bl = (std_bl/mean_bl)*100;

                    CV(h,i,j)=cv_bl;

                    clear std_bl mean_bl cv_bl

                end

            end
        end
    end    

    CV_nii.vol=CV;
    eval(sprintf('save_nifti(CV_nii,''CV_%s.nii'');',meta_short));
    eval(sprintf('command=char([''nii2mnc CV_%s.nii CV_%s.mnc'']);',meta_short,meta_short));
    system (command);
    
    clear META_full tester CV_nii CV meta meta_short command

end

% % % % % % % % % % % % % % % % % % % % % % 

CV_Cr=load_nifti('CV_Cr+PCr.nii');
CV_NAA=load_nifti('CV_NAA+NAAG.nii');
CV_Ch=load_nifti('CV_GPC+PCh.nii');

mask=CV_Cr;
mask_vol=zeros(n_h,n_i,n_j);

for h=1:n_h
for i=1:n_i
for j=1:n_j

    if CV_Cr.vol(h,i,j) <= cv_thresh && CV_Cr.vol(h,i,j) > 0
    if CV_NAA.vol(h,i,j) <= cv_thresh && CV_NAA.vol(h,i,j) > 0
    if CV_Ch.vol(h,i,j) <= cv_thresh && CV_Ch.vol(h,i,j) > 0
        
        mask_vol(h,i,j)=1;
        
    end
    end
    end

end
end
end    

mask.vol=mask_vol;
name=sprintf('mask_CV%d.nii',cv_thresh);
save_nifti(mask,name);
command=char([sprintf('nii2mnc mask_CV%d.nii mask_CV%d.mnc',cv_thresh,cv_thresh)]);
system (command);

command='rm -rf *.nii';
system (command); 

clear mask mask_vol name counter command ans a CV_Ch CV_Cr CV_NAA path

% % % % % % % % % % % % % % % % % % % % 

cd(main);

for i=1:n_blocks
    
    command=char([sprintf('mkdir QualityControlMask/BL%d',i)]);
    system (command);
      
    system(sprintf('cp %s/BLOCK%d_withDeuBasis/maps/Extra/SNR_map.mnc %s/BL%d/SNR_map.mnc',main,i,stat,i));
    system(sprintf('cp %s/BLOCK%d_withDeuBasis/maps/Extra/0_pha_map.mnc %s/BL%d/0Phase_map.mnc',main,i,stat,i));  
    system(sprintf('cp %s/BLOCK%d_withDeuBasis/maps/Extra/FWHM_map.mnc %s/BL%d/FWHM_map.mnc',main,i,stat,i));  
    
    command=char([sprintf('mnc2nii %s/BL%d/SNR_map.mnc %s/BL%d/SNR_map.nii',stat,i,stat,i)]);
    system (command);
    command=char([sprintf('mnc2nii %s/BL%d/0Phase_map.mnc %s/BL%d/0Phase_map.nii',stat,i,stat,i)]);
    system (command);
    command=char([sprintf('mnc2nii %s/BL%d/FWHM_map.mnc %s/BL%d/FWHM_map.nii',stat,i,stat,i)]);
    system (command);
    
end

fwhm_all=zeros(n_h,n_i,n_j,n_blocks);
snr_all=zeros(n_h,n_i,n_j,n_blocks);
phase_all=zeros(n_h,n_i,n_j,n_blocks);

for a=1:n_blocks
    
    cd(sprintf('%s/BL%d/',stat,a));
    
    fwhm=load_nifti('FWHM_map.nii');    
    snr=load_nifti('SNR_map.nii');
    phase=load_nifti('0Phase_map.nii');

    fwhm_all(:,:,:,a)=fwhm.vol(:,:,:);
    snr_all(:,:,:,a)=snr.vol(:,:,:);
    phase_all(:,:,:,a)=phase.vol(:,:,:);
           
    mask_fwhm=fwhm;
    mask_fwhm=rmfield(mask_fwhm,'vol');
    mask_snr=snr;
    mask_snr=rmfield(mask_snr,'vol');
    mask_phase=phase;
    mask_phase=rmfield(mask_phase,'vol');
    
    vol_fwhm=zeros(n_h,n_i,n_j);
    vol_snr=zeros(n_h,n_i,n_j);
    vol_phase=zeros(n_h,n_i,n_j);
    
    fmask=fwhm;
    smask=snr;
    pmask=phase;
    
    for h=1:n_h
    for i=1:n_i
    for j=1:n_j

        if fwhm.vol(h,i,j) <= fwhm_thresh && fwhm.vol(h,i,j) > 0               
            vol_fwhm(h,i,j)=1;                    
        end
        
        if snr.vol(h,i,j) > snr_thresh               
            vol_snr(h,i,j)=1;                    
        end
        
        if phase.vol(h,i,j) <= phase_thresh && phase.vol(h,i,j) > 0               
            vol_phase(h,i,j)=1;                    
        end

    end
    end
    end    
    
    mask_fwhm.vol=vol_fwhm;
    mask_snr.vol=vol_snr;
    mask_phase.vol=vol_phase;          
       
    save_nifti(mask_fwhm,'mask_fwhm.nii');
    system ('nii2mnc mask_fwhm.nii mask_fwhm.mnc');
    save_nifti(mask_snr,'mask_snr.nii');
    system ('nii2mnc mask_snr.nii mask_snr.mnc');
    save_nifti(mask_phase,'mask_phase.nii');
    system ('nii2mnc mask_phase.nii mask_phase.mnc');
    
    clear fwhm snr phase mask_fwhm mask_snr mask_phase vol_fwhm vol_snr vol_phase
    
    command=char(['rm -rf *.nii']);
    system (command);

end

cd(stat);

combmask=fmask;

mask_fwhm=zeros(n_h,n_i,n_j);
mask_snr=zeros(n_h,n_i,n_j);
mask_phase=zeros(n_h,n_i,n_j);
mask_comb=zeros(n_h,n_i,n_j);
    
for h=1:n_h
for i=1:n_i
for j=1:n_j

    for a=1:n_blocks

        if fwhm_all(h,i,j,a) <= fwhm_thresh && fwhm_all(h,i,j,a) > 0
            mask_fwhm(h,i,j)=mask_fwhm(h,i,j)+1;
        end

        if snr_all(h,i,j,a) > snr_thresh  
            mask_snr(h,i,j)=mask_snr(h,i,j)+1;
        end

        if phase_all(h,i,j,a) <= phase_thresh && phase_all(h,i,j,a) > 0
            mask_phase(h,i,j)=mask_phase(h,i,j)+1;
        end

    end

    if mask_fwhm(h,i,j) == n_blocks
        mask_fwhm(h,i,j)=1;
    else
        mask_fwhm(h,i,j)=0;
    end
    
    if mask_snr(h,i,j) == n_blocks
        mask_snr(h,i,j)=1;
    else
        mask_snr(h,i,j)=0;
    end
    
    if mask_phase(h,i,j) == n_blocks
        mask_phase(h,i,j)=1;
    else
        mask_phase(h,i,j)=0;
    end

end
end
end

for h=1:n_h
for i=1:n_i
for j=1:n_j

    if mask_fwhm(h,i,j) == 1 && mask_snr(h,i,j) == 1 && mask_phase(h,i,j) == 1
        mask_comb(n_h,n_i,n_j)=1;
    end

end
end
end 


fmask=rmfield(fmask,'vol');
smask=rmfield(smask,'vol');
pmask=rmfield(pmask,'vol');
combmask=rmfield(combmask,'vol');

fmask.vol=mask_fwhm;
smask.vol=mask_snr;
pmask.vol=mask_phase;
combmask.vol=mask_comb;

name=sprintf('mask_FWHM%.2f.nii',fwhm_thresh);
save_nifti(fmask,name);
command=char([sprintf('nii2mnc mask_FWHM%.2f.nii mask_FWHM%.2f.mnc',fwhm_thresh,fwhm_thresh)]);
system (command);

name=sprintf('mask_SNR%d.nii',snr_thresh);
save_nifti(smask,name);
command=char([sprintf('nii2mnc mask_SNR%d.nii mask_SNR%d.mnc',snr_thresh,snr_thresh)]);
system (command);

name=sprintf('mask_PHASE%d.nii',phase_thresh);
save_nifti(pmask,name);
command=char([sprintf('nii2mnc mask_PHASE%d.nii mask_PHASE%d.mnc',phase_thresh,phase_thresh)]);
system (command);

name=sprintf('mask_FWHM%.2f_SNR%d_PHASE%d.nii',fwhm_thresh,snr_thresh,phase_thresh);
save_nifti(pmask,name);
command=char([sprintf('nii2mnc mask_FWHM%.2f_SNR%d_PHASE%d.nii mask_FWHM%.2f_SNR%d_PHASE%d.mnc',fwhm_thresh,snr_thresh,phase_thresh,fwhm_thresh,snr_thresh,phase_thresh)]);
system (command);

command=char(['rm -rf *.nii']);                                          
system (command);

% % % % % % % % % % % % % % % % % % % % 

command=char([sprintf('mincmath -mult -nocheck_dimensions mask_FWHM%.2f_SNR%d_PHASE%d.mnc mask_CV%d.mnc mask_CV%d_FWHM%.2f_SNR%d_PHASE%d.mnc',fwhm_thresh,snr_thresh,phase_thresh,cv_thresh,cv_thresh,fwhm_thresh,snr_thresh,phase_thresh)]);
system (command);

% % % % % % % % % % % % % % % % % % % % 








