clear all
close all
clear clc

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

subj='SUBJ1_S2';

nb_blocks=7;

source='QualityAndOutlier_Clip';

mask_quali_name='mask_CV12_FWHM0.10_SNR5_PHASE40';

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

main=sprintf('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/MRSI_temporal_stability/DEU/NewBasis_%s',subj);
stat=sprintf('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/MRSI_temporal_stability/DEU/NewBasis_%s/statistics/Masks_FIN',subj);
stat2=sprintf('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/MRSI_temporal_stability/DEU/NewBasis_%s/statistics/ExtractNumb_FIN',subj);            
quali=sprintf('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/MRSI_temporal_stability/DEU/NewBasis_%s/QC',subj);
fin_fsmasks=sprintf('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/MRSI_temporal_stability/DEU/%s/NEW_MASKS/MASKS/TEST',subj);

cd(main);
command=char(['mkdir statistics']);
system (command);
command=char(['mkdir statistics/Masks_FIN']);
system (command);
command=char(['mkdir statistics/ExtractNumb_FIN']);
system (command);

for counter=1:2
    
    switch counter
        case 1
            mask='whole_GM_2csi_swapped_BIN_FIN';
        case 2
            mask='whole_WM_2csi_swapped_BIN';
        case 3
            mask='whole_cingulate_CTX';
    end
    
    cd(stat);
    path=sprintf('%s/%s.mnc',stat,mask);
    if isfile(path) 
    else
        system(sprintf('cp %s/%s.mnc %s/%s.mnc',fin_fsmasks,mask,stat,mask));
        system(sprintf('cp %s/%s.nii %s/%s.nii',fin_fsmasks,mask,stat,mask));
    end    
    system(sprintf('mnc2nii %s.mnc %s.nii',mask,mask));
    mask_vol=load_nifti(sprintf('%s.nii',mask));
    
    cd(quali);
    system(sprintf('mnc2nii %s.mnc mask.nii',mask_quali_name));
    mask_quali=load_nifti('mask.nii');
    command='rm -rf *.nii';
    system (command); 
    
    big_met_count=0;
    big_holder=zeros(nb_blocks,12*3);

for met_counter=1:14
    
    switch met_counter
        case 1
            metabol='Asp';
        case 2
            metabol='GABA';
        case 3
            metabol='Gln';
        case 4
            metabol='Glu_GD';
        case 5
            metabol='Glu_GN';
        case 6
            metabol='Ins'; 
        case 7
            metabol='Ins+Gly';
        case 8
            metabol='MM_mea';
        case 9
            metabol='NAA+NAAG';
        case 10
            metabol='PCh+GPC';
        case 11
            metabol='GPC+PCh';
        case 12
            metabol='Cr+PCr';
        case 13
            metabol='PCr+Cr';    
        case 14
            metabol='Tau';
    end
    
    path=sprintf('%s/BLOCK1_withDeuBasis/maps/%s/%s_amp_map.mnc',main,source,metabol);
    if isfile(path)
        big_met_count=big_met_count+1;
    else
        continue
    end
    
    for i=1:nb_blocks
        
        cd(sprintf('%s/BLOCK%d_withDeuBasis/maps/%s/',main,i,source))
        
        system(sprintf('mnc2nii %s_amp_map.mnc %s_amp_map.nii',metabol,metabol));
        eval(sprintf(' BL%d=load_nifti(''%s_amp_map.nii''); ',i,metabol)); 
        eval(sprintf('BL(:,:,:,i)=BL%d.vol(:,:,:);',i));

        path=sprintf('%s/Cr+PCr_amp_map.mnc',stat);
        if isfile(path) 
            system(sprintf('mnc2nii Cr+PCr_amp_map.mnc Cr+PCr_amp_map.nii'));        
            eval(sprintf(' BL%d_ratio=load_nifti(''Cr+PCr_amp_map.nii''); ',i));
            eval(sprintf('BL_rat(:,:,:,i)=BL%d_ratio.vol(:,:,:);',i));
        else
            system(sprintf('mnc2nii PCr+Cr_amp_map.mnc Cr+PCr_amp_map.nii'));        
            eval(sprintf(' BL%d_ratio=load_nifti(''Cr+PCr_amp_map.nii''); ',i));
            eval(sprintf('BL_rat(:,:,:,i)=BL%d_ratio.vol(:,:,:);',i));
        end  
        
        command='rm -rf *.nii';
        system (command); 
        
    end

    [n_h, n_i, n_j]=size(BL1.vol);

    sum_BL=zeros(1,nb_blocks);
    count=zeros(1,nb_blocks);
    avg=zeros(1,nb_blocks);
    stderror=zeros(1,nb_blocks);
    stddev=zeros(1,nb_blocks);
    variance=zeros(1,nb_blocks);

    if strcmp(mask,'whole_GM_2csi_swapped_BIN_FIN') == 1   
        surv_gm=zeros(1,nb_blocks);
    elseif strcmp(mask,'whole_WM_2csi_swapped_BIN') == 1
        surv_wm=zeros(1,nb_blocks);
    end 
    
    for a=1:nb_blocks
        
    holder(1)=0;

        for h=1:n_h
            for i=1:n_i
                for j=1:n_j
                                
                    if mask_vol.vol(h,i,j) > 0  && mask_quali.vol(h,i,j) > 0 
                        
                        if strcmp(mask,'whole_GM_2csi_swapped_BIN_FIN') == 1   
                            surv_gm(a)=surv_gm(a)+1;
                        elseif strcmp(mask,'whole_WM_2csi_swapped_BIN') == 1
                            surv_wm(a)=surv_wm(a)+1;
                        end
                        
                        if isnan(BL(h,i,j,a))==0 
                            if BL(h,i,j,a) > 0 

                                if BL(h,i,j,a) > 100
                                    BL(h,i,j,a)=100;
                                end

                                count(a)=count(a)+1;
                                sum_BL(a)=sum_BL(a)+BL(h,i,j,a);
                                holder(count(a))=BL(h,i,j,a);

                            end
                        end
                    end

                end
            end
        end 
      
        avg(a)=sum_BL(a)/count(a);       
        stddev(a)=std(holder);      
        variance(a)=var(holder);
        stderror(a)=stddev(a)/sqrt(count(a));
        
        clear holder

    end   

    cd(stat2);

    name=sprintf('Variance_AddCrit_Average_%s_%s_%s.txt',metabol,mask,source);
    fileID = fopen(name,'w');
    fprintf(fileID,'mean stddev stderror %s %s\n',mask,source);
    for i=1:nb_blocks
        fprintf(fileID,'%d %d %d\n',avg(i),stddev(i),stderror(i));
    end
    fclose(fileID);
        
    big_holder(:,(big_met_count-1)*3+1)=avg(:);
    big_holder(:,(big_met_count-1)*3+2)=stddev(:);
    big_holder(:,(big_met_count-1)*3+3)=stderror(:);
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    
    sum_BL=zeros(1,nb_blocks);
    count=zeros(1,nb_blocks);
    avg=zeros(1,nb_blocks);
    stderror=zeros(1,nb_blocks);
    stddev=zeros(1,nb_blocks);
    variance=zeros(1,nb_blocks);

    for a=1:nb_blocks
        
        holder(1)=0;

        for h=1:n_h
            for i=1:n_i
                for j=1:n_j

                    if mask_vol.vol(h,i,j) > 0  && mask_quali.vol(h,i,j) > 0 
                        if isnan(BL(h,i,j,a))==0 && isnan(BL_rat(h,i,j,a))==0
                            if BL(h,i,j,a) > 0 && BL_rat(h,i,j,a) > 0

                                if BL_rat(h,i,j,a) > 100
                                    BL_rat(h,i,j,a)=100;
                                end
                                
                                if BL(h,i,j,a) > 100
                                    BL(h,i,j,a)=100;
                                end
                                
                                count(a)=count(a)+1;
                                sum_BL(a)=sum_BL(a)+BL(h,i,j,a)/BL_rat(h,i,j,a);
                    
                            end
                        end
                    end

                end
            end
        end 
       
        avg(a)=sum_BL(a)/count(a);       
        stddev(a)=std(holder);      
        variance(a)=var(holder);
        stderror(a)=stddev(a)/sqrt(count(a));
        
        clear holder

    end

    cd(stat2);

    name=sprintf('Variance_AddCrit_Average_%s_%s_Ratios%s.txt',metabol,mask,source);

    fileID = fopen(name,'w');
    fprintf(fileID,'mean stddev stderror %s-RatiostCR %s\n',mask,source);
    for i=1:nb_blocks
        fprintf(fileID,'%d %d %d\n',avg(i),stddev(i),stderror(i));
    end
    fclose(fileID);

end

    name=sprintf('X_ALLDATA_%s_%s.txt',mask,source);
    fileID = fopen(name,'w');
    dlmwrite(name,big_holder)
    fclose(fileID);

end







