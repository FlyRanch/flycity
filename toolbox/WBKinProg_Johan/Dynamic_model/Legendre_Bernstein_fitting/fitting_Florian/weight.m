function [ WL1, WL2, WR1, WR2, WLR1, WLR2 ] = weight(a_fit_L1, a_fit_L2, a_fit_R1, a_fit_R2, n_pol)


    nr_wb = size(a_fit_L1,1);

    var_L1 = zeros(nr_wb,nr_wb,n_pol+1);
    
    weightL1 = zeros(n_pol+1,nr_wb);

    var_L2 = zeros(nr_wb,nr_wb,n_pol+1);
    
    weightL2 = zeros(n_pol+1,nr_wb);
    
    var_R1 = zeros(nr_wb,nr_wb,n_pol+1);
    
    weightR1 = zeros(n_pol+1,nr_wb);

    var_R2 = zeros(nr_wb,nr_wb,n_pol+1);
    
    weightR2 = zeros(n_pol+1,nr_wb);
    
    var_LR1 = zeros(nr_wb,nr_wb,n_pol+1);
    
    weightLR1 = zeros(n_pol+1,2*(nr_wb));
    
    var_LR2 = zeros(nr_wb,nr_wb,n_pol+1);
    
    weightLR2 = zeros(n_pol+1,2*(nr_wb));    
    
    
    for m = 1:(nr_wb)
            
        for n = 1:(nr_wb)
            
            if m ~= n
            
                var_L1(m,n,:) = var([a_fit_L1(:,m)' ; a_fit_L1(:,n)']);
                
                var_L2(m,n,:) = var([a_fit_L2(:,m)' ; a_fit_L2(:,n)']);

                var_R1(m,n,:) = var([a_fit_R1(:,m)' ; a_fit_R1(:,n)']);
                
                var_R2(m,n,:) = var([a_fit_R2(:,m)' ; a_fit_R2(:,n)']);
                
                var_LR1(m,n,:) = var([a_fit_L1(:,m)' ; a_fit_R1(:,n)']);
                
                var_LR2(m,n,:) = var([a_fit_L2(:,m)' ; a_fit_R2(:,n)']);
            end
            
        end
        
    end  

    for i = 1:(n_pol+1)
        
        weightL1(i,:) = sum(var_L1(:,:,i));
        
        weightL2(i,:) = sum(var_L2(:,:,i));
        
        weightR1(i,:) = sum(var_R1(:,:,i));
        
        weightR2(i,:) = sum(var_R2(:,:,i));    
        
        weightLR1(i,:) = sum([var_L1(:,:,i) var_LR1(:,:,i); var_LR1(:,:,i) var_R1(:,:,i)]);
        
        weightLR2(i,:) = sum([var_L2(:,:,i) var_LR2(:,:,i); var_LR2(:,:,i) var_R2(:,:,i)]);
        
    end
    
    WL1 = (1/sum(1./sum(weightL1.^2)))*(1./sum(weightL1.^2));
    
    WL2 = (1/sum(1./sum(weightL2.^2)))*(1./sum(weightL2.^2));
    
    WR1 = (1/sum(1./sum(weightR1.^2)))*(1./sum(weightR1.^2));
    
    WR2 = (1/sum(1./sum(weightR2.^2)))*(1./sum(weightR2.^2));
    
    WLR1 = (1/sum(1./sum(weightLR1.^2)))*(1./sum(weightLR1.^2));
    
    WLR2 = (1/sum(1./sum(weightLR2.^2)))*(1./sum(weightLR2.^2));
    
    
    
   
    
%     WL1 = (1/sum(1./sum(weightL1)))*(1./sum(weightL1));
%     
%     WL2 = (1/sum(1./sum(weightL2)))*(1./sum(weightL2));
%     
%     WR1 = (1/sum(1./sum(weightR1)))*(1./sum(weightR1));
%     
%     WR2 = (1/sum(1./sum(weightR2)))*(1./sum(weightR2));
%     
%     WLR1 = (1/sum(1./sum(weightLR1)))*(1./sum(weightLR1));
%     
%     WLR2 = (1/sum(1./sum(weightLR2)))*(1./sum(weightLR2));
    
end

