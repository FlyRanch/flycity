function [ p ] = test_for_clap(p,PAR)
% [ p ] = test_for_clap( p,PAR ) 
% p is the state vector and PAR contains the body morphology params
% The function checks the state to see if the wings are in a clap pose. 
% If so rotate the wings so that the wingtips are in the
% midsagital plane with the elevation taken from the state estimate.
% --TODO: edit the wing rotation so that the wing chord is parallel to the 
% plane.

    %First get the translation vector for the wings, to define the 
    %mid-sagital plane
    BL = PAR.params.bodyscale*(PAR.params.bodylen+PAR.params.headlen);
    RJTrans = BL.*([0.2021 0.1055 -0.1477]);
    LJTrans = BL.*([0.2021 -0.1055 -0.1477]);
    RJTrans = RJTrans .* 0.001*1.2;
    LJTrans = LJTrans .* 0.001*1.2;
    %now test for a clap pose - at this point just use the condition
    %that a wing has crossed the mid-sagital plane. Edit to apply a more 
    %rigorous condition.
    
    %left and right quarternions
    q_right = p(12:12+3);
    q_left = p(8:8+3);

    %Make vectors for left and right wingtip
    ltip = [0,PAR.params.winglen*PAR.params.wingscale*-1,0]';
    rtip = [0,PAR.params.winglen*PAR.params.wingscale,0]';
    
    %test for clap
    rmat_rw = quat2matNEW(q_right);
    rmat_lw = quat2matNEW(q_left);
    
    disp(rmat_lw*ltip);
    disp(LJTrans(2)*-1);
    disp('');
    disp(rmat_rw*rtip);
    disp(RJTrans(2)*-1);
    disp('');
    condition_left = ([0,1,0]*rmat_lw*ltip > LJTrans(2)*-1)
    condition_right = ([0,1,0]*rmat_rw*rtip < RJTrans(2)*-1)
    
    if condition_left
        ea = SpinConv('DCMtoEA213',rmat_lw,'1',1);
        %a = ea(1);
        a = -90;
        b = ea(2);
        c = 90;
        eal_prime = [a b c];
        rmatsl_prime = SpinConv('EA213toDCM',eal_prime,'1',1);
        ql_prime = quat2matNEW(rmatsl_prime);
        p(8:8+3) = ql_prime;
        disp('cond_left')
    end
    
    if condition_right
        ea = SpinConv('DCMtoEA213',rmat_rw,'1',1);
        %a = ea(1);
        a = -90;
        b = ea(2);
        c = -1*90;
        ear_prime = [a b c];
        rmatsr_prime = SpinConv('EA213toDCM',ear_prime,'1',1);
        qr_prime = quat2matNEW(rmatsr_prime);
        p(12:12+3) = qr_prime;
        disp('cond_right');
    end
end
