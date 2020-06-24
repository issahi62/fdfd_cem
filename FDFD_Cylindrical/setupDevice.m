function [ device] = setupDevice(fn,device,pml)   
    ER = importdata(strcat(fn,'ER_Layer1.dat'));
    
    %Apply PML ER
    device.ERxx = ER;
    device.ERyy = ER;
    device.ERzz = ER;
    
    %Form N^2 x M^2 matrices
    device.ERxx = diag(sparse(device.ERxx(:)));
    device.ERyy = diag(sparse(device.ERyy(:)));
    device.ERzz = diag(sparse(device.ERzz(:)));
    
    %Apply PML UR
    UR = importdata(strcat(fn,'UR_Layer1.dat'));
    device.URxx = UR;
    device.URyy = UR;
    device.URzz = UR;
    
    device.URxx = diag(sparse(device.URxx(:)));
    device.URyy = diag(sparse(device.URyy(:)));
    device.URzz = diag(sparse(device.URzz(:)));
end

