function Z = branch_impedance_per_length(config)
% input 
% - config: code of the type of branch
% output
% - Z: upper triangular representation of the impedance matrix
%  Z (R +jX) in ohms per mile

Z = [];
switch config
    case 1
        Z = [[0.4576 + 1i * 1.0780, 0.1560 + 1i * 0.5017, 0.1535 + 1i * 0.3849];
                                     [                   0, 0.4666 + 1i * 1.0482, 0.1580 + 1i * 0.4236];
                                     [                   0,                    0, 0.4615 + 1i * 1.0651]];
    
    case 2
        Z = [[0.4666 + 1i * 1.0482, 0.1580 + 1i * 0.4236, 0.1560 + 1i * 0.5017];
                                     [                   0, 0.4615 + 1i * 1.0651, 0.1535 + 1i * 0.3849];
                                     [                   0,                    0, 0.4576 + 1i * 1.0780]];
    
    case 3
        Z = [[0.4615 + 1i * 1.0651, 0.1535 + 1i * 0.3849, 0.1580 + 1i * 0.4236];
                                     [                   0, 0.4576 + 1i * 1.0780, 0.1560 + 1i * 0.5017];
                                     [                   0,                    0, 0.4666 + 1i * 1.0482]];
    
    case 4
        Z = [[0.4615 + 1i * 1.0651, 0.1580 + 1i * 0.4236, 0.1535 + 1i * 0.3849];
            [                   0, 0.4666 + 1i * 1.0482, 0.1560 + 1i * 0.5017];
            [                   0,                    0, 0.4576 + 1i * 1.0780]];
    
    case 5
        Z = [[0.4666 + 1i * 1.0482, 0.1560 + 1i * 0.5017, 0.1580 + 1i * 0.4236];
                                     [                   0, 0.4576 + 1i * 1.0780, 0.1535 + 1i * 0.3849];
                                     [                   0,                    0, 0.4615 + 1i * 1.0651]];
    
    case 6
        Z = [[0.4576 + 1i * 1.0780, 0.1535 + 1i * 0.3849, 0.1560 + 1i * 0.5017];
                                     [                   0, 0.4615 + 1i * 1.0651, 0.1580 + 1i * 0.4236];
                                     [                   0,                    0, 0.4666 + 1i * 1.0482]];
        
    case 7
        Z = [[0.4576 + 1i * 1.0780,                    0, 0.1535 + 1i * 0.3849];
                                     [                   0,                    0,                    0];
                                     [                   0,                    0, 0.4615 + 1i * 1.0651]];
    
    case 8
        Z = [[0.4576 + 1i * 1.0780, 0.1535 + 1i * 0.3849,                    0];
                                     [                   0, 0.4615 + 1i * 1.0651,                    0];
                                     [                   0,                    0,                    0]];
    
    case 9
        Z = [[1.3292 + 1i * 1.3475,                    0,                    0];
                                     [                   0,                    0,                    0];
                                     [                   0,                    0,                    0]];
    
    case 10
        Z = [[                   0,                    0,                    0];
                                     [                   0, 1.3292 + 1i * 1.3475,                    0];
                                     [                   0,                    0,                    0]];
    
    case 11
        Z = [[                   0,                    0,                    0];
                                     [                   0,                    0,                    0];
                                     [                   0,                    0, 1.3292 + 1i * 1.3475]];
    
    case 12
        Z = [[1.5209 + 1i * 0.7521, 0.5198 + 1i * 0.2775, 0.4924 + 1i * 0.2157];
                                     [                   0, 1.5329 + 1i * 0.7162, 0.5198 + 1i * 0.2775];
                                     [                   0,                    0, 1.5209 + 1i * 0.7521]];
                                
    
    %%%%%%%%
    case 601
        Z = [[0.3465 + 1i * 1.0179, 0.1560 + 1i * 0.5017, 0.1580 + 1i * 0.4236];
                                     [                    0, 0.3375 + 1i * 1.0478, 0.1535 + 1i * 0.3849];
                                     [                    0,                    0, 0.3414 + 1i * 1.0348]];
    
    case  602
        Z = [[0.7526 + 1i * 1.1814, 0.1580 + 1i * 0.4236, 0.1560 + 1i * 0.5017];
                                     [                   0, 0.7475 + 1i * 1.1983, 0.1535 + 1i * 0.3849];
                                     [                   0,                    0, 0.7436 + 1i * 1.2112]];

     
    case  603
        Z = [[                   0,                    0,                    0];
                                     [                   0, 1.3294 + 1i * 1.3471, 0.2066 + 1i * 0.4591];
                                     [                   0,                    0, 1.3238 + 1i * 1.3569]];
        
    
    case  604
        Z = [[1.3238 + 1i * 1.3569,                   0, 0.2066 + 1i * 0.4591];
                                     [                   0,                   0,                    0];
                                     [                   0,                   0, 1.3294 + 1i * 1.3471]];
    
    case  605
        Z = [[                   0,                   0,                    0];
                                     [                   0,                   0,                    0];
                                     [                   0,                   0, 1.3292 + 1i * 1.3475]];
        
    
    case  606
        Z = [[0.7982 + 1i * 0.4463, 0.3192 + 1i * 0.0328, 0.2849 + 1i * 0.0143];
                                     [                   0, 0.7891 + 1i * 0.4041, 0.3192 + 1i * 0.0328];
                                     [                   0,                    0, 0.7982 + 1i * 0.4463]];        
    
    case  607
        Z = [[1.3425 + 1i * 0.5124,                    0,                    0];
                                     [                   0,                    0,                    0];
                                     [                   0,                    0,                    0]];
                                

    otherwise
        disp([mfilename, ': config code ', num2str(config), ' not found']);
        
end