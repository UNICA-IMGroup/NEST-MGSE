function [opt_mesh, num_nod_path] = Mesh_optimal(path1, path2)

m = size(path1,1);
n = size(path2,1);
num_nod_path = 10^6;

for i=1:m
    a = 1;
    for j=1:n
        path_1st = path1(i,:);
        path_2nd = path2(j,:);
        var = 0;
        while var == 0
            com_node = find(path_2nd == path_1st(a));
            if size(com_node) == [1,0]
                var = 0;
                a = a+1;
            else
                num_nod_path2 = a + com_node;
                var =1;
            end
        end
        if num_nod_path2 < num_nod_path
            num_nod_path = num_nod_path2;
            way1 = path_1st(1:a);
            way2 = path_2nd(1:com_node);
        end
    end
end

opt_mesh = zeros(2,size(path1,2));
opt_mesh(1,1:length(way1)) = way1;
opt_mesh(2,1:length(way2)) = way2;

        