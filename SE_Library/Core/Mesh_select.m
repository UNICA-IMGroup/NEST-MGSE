function [best_mesh] = Mesh_select(path_mesh, br_to_nod)

b = 0;
c = 0;
for i=1:length(br_to_nod)
    a = find(path_mesh(:,1) == br_to_nod(i));
    a = max(a);
    c = [c; a];
    b = [b; a-c(end-1)];
end
b = b(2:end);
c = c(2:end);
num = length(b);
vect = 1:1:num-1;
somma = sum(vect);

best_mesh = zeros(2*(num-1),size(path_mesh,2));
num_nod_vect = 10^6*ones(num-1,1);
for i=1:num-1
    for j=i+1:num
        path1 = path_mesh(c(i)-b(i)+1:c(i),:);
        path2 = path_mesh(c(j)-b(j)+1:c(j),:);
        [opt_mesh, num_nod_path] = Mesh_optimal(path1, path2);
        idx = find(num_nod_path >= num_nod_vect);
        min_idx = find(num_nod_vect > num_nod_path,1);
        if min_idx < num-1
           num_nod_vect(min_idx+1:end) = num_nod_vect(min_idx:end-1);
           best_mesh(2*min_idx+1:end,:) = best_mesh(2*min_idx-1:end-2,:);
        end
        if length(idx) < length(num_nod_vect)
            num_nod_vect(min_idx) = num_nod_path;
            best_mesh(2*min_idx-1:2*min_idx,:) = opt_mesh;
        end
    end
end


        
















     