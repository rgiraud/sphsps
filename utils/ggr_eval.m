


function [ggr] = ggr_eval(lab_map)

lab_map = int32(lab_map);
[h,w]  = size(lab_map);
sp_ind = unique(lab_map(:))';
sp_nbr = length(sp_ind);

%Smf matrix
S_tab = int32(zeros(3*h+1,3*w+1));

% Size 3D to 2D
ss = round(sqrt(h*w/sp_nbr)/1.5);

src = 0;
sum_src = 0;
for k=sp_ind
    
    %S_k current superpixel
    S_k      = lab_map == k;
    
    %PANORAMIC
    %Handle split regions
    ll = bwlabel(S_k, 8);
    if (length(unique(ll))>=3) %Two separate blocks?
        
        %copy paste lab_map at frontier
        S_k(:,end+1:2*w) = S_k;
        
        %Take largest blob
        ll = bwlabel(S_k, 8);
        size_n = zeros(length(unique(ll(:))),1);
        for n=0:length(unique(ll(:)))-1
            size_n(n+1) = sum(sum(ll==n));
        end
        [~,max_n] = sort(size_n,'descend');
        S_k = ll==(max_n(2)-1);
        
    end
    
    [y,x] = find(S_k);
    
    if (length(unique(y(:)))>2 && length(unique(x(:)))>2)
        
        %Conversion 2D -> 3D -> 2D
        r = 1;
        xs = r*sin(y*pi/h).*cos(x*2*pi/w);
        ys = r*sin(y*pi/h).*sin(x*2*pi/w);
        zs = r*cos(y*pi/h);
        data=[xs';ys';zs'];
        %     data = repmat(data,[1 3]);
        
        x = my_own_pca(data);
        S_k = zeros(ss);
        x(:,1) = (x(:,1) - min(x(:,1)));
        x(:,2) = (x(:,2) - min(x(:,2)));
        max_1 = max(x(:,1));
        max_2 = max(x(:,2));
        maxx = max(max_1,max_2);
        for j=1:length(x(:,1))
            S_k(round(x(j,1)/maxx*ss+1),round(x(j,2)/maxx*ss+1)) = 1;
        end
        
        [yk,xk]  = find(S_k);
        size_S_k = sum(S_k(:));
        
        %Convex hull of S_k
        hull       = regionprops(S_k,'ConvexImage');
        hull       = hull.ConvexImage;
        perim_hull = regionprops(hull,'Perimeter');
        perim_hull = perim_hull.Perimeter;
        cc_hull    = perim_hull/sum(hull(:));
        
        perim_S_k = regionprops(S_k,'Perimeter');
        perim_S_k = perim_S_k.Perimeter;
        cc_S_k    = perim_S_k/size_S_k;
        
        if (cc_S_k ~= 0)
            
            %Evaluates the convexity of S_k
            cr_k = cc_hull/cc_S_k;
            
            %Evaluates the balanced repartition of S_k
            sigma_x = std(xk(:));
            sigma_y = std(yk(:));
            vxy_k   = sqrt(min(sigma_x,sigma_y)/max(sigma_x,sigma_y));
            
            %Shape Regularity Criteria (SRC)
            src_k = cr_k*vxy_k;
            src   = src + src_k*size_S_k;
            
            %SMF
            %Barycenter
            my = round(mean(yk(:)));
            mx = round(mean(xk(:)));
            
            for l=1:length(yk)
                %Registered position
                yk_r               = yk(l)+(1.5*h+1-my);
                xk_r               = xk(l)+(1.5*w+1-mx);
                S_tab(yk_r,xk_r) = S_tab(yk_r,xk_r) + 1;
            end
            
            sum_src = sum_src + size_S_k;
            
        end
        
    end
    
end
src = src/sum_src;

%Smooth Matching Factor (SMF)
S_mean = double(S_tab)/(sp_nbr);
S_mean = S_mean/sum(S_mean(:));

smf = 0;
sum_smf = 0;
for k=sp_ind
    
    %S_k current superpixel
    S_k      = lab_map == k;
    
    %PANORAMIC
    %Handle split regions
    ll = bwlabel(S_k, 8);
    if (length(unique(ll))>=3) %Two separate blocks?
        
        %copy paste lab_map at frontier
        S_k(:,end+1:2*w) = S_k;
        
        %Take largest blob
        ll = bwlabel(S_k, 8);
        size_n = zeros(length(unique(ll(:))),1);
        for n=0:length(unique(ll(:)))-1
            size_n(n+1) = sum(sum(ll==n));
        end
        [~,max_n] = sort(size_n,'descend');
        S_k = ll==(max_n(2)-1);
        
    end
    
    [y,x] = find(S_k);
    
    if (length(unique(y(:)))>2 && length(unique(x(:)))>2)
        
        %Conversion 2D -> 3D -> 2D
        r = 1;
        xs = r*sin(y*pi/h).*cos(x*2*pi/w);
        ys = r*sin(y*pi/h).*sin(x*2*pi/w);
        zs = r*cos(y*pi/h);
        data=[xs';ys';zs'];
        
        x = my_own_pca(data);
        
        S_k = zeros(ss);
        x(:,1) = (x(:,1) - min(x(:,1)));
        x(:,2) = (x(:,2) - min(x(:,2)));
        max_1 = max(x(:,1));
        max_2 = max(x(:,2));
        maxx = max(max_1,max_2);
        for j=1:length(x(:,1))
            S_k(round(x(j,1)/maxx*ss+1),round(x(j,2)/maxx*ss+1)) = 1;
        end
        
        [yk,xk]  = find(S_k);
        size_S_k = sum(S_k(:));
        
        %SMF
        %Barycenter
        my = round(mean(yk(:)));
        mx = round(mean(xk(:)));
        S_kk = zeros(3*h+1,3*w+1);
        
        for l=1:length(yk)
            %Registered position
            yk_r               = yk(l)+(1.5*h+1-my);
            xk_r               = xk(l)+(1.5*w+1-mx);
            S_kk(yk_r,xk_r) = S_kk(yk_r,xk_r) + 1;
        end
        
        size_S_k = sum(S_kk(:));
        S_kk      = S_kk/size_S_k;
        smf      = smf + size_S_k*sum(sum(abs(S_mean-S_kk)));
        
        sum_smf = sum_smf + size_S_k;
    end
end

smf = 1 - smf/(2*(sum_smf));


%%Generalizerd Global Regularity (G-GR) measure
ggr = src*smf;




function x = my_own_pca(data)


% remove the mean variable-wise (row-wise)
data=data-repmat(mean(data,2),1,size(data,2));

% calculate eigenvectors (loadings) W, and eigenvalues of the covariance matrix
[W, EvalueMatrix] = eig(cov(data'));
Evalues = diag(EvalueMatrix);

% order by largest eigenvalue
Evalues = Evalues(end:-1:1);
W = W(:,end:-1:1); W=W';

% generate PCA component space (PCA scores)
pc = W * data;

x = -pc';

end



end